from ase.calculators.interface import Calculator
import os, sys
import shutil, shlex
from subprocess import Popen, PIPE
from tempfile import mkdtemp, NamedTemporaryFile
import numpy as np
import unitconversion
from bonds import Bonds
from ffdata import FFData, SequenceType, ImproperType
from ase.io.lammpsrun import read_lammps_dump
from lammpspython import lammps
import warnings

from itertools import combinations, permutations

# "End mark" used to indicate that the calculation is done
CALCULATION_END_MARK = '__end_of_ase_invoked_calculation__'

class LAMMPSParameters:
	def __init__(self, **pars):
		self.atom_style     = pars.get('atom_style', 'full')
		self.units          = pars.get('units', 'real')
		self.neighbor       = pars.get('neighbor')
		self.newton         = pars.get('newton')
		
		self.pair_style     = pars.get('pair_style')
		self.bond_style     = pars.get('bond_style')
		self.angle_style    = pars.get('angle_style')        
		self.dihedral_style = pars.get('dihedral_style')
		self.improper_style = pars.get('improper_style')
		self.kspace_style   = pars.get('kspace_style')
		self.special_bonds  = pars.get('special_bonds')
		self.pair_modify    = pars.get('pair_modify')
		
		# Pair coeffs to be specified in the input file
		self.pair_coeffs    = pars.get('pair_coeffs', [])
		
		self.extra_cmds     = pars.get('extra_cmds', [])
		
		# Atoms groups: tuples (group_id, list_of_indices)
		self.groups         = pars.get('groups', [])

class LAMMPSData:
	def __init__(self):
		self.clear()
		
	def clear(self):
		self.tables = None
		
		self.atom_types = None
		self.bonds = None
		self.angles = None
		self.dihedrals = None
		self.impropers = None
		
		self.atom_typeorder = None
		self.bond_typeorder = None
		self.angle_typeorder = None
		self.dihedral_typeorder = None
		self.improper_typeorder = None


class LAMMPSBase(Calculator):

	def __init__(self, label='lammps', tmp_dir=None, parameters={}, 
		     update_charges=False, force_triclinic=False,
		     keep_alive=False, debug=False, mpi_comm=None):
		"""The LAMMPS calculators object """

		self.label = label
		self.parameters = LAMMPSParameters(**parameters)
		self.data = LAMMPSData()
		self.ff_data = None
		self.forces = None
		self.atoms = None
		self.atoms_after_last_calc = None
		self.update_charges = update_charges
		self.force_triclinic = force_triclinic
		self.keep_alive = keep_alive
		self.debug = debug or tmp_dir     
		self.lammps_process = LammpsLibrary(log=self.debug, mpi_comm=mpi_comm)
		self.calls = 0
		self.mpi_comm = mpi_comm
		
		self._custom_thermo_args = [
			'step', 'temp', 'press', 'cpu', 
			'pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz',
			'ke', 'pe', 'etotal',
			'vol', 'lx', 'ly', 'lz', 'atoms']
		
		self._dump_fields = [
			'id', 'type', 'x', 'y', 'z', 'vx',
			'vy', 'vz', 'fx', 'fy', 'fz', 'q']
		
		if tmp_dir is None:
			self.tmp_dir = self.single_call(mkdtemp, prefix='LAMMPS-')
		else:
			# If tmp_dir is pointing somewhere, don't remove stuff!
			self.tmp_dir=os.path.realpath(tmp_dir)
			if not os.path.isdir(self.tmp_dir):
				self.single_call(os.mkdir, self.tmp_dir, 0755)
		
		if self.debug:
			print 'LAMMPS (label: %s) running at %s' % (self.label, self.tmp_dir)
	
	def single_call(self, func, *args, **kwargs):
		if self.mpi_comm:
			if self.mpi_comm.rank == 0:
				val = func(*args, **kwargs)
			else:
				val = None
			return self.mpi_comm.bcast(val, root=0)
		else:
			return func(*args, **kwargs)
	
	def __del__(self):
		if not self.debug:
			shutil.rmtree(self.tmp_dir)
	
	def atom_types(self, atoms):
		""" Implement this method in subclasses"""
		raise NotImplementedError()
	
	def determine_charges(self, atoms, atom_types):
		""" Reimplement this method in subclasses if needed"""
		return atoms.get_initial_charges()
	
	def prepare_calculation(self, atoms, data):
		""" Implement this method in subclasses if needed"""
		pass

	def get_potential_energy(self, atoms):
		self.update(atoms)
		energy = self.lammps_process.get_thermo('pe')
		return self.to_ase_units(energy, 'energy')

	def get_forces(self, atoms):
		self.update(atoms)
		return self.forces

	def get_stress(self, atoms):
		self.update(atoms)
		thermo = self.lammps_process.get_thermo
		stress = [thermo(key) for key in ('pxx','pyy','pzz', 'pyz','pxz','pxy')]
		return self.to_ase_units(np.array(stress), 'stress')  
	
	def calculation_required(self, atoms, quantities=None):
		return atoms != self.atoms_after_last_calc or \
			(atoms.get_initial_charges() != self.atoms_after_last_calc.get_initial_charges()).any()
		
	def update(self, atoms):
		if not self.calculation_required(atoms): return
		self.setup_calculation(atoms)
		self.evaluate_forces()
		self.close_calculation()
	
	def minimize(self, atoms, etol=0, ftol=0, maxeval=100000, maxiter=100000000, min_style=None, relax_cell=False):
		if etol == 0 and ftol == 0: 
			raise RuntimeError('Specify at least one tolerance value!')
		ftol = self.from_ase_units(ftol, 'force')
		minimize_params = '%g %g %i %i' % (etol, ftol, maxeval, maxiter)
		
		self.setup_calculation(atoms)
		self.evaluate_forces()
		
		f = self.lammps_process
		if relax_cell:
			f.write('fix relax_cell all box/relax tri 0.0 nreset 20\n')
		if min_style:
			f.write('min_style %s\n' % min_style)
		f.write('minimize %s\n'  % minimize_params)
		if relax_cell:
			f.write('unfix relax_cell\n')
		self.update_forces()
		self.update_atoms()
		self.close_calculation()
		
	def molecular_dynamics(self, atoms, timestep, fix, step_iterator, update_cell, total_steps, constraints):
		timestep = str(self.from_ase_units(timestep, 'time'))
		self.setup_calculation(atoms)
		self.evaluate_forces()
		cur_step = 0
		
		f = self.lammps_process
		f.write('fix mdfix %s\n' % fix)
		f.write('timestep %s\n' % timestep)
		for c in constraints:
			for cmd in c.get_commands(self.atoms):
				f.write(cmd + '\n')
		
		for nsteps in step_iterator:
			if total_steps:
				f.write('run %s start 0 stop %s\n' % (nsteps, total_steps))
			else:
				f.write('run %s\n' % nsteps)
			self.update_forces()
			self.update_atoms(update_cell=update_cell)
			
			cur_step += nsteps
			
		self.close_calculation()
	
	def evaluate_forces(self):
		f = self.lammps_process
		f.write('run 0\n')
		self.update_forces()
		self.atoms_after_last_calc = self.atoms.copy()
	
	def setup_calculation(self, atoms):
		self.filelabel = '%s%06d' % (self.label, self.calls)
		self.lammps_process.start(self.tmp_dir, self.filelabel)
				
		if np.all(atoms.pbc == False):
			# Make sure the atoms are inside the cell
			inv_cell = np.linalg.inv(atoms.cell)
			frac_positions = np.dot(inv_cell, atoms.positions.T)
			if np.any(frac_positions < 0) or np.any(frac_positions > 1):
				atoms.center(vacuum=1)

		self.atoms = atoms
		self.prism = Prism(self.atoms.cell)
		self.prepare_data()
		self.prepare_calculation(self.atoms, self.data)
		self.write_lammps_input()
		if self.debug: print 'Calculation initialized.'
		self.calls += 1
	
	def close_calculation(self):
		if self.debug == True: self.lammps_process.close_logs()
	
	def prepare_data(self):
		""" Prepare self.data for write_lammps_data() using self.ff_data """
		
		atoms = self.atoms
		tables = []
		ff_data = self.ff_data
		
		def status_message(msg):
			if not self.debug: return
			print msg,
			sys.stdout.flush()
			
		def status_done():
			if not self.debug: return
			print 'Done.'
		
		if not 'bonds' in atoms.info:
			status_message('Detecting bonds...')
			atoms.info['bonds'] = Bonds(atoms, autodetect=True)
			status_done()
		
		status_message('Detecting atom types...')
		atom_types = self.atom_types(atoms)
		status_done()
		
		# Note: while often there is 1-1 mapping between formal atom types and
		# the nonbonded (actual) types, this is not the case for all FF's.
		# In FFData, the category 'atom' has data for each nonbonded (actual) type, instead of
		# for each formal type like one might think.
		atom_typeorder = list(set(atom_types))
		nonbonded_types = [ff_data.get_actual_type('atom', tp) for tp in atom_types]
		nonbonded_typeorder = [ff_data.get_actual_type('atom', tp) for tp in atom_typeorder]
		
		# A helper to run FFData.find() on a given set of objects of a given category ('bond', 'angle', etc...)
		# If parameters in FFData are not found, discard the object.
		def identify_objects(objects, category):
			result = []
			discarded_objects = set()
			for indices in objects:
				if category != 'improper':
					type = SequenceType([atom_types[i] for i in indices])
				else:
					a, c, b, d = (atom_types[ind] for ind in indices)
					type = ImproperType(central_type=a, other_types=(c,b,d))
					
				actual_indices, actual_type = ff_data.find(category, indices, type)
				if actual_indices == None:
					discarded_objects.add(actual_type)
					continue
				if category == 'improper' and ff_data.class2:
					actual_indices = [actual_indices[1], actual_indices[0], actual_indices[2], actual_indices[3]]
				result.append(dict(indices=actual_indices, type=actual_type))
				
			if self.debug:
				for tp in discarded_objects:
					print 'No parameters for %s. Skipping.' % tp
			return result
		
		b = atoms.info['bonds']
		parameters = self.parameters
		if parameters.bond_style:    
			bonds     = identify_objects(b, 'bond')
		else: bonds = []
		if parameters.angle_style:
			status_message('Detecting angles...')
			angles    = identify_objects(b.find_angles(), 'angle')
			status_done()
		else: angles = []
		if parameters.dihedral_style:
			status_message('Detecting dihedrals...')
			dihedrals = identify_objects(b.find_dihedrals(), 'dihedral')
			status_done()
		else: dihedrals = []
		if parameters.improper_style:
			status_message('Detecting impropers...')
			impropers = identify_objects(b.find_impropers(), 'improper')
			status_done()
		else: impropers = []
		
		# Compile a list of force field parameters for a set of objects of a given category
		# each set of parameters is added to table as a list. The return value is an ordered
		# list of all types present in objects.
		def add_coeff_tables(category, objects, typeorder=None, warn_missing=True):
			if not objects: return
			if typeorder:
				used_types = typeorder
			else:
				used_types = set(object['type'] for object in objects)
			
			available_tables = ff_data.available_tables(category)
			new_tables = {}
			for type in used_types:
				params = ff_data.get_params(category, type)
				for title, ncols in available_tables:
					try:
						values = params[title]
					except KeyError:
						if warn_missing: print 'No %s for %s!' % (title, type)
						values = [0]*ncols
					table = new_tables.setdefault(title, [])
					if self.debug:
						comment = '       # %s' % type
						values = values + [comment]
					table.append(values)
			tables.extend(new_tables.items())
			return list(used_types)
		
		# Add masses to ff_data
		masses = dict(zip(nonbonded_types, self.atoms.get_masses()))
		for type in nonbonded_typeorder:
			ff_data.add('atom', type, 'Masses', [masses[type]])
			
		add_coeff_tables('atom', nonbonded_types, nonbonded_typeorder)
		
		bond_typeorder     = add_coeff_tables('bond', bonds, warn_missing=self.debug)
		angle_typeorder    = add_coeff_tables('angle', angles, warn_missing=self.debug)
		dihedral_typeorder = add_coeff_tables('dihedral', dihedrals, warn_missing=self.debug)
		improper_typeorder = add_coeff_tables('improper', impropers, warn_missing=self.debug)
		
		# Atoms
		charges = self.determine_charges(atoms, atom_types)
		atom_typeids = [nonbonded_typeorder.index(at)+1 for at in nonbonded_types]
		positions = self.prism.vector_to_lammps(self.atoms.positions)
		positions = self.from_ase_units(positions, 'distance')
		columns = [atom_typeids, charges, positions[:,0], positions[:,1], positions[:,2]]
		if self.debug:
			comments = ['   # %s' % tp for tp in atom_types]
			columns += [comments]
		
		if self.parameters.atom_style == 'full':
			columns.insert(0, ['1']*len(self.atoms))
		elif not self.parameters.atom_style == 'charge':
			raise RuntimeError('Unsupported atom_style: %s' % self.parameters.atom_style)
			
		tables.append(('Atoms', zip(*columns)))
		
		# Bonds, Angles, etc.
		def add_object_table(title, objects, typeorder):
			if not objects or not typeorder: return
			table = []
			for obj in objects:
				typeid = typeorder.index(obj['type'])+1
				atoms = [idx+1 for idx in obj['indices']]
				values = [typeid] + atoms
				if self.debug:
					comment = '    # %s' % obj['type']
					values += [comment]
				table.append(values)
			tables.append((title, table))
		
		add_object_table('Bonds', bonds, bond_typeorder)
		add_object_table('Angles', angles, angle_typeorder)
		add_object_table('Dihedrals', dihedrals, dihedral_typeorder)
		add_object_table('Impropers', impropers, improper_typeorder)
		
		if self.atoms.has('momenta'):
			vel = self.prism.vector_to_lammps(self.atoms.get_velocities())
			lammps_velocities = self.from_ase_units(vel, 'velocity')
			tables.append(('Velocities', lammps_velocities))
		
		data = self.data
		data.tables = tables
		data.atom_types = atom_types
		data.bonds      = bonds
		data.angles     = angles
		data.dihedrals  = dihedrals
		data.impropers  = impropers
		data.atom_typeorder     = atom_typeorder
		data.bond_typeorder     = bond_typeorder
		data.angle_typeorder    = angle_typeorder
		data.dihedral_typeorder = dihedral_typeorder
		data.improper_typeorder = improper_typeorder
		
	
	def write_lammps_input(self):
		"""Write LAMMPS parameters and input data """
		f = self.lammps_process
		parameters = self.parameters
		
		# Write datafile
		datafile = self.single_call(self.write_lammps_data)
		
		f.write('# (written by ASE)\n')
		f.write('clear\n')
		f.write('atom_style %s \n' % parameters.atom_style)
		
		pbc = self.atoms.get_pbc()
		f.write('units %s \n' % parameters.units)
		f.write('boundary %s %s %s \n' % tuple('mp'[int(x)] for x in pbc))
		if parameters.neighbor:
			f.write('neighbor %s \n' % (parameters.neighbor))
		if parameters.newton:
			f.write('newton %s \n' % (parameters.newton))

		# Write interaction stuff
		f.write('\n### interactions \n')
		if parameters.pair_style:
			f.write('pair_style %s \n' % parameters.pair_style)
			
		if parameters.bond_style and self.data.bonds:
			f.write('bond_style %s \n' % parameters.bond_style)
			
		if parameters.angle_style and self.data.angles:
			f.write('angle_style %s \n' % parameters.angle_style)
			
		if parameters.dihedral_style and self.data.dihedrals:
			f.write('dihedral_style %s \n' % parameters.dihedral_style)
			
		if parameters.improper_style and self.data.impropers:
			f.write('improper_style %s \n' % parameters.improper_style)
		
		if parameters.kspace_style:
			f.write('kspace_style %s \n' % parameters.kspace_style)
		
		if parameters.special_bonds:
			f.write('special_bonds %s \n' % parameters.special_bonds)
			
		if parameters.pair_modify:
			f.write('pair_modify %s \n' % parameters.pair_modify)
		
		print 'read data'
		f.write('\n### read data \n')
		f.write('read_data %s\n' % datafile)
		print 'data read'
		
		# Extra pair coeffs
		for line in parameters.pair_coeffs:
			f.write('pair_coeff %s \n' % line)
		
		# Create groups
		for group_id, indices in parameters.groups:
			indices_str = ' '.join([str(i+1) for i in indices])
			f.write('group %s id %s\n' % (group_id, indices_str))
		
		for cmd in parameters.extra_cmds:
			f.write(cmd + '\n')
		
		
	def write_lammps_data(self):
		"""Write system configuration and force field parameters to file to be read
		with read_data by LAMMPS."""
		prefix = 'data_%s' % self.filelabel
		f = NamedTemporaryFile(prefix=prefix, dir=self.tmp_dir, delete=(not self.debug))
		filename = f.name
		data = self.data
		
		f.write(filename + ' (written by ASE) \n\n')

		f.write('%d \t atoms \n' % len(data.atom_types))
		if data.bonds:     f.write('%d \t bonds \n' % len(data.bonds))
		if data.angles:    f.write('%d \t angles \n' % len(data.angles))
		if data.dihedrals: f.write('%d \t dihedrals \n' % len(data.dihedrals))
		if data.impropers: f.write('%d \t impropers \n' % len(data.impropers))
		
		f.write('%d  atom types\n' % len(data.atom_typeorder))
		if data.bonds:     f.write('%d  bond types\n' % len(data.bond_typeorder))
		if data.angles:    f.write('%d  angle types\n' % len(data.angle_typeorder))
		if data.dihedrals: f.write('%d  dihedral types\n' % len(data.dihedral_typeorder))
		if data.impropers: f.write('%d  improper types\n' % len(data.improper_typeorder))

		xhi, yhi, zhi, xy, xz, yz = self.prism.get_lammps_prism()
		f.write('0.0 %f  xlo xhi\n' % xhi)
		f.write('0.0 %f  ylo yhi\n' % yhi)
		f.write('0.0 %f  zlo zhi\n' % zhi)
		
		if self.force_triclinic or self.prism.is_skewed():
			f.write('%f %f %f  xy xz yz\n' % (xy, xz, yz))
		f.write('\n\n')
		
		for title, table in data.tables:
			if len(table) == 0:  continue
			f.write('%s \n\n' % title)
			for index, row in enumerate(table):
				f.write(('%d'+' %s'*len(row) +'\n') % ((index+1,) + tuple(row)))
			f.write('\n\n')
		
		f.close()
		return filename

	def update_forces(self):
		f = self.lammps_process.get_forces()
		self.forces = self.prism.vector_to_ase(self.to_ase_units(f, 'force'))
		
	def update_atoms(self, update_cell=False):
		p = self.lammps_process.get_positions()
		v = self.lammps_process.get_velocities()
		
		if update_cell:
			cell, celldisp = self.lammps_process.get_cell()
			p -= celldisp
			asecell = self.prism.vector_to_ase(self.to_ase_units(cell, 'distance'))
			self.atoms.cell = asecell
		
		pos = self.prism.vector_to_ase(self.to_ase_units(p, 'distance'))
		vel = self.prism.vector_to_ase(self.to_ase_units(v, 'velocity'))	
		self.atoms.set_positions(pos)
		self.atoms.set_velocities(vel)
		
		if self.update_charges:
			self.atoms.set_initial_charges(self.lammps_process.get_charges())
			
		self.atoms_after_last_calc = self.atoms.copy()
	
	def to_ase_units(self, value, quantity):
		return unitconversion.convert(value, quantity, self.parameters.units, 'ASE')

	def from_ase_units(self, value, quantity):
		return unitconversion.convert(value, quantity, 'ASE', self.parameters.units)


class LammpsLibrary:
	''' Interface to the LAMMPS library '''
	def __init__(self, log=False, mpi_comm=None):
		self.lammps = None
		self.inlog = None
		self.log = log
		
		if mpi_comm and mpi_comm.rank > 0:
			self.log = False
	
	def start(self, tmp_dir, filelabel=''):
		if self.lammps: self.lammps.close()
		self.close_logs()
		
		if self.log == True:
			outlog = NamedTemporaryFile(prefix='log_'+filelabel, dir=tmp_dir, delete=False)
			self.inlog  = NamedTemporaryFile(prefix='in_'+filelabel, dir=tmp_dir, delete=False)
			outlogpath = outlog.name
			outlog.close()
		else:
			outlogpath = 'none'
		
		self.lammps = lammps.lammps(cmdargs=['-log', outlogpath, '-screen', 'none'])

	def running(self):
		return self.lammps != None
	
	def write(self, command):
		if self.inlog: self.inlog.write(command)
		self.lammps.command(command)
		
	def close_logs(self):
		if self.inlog: self.inlog.close()
		self.inlog = None
		
	def get_thermo(self, key):
		stress_components = ['pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz']
		if key in stress_components:
			stress = self.lammps.extract_compute('thermo_press', 0, 1)
			return stress[stress_components.index(key)]
		else:
			return self.lammps.extract_compute('thermo_%s' % key, 0, 0)

	def get_positions(self):
		pos = self.lammps.gather_atoms('x', 1, 3)
		return np.reshape(list(pos), (-1,3))
		
	def get_velocities(self):
		vel = self.lammps.gather_atoms('v', 1, 3)
		return np.reshape(list(vel), (-1,3))
		
	def get_charges(self):
		q = self.lammps.gather_atoms('q', 1, 1)
		return np.array(list(q))
	
	def get_cell(self):
		names = 'boxxlo','boxxhi','boxylo','boxyhi','boxzlo','boxzhi'
		values = [self.lammps.extract_global(name, 1) for name in names]
		xlo, xhi, ylo, yhi, zlo, zhi = values
		
		# Hack: extract_global() doesn't support tilt factors
		self.write('variable boxxy equal xy\n')
		self.write('variable boxyz equal yz\n')
		self.write('variable boxxz equal xz\n')
		xy = self.lammps.extract_variable('boxxy', None, 0)
		yz = self.lammps.extract_variable('boxyz', None, 0)
		xz = self.lammps.extract_variable('boxxz', None, 0)
		
		xhilo = (xhi - xlo) - abs(xy) - abs(xz)
		yhilo = (yhi - ylo) - abs(yz)
		zhilo = (zhi - zlo)
		celldispx = xlo - min(0, xy) - min(0, xz)
		celldispy = ylo - min(0, yz)
		celldispz = zlo

		cell = [[xhilo, 0, 0], [xy, yhilo, 0], [xz, yz, zhilo]]
		celldisp = [[celldispx, celldispy, celldispz]]
		
		return cell, celldisp

	def get_forces(self):
		force = self.lammps.gather_atoms('f', 1, 3)
		return np.reshape(list(force), (-1, 3))

class Prism:
	"""The representation of the unit cell in LAMMPS"""

	def __init__(self, cell, pbc=(True,True,True), digits=10):
		# Use LQ decomposition to get the lammps cell
		# ase_cell * R = lammps_cell
		Qtrans, Ltrans = np.linalg.qr(cell.T)
		self.R = Qtrans
		self.lammps_cell = Ltrans.T
	
		if self.is_skewed() and not np.all(pbc):
			raise RuntimeError('Skewed lammps cells MUST have '
							'PBC == True in all directions!')

	def get_lammps_prism(self):
		return self.lammps_cell[(0,1,2,1,2,2), (0,1,2,0,0,1)]
		
	def update_cell(self, xyz, offdiag):
		self.lammps_cell = self.to_cell_matrix(xyz, offdiag)
		return np.dot(self.lammps_cell, self.R.T)
	
	def to_cell_matrix(self, xyz, offdiag):
		x, y, z = xyz
		xy, xz, yz = offdiag
		return np.array([[x,0,0], [xy,y,0], [xz, yz, z]])

	def vector_to_lammps(self, vec):
		return np.dot(vec, self.R)
		
	def vector_to_ase(self, vec):
		return np.dot(vec, self.R.T)
	
	def is_skewed(self):
		tolerance = 1e-6
		cell_sq = self.lammps_cell**2
		return np.sum(np.tril(cell_sq, -1)) / np.sum(np.diag(cell_sq)) > tolerance


