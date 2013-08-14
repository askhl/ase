from lammpsbase import LAMMPSBase
import read_charmm_params
import typing, charmmtypes
import warnings

class CHARMM(LAMMPSBase):
	
	def __init__(self, ff_file_path, label='charmm', pair_cutoff=10.0, auto_charges=False, **kwargs):
		LAMMPSBase.__init__(self, label, **kwargs)
		
		self.parameters.units          = 'real'
		self.parameters.pair_style     = 'lj/charmm/coul/charmm 8.0 10.0'
		self.parameters.bond_style     = 'harmonic'
		self.parameters.angle_style    = 'charmm'
		self.parameters.dihedral_style = 'charmm'
		self.parameters.improper_style = 'harmonic'
		self.parameters.special_bonds  = 'charmm'
		
		# Read force field file
		self.ff_data = read_charmm_params.read(open(ff_file_path))
		self.type_resolver = typing.TypeResolver(charmmtypes.data)
		
		self.auto_charges = auto_charges
	
	def atom_types(self, atoms):
		self.type_resolver.resolve_atoms(atoms)
		return atoms.info['atom_types']
		
	def determine_charges(self, atoms, atom_types):
		if self.auto_charges:
			return charges_by_increments(atoms)
		else:
			return atoms.get_initial_charges()
		
	def prepare_calculation(self, atoms, data):
		if not self.auto_charges and (atoms.get_initial_charges() == 0).all():
			warnings.warn("No partial charges set! There won't be any Coulomb interactions.")


# CHARMM charges taken from top_all36_cgenff.rtf. The bond increment model seems
# to work with molecules containing C, H and O. With nitrogen it gets tricky.

def charges_by_increments(atoms):
	
	bond_increments = {
		('HGP1', 'OG311'):  .42,
		('HGA1', 'CG311'):  .09,
		('HGA2', 'CG321'):  .09,
		('HGA3', 'CG331'):  .09,
		('CG311', 'OG311'): .23,
		('CG321', 'OG311'): .23,
		('HGTIP3', 'OGTIP3'): 0.417,   # water
		('HGR61', 'CG2R61'): 0.115,    # benzene
		('CG2R61', 'OG311'): .11       # phenol
		}
	
	charges = atoms.get_initial_charges()
	types = atoms.info['atom_types']
	for i,j in atoms.info['bonds']:
		if types[i] == types[j]: continue
		if not (types[i], types[j]) in bond_increments:
			i,j = j,i
		incr = bond_increments.get((types[i], types[j]))
		if not incr:
			return atoms.get_initial_charges()
		charges[i] += incr
		charges[j] -= incr
	return charges