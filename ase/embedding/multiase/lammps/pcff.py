from lammpsbase import LAMMPSBase
import read_compass_params
import typing, pcfftypes
import numpy as np

class PCFF(LAMMPSBase):
	
	def __init__(self, ff_file_path, label='pcff', pair_cutoff=10.0, kspace=False, **kwargs):
		LAMMPSBase.__init__(self, label, **kwargs)
		
		self.parameters.units          = 'real'
		
		if kspace:
			self.parameters.pair_style     = 'lj/class2/coul/long %f' % pair_cutoff
			self.parameters.kspace_style   = 'ewald 0.0001'
		else:
			self.parameters.pair_style     = 'lj/class2/coul/cut %f' % pair_cutoff
			
		self.parameters.bond_style     = 'class2'
		self.parameters.angle_style    = 'class2'
		self.parameters.dihedral_style = 'class2'
		self.parameters.improper_style = 'class2'
		self.parameters.special_bonds  =  'lj/coul 0.0 0.0 1.0 dihedral yes'
		
		self.ff_data, self.bond_increments = read_compass_params.read(open(ff_file_path))
		self.type_resolver = typing.TypeResolver(pcfftypes.data)
	
	def atom_types(self, atoms):
		self.type_resolver.resolve_atoms(atoms)
		return atoms.info['atom_types']
	
	def determine_charges(self, atoms, atom_types):
		eq = self.ff_data.equivalence
		charges = np.zeros(len(atoms))
		
		for i, j in atoms.info['bonds']:
			t1 = eq[atom_types[i]]['bond']
			t2 = eq[atom_types[j]]['bond']
			if not (t1, t2) in self.bond_increments:
				i,j = j,i
				t1,t2 = t2,t1
			increment = self.bond_increments.get((t1, t2))
			if not increment:
				raise RuntimeError('No bond increment for (%s, %s)' % (t1, t2))
			charges[i] += increment[0]
			charges[j] += increment[1]
		return charges
		