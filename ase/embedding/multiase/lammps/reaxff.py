from lammpsbase import LAMMPSBase, FFData
import numpy as np

class ReaxFF(LAMMPSBase):
	
	def __init__(self, ff_file_path, label='reaxff', update_charges=True, **kwargs):
		LAMMPSBase.__init__(self, label, update_charges=update_charges, **kwargs)
		
		self.parameters.atom_style = 'charge'
		self.parameters.pair_style = 'reax/c NULL'
		self.parameters.extra_cmds.append('fix qeq all qeq/reax 1 0.0 10.0 1.0e-6 reax/c')
		self.parameters.units      = 'real'
				
		self.ff_file = ff_file_path
		self.ff_data = FFData()
	
	def atom_types(self, atoms):
		return atoms.get_chemical_symbols()
	
	def prepare_calculation(self, atoms, data):
		typeorder = data.atom_typeorder
		self.parameters.pair_coeffs = ['* * %s %s' % (self.ff_file, ' '.join(typeorder))]
	