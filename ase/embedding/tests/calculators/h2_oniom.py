from ase import Atoms
from ase.embedding.calculators.oniom import ONIOM
#import ase calculators for oniom type of embedding
from gpaw import GPAW
from ase.calculators.emt import EMT
from ase.embedding.multiase.lammps.reaxff import ReaxFF
from ase.embedding.multiase.utils import get_datafile


import numpy as np
d = 0.76470
a = 6.0

#create atoms
atoms = Atoms("H2",
              positions = [(0, 0, 0),
                           (0, 0, d)],
              cell = (10*a, 10*a, 10*a))

#Full region calculator and cell (also noted as MM region)
calc_emt = EMT()
cell_emt = (10*a,10*a,10*a)

#Alternative full region calculator and cell (also noted as MM region)
#to use them, replace replace the respective parameters in ONIOM()
calc_reaxff = ReaxFF(ff_file_path=get_datafile("ffield.reax.new"),
                     implementation="C")
cell_reaxff = (10*a, 10*a, 10*a)

#Local region calculator and cell (also noted as QM region)
calc_gpaw = GPAW(nbands=2, txt="h2_1.txt")
cell_gpaw = (a, a, a)
#list of atoms belonging to the localized region
atoms_gpaw = [0]

#create oniom calculator
oniom = ONIOM(atoms, calc_full=calc_emt, cell_full=cell_emt,
              calc_loc=calc_gpaw, cell_loc=cell_gpaw, atoms_loc=atoms_gpaw)

#the created calculator can be used in the same way as any other ase calculator
atoms.set_calculator(oniom)
print(atoms.get_forces())
print(atoms.get_potential_energy())
