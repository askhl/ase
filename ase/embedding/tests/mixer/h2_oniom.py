from ase import Atoms
from gpaw import GPAW
from ase.embedding.multiase.lammps.reaxff import ReaxFF
from ase.embedding.multiase.utils import get_datafile


import numpy as np

d = 0.76470
a = 6.0

atoms = Atoms("H2",
                positions = [(0, 0, 0),
                (0, 0, d)],
                cell = (10*a, 10*a, 10*a))

calc_gpaw = GPAW(nbands=2, txt="h2_1.txt")
calc_reaxff = ReaxFF(ff_file_path=get_datafile("ffield.reax.new"),
                implementation="C")
reaxff_cell = (100*a, 100*a, 100*a)
gpaw_cell = (a, a, a)


from ase.embedding.oniom import ONIOM
oniom = ONIOM(atoms, calc_full=calc_reaxff, cell_full=reaxff_cell,
              calc_loc=calc_gpaw, cell_loc=gpaw_cell)

atoms.set_calculator(oniom)
print(atoms.get_forces())
print(atoms.get_potential_energy())
