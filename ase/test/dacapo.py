from ase import *
from ase.calculators.dacapo  import Dacapo

if 1:
    h2 = Atoms('H2', positions=[(0, 0, 0), (0, 0, 1.1)],
               calculator=Dacapo(), pbc=1)
    h2.center(vacuum=2.0)
else:
    h2 = Dacapo().get_atoms()

if 0:
    print h2.get_potential_energy()
