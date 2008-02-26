import numpy as npy
from ase import *

a = 3.6
b = a / 2
cu = Atoms('Cu2Ag',
           positions=[(0, 0, 0),
                      (b, b, 0),
                      (a, a, b)],
           calculator=EMT())
e0 = cu.get_potential_energy()
print e0
