import os

from ase.test import NotAvailable

try:
    from ase.test.jacapo_installed import jacapo_installed
    jacapo_installed()
except NotAvailable:
    raise NotAvailable('Jacapo required')

from ase import Atoms, Atom
from ase.calculators.jacapo import Jacapo

atoms = Atoms([Atom('H',[0,0,0])],
            cell=(2,2,2))

calc = Jacapo('Jacapo-test.nc',
              pw=200,
              nbands=2,
              kpts=(1,1,1),
              spinpol=False,
              dipole=False,
              symmetry=False,
              ft=0.01)

atoms.set_calculator(calc)

print atoms.get_potential_energy()
os.system('rm -f Jacapo-test.nc Jacapo-test.txt')
