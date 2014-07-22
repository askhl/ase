from ase.test import cli, require
from ase.db import connect
from ase.db.jsondb import read_json
from ase.calculators.gaussian import Gaussian

require('gaussian')
cli("""\
ase-build O | ase-run gaussian -d oxygen.json &&
ase-build O2 | ase-run gaussian -d oxygen.json""")
c = connect('oxygen.json')
dct = read_json('oxygen.json')
for name in ['O2', 'O']:
    d = c.get([('name', '=', name)])
    id = d.id
    e1 = d.energy
    e2 = c.get_atoms(id).get_potential_energy()
    e3 = Gaussian.read_atoms(name).get_potential_energy()
    e4 = dct[id]['energy']
    assert e1 == e2 == e3 == e4
    print(e1)
ae = 2 * c.get('name=O').energy - c.get('name=O2').energy
assert abs(ae - 5.664) < 1e-3
