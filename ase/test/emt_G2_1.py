import sys

import numpy as np

from ase.utils.compound_test.molecule_emt import main

# http://stackoverflow.com/questions/2158395/flatten-an-irregular-list-of-lists-in-python
def flatten(x):
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

def array_almost_equal(a1, a2, tol=np.finfo(type(1.0)).eps):
    return (np.abs(a1 - a2) < tol).all()

energies_ref = {
    'C2H2': (11.686151320322335, None),
    'C2H4': (17.415147328268556, None),
    'C2H6': (23.181060638046283, None),
    'CH': (4.5774752309526692, 3.6378069042612404),
    'CH2_s1A1d': (8.0249469958387678, None),
    'CH2_s3B1d': (7.9725918072023827, None),
    'CH3': (11.156535835199753, None),
    'CH3OH': (17.976338739425081, None),
    'CH4': (14.349422888656331, None),
    'CN': (7.8355973663682041, 7.8667794377531104),
    'CO': (7.3137926543566278, 11.242584697775644),
    'CO2': (11.813727761098084, None),
    'H2CO': (13.120280586356913, None),
    'H2O': (8.4001892177132405, None),
    'H2O2': (12.517207489925241, None),
    'HCN': (10.55339402193094, None),
    'HCO': (10.243137097731383, None),
    'N2': (9.6512403798077138, 9.9077288388938562),
    'N2H4': (18.364114597251039, None),
    'NH': (5.4242942047339024, 3.6234100182010183),
    'NH2': (8.6424786489844685, None),
    'NH3': (11.363494365881161, None),
    'NO': (8.9975876653724285, 6.6222163383731045),
    'O2': (8.2773240044713159, 5.217552927745098),
    'OH': (5.5639474941954337, 4.6160146022927657),
    }

bonds_ref = {
    'CH': (1.0706867264687667, 1.0704290729529695, 1.12052),
    'CN': (0.98337140336546236, 0.97798359582285044, 1.1347990000000001),
    'CO': (1.0563275075110994, 1.0547109735465099, 1.1503399999999999),
    'N2': (0.99808780893191384, 0.994513209874611, 1.12998),
    'NH': (1.0578946697335252, 1.0576946966818899, 1.0394300000000001),
    'NO': (1.059474093772157, 1.0585972546765832, 1.142703),
    'O2': (1.1003327119060096, 1.0955369395158239, 1.2459560000000001),
    'OH': (1.0321150670542552, 1.0328365109000035, 0.97906999999999988),
    }

# The G2 module may be already loaded (by using ase.data.molecule method),
# so it has to be deleted (together with it's submodules)
# in order to access G2_1 module only!
# See ase.utils.compound_test
modulepathroot = 'ase.data.'
for module in ['G2', 'G2_2', 'G2_1']:
    modulepath = modulepathroot + module
    if modulepath in sys.modules:
        del(sys.modules[modulepath])

energies, bonds = main()

assert len(energies.keys()) == len(energies_ref.keys()), str(energies.keys())+ str(energies_ref.keys())
assert len(bonds.keys()) == len(bonds_ref.keys())
# test energies
keys = energies_ref.keys()
keys.sort()
# skip None and complex numbers
ref = [[e for e in energies_ref[f] if ((not e is None) and np.isreal(e))] for f in keys]
# skip None and complex numbers
test = [[e for e in energies[f] if ((not e is None) and np.isreal(e))] for f in keys]
assert array_almost_equal(np.array(flatten(test)), np.array(flatten(ref)), tol=1e-5)

# test bonds
keys = bonds_ref.keys()
keys.sort()
# skip None and complex numbers
ref = [[e for e in bonds_ref[f] if ((not e is None) and np.isreal(e))] for f in keys]
# skip None and complex numbers
test = [[e for e in bonds[f] if ((not e is None) and np.isreal(e))] for f in keys]
assert array_almost_equal(np.array(flatten(test)), np.array(flatten(ref)), tol=1e-3)
