# test (if cmr present) extracting of eigenvalues/occupations
# for spin paried/polarized systems with small number of electrons
# Warning: abinit converges slowly with vacuum and ecut:
# parameters used in this run are nonsense!

import os

from ase.test import NotAvailable

try:
    from ase.test.abinit_installed import abinit_installed
    abinit_installed()
except NotAvailable:
    raise NotAvailable('Abinit required')

from ase.units import Hartree

from ase.utils.compound_test import get_atomization_reactions, \
     calculate_reaction_energies, write_csv

from ase.utils.compound_test import import_module

from ase.utils.compound_test import BatchTest

from ase.utils.compound_test.molecule_abinit import \
     ABINITEnergyMoleculeTest as Test


ea_ref = {
    'N2': 5.2980191339718203,
    'O2': 1.2575288252661494,
    }

ref = {
    'N2':
    {
    'energy': -501.8774336939382,
    'direct_bandgap_0': 8.3275034124579612,
    },
    'O2':
    {
    'energy': -779.53098872176645,
    'direct_bandgap_0': 4.1559964584671585,
    'direct_bandgap_1': 8.1642350385248577,
    },
    'O':
    {
    'energy': -389.13672994825015,
    'direct_bandgap_0': 9.4380004691596238,
    'direct_bandgap_1': 0.82668220001461634,
    },
    'N':
    {
    'energy': -248.28970727998319,
    'direct_bandgap_0': 9.0390814088497393,
    'direct_bandgap_1': 17.997617086559146,
    },
    }

kwargs = dict(
    vacuum=2.5,
    xc='PBE',
    width=0.001, # eV abinit does not accept 0.0
    ecut=150,
    nbands=8,
    nstep=120,
    )

database = 'G2_1'

modulepath = 'ase.data.' + database
module = import_module(modulepath)

dir = database + '_ABINIT'

identifier = 'energy'

etest = Test(dir + '/' + identifier, data=module.data, **kwargs)

betest = BatchTest(etest)

molecules = ['N2', 'O2']
formulas = []
for m in molecules:
    assert m in etest.get_formulas()
    formulas.extend([a.symbol for a in etest.compound(m)]) # constituing atoms
formulas = list(set(formulas)) # unique
formulas.extend(molecules)

betest.run(formulas)

# analyse results

reactions = get_atomization_reactions(data=etest.data)
results = etest.get_results(formulas, dir=dir)
energies = {}
for r in results:
    # compound: energy
    energies[r[0]] = r[1]

db_data = {}
for formula in formulas:
    try:
        db_data[formula] = {}
        import cmr
        filename = os.path.join(dir, identifier + '.' + formula)
        db = cmr.read(filename + '.db')
        for item in ['iter', 'energy',
                     'direct_bandgap_0',
                     'direct_bandgap_1',
                     ]:
            db_data[formula][item] = db[item]
        # test
        for k in ref[formula].keys():
            diff = db_data[formula][k] - ref[formula][k]
            assert abs(diff) < 1.0e-2, k + ': ' + str(diff)
    except (ImportError, IOError, KeyError): # no cmr or no db data
        pass

ea = calculate_reaction_energies(energies, reactions)
write_csv(ea, dir=dir, outfilename='ea.csv')

for k in ea.keys():
   diff = ea[k] - ea_ref[k] # eV
   assert abs(diff) < 0.01, k + ': ' + str(diff)
