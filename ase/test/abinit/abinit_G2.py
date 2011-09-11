# test (if cmr present) extracting of eigenvalues/occupations
# for spin paried/polarized systems with large number of electrons
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
    'C3H7': 28.697913584416938,
    'C3H6_Cs': 27.653685158479959,
    }

ref = {
    'H':
    {
    'energy': -13.021261138720101,
    'direct_bandgap_0': 8.8268325227367033,
    'direct_bandgap_1': -11.999681256170023,
    },
    'C':
    {
    'energy': -141.6508069244357,
    'direct_bandgap_0': 0.36844229717570443,
    'direct_bandgap_1': 9.7209989839770028,
    },
    'C3H6_Cs':
    {
    'energy': -530.73367276410761,
    'direct_bandgap_0': 5.5592881324221857,
    },
    'C3H7':
    {
    'energy': -544.79916232876474,
    'direct_bandgap_0': 3.8237453175132918,
    'direct_bandgap_1': 5.8825595128097312,
    },
    }

kwargs = dict(
    vacuum=2.5,
    xc='PBE',
    width=0.001, # eV abinit does not accept 0.0
    ecut=150,
    nbands=12,
    nstep=120,
    )

database = 'G2'

modulepath = 'ase.data.' + database
module = import_module(modulepath)

dir = database + '_ABINIT'

identifier = 'energy'

etest = Test(dir + '/' + identifier, data=module.data, **kwargs)

betest = BatchTest(etest)

molecules = ['C3H6_Cs', 'C3H7',]
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
