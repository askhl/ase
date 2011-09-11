import os

from ase.test import NotAvailable

try:
    from ase.test.jacapo_installed import jacapo_installed
    jacapo_installed()
except NotAvailable:
    raise NotAvailable('Jacapo required')

from ase.utils.compound_test import get_atomization_reactions, \
     calculate_reaction_energies, write_csv

from ase.utils.compound_test import import_module

from ase.utils.compound_test import BatchTest

from ase.utils.compound_test.molecule_jacapo import \
     JACAPOEnergyMoleculeTest as Test

ref = {
    'O2':
    {
    'energy': -868.87662296493772,
    'direct_bandgap_0': 5.7791952308747101,
    'direct_bandgap_1': 6.5933598183231066,
    },
    'O':
    {
    'energy': -431.62602720861776,
    'direct_bandgap_0': 7.000288244962876,
    'direct_bandgap_1': 1.5365786012024945,
    },
    }

ea_ref = {
    # this shows how bad is O in dacapo (fully converged result is very close!)
    # but the "right" answer is ~6.2 eV (DOI: 10.1063/1.1926272)
    'O2': 5.6245685477068719,
    }

kwargs = dict(
    vacuum=2.2,
    xc='PBE',
    ft=0.0,
    pw=300,
    dipole=False,
    symmetry=False,
    )

database = 'G2_1'

modulepath = 'ase.data.' + database
module = import_module(modulepath)

dir = database + '_JACAPO'

identifier = 'energy'

etest = Test(dir + '/' + identifier, data=module.data, **kwargs)

betest = BatchTest(etest)

dimers = etest.get_formulas(natoms=2)

molecules = ['O2']
formulas = []
for m in molecules:
    assert m in dimers
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
