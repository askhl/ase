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
    # abinit converges slowly with vacuum and ecut: this is very far from converged!
    # abinit gives in fact quite good converged result
    'O2': 5.0755627357261801,
    }

ref = {
    'O2':
    {
    'energy': -828.26649603829867,
    'direct_bandgap_0': 6.2050145513276123,
    'direct_bandgap_1': 7.4472147630019769,
    },
    'O':
    {
    'energy': -411.59545654501886,
    'direct_bandgap_0': 10.469040250547174,
    'direct_bandgap_1': 1.0013793601230372,
    },
    }

kwargs = dict(
    vacuum=2.5,
    xc='PBE',
    width=0.001, # eV abinit does not accept 0.0
    ecut=300,
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
