from ase.test import NotAvailable

try:
    from ase.test.fleur_installed import fleur_installed
    fleur_installed()
except NotAvailable:
    raise NotAvailable('FLEUR required')

from ase.units import Hartree

from ase.utils.compound_test import get_atomization_reactions, \
     calculate_reaction_energies, write_csv

from ase.utils.compound_test import import_module

from ase.utils.compound_test import BatchTest

from ase.utils.compound_test.molecule_fleur import \
     FLEUREnergyMoleculeTest as Test


ea_ref = {
    'O2': None,
    }

kwargs = dict(
    vacuum=1.5,
    xc='PBE',
    # fleur is unable to converge most atoms and spin-polarized molecules
    # Crazy settings below will prevent a possible crash (with famous 'inwint')
    # and stop fleur in the middle of SCF - and result is a nonsense!
    kmax=2.0,
    width=0.0001,
    mixer={'imix' : 0, 'alpha' : 0.01, 'spinf' : 2.00},
    convergence={'energy': 2.0 * Hartree}, # eV
    )

database = 'G2_1'

modulepath = 'ase.data.' + database
module = import_module(modulepath)

dir = database + '_FLEUR'

etest = Test(dir + '/energy', data=module.data, **kwargs)

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

ea = calculate_reaction_energies(energies, reactions)
write_csv(ea, dir=dir, outfilename='ea.csv')

if 0: # don't test nonsense result
    for k in ea.keys():
       diff = ea[k] - ea_ref[k] # eV
       assert abs(diff) < 0.01, k + ': ' + str(diff)
