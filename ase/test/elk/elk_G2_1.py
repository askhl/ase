from ase.test import NotAvailable

try:
    from ase.test.elk_installed import elk_installed
    elk_installed()
except NotAvailable:
    raise NotAvailable('ELK required')

from ase.units import Hartree

from ase.utils.compound_test import get_atomization_reactions, \
     calculate_reaction_energies, write_csv

from ase.utils.compound_test import import_module

from ase.utils.compound_test import BatchTest

from ase.utils.compound_test.molecule_elk import \
     ELKEnergyMoleculeTest as Test


ea_ref = {
    'O2': None,
    }

kwargs = dict(
    vacuum=2.5,
    xc='PBE',
    autormt=True,
    autokpt=False,
    swidth=0.001, # eV elk does not accept 0.0
    stype=3, # Fermi
    # elk 1.3.31 is unable to converge (epspot) most molecules with PBE
    # and gives wrong magentic moments for atoms
    # see http://sourceforge.net/projects/elk/forums/forum/897820/topic/4058497
    # Crazy settings below will stop elk in the middle of SCF - result is a nonsense!
    rgkmax=3.5,
    gmaxvr=6,
    epspot=5.0e-1,
    epsengy=5.0e-1, # Hartree
    )

database = 'G2_1'

modulepath = 'ase.data.' + database
module = import_module(modulepath)

dir = database + '_ELK'

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
