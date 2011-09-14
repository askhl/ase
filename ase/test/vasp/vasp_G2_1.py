import os

from ase.test import NotAvailable

try:
    from ase.test.vasp_installed import vasp_installed
    vasp_installed()
except NotAvailable:
    raise NotAvailable('Vasp required')

from ase.utils.compound_test import get_atomization_reactions, \
     calculate_reaction_energies, write_csv

from ase.utils.compound_test import import_module

from ase.utils.compound_test import BatchTest

from ase.utils.compound_test.molecule_vasp import \
     VASPEnergyMoleculeTest as Test

ref = {
    'O2':
    {
    'energy': -9.4383809999999997,
    'direct_bandgap_0': 5.5442779999999994,
    'direct_bandgap_1': 6.6358930000000012,
    },
    'O':
    {
    'energy': -1.588541,
    'direct_bandgap_0': 7.1638999999999999,
    'direct_bandgap_1': 1.0761020000000006,
    },
    }

ea_ref = {
    'O2': 6.2612989999999993,
    }

kwargs = dict(
    vacuum=2.5,
    xc='PBE',
    algo='Fast', # may give wrong energies
    ediff=1.0e-3,
    prec='low',
    nbands=-3,
    ismear=-1,
    sigma=0.0,
    addgrid=True, # http://cms.mpi.univie.ac.at/vasp-forum/forum_viewtopic.php?4.383
    lasph=True, # DOI: 10.1063/1.1926272: Section IV
    # band parallelization:
    # http://cms.mpi.univie.ac.at/vasp/guide/node138.html
    # without npar=1 you get wrong number of bands in parallel!
    lplane=False,
    npar=1,
    )

database = 'G2_1'

modulepath = 'ase.data.' + database
module = import_module(modulepath)

dir = database + '_VASP'

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
