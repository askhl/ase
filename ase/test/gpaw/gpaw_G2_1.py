import os

from ase.test import NotAvailable

try:
    import gpaw
except (ImportError, NotAvailable):
    raise NotAvailable('GPAW required')

from ase.units import kcal, mol

from ase.utils.compound_test import get_atomization_reactions, \
     calculate_reaction_energies, write_csv

#from ase.utils.compound_test import import_module

from ase.utils.compound_test import BatchTest

from ase.utils.compound_test import EnergyMoleculeTest as Test

from gpaw import GPAW

# This is how one may setup a special calculator

ref = {
    'O2':
    {
    'energy': -7.8414743877758557,
    'direct_bandgap_0': 9.7337968681754159,
    'direct_bandgap_1': 6.5925575440115285,
    },
    'O':
    {
    'energy': -0.82310543098329125,
    'direct_bandgap_0': 23.076514407003771,
    'direct_bandgap_1': 0.89958249370291998,
    }
    }

ea_ref = {
    # far from converged but atomization energy is close to converged
    'O2': 6.1952635258092732,
    }

kwargs = dict(
    xc='PBE',
    h=0.3,
    mode='lcao',
    basis='dzp',
    width=0.0,
    nbands=-3,
    fixmom=True,
    )

def gpaw(system, formula, test):
    #
    # http://cms.mpi.univie.ac.at/vasp-forum/forum_viewtopic.php?4.8980.0#9073
    system.cell[1, 1] += 0.01
    system.cell[2, 2] += 0.02
    system.center()
    
    return GPAW(txt=test.get_filename(formula, extension='txt'),
                hund=(len(system) == 1),
                charge=sum(system.get_charges()),
                **kwargs)

database = 'G2_1'

#modulepath = 'ase.data.' + database
#module = import_module(modulepath)

dir = database + '_GPAWb'

identifier = 'energy'

etest = Test(gpaw, dir + '/' + identifier, vacuum=2.5)

betest = BatchTest(etest)

#dimers = etest.get_formulas(natoms=2)

#molecule = 'O2'
#assert molecule in dimers

#formulas = [a.symbol for a in etest.compound(molecule)] # constituing atoms
#formulas = list(set(formulas)) # unique
#formulas.append(molecule)
formulas = ['O','O2']
betest.run(formulas, fraction='1/1')

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

from ase.data.G2_1_vasp_ref import info

for k in ea.keys():
    diff = ea[k] - ea_ref[k] # eV
    assert abs(diff) < 0.01, k + ': ' + str(diff)
    print ea[k], info['atomization energy']['PBE'][k]*(kcal/mol), '(WIEN97)'
