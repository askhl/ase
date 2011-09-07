import os
import sys

import urllib2

import __builtin__

import numpy as np

from ase.test import NotAvailable

from ase.units import kcal, mol

from ase.utils.compound_test import get_atomization_reactions, \
     calculate_reaction_energies, write_csv

from ase.utils.compound_test import import_module

from ase.utils.compound_test import BatchTest

from ase.utils.compound_test.molecule_emt import EMTEnergyMoleculeTest

ea_ref = {
    'G2RC':
    {
    2: 0.39041398166621666,
    3: 0.25281952969937116,
    5: -1.4624731826544464,
    7: -0.40658621500291314,
    12: 0.57521101260687146,
    13: -0.46193828701794626,
    15: 3.0307018972588677,
    17: -0.41693897111701528,
    19: 0.57589628171724749,
    21: 1.3489078861482464,
    22: -0.8345735215176564,
    24: -0.82939132072898492,
    },
    'WATER27':
    {
    1: -0.0030247833463983298,
    2: 0.066159317756866898,
    3: 0.10477620488969031,
    4: 0.10958991182342537,
    5: 0.2416448305136889,
    6: 0.22225614213389377,
    7: 0.17066703383794746,
    8: 0.10063651602143509,
    9: 0.36201395757131039,
    10: 0.3592482723961048,
    11: 0.81778460324132851,
    12: 0.98519685961929326,
    13: 0.90707840977527354,
    14: 0.90558656356087397,
    15: 0.14832793640099728,
    16: 0.064897074169329549,
    17: 0.018031966705603253,
    18: 0.23173121989958645,
    19: 0.17162545836515974,
    20: 0.35931661498829115,
    21: 0.38171648158166338,
    22: 0.42596482435230953,
    23: 0.38035968187663904,
    24: 0.47895171770785439,
    25: 0.50302038543340899,
    26: 0.59109462022364667,
    27: -0.57171144822327946,
    },
    }

project = 'GMTKN30'

# download and create the project databases
try:
    resp = urllib2.urlopen("http://toc.uni-muenster.de/GMTKN/GMTKN30/GMTKN30main.html")
    from ase.data.GMTKN30 import main
    main()
except urllib2.HTTPError:
    raise NotAvailable('Retrieval of GMTKN30 failed')

# database will be read from the current directory
cwd =  os.getcwd()
sys.path.insert(0, cwd)

for database in ['G2RC', 'WATER27']:
    modulepath = project + '_' + database

    module = import_module(modulepath)

    name = database + '/energy'
    # etest object provides all data and methods needed
    # calculate energies on fixed geometries
    etest = EMTEnergyMoleculeTest(name, vacuum=3.0, data=module.data)
    test = BatchTest(etest)

    formulas = etest.get_formulas()

    atoms = etest.get_formulas(natoms=1)

    supported_elements = 'Ni, C, Pt, Ag, H, Al, O, N, Au, Pd, Cu'.split(', ')

    formulas = [formula for formula in formulas
                if np.all([symbol in supported_elements
                           for symbol
                           in etest.compound(formula).get_chemical_symbols()])]

    atoms = [symbol for symbol in atoms if symbol in supported_elements]
    dimers = [formula for formula in formulas if len(etest.compound(formula)) == 2]


    print 'Energy test'
    print '-----------'
    test.run(formulas)

    # use defined reactions
    reactions = module.info['reactions']

    results = etest.get_results(formulas, dir=database)
    write_csv(results, dir=database, outfilename='energy_' + database+ '.csv')
    energies = {}
    for r in results:
        # compound: energy
        energies[r[0]] = r[1]
    ea = calculate_reaction_energies(energies, reactions)
    write_csv(ea, dir=database, outfilename='ea_' + database+ '.csv')

    # print results comparing to the reference
    keys = ea_ref[database].keys()
    keys.sort()
    for k in keys:
        print k, ea[k]*1/(kcal/mol), \
              module.info['reaction energy']['reference'][k]

    # test reaction energies
    ref = ea_ref[database]
    assert len(ea.keys()) == len(ref.keys())
    keys = ref.keys()
    keys.sort()
    for k in keys:
        diff = abs(ea[k] - ref[k])
        assert diff < 1e-5, k + ': ' + str(diff)
