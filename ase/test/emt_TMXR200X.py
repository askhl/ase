import os
import sys

import __builtin__

import numpy as np

from ase.utils.compound_test import get_atomization_reactions, \
     calculate_reaction_energies, write_csv

from ase.utils.compound_test import import_module

from ase.utils.compound_test import BatchTest

from ase.utils.compound_test.molecule_emt import EMTEnergyMoleculeTest

ea_ref = {
    'TM1R2006':
    {
    'CuMe': 13.212914659644694,
    'Ni(CO)4': 33.109694955280631,
    'Cu(acac)2': 91.402081882891125,
    'Ni(acac)2': 91.673114346973733,
    'CuCN': 9.4891898430752821
    },
    }

project = 'TMXR200X'

for database in ['TM1R2006']:

    modulepath = 'ase.data.' + project + '_' + database

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

    # assume atomization energies
    reactions = get_atomization_reactions(data=etest.data)
    results = etest.get_results(formulas, dir=database)
    write_csv(results, dir=database, outfilename='energy_' + database + '.csv')
    energies = {}
    for r in results:
        # compound: energy
        energies[r[0]] = r[1]
    ea = calculate_reaction_energies(energies, reactions)
    write_csv(ea, dir=database, outfilename='ea_' + database + '.csv')

    # test atomization energies
    ref = ea_ref[database]
    assert len(ea.keys()) == len(ref.keys())
    keys = ref.keys()
    keys.sort()
    for k in keys:
        diff = ea[k] - ref[k]
        assert abs(diff) < 1e-5, k + ': ' + str(diff)
