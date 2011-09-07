"""This module defines extensible classes for running tests on molecules.

Use this to compare different calculators, XC functionals and so on by
calculating e.g. atomization energies, bond lengths across databases
of molecules.
"""

import os
import sys

import __builtin__

import numpy as np

from ase.parallel import rank

from ase.calculators.emt import EMT

from ase.data.molecules import latex

from ase.utils.compound_test import get_atomization_reactions, \
     calculate_reaction_energies, write_csv

from ase.utils.compound_test import import_module

from ase.utils.compound_test import BatchTest
from ase.utils.compound_test import MoleculeTest
from ase.utils.compound_test import EnergyMoleculeTest, \
     GeometryMoleculeTest, BondLengthDimerTest, \
     BondLengthFitDimerTest

from ase.utils.compound_test import BondLengthTest

class EMTMoleculeTest(MoleculeTest):

    def setup_calculator(self, system, calculator):
        return EMT()

class EMTEnergyMoleculeCMRTest(EnergyMoleculeTest, EMTMoleculeTest):
    """This uses init from EMTMoleculeTest and run from EnergyMoleculeTest.  """

    def check_system(self, system, formula):
        pass

    def setup_calculator(self, system, calculator):
        return EMTMoleculeTest.setup_calculator(self, system, calculator)

    def write_db(self, filename, system, **kwargs):
        write_db(filename, system, **kwargs)

class EMTEnergyMoleculeTest(EnergyMoleculeTest, EMTMoleculeTest):
    """This uses init from EMTMoleculeTest and run from EnergyMoleculeTest.  """

    def check_system(self, system, formula):
        pass

    def setup_calculator(self, system, calculator):
        return EMTMoleculeTest.setup_calculator(self, system, calculator)

    def write_db(self, filename, system, **kwargs):
        pass


class EMTBondLengthFitDimerTest(BondLengthFitDimerTest, EMTMoleculeTest):

    def get_linspace_parameters(self):
        """Owervrite the method to use different sampling.  """
        return 0.96, 1.04, 5

    def check_system(self, system, formula):
        pass

    def setup_calculator(self, system, calculator):
        return EMTMoleculeTest.setup_calculator(self, system, calculator)

    def write_db(self, filename, system, **kwargs):
        pass


class EMTBondLengthDimerTest(BondLengthDimerTest, EMTMoleculeTest):

    def check_system(self, system, formula):
        pass

    def setup_calculator(self, system, calculator):
        return EMTMoleculeTest.setup_calculator(self, system, calculator)

    def write_db(self, filename, system, **kwargs):
        pass


def plot_opt(systems, opt_path, opt_values, references, mode='fit',
             xlabel=u'Bond length [A]', ylabel='Energy [eV]',
             figname='bond.png', dir='.'):
    """Plot results of optimization around reference values.  """
    import numpy as np
    import os
    try:
        import matplotlib
        matplotlib.use('Agg')
        import pylab as plt
    except ImportError:
        pass
    #
    assert mode in ['fit', 'opt']
    B = []
    E0 = []
    opt_path_keys = opt_path.keys()
    opt_path_keys.sort()
    for formula in opt_path_keys:
        bonds = [a[0] for a in opt_path[formula]]
        energies = [a[1] for a in opt_path[formula]]
        if not formula in references.keys():
            continue
        earef = references[formula]['ea']
        bref = references[formula]['d0']
        ea = opt_values[formula]
        if mode == 'opt':
            e0 = energies[-1]
            b0 = bonds[-1]
            # need at least 2 optimizer steps to fit a fake parabola
            if len(energies) < 2:
                if formula == opt_path_keys[0]:
                    plt.plot([b0], [e0-ea], 'g.', color='0.7', label=mode)
                else:
                    plt.plot([b0], [e0-ea], 'g.', color='0.7', label='_nolegend_')
            else:
                # a fake parabola centered at b0, using e0 and max(energies) points
                maxe = max([e for e in energies])
                efit = [maxe, e0, maxe]
                maxb = bonds[energies.index(maxe)]
                if maxb > b0:
                    bfit = [b0 - (maxb - b0), b0, maxb]
                else:
                    bfit = [maxb, b0, b0 + (b0 - maxb)]
                e = np.polyfit(bfit, efit, 2)
                b = np.linspace(bfit[0], bfit[-1], 20)
        else: # manual fit
            # conssistency check
            from ase.utils.compound_test.molecule_emt import EMTBondLengthFitDimerTest
            dummy = EMTBondLengthFitDimerTest()
            start, stop, num = dummy.get_linspace_parameters()
            del dummy
            btest = np.linspace(start * bref, stop * bref, num)
            for n in range(num):
                assert abs(btest[n] - bonds[n]) < 1e-10, 'Error: inconsistent bond length sampling'
            e0 =  energies[len(energies)/2] # use ~middle value for labels
            assert len(energies) == len(bonds), 'data missing for: '+formula
            e = np.polyfit(bonds, energies, 3)
            dedb = np.polyder(e, 1)
            b0 = np.roots(dedb)[1]
            b = np.linspace(bonds[0], bonds[-1], 20)
        #
        e = np.polyval(e, b) - ea
        if formula == opt_path_keys[0]:
            plt.plot(b, e, '-', color='0.7', label=mode)
        else:
            plt.plot(b, e, '-', color='0.7', label='_nolegend_')
        name = latex(systems[formula]['name'])
        plt.text(b0, e0 - ea + 0.0, name, color='0.7')
        plt.text(bref, - earef + 0.0, name)
        B.append(bref)
        E0.append(-earef)

    plt.plot(B, E0, 'g.', label='reference')
    plt.legend(loc='lower right')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(os.path.join(dir, figname))
    plt.clf()

def write_db(filename, system, **kwargs):
    """Write a minimal db file.  """
    import os
    from ase.dft.kpoints import monkhorst_pack, get_monkhorst_shape
    from ase.calculators.eigenvalues import get_bandgap
    try:
        import cmr
        from cmr.base.converter import Converter
        from cmr.static import CALCULATOR_SPREADSHEET
        import datetime
        os.environ["CMR_REPOSITORY"]= "."
        symbols = system.get_chemical_symbols()
        calc = system.get_calculator()
        params = {
            "user": os.environ["USER"],
            "location": "CAMd/DTU/Physics",
            "date": datetime.datetime.today(),
            "description": filename.split('.')[0],
            # atoms
            "charges": system.get_charges(),
            "magmoms": system.get_magnetic_moments(),
            "name": ''.join(symbols),
            "positions": system.get_positions(),
            "symbols": symbols,
            "cell": system.get_cell(),
            "volume": system.get_volume(),
            # calculator
            "energy": system.get_potential_energy(), # eV
            }
        # kwargs: time, number of iterations, bulk modulus, etc.
        for k, v in kwargs.items():
            params.update({k: v})
        spreadsheet = True
        #add output parameters:
        params["output"] = filename
        params["db"] = False #add to CMR_REPOSITORY
        params["private"] = 600
        if spreadsheet:
            data = Converter.get_xml_writer(CALCULATOR_SPREADSHEET)
            #
            for p in params.keys():
                data.set_user_variable(p, params[p])
            #
            if rank == 0:
                data.write(params)
    except ImportError: # no cmr
        pass


def main():

    import os
    import sys

    import __builtin__

    from ase.units import kcal, mol
    from ase.io import read

    database = 'G2_1'

    modulepath = 'ase.data.' + database
    module = import_module(modulepath)

    get_atomization_energy = module.get_atomization_energy

    dir = database + '_EMT'
    name1 = dir + '/energy'
    name2 = dir + '/bond'
    name3 = dir + '/bondfit'
    # etest object provides all data and methods needed
    # calculate energies on fixed geometries
    try:
        import cmr
        etest = EMTEnergyMoleculeCMRTest(name1, vacuum=3.0, data=module.data)
    except ImportError:
        etest = EMTEnergyMoleculeTest(name1, vacuum=3.0, data=module.data)
    test1 = BatchTest(etest)
    # calculate energies on relaxed geometries
    btest1 = EMTBondLengthDimerTest(name2, vacuum=3.0, data=module.data)
    test2 = BatchTest(btest1)
    # calculate energies on manually fitted geometries
    btest2 = EMTBondLengthFitDimerTest(name3, vacuum=3.0, data=module.data)
    test3 = BatchTest(btest2)

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
    test1.run(formulas)

    # assume atomization energies
    reactions = get_atomization_reactions(data=etest.data)
    results = etest.get_results(formulas, dir=dir)
    write_csv(results, dir=dir, outfilename='energy.csv')
    energies = {}
    for r in results:
        # compound: energy
        energies[r[0]] = r[1]
    ea = calculate_reaction_energies(energies, reactions)
    write_csv(ea, dir=dir, outfilename='ea.csv')

    print
    print 'Bond length test'
    print '----------------'
    test2.run(dimers + atoms)

    # assume atomization energies
    reactions = get_atomization_reactions(data=etest.data)
    results = etest.get_results(dimers + atoms, dir=dir, identifier='bond')
    write_csv(results, dir=dir, outfilename='bond.csv')
    energies = {}
    for r in results:
        # compound: energy
        energies[r[0]] = r[1]
    ea_opt = calculate_reaction_energies(energies, reactions)
    write_csv(ea_opt, dir=dir, outfilename='ea_opt.csv')
    # optimized bond length
    results = btest1.get_results(dimers, dir=dir)
    write_csv(results, dir=dir, outfilename='energy_opt.csv')
    bond_opt = {}
    for r in results:
        # compound: bond
        bond_opt[r[0]] = r[1]
    for formula, bond in test2.collect(dimers, verbose=False):
        pass # use it for test coverage
    write_csv(bond_opt, dir=dir, outfilename='bond_opt.csv')
    # the remaining part is needed to make a plot
    # store formula and (bond, energy) values
    opt_path = {}
    for formula in dimers:
        filename = btest1.get_filename(formula, extension='traj')
        nsteps = len(read(filename))
        data = []
        for n in range(nsteps):
            e = btest1.retrieve_energy(formula, filename, index=n)
            d = btest1.retrieve_distance(formula, filename, mic=False, index=n)
            data.append([d, e])
        opt_path[formula] = data
    # plot bond lengths
    references = {}
    for f in dimers:
        references[f] = {
        # experimental enthalpy of formation
        'ea': get_atomization_energy(f) * kcal/mol,
        # reference bond length
        'd0': btest1.compound(f).get_distance(0, 1, mic=False),
        }
    # plot only reasonable EMT results
    opt_path_plot = {}
    for f in ['CN', 'N2', 'NO', 'O2']:
        opt_path_plot[f] = opt_path[f]
    plot_opt(btest1.data, opt_path_plot, ea_opt, references, mode='opt',
             figname='bond_opt.png', dir=dir)

    print
    print 'Bond length test fit'
    print '--------------------'
    test3.run(dimers + atoms)

    # assume atomization energies
    reactions = get_atomization_reactions(data=etest.data)
    results = etest.get_results(dimers + atoms, dir=dir, identifier='bondfit')
    write_csv(results, dir=dir, outfilename='bondfit.csv')
    # we don't have atomization energies - we fit them
    # optimized energy
    energies = {}
    for r in results:
        # compound: energy
        energies[r[0]] = r[1]
    # optimized bond length
    bond_fit = {}
    # store formula and (bond, energy) values
    opt_path = {}
    for formula in dimers:
        filename = btest2.get_filename(formula, extension='traj')
        d, e, d0, e0, polynomial = btest2.retrieve(formula, filename)
        # we have fitted energy and bond now
        energies[formula] = e0 # overwrite dimers energies
        bond_fit[formula] = d0
        nsteps = len(d)
        data = []
        for n in range(nsteps):
            data.append([d[n], e[n]])
        opt_path[formula] = data
    #
    ea_fit = calculate_reaction_energies(energies, reactions)
    write_csv(ea_fit, dir=dir, outfilename='ea_fit.csv')
    write_csv(bond_fit, dir=dir, outfilename='bond_fit.csv')
    # the remaining part is needed to make a plot
    # plot bond lengths
    references = {}
    for f in dimers:
        references[f] = {
        # experimental enthalpy of formation
        'ea': get_atomization_energy(f) * kcal/mol,
        # reference bond length
        'd0': btest2.compound(f).get_distance(0, 1, mic=False),
        }
    # plot only reasonable EMT results
    opt_path_plot = {}
    for f in ['CN', 'N2', 'NO', 'O2']:
        opt_path_plot[f] = opt_path[f]
    plot_opt(btest2.data, opt_path_plot, ea_fit, references, mode='fit',
             figname='bond_fit.png', dir=dir)
    print
    print 'Atomization energies: EMT vs reference [eV]'
    print '-------------------------------------------'
    energies = {}
    for f in formulas:
        system = etest.compound(f)
        if len(system) > 1: # molecule, ea (ref), ea (fixed)
            try:
                ref = references[f]['ea']
                values = (f, ea[f], ref)
                print '%10s %5.1f %5.1f' % values
            except KeyError:
                ref = None # no reference bond length
                values = (f, ea[f], ref)
                print '%10s %5.1f %s' % values
            energies[f] = values[1:]

    print
    print 'Bond lengths: EMT optimized and fitted vs reference [A]'
    print '-------------------------------------------------------'
    bonds = {}
    for f in dimers:
        # molecule, bond (ref), bond (opt), bond (fit)
        values = (f, bond_opt[f], bond_fit[f], references[f]['d0'])
        print '%10s %6.3f %6.3f %6.3f' % values
        bonds[f] = values[1:]
    print 'Note: the fit method is wrong due to fitting outside of range!'

    return energies, bonds

if __name__ == '__main__':
    energies, bonds = main()
