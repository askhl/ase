"""This module defines extensible classes for running batch tests.

Use this to compare different calculators, XC functionals and so on by
calculating e.g. atomization energies,
bond lengths across databases of molecules, bulk systems.

One compound test to rule them all
One compound test to run them
One compound test to save them all
And on the webpage plot them (implementation pending)
"""

import os
import sys
import traceback
import time

import __builtin__

import numpy as np

from ase import Atoms
from ase.atoms import string2symbols
from ase.io import read
from ase.io.trajectory import PickleTrajectory
from ase.optimize import QuasiNewton
from ase.optimize.sciopt import Converged
if 0: # other optimizers
    from ase.optimize.oldqn import GoodOldQuasiNewton as QuasiNewton
    from ase.optimize.sciopt import SciPyFminBFGS as QuasiNewton
    from ase.optimize.sciopt import SciPyFminCG as QuasiNewton
    from ase.optimize.bfgslinesearch import BFGSLineSearch as QuasiNewton
    from ase.optimize.bfgs import BFGS as QuasiNewton
from ase.parallel import barrier, rank, paropen
from ase.parallel import paropen as open
#from ase.data.molecules import molecule


class BatchTest:
    """Contains logic for looping over tests and file management."""
    def __init__(self, test):
        self.test = test
        self.txt = sys.stdout # ?

    def run_single_test(self, formula):
        filename = self.test.get_filename(formula, extension='traj')
        if rank == 0:
            print >> self.txt, self.test.name, formula, '...',
        barrier()
        self.txt.flush()
        if os.path.exists(filename):
            if rank == 0:
                print >> self.txt, 'Skipped.'
            return
        barrier()
        try:
            open(filename, 'w').close() # Empty file
            system = self.test.setup(formula)
            self.test.run(formula, system, filename)
            if rank == 0:
                print >> self.txt, 'OK!'
            self.txt.flush()
        except 'asdfg':
            if rank == 0:
                print >> self.txt, 'Failed!'
            traceback.print_exc(file=self.txt)
            print >> self.txt
            self.txt.flush()

    def run(self, formulas, fraction='1/1'):
        """Run a batch of tests.

        This will invoke the run_single_test method on each formula, printing
        status to stdout.

        The formulas that already have *.traj files are skipped.

        Keyword fraction allows one to split the run over several sets,
        for example:
        fraction='1/4' will appox. run the first one-fourth of the systems,
        fraction='4/4' will appox. run the last one-fourth of the systems
        (including the remaining systems if not divisible by 4),
        fraction=1/1 runs all the systems.

        """

        if rank == 0:
            # Create directories if necessary
            if self.test.dir and not os.path.isdir(self.test.dir):
                os.mkdir(self.test.dir) # Won't work on 'dir1/dir2', but oh well
        barrier()

        # split the formulas set into subsets
        nominator, denominator = fraction.split('/')
        nominator, denominator = int(nominator), int(denominator)
        assert nominator> 0
        assert nominator<= denominator
        reminder = len(formulas) % denominator
        quotient = int(len(formulas)/denominator)
        start = (nominator-1)*quotient
        if nominator == denominator:
            formulas_set = formulas[start:]
        else:
            stop = nominator*quotient
            formulas_set = formulas[start:stop]

        for formula in formulas_set:
            self.run_single_test(formula)

    def collect(self, formulas, verbose=False):
        """Yield results of previous calculations."""
        for formula in formulas:
            try:
                filename = self.test.get_filename(formula, extension='traj')
                results = self.test.retrieve(formula, filename)
                if verbose:
                    if rank == 0:
                        print >> self.txt, 'Loaded:', formula, filename
                yield formula, results
            except (IOError, RuntimeError, TypeError):
                # XXX which errors should we actually catch?
                if verbose:
                    if rank == 0:
                        print >> self.txt, 'Error:', formula, '[%s]' % filename
                    traceback.print_exc(file=self.txt)

class CompoundTest:
    """Generic class for runnings various tests.

    Usage: instantiate CompoundTest with desired test settings and
    invoke its run() method on the desired formulas.

    You can create a subclass using an arbitrary calculator by overriding the
    setup_calculator method.  Most methods can be overridden to
    provide highly customized behaviour.  """

    def __init__(self, calculator, name='test'):
        """Create a test.

        The name parameter will be part of all output files generated
        by this test.  If name contains a '/' character, the
        preceding part will be interpreted as a directory in which to
        put files.

        The vacuum parameter is used to set the cell size
        (along z-axis for PeriodicSystem), or cell is set directly to cell.

        A tuple of exception types can be provided which will be
        caught during a batch of calculations.  Types not specified
        will be considered fatal.

        data is a dictionary of full definitions of compounds
        (sufficient to create an Atoms object instance).
        If data is None then ase.data.G2 is used.  """

        self.calculator = calculator

        dir, path = os.path.split(name)
        self.dir = dir
        self.identifier = path
        self.name = name


    def get_formulas(self, natoms=None):
        """Get sorted list of systems with given number of atoms in database.  """
        if natoms is None:
            l = [s for s in self.data.keys()]
        else:
            assert natoms > 0, 'Error: number of atoms must be positive: ' + str(natoms)
            l = [s for s in self.data.keys() if (len(self.compound(s)) == natoms)]
        l = list(set(l)) # unique
        l.sort()
        return l

    def setup_calculator(self, system, formula):
        """Create a new calculator.

        Every implementation has to provide this method."""
        return self.calculator(system, formula, self)

    def setup_system(self, formula):
        """Create an Atoms object from the given formula.

        Every implementation has to provide this method."""
        raise NotImplementedError

    def check_system(self, system, formula):
        """Check the electronic state.

        Every implementation has to provide this method."""
        raise NotImplementedError

    def setup(self, formula):
        """Build calculator and atoms objects.

        This will invoke the setup_calculator and setup_system methods."""
        system = self.setup_system(formula)
        calc = self.setup_calculator(system, formula)
        system.set_calculator(calc)
        return system

    def get_filename(self, formula, extension=''):
        """Returns the filename for a test result file.

        Default format is <name>.<formula>[.<extension>]

        The test may write other files, but this filename is used as a
        flag denoting whether the calculation has been done
        already."""
        if extension == '':
            return '.'.join([self.name, formula])
        else:
            return '.'.join([self.name, formula, extension])

    def run(self, formula, system, filename):
        """Calculate energy of specified system and save to file.

        Every implementation has to provide this method."""
        raise NotImplementedError

    def retrieve(self, formula, filename):
        """Retrieve results of previous calculation from file.

        Default implementation returns the total energy.

        This method should be overridden whenever the test method is
        overridden to calculate something else than the total energy."""
        raise NotImplementedError

    def get_results(self, formulas, dir='test', identifier=None):
        """Return the formula and value (given by self.retrieve),
        number of iterations and time per iteration (if found in db file).
        Return results as a list.  """
        if identifier is None:
            identifier = self.identifier
        results = []
        for formula in formulas:
            value = t = iter = 0
            titer = np.nan
            filename = os.path.join(dir, identifier + '.' + formula)
            # read time information from db file
            try:
                import cmr
                data = cmr.read(filename + '.db')
                t = data['time']
                iter = data['iter']
                if iter != 0:
                    titer = t / iter
            except (ImportError, IOError, KeyError): # no cmr or no db data
                pass
            try:
                trajectory = filename + '.traj'
                value = self.retrieve(formula, trajectory)
                # formula, value, number of iterations, time per iteration
                row = [formula, value]
                row.extend([iter, titer])
                results.append(row)
            except IOError:
                continue
        return results

    def write_db(self, filename, system, **kwargs):
        write_db(filename, system, **kwargs)

class MoleculeTest(CompoundTest):
    """Generic class for runnings various molecule tests.

    Usage: instantiate MoleculeTest with desired test settings and
    invoke its run() method on the desired formulas.

    You can create a subclass using an arbitrary calculator by overriding the
    setup_calculator method.  Most methods can be overridden to
    provide highly customized behaviour.  """

    def __init__(self, calculator, name='test',
                 vacuum=8.0, cell=None,
                 data=None):
        CompoundTest.__init__(self, calculator, name)

        self.vacuum = vacuum
        self.cell = cell

        if data is None:
            from ase.data.G2_1 import data
        self.data = data

    def compound(self, name, **kwargs):
        """Create formula from the database."""
        data = self.data
        if name not in data.keys():
            raise NotImplementedError('System %s not in database.'
                                      % (name))
        d = data[name]
        if 'magmoms' not in kwargs:
            kwargs['magmoms'] = d['magmoms']
        if 'charges' not in kwargs:
            kwargs['charges'] = d['charges']
        if 'cell' not in kwargs:
            try:
                kwargs['cell'] = d['cell']
            except KeyError:
                pass # cell not specified
        return Atoms(symbols=d['symbols'],
                     positions=d['positions'],
                     **kwargs)

    def setup_system(self, formula):
        """Create an Atoms object from the given formula.

        By default this will be loaded from the database, setting
        the cell size by means of the molecule test's vacuum parameter."""
        #system = molecule(formula)
        system = self.compound(formula)
        if self.vacuum is not None:
            system.center(vacuum=self.vacuum)
        if self.cell is not None:
            system.set_cell(self.cell)
        return system

    def check_system(self, system, formula):
        """Check the electronic state.  """
        calc = system.get_calculator()
        nelectrons = calc.get_number_of_electrons()
        nspins = calc.get_number_of_spins()
        fa = calc.get_occupation_numbers(spin=0)
        width = calc.get_electronic_temperature()
        if nspins == 1:
            assert (fa.sum() - nelectrons) < 1e-6, formula+ ' :wrong number of electrons: ' + str(fa.sum())
        if not (width > 0.0):
            assert ((fa.round() - fa)**2).sum() < 1e-14, formula + ' :large fractional occupancies: '  + str(fa)
        if nspins == 2:
            fb = calc.get_occupation_numbers(spin=1)
            # program may prefer to occupy spin=0 or spin=1 first
            if len([a for a in fa if abs(a) > 0.0]) > len([a for a in fb if abs(a) > 0.0]):
                ftmp = fa
                fa = fb
                fb = ftmp
            if not (width > 0.0):
                assert ((fb.round() - fb)**2).sum() < 1e-9, formula+ ' :large fractional occupancies: ' + str(fb)
            M = 0.0
            if (not self.data.get(formula, None) is None): # formula in data
                if (not self.data[formula]['magmoms'] is None): # formula has non-zero magmoms
                    M = sum(self.data[formula]['magmoms'])
            assert abs(abs((fa-fb).sum()) - M) < 1e-9, formula + ' :incorrect magnetic moment: ' + str((fa-fb).sum())


class PeriodicSystemTest(CompoundTest):
    """Generic class for runnings various tests.

    Usage: instantiate PeriodicSystemTest with desired test settings and
    invoke its run() method on the desired formulas.

    You can create a subclass using an arbitrary calculator by overriding the
    setup_calculator method.  Most methods can be overridden to
    provide highly customized behaviour.  """

    def __init__(self, name='test',
                 vacuum=0.0, cell=None,
                 exceptions=None, data=None):
        CompoundTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                             exceptions=exceptions, data=data)

    def setup_system(self, formula):
        """Create an Atoms object from the given formula.

        By default this will be loaded from the database, setting
        the cell size by means of the molecule test's vacuum parameter."""
        system = self.compound(formula)
        if not self.vacuum is None:
            cell = system.get_cell()
            shape = cell.shape
            cell = np.ravel(cell)
            cell[-1] += self.vacuum # z-axis
            cell.resize(shape)
            system.set_cell(cell)
            # MDTMP: system.center(vacuum=self.vacuum, axis=2)
            # cannot be used as vacuum is added related to atoms not the cell
        if not self.cell is None:
            system.set_cell(self.cell)
        system.set_pbc([1, 1, 1])
        return system

    def check_system(self, system, formula):
        """Check.  """
        pass

class EnergyMoleculeTest(MoleculeTest):
    """This uses init from MoleculeTest.  """

    def add2db(self, system, iterable):
        kwargs = {
            'iter': iterable.get_number_of_iterations(),
            }
        return kwargs

    def run(self, formula, system, filename):
        """Calculate energy of specified system and save to file."""
        t0 = time.time()
        system.get_potential_energy()
        self.check_system(system, formula)
        kwargs = {
            't': time.time() - t0,
            }
        kwargsup = self.add2db(system, system.get_calculator())
        if kwargsup: kwargs.update(kwargsup)
        self.write_db(os.path.splitext(filename)[0] + '.db', system, **kwargs)
        # Won't create .bak file:
        traj = PickleTrajectory(open(filename, 'w'), 'w')
        traj.write(system)
        traj.close()

    def retrieve_energy(self, formula, filename):
        system = read(filename)
        energy = system.get_potential_energy()
        return energy

    def retrieve(self, formula, filename):
        return self.retrieve_energy(formula, filename)

    def calculate_atomization_energies(self, molecular_energies,
                                       atomic_energies):
        atomic_energy_dict = dict(atomic_energies)
        for formula, molecular_energy in molecular_energies:
            try:
                system = self.compound(formula)
                atomic = [atomic_energy_dict[s]
                          for s in system.get_chemical_symbols()]
                atomization_energy = molecular_energy - sum(atomic)
                yield formula, atomization_energy
            except KeyError:
                pass


class EnergyPeriodicSystemTest(PeriodicSystemTest):
    """This uses init from PeriodicSystemTest.  """

    def add2db(self, system, iterable):
        kwargs = {
            'iter': iterable.get_number_of_iterations(),
            }
        return kwargs

    def run(self, formula, system, filename):
        """Calculate energy of specified system and save to file."""
        t0 = time.time()
        system.get_potential_energy()
        self.check_system(system, formula)
        kwargs = {
            't': time.time() - t0,
            }
        kwargsup = self.add2db(system, system.get_calculator())
        if kwargsup: kwargs.update(kwargsup)
        self.write_db(os.path.splitext(filename)[0] + '.db', system, **kwargs)
        # Won't create .bak file:
        traj = PickleTrajectory(open(filename, 'w'), 'w')
        traj.write(system)
        traj.close()

    def retrieve_energy(self, formula, filename):
        system = read(filename)
        energy = system.get_potential_energy()
        return energy

    def retrieve(self, formula, filename):
        return self.retrieve_energy(formula, filename)

    def calculate_atomization_energies(self, molecular_energies,
                                       atomic_energies):
        raise NotImplementedError


class GeometryMoleculeTest(MoleculeTest):
    """This uses init from MoleculeTest.  """

    def add2db(self, system, iterable):
        kwargs = {
            'iter': iterable.get_number_of_steps() + 1,
            }
        return kwargs

    def run(self, formula, system, filename):
        """Calculate geometry of a molecule using the default optimizer.

        Warning: determination of total energies by varying atomic
        separations close to the bond length, allowing
        determination of bond length by fitting is unreliable when
        sampled around a bond length far from equilibrium!
        """
        dyn = QuasiNewton(system, logfile=filename+'.log' ,trajectory=filename)
        t0 = time.time()
        try:
            dyn.run(fmax=0.02, steps=300)
            #dyn.run(fmax=0.003, steps=300) # needed for VdW
        except Converged:
            pass # MDTMP: sciopt raises Converged when converged
        self.check_system(system, formula)
        # ase reports one less optimization step
        kwargs = {
            't': time.time() - t0,
            }
        kwargsup = self.add2db(system, dyn)
        if kwargsup: kwargs.update(kwargsup)
        self.write_db(os.path.splitext(filename)[0] + '.db', system, **kwargs)

    def retrieve(self, formula, filename):
        return self.retrieve_energy(formula, filename, index=-1)

    def retrieve_energy(self, formula, filename, index=-1):
        system = read(filename, index=index)
        energy = system.get_potential_energy()
        return energy

    def retrieve_distance(self, formula, filename, a0, a1, mic=False, index=-1):
        system = read(filename, index=index)
        return system.get_distance(a0, a1, mic)

    def retrieve_dihedral(self, list, index=-1):
        system = read(filename, index=index)
        return system.get_dihedral(list)

class GeometryPeriodicSystemTest(PeriodicSystemTest):
    """This uses init from PeriodicSystemTest.  """

    def add2db(self, system, iterable):
        kwargs = {
            'iter': iterable.get_number_of_steps() + 1,
            }
        return kwargs

    def run(self, formula, system, filename):
        """Calculate geometry using the default optimizer.

        Warning: determination of total energies by varying atomic
        separations close to the bond length, allowing
        determination of bond length by fitting is unreliable when
        sampled around a bond length far from equilibrium!
        """
        dyn = QuasiNewton(system, logfile=filename+'.log' ,trajectory=filename)
        t0 = time.time()
        try:
            dyn.run(fmax=0.02, steps=300)
            #dyn.run(fmax=0.003, steps=300) # needed for VdW
        except Converged:
            pass # MDTMP: sciopt raises Converged when converged
        self.check_system(system, formula)
        # ase reports one less optimization step
        kwargs = {
            't': time.time() - t0,
            }
        kwargsup = self.add2db(system, dyn)
        if kwargsup: kwargs.update(kwargsup)
        self.write_db(os.path.splitext(filename)[0] + '.db', system, **kwargs)

    def retrieve(self, formula, filename):
        return self.retrieve_energy(formula, filename, index=-1)

    def retrieve_energy(self, formula, filename, index=-1):
        system = read(filename, index=index)
        energy = system.get_potential_energy()
        return energy

    def retrieve_distance(self, formula, filename, a0, a1, mic=False, index=-1):
        system = read(filename, index=index)
        return system.get_distance(a0, a1, mic)

    def retrieve_dihedral(self, list, index=-1):
        system = read(filename, index=index)
        return system.get_dihedral(list)


class BondLengthDimerTest(GeometryMoleculeTest):
    """This uses init from MoleculeTest.  """

    def retrieve(self, formula, filename):
        return self.retrieve_distance(formula, filename, mic=False, index=-1)

    def retrieve_energy(self, formula, filename, index=-1):
        system = read(filename, index=index)
        if len(system) != 2:
            raise ValueError('Not a dimer')
        energy = system.get_potential_energy()
        return energy

    def retrieve_distance(self, formula, filename, mic=False, index=-1):
        system = read(filename, index=index)
        if len(system) != 2:
            raise ValueError('Not a dimer')
        return system.get_distance(0, 1, mic)

class BondLengthFitDimerTest(GeometryMoleculeTest):
    """This uses init from MoleculeTest.  """

    def get_linspace_parameters(self):
        """Owervrite the method to use different sampling.  """
        return 0.96, 1.04, 5

    def add2db(self, system, iterable):
        kwargs = {
            'iter': iterable.get_number_of_iterations(),
            }
        return kwargs

    def run(self, formula, system, filename):
        """Calculate bond length of a dimer.

        This will calculate total energies for varying atomic
        separations close to the starting bond length, allowing
        determination of bond length by fitting.
        Warning: see the GeometryMoleculeTest class.
        """
        start, stop, num = self.get_linspace_parameters()
        interval = np.linspace(start, stop, num)
        traj = PickleTrajectory(open(filename, 'w'), 'w')
        if len(system) != 2:
            t0 = time.time()
            system.get_potential_energy()
            self.check_system(system, formula)
            kwargs = {
                't': time.time() - t0,
                }
            kwargsup = self.add2db(system, system.get_calculator())
            if kwargsup: kwargs.update(kwargsup)
            self.write_db(os.path.splitext(filename)[0] + '.db', system, **kwargs)
            # Won't create .bak file:
            traj.write(system)
            traj.close()
        else:
            d = system.get_distance(0, 1)
            for x in interval:
                system.set_distance(0, 1, d * x)
                system.get_potential_energy()
                self.check_system(system, formula)
                traj.write(system)
            traj.close()

    def retrieve(self, formula, filename):
        traj = PickleTrajectory(filename, 'r')
        distances = np.array([a.get_distance(0, 1) for a in traj])
        energies = np.array([a.get_potential_energy() for a in traj])
        polynomial = np.polyfit(distances, energies, 3) # or maybe 3rd order?
        # With 3rd order it is not always obvious which root is right
        dedb = np.polyder(polynomial, 1)
        d0 = np.roots(dedb)[1]
        e0 = np.polyval(energies, d0)
        return distances, energies, d0, e0, polynomial

class BondLengthTest(BondLengthFitDimerTest):
    """This uses init from MoleculeTest.  """

    def get_linspace_parameters(self):
        "Deprecated class, may give wrong results! Use BondLengthFitDimerTest instead."
        import warnings
        warnings.warn('BondLengthTest is deprecated. '
                      ' Please use BondLengthFitDimerTest' \
                      ' instead.', DeprecationWarning, stacklevel=2)
        return BondLengthFitDimerTest.get_linspace_parameters(self)

    def run(self, formula, system, filename):
        "Deprecated class, may give wrong results! Use BondLengthFitDimerTest instead."
        import warnings
        warnings.warn('BondLengthTest is deprecated. '
                      ' Please use BondLengthFitDimerTest' \
                      ' instead.', DeprecationWarning, stacklevel=2)
        return BondLengthFitDimerTest.run(self, formula, system, filename)

    def retrieve(self, formula, filename):
        "Deprecated class, may give wrong results! Use BondLengthFitDimerTest instead."
        import warnings
        warnings.warn('BondLengthTest is deprecated. '
                      ' Please use BondLengthFitDimerTest' \
                      ' instead.', DeprecationWarning, stacklevel=2)
        return BondLengthFitDimerTest.retrieve(self, formula, filename)

def write_db(filename, system, **kwargs):
    """Write db file.  """
    import os
    from ase.dft.kpoints import get_monkhorst_pack_size_and_offset
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
        # the ones below may not exist in all calculators
        params.update({
            # calculator
            #"kpts": get_monkhorst_pack_size_and_offset(calc.get_ibz_k_points())[0], # MDTMP - fails sometimes!
            "width": calc.get_electronic_temperature(), # eV
            "niter": calc.get_number_of_iterations(),
            "nelect": calc.get_number_of_electrons(),
            "nbands": calc.get_number_of_bands(),
            })
        # bandgap in case of no fractional occupancies
        occ_0 = calc.get_occupation_numbers(spin=0)
        if ((occ_0.round() - occ_0)**2).sum() < 1e-9:
            params.update({
                # calculator
                "direct_bandgap_0":
                get_bandgap(system, mode='direct', spin=0),
                "indirect_bandgap_0":
                get_bandgap(system, mode='indirect', spin=0),
                })
        if calc.get_spin_polarized():
            occ_1 = calc.get_occupation_numbers(spin=1)
            if ((occ_1.round() - occ_1)**2).sum() < 1e-9:
                params.update({
                    # calculator
                    "direct_bandgap_1":
                    get_bandgap(system, mode='direct', spin=1),
                    "indirect_bandgap_1":
                    get_bandgap(system, mode='indirect', spin=1),
                    })
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

def get_atomization_reactions(data):
    """Get definitions of atomization reactions in the database data.  """
    from ase.atoms import string2symbols
    reactions = []
    rno = 0
    systemkeys = data.keys()
    systemkeys.sort()
    for system in systemkeys:
        if len(data[system]['positions']) == 1:
            continue # skip atom
        else:
            rno += 1
            energy = None
            stoich = []
            # find atoms names in the database
            for a in string2symbols(data[system]['symbols']):
                for s in systemkeys:
                    if data[s]['symbols'] == a:
                        stoich.append((s, 1))
            stoich.append((system, - 1))
            # for atomization reactions without reference data
            # use system name as reaction identifier
            stoich.append(('reaction_id', system))
            reactions.append(stoich)
    return reactions

def calculate_reaction_energies(energies, reactions):
    """Calculate reaction energies according to reactions definitions.  """
    re = {}
    for rno, r in enumerate(reactions):
        assert r[-1][0] == 'reaction_id'
        reaction_id = r[-1][1]
        stoich = r[:-1]
        try:
            energy = 0.0
            for compound, weight in stoich:
                energy += float(energies[compound])*float(weight)
            re[reaction_id] = energy
        except KeyError:
            continue # compound not found in energies
    return re

def write_csv(data, dir='.', outfilename='data.csv'):
    """Write csv file using data in a list or dictionary.   """
    import os
    import csv
    if type(data) == type({}):
        keys = data.keys()
        keys.sort()
        d = []
        for k in keys:
           d.append([k, data[k]])
    elif type(data) == type([]):
        d = [l for l in data]
    csvwriter = csv.writer(open(os.path.join(dir, outfilename), 'wb'))
    for row in d:
        csvwriter.writerow(row)

def write_dict(results, key, name='data', outfilename='out.py'):
    """Write dictionary name of results of key into python file.  """
    import datetime
    import pprint
    try:
        from gpaw.output import initialize_text_stream
        fh, firsttime = initialize_text_stream(outfilename, rank=0)
    except ImportError:
        fh = open(outfilename, 'w')
    fh.write('# Computer generated code! Hands off!\n')
    fh.write('# Generated: ' + str(datetime.date.today()) + '\n')
    fh.write('from numpy import array\n')
    fh.write(name + ' = ')
    data = {}
    data[key] = results
    pprint.pprint(data, stream=fh)
    fh.close()

def read_dict(modulepath, key, name='data'):
    """Read dictionary name from a python module.  """
    import sys
    import __builtin__
    __builtin__.__import__(modulepath)
    f = 'sys.modules[modulepath].' + name + '[key]'
    return eval(f)

def import_module(modulepath):
    """Import a python module.  """
    __builtin__.__import__(modulepath)
    return sys.modules[modulepath]

def break_cell_symmetry(atoms):
    import numpy as np
    cell = atoms.get_cell()
    if len(cell.shape) == 2: # general cell
        cell = cell + np.array([
            [0.00, 0.00, 0.00],
            [0.00, 0.01, 0.00],
            [0.00, 0.00, 0.02]])
    else: # diagonal cell
        cell = cell + np.array([0.0, 0.01, 0.02])
    return cell
