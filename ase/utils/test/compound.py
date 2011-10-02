"""This module defines extensible classes for running tests.

Use this to compare different calculators, XC functionals and so on by
calculating e.g. atomization energies,
bond lengths across databases of molecules, bulk systems.

One compound test to rule them all
One compound test to run them
One compound test to save them all
And on the webpage plot them (implementation pending)
"""

import os

import time

from ase import Atoms
import ase.io

from ase.structure import molecule
from ase.data.g2 import data as data_g2

from ase.parallel import rank

class CompoundTest:
    """Base class for runnings various tests.  """

    def __init__(self, calculator, calculate, retrieve,
                 name='test', data=None):

        """Create a test which will run all run methods using the calculator.

        The name parameter will be part of all output files generated
        by this test.  If name contains a '/' character, the
        preceding part will be interpreted as a directory in which to
        put files.

        data is a dictionary of formula keys with corresponding atoms objects.
        If data is None then ase.data.g2 is used.  """

        assert name.count(os.path.sep) <= 1, 'More than 1 separator in name'

        self.calculator = calculator

        self.calculate = calculate

        self.retrieve = retrieve

        dir, path = os.path.split(name)
        self.dir = dir
        self.identifier = path
        self.name = name

        self.data = {}
        if data is None:
            for k in data_g2.keys():
                self.data[k] = molecule(k, data=data)
        else:
            for k in data.keys():
                if not isinstance(data[k], Atoms):
                    self.data[k] = molecule(k, data=data)
                else:
                    self.data[k] = data[k]

    def get_filename(self, formula, extension=''):
        """Returns the filename for a test result file.

        Default format is <name>.<formula>[.<extension>]
        """

        if extension == '':
            return '.'.join([self.name, formula])
        else:
            return '.'.join([self.name, formula, extension])

    def get_lock_filename(self, formula):
        """Returns the filename for the main test result file.
        For example a trajectory filename for single-point calculations.

        The test may write other files, but this filename is used as a
        lock file whether the calculation has been done already,
        and it's used to retrieve the results.

        Every implementation has to provide this method."""
        raise NotImplementedError

    def get_formulas(self, natoms=None):
        """Get sorted list of systems with given number of atoms in database.  """
        if natoms is None:
            l = [s for s in self.data.keys()]
        else:
            assert natoms > 0, 'Error: number of atoms must be positive: ' + str(natoms)
            l = [s for s in self.data.keys() if (len(self.data(s)) == natoms)]
        l = list(set(l)) # unique
        l.sort()
        return l

    def setup_calculator(self, formula, system):
        """Return a requested calculator for the given formula.  """
        return self.calculator(formula, system, self)

    def setup_system(self, formula):
        """Create an Atoms object from the given formula.

        Every implementation has to provide this method."""
        raise NotImplementedError

    def setup(self, formula):
        """Build calculator and atoms objects.

        This will invoke the setup_calculator and setup_system methods."""
        system = self.setup_system(formula)
        calc = self.setup_calculator(formula, system)
        system.set_calculator(calc)
        return system

    def run(self, formula, system):
        """Apply all the calculate method to the formula.  """

        self.calculate(formula, system, self)

    def retrieve_results(self, formula):
        """Retrieve results of a calculation.  """
        return self.retrieve(formula, self)

class CompoundTestEnergy(CompoundTest):
    """Perform single point energy calculation. """

    def __init__(self, calculator, calculate=None, retrieve=None,
                 name='test', data=None):

        # provide default methods
        if calculate is None:
            calculate = calculate_energy

        if retrieve is None:
            retrieve = retrieve_energy

        CompoundTest.__init__(self, calculator, calculate, retrieve,
                              name=name, data=data)

    def get_lock_filename(self, formula):
        return self.get_filename(formula, extension='traj')

    def setup_system(self, formula):
        """Return an Atoms object for the given formula.  """
        return self.data[formula]

def calculate_energy(formula, system, test):
    t = time.time()
    system.get_potential_energy()
    t = time.time() - t
    M = sum(test.data[formula].get_initial_magnetic_moments())
    check_occupations(formula, system, M=M)
    ase.io.write(test.get_lock_filename(formula), system, format='traj')
    kwargs = {'t': t}
    kwargs.update({'iter': system.get_calculator().get_number_of_iterations()})
    write_db(test.get_filename(formula, extension='db'), system, **kwargs)

def retrieve_energy(formula, test):
    system = ase.io.read(test.get_lock_filename(formula))
    return system.get_potential_energy()

def check_occupations(formula, system, M=0.0):
    """Check the electronic state (compare occupations with magnetic moment).  """
    calculator = system.get_calculator()
    nelectrons = calculator.get_number_of_electrons()
    nspins = calculator.get_number_of_spins()
    fa = calculator.get_occupation_numbers(spin=0)
    width = calculator.get_electronic_temperature()
    if nspins == 1:
        assert (fa.sum() - nelectrons) < 1e-6, formula+ ' :wrong number of electrons: ' + str(fa.sum())
    if not (width > 0.0):
        assert ((fa.round() - fa)**2).sum() < 1e-14, formula + ' :large fractional occupancies: '  + str(fa)
    if nspins == 2:
        fb = calculator.get_occupation_numbers(spin=1)
        # program may prefer to occupy spin=0 or spin=1 first
        if len([a for a in fa if abs(a) > 0.0]) > len([a for a in fb if abs(a) > 0.0]):
            ftmp = fa
            fa = fb
            fb = ftmp
        if not (width > 0.0):
            assert ((fb.round() - fb)**2).sum() < 1e-9, formula+ ' :large fractional occupancies: ' + str(fb)
        assert abs(abs((fa-fb).sum()) - M) < 1e-9, formula + ' :incorrect magnetic moment: ' + str((fa-fb).sum())

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

def read_db(filename):
    import cmr
    return cmr.read(filename)

def retrieve_db(formula, test):
    return read_db(test.get_filename(formula, extension='db'))

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
