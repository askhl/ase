from ase.utils.compound_test import BatchTest

from ase.utils.compound_test import break_cell_symmetry

from ase.utils.compound_test import MoleculeTest
from ase.utils.compound_test import EnergyMoleculeTest, \
     GeometryMoleculeTest, BondLengthDimerTest, \
     BondLengthFitDimerTest

from gpaw import GPAW as Calculator

from gpaw import ConvergenceError
from gpaw.mixer import Mixer

from ase.dft.kpoints import get_kpoints_guess

class GPAWMoleculeTest(MoleculeTest):
    def __init__(self, name='gpaw', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError, ConvergenceError), data=None,
                 **kwargs):
        MoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                              exceptions=exceptions, data=data)
        self.parameters = kwargs
        if not kwargs.get('kpts', None): # always set kpts
            self.parameters['kpts'] = [1, 1, 1]

    def setup_calculator(self, system, formula):
        assert self.parameters.get('kpts', None), 'set kpts'
        # set k-points
        self.parameters['kpts'] = get_kpoints_guess(system,
                                                    self.parameters['kpts'])
        # set output file
        self.parameters['txt'] = self.get_filename(formula, extension='txt')
        #
        self.parameters['hund'] = (len(system) == 1)
        #
        # http://cms.mpi.univie.ac.at/vasp-forum/forum_viewtopic.php?4.8980.0#9073
        if len(system) == 1: # break cell symmetry for atoms
            cell = break_cell_symmetry(system)
            system.set_cell(cell)
        system.center()
        #
        # charge
        if not self.parameters.get('charge', None):
            self.parameters['charge'] = sum(system.get_charges())
        #
        calc = Calculator(**self.parameters)

        return calc


class GPAWEnergyMoleculeTest(EnergyMoleculeTest, GPAWMoleculeTest):
    """This uses init from GPAWMoleculeTest and run from EnergyMoleculeTest.  """

    def __init__(self, name='gpaw', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError, ConvergenceError), data=None,
                 **kwargs):
        GPAWMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                  exceptions=exceptions, data=data,
                                  **kwargs)

    def setup_calculator(self, system, formula):
        return GPAWMoleculeTest.setup_calculator(self, system, formula)

    def run(self, formula, system, filename):
        EnergyMoleculeTest.run(self, formula, system, filename)
        if 0:
            # write also gpw file
            gpwf = GPAWMoleculeTest.get_filename(self,
                                                 formula, extension='gpw')
            calc = system.get_calculator()
            calc.write(gpwf)


class GPAWGeometryMoleculeTest(GeometryMoleculeTest, GPAWMoleculeTest):
    """This uses init from GPAWMoleculeTest and run from GeometryMoleculeTest.  """

    def __init__(self, name='gpaw', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError, ConvergenceError), data=None,
                 **kwargs):
        GPAWMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                  exceptions=exceptions, data=data,
                                  **kwargs)

    def setup_calculator(self, system, formula):
        return GPAWMoleculeTest.setup_calculator(self, system, formula)

    def run(self, formula, system, filename):
        GeometryMoleculeTest.run(self, formula, system, filename)
        if 0:
            # write also gpw file
            gpwf = GPAWMoleculeTest.get_filename(self,
                                                 formula, extension='gpw')
            calc = system.get_calculator()
            calc.write(gpwf)


class GPAWBondLengthFitDimerTest(BondLengthFitDimerTest, GPAWMoleculeTest):
    """This uses init from GPAWMoleculeTest and run from BondLengthFitDimerTest.  """

    def __init__(self, name='gpaw', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError, ConvergenceError), data=None,
                 **kwargs):
        GPAWMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                  exceptions=exceptions, data=data,
                                  **kwargs)

    def setup_calculator(self, system, formula):
        return GPAWMoleculeTest.setup_calculator(self, system, formula)

    def run(self, formula, system, filename):
        BondLengthFitTest.run(self, formula, system, filename)
        if 0:
            # write also gpw file
            gpwf = GPAWMoleculeTest.get_filename(self,
                                                 formula, extension='gpw')
            calc = system.get_calculator()
            calc.write(gpwf)


class GPAWBondLengthDimerTest(BondLengthDimerTest, GPAWMoleculeTest):
    """This uses init from GPAWMoleculeTest and run from BondLengthDimerTest.  """

    def __init__(self, name='gpaw', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError, ConvergenceError), data=None,
                 **kwargs):
        GPAWMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                  exceptions=exceptions, data=data,
                                  **kwargs)

    def setup_calculator(self, system, formula):
        return GPAWMoleculeTest.setup_calculator(self, system, formula)

    def run(self, formula, system, filename):
        BondLengthTest.run(self, formula, system, filename)
        if 0:
            # write also gpw file
            gpwf = GPAWMoleculeTest.get_filename(self,
                                                 formula, extension='gpw')
            calc = system.get_calculator()
            calc.write(gpwf)
