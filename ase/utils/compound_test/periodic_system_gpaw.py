from ase.utils.compound_test import BatchTest

from ase.utils.compound_test import break_cell_symmetry

from ase.utils.compound_test import PeriodicSystemTest
from ase.utils.compound_test import EnergyPeriodicSystemTest, \
     GeometryPeriodicSystemTest

from gpaw import GPAW as Calculator

from gpaw import ConvergenceError
from gpaw.mixer import Mixer

from ase.dft.kpoints import get_kpoints_guess

class GPAWPeriodicSystemTest(PeriodicSystemTest):
    def __init__(self, name='gpaw', vacuum=0.0, cell=None,
                 exceptions=(RuntimeError, ConvergenceError), data=None,
                 **kwargs):
        PeriodicSystemTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                    exceptions=exceptions, data=data)
        self.parameters = kwargs
        assert kwargs.get('kpts', None), 'set kpts'

    def setup_calculator(self, system, formula):
        assert self.parameters.get('kpts', None), 'set kpts'
        # set k-points
        self.parameters['kpts'] = get_kpoints_guess(system,
                                                    self.parameters['kpts'])
        # set output file
        self.parameters['txt'] = self.get_filename(formula, extension='txt')
        #
        # charge
        if not self.parameters.get('charge', None):
            self.parameters['charge'] = sum(system.get_charges())
        #
        calc = Calculator(**self.parameters)

        return calc

class GPAWEnergyPeriodicSystemTest(EnergyPeriodicSystemTest, GPAWPeriodicSystemTest):
    """This uses init from GPAWPeriodicSystemTest and run from EnergyPeriodicSystemTest.  """

    def __init__(self, name='gpaw', vacuum=0.0, cell=None,
                 exceptions=(RuntimeError, ConvergenceError), data=None,
                 **kwargs):
        GPAWPeriodicSystemTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                        exceptions=exceptions, data=data,
                                        **kwargs)

    def setup_calculator(self, system, formula):
        return GPAWPeriodicSystemTest.setup_calculator(self, system, formula)

    def run(self, formula, system, filename):
        EnergyPeriodicSystemTest.run(self, formula, system, filename)
        if 0:
            # write also gpw file
            gpwf = GPAWPeriodicSystemTest.get_filename(self,
                                                       formula, extension='gpw')
            calc = system.get_calculator()
            calc.write(gpwf)


class GPAWGeometryPeriodicSystemTest(GeometryPeriodicSystemTest, GPAWPeriodicSystemTest):
    """This uses init from GPAWPeriodicSystemTest and run from GeometryPeriodicSystemTest.  """

    def __init__(self, name='gpaw', vacuum=0.0, cell=None,
                 exceptions=(RuntimeError, ConvergenceError), data=None,
                 **kwargs):
        GPAWPeriodicSystemTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                        exceptions=exceptions, data=data,
                                        **kwargs)

    def setup_calculator(self, system, formula):
        return GPAWPeriodicSystemTest.setup_calculator(self, system, formula)

    def run(self, formula, system, filename):
        GeometryPeriodicSystemTest.run(self, formula, system, filename)
        if 0:
            # write also gpw file
            gpwf = GPAWPeriodicSystemTest.get_filename(self,
                                                       formula, extension='gpw')
            calc = system.get_calculator()
            calc.write(gpwf)
