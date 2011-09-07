from ase.utils.compound_test import BatchTest

from ase.utils.compound_test import break_cell_symmetry

from ase.utils.compound_test import PeriodicSystemTest
from ase.utils.compound_test import EnergyPeriodicSystemTest, \
     GeometryPeriodicSystemTest

from ase.calculators.jacapo import Jacapo

from ase.dft.kpoints import get_kpoints_guess

class JACAPOPeriodicSystemTest(PeriodicSystemTest):
    def __init__(self, name='jacapo', vacuum=0.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        PeriodicSystemTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                    exceptions=exceptions, data=data)
        self.parameters = kwargs
        assert kwargs.get('kpts', None), 'set kpts'

    def setup_calculator(self, system, formula):
        import os
        from glob import glob
        assert self.parameters.get('kpts', None), 'set kpts'
        # clean old restarts to save disk space
        if 0:
            for file in glob(os.path.join(self.dir, '*.nc')):
                if os.path.isfile(file): os.remove(file)
        # set k-points
        self.parameters['kpts'] = get_kpoints_guess(system,
                                                    self.parameters['kpts'])
        # set output file
        self.parameters['nc'] = self.get_filename(formula, extension='nc')
        #
        # spin-polarized system
        if system.get_initial_magnetic_moments().any():
            if not 'spinpol' in self.parameters.keys():
                self.parameters['spinpol'] = True
        #
        calc = Calculator(**self.parameters)

        return calc


class JACAPOEnergyPeriodicSystemTest(EnergyPeriodicSystemTest, JACAPOPeriodicSystemTest):
    """This uses init from JACAPOPeriodicSystemTest.  """

    def __init__(self, name='jacapo', vacuum=0.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        JACAPOPeriodicSystemTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                          exceptions=exceptions, data=data,
                                          **kwargs)

    def setup_calculator(self, system, formula):
        return JACAPOPeriodicSystemTest.setup_calculator(self, system, formula)


class JACAPOGeometryPeriodicSystemTest(GeometryPeriodicSystemTest, JACAPOPeriodicSystemTest):
    """This uses init from JACAPOPeriodicSystemTest.  """

    def __init__(self, name='jacapo', vacuum=0.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        JACAPOPeriodicSystemTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                          exceptions=exceptions, data=data,
                                          **kwargs)

    def setup_calculator(self, system, formula):
        return JACAPOPeriodicSystemTest.setup_calculator(self, system, formula)
