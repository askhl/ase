from ase.utils.compound_test import BatchTest

from ase.utils.compound_test import break_cell_symmetry

from ase.utils.compound_test import PeriodicSystemTest
from ase.utils.compound_test import EnergyPeriodicSystemTest, \
     GeometryPeriodicSystemTest

from ase.calculators.abinit import Abinit as Calculator

from ase.dft.kpoints import get_kpoints_guess

class ABINITPeriodicSystemTest(PeriodicSystemTest):
    def __init__(self, name='abinit', vacuum=0.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        PeriodicSystemTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                              exceptions=exceptions, data=data)
        self.parameters = kwargs
        assert kwargs.get('kpts', None), 'set kpts'
        assert kwargs.get('ecut', None), 'set ecut'

    def setup_calculator(self, system, formula):
        import os
        from glob import glob
        assert self.parameters.get('kpts', None), 'set kpts'
        assert self.parameters.get('ecut', None), 'set ecut'
        # clean old restarts to save disk space
        if 0:
            for ext in ['o_DDB', 'o_DEN', 'o_EIG', 'o_WFK']:
                for file in glob(os.path.join(self.dir, '*.' + ext)):
                    if os.path.isfile(file): os.remove(file)
        # set k-points
        self.parameters['kpts'] = get_kpoints_guess(system,
                                                    self.parameters['kpts'])
        # set output file
        self.parameters['label'] = self.get_filename(formula, extension='')
        #
        calc = Calculator(**self.parameters)

        return calc


class ABINITEnergyPeriodicSystemTest(EnergyPeriodicSystemTest, ABINITPeriodicSystemTest):
    """This uses init from ABINITPeriodicSystemTest.  """

    def __init__(self, name='abinit', vacuum=0.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        ABINITPeriodicSystemTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                          exceptions=exceptions, data=data,
                                          **kwargs)

    def setup_calculator(self, system, formula):
        return ABINITPeriodicSystemTest.setup_calculator(self, system, formula)

    def check_system(self, system, formula):
        pass


class ABINITGeometryPeriodicSystemTest(GeometryPeriodicSystemTest, ABINITPeriodicSystemTest):
    """This uses init from ABINITPeriodicSystemTest.  """

    def __init__(self, name='abinit', vacuum=0.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        ABINITPeriodicSystemTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                          exceptions=exceptions, data=data,
                                          **kwargs)

    def setup_calculator(self, system, formula):
        return ABINITPeriodicSystemTest.setup_calculator(self, system, formula)

    def check_system(self, system, formula):
        pass
