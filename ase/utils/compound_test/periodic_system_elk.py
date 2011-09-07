from ase.utils.compound_test import BatchTest

from ase.utils.compound_test import break_cell_symmetry

from ase.utils.compound_test import PeriodicSystemTest
from ase.utils.compound_test import EnergyPeriodicSystemTest, \
     GeometryPeriodicSystemTest

from ase.calculators.elk import ELK as Calculator

from ase.dft.kpoints import get_kpoints_guess

class ELKPeriodicSystemTest(PeriodicSystemTest):
    def __init__(self, name='elk', vacuum=0.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        PeriodicSystemTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                              exceptions=exceptions, data=data)
        self.parameters = kwargs
        assert (kwargs.get('kpts', None) or
                kwargs.get('autokpt', None)), 'set kpts or autokpt'
        assert not (kwargs.get('kpts', None) and
                    kwargs.get('autokpt', None)), 'set kpts or autokpt'

    def setup_calculator(self, system, formula):
        import os
        from glob import glob
        assert (self.parameters.get('kpts', None) or
                self.parameters.get('autokpt', None)), 'set kpts or autokpt'
        # clean old restarts to avoid conflicts
        files = glob('EVECFV.OUT')
        files.extend(glob('EVECSV.OUT'))
        files.extend(glob('STATE.OUT'))
        for file in files:
            if os.path.isfile(file): os.remove(file)
        # set k-points
        if self.parameters.get('kpts', None):
            self.parameters['kpts'] = get_kpoints_guess(system,
                                                        self.parameters['kpts'])
        # set output file
        self.parameters['dir'] = os.path.join(self.dir, formula)
        #
        # spin-polarized system
        # http://sourceforge.net/projects/elk/forums/forum/897820/topic/4058497
        if system.get_initial_magnetic_moments().any():
            if not 'bfieldc' in self.parameters.keys():
                self.parameters['bfieldc'] = (0.0, 0.0, 0.01)
            if not 'reducebf' in self.parameters.keys():
                self.parameters['reducebf'] = 0.8
            if not 'spinpol' in self.parameters.keys():
                self.parameters['spinpol'] = True
        #
        calc = Calculator(**self.parameters)

        return calc

class ELKEnergyPeriodicSystemTest(EnergyPeriodicSystemTest, ELKPeriodicSystemTest):
    """This uses init from ELKPeriodicSystemTest and run from EnergyPeriodicSystemTest.  """

    def __init__(self, name='elk', vacuum=0.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        ELKPeriodicSystemTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                       exceptions=exceptions, data=data,
                                       **kwargs)

    def setup_calculator(self, system, formula):
        return ELKPeriodicSystemTest.setup_calculator(self, system, formula)

    def check_system(self, system, formula):
        pass

    def add2db(self, system, iterable):
        pass

    def write_db(self, filename, system, **kwargs):
        pass


class ELKGeometryPeriodicSystemTest(GeometryPeriodicSystemTest, ELKPeriodicSystemTest):
    """This uses init from ELKPeriodicSystemTest and run from GeometryPeriodicSystemTest.  """

    def __init__(self, name='elk', vacuum=0.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        ELKPeriodicSystemTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                       exceptions=exceptions, data=data,
                                       **kwargs)

    def setup_calculator(self, system, formula):
        return ELKPeriodicSystemTest.setup_calculator(self, system, formula)

    def check_system(self, system, formula):
        pass

    def add2db(self, system, iterable):
        pass

    def write_db(self, filename, system, **kwargs):
        pass
