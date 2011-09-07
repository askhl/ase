from ase.utils.compound_test import BatchTest

from ase.utils.compound_test import break_cell_symmetry

from ase.utils.compound_test import PeriodicSystemTest
from ase.utils.compound_test import EnergyPeriodicSystemTest, \
     GeometryPeriodicSystemTest

from ase.calculators.fleur import FLEUR as Calculator

from ase.dft.kpoints import get_kpoints_guess

class FLEURPeriodicSystemTest(PeriodicSystemTest):
    def __init__(self, name='fleur', vacuum=0.0, cell=None,
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
        # clean old restarts to avoid conflicts
        files = glob('cdn*')
        files.extend(glob('eig*'))
        files.extend(glob('broyd*'))
        files.extend(glob('pot*'))
        for file in files:
            if os.path.isfile(file): os.remove(file)
        # set k-points
        self.parameters['kpts'] = get_kpoints_guess(system,
                                                    self.parameters['kpts'])
        # set output file
        self.parameters['workdir'] = os.path.join(self.dir, formula)
        #
        # avoid ptsym: STOP rw_symfile: complex phases not implemented!
        # inp file supports only 6 digits for the position
        if 0: # this was required in fleur v25
            system.positions[0, 0] += 0.000002
        # equivatoms=False in spin-polarized system
        if system.get_initial_magnetic_moments().any():
            if not 'equivatoms' in self.parameters.keys():
                self.parameters['equivatoms'] = False
        #
        calc = Calculator(**self.parameters)

        return calc

class FLEUREnergyPeriodicSystemTest(EnergyPeriodicSystemTest, FLEURPeriodicSystemTest):
    """This uses init from FLEURPeriodicSystemTest and run from EnergyPeriodicSystemTest.  """

    def __init__(self, name='fleur', vacuum=0.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        FLEURPeriodicSystemTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                         exceptions=exceptions, data=data,
                                         **kwargs)

    def setup_calculator(self, system, formula):
        return FLEURPeriodicSystemTest.setup_calculator(self, system, formula)

    def check_system(self, system, formula):
        pass

    def add2db(self, system, iterable):
        pass

    def write_db(self, filename, system, **kwargs):
        pass


class FLEURGeometryPeriodicSystemTest(GeometryPeriodicSystemTest, FLEURPeriodicSystemTest):
    """This uses init from FLEURPeriodicSystemTest and run from GeometryPeriodicSystemTest.  """

    def __init__(self, name='fleur', vacuum=0.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        FLEURPeriodicSystemTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                         exceptions=exceptions, data=data,
                                         **kwargs)

    def setup_calculator(self, system, formula):
        return FLEURPeriodicSystemTest.setup_calculator(self, system, formula)

    def check_system(self, system, formula):
        pass

    def add2db(self, system, iterable):
        pass

    def write_db(self, filename, system, **kwargs):
        pass
