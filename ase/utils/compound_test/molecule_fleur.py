from ase.utils.compound_test import BatchTest

from ase.utils.compound_test import break_cell_symmetry

from ase.utils.compound_test import MoleculeTest
from ase.utils.compound_test import EnergyMoleculeTest, \
     GeometryMoleculeTest, BondLengthDimerTest, \
     BondLengthTest

from ase.calculators.fleur import FLEUR as Calculator

from ase.dft.kpoints import get_kpoints_guess

class FLEURMoleculeTest(MoleculeTest):
    def __init__(self, name='fleur', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        MoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                              exceptions=exceptions, data=data)
        self.parameters = kwargs
        if not kwargs.get('kpts', None): # always set kpts
            self.parameters['kpts'] = [1, 1, 1]

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
        # http://cms.mpi.univie.ac.at/vasp-forum/forum_viewtopic.php?4.8980.0#9073
        if len(system) == 1: # break cell symmetry for atoms
            cell = break_cell_symmetry(system)
            system.set_cell(cell)
        system.center()
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

class FLEUREnergyMoleculeTest(EnergyMoleculeTest, FLEURMoleculeTest):
    """This uses init from FLEURMoleculeTest and run from EnergyMoleculeTest.  """

    def __init__(self, name='fleur', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        FLEURMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                   exceptions=exceptions, data=data,
                                   **kwargs)

    def setup_calculator(self, system, formula):
        return FLEURMoleculeTest.setup_calculator(self, system, formula)

    def check_system(self, system, formula):
        pass

    def add2db(self, system, iterable):
        pass

    def write_db(self, filename, system, **kwargs):
        pass


class FLEURGeometryMoleculeTest(GeometryMoleculeTest, FLEURMoleculeTest):
    """This uses init from FLEURMoleculeTest and run from GeometryMoleculeTest.  """

    def __init__(self, name='fleur', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        FLEURMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                   exceptions=exceptions, data=data,
                                   **kwargs)

    def setup_calculator(self, system, formula):
        return FLEURMoleculeTest.setup_calculator(self, system, formula)

    def check_system(self, system, formula):
        pass

    def add2db(self, system, iterable):
        pass

    def write_db(self, filename, system, **kwargs):
        pass


class FLEURBondLengthFitDimerTest(BondLengthTest, FLEURMoleculeTest):
    """This uses init from FLEURMoleculeTest and run from BondLengthTest.  """

    def __init__(self, name='fleur', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        FLEURMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                   exceptions=exceptions, data=data,
                                   **kwargs)

    def setup_calculator(self, system, formula):
        return FLEURMoleculeTest.setup_calculator(self, system, formula)

    def check_system(self, system, formula):
        pass

    def add2db(self, system, iterable):
        pass

    def write_db(self, filename, system, **kwargs):
        pass


class FLEURBondLengthDimerTest(BondLengthDimerTest, FLEURMoleculeTest):
    """This uses init from FLEURMoleculeTest and run from BondLengthDimerTest.  """

    def __init__(self, name='fleur', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        FLEURMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                   exceptions=exceptions, data=data,
                                   **kwargs)

    def setup_calculator(self, system, formula):
        return FLEURMoleculeTest.setup_calculator(self, system, formula)

    def check_system(self, system, formula):
        pass

    def add2db(self, system, iterable):
        pass

    def write_db(self, filename, system, **kwargs):
        pass
