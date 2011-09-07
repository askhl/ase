from ase.utils.compound_test import BatchTest

from ase.utils.compound_test import break_cell_symmetry

from ase.utils.compound_test import MoleculeTest
from ase.utils.compound_test import EnergyMoleculeTest, \
     GeometryMoleculeTest, BondLengthDimerTest, \
     BondLengthTest

from ase.calculators.elk import ELK as Calculator

from ase.dft.kpoints import get_kpoints_guess

class ELKMoleculeTest(MoleculeTest):
    def __init__(self, name='elk', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        MoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                              exceptions=exceptions, data=data)
        self.parameters = kwargs
        if not kwargs.get('kpts', None): # always set kpts
            self.parameters['kpts'] = [1, 1, 1]
            self.parameters['autokpt'] = False

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
        # http://cms.mpi.univie.ac.at/vasp-forum/forum_viewtopic.php?4.8980.0#9073
        if len(system) == 1: # break cell symmetry for atoms
            cell = break_cell_symmetry(system)
            system.set_cell(cell)
        system.center()
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

class ELKEnergyMoleculeTest(EnergyMoleculeTest, ELKMoleculeTest):
    """This uses init from ELKMoleculeTest and run from EnergyMoleculeTest.  """

    def __init__(self, name='elk', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        ELKMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                 exceptions=exceptions, data=data,
                                 **kwargs)

    def setup_calculator(self, system, formula):
        return ELKMoleculeTest.setup_calculator(self, system, formula)

    def check_system(self, system, formula):
        pass

    def add2db(self, system, iterable):
        pass

    def write_db(self, filename, system, **kwargs):
        pass


class ELKGeometryMoleculeTest(GeometryMoleculeTest, ELKMoleculeTest):
    """This uses init from ELKMoleculeTest and run from GeometryMoleculeTest.  """

    def __init__(self, name='elk', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        ELKMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                 exceptions=exceptions, data=data,
                                 **kwargs)

    def setup_calculator(self, system, formula):
        return ELKMoleculeTest.setup_calculator(self, system, formula)

    def check_system(self, system, formula):
        pass

    def add2db(self, system, iterable):
        pass

    def write_db(self, filename, system, **kwargs):
        pass


class ELKBondLengthFitDimerTest(BondLengthTest, ELKMoleculeTest):
    """This uses init from ELKMoleculeTest and run from BondLengthTest.  """

    def __init__(self, name='elk', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        ELKMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                 exceptions=exceptions, data=data,
                                 **kwargs)

    def setup_calculator(self, system, formula):
        return ELKMoleculeTest.setup_calculator(self, system, formula)

    def check_system(self, system, formula):
        pass

    def add2db(self, system, iterable):
        pass

    def write_db(self, filename, system, **kwargs):
        pass


class ELKBondLengthDimerTest(BondLengthDimerTest, ELKMoleculeTest):
    """This uses init from ELKMoleculeTest and run from BondLengthDimerTest.  """

    def __init__(self, name='elk', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        ELKMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                 exceptions=exceptions, data=data,
                                 **kwargs)

    def setup_calculator(self, system, formula):
        return ELKMoleculeTest.setup_calculator(self, system, formula)

    def check_system(self, system, formula):
        pass

    def add2db(self, system, iterable):
        pass

    def write_db(self, filename, system, **kwargs):
        pass
