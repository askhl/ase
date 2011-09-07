from ase.utils.compound_test import BatchTest

from ase.utils.compound_test import break_cell_symmetry

from ase.utils.compound_test import MoleculeTest
from ase.utils.compound_test import EnergyMoleculeTest, \
     GeometryMoleculeTest, BondLengthDimerTest

from ase.calculators.abinit import Abinit as Calculator

from ase.dft.kpoints import get_kpoints_guess

class ABINITMoleculeTest(MoleculeTest):
    def __init__(self, name='abinit', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        MoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                              exceptions=exceptions, data=data)
        self.parameters = kwargs
        if not kwargs.get('kpts', None): # always set kpts
            self.parameters['kpts'] = [1, 1, 1]
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
        # http://cms.mpi.univie.ac.at/vasp-forum/forum_viewtopic.php?4.8980.0#9073
        if len(system) == 1: # break cell symmetry for atoms
            cell = break_cell_symmetry(system)
            system.set_cell(cell)
        system.center()
        #
        calc = Calculator(**self.parameters)

        return calc


class ABINITEnergyMoleculeTest(EnergyMoleculeTest, ABINITMoleculeTest):
    """This uses init from ABINITMoleculeTest.  """

    def __init__(self, name='abinit', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        ABINITMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                    exceptions=exceptions, data=data,
                                    **kwargs)

    def setup_calculator(self, system, formula):
        return ABINITMoleculeTest.setup_calculator(self, system, formula)

#    def check_system(self, system, formula):
#        pass
#
#    def add2db(self, system, iterable):
#        pass
#
#    def write_db(self, filename, system, **kwargs):
#        pass


class ABINITGeometryMoleculeTest(GeometryMoleculeTest, ABINITMoleculeTest):
    """This uses init from ABINITMoleculeTest.  """

    def __init__(self, name='abinit', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        ABINITMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                    exceptions=exceptions, data=data,
                                    **kwargs)

    def setup_calculator(self, system, formula):
        return ABINITMoleculeTest.setup_calculator(self, system, formula)

    def check_system(self, system, formula):
        pass

    def add2db(self, system, iterable):
        pass

    def write_db(self, filename, system, **kwargs):
        pass


class ABINITBondLengthDimerTest(BondLengthDimerTest, ABINITMoleculeTest):
    """This uses init from ABINITMoleculeTest.  """

    def __init__(self, name='abinit', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        ABINITMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                    exceptions=exceptions, data=data,
                                    **kwargs)

    def setup_calculator(self, system, formula):
        return ABINITMoleculeTest.setup_calculator(self, system, formula)

    def check_system(self, system, formula):
        pass

    def add2db(self, system, iterable):
        pass

    def write_db(self, filename, system, **kwargs):
        pass
