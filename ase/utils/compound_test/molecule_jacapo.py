from ase.utils.compound_test import BatchTest

from ase.utils.compound_test import break_cell_symmetry

from ase.utils.compound_test import MoleculeTest
from ase.utils.compound_test import EnergyMoleculeTest, \
     GeometryMoleculeTest, BondLengthDimerTest

from ase.calculators.jacapo import Jacapo as Calculator

from ase.dft.kpoints import get_kpoints_guess

class JACAPOMoleculeTest(MoleculeTest):
    def __init__(self, name='jacapo', vacuum=8.0, cell=None,
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
        # http://cms.mpi.univie.ac.at/vasp-forum/forum_viewtopic.php?4.8980.0#9073
        if len(system) == 1: # break cell symmetry for atoms
            cell = break_cell_symmetry(system)
            system.set_cell(cell)
        system.center()
        #
        # spin-polarized system
        if system.get_initial_magnetic_moments().any():
            if not 'spinpol' in self.parameters.keys():
                self.parameters['spinpol'] = True
        #
        calc = Calculator(**self.parameters)

        return calc


class JACAPOEnergyMoleculeTest(EnergyMoleculeTest, JACAPOMoleculeTest):
    """This uses init from JACAPOMoleculeTest.  """

    def __init__(self, name='jacapo', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        JACAPOMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                    exceptions=exceptions, data=data,
                                    **kwargs)

    def setup_calculator(self, system, formula):
        return JACAPOMoleculeTest.setup_calculator(self, system, formula)


class JACAPOGeometryMoleculeTest(GeometryMoleculeTest, JACAPOMoleculeTest):
    """This uses init from JACAPOMoleculeTest.  """

    def __init__(self, name='jacapo', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        JACAPOMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                    exceptions=exceptions, data=data,
                                    **kwargs)

    def setup_calculator(self, system, formula):
        return JACAPOMoleculeTest.setup_calculator(self, system, formula)

class JACAPOBondLengthDimerTest(BondLengthDimerTest, JACAPOMoleculeTest):
    """This uses init from JACAPOMoleculeTest.  """

    def __init__(self, name='jacapo', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        JACAPOMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                    exceptions=exceptions, data=data,
                                    **kwargs)

    def setup_calculator(self, system, formula):
        return JACAPOMoleculeTest.setup_calculator(self, system, formula)
