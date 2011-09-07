import os

import shutil

from ase.utils.compound_test import break_cell_symmetry

from ase.utils.compound_test import PeriodicSystemTest
from ase.utils.compound_test import EnergyPeriodicSystemTest, \
     GeometryPeriodicSystemTest

from ase.calculators.vasp import Vasp as Calculator

from ase.dft.kpoints import get_kpoints_guess

class VASPPeriodicSystemTest(PeriodicSystemTest):
    def __init__(self, name='vasp', vacuum=0.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        PeriodicSystemTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                    exceptions=exceptions, data=data)
        self.parameters = kwargs
        assert kwargs.get('kpts', None), 'set kpts'
        if not kwargs.get('setups', None): # set better setups
            # http://cms.mpi.univie.ac.at/vasp/vasp/node243.html
            # http://cms.mpi.univie.ac.at/vasp/vasp/PAW_potentials.html
            setups = {
                'B':'_h', 'C':'_h', 'N':'_h', 'O':'_h', 'F':'_h',
                'H':'_h',
                'Li': '_sv', 'Be': '_sv',
                'Na': '_sv', 'Mg': '_pv',
                'K': '_sv', 'Ca': '_sv',
                'Rb': '_sv', 'Sr': '_sv',
                'Cs': '_sv', 'Ba': '_sv',
                'Sc': '_sv',
                'Ti': '_pv', 'V': '_pv', 'Cr': '_pv', 'Mn': '_pv',
                'Fe': '_pv', 'Ni': '_pv', 'Cu': '_pv',
                'Y': '_sv',
                'Zr': '_sv', 'Nb': '_pv', 'Mo': '_pv', 'Tc': '_pv',
                'Ru': '_pv', 'Rh': '_pv', 'Pd': '_pv',
                'Hf': '_pv', 'Ta': '_pv', 'W': '_pv', 'Re': '_pv',
                'Os': '_pv',
                # https://listserv.fysik.dtu.dk/pipermail/ase-users/2011-January/000906.html
                #'Al': '_h',  'Si': '_h', 'P': '_h', 'S': '_h', 'Cl': '_h',
                'Si': '_h', 'P': '_h', 'S': '_h', 'Cl': '_h',
                'Ga': '_h', 'Ge': '_h',
                'In': '_d', 'Sn': '_d',
                'Tl': '_d', 'Pb': '_d', 'Bi': '_d',
                }
            self.parameters['setups'] = setups

    def setup_calculator(self, system, formula):
        import os
        assert self.parameters.get('kpts', None), 'set kpts'
        # clean old vasp restarts to avoid conflicts
        for file in ['CHG', 'CHGCAR', 'WAVECAR']:
            if os.path.isfile(file): os.remove(file)
        # set k-points
        self.parameters['kpts'] = get_kpoints_guess(system,
                                                    self.parameters['kpts'])
        calc = Calculator(**self.parameters)
        # calculate number of electrons and bands
        nupdown = sum(system.get_initial_magnetic_moments())
        self.parameters['nupdown'] = nupdown

        return calc


class VASPEnergyPeriodicSystemTest(EnergyPeriodicSystemTest, VASPPeriodicSystemTest):
    """This uses init from VASPPeriodicSystemTest and run from EnergyPeriodicSystemTest.  """

    def __init__(self, name='vasp', vacuum=0.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        VASPPeriodicSystemTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                        exceptions=exceptions, data=data,
                                        **kwargs)

    def setup_calculator(self, system, formula):
        return VASPPeriodicSystemTest.setup_calculator(self, system, formula)

    def run(self, formula, system, filename):
        EnergyPeriodicSystemTest.run(self, formula, system, filename)
        # save CAR files
        for f in ['OUTCAR', 'INCAR', 'POSCAR', 'POTCAR', 'KPOINTS']:
            fname = VASPPeriodicSystemTest.get_filename(self,
                                                        formula, extension=f)
            if os.path.isfile(f):
                shutil.move(f, fname)

class VASPGeometryPeriodicSystemTest(GeometryPeriodicSystemTest, VASPPeriodicSystemTest):
    """This uses init from VASPPeriodicSystemTest and run from GeometryPeriodicSystemTest.  """

    def __init__(self, name='vasp', vacuum=0.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        VASPPeriodicSystemTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                        exceptions=exceptions, data=data,
                                        **kwargs)

    def setup_calculator(self, system, formula):
        return VASPPeriodicSystemTest.setup_calculator(self, system, formula)

    def run(self, formula, system, filename):
        GeometryPeriodicSystemTest.run(self, formula, system, filename)
        # save CAR files
        for f in ['OUTCAR', 'INCAR', 'POSCAR', 'POTCAR', 'KPOINTS']:
            fname = VASPPeriodicSystemTest.get_filename(self,
                                                        formula, extension=f)
            if os.path.isfile(f):
                shutil.move(f, fname)
