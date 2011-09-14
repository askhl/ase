import os

import shutil

from ase.utils.compound_test import break_cell_symmetry

from ase.utils.compound_test import MoleculeTest
from ase.utils.compound_test import EnergyMoleculeTest, \
     GeometryMoleculeTest, BondLengthDimerTest

from ase.calculators.vasp import Vasp as Calculator

from ase.dft.kpoints import get_kpoints_guess

class VASPMoleculeTest(MoleculeTest):
    def __init__(self, name='vasp', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        MoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                              exceptions=exceptions, data=data)
        self.parameters = kwargs
        if not kwargs.get('kpts', None): # always set kpts
            self.parameters['kpts'] = [1, 1, 1]
        # leaking electons problem for atoms
        # http://cms.mpi.univie.ac.at/vasp-forum/forum_viewtopic.php?3.92
        # ferwe/ferdo is the solution, but this forces one to set nbands
        #assert kwargs.get('nbands', None), 'set nbands' # MDTMP
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
        # http://cms.mpi.univie.ac.at/vasp-forum/forum_viewtopic.php?4.8980.0#9073
        if len(system) == 1: # break cell symmetry for atoms
            cell = break_cell_symmetry(system)
            system.set_cell(cell)
        system.center()
        #
        if len(system) > 1:
            # http://cms.mpi.univie.ac.at/vasp/guide/node143.html
            self.parameters['ldipol'] = True
            self.parameters['idipol'] = 4
            self.parameters['dipol'] = [0.5, 0.5, 0.5]
        calc = Calculator(**self.parameters)
        # calculate number of electrons and bands
        nupdown = sum(system.get_initial_magnetic_moments())
        self.parameters['nupdown'] = nupdown
        calc.initialize(system)
        potcar = 'POTCAR'
        calc.write_potcar()
        nelect = calc.get_default_number_of_electrons(filename=potcar)
        if os.path.isfile(potcar): os.remove(potcar)
        nelect_sum = 0.0
        for atom in system:
            for s, n in nelect:
                if atom.symbol == s:
                    nelect_sum += n
        assert nelect_sum > 0.0
        # set the number of electrons
        nelect = nelect_sum-sum(system.get_charges())
        calc.set(nelect=nelect)
        # set the number of bands
        nbands = self.parameters.get('nbands', None)
        if not nbands is None:
            assert not nbands is None, 'set nbands'
            if nbands < 0:
                nbands = int(0.5*nelect) - nbands
            calc.set(nbands=nbands)
            # set ferwe/ferdo for spin polarized systems to avoid missing electrons
            # http://cms.mpi.univie.ac.at/vasp-forum/forum_viewtopic.php?3.92
            if (system.get_initial_magnetic_moments().any() and
                self.parameters.get('sigma', 0.2) == 0.0): # only with width of 0.0
                # used in hack to set ferwe/ferdo to assure occupations
                assert abs(nbands) > 0
                # http://cms.mpi.univie.ac.at/vasp-forum/forum_viewtopic.php?4.1434
                nferwe = int((nelect + nupdown)/2.)
                ferwe = [1 for a in range(nferwe)]
                for n in range(nbands - nferwe):
                    ferwe.append(0)
                ferdo = [1 for a in range(int(nferwe - nupdown))]
                for n in range(nbands - int(nferwe - nupdown)):
                    ferdo.append(0)
                calc.set(ferwe=ferwe)
                calc.set(ferdo=ferdo)

        return calc


class VASPEnergyMoleculeTest(EnergyMoleculeTest, VASPMoleculeTest):
    """This uses init from VASPMoleculeTest and run from EnergyMoleculeTest.  """

    def __init__(self, name='vasp', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        VASPMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                  exceptions=exceptions, data=data,
                                  **kwargs)

    def setup_calculator(self, system, formula):
        return VASPMoleculeTest.setup_calculator(self, system, formula)

    def run(self, formula, system, filename):
        EnergyMoleculeTest.run(self, formula, system, filename)
        # save CAR files
        for f in ['OUTCAR', 'INCAR', 'POSCAR', 'POTCAR', 'KPOINTS']:
            fname = VASPMoleculeTest.get_filename(self,
                                                  formula, extension=f)
            if os.path.isfile(f):
                shutil.move(f, fname)

class VASPGeometryMoleculeTest(GeometryMoleculeTest, VASPMoleculeTest):
    """This uses init from VASPMoleculeTest and run from GeometryMoleculeTest.  """

    def __init__(self, name='vasp', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        VASPMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                  exceptions=exceptions, data=data,
                                  **kwargs)

    def setup_calculator(self, system, formula):
        return VASPMoleculeTest.setup_calculator(self, system, formula)

    def run(self, formula, system, filename):
        GeometryMoleculeTest.run(self, formula, system, filename)
        # save CAR files
        for f in ['OUTCAR', 'INCAR', 'POSCAR', 'POTCAR', 'KPOINTS']:
            fname = VASPMoleculeTest.get_filename(self,
                                                  formula, extension=f)
            if os.path.isfile(f):
                shutil.move(f, fname)

class VASPBondLengthDimerTest(BondLengthDimerTest, VASPMoleculeTest):
    """This uses init from VASPMoleculeTest and run from BondLengthDimerTest.  """

    def __init__(self, name='vasp', vacuum=8.0, cell=None,
                 exceptions=(RuntimeError), data=None,
                 **kwargs):
        VASPMoleculeTest.__init__(self, name=name, vacuum=vacuum, cell=cell,
                                  exceptions=exceptions, data=data,
                                  **kwargs)

    def setup_calculator(self, system, formula):
        return VASPMoleculeTest.setup_calculator(self, system, formula)

    def run(self, formula, system, filename):
        BondLengthDimerTest.run(self, formula, system, filename)
        # save CAR files
        for f in ['OUTCAR', 'INCAR', 'POSCAR', 'POTCAR', 'KPOINTS']:
            fname = VASPMoleculeTest.get_filename(self,
                                                  formula, extension=f)
            if os.path.isfile(f):
                shutil.move(f, fname)
