import os
import time

import numpy as np

from ase.io.pupynere import NetCDFFile


def restart(filename, **kwargs):
    calc = Dacapo(filename, **kwargs)
    atoms = calc.get_atoms()
    return atoms, calc


class Dacapo:
    def __init__(self, filename='out.nc', 
                 stay_alive=False, stress=False,
                 atoms=None, **kwargs):
        self.filename = filename
        self.stay_alive = stay_alive
        self.calculate_stress = stress

        try:
            self.nc = NetCDFFile(filename)
        except IOError:
            self.nc = None
            self.parameters = {'ecut': 350, 'xc': 'PW91'}
            self.set_atoms(atoms)
        else:
            self.parameters = self._read_parameters_from_nc_file()
            from ase.io.dacapo import read_dacapo
            self.atoms = read_dacapo(filename)

        self.ncvars = {}
        self.pps = []

        self.set(**kwargs)

    def set_atoms(self, atoms):
        if atoms is None:
            self.atoms = None
        else:
            self.atoms = atoms.copy()

    def set_pp(self, Z, path):
        self.pps.append((Z, path))

    def set(self, name=None, value=None, **kwargs):
        if name:
            kwargs[name] = value
        for name, value in kwargs.items():
            if name in ['ecut', 'xc']:
                self.parameters[name] = value
            else:
                self.ncvars[name] = value

    def calculate(self, atoms):
        old = self.atoms
        new = atoms
        if (self.nc is None or
            old is None or
            len(old) != len(new) or
            (old.positions != new.positions).any() or
            (old.numbers != new.numbers).any() or
            (old.cell != new.cell).any()):
            self.write_input_file(atoms)
            self.run()
            self.nc = NetCDFFile('out.nc')
            self.atoms = atoms.copy()

    def _read_parameters_from_nc_file(self):
        return {'ecut': 350, 'xc': 'PW91'}

    def write_input_file(self, atoms):
        input = NetCDFFile('input.nc', 'w')

        par = self.parameters
        if 'nbands' not in par:
            n = sum([valence[atom.symbol] for atom in atoms])
            par['nbands'] = int(n * 0.65) + 4

        magmoms = atoms.get_initial_magnetic_moments()
        if magmoms.any():
            par['spinpol'] = True
        nspins = 1
        if par.get('spinpol'):
            nspins = 2

        assert 'kpts' not in par

        for name, value in [
            ('dim3', 3),
            ('number_ionic_steps', 1),
            ('number_of_dynamic_atoms', len(atoms)),
            ('dim24', 24),
            ('dim2', 2),
            ('dim5', 5),
            ('dim4', 4),
            ('number_of_spin', nspins),
            ('number_of_bands', par['nbands']),
            ('number_BZ_kpoints', 1)]:
            input.createDimension(name, value)

        ncvars = {
            'ChargeMixing:Pulay_KerkerPrecondition': 'No',
            'ChargeMixing:UpdateCharge': 'Yes',
            'ElectronicMinimization:Method': 'eigsolve',
            'ConvergenceControl:DensityConvergence': 0.0001,
            'ConvergenceControl:AbsoluteEnergyConvergence': 1.0e-05,
            'ConvergenceControl:OccupationConvergence': 0.001,
            'ElectronicBands:OccupationStatistics': 'FermiDirac',
            'ElectronicBands:OccupationStatistics_FermiTemperature': 0.1,
            'ElectronicBands:SpinPolarization': 1,
            'ElectronicBands:NumberOfBands': 5,
            'UseSymmetry(dim3)': 'Off',
            'Date(dim24)': time.asctime(),
            'PlaneWaveCutoff()': par['ecut'],
            'ExcFunctional(dim4)': par['xc'],
            'BZKpoints(number_BZ_kpoints, dim3)': np.zeros((1, 3)),
            'DynamicAtomPositions' +
            '(number_ionic_steps, number_of_dynamic_atoms, dim3)':
                atoms.get_scaled_positions(),
            'DynamicAtomSpecies(number_of_dynamic_atoms, dim2)':
                ['%2s' % symbol
                 for symbol in atoms.get_chemical_symbols()],
            'UnitCell(number_ionic_steps, dim3, dim3)': atoms.get_cell(),
            'InitialAtomicMagneticMoment(number_of_dynamic_atoms)': magmoms,
            'AtomProperty_H:PspotFile':
                '/home/jjmo/dacapo/psp/H/PW91/ch_e9g4.pseudo'}
        
        if self.calculate_stress:
            self.calc.CalculateStress()

        for Z, path in self.pps:
            self.calc.SetPseudoPotential(Z, path)

        ncvars.update(self.ncvars)

        for name, value in ncvars.items():
            i = name.find('(')
            if i >= 0:
                dims = name[i + 1:-1].replace(',', ' ').split()
                if isinstance(value, str):
                    value = list(value)
                value = np.asarray(value)
                var = input.createVariable(name[:i], value.dtype, dims)
                if len(dims) == 0:
                    var.assignValue(value)
                else:
                    var[:] = value
            else:
                name, attr = name.split(':')
                if name not in input.variables:
                    input.createVariable(name, 'c', ())
                input.variables[name].__setattr__(attr, value)
        input.close()

    def run(self):
        if self.filename.endswith('.nc'):
            txt = self.filename[:-2] + 'txt'
        else:
            txt = self.filename + 'txt'
        dacapo = os.environ['DACAPO_SCRIPT']
        locals = {'nc': self.filename, 'txt': txt}
        execfile(dacapo, {}, locals)
        exitcode = locals['exitcode']
        if exitcode != 0:
            raise RuntimeError(('Dacapo exited with exit code: %d.  ' +
                                'Check %s for more information.') %
                               (exitcode, txt))

    def get_atoms(self):
        atoms = self.atoms.copy()
        atoms.set_calculator(self)
        return atoms
    
    def get_potential_energy(self, atoms):
        self.calculate(atoms)
        return self.nc.variables['TotalEnergy'][-1]

    def get_forces(self, atoms):
        self.calculate(atoms)
        return self.nc.variables['DynamicAtomForces'][-1]

    def get_stress(self, atoms):
        #check for stress flag ...
        self.update(atoms)
        stress = np.array(self.calc.GetStress())
        if stress.ndim == 2:
            return stress.ravel()[[0, 4, 8, 5, 2, 1]]
        else:
            return stress

    def calculation_required(self, atoms, quantities):
        if self.calc is None:
            return True

        if atoms != self.get_atoms():
            return True

        return False


valence = {
'H':   1,
'B':   3,
'C':   4,
'N':   5,
'O':   6,
'Li':  1,
'Na':  1,
'K':   9,
'Mg':  8,
'Ca': 10,
'Sr': 10,
'Al':  3,
'Ga': 13,
'Sc': 11,
'Ti': 12,
'V':  13,
'Cr': 14,
'Mn':  7,
'Fe':  8,
'Co':  9,
'Ni': 10,
'Cu': 11,
'Zn': 12,
'Y':  11,
'Zr': 12,
'Nb': 13,
'Mo':  6,
'Ru':  8,
'Rh':  9,
'Pd': 10,
'Ag': 11,
'Cd': 12,
'Pt': 10,
}


if 'nc' in locals():
    import os
    path = '~/dacapo/src/gfortran_fnosecond_underscore_serial/'
    exitcode = os.system('%s/dacapo.run input.nc %s -out %s' % (path, nc, txt))
