import optparse
from math import sqrt, pi

import numpy as np


class CalculatorWrapper:
    def __init__(self, Class, name, **kwargs):
        self.Class = Class
        self.name = name
        self.kwargs = kwargs

    def add_options(self, parser):
        pass

    def parse(self, opts):
        pass

    def __call__(self, name, atoms):
        return self.Class(**self.kwargs)


class ElectronicStructureCalculatorWrapper(CalculatorWrapper):
    def __init__(self, name, xc='LDA', kpts=None, kptdensity=3.0, **kwargs):
        self.xc = xc
        self.kpts = kpts
        self.kptdensity = kptdensity

        CalculatorWrapper.__init__(self, None, name, **kwargs)

    def add_options(self, parser):
        calc = optparse.OptionGroup(parser, 'Calculator')
        calc.add_option('-f', '--xc', default='LDA',
                        help='XC functional.',
                        choices=['LDA', 'PBE'])
        calc.add_option('-k', '--monkhorst-pack',
                        help='Monkhorst-Pack sampling of BZ.  Example: ' +
                        '"4,4,4": 4x4x4 k-points, "4,4,4g": same set of ' +
                        'k-points shited to include the Gamma point.')
        calc.add_option('--k-point-density', type='float', default=3.0,
                        help='Density of k-points in Angstrom.')
        parser.add_option_group(calc)

    def parse(self, opts):
        self.xc = opts.xc
        
        mp = opts.monkhorst_pack
        if mp is not None:
            if mp[-1].lower() == 'g':
                kpts = np.array([int(k) for k in mp[:-1].split(',')])
                shift = 0.5 * ((kpts + 1) % 2) / kpts
                self.kpts = monkhorst_pack(kpts) + shift
            else:
                self.kpts = [int(k) for k in mp.split(',')]

        self.kptdensity = opts.k_point_density

    def calculate_kpts(self, atoms):
        if self.kpts is not None:
            return self.kpts
        
        recipcell = atoms.get_reciprocal_cell()
        kpts = []
        for i in range(3):
            if atoms.pbc[i]:
                k = 2 * pi * sqrt((recipcell[i]**2).sum()) * self.kptdensity
                kpts.append(max(1, 2 * int(round(k / 2))))
            else:
                kpts.append(1)

        return kpts


# Recognized names of calculators sorted alphabetically:
calcnames = ['abinit', 'emt', 'gpaw']


def get_calculator_wrapper(name):
    if hasattr(name, 'get_potential_energy'):
        return CalculatorWrapper(name)

    if name == 'emt':
        from ase.calculators.emt import EMT
        return CalculatorWrapper(EMT, 'EMT')

    if name == 'abinit':
        from ase.calculators.abinit import AbinitWrapper
        return AbinitWrapper()

    if name == 'gpaw':
        from gpaw.wrapper import GPAWWrapper
        return GPAWWrapper()

    2 / 0
