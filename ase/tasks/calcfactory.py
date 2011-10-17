import optparse
from math import sqrt, pi

import numpy as np


def str2dict(s, namespace={}):
    dct = {}
    s = (s + ',').split('=')
    for i in range(len(s) - 1):
        key = s[i]
        m = s[i + 1].rfind(',')
        value = s[i + 1][:m]
        try:
            value = eval(value, namespace)
        except (NameError, SyntaxError):
            pass
        dct[key] = value
        s[i + 1] = s[i + 1][m + 1:]
    return dct


class CalculatorFactory:
    def __init__(self, Class, name, label='label',
                 kpts=None, kptdensity=3.0,
                 **kwargs):
        self.Class = Class
        self.name = name
        self.label = label
        self.kpts = kpts
        self.kptdensity = kptdensity
        self.kwargs = kwargs

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

    def __call__(self, name, atoms):
        kpts = self.calculate_kpts(atoms)
        if kpts != 'no k-points':
            self.kwargs['kpts'] = kpts

        if self.label is not None:
            self.kwargs[self.label] = name

        return self.Class(**self.kwargs)

    def add_options(self, parser):
        calc = optparse.OptionGroup(parser, 'Calculator')
        calc.add_option('-k', '--monkhorst-pack',
                        metavar='K1,K2,K3',
                        help='Monkhorst-Pack sampling of BZ.  Example: ' +
                        '"4,4,4": 4x4x4 k-points, "4,4,4g": same set of ' +
                        'k-points shifted to include the Gamma point.')
        calc.add_option('--k-point-density', type='float', default=3.0,
                        help='Density of k-points in Angstrom.')
        calc.add_option('-p', '--parameters', metavar='key=value,...',
                        help='Comma-separated key=value pairs of ' +
                        'calculator specific parameters.')
        parser.add_option_group(calc)

    def parse(self, opts, args):
        mp = opts.monkhorst_pack
        if mp is not None:
            if mp[-1].lower() == 'g':
                kpts = np.array([int(k) for k in mp[:-1].split(',')])
                shift = 0.5 * ((kpts + 1) % 2) / kpts
                self.kpts = monkhorst_pack(kpts) + shift
            else:
                self.kpts = [int(k) for k in mp.split(',')]

        self.kptdensity = opts.k_point_density

        if opts.parameters:
            self.kwargs.update(str2dict(opts.parameters))


# Recognized names of calculators sorted alphabetically:
calcnames = ['abinit', 'emt', 'gpaw', 'nwchem', 'vasp']


def calculator_factory(name, **kwargs):
    if name == 'gpaw':
        from gpaw.factory import GPAWFactory
        return GPAWFactory(**kwargs)

    classname = {'emt': 'EMT', 'nwchem': 'NWchem'}.get(name, name.title())
    module = __import__('ase.calculators.' + name, fromlist=[classname])
    Class = getattr(module, classname)

    if name in ['emt']:
        kpts = 'no k-points'
    else:
        kpts = None

    if name in ['emt']:
        label = None
    else:
        label = 'label'

    return CalculatorFactory(Class, classname, label, kpts, **kwargs)

