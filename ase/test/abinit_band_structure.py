import numpy as np

from ase.test import NotAvailable

try:
    from ase.test.abinit_installed import abinit_installed
    abinit_installed()
except NotAvailable:
    raise NotAvailable('Abinit required')

from ase import Atoms
from ase.units import Ry

from ase.calculators.abinit import Abinit

# http://stackoverflow.com/questions/2158395/flatten-an-irregular-list-of-lists-in-python
def flatten(x):
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

def array_almost_equal(a1, a2, tol=np.finfo(type(1.0)).eps):
    return (np.abs(a1 - a2) < tol).all()

def plot_save(directory_name, out_prefix):
    from os.path import exists, sep
    assert exists(directory_name)
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pylab

    pylab.savefig(directory_name + sep + out_prefix +'.png')


ref = {
    'e_fermi': 0.555928813242,
    'eigs':
       [[-1.30859602, -1.5129536 , -1.70778719, -1.89282468, -2.06752184,
        -2.2324229 , -2.38698363, -2.53120402, -2.66508409, -2.78862383,
        -2.90155112, -3.00386597, -3.09584048, -3.17665833, -3.24713584,
        -3.3067288 , -3.3554372 , -3.39353315, -3.42074455, -3.43707139,
        -3.44251366],
       [-1.30859602, -1.09417022, -0.87049255, -0.63701877, -0.39483735,
        -0.14313194,  0.11728112,  0.38640182,  0.66395805,  0.94967771,
         1.24301655,  1.54397459,  1.85200759,  2.16657132,  2.48657734,
         2.81120929,  3.13937872,  3.46890872,  3.79680604,  4.11545148,
         4.34484354],
       [ 6.37671846,  6.19902804,  6.03031739,  5.8569808 ,  5.68446055,
         5.52228063,  5.37044105,  5.2292139 ,  5.09832709,  4.97805272,
         4.86811868,  4.76879709,  4.67981582,  4.60144701,  4.53341852,
         4.47573036,  4.42865464,  4.39191926,  4.36579632,  4.35001371,
         4.34484354],
       [ 6.37671846,  6.19902804,  6.03031739,  5.87113073,  5.72092382,
         5.58024091,  5.44880987,  5.32635859,  5.21315918,  5.10921165,
         5.01424388,  4.92852798,  4.85179185,  4.78430758,  4.72553097,
         4.67600623,  4.63573336,  4.60416814,  4.58158269,  4.5682491 ,
         4.56362317],
       [ 6.37671846,  6.23358652,  6.04011349,  5.87113073,  5.72092382,
         5.58024091,  5.44880987,  5.32635859,  5.21315918,  5.10921165,
         5.01424388,  4.92852798,  4.85179185,  4.78430758,  4.72553097,
         4.67600623,  4.63573336,  4.60416814,  4.58158269,  4.5682491 ,
         4.56362317],
       [ 6.3794396 ,  6.56420497,  6.76067125,  6.96666152,  7.18190366,
         7.25047637,  7.08884068,  6.9350963 ,  6.78679419,  6.64366225,
         6.50406779,  6.36610601,  6.22678367,  6.08174693,  5.9252814 ,
         5.75004002,  5.54949203,  5.31901151,  5.05968691,  4.78049799,
         4.56362317],
       [ 6.43739987,  6.56420497,  6.76067125,  6.96666152,  7.18489691,
         7.40666978,  7.64041567,  7.88368555,  8.1362073 ,  8.39770881,
         8.66791797,  8.94737901,  8.77104916,  8.40369532,  8.0679067 ,
         7.77211883,  7.52476724,  7.33238267,  7.19686992,  7.11714053,
         7.09101759]],
    }


a = 4.23
atoms = Atoms('Na2', cell=(a, a, a), pbc=True,
              scaled_positions=[[0, 0, 0], [.5, .5, .5]])

# the calculator parameters illustrate only the way of calculating
# a band structure - they are far from converged!

nbands = 3
label = 'Na_sc'
# Make self-consistent calculation and save results
calc = Abinit(label=label,
              nbands=nbands,
              xc='LDA',
              ecut=20 * Ry,
              width=0.1,
              kpts=[4, 4, 4])

atoms.set_calculator(calc)
atoms.get_potential_energy()

# Extract Fermi level from the self-consistent calculation
e_fermi = calc.get_fermi_level()
assert nbands == calc.get_number_of_bands()

# parameters for calculation of band structure
# see http://www.abinit.org/Infos_v5.6/tutorial/lesson_3.html#35

# perform band structure calculation
label = 'Na_sc_band'
# Make self-consistent calculation and save results
calc = Abinit(label=label,
              nbands=nbands,
              xc='LDA',
              ecut=20 * Ry,
              width=0.1,
              kpts=[4, 4, 4])

calc.set_inp('ndtset', 2) # two datasets are used
calc.set_inp('iscf2', -2) # make a non-self-consistent calculation ;
calc.set_inp('getden2', -1) # to take the output density of dataset 1
calc.set_inp('kptopt2', -1) # to define one segment in the brillouin Zone
nband2 = 7
calc.set_inp('nband2', nband2) # use 7 bands in band structure calculation
calc.set_inp('ndivk2', 20) # with 21 divisions of the first segment
calc.set_inp('kptbounds2', "\n0.5  0.0  0.0\n0.0  0.0  0.0\n0.0  0.5  0.5\n1.0  1.0  1.0\n")
calc.set_inp('tolwfr2', 1.0e-12) #

atoms.set_calculator(calc)
atoms.get_potential_energy()

# interface returns always the values from the last dataset
assert nband2 == calc.get_number_of_bands()

# Calculate band structure along Gamma-X i.e. from 0 to 0.5

kpts2 = calc.get_ibz_k_points()
nkpts2 = len(kpts2)

eigs = np.empty((nband2, nkpts2), float)

for k in range(nkpts2):
    eigs[:, k] = calc.get_eigenvalues(kpt=k)

eigs -= e_fermi

results = {
    'e_fermi': e_fermi,
    'eigs': eigs,
    }

keys = ref.keys()
keys.sort()
for k in keys:
    if isinstance(ref[k], list):
        array_almost_equal(
            np.array(flatten(results[k])),
            np.array(flatten(ref[k])),
            tol=1e-4)
    else:
        diff = abs(results[k] - ref[k])
        assert diff < 1e-4, k + ': ' + str(diff)

try:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pylab

    for n in range(nband2):
        pylab.plot(kpts2[:, 0], eigs[n], '.m')
    plot_save(".", label)
except ImportError:
    raise NotAvailable('matplotlib required')
