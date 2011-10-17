# -*- coding: utf-8 -*-
# test ASE3 eos vs ASE2' on EMT Al bulk

import numpy as np

from ase.structure import bulk

from ase.io.trajectory import PickleTrajectory

from ase.calculators.emt import EMT

# old ASE2 conversion factor
eVA3ToGPA = 160.21773

ref = {
    'volumes': [28.47776421093749, 30.785692069335923, 33.215062499999988, 35.768989415039037, 38.450586726562491],
    'energies': [0.035451988788178568, -0.0050928246886510209, -0.0030041812880963192, 0.041478738539716176, 0.12398550538752318],
    # name: (V0 A**3, E0 eV, B eV/A**3)
    # ASE2: ScientificPython 2.6.2/Numeric 24.2
    'Taylor': (31.840250671574669, -0.0093696499506199793, 37.32864897185673/eVA3ToGPA),
    'Murnaghan': (31.811886131656994, -0.0092309977511284558, 37.68507176723223/eVA3ToGPA),
    'Birch': (31.816150251684537, -0.0092578197283539661, 37.784548094764467/eVA3ToGPA),
    'BirchMurnaghan': (31.815036532251092, -0.0091856404643270734, 36.636139353078576/eVA3ToGPA),
    'PourierTarantola': (31.801306306064337, -0.0091924772234127095, 37.469831797604471/eVA3ToGPA),
    'Vinet': (31.812091857456515, -0.0092375570142525093, 37.679714432047632/eVA3ToGPA),
    'AntonSchmidt': (31.209372128872122, 0.030395802834144892, 29.565499969260323/eVA3ToGPA),
    # ASE3: scipy 0.7.0/numpy 1.3.0
    'poly': (31.830028344128394, -0.0093740142335718701, 0.23105595482095204),
    #
    'taylor': (31.846130632992622, -0.0094419960571400958, 0.23284801776858643),
    'murnaghan': (31.82352580694722, -0.0093556176994294381, 0.23032026740510439),
    'birch': (31.840266576530414, -0.0094115178983234758, 0.23247814615565363),
    'birchmurnaghan': (31.833413440234892, -0.0094036485972671252, 0.23163635940883145),
    'pouriertarantola': (31.82532984085244, -0.0093642551073028057, 0.2305735519625095),
    'vinet': (31.825762404500544, -0.0093694617768735004, 0.23068435951129596),
    'oldpoly': (31.846130632074892, -0.0094419960411316062, 0.23284801766725688),
    }

# original ASE2 methods

eos_strl = [
    'Taylor',
    'Murnaghan',
    'Birch',
    'BirchMurnaghan',
    'PourierTarantola',
    'Vinet',
    'AntonSchmidt',
    ]

# RuntimeError: Optimal parameters not found:
# Number of calls to function has reached maxfev = 1000.
eos_strl3 = [m for m in eos_strl]
eos_strl3.remove('AntonSchmidt')

results = {}

# prepare energies and volumes

b = bulk('Al', 'fcc', a=4.05, orthorhombic=True)
b.set_calculator(EMT())
cell = b.get_cell()

volumes = []
energies = []
traj = PickleTrajectory('eos.traj', 'w')
for x in np.linspace(0.95, 1.05, 5):
    b.set_cell(cell * x, scale_atoms=True)
    volumes.append(b.get_volume())
    energies.append(b.get_potential_energy())
    traj.write(b)

for n, (v, e) in enumerate(zip(volumes, energies)):
    vref = ref['volumes'][n]
    eref = ref['energies'][n]
    vabserr = abs((v - vref) / vref)
    vstrerr = str(n) + ' volume: ' + str(v) + ': ' + str(vref) + ': ' + str(vabserr)
    assert vabserr < 1.e-6, vstrerr
    eabserr = abs((e - eref) / eref)
    estrerr = str(n) + ' energy: ' + str(e) + ': ' + str(eref) + ': ' + str(eabserr)
    assert eabserr < 1.e-4, estrerr

# ASE2

try:
    from ASE.Utilities.EquationOfState import EquationOfState

    import sys

    from ase.utils import devnull

    for e in eos_strl:
        eos = EquationOfState(e, volumes, energies)
        sys.stdout = devnull
        print eos
        sys.stdout = sys.__stdout__
        results[e] = (eos.GetV0(), eos.GetEnergy(), eos.GetBulkModulus()/eVA3ToGPA)
except ImportError:
    pass

# ASE3

from ase.utils.eos import EquationOfState
from ase.utils.eosase2 import EquationOfState as EquationOfState2

for e in eos_strl3 + ['poly', 'oldpoly']:
    if e in eos_strl3 + ['oldpoly']:
        eos = EquationOfState2(volumes, energies, eos=e.lower())
    else:
        eos = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    results[e.lower()] = (v0, e0, B)


# test ASE2 vs ASE2 (if available)

for e in eos_strl:
    for n, v2 in enumerate(ref[e]):
        if n in [0, 2]: # only test volume and bulk modulus
            try:
                v3 = results[e][n]
                abserr = abs((v3 - v2) / v2)
                #print e, abserr
                strerr = e + ' 2 vs 2: ' + str(v3) + ': ' + str(v2) + ': ' + str(abserr)
                assert abserr < 1.e-6, strerr
            except KeyError:
                pass

# test ASE3 vs ASE2 reference

for e in eos_strl3:
    for n, v2 in enumerate(ref[e]):
        if n in [0, 2]: # only test volume and bulk modulus
            v3 = results[e.lower()][n]
            abserr = abs((v3 - v2) / v2)
            #print e, abserr
            strerr = e + ' 3 vs 2: ' + str(v3) + ': ' + str(v2) + ': ' + str(abserr)
            # ASE2/ScientificPython/Numeric vs ASE2 methods/scipy/numpy error ~ 2% for B!
            assert abserr < 3.e-2, strerr

# test ASE3 vs ASE3

for e in eos_strl3:
    for n, v2 in enumerate(ref[e.lower()]):
        if n in [0, 2]: # only test volume and bulk modulus
            v3 = results[e.lower()][n]
            abserr = abs((v3 - v2) / v2)
            #print e, abserr
            strerr = e + ' 3 vs 3: ' + str(v3) + ': ' + str(v2) + ': ' + str(abserr)
            assert abserr < 1.e-6, strerr

# test poly vs oldpoly

for e in ['poly']:
    for n, v2 in enumerate(ref['oldpoly']):
        if n in [0, 2]: # only test volume and bulk modulus
            v3 = results[e.lower()][n]
            abserr = abs((v3 - v2) / v2)
            #print e, abserr
            strerr = e + ' 3 vs 3: ' + str(v3) + ': ' + str(v2) + ': ' + str(abserr)
            assert abserr < 1.e-2, strerr
