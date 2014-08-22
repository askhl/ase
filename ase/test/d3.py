from ase.lattice.surface import fcc111
from ase.data.s22 import create_s22_system
from ase.calculators.d3 import D3

fortran = True

try:
    from ase.calculators.d3.d3ef import d3ef
except ImportError:
    fortran = False

def testforce(system, fortran, bj, threebody):
    system.set_calculator(D3(fortran=fortran, bj=bj, threebody=threebody))
    f1 = system.get_forces()
    f2 = system.calc.calculate_numerical_forces(system, 0.0001)
    
    assert abs(f1 - f2).max() < 1e-6

waterdim = create_s22_system('Water_dimer')
slab = fcc111('Au', (1, 1, 4), vacuum=7.5)

testforce(waterdim, fortran=False, bj=False, threebody=False)
testforce(waterdim, fortran=False, bj=True, threebody=False)
testforce(waterdim, fortran=False, bj=False, threebody=True)
testforce(waterdim, fortran=False, bj=True, threebody=True)

testforce(slab, fortran=False, bj=False, threebody=False)
testforce(slab, fortran=False, bj=True, threebody=False)
#testforce(slab, fortran=False, bj=False, threebody=True) # This is too slow!
#testforce(slab, fortran=False, bj=True, threebody=True) # This is too slow!

if fortran:
    testforce(waterdim, fortran=True, bj=False, threebody=False)
    testforce(waterdim, fortran=True, bj=True, threebody=False)
    testforce(waterdim, fortran=True, bj=False, threebody=True)
    testforce(waterdim, fortran=True, bj=True, threebody=True)

    testforce(slab, fortran=True, bj=False, threebody=False)
    testforce(slab, fortran=True, bj=True, threebody=False)
    testforce(slab, fortran=True, bj=False, threebody=True)
    testforce(slab, fortran=True, bj=True, threebody=True)
