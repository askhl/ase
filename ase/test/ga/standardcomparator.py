from ase.optimize.genetic_algorithm.standardcomparator import StandardComparator
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

a1 = Atoms('AgAgAg', positions = [[0,0,0], [1.5,0,0], [1.5, 1.5, 0]])
a2 = Atoms('AgAgAg', positions = [[0,0,0], [1.4,0,0], [1.5, 1.5, 0]])

e1 = 1.0
e2 = 0.8

a1.set_calculator(SinglePointCalculator(e1, None, None, None, a1))
a2.set_calculator(SinglePointCalculator(e2, None, None, None, a2))

comp1 = StandardComparator(n_top=3,
                           pair_cor_cum_diff=0.03,
                           pair_cor_max=0.7,
                           dE=0.3)
assert comp1.looks_like(a1, a2)


comp2 = StandardComparator(n_top=3,
                           pair_cor_cum_diff=0.03,
                           pair_cor_max=0.7,
                           dE=0.15)
assert not comp2.looks_like(a1, a2)


comp3 = StandardComparator(n_top=3,
                           pair_cor_cum_diff=0.02,
                           pair_cor_max=0.7,
                           dE=0.3)
assert not comp3.looks_like(a1, a2)
