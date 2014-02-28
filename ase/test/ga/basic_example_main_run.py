from random import random
from ase.io import write
from ase.optimize import BFGS
from ase.calculators.emt import EMT

from ase.optimize.genetic_algorithm.data import DataConnection
from ase.optimize.genetic_algorithm.population import Population
from ase.optimize.genetic_algorithm.standardcomparator import StandardComparator
from ase.optimize.genetic_algorithm.cutandsplicepairing import CutAndSplicePairing
from ase.optimize.genetic_algorithm.utilities import closest_distances_generator
from ase.optimize.genetic_algorithm.utilities import get_all_atom_types
from ase.optimize.genetic_algorithm.standardmutations import MutationSelector
from ase.optimize.genetic_algorithm.standardmutations import MirrorMutation
from ase.optimize.genetic_algorithm.standardmutations import RattleMutation
from ase.optimize.genetic_algorithm.standardmutations import PermutationMutation

# Change the following three parameters to suit your needs
population_size = 20
mutation_probability = 0.3
n_to_test = 20

# Initialize the different components of the GA
da = DataConnection('ga_db.sql')
tmp_folder = da.get_tmp_folder()
atom_numbers_to_optimize = da.get_atom_numbers_to_optimize()
n_to_optimize = len(atom_numbers_to_optimize)
slab = da.get_slab()
blmin = closest_distances_generator(get_all_atom_types(slab, atom_numbers_to_optimize), 
                                    ratio_of_covalent_radii=0.7)

comp = StandardComparator(n_top=n_to_optimize,
                          pair_cor_cum_diff=0.015,
                          pair_cor_max=0.7,
                          dE=0.02,
                          mic=False)

pairing = CutAndSplicePairing(slab, n_to_optimize, blmin)
mutations = MutationSelector([1., 1., 1.],
                             [MirrorMutation(blmin, n_to_optimize),
                              RattleMutation(blmin, n_to_optimize),
                              PermutationMutation(n_to_optimize)])

# Relax all unrelaxed structures (e.g. the starting population)
while da.get_number_of_unrelaxed_candidates() > 0:
    a = da.get_an_unrelaxed_candidate()
    a.set_calculator(EMT())
    print 'Relaxing starting candidate {0}'.format(a.info['confid'])
    dyn = BFGS(a, trajectory=None, logfile=None)
    dyn.run(fmax=0.05, steps=100)
    da.add_relaxed_step(a)

# create the population
population = Population(data_connection=da,
                        population_size=population_size,
                        comparator=comp)

# test n_to_test new candidates
for i in xrange(n_to_test):
    print 'Now starting configuration number {0}'.format(i)
    a1, a2 = population.get_two_candidates()
    a3, desc = pairing.pair_candidates(a1, a2)
    if a3 == None:
        continue
    da.add_unrelaxed_candidate(a3, description=desc)

    # Check if we want to do a mutation
    if random() < mutation_probability:
        a3_mut, desc = mutations.mutate(a3)
        if a3_mut != None:
            da.add_unrelaxed_step(a3_mut, desc)
            a3 = a3_mut

    # Relax the new candidate
    a3.set_calculator(EMT())
    dyn = BFGS(a3, trajectory=None, logfile=None)
    dyn.run(fmax=0.05, steps=100)
    da.add_relaxed_step(a3)

write('all_candidates.traj', da.get_all_relaxed_candidates())

da.close()

