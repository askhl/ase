from random import random
from ase.io import write
import time
from ase.optimize.genetic_algorithm.data import DataConnection
from ase.optimize.genetic_algorithm.population import Population
from ase.optimize.genetic_algorithm.standardcomparator import StandardComparator
from ase.optimize.genetic_algorithm.cutandsplicepairing import CutAndSplicePairing
from ase.optimize.genetic_algorithm.standardmutations import MutationSelector
from ase.optimize.genetic_algorithm.standardmutations import MirrorMutation
from ase.optimize.genetic_algorithm.standardmutations import RattleMutation
from ase.optimize.genetic_algorithm.standardmutations import PermutationMutation
from ase.optimize.genetic_algorithm.utilities import closest_distances_generator
from ase.optimize.genetic_algorithm.utilities import get_all_atom_types
from ase.optimize.genetic_algorithm.parallellocalrun import ParallelLocalRun

population_size = 20
mutation_probability = 0.3
n_to_test = 100


# Initialize the different components of the GA
da = DataConnection('gadb.db')
tmp_folder = 'tmp_folder/'

# An extra object is needed to handle the parallel execution
parallel_local_run = ParallelLocalRun(data_connection=da,
                                      tmp_folder=tmp_folder,
                                      n_simul=4,
                                      calc_script='calc.py')

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
    parallel_local_run.relax(a)

# Wait until the starting population is relaxed
while parallel_local_run.get_number_of_jobs_running() > 0:
    time.sleep(5.)

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
    parallel_local_run.relax(a3)
    population.update()

# Wait until the last candidates are relaxed
while parallel_local_run.get_number_of_jobs_running() > 0:
    time.sleep(5.)

write('all_candidates.traj', da.get_all_relaxed_candidates())
