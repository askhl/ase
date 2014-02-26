from random import random
from ase.io import write
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
from ase.optimize.genetic_algorithm.pbs_queue_run import PBSQueueRun

def jtg(job_name, traj_file):
    s = '#!/bin/sh\n'
    s += '#PBS -l nodes=1:ppn=12\n'
    s += '#PBS -l walltime=48:00:00\n'
    s += '#PBS -N {0}\n'.format(job_name)
    s += '#PBS -q q12\n'
    s += 'cd $PBS_O_WORKDIR\n'
    s += 'python calc.py {0}\n'.format(traj_file)
    return s

population_size = 20
mutation_probability = 0.3

# Initialize the different components of the GA
da = DataConnection('ga_db.sql')

# The PBS queing interface is created
pbs_run = PBSQueueRun(da, 
                      job_prefix='Ag2Au2_opt',
                      n_simul=5,
                      job_template_generator=jtg)

tmp_folder = da.get_tmp_folder()
atom_numbers_to_optimize = da.get_atom_numbers_to_optimize()
n_to_optimize = len(atom_numbers_to_optimize)
slab = da.get_slab()
blmin = closest_distances_generator(get_all_atom_types(slab, atom_numbers_to_optimize), 
                                    ratio_of_covalent_radii=0.7)

comp = StandardComparator(n_top=n_to_optimize)
pairing = CutAndSplicePairing(slab, n_to_optimize, blmin)
mutations = MutationSelector([1., 1., 1.],
                             [MirrorMutation(blmin, n_to_optimize),
                              RattleMutation(blmin, n_to_optimize),
                              PermutationMutation(n_to_optimize)])

# Relax all unrelaxed structures (e.g. the starting population)
while da.get_number_of_unrelaxed_candidates() > 0 and not pbs_run.enough_jobs_running():
    a = da.get_an_unrelaxed_candidate()
    pbs_run.relax(a)

# create the population
population = Population(data_connection=da,
                        population_size=population_size,
                        comparator=comp)

# Submit new candidates until enough are running
while not pbs_run.enough_jobs_running() and len(population.get_current_population()) > 2:
    a1, a2 = population.get_two_candidates()
    a3, desc = pairing.pair_candidates(a1, a2)
    if a3 == None:
        continue
    da.add_unrelaxed_candidate(a3, description = desc)

    if random() < mutation_probability:
        a3_mut, desc = mutations.mutate(a3)
        if a3_mut != None:
            da.add_unrelaxed_step(a3_mut, desc)
            a3 = a3_mut
            
    pbs_run.relax(a3)

write('all_candidates.traj', da.get_all_relaxed_candidates())

da.close()
