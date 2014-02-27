from ase.optimize.genetic_algorithm.data import PrepareDB
from ase.optimize.genetic_algorithm.data import DataConnection
from ase.optimize.genetic_algorithm.startgenerator import StartGenerator
from ase.optimize.genetic_algorithm.utilities import closest_distances_generator
import os
import numpy as np
from ase.lattice.surface import fcc111
from ase.constraints import FixAtoms
from ase.calculators.singlepoint import SinglePointCalculator

db_file = 'ga_db.sql'
db_folder = 'db_folder/'
tmp_folder = 'tmp_folder/'

d = PrepareDB(db_file_name=db_file,
              db_data_folder=db_folder,
              tmp_folder=tmp_folder)

slab = fcc111('Au', size=(4,4,2), vacuum=10.0, orthogonal = True)
slab.set_constraint(FixAtoms(mask = slab.positions[:,2] <=10.))

# define the volume in which the adsorbed cluster is optimized
# the volume is defined by a corner position (p0)
# and three spanning vectors (v1, v2, v3)
pos = slab.get_positions()
cell = slab.get_cell()
p0 = np.array([0., 0., max(pos[:, 2]) + 2.])
v1 = cell[0, :] * 0.8
v2 = cell[1, :] * 0.8
v3 = cell[2, :]
v3[2] = 3.

# define the closest distance two atoms of a given species can be to each other
cd = closest_distances_generator(atom_numbers = [47, 79],
                                 ratio_of_covalent_radii = 0.7)

# Define the composition of the atoms to optimize
atom_numbers = 2 * [47] + 2 * [79]

# create the starting population
sg = StartGenerator(slab = slab,
                    atom_numbers = atom_numbers,
                    closest_allowed_distances = cd,
                    box_to_place_in = [p0, [v1, v2, v3]])

# generate the starting population
starting_population = [sg.get_new_candidate() for i in xrange(20)]

d.add_slab(slab)
d.define_atom_numbers(atom_numbers)

for a in starting_population:
    d.add_unrelaxed_candidate(a)

d.close()

# and now for the actual test
dc = DataConnection(db_file)

slab_get = dc.get_slab()
an_get = dc.get_atom_numbers_to_optimize()

assert dc.get_number_of_unrelaxed_candidates() == 20

a1 = dc.get_an_unrelaxed_candidate()
dc.mark_as_queued(a1)

assert dc.get_number_of_unrelaxed_candidates() == 19
assert len(dc.get_all_candidates_in_queue()) == 1

a1.set_calculator(SinglePointCalculator(0.0, None, None, None, a1))
dc.add_relaxed_step(a1)

assert dc.get_number_of_unrelaxed_candidates() == 19
assert len(dc.get_all_candidates_in_queue()) == 0

assert len(dc.get_all_relaxed_candidates()) == 1

a2 = dc.get_an_unrelaxed_candidate()
dc.mark_as_queued(a2)
confid = a2.info['confid']
assert dc.get_all_candidates_in_queue()[0] == confid

dc.remove_from_queue(confid)
assert len(dc.get_all_candidates_in_queue()) == 0

import shutil
for o in [db_folder, tmp_folder]:
    shutil.rmtree(o)
os.remove(db_file)
