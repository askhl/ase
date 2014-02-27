from ase.optimize.genetic_algorithm.data import PrepareDB
from ase.optimize.genetic_algorithm.startgenerator import StartGenerator
from ase.optimize.genetic_algorithm.utilities import closest_distances_generator
from ase.io import read
from ase.visualize import view
from ase.constraints import FixAtoms
import numpy as np
from ase.lattice.surface import fcc111

# create the surface
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

# Define the composition of the atoms to optimize
atom_numbers = 2 * [47] + 2 * [79]

# define the closest distance two atoms of a given species can be to each other
cd = closest_distances_generator(atom_numbers = [47, 79],
                                 ratio_of_covalent_radii = 0.7)

# create the starting population
sg = StartGenerator(slab = slab,
                    atom_numbers = atom_numbers, 
                    closest_allowed_distances = cd,
                    box_to_place_in = [p0, [v1, v2, v3]])

# generate the starting population
starting_population = [sg.get_new_candidate() for i in xrange(20)]

# view(starting_population) # uncomment this line to see the starting population

# create the database to store information in
d = PrepareDB(db_file_name = 'ga_db.sql',
              db_data_folder = './db_folder/',
              tmp_folder = './tmp_folder/')

# add the information created above
d.add_slab(slab)
d.define_atom_numbers(atom_numbers)
for a in starting_population:
    d.add_unrelaxed_candidate(a)

# close the database
d.close()
