# creates: Al_phonon.png Al_mode.gif

from ase.structure import bulk
from ase.calculators.emt import EMT
from ase.dft.kpoints import ibz_points, get_bandpath
from ase.phonons import Phonons

# Setup crystal and EMT calculator
atoms = bulk('Al', 'fcc', a=4.05)
calc = EMT()

# Phonon calculator
N = 6
ph = Phonons(atoms, calc, supercell=(N, N, N))
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True)
# ph.clean()

# High-symmetry points in the Brillouin zone
points = ibz_points['fcc']
G = points['Gamma']
X = points['X']
W = points['W']
K = points['K']
L = points['L']
U = points['U']

point_names = ['$\Gamma$', 'X', 'U', 'L', '$\Gamma$', 'K']
path = [G, X, U, L, G, K]

path_kc, q, Q = get_bandpath(path, atoms.cell, 100)
omega_kn = 1000 * ph.band_structure(path_kc)

# Plot phonon dispersion
import pylab as plt

for n in range(len(omega_kn[0])):
    omega_n = omega_kn[:, n]
    plt.plot(q, omega_n, 'k-', lw=2)

plt.xticks(Q, point_names, fontsize=18)
plt.yticks(fontsize=18)
plt.xlim(q[0], q[-1])
plt.ylabel("Frequency ($\mathrm{meV}$)", fontsize=22)
plt.grid('on')
plt.savefig('Al_phonon.png')


# Write modes for specific q-vector to trajectory files
ph.write_modes(L, branches=[2], repeat=(5, 5, 5), kT=2e-4)

# Generate png animation
from subprocess import call
from ase.io import PickleTrajectory, write

trajfile = 'phonon.mode.2.traj'
trajectory = PickleTrajectory(trajfile, 'r')

for i, atoms in enumerate(trajectory):
    write('picture%02i.png' %i, atoms, show_unit_cell=2, rotation='-36x,26.5y,-25z')

call(['mplayer', '-vo', 'gif89a:fps=10.0:output=Al_mode.gif',
      '-mf', 'w=512:h=512:type=png:fps=10.0', 'mf://picture*.png'])



