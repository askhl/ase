from ase.optimize import BFGS
from ase.io import read, write
from ase.calculators.emt import EMT
from ase.optimize.genetic_algorithm.relax_attaches import VariansBreak
import sys


fname = sys.argv[1]

print 'Now relaxing {}'.format(fname)
a = read(fname)

a.set_calculator(EMT())
dyn = BFGS(a, trajectory=None, logfile=None)
vb = VariansBreak(a, dyn)
dyn.attach(vb.write)
dyn.run(fmax = 0.05)

write(fname[:-5] + '_relax.traj', a)

print 'Done relaxing {}'.format(fname)
