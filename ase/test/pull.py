from ase import *

Cu = Atoms('Cu',
           positions=[(0, 0, 0)],
           pbc=(1, 0, 0),
           calculator=EMT())
traj = PickleTrajectory('Cu.traj', 'w')
for a in linspace(2.0, 4.0, 20):
    Cu.set_cell([a, 1, 1])
    traj.write(Cu)
