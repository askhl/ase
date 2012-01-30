import os
from ase import Atom, Atoms
from ase.io import PickleTrajectory

co = Atoms([Atom('C', (0, 0, 0)),
            Atom('O', (0, 0, 1.2))])
traj = PickleTrajectory('1.traj', 'w', co)
for i in range(5):
    co.positions[:, 2] += 0.1
    traj.write()
del traj
traj = PickleTrajectory('1.traj', 'a')
co = traj[-1]
print co.positions
co.positions[:] += 1
traj.write(co)
del traj
t = PickleTrajectory('1.traj', 'a')
print t[-1].positions
print '.--------'
for a in t:
    print 1, a.positions[-1,2]
co.positions[:] += 1
t.write(co)
for a in t:
    print 2, a.positions[-1,2]
assert len(t) == 7

# Change of atomic numbers and pbc not allowed
co[0].number = 1
try:
    t.write(co)
except ValueError:
    pass
else:
    assert False

co[0].number = 6
co.pbc = True
try:
    t.write(co)
except ValueError:
    pass
else:
    assert False

# Change of number of atoms not allowed
co.pbc = False
o = co.pop(1)
try:
    t.write(co)
except ValueError:
    pass
else:
    assert False

co.append(o)
t.write(co)

# Append to a nonexisting file
fname = '2.traj'
if os.path.isfile(fname):
    os.remove(fname)
traj = PickleTrajectory(fname, 'a', co)
traj.write()
del traj
os.remove(fname)

# Check offsets with changing image byte-size
traj = PickleTrajectory('1.traj', 'a')
co = traj[0]
co.set_momenta([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
for i in range(20):
    co[1].x += 0.05
    co[1].y += 0.05
    co[1].z += 0.05
    traj.write(co)
del co.arrays['momenta']
for i in range(20):
    co[1].x += 0.05
    co[1].y += 0.05
    co[1].z += 0.05
    traj.write(co)
del traj
traj = PickleTrajectory('1.traj', 'a')
for a in traj:
    print a[0].z

# Check slicing of trajectory
traj1 = traj[10:]
traj2 = traj1[10:]
assert(traj[35] == traj1[25])
assert(traj[35] == traj2[15])

# Check the length
assert(len(traj) == 48)
assert(len(traj1) == 38)
assert(len(traj2) == 28)

