from ase import *
a = 3.6
b = a / 2
fcc = Atoms(symbols='Cu', positions=[(0, 0, 0)],
            cell=[(0, b, b), (b, 0, b), (b, b, 0)],
            pbc=1)
fcc *= (2, 1, 1)
fcc.set_calculator(EMT())
fcc.set_momenta([(0.9, 0.0, 0.0), (-0.9, 0, 0)])
md = VelocityVerlet(fcc)
def f():
    print fcc.get_potential_energy(), fcc.get_total_energy()
md.attach(f)
md.attach(PickleTrajectory('Cu2.traj', 'w', fcc).write, interval=3)
md.run(dt=0.1, steps=20)
fcc2 = PickleTrajectory('Cu2.traj', 'r')[-1]

