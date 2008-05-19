"""Structure optimization. """

import sys
import pickle
from math import sqrt
from os.path import isfile

from ase.parallel import rank

import numpy as npy


class Dynamics:
    """Base-class for all MD and structure optimization classes."""
    def __init__(self, atoms, logfile):
        self.atoms = atoms

        if rank != 0:
            logfile = None
        elif isinstance(logfile, str):
            if logfile == '-':
                logfile = sys.stdout
            else:
                logfile = open(logfile, 'a')
        self.logfile = logfile
        
        self.observers = []
        self.nsteps = 0

    def get_number_of_steps(self):
        return self.nsteps

    def attach(self, function, interval=1, *args, **kwargs):
        """Attach callback function.

        At every *interval* steps, call *function* with arguments
        *args* and keyword arguments *kwargs*."""

        self.observers.append((function, interval, args, kwargs))

    def call_observers(self):
        for function, interval, args, kwargs in self.observers:
            if self.nsteps % interval == 0:
                function(*args, **kwargs)


class Optimizer(Dynamics):
    """Base-class for all structure optimization classes."""
    def __init__(self, atoms, restart, logfile):
        Dynamics.__init__(self, atoms, logfile)
        self.restart = restart

        if restart is None or not isfile(restart):
            self.initialize()
        else:
            self.read()

    def run(self, fmax=0.05, steps=100000000):
        """Run structure optimization algorithm.

        This method will return when the forces on all individual
        atoms are less than *fmax* or when the number of steps exceeds
        *steps*."""

        self.fmax = fmax
        step = 0
        f = self.atoms.get_forces()
        while step < steps:
            self.log(f)
            if self.converged(f):
                return
            self.step(f)
            self.nsteps += 1
            step += 1
            f = self.atoms.get_forces()
            self.call_observers()

    def converged(self, forces=None):
        if forces is None:
            forces = self.atoms.get_forces()
        return (forces**2).sum(axis=1).max() < self.fmax**2

    def log(self, forces):
        if self.logfile is None:
            return
        fmax = sqrt((forces**2).sum(axis=1).max())
        e = self.atoms.get_potential_energy()
        name = self.__class__.__name__
        self.logfile.write('%s: %3d %15.6f %12.4f\n' %
                           (name, self.nsteps, e, fmax))
        self.logfile.flush()
        
    def dump(self, data):
        if rank == 0 and self.restart is not None:
            pickle.dump(data, open(self.restart, 'wb'), protocol=2)

    def load(self):
        return pickle.load(open(self.restart))
