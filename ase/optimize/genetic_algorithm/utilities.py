""" Various utility methods used troughout the GA. """
from ase.data import covalent_radii
import itertools
import numpy as np
from ase.io import write, read
import sqlite3
import os
import time


def closest_distances_generator(atom_numbers, ratio_of_covalent_radii):
    """ Generates the blmin dict used across the GA.
        The distances are based on the covalent radii of the atoms.
    """
    cr = covalent_radii
    ratio = ratio_of_covalent_radii

    blmin = dict()
    for i in atom_numbers:
        blmin[(i, i)] = cr[i] * 2 * ratio
        for j in atom_numbers:
            if i == j:
                continue
            if (i, j) in blmin.keys():
                continue
            blmin[(i, j)] = blmin[(j, i)] = ratio * (cr[i] + cr[j])
    return blmin


def get_mic_distance(p1, p2, cell, pbc):
    """ This method calculates the shortest distance between p1 and p2
         through the cell boundaries defined by cell and pbc.
         This method works for reasonable unit cells, but not for extremely
         elongated ones.
    """
    ct = cell.T
    pos = np.mat((p1, p2))
    scaled = np.linalg.solve(ct, pos.T).T
    for i in xrange(3):
        if pbc[i]:
            scaled[:, i] %= 1.0
            scaled[:, i] %= 1.0
    P = np.dot(scaled, cell)

    pbc_directions = [[-1, 1] * direction + [0] for direction in pbc]
    translations = np.mat(list(itertools.product(*pbc_directions))).T
    p0r = np.tile(np.reshape(P[0, :], (3, 1)), (1, translations.shape[1]))
    p1r = np.tile(np.reshape(P[1, :], (3, 1)), (1, translations.shape[1]))
    dp_vec = p0r + ct * translations
    d = np.min(np.power(p1r - dp_vec, 2).sum(axis=0))**0.5
    return d


def db_call_with_error_tol(db_cursor, expression, args=[]):
    """ In case the GA is used on older versions of networking
         filesystems there might be some delays. For this reason
         some extra error tolerance when calling the SQLite db is
         employed.
    """
    i = 0
    while i < 10:
        try:
            db_cursor.execute(expression, args)
            return
        except sqlite3.OperationalError, e:
            print(e)
            time.sleep(2.)
        i += 1
    raise sqlite3.OperationalError(
        'Database still locked after 10 attempts (20 s)')


def save_trajectory(confid, trajectory, folder):
    """ Saves traj files to the database folder.
         This method should never be used directly,
         but only through the DataConnection object.
    """
    fname = os.path.join(folder, 'traj%05d.traj' % confid)
    write(fname, trajectory)
    return fname


def get_trajectory(fname):
    """ Extra error tolerance when loading traj files. """
    fname = str(fname)
    try:
        t = read(fname)
    except IOError, e:
        print('get_trajectory error ' + e)
    return t


def atoms_too_close(a, bl):
    """ Checks if any atoms in a are too close, as defined by
        the distances in the bl dictionary. """
    num = a.numbers
    for i in xrange(len(a)):
        for j in xrange(i + 1, len(a)):
            if a.get_distance(i, j, True) < bl[(num[i], num[j])]:
                return True
    return False


def atoms_too_close_two_sets(a, b, bl):
    """ Checks if any atoms in a are too close to an atom in b,
        as defined by the bl dictionary. """
    tot = a + b
    num = tot.numbers
    for i in xrange(len(a)):
        for j in xrange(len(a), len(tot)):
            if tot.get_distance(i, j, True) < bl[(num[i], num[j])]:
                return True
    return False


def get_all_atom_types(slab, atom_numbers_to_optimize):
    """ Utility method used to extract all unique atom types
        from the atoms object slab and the list of atomic numbers
        atom_numbers_to_optimize. """
    from_slab = list(set(slab.numbers))
    from_top = list(set(atom_numbers_to_optimize))
    from_slab.extend(from_top)
    return list(set(from_slab))
