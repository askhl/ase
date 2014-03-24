"""
    Objects which handle all communication with the SQLite database.
"""
import sqlite3
import os
from ase.optimize.genetic_algorithm.utilities import db_call_with_error_tol
from ase.optimize.genetic_algorithm.utilities import save_trajectory
from ase.optimize.genetic_algorithm.utilities import get_trajectory
import cPickle as pickle
import datetime
import ase.db

def split_description(desc):
    d = desc.split(':')
    assert len(d) == 2, desc
    return d[0], d[1]

class DataConnection(object):
    """ Class that handles all database communication.

        All data communication is collected in this class in order to
        make a decoupling of the data representation and the GA method.

        Parameters:

        db_file_name: Path to the SQLLite data file.

    """
    def __init__(self, db_file_name):
        self.db_file_name = db_file_name
        if not os.path.isfile(self.db_file_name):
            raise IOError('DB file {0} not found'.format(self.db_file_name))
        self.c = ase.db.connect(self.db_file_name)
        self.already_returned = []

    def get_number_of_unrelaxed_candidates(self):
        """ Returns the number of candidates not yet queued or relaxed. """
        return len(self.__get_ids_of_all_unrelaxed_candidates__())

    def get_an_unrelaxed_candidate(self):
        """ Returns a candidate ready for relaxation. """
        to_get = self.__get_ids_of_all_unrelaxed_candidates__()
        if len(to_get) == 0:
            raise ValueError('No unrelaxed candidate to return')
        
        a = self.c.get_atoms(gaid=to_get[0])
        a.info['confid']=to_get[0]
        return a

    def __get_ids_of_all_unrelaxed_candidates__(self):
        """ Helper method used by the two above methods. """

        all_unrelaxed_ids = [t.gaid for t in self.c.select(relaxed=0)]
        all_relaxed_ids = [t.gaid for t in self.c.select(relaxed=1)]
        all_queued_ids = [t.gaid for t in self.c.select(queued=1)]
        
        actually_unrelaxed = [gaid for gaid in all_unrelaxed_ids 
                              if (gaid not in all_relaxed_ids and
                                  gaid not in all_queued_ids)]

        return actually_unrelaxed

    def __get_latest_traj_for_confid__(self, confid):
        """ Method for obtaining the latest traj
            file for a given configuration.
            There can be several traj files for
            one configuration if it has undergone
            several changes (mutations, pairings, etc.)."""
        all = self.c.select(gaid=confid)
        all.sort(key = lambda x: x.mtime)
        return self.c.get_atoms(all[-1].gaid)
        
    def mark_as_queued(self, a):
        """ Marks a configuration as queued for relaxation. """
        gaid = a.info['confid']
        self.c.write(None, gaid=gaid, queued=1)

    def add_relaxed_step(self, a):
        """ After a candidate is relaxed it must be marked as such. """
        a.get_potential_energy()  # test that energy can be extracted
        gaid = a.info['confid']
        self.c.write(a, gaid=gaid, relaxed=1)
        
    def add_unrelaxed_candidate(self, candidate, description):
        """ Adds a new candidate which needs to be relaxed. """
        t, desc = split_description(description)
        kwargs = {'relaxed': 0,
                  t: 1,
                  'description': desc}

        gaid = self.c.write(candidate, **kwargs)
        self.c.update(gaid, gaid=gaid)
        candidate.info['confid'] = gaid

    def add_unrelaxed_step(self, candidate, description):
        """ Add a change to a candidate without it having been relaxed.
            This method is typically used when a
            candidate has been mutated. """
        gaid = candidate.info['confid']
        t, desc = split_description(description)
        kwargs = {'relaxed': 0,
                  t: 1,
                  'description': desc,
                  'gaid': gaid}
        self.c.write(candidate, **kwargs)

    def get_number_of_atoms_to_optimize(self):
        """ Get the number of atoms being optimized. """
        v = self.c.get(simulation_cell=True)
        return len(v.data.optimization_stoichiometry)

    def get_atom_numbers_to_optimize(self):
        """ Get the list of atom numbers being optimized. """
        v = self.c.get(simulation_cell=True)
        
        return v.data.optimization_stoichiometry

    def get_slab(self):
        """ Get the super cell, including stationary atoms, in which
            the structure is being optimized. """
        return self.c.get_atoms(simulation_cell=True)

    def get_participation_in_pairing(self):
        """ Get information about how many direct
            offsprings each candidate has, and which specific
            pairings have been made. This information is used
            for the extended fitness calculation described in
            L.B. Vilhelmsen et al., JACS, 2012, 134 (30), pp 12807-12816
        """
        entries = self.c.select(pairing=1)

        ntot = 0
        frequency = dict()
        pairs = []
        for e in entries:
            ntot += 1
            txt = e.description
            tsplit = txt.split(' ')
            c1 = int(tsplit[1])
            c2 = int(tsplit[2])
            pairs.append((min(c1, c2), max(c1, c2)))
            if c1 not in frequency.keys():
                frequency[c1] = 0
            frequency[c1] += 1
            if c2 not in frequency.keys():
                frequency[c2] = 0
            frequency[c2] += 1
        return (frequency, pairs)

    def get_all_relaxed_candidates(self, only_new=False):
        """ Returns all candidates that have been relaxed. The optional
            parameter only_new can be used to specify only to get
            candidates relaxed since last time this function was
            invoced. """

        entries = self.c.select(relaxed=1)
        
        trajs = []
        for v in entries:
            if only_new and v.gaid in self.already_returned:
                continue
            t = self.c.get_atoms(id=v.id)
            t.info['confid'] = v.gaid
            trajs.append(t)
            self.already_returned.append(v.gaid)
        trajs.sort(key=lambda x: x.get_potential_energy())
        return trajs

    def get_all_candidates_in_queue(self):
        """ Returns all structures that are queued, but have not yet
            been relaxed. """
        all_queued_ids = [t.gaid for t in self.c.select(queued=1)]
        all_relaxed_ids = [t.gaid for t in self.c.select(relaxed=1)]
        
        in_queue = [qid for qid in all_queued_ids if qid not in all_relaxed_ids]
        return in_queue

    def remove_from_queue(self, confid):
        """ Removes the candidate confid from the queue. """
        
        queued_ids = self.c.select(queued=1, gaid=confid)
        ids = [q.id for q in queued_ids]
        self.c.delete(ids)

class PrepareDB(object):

    """ Class used to initialize a database.

        This class is used once to setup the database and create
        working directories.

        Parameters:

        db_file_name: Database file to use
        db_data_folder: Directory the database can use for storage.
        tmp_folder: Temporary folder used during relaxations.

    """
    def __init__(self, db_file_name, simulation_cell, stoichiometry):
        if os.path.exists(db_file_name):
            raise IOError('DB file {0} already exists'.format(db_file_name))
        self.db_file_name = db_file_name

        self.c = ase.db.connect(self.db_file_name)

        self.c.write(simulation_cell, data={'optimization_stoichiometry': stoichiometry},
                     simulation_cell=True)

    def add_unrelaxed_candidate(self, candidate):
        """ Add an unrelaxed starting candidate. """
        gaid = self.c.write(candidate, origin='StartingCandidateUnrelaxed', relaxed=False)
        self.c.update(gaid, gaid=gaid)

    def add_relaxed_candidate(self, candidate):
        """ Add a relaxed starting candidate. """
        gaid = self.c.write(candidate, origin='StartingCandidateRelaxed', relaxed=True)
        self.c.update(gaid, gaid=gaid)
