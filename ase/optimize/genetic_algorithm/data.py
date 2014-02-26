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

        self.conn = sqlite3.connect(db_file_name)
        self.cur = self.conn.cursor()

    def get_number_of_unrelaxed_candidates(self):
        """ Returns the number of candidates not yet queued or relaxed. """
        return len(self.__get_ids_of_all_unrelaxed_candidates__())

    def get_an_unrelaxed_candidate(self):
        """ Returns a candidate ready for relaxation. """
        to_get = self.__get_ids_of_all_unrelaxed_candidates__()
        return self.__get_latest_traj_for_confid__(to_get[0])

    def __get_ids_of_all_unrelaxed_candidates__(self):
        """ Helper method used by the two above methods. """
        sql = 'SELECT ConfID, HistType, TrajFile FROM History'
        db_call_with_error_tol(self.cur, sql)
        s = list(self.cur.fetchall())
        groups = dict()
        for l in s:
            if l[0] not in groups.keys():
                groups[l[0]] = []
            groups[l[0]].append(l)
        pending = []
        for k in groups.keys():
            relaxed = False
            for p in groups[k]:
                if p[1].find('Relaxed') != -1 or p[1].find('Queued') != -1:
                    relaxed = True
                    break
            if not relaxed:
                pending.append(k)
        return pending

    def __get_latest_traj_for_confid__(self, confid):
        """ Method for obtaining the latest traj
            file for a given configuration.
            There can be several traj files for
            one configuration if it has undergone
            several changes (mutations, pairings, etc.)."""
        sql = 'SELECT TrajFile FROM History WHERE ConfID = ?'
        sql += 'AND HistType != "Queued" ORDER BY HistTime DESC LIMIT 1'
        db_call_with_error_tol(self.cur, sql, [confid])
        s = list(self.cur.fetchall())
        ls = len(s)
        if ls != 1:
            raise ValueError('__get_latest_traj_for_confid__ ' +
                             'returned {0} possible traj files'.format(ls))
        t = get_trajectory(s[0][0])
        t.info['confid'] = confid
        return t

    def mark_as_queued(self, a):
        """ Marks a configuration as queued for relaxation. """
        confid = a.info['confid']
        sql = 'INSERT INTO History(ConfID, HistTime, ' + \
        'HistType, Description) VALUES(?, ?, ?, ?)'
        t = datetime.datetime.now()
        db_call_with_error_tol(self.cur, sql,
                               [confid, t, 'Queued', ''])
        self.conn.commit()

    def add_relaxed_step(self, a):
        """ After a candidate is relaxed it must be marked as such. """
        a.get_potential_energy()  # test that energy can be extracted
        confid = a.info['confid']

        sql = 'INSERT INTO History(ConfID,HistTime,HistType,Description) ' + \
          'VALUES(?, ?, ?, ?)'
        db_call_with_error_tol(self.cur, sql,
                               [confid, datetime.datetime.now(),
                                'Relaxed', ''])
        histid = self.cur.lastrowid
        folder = self.get_db_folder()
        fname = save_trajectory(histid, a, folder)
        sql = 'UPDATE History SET TrajFile = ? WHERE HistID = ?'
        db_call_with_error_tol(self.cur, sql, [fname, histid])
        self.conn.commit()

    def add_unrelaxed_candidate(self, candidate, description):
        """ Adds a new candidate which needs to be relaxed. """
        db_cursor = self.cur
        sql = 'INSERT INTO Configuration(CreationTime) VALUES(?)'
        db_call_with_error_tol(db_cursor, sql, [datetime.datetime.now()])
        confid = db_cursor.lastrowid
        candidate.info['confid'] = confid
        sql = 'INSERT INTO History(ConfID, ' + \
            'HistTime, HistType, Description) VALUES(?, ?, ?, ?)'
        db_call_with_error_tol(db_cursor, sql,
                               [confid, datetime.datetime.now(),
                                'NewCandNoEnergy', description])
        histid = db_cursor.lastrowid
        folder = self.get_db_folder()
        fname = save_trajectory(histid, candidate, folder)
        sql = 'UPDATE History SET TrajFile = ? WHERE HistID = ?'
        db_call_with_error_tol(db_cursor, sql, [fname, histid])
        self.conn.commit()

    def add_unrelaxed_step(self, candidate, description):
        """ Add a change to a candidate without it having been relaxed.
            This method is typically used when a
            candidate has been mutated. """
        confid = candidate.info['confid']
        sql = 'INSERT INTO History(ConfID, HistTime, HistType, ' + \
            'Description) VALUES(?, ?, ?, ?)'
        db_call_with_error_tol(self.cur, sql, [confid,
                                               datetime.datetime.now(),
                                               'StepNoEnergy', description])
        id2 = self.cur.lastrowid
        folder = self.get_db_folder()
        fname = save_trajectory(id2, candidate, folder)
        sql = 'UPDATE History SET TrajFile = ? WHERE HistID = ?'
        db_call_with_error_tol(self.cur, sql, [fname, id2])
        self.conn.commit()

    def get_tmp_folder(self):
        """ Get the temporary folder used for active relaxations. """
        if hasattr(self, 'tmp_folder'):
            return self.tmp_folder
        sql = 'SELECT Description FROM MetaData WHERE Name="tmp_folder"'
        db_call_with_error_tol(self.cur, sql)
        l = self.cur.fetchall()
        return str(l[0][0])

    def get_db_folder(self):
        """ Get the folder where the database stores all traj files. """
        if hasattr(self, 'db_data_folder'):
            return self.db_folder
        sql = 'SELECT Description FROM MetaData WHERE Name="db_data_folder"'
        db_call_with_error_tol(self.cur, sql)
        l = self.cur.fetchall()
        return str(l[0][0])

    def get_number_of_atoms_to_optimize(self):
        """ Get the number of atoms being optimized. """
        return len(self.get_atom_numbers_to_optimize())

    def get_atom_numbers_to_optimize(self):
        """ Get the list of atom numbers being optimized. """
        if hasattr(self, 'atom_numbers_to_optimize'):
            return len(self.atom_numbers_to_optimize)
        sql = 'SELECT Data FROM MetaData WHERE Name="atom_numbers"'
        db_call_with_error_tol(self.cur, sql)
        s = self.cur.fetchall()
        self.atom_numbers_to_optimize = pickle.loads(str(s[0][0]))
        return self.atom_numbers_to_optimize[:]

    def get_slab(self):
        """ Get the super cell, including stationary atoms, in which
            the structure is being optimized. """
        if hasattr(self, 'slab'):
            return self.slab.copy()
        sql = 'SELECT Data FROM MetaData WHERE Name="slab"'
        db_call_with_error_tol(self.cur, sql)
        s = self.cur.fetchall()
        self.slab = pickle.loads(str(s[0][0]))
        return self.slab.copy()

    def get_participation_in_pairing(self):
        """ Get information about how many direct
            offsprings each candidate has, and which specific
            pairings have been made. This information is used
            for the extended fitness calculation described in
            L.B. Vilhelmsen et al., JACS, 2012, 134 (30), pp 12807-12816
        """
        sql = 'SELECT ConfID, Description FROM ' + \
            'History WHERE Description LIKE ?'
        db_call_with_error_tol(self.cur, sql, ['pairing:%'])

        entries = self.cur.fetchall()

        frequency = dict()
        pairs = []
        for e in entries:
            txt = str(e[1])
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

    def get_all_relaxed_candidates(self, since=None):
        """ Returns all candidates that have been relaxed.
            The optional parameter since can be used to specify only
            to get candidates relaxed since this time. """
        if since == None:
            sql = 'SELECT ConfID, TrajFile FROM History WHERE ' + \
                'HistType = "Relaxed"'
            db_call_with_error_tol(self.cur, sql)
        else:
            sql = 'SELECT ConfID, TrajFile FROM History WHERE ' + \
                'HistType = "Relaxed" AND HistTime > ?'
            db_call_with_error_tol(self.cur, sql, [since])
        s = list(self.cur.fetchall())
        trajs = []
        for (confid, trajfile) in s:
            t = get_trajectory(trajfile)
            t.info['confid'] = confid
            trajs.append(t)
        trajs.sort(key=lambda x: x.get_potential_energy())
        return trajs

    def get_all_candidates_in_queue(self):
        """ Returns all structures that are queued, but have not yet
            been relaxed. """
        sql = ''' SELECT DISTINCT ConfID
                  FROM
                    History
                  WHERE
                    HistType = "Queued"
                   AND ConfID NOT IN
                    (SELECT ConfID FROM History WHERE
                         HistType LIKE "%Relaxed")'''
        db_call_with_error_tol(self.cur, sql)
        s = [l[0] for l in self.cur.fetchall()]
        return s

    def remove_from_queue(self, confid):
        """ Removes the candidate confid from the queue. """
        sql = ''' DELETE FROM History
                   WHERE ConfID = ? AND HistType="Queued" '''
        db_call_with_error_tol(self.cur, sql, [confid])
        self.conn.commit()

    def close(self):
        """ Close the database. """
        self.conn.commit()
        self.conn.close()


class PrepareDB(object):

    """ Class used to initialize a database.

        This class is used once to setup the database and create
        working directories.

        Parameters:

        db_file_name: Database file to use
        db_data_folder: Directory the database can use for storage.
        tmp_folder: Temporary folder used during relaxations.

    """
    def __init__(self, db_file_name, db_data_folder, tmp_folder):
        if os.path.exists(db_file_name):
            raise IOError('DB file {0} already exists'.format(db_file_name))
        self.db_file_name = db_file_name

        if os.path.exists(db_data_folder):
            txt = 'DB data folder {0} already exists'.format(db_data_folder)
            raise IOError(txt)
        self.db_data_folder = db_data_folder

        if os.path.exists(tmp_folder):
            txt = 'tmp folder {0} already exists'.format(db_data_folder)
            raise IOError(txt)
        self.tmp_folder = tmp_folder

        os.mkdir(db_data_folder)
        os.mkdir(tmp_folder)

        self.conn = sqlite3.connect(db_file_name)
        self.cur = self.conn.cursor()
        self.__create_tables__(self.cur)
        sql = 'INSERT INTO MetaData(Name, Description) VALUES(?, ?)'
        db_call_with_error_tol(self.cur, sql, ['tmp_folder', tmp_folder])
        db_call_with_error_tol(self.cur, sql,
                               ['db_data_folder', db_data_folder])
        db_call_with_error_tol(self.cur, sql, ['db_file', db_file_name])

        self.conn.commit()

    def __create_tables__(self, cur):
        """ Private method used to create all tables. """
        c = cur
        c.execute('''CREATE TABLE MetaData(
                              Name TEXT PRIMARY KEY,
                              Description TEXT,
                              Data LONGBLOG)''')

        c.execute('''CREATE TABLE Configuration(
                              ConfID INTEGER PRIMARY KEY AUTOINCREMENT,
                              CreationTime timestamp)''')

        c.execute('''CREATE TABLE History(
                              HistID INTEGER PRIMARY KEY AUTOINCREMENT,
                              ConfID INTEGER,
                              HistTime timestamp,
                              HistType TEXT,
                              Description TEXT,
                              TrajFile TEXT)''')

        c.execute('''CREATE TABLE HistoryType(
                              HistTypeID INTEGER PRIMARY KEY AUTOINCREMENT,
                              Title TEXT)''')

        c.executemany('INSERT INTO HistoryType(Title) values(?)',
                      [('StartCandNoEnergy', ),
                       ('StartCandRelaxed', ),
                       ('NewCandNoEnergy', ),
                       ('NewCandRelaxed', ),
                       ('StepNoEnergy', ),
                       ('Queued', ),
                       ('Relaxed', )])

    def add_slab(self, slab):
        """ Add the super cell, including all
            stationary atoms, to the database. """
        sql = 'SELECT Description FROM MetaData WHERE Name="slab"'
        db_call_with_error_tol(self.cur, sql)
        slabs = self.cur.fetchall()
        if len(slabs) > 0:
            raise ValueError('There is already a slab in the database')
        sql = '''INSERT INTO MetaData(Name, Data) VALUES("slab", ?)'''
        db_call_with_error_tol(self.cur, sql, [pickle.dumps(slab)])

    def define_atom_numbers(self, atom_numbers):
        """ Set the composition of the atoms to optimize. """
        sql = 'SELECT Description FROM MetaData WHERE Name="atom_numbers"'
        db_call_with_error_tol(self.cur, sql)
        atom_numbers_db = self.cur.fetchall()
        if len(atom_numbers_db) > 0:
            err_text = 'There is already atomic numbers in the database'
            raise ValueError(err_text)
        sql = 'INSERT INTO MetaData(Name, Data) VALUES("atom_numbers", ?)'
        db_call_with_error_tol(self.cur, sql, [pickle.dumps(atom_numbers)])

    def add_unrelaxed_candidate(self, candidate):
        """ Add an unrelaxed starting candidate. """
        self.add_candidate(candidate, 'StartCandNoEnergy')

    def add_relaxed_candidate(self, candidate):
        """ Add a relaxed starting candidate. """
        self.add_candidate(candidate, 'StartCandRelaxed')

    def add_candidate(self, candidate, cand_type):
        """ Helper method used by the two above methods. """
        db_cursor = self.cur
        sql = 'INSERT INTO Configuration(CreationTime) VALUES(?)'
        db_call_with_error_tol(db_cursor, sql, [datetime.datetime.now()])
        confid = db_cursor.lastrowid
        sql = 'INSERT INTO History(ConfID, HistTime, ' + \
            'HistType, Description) VALUES(?, ?, ?, ?)'
        db_call_with_error_tol(db_cursor,
                               sql,
                               [confid,
                                datetime.datetime.now(), cand_type, ''])
        histid = db_cursor.lastrowid
        folder = self.db_data_folder
        candidate.info['confid'] = confid
        fname = save_trajectory(histid, candidate, folder)
        sql = 'UPDATE History SET TrajFile = ? WHERE HistID = ?'
        db_call_with_error_tol(db_cursor, sql, [fname, histid])

    def close(self):
        """ Commit changes to the database and close the connection. """
        self.conn.commit()
        self.conn.close()
