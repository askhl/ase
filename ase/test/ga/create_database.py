from ase.optimize.genetic_algorithm.data import PrepareDB
from ase.optimize.genetic_algorithm.data import DataConnection
import os
import sqlite3
import numpy as np

db_file = 'ga_db.sql'
db_folder = 'db_folder/'
tmp_folder = 'tmp_folder/'
import shutil
for o in [db_folder, tmp_folder]:
    shutil.rmtree(o, ignore_errors=True)
if os.path.isfile(db_file):
    os.remove(db_file)

d = PrepareDB(db_file_name=db_file,
              db_data_folder=db_folder,
              tmp_folder=tmp_folder)

assert os.path.isfile(db_file)
assert os.path.isdir(db_folder)
assert os.path.isdir(tmp_folder)

conn = sqlite3.connect(db_file)
cur = conn.cursor()

cur.execute('SELECT Name FROM sqlite_master WHERE type="table"')

tables = [t[0] for t in cur.fetchall()]
assert 'Configuration' in tables
assert 'History' in tables
assert 'MetaData' in tables
assert 'HistoryType' in tables

from ase.lattice.surface import fcc111

atom_numbers = np.array([78, 78, 79, 79])
slab = fcc111('Ag', size=(4,4,2), vacuum=10.)
d.add_slab(slab)
d.define_atom_numbers(atom_numbers)

conn.close()
d.close()

dc = DataConnection(db_file)

slab_get = dc.get_slab()
an_get = dc.get_atom_numbers_to_optimize()

assert len(slab) == len(slab_get)
assert np.all(slab.numbers == slab_get.numbers)
assert np.all(slab.get_positions() == slab_get.get_positions())
assert np.all(an_get == atom_numbers)
assert db_folder == dc.get_db_folder()
assert tmp_folder == dc.get_tmp_folder()


import shutil
for o in [db_folder, tmp_folder]:
    shutil.rmtree(o)
os.remove(db_file)
