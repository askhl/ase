import cPickle as cp
import os.path
import shutil
def save(obj, name='WHAM.pckl'):
    if os.path.isfile(name):
        shutil.move(name, name+'.bak')
    cp.dump((obj), open(name, 'wb'), 2)
    print 'object has been saved'	

def load(name='WHAM.pckl'):
    obj = cp.load(file(name))
    return obj 

