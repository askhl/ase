#from __future__ import absolute_import  # PY24
import os
import operator

import numpy as np

if 1:
    def encode(obj):
        if isinstance(obj, str):
            return '"' + obj + '"'
        if isinstance(obj, (bool, np.bool_)):
            return repr(obj).lower()
        if isinstance(obj, (int, float)):
            return repr(obj)
        if isinstance(obj, dict):
            return '{' + ', '.join(['"' + key + '": ' + encode(value)
                                    for key, value in obj.items()]) + '}'
        if obj is None:
            return 'null'
        return '[' + ','.join([encode(value) for value in obj]) + ']'
    def loads(txt):
        return eval(txt, {'false': False, 'true': True, 'null': None})
else:
    from json import JSONEncoder, loads
    class NDArrayEncoder(JSONEncoder):
        def default(self, obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return JSONEncoder.default(self, obj)
    encode = NDArrayEncoder().encode
        
from ase.parallel import world
from ase.db import IdCollisionError
from ase.db.core import NoDatabase, dict2atoms


ops = {'<': operator.lt,
       '<=': operator.le,
       '=': operator.eq,
       '>=': operator.ge,
       '>': operator.gt}


def numpyfy(obj):
    if isinstance(obj, dict):
        return dict((key, numpyfy(value)) for key, value in obj.items())
    if isinstance(obj, list):
        try:
            a = np.array(obj)
        except ValueError:
            obj = [numpyfy(value) for value in obj]
        else:
            if a.dtype != object:
                obj = a
    return obj


def write_json(name, results):
    print results
    if world.rank == 0:
        fd = open(name, 'w')
        fd.write(encode(results))
        fd.close()


def read_json(name):
    fd = open(name, 'r')
    results = loads(fd.read())
    fd.close()
    world.barrier()
    return numpyfy(results)


class JSONDatabase(NoDatabase):
    def _write(self, id, atoms, keywords, key_value_pairs, data, replace):
        if os.path.isfile(self.filename):
            bigdct = read_json(self.filename)
            if not replace and id in bigdct:
                raise IdCollisionError
        else:
            bigdct = {}

        dct = self.collect_data(atoms)
        dct['id'] = id
        dct['keywords'] = keywords
        dct['key_value_pairs'] = key_value_pairs
        dct['data'] = data
        bigdct[id] = dct
        write_json(self.filename, bigdct)

    def get_dict(self, id):
        bigdct = read_json(self.filename)
        if id in [-1, 0]:
            assert len(bigdct) == 1
            id = bigdct.keys()[0]
        return bigdct[id]

    def _iselect(self, keywords, cmps):
        bigdct = read_json(self.filename)
        cmps = [(key, ops[op], val) for key, op, val in cmps]
        for id, dct in bigdct.items():
            for keyword in keywords:
                if keyword not in dct['keywords']:
                    break
            else:
                for key, op, val in cmps:
                    value = get_value(dct, key)
                    print key,value,op,val
                    if value is None or not op(value,val):
                        break
                else:
                    yield dct


def get_value(dct, key):
    value = dct['key_value_pairs'].get(key)
    if value is not None:
        return value
    if key in ['energy', 'magmom']:
        return dct.get('results', {}).get(key)
    if key in ['timestamp', 'username', 'calculator']:
        return dct.get(key)
    if isinstance(key, int):
        return (dct['numbers'] == key).sum()