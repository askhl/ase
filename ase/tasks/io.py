import numpy as np

from ase.parallel import world


def dumps(obj):
    if isinstance(obj, str):
        return '"' + obj + '"'
    if isinstance(obj, (int, float)):
        return repr(obj)
    if isinstance(obj, dict):
        return '{' + ','.join(dumps(key) + ':' + dumps(value)
                               for key, value in obj.items()) + '}'
    return '[' + ','.join(dumps(value) for value in obj) + ']'
      

def loads(s):
    obj = eval(s)
    return npify(obj)


def npify(obj):
    if isinstance(obj, dict):
        return dict((key, npify(value)) for key, value in obj.items())
    if isinstance(obj, list):
        try:
            obj = np.array(obj)
        except ValueError:
            obj = [npify(value) for value in obj]
    return obj


class JSONWriter:
    def write(self, name, atoms, results):
        if world.rank == 0:
            fd = open(name + '.json', 'w')
            fd.write(dumps(results))
            fd.close()

class JSONReader:
    def read(self, name):
        fd = open(name + '.json', 'r')
        results = loads(fd.read())
        fd.close()
        return results
