import numpy as np

from ase.test import NotAvailable

try:
    import gpaw
except (ImportError, NotAvailable):
    raise NotAvailable('GPAW required')

from ase.structure import bulk

from ase.utils.compound_test import get_atomization_reactions, \
     calculate_reaction_energies, write_csv

from ase.utils.compound_test import BatchTest

from ase.utils.compound_test.periodic_system_gpaw import \
     GPAWEnergyPeriodicSystemTest

def atoms2data(atoms):
    data = {}
    for k in [
        'positions',
        'charges',
        'cell',
        ]:
        f = 'atoms.get_' + k + '()'
        data[k] = eval(f)
    data['magmoms'] = atoms.get_initial_magnetic_moments()
    data['symbols'] = atoms.get_chemical_symbols()
    return data

ref = {
    'Al': {
    'energy': -13.619147756299263,
    },
    }


data = {
    'Al': atoms2data(bulk('Al', 'fcc', 4.05, cubic=True)),
    }

def main():

    import os

    kwargs = dict(
        vacuum=0.0,
        xc='PBE',
        h=0.3,
        mode='lcao',
        basis='szp(dzp)',
        width=0.1,
        nbands=-5,
        kpts=[2, 2, 2],
        )

    dir = 'GPAW_Al_bulk'
    identifier = 'energy'

    etest = GPAWEnergyPeriodicSystemTest(dir + '/' + identifier,
                                         data=data, **kwargs)

    betest = BatchTest(etest)

    formulas = etest.get_formulas()

    betest.run(formulas)

    # analyse results

    results = etest.get_results(formulas, dir=dir)
    energies = {}
    for r in results:
        # compound: energy
        energies[r[0]] = r[1]

    db_data = {}
    for formula in formulas:
        try:
            db_data[formula] = {}
            import cmr
            filename = os.path.join(dir, identifier + '.' + formula)
            db = cmr.read(filename + '.db')
            for item in ['iter', 'energy']:
                db_data[formula][item] = db[item]
            # test
            for k in ref[formula].keys():
                diff = db_data[formula][k] - ref[formula][k]
                assert abs(diff) < 1.0e-2, k + ': ' + str(diff)
        except (ImportError, IOError, KeyError): # no cmr or no db data
            pass

        # test without cmr
        diff = energies[formula] - ref[formula]['energy']
        assert abs(diff) < 1.0e-2, formula + ': ' + str(diff)

if __name__ == '__main__':
    main()
