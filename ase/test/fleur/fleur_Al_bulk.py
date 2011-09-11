from ase.test import NotAvailable

try:
    from ase.test.fleur_installed import fleur_installed
    fleur_installed()
except NotAvailable:
    raise NotAvailable('FLEUR required')

from ase.utils.compound_test import get_atomization_reactions, \
     calculate_reaction_energies, write_csv

from ase.utils.compound_test import BatchTest

from ase.utils.compound_test.periodic_system_fleur import \
     FLEUREnergyPeriodicSystemTest

from ase.test.gpaw.gpaw_Al_bulk import data

ref = {
    'Al': {
    'energy': -26429.277732777817,
    },
    }


def main():

    import os

    kwargs = dict(
        vacuum=0.0,
        xc='PBE',
        kpts=[2, 2, 2],
        )

    dir = 'FLEUR_Al_bulk'
    identifier = 'energy'

    etest = FLEUREnergyPeriodicSystemTest(dir + '/' + identifier,
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
            # currently ase.utils.compound_test.write_db
            # cannot be used for fleur due to different arguments
            # in calculator methods compared to other calculators (gpaw, ...)
            # See point 7. of https://trac.fysik.dtu.dk/projects/ase/ticket/27
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
