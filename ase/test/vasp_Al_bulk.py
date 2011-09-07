from ase.test import NotAvailable

try:
    from ase.test.vasp_installed import vasp_installed
    vasp_installed()
except NotAvailable:
    raise NotAvailable('Vasp required')

from ase.utils.compound_test import get_atomization_reactions, \
     calculate_reaction_energies, write_csv

from ase.utils.compound_test import BatchTest

from ase.utils.compound_test.periodic_system_vasp import \
     VASPEnergyPeriodicSystemTest

from ase.test.gpaw_Al_bulk import data

ref = {
    'Al': {
    'energy': -14.158018999999999,
    },
    }


def main():

    import os

    kwargs = dict(
        vacuum=0.0,
        xc='PBE',
        algo='Fast', # may give wrong energies
        ediff=1.0e-3,
        prec='low',
        ismear=-1,
        sigma=0.1,
        # band parallelization: http://cms.mpi.univie.ac.at/vasp/guide/node138.html
        lplane=False,
        npar=1,
        kpts=[2, 2, 2],
        )

    dir = 'VASP_Al_bulk'
    identifier = 'energy'

    etest = VASPEnergyPeriodicSystemTest(dir + '/' + identifier,
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
