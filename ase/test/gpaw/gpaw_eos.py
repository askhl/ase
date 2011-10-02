import os

import time

from ase.test import NotAvailable

try:
    import gpaw
except (ImportError, NotAvailable):
    raise NotAvailable('GPAW required')

from ase.structure import bulk

from ase.units import kcal, mol

from ase.utils.test.batch import BatchTest

from ase.utils.test.bulk import BulkTestEOS as Test

from gpaw import GPAW

kwargs = dict(
    xc='LDA',
    h=0.3,
    mode='lcao',
    basis='szp(dzp)',
    width=0.1,
    nbands=-5,
    kpts=[2, 2, 2],
    )

def gpaw(formula, system, test):
    return GPAW(txt=test.get_filename(formula, extension='txt'),
                charge=sum(system.get_charges()),
                **kwargs)

data = {
    'Al': bulk('Al', 'fcc', 4.05, cubic=True),
    }

dir = 'bulk_GPAW'

identifier = 'eos'

etest = Test(gpaw, name=dir + '/' + identifier, data=data)

betest = BatchTest(etest)

formulas = etest.get_formulas()

betest.run(formulas, fraction='1/1')

