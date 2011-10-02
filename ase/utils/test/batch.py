import os
import sys
import traceback

from ase.utils import opencew

from ase.parallel import barrier, rank

class BatchTest:
    """Contains logic for looping over tests and file management."""
    def __init__(self, test):
        self.test = test
        self.txt = sys.stdout # ?

    def run_single_test(self, formula):
        filename = self.test.get_lock_filename(formula)
        if rank == 0:
            print >> self.txt, self.test.name, formula, '...',
        barrier()
        self.txt.flush()
        if os.path.exists(filename):
            if rank == 0:
                print >> self.txt, 'Skipped.'
            return
        barrier()
        try:
            fd = opencew(filename) # Empty file
            if fd is not None:
                fd.close()
            system = self.test.setup(formula)
            self.test.run(formula, system)
            if rank == 0:
                print >> self.txt, 'OK!'
            self.txt.flush()
        except 'asdfg':
            if rank == 0:
                print >> self.txt, 'Failed!'
            traceback.print_exc(file=self.txt)
            print >> self.txt
            self.txt.flush()

    def run(self, formulas, fraction='1/1'):
        """Run a batch of tests.

        This will invoke the run_single_test method on each formula, printing
        status to stdout.

        The formulas that already have *.traj files are skipped.

        Keyword fraction allows one to split the run over several sets,
        for example:
        fraction='1/4' will appox. run the first one-fourth of the systems,
        fraction='4/4' will appox. run the last one-fourth of the systems
        (including the remaining systems if not divisible by 4),
        fraction=1/1 runs all the systems.

        """

        if rank == 0:
            # Create directories if necessary
            if self.test.dir and not os.path.isdir(self.test.dir):
                os.makedirs(self.test.dir)
        barrier()

        # split the formulas set into subsets
        nominator, denominator = fraction.split('/')
        nominator, denominator = int(nominator), int(denominator)
        assert nominator> 0
        assert nominator<= denominator
        reminder = len(formulas) % denominator
        quotient = int(len(formulas)/denominator)
        start = (nominator-1)*quotient
        if nominator == denominator:
            formulas_set = formulas[start:]
        else:
            stop = nominator*quotient
            formulas_set = formulas[start:stop]

        for formula in formulas_set:
            self.run_single_test(formula)

    def collect_results(self, formulas, verbose=False):
        """Yield results of calculations.  """
        for formula in formulas:
            try:
                results = self.test.retrieve_results(formula)
                if verbose:
                    if rank == 0:
                        print >> self.txt, 'Loaded:', formula
                yield formula, results
            except (IOError, RuntimeError, TypeError):
                # XXX which errors should we actually catch?
                if verbose:
                    if rank == 0:
                        print >> self.txt, 'Error:', formula, '[%s]' % filename
                    traceback.print_exc(file=self.txt)
