from ase.cli.run import RunCommand
from ase.constraints import FixAtoms, UnitCellFilter
from ase.optimize import LBFGS
from ase.io.trajectory import PickleTrajectory


class OptimizeCommand(RunCommand):
    def add_parser(self, subparser):
        parser = subparser.add_parser('optimize', help='relax ...')
        self.add_arguments(parser)
        
    def add_arguments(self, parser):
        RunCommand.add_arguments(self, parser)
        add = parser.add_argument
        add('-f', '--maximum-force', default=0.05, type=float,
            help='Relax internal coordinates.')
        add('--constrain-tags',
            metavar='T1,T2,...',
            help='Constrain atoms with tags T1, T2, ...')
        add('-s', '--maximum-stress', type=float,
            help='Relax unit-cell and internal coordinates.')

    def calculate(self, atoms, name):
        args = self.args
        if args.constrain_tags:
            tags = [int(t) for t in args.constrain_tags.split(',')]
            mask = [t in tags for t in atoms.get_tags()]
            atoms.constraints = FixAtoms(mask=mask)
        
        trajectory = PickleTrajectory(self.get_filename(name, 'traj'), 'w',
                                      atoms)
        if args.maximum_stress:
            optimizer = LBFGS(UnitCellFilter(atoms), logfile=self.logfile)
            fmax = args.maximum_stress
        else:
            optimizer = LBFGS(atoms, logfile=self.logfile)
            fmax = args.maximum_force

        optimizer.attach(trajectory)
        optimizer.run(fmax=fmax)

        data = RunCommand.calculate(self, atoms, name)

        if hasattr(optimizer, 'force_calls'):
            data['force calls'] = optimizer.force_calls

        return data
