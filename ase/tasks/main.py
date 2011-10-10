import os
import sys
import tempfile
import traceback

from ase.tasks.calcwrapper import calcnames, get_calculator_wrapper
from ase.tasks.molecule import MoleculeTask
from ase.tasks.bulk import BulkTask


usage = """\
Usage: ase [task] [calculator] [options] system(s)

task:       mol(ecule), bulk or the name of Python script that instantiates
            a Task object.
calculator: %s.
systems:    chemical formulas or filenames of files containing the atomic
            structure.

Try "ase mol --help" or "ase bulk --help".
""" % (', '.join(calcnames[:-1]) + ' or ' + calcnames[-1])


def run(args=sys.argv[1:], calcname='emt'):
    if isinstance(args, str):
        args = args.split(' ')

    argsoriginal = args[:]

    if (len(args) > 0 and
        (args[0] in ['mol', 'molecule', 'bulk'] or args[0].endswith('.py'))):
        taskname = args.pop(0)
    else:
        sys.stderr.write(usage)
        return

    if len(args) > 0:
        if args[0] in calcnames:
            calcname = args.pop(0)
    
    calcwrapper = get_calculator_wrapper(calcname)
    
    if taskname.endswith('.py'):
        locals = {}
        execfile(taskname, {}, locals)
        tasks = [task for task in locals.values()]
        assert len(tasks) == 1
        task = tasks[0]
    elif taskname == 'bulk':
        task = BulkTask(calcwrapper=calcwrapper)
    else:
        task = MoleculeTask(calcwrapper=calcwrapper)

    parser = task.create_parser()
    calcwrapper.add_options(parser)
    opts, args = task.parse(parser, args)
    calcwrapper.parse(opts)

    if opts.interactive_python_session:
        if '-i' in argsoriginal:
            argsoriginal.remove('-i')
        if '--interactive-python-session' in argsoriginal:
            argsoriginal.remove('--interactive-python-session')
        file = tempfile.NamedTemporaryFile()
        file.write('import os\n')
        file.write('if "PYTHONSTARTUP" in os.environ:\n')
        file.write('    execfile(os.environ["PYTHONSTARTUP"])\n')
        file.write('from ase.tasks.main import run\n')
        file.write('atoms, task = run(%r)\n' % argsoriginal)
        file.flush()
        os.system('python -i %s' % file.name)
        return

    atoms = task.run(args)

    return atoms, task


def main():
    try:
        run()
    except KeyboardInterrupt:
        raise
    except Exception:
        traceback.print_exc()
        sys.stderr.write("""
An exception occurred!  Please report the issue to
ase-developer@listserv.fysik.dtu.dk - thanks!  Please also report this
if it was a user error, so that a better error message can be provided
next time.""")
        raise
