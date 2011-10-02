import time

import numpy as np

from ase.io.trajectory import PickleTrajectory

from ase.utils.eos import EquationOfState

from ase.utils.test.compound import CompoundTestEnergy

from ase.utils.test.compound import write_db

class BulkTestEOS(CompoundTestEnergy):
    """Perform equation of state. """

    def __init__(self, calculator, calculate=None, retrieve=None,
                 name='test', data=None):

        # provide default methods
        if calculate is None:
            calculate = calculate_eos1

        if retrieve is None:
            retrieve = retrieve_bulk_properties

        CompoundTestEnergy.__init__(self, calculator, calculate, retrieve,
                                    name=name, data=data)

def calculate_eos1(formula, system, test):
    # standard equation of state (scale all cell dimension in the same way)
    t = time.time()
    nstep, npoints = 0.02, 5
    assert npoints >= 3
    cell = system.get_cell()
    assert not (np.array(cell - np.transpose(cell)).any()), 'general cells unsuported'
    if cell.ndim == 2: # find diagonal cell parameters
        a0, b0, c0 = cell.diagonal()
    else:
        a0, b0, c0 = cell
    a, b, c = a0, b0, c0
    # prefer right side
    strains = nstep * np.array(range(- npoints/2 + 1, npoints/2 + 1))
    filename = test.get_lock_filename(formula)
    traj = PickleTrajectory(open(filename, 'w'), 'w')
    #
    calculator = system.get_calculator()
    energies = []
    volumes = []
    for ns, strain in enumerate(strains):
        a = a0 + strain
        b = b0 * a / a0
        c = c0 * a / a0
        system_abc = system.copy()
        system_abc.set_cell([a, b, c], scale_atoms=True)
        system_abc.set_calculator(calculator)
        energies.append(system_abc.get_potential_energy())
        volumes.append(system_abc.get_volume())
        traj.write(system_abc)
    traj.close()
    eos = EquationOfState(volumes, energies)
    v, e, B = eos.fit()
    x = (v / system.get_volume())**(1.0 / 3)
    # found cell
    a = a0 * x
    b = b0 * x
    c = c0 * x
    t = time.time() - t
    kwargs = {'t': t}
    kwargs.update({'iter': npoints})
    kwargs.update({
        'B': B,
        'v': v,
        'a': a,
        'c/a': c/a,
        'eos_energy': e,
        'eos_nstep': nstep,
        'eos_npoints': npoints,
        })
    write_db(test.get_filename(formula, extension='db'), system, **kwargs)

def retrieve_bulk_properties(formula, test):
    # this is wrong!
    system = ase.io.read(test.get_lock_filename(formula))
    return system.get_potential_energy()
