
from ase.embedding.multiase.mixer.selector import AtomListSelector
from ase.embedding.multiase.mixer.mixer import Mixer, EnergyCalculation, ForceCalculation
#from ase.embedding.multiase.utils import get_datafile

import numpy as np


class ONIOM(Mixer):
    def __init__(self, atoms, calc_full=None, cell_full=None,
                 calc_loc=None, cell_loc=None):
        """ Parameters in oniom
        atoms: ASE object with all atoms
        calc_full: ase-calculator for the full region
        cell_full: 
        calc_loc: ase-calculator for the embedded (localized) region
        cell_loc:
        """
        Mixer.set_atom_ids(atoms) # sequence of numbers in same order as positions
                          # were given above, index starts from 0


        calc_gpaw = calc_loc
        calc_reaxff = calc_full
        reaxff_cell = cell_full 
        gpaw_cell = cell_loc

        filter_full_system = AtomListSelector((0, 1),
                                              {0: 1.0,
                                               1: 1.0})
        
        filter_qm_region = AtomListSelector((0, 1), {0: 1.0, 1: 0.0})
        
        forces_full_system = ForceCalculation("forces_full_sys",
                                              selector=filter_full_system,
                                              calculator=calc_reaxff,
                                              cell=reaxff_cell)

        forces_qm_gpaw = ForceCalculation("forces_qm_gpaw",
                                          selector=filter_qm_region,
                                          calculator=calc_gpaw,
                                          cell=gpaw_cell)

        forces_qm_reaxff = ForceCalculation("forces_qm_reaxff",
                                            selector=filter_qm_region,
                                            calculator=calc_reaxff,
                                            cell=reaxff_cell)
        
        energy_full_system_reaxff = EnergyCalculation("energy_full_sys",
                                                      selector=filter_full_system,
                                                      calculator=calc_reaxff,
                                                      cell=reaxff_cell)

        energy_qm_region_reaxff = EnergyCalculation("energy_qm_reaxff",
                                                    selector=filter_qm_region,
                                                    calculator=calc_reaxff,
                                                    cell=reaxff_cell,
                                                    coeff=-1.0)

        energy_qm_region_gpaw = EnergyCalculation("energy_qm_gpaw",
                                                  selector=filter_qm_region,
                                                  calculator=calc_gpaw,
                                                  cell=gpaw_cell)

        mixer_forces = (forces_full_system, forces_qm_gpaw, forces_qm_reaxff)
        mixer_energies = (energy_full_system_reaxff,
                          energy_qm_region_reaxff,
                          energy_qm_region_gpaw)


        atoms.center()

        Mixer.__init__(self, name="H2_mixer", forces=mixer_forces,
                       energies=mixer_energies)

        
