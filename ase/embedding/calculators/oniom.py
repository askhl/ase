
from ase.embedding.multiase.mixer.selector import AtomListSelector
from ase.embedding.multiase.mixer.mixer import Mixer, EnergyCalculation, ForceCalculation

import numpy as np


class ONIOM(Mixer):
    def __init__(self, atoms, calc_full=None, cell_full=None,
                 calc_loc=None, cell_loc=None, atoms_loc=None):
        """ The total energy and forces are a combination: 
        H = H_calc_full(cell_full) + ( H_calc_loc(cell_loc) 
                                       - H_calc_full(cell_loc))
          = H_1(A) + ( H_2(B) - H_1(B) )
        Parameters in oniom
        atoms: ASE object with all atoms
        calc_full: ase-calculator for the full region
        cell_full: defines de size of the full region
        calc_loc: ase-calculator for the localized region
        cell_loc: defines the size of the localized region
        atoms_loc: tuple with atom indices belonging to the localized region
        """

        Mixer.set_atom_ids(atoms) # sequence of numbers in same order as 
                               #positions were given above, index starts from 0

        #Create filter for the regions full (A) and localized (B)
        filter_A = AtomListSelector((0, 1), {0: 1.0, 1: 1.0})        
        filter_B = AtomListSelector((0, 1), {0:1.0, 1: 0.0})
        
        #Create energies and forces, default coefficient is 1.0
        #only H_1_B and F_1_B have negative coefficient
        forces_1_A = ForceCalculation("forces_full",
                                       selector=filter_A,
                                       calculator=calc_full,
                                       cell=cell_full)

        forces_2_B = ForceCalculation("forces_loc",
                                      selector=filter_B,
                                      calculator=calc_loc,
                                      cell=cell_loc)

        forces_1_B = ForceCalculation("forces_calc_full_cell_loc",
                                      selector=filter_B,
                                      calculator=calc_full,
                                      cell=cell_loc,
                                      coeff=-1.0)
        
        energy_1_A = EnergyCalculation("energy_full",
                                       selector=filter_A,
                                       calculator=calc_full,
                                       cell=cell_full)
        energy_2_B = EnergyCalculation("energy_loc",
                                       selector=filter_B,
                                       calculator=calc_loc,
                                       cell=cell_loc)

        energy_1_B = EnergyCalculation("energy_calc_full_cell_loc",
                                       selector=filter_B,
                                       calculator=calc_full,
                                       cell=cell_loc,
                                       coeff=-1.0)


        mixer_forces = (forces_1_A, forces_2_B, forces_1_B)
        mixer_energies = (energy_1_A, energy_2_B, energy_1_B)


        atoms.center()

        Mixer.__init__(self, name="oniom", forces=mixer_forces,
                       energies=mixer_energies)

        
