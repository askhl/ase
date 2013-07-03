
from ase.embedding.multiase.mixer.selector import AtomListSelector
from ase.embedding.multiase.mixer.mixer import Mixer, EnergyCalculation 
from ase.embedding.multiase.mixer.mixer import ForceCalculation

class ONIOM(Mixer):
    def __init__(self, atoms, calc_full=None, cell_full=None,
                 calc_loc=None, cell_loc=None, atoms_loc=None,
                 selector='AtomListSelector'):
        """ The total energy and forces follow the linear combination: 
        H = H_calc_full(cell_full) + ( H_calc_loc(cell_loc) 
                                       - H_calc_full(cell_loc))
          = H_1(A+B) + ( H_2(B) - H_1(B) )
        Parameters in oniom
        atoms: ASE object with all atoms
        calc_full: ase-calculator for the full region
        cell_full: defines de size of the full region
        calc_loc: ase-calculator for the local region
        cell_loc: defines the size of the local region
        atoms_loc: tuple with atom indices belonging to the local region
        """

        Mixer.set_atom_ids(atoms) # sequence of numbers in same order as 
                               #positions were given above, index starts from 0

        #Create filter for the regions full (A+B) and local (B)
        if selector == 'AtomListSelector':
            atoms_full = tuple(range(len(atoms)))
            weights_full = {a:1.0 for a in atoms_full}
            weights_loc = {a:0.0 for a in atoms_full}
            for a in atoms_loc:
                weights_loc[a] = 1.0
            filter_AB = AtomListSelector(atoms_full, weights_full)        
            filter_B = AtomListSelector(atoms_full, weights_loc)
        else: # TODO use the box selector implemented in mixer
            raise ValueError("Invalid selector given")

        #Create energies and forces, default coefficient is 1.0
        #only H_1_B and F_1_B have negative coefficient
        forces_1_AB = ForceCalculation("forces_full",
                                       selector=filter_AB,
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
        
        energy_1_AB = EnergyCalculation("energy_full",
                                       selector=filter_AB,
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


        mixer_forces = (forces_1_AB, forces_2_B, forces_1_B)
        mixer_energies = (energy_1_AB, energy_2_B, energy_1_B)


        atoms.center()

        Mixer.__init__(self, name="oniom", forces=mixer_forces,
                       energies=mixer_energies)

        
