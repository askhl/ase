"""
    Creating of all the chemical elements for the periodic table.
"""
# Import ASE modules.
from ase.data import chemical_symbols, atomic_numbers, atomic_names, \
                     atomic_masses, covalent_radii
from ase.data.colors import jmol_colors                        
from ase.chemicalelement import ChemicalElement
# Import other modules.
import numpy as npy


# Create a list of all the chemical elements.
chemical_elements = [None]
for i in range(1,len(chemical_symbols)):
    # Fetch all the relevant data.
    name = atomic_names[i]
    symbol = chemical_symbols[i]
    atomic_number = atomic_numbers[symbol]
    # Create an variable with the correct name for easy importing.
    vars()[name] = ChemicalElement(None,
                                   symbol=symbol,
                                   name=name,
                                   atomic_number=atomic_number)
