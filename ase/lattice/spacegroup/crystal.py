# Copyright (C) 2010, Jesper Friis
# (see accompanying license files for details).

"""
A module for ASE for simple creation of crystalline structures from
knowledge of the space group.

"""

import numpy as np

import ase
from ase.atoms import string2symbols
from spacegroup import Spacegroup
from cell import cellpar_to_cell

__all__ = ['crystal']



def crystal(symbols=None, basis=None, spacegroup=1, setting=1, 
            cell=None, cellpar=None, 
            ab_normal=(0,0,1), a_direction=None, size=(1,1,1),
            ondublicates='warn', symprec=0.001, 
            pbc=True, **kwargs):
    """Create an Atoms instance for a conventional unit cell of a
    space group.

    Parameters:

    symbols : string | sequence of strings
        Either a string formula or a sequence of element
        symbols. E.g. ('Na', 'Cl') and 'NaCl' are equivalent.
    basis : list of scaled coordinates | atoms instance
        Positions of the non-equivalent sites given either as
        scaled positions or through an atoms instance.
    spacegroup : int | string | Spacegroup instance
        Space group given either as its number in International Tables
        or as its Hermann-Mauguin symbol.
    setting : 1 | 2
        Space group setting.
    cell : 3x3 matrix
        Unit cell vectors.
    cellpar : [a, b, c, alpha, beta, gamma]
        Cell parameters with angles in degree. Is not used when `cell`
        is given. 
    ab_normal : vector
        Is used to define the orientation of the unit cell relative
        to the Cartesian system when `cell` is not given. It is the
        normal vector of the plane spanned by a and b.
    a_direction : vector
        Defines the orientation of the unit cell a vector. a will be 
        parallel to the projection of `a_direction` onto the a-b plane.
    size : 3 positive integers
        How many times the conventional unit cell should be repeated
        in each direction.
    ondublicates : 'keep' | 'replace' | 'warn' | 'error'
        Action if `basis` contain symmetry-equivalent positions:
            'keep'    - ignore additional symmetry-equivalent positions
            'replace' - reolace
            'warn'    - like 'keep', but issue an UserWarning
            'error'   - raises a SpacegroupValueError
    symprec : float
        Minimum "distance" betweed two sites in scaled coordinates
        before they are counted as the same site.
    pbc : one or three bools
        Periodic boundary conditions flags.  Examples: True,
        False, 0, 1, (1, 1, 0), (True, False, False).  Default
        is True.

    Keyword arguments:

    All additional keyword arguments are passed on to the Atoms constructor. 
    Currently, probably the most useful additional keyword arguments are
    `constraint` and `calculator`.

    Examples:

    Two diamond unit cells (space group number 227)

    >>> diamond = crystal('C', [(0,0,0)], spacegroup=227, 
    ...     cellpar=[3.57, 3.57, 3.57, 90, 90, 90], size=(2,1,1))
    >>> ase.view(diamond)  # doctest: +SKIP

    A CoSb3 skutterudite unit cell containing 32 atoms

    >>> skutterudite = crystal(('Co', 'Sb'), 
    ...     basis=[(0.25,0.25,0.25), (0.0, 0.335, 0.158)], 
    ...     spacegroup=204, cellpar=[9.04, 9.04, 9.04, 90, 90, 90])
    >>> len(skutterudite)
    32
    """
    sg = Spacegroup(spacegroup, setting)
    if isinstance(symbols, ase.Atoms):
        basis = symbols
        symbols = basis.get_chemical_symbols()
    if isinstance(basis, ase.Atoms):
        basis_coords = basis.get_scaled_positions()
        if cell is None and cellpar is None:
            cell = basis.cell
        if symbols is None:
            symbols = basis.get_chemical_symbols()
    else:
        basis_coords = np.array(basis, dtype=float, copy=False, ndmin=2)
    sites, kinds = sg.equivalent_sites(basis_coords, 
                                       ondublicates=ondublicates, 
                                       symprec=symprec)
    symbols = parse_symbols(symbols)
    symbols = [symbols[i] for i in kinds]
    if cell is None:
        cell = cellpar_to_cell(cellpar, ab_normal, a_direction)

    atoms = ase.Atoms(symbols, 
                      scaled_positions=sites, 
                      cell=cell,
                      pbc=pbc,
                      **kwargs)

    if isinstance(basis, ase.Atoms):
        for name in basis.arrays:
            if not atoms.has(name):
                array = basis.get_array(name)
                atoms.new_array(name, [array[i] for i in kinds], 
                                dtype=array.dtype, shape=array.shape[1:])

    if size != (1, 1, 1):
        atoms = atoms.repeat(size)
    return atoms

def parse_symbols(symbols):
    """Return `sumbols` as a sequence of element symbols."""
    if isinstance(symbols, basestring):
        symbols = string2symbols(symbols)
    return symbols






#-----------------------------------------------------------------
# Self test
if __name__ == '__main__':
    import doctest
    print 'doctest: ', doctest.testmod()
    
