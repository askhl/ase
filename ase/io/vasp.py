"""
This module contains the class Atoms which extends ase Atoms
to handle VASP simulations. It will only work with ASE-3.0 and later
versions.
"""

import sys
import os
import re

def read_vasp(filename='CONTCAR'):
    """Import POSCAR/CONTCAR type file.

    Reads unitcell, atom positions and constraints from the POSCAR/CONTCAR
    file and tries to read atom types from POSCAR/CONTCAR header, if this fails
    the atom types are read from OUTCAR or POTCAR file.
    """
 
    import os
    from ase import Atoms, Atom
    from ase.constraints import FixAtoms, fix_scaled
    from ase.data import chemical_symbols
    import numpy as np

    if type(filename) == str:
        f = open(filename)
    elif type(filename) == file:
        f = filename
    else:
        raise TypeError("filename argument must be a string or a file object.")

    # First line should contain the atom symbols , eg. "Ag Ge" in
    # the same order
    # as later in the file (and POTCAR for the full vasp run)
    atomtypes = f.readline().split()
    try:
        for atype in atomtypes:
           if not atype in chemical_symbols:
              raise KeyError
    except KeyError:
        try:
            file_outcar = ReadOUTCAR(dir)
            atomtypes = file_outcar.atom_types()
        except IOError:
            file_potcar = ReadPOTCAR(dir)
            atomtypes = file_potcar.atom_types()

    lattice_constant = float(f.readline())

    # Now the lattice vectors
    a = []
    for ii in range(3):
        s = f.readline().split()
        floatvect = float(s[0]), float(s[1]), float(s[2])
        a.append(floatvect)

    basis_vectors = np.array(a) * lattice_constant

    # Number of atoms. Again this must be in the same order as
    # in the first line
    # or in the POTCAR or OUTCAR file
    atom_symbols = []
    numofatoms = f.readline().split()
    for i, num in enumerate(numofatoms):
        numofatoms[i] = int(num)
        if (len(atomtypes) < i + 1):
            atomtypes.append("Unknown")
        [atom_symbols.append(atomtypes[i]) for na in xrange(numofatoms[i])]

    # Check if Selective dynamics is switched on
    sdyn = f.readline()
    selective_dynamics = sdyn[0].lower() == "s"

    # Check if atom coordinates are cartesian or direct
    if selective_dynamics:
        ac_type = f.readline()
    else:
        ac_type = sdyn
    cartesian = ac_type[0].lower() == "c" or ac_type[0].lower() == "k"
    tot_natoms = sum(numofatoms)
    atoms_pos = np.empty((tot_natoms, 3))
    if selective_dynamics:
        selective_flags = np.empty((tot_natoms, 3), dtype=bool)
    for atom in xrange(tot_natoms):
        ac = f.readline().split()
        atoms_pos[atom] = (float(ac[0]), float(ac[1]), float(ac[2]))
        if selective_dynamics:
            curflag = []
            for flag in ac[3:6]:
                curflag.append(flag == 'F')
            selective_flags[atom] = curflag
    # Done with all reading
    if type(filename) == str:
        f.close()
    if cartesian:
        atoms_pos *= lattice_constant
    atoms = Atoms(symbols = atom_symbols, cell = basis_vectors, pbc = True)
    if cartesian:
        atoms.set_positions(atoms_pos)
    else:
        atoms.set_scaled_positions(atoms_pos)
    if selective_dynamics:
        constraints = []
        indices = []
        for ind, sflags in enumerate(selective_flags):
            if sflags.any() and not sflags.all():
                constraints.append(fix_scaled(atoms.get_cell(), ind, sflags))
            elif sflags.all():
                indices.append(ind)
        if indices:
            constraints.append(FixAtoms(indices))
        if constraints:
            atoms.set_constraint(constraints)
    return atoms

def write_vasp(filename, atoms, label='', direct=False, sort=None):
    """Method to write VASP position (POSCAR/CONTCAR) files.

    Writes label, scalefactor, unitcell, # of various kinds of atoms,
    positions in cartesian or scaled coordinates (Direct), and constraints
    to file. Cartesian coordiantes is default and default label is the 
    atomic species, e.g. 'C N H Cu'.
    """
    
    import numpy as np
    from ase.constraints import FixAtoms, fix_scaled

    if type(filename) == str:
        f = open(filename, 'w')
    elif type(filename) == file:
        f = file
    else:
        raise TypeError('filename argument must be either a string or a file \
object.')
    
    # Write atom positions in scaled or cartesian coordinates
    if direct:
        coord = atoms.get_scaled_positions()
    else:
        coord = atoms.get_positions()

    if atoms.constraints:
        sflags = np.zeros((len(atoms), 3), dtype=bool)
        for constr in atoms.constraints:
            if isinstance(constr, fix_scaled):
                sflags[constr.a] = constr.mask
            elif isinstance(constr, FixAtoms):
                sflags[constr.index] = [True, True, True]

    if sort:
        ind = np.argsort(atoms.get_chemical_symbols())
        symbols = np.array(atoms.get_chemical_symbols())[ind]
        coord = coord[ind]
        if atoms.constraints:
            sflags = sflags[ind]
    else:
        symbols = atoms.get_chemical_symbols()

    # Create a list of (symbol, count) pairs
    sc = []
    psym = symbols[0]
    count = 0
    for sym in symbols:
        if sym != psym:
            sc.append((psym, count))
            psym = sym
            count = 1
        else:
            count += 1
    sc.append((psym, count))

    # Create the label
    if label == '':
        for sym, c in sc:
            label += '%3s' % sym
    f.write(label + '\n')

    # Write unitcell in real coordinates and adapt to VASP convention 
    # for unit cell
    # ase Atoms doesn't store the lattice constant separately, so always
    # write 1.0.
    f.write('%19.16f\n' % 1.0)
    for vec in atoms.get_cell():
        f.write(' ')
        for el in vec:
            f.write('%22.16f' % el)
        f.write('\n')

    # Numbers of each atom
    for sym, count in sc:
        f.write('%4i' % count)
    f.write('\n')

    if atoms.constraints:
        f.write('Selective dynamics\n')

    if direct:
        f.write('Direct\n')
    else:
        f.write('Cartesian\n')

    for iatom, atom in enumerate(coord):
        for dcoord in atom:
            f.write('%20.16f' % dcoord)
        if atoms.constraints:
            for flag in sflags[iatom]:
                if flag:
                    s = 'F'
                else:
                    s = 'T'
                f.write('%4s' % s)
        f.write('\n')

    if type(filename) == str:
        f.close()

class ReadPOTCAR:
    """Class that read data from POTCAR file.

    Directory can be specified, default is current directory.
    """
    def __init__(self,dir='./'):
        self._file_ = os.path.join(dir, 'POTCAR')

    def atom_types(self):
        """Method that returns list of atomtypes."""
        file=open(self._file_,'r')
        lines=file.readlines()
        file.close()
        atomtypes=[]
        for line in lines:
            if re.search('TITEL',line):
                atomtypes.append(line.split()[3].split('_')[0].split('.')[0])
        return atomtypes

class ReadOUTCAR:
    """Class that read data from OUTCAR file. 

    Directory can be specified, default is current directory.
    """
    def __init__(self,dir='./'):
        
        self._file_ = os.path.join(dir, 'OUTCAR')

    def atom_types(self):
        """Method that returns list of atomtypes."""
        file=open(self._file_,'r')
        lines=file.readlines()
        file.close()
        atomtypes=[]
        for line in lines:
            if re.search('TITEL',line):
                atomtypes.append(line.split()[3].split('_')[0].split('.')[0])
        return atomtypes

