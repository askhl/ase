"""
This module contains functionality for reading and writing an ASE
Atoms object in VASP POSCAR format.

"""

import os

def get_atomtypes(fname):
    """Given a file name, get the atomic symbols. 

    The function can get this information from OUTCAR and POTCAR
    format files.  The files can also be compressed with gzip or
    bzip2.

    """
    atomtypes=[]
    if fname.find('.gz') != -1:
        import gzip
        f = gzip.open(fname)
    elif fname.find('.bz2') != -1:
        import bz2
        f = bz2.BZ2File(fname)
    else:
        f = open(fname)
    for line in f:
        if line.find('TITEL') != -1:
            atomtypes.append(line.split()[3].split('_')[0].split('.')[0])
    return atomtypes

def atomtypes_outpot(posfname, numsyms):
    """Try to retreive chemical symbols from OUTCAR or POTCAR
    
    If getting atomtypes from the first line in POSCAR/CONTCAR fails, it might
    be possible to find the data in OUTCAR or POTCAR, if these files exist.

    posfname -- The filename of the POSCAR/CONTCAR file we're trying to read
    
    numsyms -- The number of symbols we must find

    """
    import os.path as op
    import glob

    # First check files with exactly same name except POTCAR/OUTCAR instead
    # of POSCAR/CONTCAR.
    fnames = [posfname.replace('POSCAR', 'POTCAR').replace('CONTCAR', 
                                                           'POTCAR')]
    fnames.append(posfname.replace('POSCAR', 'OUTCAR').replace('CONTCAR',
                                                               'OUTCAR'))
    # Try the same but with compressed files
    fsc = []
    for fn in fnames:
        fsc.append(fn + '.gz')
        fsc.append(fn + '.bz2')
    for f in fsc:
        fnames.append(f)
    # Finally try anything with POTCAR or OUTCAR in the name
    vaspdir = op.dirname(posfname)
    fs = glob.glob(vaspdir + '*POTCAR*')
    for f in fs:
        fnames.append(f)
    fs = glob.glob(vaspdir + '*OUTCAR*')
    for f in fs:
        fnames.append(f)

    tried = []
    files_in_dir = os.listdir('.')
    for fn in fnames:
        if fn in files_in_dir:
            tried.append(fn)
            at = get_atomtypes(fn)
            if len(at) == numsyms:
                return at

    raise IOError('Could not determine chemical symbols. Tried files ' 
                  + str(tried))


def get_atomtypes_from_formula(formula):
    """Return atom types from chemical formula (optionally prepended
    with and underscore).
    """
    from ase.atoms import string2symbols
    symbols = string2symbols(formula.split('_')[0])
    atomtypes = [symbols[0]]
    for s in symbols[1:]:
        if s != atomtypes[-1]: atomtypes.append(s)
    return atomtypes

def _get_nlines(f, n, skip=0, byte_offset=None):
    """Get the content of n lines"""
    import numpy as np
    # Consider using itertools.islice if it raises problems
    # Avoid f.readline().split()
    # Change offset if needed
    if byte_offset is not None:
        f.seek(byte_offset)
    # Skip n lines
    for _ in xrange(skip):
        f.readline()
    # Get the data block borders
    pos = f.tell()
    for _ in xrange(n):
        f.readline()
    end = f.tell()
    f.seek(pos)
    # Return the data using only one call to ".split()" (faster)
    return np.array(f.read(end-pos).split())

def read_vasp(filename='CONTCAR'):
    """Import POSCAR/CONTCAR type file.

    Reads unitcell, atom positions and constraints from the POSCAR/CONTCAR
    file and tries to read atom types from POSCAR/CONTCAR header, if this 
    fails the atom types are read from OUTCAR or POTCAR file.
    """
 
    from ase import Atoms, Atom
    from ase.constraints import FixAtoms, FixScaled
    from ase.data import chemical_symbols
    from itertools import repeat, chain
    import numpy as np

    # Open the file
    if isinstance(filename, str):
        # Assume filename is the name of a file
        f = open(filename, 'rb')
    else:
        # Falling back assuming filename is a file object
        f = filename

    ## Reading the file
    basis_vectors, atom_symbols = _get_pos_hdrs(f)
    tot_natoms = len(atom_symbols)

    # Check if 'selective dynamics' is switched on
    sdyn = f.readline()
    has_selective_dynamics = (sdyn[0].lower() == "s")

    if has_selective_dynamics:
        ac_type = f.readline()
        _data = _get_nlines(f, tot_natoms).reshape(tot_natoms, 6)
        atoms_pos = _data[0:tot_natoms, 0:3].astype(float)
        # 'F': Fixed and 'T': Free
        # selective_flags: True if coordinate is FIXED ("F")
        selective_flags = (_data[0:tot_natoms, 3:6] == "F")
    else:
        ac_type = sdyn
        atoms_pos = \
            _get_nlines(f, tot_natoms).reshape(tot_natoms, 3).astype(float)
        
    try:
        time_step = None
        ac_v_type = f.readline()
        # VASP velocities are A/fs, while ASE want m/s
        velocities = \
    _get_nlines(f, tot_natoms).reshape(tot_natoms, 3).astype(float)
        try:
            time_step = float(_get_nlines(f, 1, skip=2)[0])
        except (ValueError, IndexError):
            pass
    except ValueError: 
        velocities = None
    
    ## Done with all reading
    try:
        f.close()
    except AttributeError:
        pass
    finally:
        del f
    
    # Check if coordinate sets are cartesian or direct
    is_cartesian = (ac_type[0].lower() in "ck")
  
    # Create the Atoms object, set positions and velocities
    atoms = Atoms(symbols = atom_symbols, cell = basis_vectors, pbc = True)
    if is_cartesian:
        atoms.set_positions(atoms_pos*scaling_factor)
    else:
        atoms.set_scaled_positions(atoms_pos)
    del atoms_pos
    if velocities is not None:
        # "Cartesian" is default if nothing is given
        is_v_direct = (ac_v_type[0].lower() in "d")
        if is_v_direct:
            # Transform to Cartesian and scale (vectors/time_step -> A/fs)
            atoms.set_velocities(velocities.dot(basis_vectors)/time_step)
        else:
            atoms.set_velocities(velocities)
    
    atoms.info["time_step"] = time_step or 1

    
    if has_selective_dynamics:
        constraints = []
        indices = []

        for ind, sflags in enumerate(selective_flags):
            if sflags.any() and not sflags.all():
                constraints.append(FixScaled(atoms.get_cell(), ind, sflags))
            elif sflags.all():
                indices.append(ind)
        if indices:
            constraints.append(FixAtoms(indices))
        if constraints:
            atoms.set_constraint(constraints)
    return atoms

def _get_pos_hdrs(f):
    """Gets through POSCAR/CONTCAR/XDATCAR headers
    and returns important informations"""
    
    from itertools import repeat, chain
    from numpy.linalg import det
    
    # Headers are at the begining of the file
    f.seek(0)

    # Assume line 1 gives atom types (VASP 4.x format, check for that later)
    atomtypes = f.readline().split()

    # Get the "lattice constant" (aka scaling factor)
    scaling_factor = float(f.readline())
    basis_vectors =  _get_nlines(f, 3).astype(float).reshape(3,3)
    
    if scaling_factor < 0:
        # The "lattice constant" is the cell volume
        basis_vectors *= (-scaling_factor/abs(det(lattice)))**(1./3)
    else:
        # The "lattice constant" is a scaling factor
        basis_vectors *= scaling_factor
    
    # Get the respective types and numbers of atoms
    try:
        # Assuming VASP 4.x format:
        numofatoms = f.readline()
        int(numofatoms.split()[0])
        # In case 1st line does not contain all symbols
        # Check for "CoP3_In-3.pos"-like syntax
        numsyms = len(numofatoms)
        numtypes = len(atomtypes)
        if numtypes < numsyms:
            if numtypes == 1 and '_' in atomtypes[0]:
                atomtypes = get_atomtypes_from_formula(atomtypes[0])
            else:
                #Nothing was found: check OUTCAR and POTCAR
                atomtypes = atomtypes_outpot(f.name, numsyms)
        else:
            # Too many types, too few atoms
            try:
                for atype in atomtypes[:numsyms]:
                    if not atype in chemical_symbols:
                        raise KeyError
            except KeyError:
                atomtypes = atomtypes_outpot(f.name, numsyms)
    except ValueError:
        # Falling back to VASP 5.x format (5th line: atom symbols):
        atomtypes = numofatoms.split()
        numofatoms = f.readline()
    finally:
        # Remove commented part of the line (if any)
        numofatoms = numofatoms.split("!")[0]
        # Get the actual number of atoms
        numofatoms = [int(num) for num in numofatoms.split()]

    # Build atom list: ["Si", "Fe"] + [2, 3] --> ["Si", "Si", "Fe", "Fe", "Fe"]
    # itertools.chain([[1, 2, 3], [4, 5]]) --> [1, 2, 3, 4, 5]
    atom_symbols = list(
            chain(*[[el[0]]*el[1] for el in zip(atomtypes, numofatoms)])
            )
    
    return basis_vectors, atom_symbols

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
                constraints.append(FixScaled(atoms.get_cell(), ind, sflags))
            elif sflags.all():
                indices.append(ind)
        if indices:
            constraints.append(FixAtoms(indices))
        if constraints:
            atoms.set_constraint(constraints)
    return atoms

def read_vasp_out(filename='OUTCAR',index = -1):
    """Import OUTCAR type file.

    Reads unitcell, atom positions, energies, and forces from the OUTCAR file
    and attempts to read constraints (if any) from CONTCAR/POSCAR, if present. 
    """
    import os
    import numpy as np
    from ase.calculators.singlepoint import SinglePointCalculator
    from ase import Atoms, Atom

    try:          # try to read constraints, first from CONTCAR, then from POSCAR
        constr = read_vasp('CONTCAR').constraints
    except:
        try:
            constr = read_vasp('POSCAR').constraints
        except:
            constr = None

    if isinstance(filename, str):
        f = open(filename)
    else: # Assume it's a file-like object
        f = filename
    data    = f.readlines()
    natoms  = 0
    images  = []
    atoms   = Atoms(pbc = True, constraint = constr)
    energy  = 0
    species = []
    species_num = []
    symbols = []
    ecount = 0
    poscount = 0
    magnetization = []

    for n,line in enumerate(data):
        if 'POTCAR:' in line:
            temp = line.split()[2]
            for c in ['.','_','1']:
                if c in temp:
                    temp = temp[0:temp.find(c)]
            species += [temp]
        if 'ions per type' in line:
            species = species[:len(species)/2]
            temp = line.split()
            for ispecies in range(len(species)):
                species_num += [int(temp[ispecies+4])]
                natoms += species_num[-1]
                for iatom in range(species_num[-1]): symbols += [species[ispecies]]
        if 'direct lattice vectors' in line:
            cell = []
            for i in range(3):
                temp = data[n+1+i].split()
                cell += [[float(temp[0]), float(temp[1]), float(temp[2])]]
            atoms.set_cell(cell)
        if 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM' in line:
            energy = float(data[n+4].split()[6])
            if ecount < poscount:
                # reset energy for LAST set of atoms, not current one - VASP 5.11? and up
                images[-1].calc.energy = energy
            ecount += 1
        if 'magnetization (x)' in line:
            magnetization = []
            for i in range(natoms):
                magnetization += [float(data[n + 4 + i].split()[4])]
        if 'POSITION          ' in line:
            forces = []
            for iatom in range(natoms):
                temp    = data[n+2+iatom].split()
                atoms  += Atom(symbols[iatom],[float(temp[0]),float(temp[1]),float(temp[2])])
                forces += [[float(temp[3]),float(temp[4]),float(temp[5])]]
                atoms.set_calculator(SinglePointCalculator(energy,forces,None,None,atoms))
            images += [atoms]
            if len(magnetization) > 0:
                images[-1].calc.magmoms = np.array(magnetization, float)
            atoms = Atoms(pbc = True, constraint = constr)
            poscount += 1


    # return requested images, code borrowed from ase/io/trajectory.py
    if isinstance(index, int):
        return images[index]
    else:
        step = index.step or 1
        if step > 0:
            start = index.start or 0
            if start < 0:
                start += len(images)
            stop = index.stop or len(images)
            if stop < 0:
                stop += len(images)
        else:
            if index.start is None:
                start = len(images) - 1
            else:
                start = index.start
                if start < 0:
                    start += len(images)
            if index.stop is None:
                stop = -1
            else:
                stop = index.stop
                if stop < 0:
                    stop += len(images)
        return [images[i] for i in range(start, stop, step)]

def write_vasp(filename, atoms, label='', direct=False, sort=None, symbol_count = None, long_format=True, vasp5=False):
    """Method to write VASP position (POSCAR/CONTCAR) files.

    Writes label, scalefactor, unitcell, # of various kinds of atoms,
    positions in cartesian or scaled coordinates (Direct), and constraints
    to file. Cartesian coordiantes is default and default label is the 
    atomic species, e.g. 'C N H Cu'.
    """
    
    import numpy as np
    from ase.constraints import FixAtoms, FixScaled

    if isinstance(filename, str):
        f = open(filename, 'w')
    else: # Assume it's a 'file-like object'
        f = filename
    
    if isinstance(atoms, (list, tuple)):
        if len(atoms) > 1:
            raise RuntimeError("Don't know how to save more than "+
                               "one image to VASP input")
        else:
            atoms = atoms[0]

    # Write atom positions in scaled or cartesian coordinates
    if direct:
        coord = atoms.get_scaled_positions()
    else:
        coord = atoms.get_positions()

    if atoms.constraints:
        sflags = np.zeros((len(atoms), 3), dtype=bool)
        for constr in atoms.constraints:
            if isinstance(constr, FixScaled):
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

    # Create a list sc of (symbol, count) pairs
    if symbol_count:
        sc = symbol_count
    else:
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
            label += '%2s ' % sym
    f.write(label + '\n')

    # Write unitcell in real coordinates and adapt to VASP convention 
    # for unit cell
    # ase Atoms doesn't store the lattice constant separately, so always
    # write 1.0.
    f.write('%19.16f\n' % 1.0)
    if long_format:
        latt_form = ' %21.16f'
    else:
        latt_form = ' %11.6f'
    for vec in atoms.get_cell():
        f.write(' ')
        for el in vec:
            f.write(latt_form % el)
        f.write('\n')

    # If we're writing a VASP 5.x format POSCAR file, write out the
    # atomic symbols
    if vasp5:
        for sym, c in sc:
            f.write(' %3s' % sym)
        f.write('\n')

    # Numbers of each atom
    for sym, count in sc:
        f.write(' %3i' % count)
    f.write('\n')

    if atoms.constraints:
        f.write('Selective dynamics\n')

    if direct:
        f.write('Direct\n')
    else:
        f.write('Cartesian\n')

    if long_format:
        cform = ' %19.16f'
    else:
        cform = ' %9.6f'
    for iatom, atom in enumerate(coord):
        for dcoord in atom:
            f.write(cform % dcoord)
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
