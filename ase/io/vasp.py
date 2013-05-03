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
    
    b_vel_ptr = f.tell()
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
        f.seek(b_vel_ptr)
    
    ## Done with all reading
    
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

def _get_read_positions_from_xdat(filename='XDATCAR', index=-1):
    
    from os.path import getsize
    
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
    # NO f.readlines() !!!
    # Too memory consumming and slow random access
    
    # Get the offset of the 1st block
    beg_block = f.tell()
    # +1 is for the "Direct configuration=    1" line
    _get_nlines(f, 0, skip=tot_natoms+1)
    # Get the offset of the 2nd block
    end_block = f.tell()
    # Get position data block size in bytes
    block_size = end_block - beg_block
    
    #Check boundaries to allow negative indexes
    f_size = getsize(f.name)
    # Get the total number of structures of file
    maxnumofstruct = (f_size - beg_block) // block_size
    
    # Check if a slice object was given
    try:
        start = index.start or 0
        step = index.step or 1
        stop = index.stop or maxnumofstruct 
    except AttributeError:
        start = index
        step = 1
        stop = index + 1 if start != -1 else maxnumofstruct
    
    if start < 0:
        start += maxnumofstruct
    if (step is not None) and (step < 0):
        # Don't want to bother with negative steps
        step = - step
    if (stop is not None) and (stop < 0):
        stop += maxnumofstruct
    
    # Returns the positions asked for
    return [_get_nlines(
            f, tot_natoms, skip=1, byte_offset=i*block_size+beg_block
            ).reshape(tot_natoms,3).astype("float") 
            for i in xrange(start, stop, step)]

def read_vasp_xdat(filename='XDATCAR', index=-1):
    """Import POSCAR/CONTCAR type file.

    Reads unitcell, and atom positions from the XDATCAR file
    and read contrains from CONTCAR/POSCAR file if any.
    """
    
    from ase import Atoms, Atom
    from ase.constraints import FixAtoms, FixScaled
    from ase.data import chemical_symbols
    from itertools import repeat, chain, izip
    import numpy as np

    # Open the file
    if isinstance(filename, str):
        # Assume filename is the name of a file
        f = open(filename, 'rb')
    else:
        # Falling back assuming filename is a file object
        f = filename
    
    basis_vectors, atom_symbols = _get_pos_hdrs(f)
    tot_natoms = len(atom_symbols)
    
    # Try to read constraints and last velocities 
    # first from CONTCAR, then from POSCAR
    try:      
        _file = read_vasp('CONTCAR')
        constr = _file.constraints
        veloci = _file.get_velocities()
    except:
        try:
            _file = read_vasp('POSCAR')
            constr = _file.constraints
            veloci = _file.get_velocities()
        except:
            constr = None
            veloci = None
    finally:
        time_step = _file.info.get("time_step") or 1
        del _file
    
    pos = _get_read_positions_from_xdat(f, index)
    images = [Atoms(symbols = atom_symbols, cell = basis_vectors, pbc = True)
              for i in repeat(None, len(pos))]
        
    for image, p in izip(images, pos):
        image.set_scaled_positions(p)
    
    ## Done with all reading
    try:
        f.close()
    except AttributeError:
        pass
    finally:
        del f
    
    if constr is not None:
        constraints = []
        indices = []
        for ind, sflags in enumerate(constr):
            if sflags.any() and not sflags.all():
                constraints.append(FixScaled(atoms.get_cell(), ind, sflags))
            elif sflags.all():
                indices.append(ind)
        if indices:
            constraints.append(FixAtoms(indices))
        if constraints:
            atoms.set_constraint(constraints)
    return images

def _order_of_appearance(f, search_for=[], breaker=None, chunk_size=32768):
    """Find the order of appearance of the terms searched for.
    The file pointer is *NOT* moved after using this function. 

    Parameters:

    f: file object
    search_for: list of strings/bytes in whatever order that are searched for
    breaker: stop the search after the current chunk
    chunk_size: size of a chunk of data in bytes (Defaults: 32768)

    Returns:

    list: searched terms sorted in order of appearance
    terms not encountered are removed
    """
    mxl = max(map(len, search_for))
    if chunk_size < mxl:
        raise ValueError("Don't choose chunk_size smaller than the search term")
    idx = range(len(search_for))
    origin = f.tell()
    start = origin
    out = [None]*len(search_for)
    while None in out:
        _data = f.read(chunk_size)
        if not _data:
            break
        # zip is needed. izip is not [].remove-safe
        for i, term in zip(idx, search_for):
            try:
                out[i] = (term, _data.index(term) + start)
                search_for.remove(term)
                idx.remove(i)
            except ValueError:
                pass
        if breaker is not None:
            try:
                _data.index(breaker)
                break
            except ValueError:
                pass
        start += chunk_size - mxl + 1
        f.seek(start)
    f.seek(origin)
    return [i[0] for i in sorted([i for i in out if i is not None], key=lambda el: el[1])]

def _go_next(f, search_for="", chunk_size=32768):
    """Find the next occurence of a string/bytes in a file.
    The file pointer is *MOVED* to the location of the searched term.

    Parameters:

    f: file object
    search_for: string/array of bytes to be searched for
    chunk_size: size of a chunk of data in bytes (Defaults: 32768)

    Returns:

    None: if 'search_for' has not been found
    Current pointer location: if 'search_for' has been found
    """
    ml = len(search_for)
    if chunk_size < ml:
        raise ValueError("Don't choose chunk_size smaller than the search term")
    start = f.tell() + 1
    f.seek(start)
    while True:
        # Get the data
        _data = f.read(chunk_size)
        if not _data:
            return None
        # Search for term
        try:
            _pos = _data.index(search_for)
            f.seek(_pos + start)
            return _pos + start
        except ValueError:
            pass
        # +1 : avoid double counting: 
        # "... seekable" --> "eekable ..."
        # "... seekabl" --> "seekable ..."
        start += chunk_size - ml + 1 
        # Change the position for the next step
        f.seek(start)

def _count_and_go(f, search_for="", n=1, chunk_size=32768):
    """Go to the nth occurence of the term searched for.
    The file pointer is *MOVED* to the location of the searched term.

    Parameters:

    f: file object
    search_for: string/array of bytes to be searched for
    n: counter
    chunk_size: size of a chunk of data in bytes (Defaults: 32768)

    Returns:

    None: if 'search_for' has not been found
    Current pointer location: if 'search_for' has been found
    """
    ml = len(search_for)
    if chunk_size < ml:
        raise ValueError("Don't choose chunk_size smaller than the search term")
    start = f.tell()
    f.seek(start)
    while n > 0:
        # Get the data
        _data = f.read(chunk_size)
        if not _data:
            return None
        # Search for term
        n -= _data.count(search_for)
        if n > 0:
            start += chunk_size - ml + 1
        f.seek(start)
    _offset = -1
    for _ in xrange(abs(n)+1):
        _offset = _data.index(search_for, _offset+1)
    f.seek(start+_offset)
    return start+_offset
        
def _count(f, search_for="", chunk_size=32768*4):
    """Count the number of times the searched term appear in a file from 
    current file pointer location to the end.
    The file pointer is *NOT* moved to the location of the searched term.

    Parameters:

    f: file object
    search_for: string/array of bytes to be searched for
    chunk_size: size of a chunk of data in bytes (Defaults: 32768*4)

    Returns:

    Number of times 'search_for' have been counted
    """
    ml = len(search_for)
    if chunk_size < ml:
        raise ValueError("Don't choose chunk_size smaller than the search term")
    origin = f.tell()
    start = origin
    f.seek(start)
    _counter = 0
    while True:
        # Get the data
        _data = f.read(chunk_size)
        if not _data:
            f.seek(origin)
            return _counter
        # Search for term
        _counter += _data.count(search_for)
        start += chunk_size - ml + 1
        f.seek(start)

def read_vasp_out(filename='OUTCAR', index=-1):
    """Import OUTCAR type file.

    Reads unitcell, atom positions, energies, and forces from the OUTCAR file
    and attempts to read constraints (if any) from CONTCAR/POSCAR, if present. 
    """
    import os
    import numpy as np
    from ase.calculators.singlepoint import SinglePointCalculator
    from ase import Atoms, Atom
    from itertools import chain, repeat, izip
    
    #1) Read contraints from CONTCAR/POSCAR if present
    #2) Use the keys in 'parsing_dct' as search terms to 
    # get the search terms ordering in the current OUTCAR
    # (for instance, do we have Energy -> Position or
    # Magnetization -> Position -> Energy or ...?)
    #3) Use 'parsing_dct' to create a new list containing
    # the ordered set of parsing functions ('step_func')
    #4) Get information from the headers of the OUTCAR
    # (atom names, multiplicity, unit cell)
    #5) Get boundaries to allow negative index searches
    #6) Create a list containing intialized "Atoms" objects
    #7) Iterate over 'step_func' and update "Atoms"
    
    # The beginning of each loop is not reliably placed in the file
    # (no contant offset unfortunately)
    
    # A couple of functions were created to allow fast reading/parsing
    
    
    # Define parsing functions
    def get_pos(f, atoms):
        if _go_next(f, "POSITION          ") is not None:
            _data = _get_nlines(f, natoms, skip=2).reshape(natoms, 6).astype(float)
            atoms.set_positions(_data[:natoms, :3])
            atoms.calc.forces = _data[:natoms, 3:]
    def get_ene(f, atoms):
        if _go_next(f, "FREE ENERGIE OF THE ION-ELECTRON SYSTEM") is not None:
            _get_nlines(f, 0, skip=4)
            # return energy
            atoms.calc.energy = float(f.readline().split()[-1])
    def get_mag(f, atoms):
        if _go_next(f, "magnetization (x)") is not None:
            _t = _get_nlines(f, natoms, skip=4)
            atoms.calc.magmoms = _t.reshape(natoms, 6)[:,5].astype(float)
    
    # Define parsing functions in a dict:
    parsing_dct = {
           "POSITION          ": get_pos,
           "FREE ENERGIE OF THE ION-ELECTRON SYSTEM": get_ene,
           "magnetization (x)" : get_mag
           }
           
    search_terms = list(parsing_dct.keys())

    # Try to read constraints, first from CONTCAR, then from POSCAR
    try:          
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
        
    step_func = [parsing_dct[key] for key in _order_of_appearance(f, search_for=search_terms, breaker="Iteration    2")]
    
    # Get the atom names
    _pos = _go_next(f, "POTCAR:")
    if _pos is not None:
        line = f.readline()
        counter = 0
        while "POTCAR:" in f.readline():
            counter +=1
        f.seek(_pos)  
        species = [f.readline().split()[2][:2] for i in repeat(None, counter)]
        species = [i if i[-1] not in "_.1 " else i[0] for i in species]
    
    # Get the associated multiplicity
    if _go_next(f, "ions per type =") is not None:
        species_num = [int(i) for i in f.readline()[len("ions per type ="):].split()]
        natoms = sum(species_num)
    
        # Construct the list of atom symbols 
        symbols = list(
           chain(*[[el[0]]*el[1] for el in zip(species, species_num)])
               )
    
    # Get the unit cell
    if _go_next(f, "direct lattice vectors") is not None:
        cell = _get_nlines(f, 3, skip=1).reshape(3,6)[:3,:3].astype(float)
    
    # Get bounds
    num_calc = _count(f, "POSITION      ")
    
    try:
        # Steps are ignored for now
        start = index.start or 0
        stop = index.stop or num_calc
    except AttributeError:
        start = index
        if start == -1:
            stop = num_calc
        else:
            stop = start + 1
    except TypeError:
        start = 0
        stop = num_calc
        
    if start < 0:
        start += num_calc
    if stop < 0:
        stop += num_calc
    
    # No step for now...
    l_slice = (stop-start)#//step == 1
    
    # Initialize the list of atoms. [].append() are expensive!
    images = [Atoms(symbols=symbols, cell=cell, pbc=True, constraint=constr)
              for i in repeat(None, l_slice)] 
    
    # Get started
    if start < 9999:
        _pos = _go_next(f, "Iteration%5d"%(start+1))
    else:
        # Then we got "Iteration ****(   1)" and can't use it
        _pos = _count_and_go(f, "(   1)  -----", start+1)
        
    if _pos is None:
        raise IndexError("Index %s do not exist (max: %s)"%(start, num_calc-1))
    
    for atoms in images:
        atoms.set_calculator(SinglePointCalculator(None,None,None,None,atoms))
        for function in step_func:
            function(f, atoms)
        
        atoms.calc.atoms = atoms

    return images

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

    velo = atoms.get_velocities()

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

    # Write velocities from atoms.get_velocities()
    if velo is not None:
        f.write(" \n")
        np.savetxt(f, velo, fmt="%16.8e%16.8e%16.8e")
        f.write(' \n')

    if type(filename) == str:
        f.close()
