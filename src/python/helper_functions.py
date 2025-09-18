import numpy as np
import ctypes
import pandas as pd
import ctypes
from ctypes import c_char, POINTER, create_string_buffer

angtobohr = 1.8897259886

# Periodic table mapping
periodic_table = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
    'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16,
    'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23,
    'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
    'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
    'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44,
    'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51,
    'Te': 52, 'I': 53, 'Xe': 54
}

def symbols_to_z(species):
    return [periodic_table[s] for s in species]


def read_xyz(filename):
    """
    Reads an .xyz file and returns:
    - natoms: int, number of atoms
    - coords: np.ndarray of shape (3, natoms), atomic coordinates
    - species: np.ndarray of shape (natoms,), element symbols as strings
    """
    filename = f"./data/molecules/{filename}"
    with open(filename, 'r') as file:
        lines = file.readlines()

    natoms = int(lines[0].strip())
    # lines[1] is the optional comment line, can be ignored or returned if needed

    species = []
    coords = []

    for line in lines[2:2 + natoms]:
        parts = line.strip().split()
        if len(parts) != 4:
            raise ValueError(f"Invalid line in XYZ file: {line}")
        species.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        
    
    coords = np.array(coords).T  # Transpose to shape (3, natoms)
    coords = coords * angtobohr # transform to bohr
    species = np.array(species, dtype=str)

    return natoms, coords, species


# Optional known exceptions
known_multiplicities = {
    ('O', 2): 3,  # O2 molecule
    ('N', 2): 4   # N2 molecule
}

def guess_multiplicity(atoms, charge=0):
    # Calculate total electrons
    total_electrons = sum(periodic_table[atom] for atom in atoms) - charge

    # Look for known exceptions
    atom_counts = {atom: list(atoms).count(atom) for atom in set(atoms)}
    for (symbol, count), mult in known_multiplicities.items():
        if atom_counts.get(symbol, 0) == count and len(atom_counts) == 1:
            return mult

    # Default guess
    return 1 if total_electrons % 2 == 0 else 2

def compute_spin_electrons(atoms, charge=0):
    """Returns (n_alpha, n_beta, multiplicity)"""
    total_electrons = sum(periodic_table[atom] for atom in atoms) - charge
    multiplicity = guess_multiplicity(atoms, charge)
    S = (multiplicity - 1) / 2

    if (total_electrons + 2*S) % 2 != 0:
        raise ValueError("Inconsistent electron count and spin state.")

    n_alpha = int((total_electrons + 2*S) // 2)
    n_beta = total_electrons - n_alpha
    return n_alpha, n_beta, multiplicity


def run_scf(lib, atom_file, eps_scf, max_iter, final_energy, basis_set, charge = 0):
    from ctypes import c_int, c_double, POINTER, byref, create_string_buffer

    # Read atomic data from XYZ file
    natoms, coords, species_list = read_xyz(atom_file)

    # Compute electronic information
    nalpha, nbeta, _ = compute_spin_electrons(species_list, charge)
    nelectrons = nalpha + nbeta
    

    # Convert species to atomic numbers
    z_num_list = symbols_to_z(species_list)

    # Convert species list to C-compatible string
    species_fixed = [s.ljust(2) for s in species_list]
    flat_bytes = b''.join(s.encode('utf-8') for s in species_fixed)
    species_c_array = (c_char * (2 * natoms))(*flat_bytes)
    
    # Flatten coordinates and convert to C array
    coords_flat = np.array(coords, dtype=np.float64)
    coords_c = coords_flat.ctypes.data_as(POINTER(c_double))

    # Convert z_num list to C array
    z_array = np.array(z_num_list, dtype=np.int32)
    z_num_c = z_array.ctypes.data_as(POINTER(c_int))

    # Prepare scalar arguments
    natoms_c = c_int(natoms)
    nelectrons_c = c_int(nelectrons)
    nalpha_c = c_int(nalpha)
    nbeta_c = c_int(nbeta)
    eps_scf_c = c_double(eps_scf)
    max_iter_c = c_int(max_iter)
    charge_c = c_int(charge)
    final_energy_c = c_double(0.0)
    basis_set_c = basis_set.encode('utf-8')


    # Print all information being passed
    '''
    print("\n=== SCF CALL DEBUG INFO ===")
    print(f"File: {atom_file}")
    print(f"natoms: {natoms}")
    print(f"nelectrons: {nelectrons}")
    print(f"nalpha: {nalpha}")
    print(f"nbeta: {nbeta}")
    print(f"species: {species_list}")
    print(f"atomic numbers: {z_num_list}")
    print(f"coords (Bohr, shape {coords.shape}):\n{coords}")
    print(f"eps_scf: {eps_scf}")
    print(f"max_iter: {max_iter}")
    print(f"basis_set: {basis_set} -> n = {basis_set_n}")
    print("============================\n")'''

    # Call the C wrapper
    lib.run_scf_c(natoms_c, nelectrons_c, charge_c, nalpha_c, nbeta_c, species_c_array, coords_c,
                  eps_scf_c, max_iter_c, byref(final_energy_c), basis_set_c, z_num_c)

    return final_energy_c.value
