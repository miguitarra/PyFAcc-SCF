import numpy as np
import ctypes

def parser_input(filename):
    '''
    This function takes a filename as input and parses the file to extract the number of atoms, an array with the atomic numbers,
    the coordinates of the atoms, the number of basis functions per atom, and the basis function exponents.

    Input:
        - filename: str, path to the input file

    Output:
        - num_atoms: int, number of atoms
        - atomic_numbers: np.ndarray, array of atomic numbers
        - coords: np.ndarray, array of shape (num_atoms, 3) with the coordinates of the atoms
        - num_basisf_atom: np.ndarray, array with the number of basis functions for each atom
        - basisf_exp: np.ndarray, array with the basis function exponents
        
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()

    # First line: number of atoms, number of electrons, total basis functions
    num_atoms, num_electrons, nalpha, nbeta, total_basis_functions = map(int, lines[0].split())

    # Initialize arrays
    coords = np.zeros((num_atoms, 3), dtype=np.float64)  # Shape: num_atoms x 3
    atomic_numbers = []
    num_basisf_atom = []
    basisf_exp = []

    i = 1  # Start reading from the second line
    for atom_idx in range(num_atoms):
        # Parse atomic data (coordinates, nuclear charge, number of basis functions)
        atom_data = lines[i].split()
        x, y, z = map(float, atom_data[:3])  # Cartesian coordinates
        coords[atom_idx] = [x, y, z]

        nuclear_charge = float(atom_data[3])  # Nuclear charge
        atomic_numbers.append(int(nuclear_charge))  # Assuming nuclear charge is the atomic number

        cur_nbasisf = int(atom_data[4])  # Number of basis functions for this atom
        num_basisf_atom.append(cur_nbasisf)

        # Parse basis function exponents
        i += 1
        for _ in range(cur_nbasisf):
            basisf_exp.append(float(lines[i].strip()))
            i += 1

    return num_atoms, num_electrons, nalpha, nbeta, np.array(atomic_numbers, dtype=np.int32), coords, np.array(num_basisf_atom, dtype=np.int32), np.array(basisf_exp, dtype=np.float64)

def print_atom_info(filename):
    """
    This function reads the input file and prints out the parsed information
    about the atoms, including their coordinates, atomic numbers, and basis functions.

    Input:
        - filename: str, path to the input file
    """
    num_atoms, num_electrons, atomic_numbers, coords, num_basisf_atom, basisf_exp = parser_input(filename)

    print(f"Number of atoms: {num_atoms}")
    print(f"Number of electrons: {num_electrons}")
    print("\nAtomic Information:")
    for i in range(num_atoms):
        print(f"Atom {i + 1}:")
        print(f"  Atomic Number: {atomic_numbers[i]}")
        print(f"  Coordinates: {coords[i]}")
        print(f"  Number of Basis Functions: {num_basisf_atom[i]}")
        print(f"  Basis Function Exponents: {basisf_exp[sum(num_basisf_atom[:i]):sum(num_basisf_atom[:i + 1])]}")

def run_scf(lib, natoms, nelectrons, nalpha, nbeta, atomic_numbers, coords, eps_scf, max_iter, nbasis_atom, basis_exp, basis_set = 6):
    # Convert inputs to ctypes
    atomic_numbers = (ctypes.c_int * len(atomic_numbers))(*atomic_numbers)
    coords = (ctypes.c_double * len(coords.flatten()))(*coords.flatten())
    nbasis_atom = (ctypes.c_int * len(nbasis_atom))(*nbasis_atom)
    basis_exp = (ctypes.c_double * len(basis_exp))(*basis_exp)
    final_energy = ctypes.c_double()
    eps_scf = ctypes.c_double(eps_scf)

    # Call the C wrapper
    lib.run_scf_c(natoms, nelectrons, nalpha, nbeta, atomic_numbers, coords, eps_scf, max_iter, 
                      ctypes.byref(final_energy), nbasis_atom, basis_exp, basis_set)

    return final_energy.value
