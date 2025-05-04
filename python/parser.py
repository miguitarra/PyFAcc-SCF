import numpy as np

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
    num_atoms, num_electrons, total_basis_functions = map(int, lines[0].split())

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

    return num_atoms, np.array(atomic_numbers, dtype=np.int32), coords, np.array(num_basisf_atom, dtype=np.int32), np.array(basisf_exp, dtype=np.float64)
