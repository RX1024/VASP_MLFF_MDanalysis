import matplotlib.pyplot as plt
import numpy as np

def parse_atom_indices(input_str):
    atoms = set()
    for part in input_str.split(','):
        if '-' in part:
            start, end = map(int, part.split('-'))
            atoms.update(range(start, end + 1))
        else:
            atoms.add(int(part))
    return sorted(atoms)

def read_xdatcar(atom_indices):
    coordinates = {atom_index: [] for atom_index in atom_indices}
    a = 0
    with open("XDATCAR", "r") as file:
        lines = file.readlines()
        for i in range(len(lines)):
            line = lines[i].strip()
            if "IL_liquid" in line:
                a = float(lines[i + 2].split()[0])
            elif "Direct configuration" in line:
                if a == 0:
                    raise ValueError("Cell parameter 'a' not found before coordinates.")
                for atom_index in atom_indices:
                    coord_line = lines[i + atom_index].split()
                    coordinates[atom_index].append([float(coord) * a for coord in coord_line[:3]])
    return coordinates

def calculate_msd(coordinates):
    msd_values = []
    for step in range(1, len(next(iter(coordinates.values())))):
        displacements = [np.linalg.norm(np.array(coordinates[atom_index][step]) - np.array(coordinates[atom_index][0])) for atom_index in coordinates]
        msd = np.mean([d**2 for d in displacements])
        msd_values.append(msd)
    return msd_values

def plot_msd(msd_values):
    steps = range(1, len(msd_values) + 1)
    plt.figure()
    plt.plot(steps, msd_values, marker='o')
    plt.xlabel('Step')
    plt.ylabel('Mean Squared Displacement')
    plt.title('MSD of Atoms Over Time')
    plt.grid(True)
    plt.show()

def main():
    atom_indices_input = input("Enter the atom indices (e.g., 1,2,3 or 1-3 or 1-3,5-7): ")
    atom_indices = parse_atom_indices(atom_indices_input)
    coordinates = read_xdatcar(atom_indices)
    msd_values = calculate_msd(coordinates)
    plot_msd(msd_values)
    print("MSD plot generated.")

if __name__ == "__main__":
    main()
