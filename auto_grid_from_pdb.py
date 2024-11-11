from rdkit import Chem
import numpy as np
import tkinter as tk
from tkinter import filedialog
import os

exhaustiveness_time = 20
def ask_input_dir():
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(title="Select Input Directory")
    if directory:
        print(f"Selected directory: {directory}")
    else:
        print("No directory selected")
    return directory

def get_receptor_coordinates(mol):
    """Extract coordinates from the receptor molecule."""
    conformer = mol.GetConformer()
    coords = []
    for atom in mol.GetAtoms():
        pos = conformer.GetAtomPosition(atom.GetIdx())
        coords.append([pos.x, pos.y, pos.z])
    return np.array(coords)

def calculate_grid_box(coords):
    """Calculate the center and size of the grid box."""
    center = coords.mean(axis=0)
    size = coords.ptp(axis=0) + 2  # Adding a margin of 10 Å
    return center, size

def write_config_file(center, size, output_file='config.txt'):
    """Write the Vina configuration to a text file."""
    with open(output_file, 'w') as f:
        f.write("center_x = {:.2f}\n".format(center[0]))
        f.write("center_y = {:.2f}\n".format(center[1]))
        f.write("center_z = {:.2f}\n".format(center[2]))
        f.write("size_x = {:.2f}\n".format(size[0]))
        f.write("size_y = {:.2f}\n".format(size[1]))
        f.write("size_z = {:.2f}\n".format(size[2]))
        f.write(f"exhaustiveness = {exhaustiveness_time}")
    print("Configuration written to", output_file)
    print(f"Configuration written to {output_file}")

pdb_dir = ask_input_dir()
pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith(".pdb")]

for pdb_file in pdb_files:
    pdb_file_path = os.path.join(pdb_dir, pdb_file)
    mol = Chem.MolFromPDBFile(pdb_file_path, removeHs=False)
    coords = get_receptor_coordinates(mol)
    center, size = calculate_grid_box(coords)
    output_file = pdb_file.replace('.pdb', '.txt')
    output_file = "config_" + output_file
    output_file_path = os.path.join(pdb_dir, output_file)
    write_config_file(center, size, output_file=output_file_path)




