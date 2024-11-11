from rdkit import Chem
from rdkit.Chem import AllChem
import tkinter as tk
from tkinter import filedialog
import os
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy

def ask_input_dir():
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(title="Select Input Directory")
    if directory:
        print(f"Selected directory: {directory}")
    else:
        print("No directory selected")
    return directory

def is_3d(mol):
    if mol.GetNumConformers() > 0:
        conf = mol.GetConformer()
        pos = conf.GetPositions()
        for atom_pos in pos:
            if atom_pos[2] != 0.0:
                return True
    return False

def mol2d_to_3d(mol):
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xf00d
    AllChem.EmbedMolecule(mol, params)
    return mol

ligands_dir = ask_input_dir()
sdf_files = [f for f in os.listdir(ligands_dir) if f.endswith(".sdf")]
preparator = MoleculePreparation()

for sdf_file in sdf_files:
    sdf_file_path = os.path.join(ligands_dir, sdf_file)
    mols_from_sdf = Chem.SDMolSupplier(sdf_file_path, removeHs=False)
    for mol in mols_from_sdf:
        if mol is not None:  # If molecule is valid
            if not is_3d(mol):
                mol = mol2d_to_3d(mol)
    mol_addHs = Chem.AddHs(mol)
    mol_setups = preparator.prepare(mol_addHs)
    for setup in mol_setups:
        pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
    new_name = sdf_file.replace('.sdf', '.pdbqt')
    new_name = 'ligand_' + new_name
    new_file = os.path.join(ligands_dir, new_name)
    with open(new_file,'w') as fw:
        for line in pdbqt_string:
            fw.write(line)