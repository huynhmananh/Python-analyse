dir_vina = "C:/Program Files (x86)/The Scripps Research Institute/Vina/vina.exe"

import subprocess
import os
import time
import tkinter as tk
from tkinter import filedialog
def ask_input_dir():
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(title="Select Input Directory")
    if directory:
        print(f"Selected directory: {directory}")
    else:
        print("No directory selected")
    return directory

def filter_files_by_keyword(directory, keyword):
    filtered_files = []
    for filename in os.listdir(directory):
        if os.path.isfile(os.path.join(directory, filename)):
            if keyword in filename:
                filtered_files.append(filename)
    return filtered_files

def wait_for_file(directory, filename, wait_time=15):
    file_path = os.path.join(directory, filename)
    while not os.path.exists(file_path):
        time.sleep(wait_time)

def run_commands(command):
        try:
            process = subprocess.run(command, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while executing '{command}': {e}")


dir = ask_input_dir()
receptor_list = filter_files_by_keyword(dir, "receptor")
receptor_list = [receptor for receptor in receptor_list if receptor.endswith(".pdbqt")]
ligand_list = filter_files_by_keyword(dir, "ligand")
os.chdir(dir)
for receptor in receptor_list:
    config_file = "config_" + receptor
    config_file = config_file.replace(".pdbqt", ".txt")
    for ligand in ligand_list:
        receptor_name = receptor.replace("receptor_", "")
        receptor_name = receptor_name.replace(".pdbqt", "")
        ligand_name = ligand.replace("ligand_", "")
        ligand_name = ligand_name.replace(".pdbqt", "")
        log_file = f"log_{receptor_name}_{ligand_name}.txt"
        output_file = f"output_{receptor_name}_{ligand_name}.pdbqt"
        command = f'''"{dir_vina}" --receptor {receptor} --ligand {ligand} --config {config_file} --log {log_file} --out {output_file} --exhaustiveness {32} --cpu {8} --num_modes {10} --seed {20242510}'''
        run_commands(command)
        wait_for_file(dir, output_file)
        time.sleep(30)