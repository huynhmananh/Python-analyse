from tkinter import Tk, filedialog
import glob
import os
import pandas as pd

root = Tk()
# pointing root to Tk() to use it as Tk() in program.
root.withdraw()
# Hides small tkinter window.
root.attributes('-topmost', True)
# Opened windows will be active. above all windows despite of selection.
open_file = filedialog.askdirectory()

os.chdir(open_file)
file_list = []
for file in glob.glob("*.csv"):
    file_list.append(file)
# Returns opened path as str
print(file_list)

total_average = {}
total_average_df = pd.DataFrame(total_average)

for file_climbing in file_list:
    data_climbing = pd.read_csv(file_climbing)
    average = data_climbing.mean(axis=1)


    name = os.path.splitext(file_climbing)[0]
    new_name = name.replace("Climbing Result ", "")
    new_column = pd.DataFrame({new_name: []})
    new_column[new_name] = average

    total_average_df = pd.concat([total_average_df, new_column], axis=1, ignore_index=False)

save_name = open_file + '/summary_result.csv'
if os.path.exists(save_name):
    os.remove(save_name)
total_average_df.to_csv(save_name)
print("FINISH")
