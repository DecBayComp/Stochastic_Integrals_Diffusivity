import glob
import os
import sys

import pandas as pd
from tqdm import tqdm

# def rename_neurons_data_sets():
"""
Rename the file by using the average location of the points.
Will not process files that have '=' in their names
"""

path = r"\\atlas.pasteur.fr\@Dbc\LAB_shared_stuff\alexander_serov\Synapse_drug pour sacha\Drug_PP2_septembre_2015\CONTROLE"
columns = ['x', 'y', 't']
extension = '.trxyt'

bl_recursive = False
if bl_recursive:
    file_list = [f for f in glob.iglob(os.path.join(
        path, "/**/*" + extension), recursive=True)]
else:
    file_list = [f for f in glob.iglob(os.path.join(
        path, r"*" + extension), recursive=False)]

# folder, _ = os.path.split(file_list[0])
for file in tqdm(file_list):
    if '=' in file:
        continue

    # Split name on the last underscore
    index = file.rfind('_')
    first_part = file[:index]
    second_part = file[index + 1:]

    # Load table
    data = pd.read_csv(file, names=columns, sep='\t')
    xmean, ymean = data[['x', 'y']].mean()
    new_name = f"{first_part}_x={xmean:05.2f}_y={ymean:05.2f}_{second_part}"
    # print(new_name)

    os.rename(file, new_name)

    # sys.exit(0)

    # print(first_part, second_part)
