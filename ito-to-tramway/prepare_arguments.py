"""
Create a list of all trajectories in a folder to put them in the arguments list and process them one by one
"""


from constants import data_folder, args_file
import glob
import os


def prepare_arguments():
    file_list = glob.glob(data_folder + "*.csv")
    with open(args_file, 'w') as file_object:
        for i in range(len(file_list)):
            args_string = '-f=%s\n' % (os.path.normpath(file_list[i]))
            # args_string.replace('\\', '/')
            file_object.write(args_string)
