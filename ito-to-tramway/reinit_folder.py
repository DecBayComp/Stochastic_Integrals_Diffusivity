

import os
import shutil


def reinit_folder(folder):
    """
    Clear folder contents or create folder if necessary.
    """

    if not os.path.isdir(folder):
        os.makedirs(folder)
    else:
        for file in os.listdir(folder):
            file_path = os.path.join(folder, file)
            if os.path.isfile(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)

    print("Folder {folder} cleaned successfully!".format(folder=folder))
