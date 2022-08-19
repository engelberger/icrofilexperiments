
# This file contains all the files necesessary to run an amber experiment


# Import libraries
import os
import sys


# A function that creates a folder
def create_folder(folder_path):
    """_summary_
    Creates a folder
    
    Args:
        folder_path (_type_): _description_
    
    """
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print("Folder created: " + folder_path)
    else:
        print("Folder already exists: " + folder_path)


def create_system_folders(n_replicas, folder_path,folder_name):
    """_summary_
    Given a number of replicas creates the following folders architecture:
    .
    ├── folder_name  <- Name of the MD experiment
    │   ├── 1    <- folder for replica 1
    │       ├── in <- input files
    │       └── out <- output files
    │   ├── 2    <- folder for replica 2
    │       ├── in <- input files
    │       └── out <- output files
    │   ├── n    <- folder for replica n
    │       ├── in <- input files
    │       └── out <- output files
    
    Args:
        n_replicas (_type_): _description_
        folder_path (_type_): _description_
        folder_name (_type_): _description_
    
    """
    for i in range(1, n_replicas+1):
        create_folder(folder_path + "/" + folder_name + "/" + str(i))
        create_folder(folder_path + "/" + folder_name + "/" + str(i) + "/in")
        create_folder(folder_path + "/" + folder_name + "/" + str(i) + "/out")
   
