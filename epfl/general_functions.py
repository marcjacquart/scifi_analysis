import os

def create_folder_if_missing(path):
    if not os.path.isdir(path):
        os.mkdir(path)
        print('New directory created: '+path)
