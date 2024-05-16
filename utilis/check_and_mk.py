import os

def check_and_mk(folder):
    # if folder does not exist, make it.
    if not os.path.exists(folder):
        os.makedirs(folder)

