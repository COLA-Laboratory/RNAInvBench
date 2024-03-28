import pandas as pd
import os

directory = "data/"

# Iterate over each file in the directory
for file_name in os.listdir(directory):
    # Check if the file ends with ".plk.gz"
    if file_name.endswith(".plk.gz"):
        data = pd.read_pickle(directory+file_name)
        print(file_name)
        print(data)
        print("Max Length: ", (max(len(lst) for lst in data["sequence"])-1))
        print("Min length: ", (min(len(lst) for lst in data["sequence"])-1))
        print("Number of Samples: ", len(data))

