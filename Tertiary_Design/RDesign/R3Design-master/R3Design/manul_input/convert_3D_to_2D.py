import moderna
from get_secondary_structure import get_secondary_structure
import argparse
import os
import json

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a directory of 3D RNA into 2D RNA")
    parser.add_argument('--directory', type=str, required=True, help='The directory containing the 3D RNAs as PDB files')
    parser.add_argument("--save_file", type=str, required=True, help="The name of the saved file, no file extension.")
    args = parser.parse_args()
    data_dir = args.directory
    save_file_name = args.save_file
    files = os.listdir(data_dir)
    files = [x for x in files if x[-4:] == ".pdb"]
    parsed_data = []
    for file in files:
        ss = get_secondary_structure(data_dir+"/"+file)
        parsed_data.append({"pdb": file, "ss": ss})
    with open(save_file_name+".txt", "w") as file:
        for row in parsed_data:
            file.write(json.dumps(row)+"\n")