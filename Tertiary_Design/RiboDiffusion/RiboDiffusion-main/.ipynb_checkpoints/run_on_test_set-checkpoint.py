import numpy as np
import subprocess
import os
import torch
import re
import json
from tqdm import tqdm


DATA_PATH = "/home/jack/Raw_Tertiary_Data/"

def read_data_splits(file_path):
    ids = []
    with open(file_path, "r") as file:
        for line in file.readlines():
            line = line.strip()
            temp = [id+".pdb" for id in line.split(",")]
            ids.append(temp)
    return ids

def get_data_files(train_path, val_path, test_path):
    
    train_list = read_data_splits(train_path)
    val_list = read_data_splits(val_path)
    test_list = read_data_splits(test_path)
    
    # train_list = [x+".pdb" for x in train_list]
    # val_list = [x+".pdb" for x in val_list]
    # test_list = [x+".pdb" for x in test_list]
    
    return train_list, val_list, test_list


if __name__ == "__main__":

    
    script_path = "main.py"

    test_results = []
    pdb_files = [f for f in os.listdir("/home/jack/Tertiary_Models/gRNAde/geometric-rna-design-main/data/DAS_Split/") if f[-4:] == ".pdb"]
    # pdb_files = [f for f in os.listdir("/home/jack/Tertiary_Models/gRNAde/geometric-rna-design-main/data/CASP15/split_pdbs2") if f[-4:] == ".pdb"]
    DATA_PATH = "/home/jack/Tertiary_Models/gRNAde/geometric-rna-design-main/data/DAS_Split/"
    puzzle_name = "DAS_Split"
    # test_list = pdb_files
    seqs = []
    i = 0
    for pdb_files_list in tqdm(pdb_files):
        found_seq = ""
        t_pdb_files = [pdb_files_list]
        single_test_rec = []
        single_test_per = []
        single_test_ssc = []
        for pdb_file in t_pdb_files:
            pdb_file_path = os.path.join(DATA_PATH, pdb_file)
            
            command = [
                "python", script_path, "-mode", "inference", "-PDB_file", pdb_file_path
            ]
            
            # print(f"Running command: {' '.join(command)}")
            
            try:
                result = subprocess.run(command, check=True, capture_output=True, text=True)
                print(f"Successfully processed {pdb_file}")
                output = result.stdout
                recovery_match = re.search(r"recovery_rate (\d+\.\d+)", output)
                recovery_rate = float(recovery_match.group(1)) if recovery_match else None
                perplexity_match = re.search(r"perplexity (\d+\.\d+)", output)
                perplexity = float(perplexity_match.group(1)) if perplexity_match else None
                ssc_match = re.search(r"secondary_struct_consistency (-?\d+\.\d+)", output)
                ssc = float(ssc_match.group(1)) if ssc_match else None
                if recovery_rate != None:
                    single_test_rec.append(recovery_rate)
                if perplexity != None:
                    single_test_per.append(perplexity)
                if ssc != None:
                    single_test_ssc.append(ssc)
                # get sequence from file
                with open("exp_inf/fasta/"+pdb_file[:-4]+"_0.fasta") as file:
                    for line in file.readlines():
                        if line[0] == ">":
                            continue
                        else:
                            seqs.append(line.strip())
                            break

            except subprocess.CalledProcessError as e:
                print(f"Error processing {pdb_file}: {e}")
        if len(single_test_rec) > 0:
            print(f"Test {i} Recovery: {np.mean(single_test_rec)}")
            print(f"Test {i} Perplexity: {np.mean(single_test_per)}")
            print(f"Test {i} SSC: {np.mean(single_test_ssc)}")
            test_results.append({"Recovery": np.mean(single_test_rec), "Perplexity": np.mean(single_test_per),
                                 "SSC": np.mean(single_test_ssc)})
            test_list.append(pdb_file)
        i += 1
    print(f"Overall Recovery: {np.mean([x["Recovery"] for x in test_results])}")
    print(f"Overall Perplexity: {np.mean([x["Perplexity"] for x in test_results])}")
    print(f"Overall SSC: {np.mean([x["SSC"] for x in test_results])}")
    print(f"Writing to JSON file...")
    with open("RiboDiffusionResults"+puzzle_name+".json", "w") as file:
        for result, files, seq in zip(test_results, test_list, seqs):
            output = {"pdb_file":os.path.join(DATA_PATH, pdb_file), "result":result, "seq":seq}
            file.write(json.dumps(output)+"\n")