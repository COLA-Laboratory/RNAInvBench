import subprocess
import os
import re
import argparse
from datetime import datetime
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from torchmetrics.functional.classification import binary_matthews_corrcoef
import torch
from tqdm import tqdm




def predict_sec_struct(seq):
    eternafold_path = "/home/jack/Tertiary_Models/gRNAde/geometric-rna-design-main/tools/EternaFold/EternaFold-master 2/src/contrafold"
    current_datetime = datetime.now().strftime("%Y%m%d_%H%M%S")
    fasta_file_path = f"fastas_temp/temp_{current_datetime}.fasta"
    cmd = [eternafold_path, "predict", fasta_file_path]
    SeqIO.write(SeqRecord(Seq(seq), id="temp"),fasta_file_path, "fasta")
    output = subprocess.run(cmd, check=True, capture_output=True).stdout.decode("utf-8")
    return [x for x in output.split("\n")[-2]]

def fetch_secondary_struct_consistency(real_seq, pred_seq):
    real_sec_struct = predict_sec_struct(real_seq)
    pred_sec_struct = predict_sec_struct(pred_seq)
    pred_sec_struct = torch.tensor([1 if c in "()" else 0 for c in pred_sec_struct], dtype=torch.float32)
    real_sec_struct = torch.tensor([1 if c in "()" else 0 for c in real_sec_struct], dtype=torch.float32)
    return binary_matthews_corrcoef(pred_sec_struct,real_sec_struct).float().mean()



parser = argparse.ArgumentParser(description="Run RDesign with specified test dataset.")
parser.add_argument('--data', type=str, required=True, help='Name of the testing dataset to use')
args = parser.parse_args()

data_choice = args.data

if data_choice.lower() == "rna_puzzles":
    pdb_directory = "data/3D_Puzzles/split_pdbs2/"
elif data_choice.lower() == "casp15":
    pdb_directory = "data/CASP15/split_pdbs2/"
elif data_choice.lower() == "das_split":
    pdb_directory = "data/DAS_Split/"
elif data_choice.lower() == "structsim_v2_split":
    pdb_directory = "data/structsim_v2_split/"

pdb_id_pattern = re.compile(r'PDB_ID:\s+(\S+)')
f1_score_pattern = re.compile(r'F1 Score:\s+([\d\.]+)')
recovery_score_pattern = re.compile(r'Recovery Score:\s+([\d\.]+)')
predicted_sequence_pattern = re.compile(r'Predicted Sequence:\s+(\S+)')
true_sequence_pattern = re.compile(r'True Sequence:\s+(\S+)')
perplexity_pattern = re.compile(r'Perplexity:\s+([\d\.]+)')

results = []

items = os.listdir(pdb_directory)
# Separate into multi-state and single-state design
multi_state_design = [item for item in items if os.path.isdir(os.path.join(pdb_directory, item)) and "Puzzle" in item]
single_state_design = [item for item in items if item.endswith(".pdb") and os.path.isfile(os.path.join(pdb_directory, item))]
for filename in tqdm(single_state_design):
    if filename[-4:] == '.pdb':
        file_path = os.path.join(pdb_directory, filename)

        try:
            command = ['python', 'manul_input/sequence_design_edited.py', file_path]
            result = subprocess.run(command, capture_output=True, text=True)
        except Exception as e:
            print(e)
            
        if pdb_id_pattern.search(result.stdout) == None:
            # print(result.stdout)
            # print(f"Result: \n{result}")
            command = ['python', 'manul_input/sequence_design.py', file_path]
            result = subprocess.run(command, capture_output=True, text=True)
            if pdb_id_pattern.search(result.stdout) == None:
                print(f"Recovery of {file_path} Failed.")
                # print(result)
                one_shot_result = {"pdb_file": str(pdb_directory), "PDB_ID": str(file_path), "F1": 0, "Recovery": 0,
                                   "Pred_Seq": None, "True_Seq": None, "Perplexity": float(0)}
                results.append(one_shot_result)
                continue
            else:
                print(f"Recovery of {file_path} Succeeded.")
        else:
            print(f"{file_path} Succeeded")
        pdb_id = str(pdb_id_pattern.search(result.stdout).group(1))
        f1_score = float(f1_score_pattern.search(result.stdout).group(1))
        recovery_score = float(recovery_score_pattern.search(result.stdout).group(1))
        predicted_sequence = str(predicted_sequence_pattern.search(result.stdout).group(1))
        true_sequence = str(true_sequence_pattern.search(result.stdout).group(1))
        perplexity = float(perplexity_pattern.search(result.stdout).group(1))
        ssc = fetch_secondary_struct_consistency(true_sequence, predicted_sequence).item()
        
        print(f"PDB_ID: {pdb_id}")
        print(f"F1 Score: {f1_score}")
        print(f"Recovery Score: {recovery_score}")
        print(f"Predicted Sequence: {predicted_sequence}")
        print(f"True Sequence: {true_sequence}")
        print(f"Perplexity: {perplexity}")
        print(f"SSC: {ssc}")

        one_shot_result = {"pdb_file": pdb_directory, "PDB_ID": pdb_id, "F1": f1_score,
                           "Recovery": recovery_score, "Pred_Seq":
                           predicted_sequence, "True_Seq": true_sequence, "Perplexity": perplexity,
                           "SSC": ssc}
        
        results.append(one_shot_result)

current_time = datetime.now()
time_string = current_time.strftime("%Y-%m-%d %H:%M:%S")

with open("rdesign_result"+data_choice+time_string+".json", "w") as file:
    for result in results:
        print(type(result))
        print(result)
        file.write(json.dumps(result)+"\n")
print("Results updated!")