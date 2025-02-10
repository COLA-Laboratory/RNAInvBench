# Run any of the Tertiary Folding Algorithms using this file and the tertiary_models_env
import numpy as np
import subprocess
import os
import torch
import re
import json
from tqdm import tqdm
import argparse
from datetime import datetime
from gRNAde.geometric_rna_design_main.gRNAde import gRNAde
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from torchmetrics.functional.classification import binary_matthews_corrcoef



def run_grnade_single_state(g, pdb_filepath=None, directory_filepath=None, split="all", output_filepath=None,
               n_samples=15, temperature=1, partial_seq=None, seed=42, max_num_conformers=1,
               gpu_id=0):
    
    if pdb_filepath is not None:
        sequences, samples, perplexities, recoveries, sc_scores = g.design_from_pdb_file(
            pdb_filepath=pdb_filepath,
            output_filepath=output_filepath,
            n_samples=n_samples,
            temperature=temperature,
            partial_seq=partial_seq,
            seed=seed
        )

    return sequences, samples, perplexities, recoveries, sc_scores


def run_grnade_multi_state(g, pdb_filepath=None, directory_filepath=None, split="all", output_filepath=None,
               n_samples=15, temperature=1, partial_seq=None, seed=42, max_num_conformers=5,
               gpu_id=0):    
    if directory_filepath is not None:
        sequences, samples, perplexities, recoveries, sc_scores = g.design_from_directory(
            directory_filepath=directory_filepath,
            output_filepath=output_filepath,
            n_samples=n_samples,
            temperature=temperature,
            partial_seq=partial_seq,
            seed=seed
        )

    return sequences, samples, perplexities, recoveries, sc_scores


def run_grnade(pdb_file_path, num_conformers=1, state_design="single", 
               model_choice="all", output_file_name="gRNAdeResults"):
    items = os.listdir(pdb_directory)
    # items = [item for item in items if "cleaned" in item]
    # Separate into multi-state and single-state design
    multi_state_design = [item for item in items if os.path.isdir(os.path.join(pdb_directory, item))]
    single_state_design = [item for item in items if item.endswith(".pdb") and os.path.isfile(os.path.join(pdb_directory, item))]
    pdb_files = []
    seqs = []
    samples = []
    perplexities = []
    recoveries = []
    sc_scores = []
    
    if state_design == "single":
        
        print("----Begin Single State Design----")
        if model_choice.lower() == "das":
            g = gRNAde(
                split="das",
                max_num_conformers=1, 
                gpu_id=1
            )
        else:
            g = gRNAde(
                split="all",
                max_num_conformers=1, 
                gpu_id=1
            )
        # single_state_design = ["8hb8_A.pdb"]
        for file in tqdm(single_state_design):
            pdb_files.append(file)
            try:
                seq, sample, perplexity, recovery, sc_score = run_grnade_single_state(g, pdb_filepath=pdb_directory+file,
                                                                                      n_samples=3, gpu_id=1)
                # Find average result
                perplexities.append(np.mean(perplexity))
                recoveries.append(np.mean(recovery))
                sc_scores.append(np.mean(sc_score))
            except Exception as e:
                print(f"Error: Problem with input file {file}")
                print(e)
                #break
        print(f"Mean Perplexity: {np.mean(perplexities)}")
        print(f"Mean Recovery: {np.mean(recoveries)}")
        print(f"Mean SC Score: {np.mean(sc_scores)}")
    
    else:
        
        print("----Begin Multi State Design----")
        multi_state_design = [f for f in multi_state_design if f not in single_state_design and "checkpoint" not in f]
        if model_choice.lower() == "multi":
            g = gRNAde(
                split="multi",
                max_num_conformers=num_conformers, 
                gpu_id=1
            )
        else:
            g = gRNAde(
                split="all",
                max_num_conformers=num_conformers, 
                gpu_id=gpu_id
            )
        for folder in tqdm(multi_state_design):
            pdb_files.append(os.listdir(pdb_directory+folder))
            if not any(os.path.isfile(os.path.join(pdb_directory+folder, f)) for f in os.listdir(pdb_directory+folder)):
                continue
    
            print(f"Cluster: {folder}")
            seq, sample, perplexity, recovery, sc_score = run_grnade_multi_state(g, directory_filepath=pdb_directory+folder,
                                                                                 max_num_conformers=5, 
                                                                                 n_samples=3, gpu_id=1)
            # Find average result
            perplexities.append(np.mean(perplexity))
            recoveries.append(np.mean(recovery))
            print(sc_score)
            flat_sc_score = [item for sublist in sc_score for item in (sublist if isinstance(sublist, list) else [sublist])]
            sc_scores.append(np.mean(flat_sc_score))
        if len(multi_state_design) >= 1:
            print(f"Mean Perplexity: {np.mean(perplexities)}")
            print(f"Mean Recovery: {np.mean(recoveries)}")
            print(f"Mean SC Score: {np.mean(sc_scores)}")

    
    current_time = datetime.now()
    time_string = current_time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"Writing to JSON file...")
    
    current_time = datetime.now()
    time_string = current_time.strftime("%Y-%m-%d %H:%M:%S")
    with open(output_file_name+time_string+".json", "w") as file:
        for pdb_file, recovery, perplexity, sc_score in zip(pdb_files, recoveries, perplexities, sc_scores):
            output = {"pdb_file":pdb_file, "recovery":recovery,
                      "perplexity":perplexity, "sc_score":sc_score}
            file.write(json.dumps(output)+"\n")

    

def run_ribodiffusion(pdb_file_path, output_file_name="RiboDiffusionResults", model_choice="all"):
    script_path = "RiboDiffusion/RiboDiffusion-main/main.py"
    test_results = []
    pdb_files = [f for f in os.listdir(pdb_file_path) if f[-4:] == ".pdb" and "checkpoint" not in f]
    # pdb_files = [f for f in os.listdir("/home/jack/Tertiary_Models/gRNAde/geometric-rna-design-main/data/CASP15/split_pdbs2") if f[-4:] == ".pdb"]
    DATA_PATH = pdb_file_path
    test_list = []
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
                with open("RiboDiffusion/RiboDiffusion-main/exp_inf/fasta/"+pdb_file[:-4]+"_0.fasta") as file:
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
    
    current_time = datetime.now()
    time_string = current_time.strftime("%Y-%m-%d %H:%M:%S")
    with open(output_file_name+time_string+".json", "w") as file:
        for result, files, seq in zip(test_results, test_list, seqs):
            output = {"pdb_file":os.path.join(DATA_PATH, pdb_file), "result":result, "seq":seq}
            file.write(json.dumps(output)+"\n")


def predict_sec_struct(seq):
    eternafold_path = "gRNAde/geometric_rna_design_main/tools/EternaFold/EternaFold-master 2/src/contrafold"
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


def run_rdesign(pdb_directory, output_file_name="RiboDiffusionResults", model_choice="all"):
    pdb_id_pattern = re.compile(r'PDB_ID:\s+(\S+)')
    f1_score_pattern = re.compile(r'F1 Score:\s+([\d\.]+)')
    recovery_score_pattern = re.compile(r'Recovery Score:\s+([\d\.]+)')
    predicted_sequence_pattern = re.compile(r'Predicted Sequence:\s+(\S+)')
    true_sequence_pattern = re.compile(r'True Sequence:\s+(\S+)')
    perplexity_pattern = re.compile(r'Perplexity:\s+([\d\.]+)')
    
    results = []
    
    items = os.listdir(pdb_directory)
    
    single_state_design = [item for item in items if item.endswith(".pdb") and os.path.isfile(os.path.join(pdb_directory, item))]
    for filename in tqdm(single_state_design):
        if filename[-4:] == '.pdb':
            file_path = os.path.join(pdb_directory, filename)
    
            try:
                command = ['python', 'RDesign/R3Design-master/R3Design/manul_input/sequence_design_edited.py', file_path]
                result = subprocess.run(command, capture_output=True, text=True)
            except Exception as e:
                print(e)
                
            if pdb_id_pattern.search(result.stdout) == None:
                # print(result.stdout)
                # print(f"Result: \n{result}")
                command = ['python', 'RDesign/R3Design-master/R3Design/manul_input/sequence_design.py', file_path]
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
    
    with open(output_file_name+time_string+".json", "w") as file:
        for result in results:
            print(type(result))
            print(result)
            file.write(json.dumps(result)+"\n")
    print("Results updated!")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run gRNAde with specified test dataset.")
    parser.add_argument('--data', type=str, required=True, help='Name of the testing dataset to use')
    parser.add_argument('--model', type=str, required=True, help='Name of the model to run')
    parser.add_argument('--model_train', type=str, required=False, default="all", help="Name of the trained model to use")
    parser.add_argument('--num_conformers', type=int, required=False, default=1, help='Name of the testing dataset to use')
    args = parser.parse_args()
    model_choice = args.model
    data_choice = args.data
    model_train_choice = args.model_train
    num_conformers = args.num_conformers

    if data_choice.lower() == "rna_puzzles":
        pdb_directory = "data/3D_Puzzles/split_pdbs/"
    elif data_choice.lower() == "casp15":
        pdb_directory = "data/CASP15/split_pdbs/"
    elif data_choice.lower() == "das_split":
        pdb_directory = "data/DAS_Split/"
    elif data_choice.lower() == "structsim_v2_split":
        pdb_directory = "data/structsim_v2_split/"

    if model_choice.lower() == "grnade":
        run_grnade(pdb_directory, num_conformers=num_conformers,
                   model_choice=model_train_choice)
    elif model_choice.lower() == "ribodiffusion":
        run_ribodiffusion(pdb_directory, model_choice=model_train_choice)
    elif model_choice.lower() == "rdesign":
        run_rdesign(pdb_directory, model_choice=model_train_choice)
