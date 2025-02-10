import numpy as np
import subprocess
import os
import re
import argparse
from datetime import datetime
import json
from gRNAde import gRNAde
from tqdm import tqdm



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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run gRNAde with specified test dataset.")
    parser.add_argument('--data', type=str, required=True, help='Name of the testing dataset to use')
    parser.add_argument('--num_conformers', type=int, required=False, default=1, help='Name of the testing dataset to use')
    args = parser.parse_args()
    
    data_choice = args.data
    num_conformers = args.num_conformers
    
    items = os.listdir(pdb_directory)
    # items = [item for item in items if "cleaned" in item]
    # Separate into multi-state and single-state design
    multi_state_design = [item for item in items if os.path.isdir(os.path.join(pdb_directory, item))]
    single_state_design = [item for item in items if item.endswith(".pdb") and os.path.isfile(os.path.join(pdb_directory, item))]
    seqs = []
    samples = []
    perplexities = []
    recoveries = []
    sc_scores = []
    print("----Begin Single State Design----")
    if data_choice.lower() == "das_split" or data_choice.lower() == "rna_puzzles":
        g = gRNAde(
            split="das",
            max_num_conformers=1, 
            gpu_id=1
        )
    # single_state_design = ["8hb8_A.pdb"]
    for file in tqdm(single_state_design):
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

    print("----Begin Multi State Design----")
    multi_state_design = [f for f in multi_state_design if f not in single_state_design and "checkpoint" not in f]
    if data_choice.lower() == "structsim_v2_split":
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