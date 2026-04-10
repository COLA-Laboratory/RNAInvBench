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
from gRNAde.geometric_rna_design_main.src.evaluator import self_consistency_score_eternafold, \
    self_consistency_score_rhofold
from gRNAde.geometric_rna_design_main.src.data.data_utils import pdb_to_tensor, get_c4p_coords
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from torchmetrics.functional.classification import binary_matthews_corrcoef
import ViennaRNA
import gc
from fetch_contact_map import get_contact_map
import sys
from pathlib import Path
import pandas as pd
from collections import Counter


def run_dssr(pdb_path: Path, json_output: bool = False) -> str:
    """
    Runs x3dna-dssr on the given PDB file with default settings.

    Args:
        pdb_path (Path): Path to the input PDB file.
        json_output (bool): Whether to request JSON output via --json.

    Returns:
        str: The standard output from DSSR (text or JSON).
    """
    if not pdb_path.is_file():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    # Build the command
    cmd = ["/home/jc1417/RNAInvBench/Tertiary_Design/x3dna-dssr-2.5.4/x3dna-dssr", f"-i={pdb_path}"]
    if json_output:
        cmd.append("--json")

    # Execute DSSR
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True
    )
    # Cleanup the files
    subprocess.run(["/home/jc1417/RNAInvBench/Tertiary_Design/x3dna-dssr-2.5.4/x3dna-dssr", "--cleanup"])
    
    return result.stdout


def fetch_structural_features(res):
    lines = res.splitlines()
    
    pattern = re.compile(
        r'^\s*(\d+)\s+'         # Serial number
        r'([ACUG])\s+'          # Base letter
        r'([^\s]+)\s+'          # DBN (e.g., . or complex symbol)
        r'([^\s]+)\s+'          # ID string (e.g., A.A1)
        r'([0-9.]+|---)\s+'     # RMSD (number or '---')
        r'(.*)$'                # Features list
    )
    
    data = []
    start_parsing = False
    
    for line in lines:
        match = pattern.match(line)
        if match:
            start_parsing = True  # Enable parsing from the first match
        if start_parsing and match:
            serial = int(match.group(1))
            base = match.group(2)
            dbn = match.group(3)
            id_str = match.group(4)
            rmsd_str = match.group(5)
            features = [f.strip() for f in match.group(6).split(',') if f.strip()]
            rmsd = float(rmsd_str) if rmsd_str != '---' else None
    
            data.append({
                'Serial': serial,
                'Base': base,
                'DBN': dbn,
                'ID': id_str,
                'RMSD': rmsd,
                'Features': features
            })
    
        
    df = pd.DataFrame(data)
    # print(df.head())
    features = (df["Features"].to_list())
    flat_features = [feature for sublist in features for feature in sublist]
    not_useful_features = ["anti", "canonical", "non-pair-contact"]
    # not_useful_features = []
    secondary_motif_features = ["helix", "stem", "hairpin-loop", "internal-loop", "bulge", "junction-loop",
                        "stem-end", "helix-end", "multi-loop"]
    # all other motifs are tertiary
    counts = Counter(flat_features)
    
    secondary_motifs = {}    
    tertiary_motifs = {}
    
    for feature, count in counts.items():
        if feature in not_useful_features:
            continue
        elif feature in secondary_motif_features:
            secondary_motifs[feature] = count
            continue
        else:
            tertiary_motifs[feature] = count
    return secondary_motifs, tertiary_motifs


def get_gddt(true_pdb, pred_pdb):
    """Global Distance Deviation Test metric (GDDT).

    Credit: Arian Jamasb, graphein (https://github.com/a-r-j/graphein)

    https://en.wikipedia.org/wiki/Global_distance_test

    The GDT score is calculated as the largest set of amino acid residues'
    alpha carbon atoms in the model structure falling within a defined
    distance cutoff of their position in the experimental structure, after
    iteratively superimposing the two structures. By the original design the
    GDT algorithm calculates 20 GDT scores, i.e. for each of 20 consecutive distance
    cutoffs (``0.5 Å, 1.0 Å, 1.5 Å, ... 10.0 Å``). For structure similarity assessment
    it is intended to use the GDT scores from several cutoff distances, and scores
    generally increase with increasing cutoff. A plateau in this increase may
    indicate an extreme divergence between the experimental and predicted structures,
    such that no additional atoms are included in any cutoff of a reasonable distance.
    The conventional GDT_TS total score in CASP is the average result of cutoffs at
    ``1``, ``2``, ``4``, and ``8`` Å.

    Random predictions give around 20; getting the gross topology right gets one to ~50; 
    accurate topology is usually around 70; and when all the little bits and pieces, 
    including side-chain conformations, are correct, GDT_TS begins to climb above 90.

    Credit: Chaitanya Joshi (https://github.com/chaitjo/geometric-rna-design)

    """
    # Use rotation matrix to align the structures together
    print(pred_pdb)
    print(true_pdb)
    _, y_hat, _, _ = pdb_to_tensor(pred_pdb, return_sec_struct=False, return_sasa=False, keep_insertions=False)
    _, y, _, _ = pdb_to_tensor(true_pdb, return_sec_struct=False, return_sasa=False, keep_insertions=False)
    print(y_hat.shape)
    print(y.shape)
    true_coords = get_c4p_coords(y)
    true_coords = true_coords - true_coords.mean(dim=0)
    
    from MDAnalysis.analysis.align import rotation_matrix
    pred = get_c4p_coords(y_hat)
    pred = pred - pred.mean(dim=0)

    R_hat = rotation_matrix(
        pred,  # mobile set
        true_coords # reference set
    )[0]

    pred_coords = pred @ R_hat.T

    # Get distance between points
    # print(f"Pred Shape {pred_coords}")
    # print(f"True Shape {true_coords}")
    dist = torch.norm(true_coords - pred_coords, dim=1)
    # print(f"Calculated dist: {dist}")

    # Return mean fraction of distances below cutoff for each cutoff (1, 2, 4, 8)
    count_1 = (dist < 1).sum() / dist.numel()
    count_2 = (dist < 2).sum() / dist.numel()
    count_4 = (dist < 4).sum() / dist.numel()
    count_8 = (dist < 8).sum() / dist.numel()
    out = torch.mean(torch.tensor([count_1, count_2, count_4, count_8]))
    if torch.isnan(out):
        return torch.tensor(0.0)
    return out.item()



def align_pdb_fetch_tm_rmsd(gt_path, pred_path):
    """
    Aligns two PDB files using USalign and returns the TM-score.

    Parameters:
        gt_file (str): Path to ground truth/reference PDB.
        compared_file (str): Path to predicted/aligned PDB.
        usalign_path (str): Path to USalign binary.

    Returns:
        float: TM-score if successful, -1 if failed.
    """
    # try:
    print('/home/jc1417/RNAInvBench/Tertiary_Design/gRNAde/geometric_rna_design_main/tools/USalign/USalign',
        gt_path,
        pred_path,
        '-outfmt', '2')
    subprocess.run([
        '/home/jc1417/RNAInvBench/Tertiary_Design/gRNAde/geometric_rna_design_main/tools/USalign/USalign',
        gt_path,
        pred_path,
        '-outfmt', '2'
    ], stdout=open('./tmscore.txt', 'w'), check=True)

    with open('./tmscore.txt', 'r') as f:
        data = f.readlines()

    header = data[0].strip().split('\t')
    values = data[1].strip().split('\t')
    record = dict(zip(header, values))

    # print(record)  # Debug print
    tm_score = (float(record["TM1"]) + float(record["TM2"])) / 2
    rmsd = float(record["RMSD"])
    gddt = get_gddt(gt_path, pred_path)
        
    # except Exception as e:
    #     print(e)
    #     print(pred_path.split('/')[-1])
    #     return -1
    return tm_score, rmsd, gddt


def get_motifs_differences(predicted_seq, pdb_file, use_relax=False):
    predict_folder = "/home/jc1417/RNAInvBench/Tertiary_Design/fastas_temp/"
    predict_file = "/home/jc1417/RNAInvBench/Tertiary_Design/fastas_temp/design.pdb"
    # RhoFold already predicted the structure, just get it and analyse the motifs
    res = run_dssr(Path(predict_file))
    predicted_secondary_motifs, predicted_tertiary_motifs = fetch_structural_features(res)
    # Get Real Seq Secondary and Tertiary Motifs
    res = run_dssr(Path(pdb_file))
    ground_truth_secondary_motifs, ground_truth_tertiary_motifs = fetch_structural_features(res)
    secondary_motif_difference = {}
    for motif in set(ground_truth_secondary_motifs) | set(predicted_secondary_motifs):
        secondary_motif_difference[motif] = (
            ground_truth_secondary_motifs.get(motif, 0)
            - predicted_secondary_motifs.get(motif, 0)
        )
    
    tertiary_motif_difference = {}
    for motif in set(ground_truth_tertiary_motifs) | set(predicted_tertiary_motifs):
        tertiary_motif_difference[motif] = (
            ground_truth_tertiary_motifs.get(motif, 0)
            - predicted_tertiary_motifs.get(motif, 0)
        )
    return secondary_motif_difference, tertiary_motif_difference


def rhofold_predict(real_seq, pred_seq, output_dir, use_relax=False):
    from gRNAde.geometric_rna_design_main.tools.rhofold.rf import RhoFold
    from gRNAde.geometric_rna_design_main.tools.rhofold.config import rhofold_config
    # return - If we want to skip alphafold
    device = "cuda:1"
    rhofold = None
    # remove output file
    if os.path.isfile(os.path.join(output_dir, f"true_design.pdb")):
        os.remove(os.path.join(output_dir, f"true_design.pdb"))
    if os.path.isfile(os.path.join(output_dir, f"pred_design.pdb")):
        os.remove(os.path.join(output_dir, f"pred_design.pdb"))
    tries = 0
    while (True):
        rhofold = RhoFold(rhofold_config, device)
        rhofold_path = "/home/jc1417/RNAInvBench/Tertiary_Design/gRNAde/geometric_rna_design_main/tools/rhofold/model_20221010_params.pt"
        print(f"Loading RhoFold checkpoint: {rhofold_path}")
        rhofold.load_state_dict(torch.load(rhofold_path, map_location=torch.device('cpu'))['model'])

        rhofold = rhofold.to(device)
        rhofold.eval()
        design_fasta_path = os.path.join(output_dir, f"design.fasta")
        SeqIO.write(SeqRecord(Seq(real_seq), id="real_sequence", description=""), design_fasta_path, "fasta")
        design_pdb_path = os.path.join(output_dir, f"true_design.pdb")
        rhofold.predict(design_fasta_path, design_pdb_path, use_relax)
        design_fasta_path = os.path.join(output_dir, f"design.fasta")
        SeqIO.write(pred_seq, design_fasta_path, "fasta")

        design_pdb_path = os.path.join(output_dir, f"pred_design.pdb")
        
        rhofold.predict(design_fasta_path, design_pdb_path, use_relax)
        if not os.path.exists(design_pdb_path) and tries < 1:
            print("RhoFold failed to predict structure, trying again...")
            continue
            
        # except Exception as e:
        #     print("RhoFold failed to predict structure, trying again...")
        #     print(e)
        #     continue
            
        del rhofold
        torch.cuda.empty_cache()
        gc.collect()
        break


def get_tertiary_structure_consistency_tm_and_motifs(predicted_seq, real_seq, pdb_file, get_motifs=False):
    predict_file = "/home/jc1417/RNAInvBench/Tertiary_Design/fastas_temp/"
    rhofold_true_pdb_file = "/home/jc1417/RNAInvBench/Tertiary_Design/fastas_temp/"+"true_design.pdb"
    rhofold_predict_pdb_file = "/home/jc1417/RNAInvBench/Tertiary_Design/fastas_temp/"+"pred_design.pdb"
    rhofold_predict(real_seq, predicted_seq, predict_file)
    tm, rmsd, gddt = align_pdb_fetch_tm_rmsd(rhofold_true_pdb_file, rhofold_predict_pdb_file)
    # tm = None
    # rmsd = None
    # gddt = None
    if get_motifs:
        try:
            secondary_motif_difference, tertiary_motif_difference = get_motifs_differences(predicted_seq, pdb_file)
        except:
            secondary_motif_difference = None
            tertiary_motif_difference = None
        return tm, rmsd, gddt, secondary_motif_difference, tertiary_motif_difference
    else:
        return tm, rmsd, gddt, None, None


def run_grnade_single_state(g, pdb_filepath=None, directory_filepath=None, split="all", output_filepath=None,
                            n_samples=15, temperature=1, partial_seq=None, seed=42, max_num_conformers=1,
                            gpu_id=0):
    if pdb_filepath is not None:
        print(pdb_filepath)
        sequences, samples, perplexities, recoveries, sc_scores, true_seq = g.design_from_pdb_file(
            pdb_filepath=pdb_filepath,
            output_filepath=output_filepath,
            n_samples=n_samples,
            temperature=temperature,
            partial_seq=partial_seq,
            seed=seed
        )

    return sequences, samples, perplexities, recoveries, sc_scores, true_seq


def run_grnade_multi_state(g, pdb_filepath=None, directory_filepath=None, split="all", output_filepath=None,
                           n_samples=15, temperature=1, partial_seq=None, seed=42, max_num_conformers=5,
                           gpu_id=0):
    if directory_filepath is not None:
        sequences, samples, perplexities, recoveries, sc_scores, true_seq = g.design_from_directory(
            directory_filepath=directory_filepath,
            output_filepath=output_filepath,
            n_samples=n_samples,
            temperature=temperature,
            partial_seq=partial_seq,
            seed=seed
        )

    return sequences, samples, perplexities, recoveries, sc_scores, true_seq


def get_recovery(seqrecord):
    desc = seqrecord.description
    for part in desc.split(','):
        if 'recovery=' in part:
            return float(part.split('=')[1])
    return 0.0  # fallback if missing or malformed


def run_grnade(pdb_directory, num_conformers=1, state_design="single",
               model_choice="all", output_file_name="gRNAdeResults", get_motifs=False):
    import random
    seed = random.randint(1, 10000)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    items = os.listdir(pdb_directory)
    pdb_directory += "/"
    # items = [item for item in items if "cleaned" in item]
    # Separate into multi-state and single-state design
    multi_state_design = [item for item in items if os.path.isdir(os.path.join(pdb_directory, item))]
    single_state_design = [item for item in items if
                           item.endswith(".pdb") and os.path.isfile(os.path.join(pdb_directory, item))]
    pdb_files = []
    seqs = []
    true_seqs = []
    samples = []
    perplexities = []
    recoveries = []
    sc_scores = []
    rmsd_scores = []
    tm_scores = []
    gdt_scores = []
    secondary_differences = []
    tertiary_differences = []

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
                seq, sample, perplexity, recovery, sc_score, true_seq = run_grnade_single_state(g, pdb_filepath=pdb_directory + file,
                                                                                                n_samples=3, gpu_id=1)
                best_seq = max(seq, key=get_recovery)
                print(f"Best Seq {str(best_seq.seq)}\nBest Seq Len {len(str(best_seq.seq))}")
                print(f"Pred Seq Len {len(true_seq)}")
                # Find average result
                perplexities.append(np.mean(perplexity))
                recoveries.append(np.mean(recovery))
                sc_scores.append(np.mean(sc_score))
                print(f"Predicted Seq: {seq}")
                print(f"Real Seq: {true_seq}")
                tm, tsc, gdt, secondary_difference, tertiary_difference = get_tertiary_structure_consistency_tm_and_motifs(
                best_seq, true_seq, pdb_directory + file, get_motifs=get_motifs)
                seqs.append(best_seq)
                true_seqs.append(true_seq)
                rmsd_scores.append(tsc)
                tm_scores.append(tm)
                gdt_scores.append(gdt)
                print("rmsd: ", str(tsc))
                print("tm: ", str(tm))
                print("gdt: ", str(gdt))
                if get_motifs:
                    print(f"Secondary Difference: {secondary_difference}")
                    secondary_differences.append(secondary_difference)
                    print(f"Tertiary Difference: {tertiary_difference}")
                    tertiary_differences.append(tertiary_difference)

            except Exception as e:
                print(f"Error: Problem with input file {file}")
                continue

            # except Exception as e:
            #     print(f"Error: Problem with input file {file}")
            #     print(e)
            # break
        print(f"Mean Perplexity: {np.mean(perplexities)}")
        print(f"Mean Recovery: {np.mean(recoveries)}")
        print(f"Mean SC Score: {np.mean(sc_scores)}")
        #print(f"Mean TM Score: {np.mean(tm_scores)}")
        #print(f"Mean RMSD Score: {np.mean(rmsd_scores)}")
        #print(f"Mean GDT Score: {np.mean(gdt_scores)}")

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
            pdb_files.append(os.listdir(pdb_directory + folder))
            if not any(os.path.isfile(os.path.join(pdb_directory + folder, f)) for f in
                       os.listdir(pdb_directory + folder)):
                continue

            print(f"Cluster: {folder}")
            seq, sample, perplexity, recovery, sc_score = run_grnade_multi_state(g,
                                                                                 directory_filepath=pdb_directory + folder,
                                                                                 max_num_conformers=5,
                                                                                 n_samples=3, gpu_id=1)
            # Find average result
            perplexities.append(np.mean(perplexity))
            recoveries.append(np.mean(recovery))
            print(sc_score)
            flat_sc_score = [item for sublist in sc_score for item in
                             (sublist if isinstance(sublist, list) else [sublist])]
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
    dataset_name = pdb_directory[23:-1]
    out_dir = Path(output_file_name)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(output_file_name + dataset_name + time_string + ".json", "w") as file:
        for pdb_file, seq, true_seq, recovery, perplexity, sc_score, rmsd, tm_score, gdt_score in zip(pdb_files, seqs, true_seqs, recoveries, perplexities, sc_scores, rmsd_scores, tm_scores, gdt_scores):
        # for pdb_file, recovery, perplexity, sc_score in zip(pdb_files, recoveries, perplexities, sc_scores):
            output = {"pdb_file": str(pdb_file), "Pred_Seq": str(seq.seq), "True_Seq": true_seq, "Recovery": str(recovery),
                      "Perplexity": str(perplexity), "SSC": str(sc_score), "RMSD": rmsd,
                      "TM_Score": tm_score, "GDT_Score": gdt_score}
            print(output)
            file.write(json.dumps(output) + "\n")



def run_ribodiffusion(pdb_directory, output_file_name="RiboDiffusionResults", model_choice="all", get_motifs=False):
    script_path = "RiboDiffusion/RiboDiffusion-main/main.py"
    test_results = []
    pdb_files = [f for f in os.listdir(pdb_directory) if f[-4:] == ".pdb" and "checkpoint" not in f]
    # pdb_files = [f for f in os.listdir("/home/jack/Tertiary_Models/gRNAde/geometric-rna-design-main/data/CASP15/split_pdbs2") if f[-4:] == ".pdb"]
    DATA_PATH = pdb_directory
    test_list = []
    seqs = []
    true_pdb_files = []
    i = 0
    for pdb_files_list in tqdm(pdb_files):
        found_seq = ""
        t_pdb_files = [pdb_files_list]
        single_test_rec = []
        single_test_per = []
        single_test_ssc = []
        tm_scores = []
        rmsd_scores = []
        gdt_scores = []
        secondary_differences = []
        tertiary_differences = []
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
                true_seq_match = re.search(r"true sequence\s+([ACGU]+)", output)
                true_seq = str(true_seq_match.group(1)) if true_seq_match else None
                if recovery_rate != None:
                    single_test_rec.append(recovery_rate)
                if perplexity != None:
                    single_test_per.append(perplexity)
                if ssc != None:
                    single_test_ssc.append(ssc)
                # get sequence from file
                with open("RiboDiffusion/RiboDiffusion-main/exp_inf/fasta/" + pdb_file[:-4] + "_0.fasta") as file:
                    for line in file.readlines():
                        if line[0] == ">":
                            continue
                        else:
                            seqs.append(line.strip())
                            break
                print("RiboDiffusion Pred Seq: ", seqs[-1])
                # try:
                tm, tsc, gdt, secondary_difference, tertiary_difference = get_tertiary_structure_consistency_tm_and_motifs(
                    SeqRecord(Seq(seqs[-1]), id="predicted_sequence", description=""), true_seq, pdb_file_path, get_motifs=get_motifs)
        
                rmsd_scores.append(tsc)
                tm_scores.append(tm)
                gdt_scores.append(gdt)
                secondary_differences.append(secondary_difference)
                tertiary_differences.append(tertiary_difference)
                # except Exception as e:
                #     print(e)
                #     print("predicting tertiary structure failed")

            except subprocess.CalledProcessError as e:
                print(f"Error processing {pdb_file}: {e}")
        if len(single_test_rec) > 0:
            print(f"Test {i} Recovery: {np.mean(single_test_rec)}")
            print(f"Test {i} Perplexity: {np.mean(single_test_per)}")
            print(f"Test {i} SSC: {np.mean(single_test_ssc)}")
            #print(f"Mean TM Score: {np.mean(tm_scores)}")
            #print(f"Mean RMSD Score: {np.mean(rmsd_scores)}")
            #print(f"Mean GDT Score: {np.mean(gdt_scores)}")
            test_results.append({"Recovery": np.mean(single_test_rec), "Pred_Seq": seqs[-1],
                                 "True_Seq": true_seq, "Perplexity": np.mean(single_test_per),
                                 "SSC": np.mean(single_test_ssc), "RMSD": np.mean(rmsd_scores),
                                 "TM_Score": np.mean(tm_scores), "GDT_Score": np.mean(gdt_scores),
                                 "Secondary_Differences": secondary_differences,
                                 "Tertiary_Differences": tertiary_differences})
            test_list.append(pdb_file)
        i += 1
    print(f"Overall Recovery: {np.mean([x['Recovery'] for x in test_results])}")
    print(f"Overall Perplexity: {np.mean([x['Perplexity'] for x in test_results])}")
    print(f"Overall SSC: {np.mean([x['SSC'] for x in test_results])}")
    # print(f"Overall RMSD: {np.mean([x['RMSD'] for x in test_results])}")
    # print(f"Overall TM Score: {np.mean([x['TM_Score'] for x in test_results])}")
    print(f"Writing to JSON file...")

    current_time = datetime.now()
    time_string = datetime.now().strftime("%Y%m%d-%H%M%S")
    dataset_name = pdb_directory[23:-1]
    out_dir = Path(output_file_name)
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(output_file_name + dataset_name + time_string + ".json", "w") as file:
        for result, files, seq, pdb_file in zip(test_results, test_list, seqs, pdb_files):
            output = {"pdb_file": os.path.join(DATA_PATH, pdb_file), "result": result, "seq": seq}
            file.write(json.dumps(output) + "\n")
    print("Results updated!")



def predict_sec_struct_eternafold(seq):
    eternafold_path = "/home/jc1417/RNAInvBench/Tertiary_Design/gRNAde/geometric_rna_design_main/tools/EternaFold/EternaFold-master 2/src/contrafold"
    current_datetime = datetime.now().strftime("%Y%m%d_%H%M%S")
    fasta_file_path = f"fastas_temp/temp_{current_datetime}.fasta"
    cmd = [eternafold_path, "predict", fasta_file_path, "--params", "parameters/EternaFoldParams.v1"]
    SeqIO.write(SeqRecord(Seq(seq), id="temp"), fasta_file_path, "fasta")
    output = subprocess.run(cmd, check=True, capture_output=True).stdout.decode("utf-8")
    return [x for x in output.split("\n")[-2]]

def predict_sec_struct(seq):
    struct, mfe = ViennaRNA.fold(seq)
    return [x for x in struct]


def fetch_secondary_struct_consistency(real_seq, pred_seq):
    real_sec_struct = predict_sec_struct_eternafold(real_seq)
    pred_sec_struct = predict_sec_struct_eternafold(pred_seq)
    pred_sec_struct = torch.tensor([1 if c in "()" else 0 for c in pred_sec_struct], dtype=torch.float32)
    real_sec_struct = torch.tensor([1 if c in "()" else 0 for c in real_sec_struct], dtype=torch.float32)
    return binary_matthews_corrcoef(pred_sec_struct, real_sec_struct).float().mean()


def run_rdesign(pdb_directory, output_file_name="RDesignResults", model_choice="all", get_motifs=False):
    pdb_id_pattern = re.compile(r'PDB_ID:\s+(\S+)')
    f1_score_pattern = re.compile(r'F1 Score:\s+([\d\.]+)')
    recovery_score_pattern = re.compile(r'Recovery Score:\s+([\d\.]+)')
    predicted_sequence_pattern = re.compile(r'Predicted Sequence:\s+(\S+)')
    true_sequence_pattern = re.compile(r'True Sequence:\s+(\S+)')
    perplexity_pattern = re.compile(r'Perplexity:\s+([\d\.]+)')

    results = []

    items = os.listdir(pdb_directory)

    single_state_design = [item for item in items if
                           item.endswith(".pdb") and os.path.isfile(os.path.join(pdb_directory, item))]
    for filename in tqdm(single_state_design):
        if filename[-4:] == '.pdb':
            file_path = os.path.join(pdb_directory, filename)

            try:
                command = ['python', 'RDesign/R3Design-master/R3Design/manul_input/sequence_design_edited.py',
                           file_path]
                result = subprocess.run(command, capture_output=True, text=True)
            except Exception as e:
                print(e)
            print(result)
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
            try:
                print(predicted_sequence)
                print(file_path)
                tm, tsc, gdt, secondary_difference, tertiary_difference = get_tertiary_structure_consistency_tm_and_motifs(
                SeqRecord(Seq(predicted_sequence), id="predicted_sequence", description=""), true_sequence,
                "/home/jc1417/RNAInvBench/Tertiary_Design/" + file_path, get_motifs=get_motifs)
                rmsd_score = tsc
                tm_score = tm
            except Exception as e:
                print("predicting tertiary structure failed")
                print(e)
                rmsd_score = 0
                tm_score = 0
                gdt = 0
                secondary_difference = None
                tertiary_difference = None

            print(f"PDB_ID: {pdb_id}")
            print(f"F1 Score: {f1_score}")
            print(f"Recovery Score: {recovery_score}")
            print(f"Predicted Sequence: {predicted_sequence}")
            print(f"True Sequence: {true_sequence}")
            print(f"Perplexity: {perplexity}")
            print(f"SSC: {ssc}")
            print(f"TM Score: {tm_score}")
            print(f"RMSD Score: {rmsd_score}")
            print(f"GDT Score: {gdt}")
            if get_motifs:
                print(f"Secondary Difference: {secondary_difference}")
                print(f"Tertiary Difference: {tertiary_difference}")

            one_shot_result = {"pdb_file": pdb_directory, "PDB_ID": pdb_id, "F1": f1_score,
                               "Recovery": recovery_score, "Pred_Seq": predicted_sequence,
                               "True_Seq": true_sequence, "Perplexity": perplexity, "SSC": ssc,
                               "RMSD": rmsd_score, "TM_Score": tm_score, "GDT_Score": gdt,
                               "Secondary_Difference": secondary_difference,
                               "Tertiary_Difference": tertiary_difference}

            results.append(one_shot_result)

    current_time = datetime.now()
    time_string = datetime.now().strftime("%Y%m%d-%H%M%S")
    # data/RNA_Type_Benchmark/CRISPR_RNA_only/
    dataset_name = pdb_directory[23:-1]
    out_dir = Path(output_file_name)
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(output_file_name + dataset_name + time_string + ".json", "w") as file:
        for result in results:
            file.write(json.dumps(result) + "\n")
    print("Results updated!")
    # print(f"Failed: {failed}")

def run_rhodesign(pdb_directory, output_file_name="RhoDesignResults", model_choice="all", contact_map=False, get_motifs=False, seed=1):
    print(pdb_directory)
    # pdb_id_pattern = re.compile(r'PDB_ID:\s+(\S+)')
    # f1_score_pattern = re.compile(r'F1 Score:\s+([\d\.]+)')
    recovery_score_pattern = re.compile(r'recovery rate:\s+([\d\.]+)')
    predicted_sequence_pattern = re.compile(r'sequence:\s+(\S+)')
    true_sequence_pattern = re.compile(r'true sequence:\s+(\S+)')
    perplexity_pattern = re.compile(r'perplexity:\s+([\d\.]+)')

    results = []

    items = os.listdir(pdb_directory)

    single_state_design = [item for item in items if
                           (item.endswith(".pdb") or item.endswith(".cif")) and os.path.isfile(os.path.join(pdb_directory, item))]
    failed = 0
    for filename in tqdm(single_state_design):
        if filename[-4:] == '.pdb' or filename[-4:] == ".cif":
            file_path = os.path.join(pdb_directory, filename)
            print(file_path)
            if contact_map:
                print(f"Using Contact Map...")
                try:
                    get_contact_map(file_path)
                    command = ['python', 'RhoDesign/RhoDesign_main/src/inference.py',
                               "-pdb", file_path, "-ss", "/home/jc1417/RNAInvBench/Tertiary_Design/ss_contact_map.npy", "-save", "fastas_temp/", "-temp", "1"]
                    result = subprocess.run(command, capture_output=True, text=True)
                except Exception as e:
                    print(e)
                    continue
            else:
                print(f"No Contact Map Supplied...")
                try:
                    command = ['python', 'RhoDesign/RhoDesign_main/src/inference_without2d.py',
                               "-pdb", file_path, "-save", "fastas_temp/"]
                    result = subprocess.run(command, capture_output=True, text=True)
                except Exception as e:
                    print(e)
                    continue
            print(result)
            try:
                recovery_score = float(recovery_score_pattern.search(result.stdout).group(1))
                predicted_sequence = str(predicted_sequence_pattern.search(result.stdout).group(1))
                perplexity = float(perplexity_pattern.search(result.stdout).group(1))
                true_sequence = str(true_sequence_pattern.search(result.stdout).group(1))
                print("true sequence: ", true_sequence)
                print("pred sequence: ", predicted_sequence)
                ssc = fetch_secondary_struct_consistency(true_sequence, predicted_sequence).item()
                tm, tsc, gdt, secondary_difference, tertiary_difference = get_tertiary_structure_consistency_tm_and_motifs(
                    SeqRecord(Seq(predicted_sequence), id="predicted_sequence", description=""), true_sequence,
                    "/home/jc1417/RNAInvBench/Tertiary_Design/" + file_path, get_motifs=get_motifs)
            except Exception as e:
                print(e)
                continue
            rmsd_score = tsc
            tm_score = tm
            # except Exception as e:
            #     print(e)
            #     failed += 1
            #     continue

            print(f"PDB: {filename[:-4]}")
            print(f"Recovery Score: {recovery_score}")
            print(f"Predicted Sequence: {predicted_sequence}")
            print(f"Perplexity: {perplexity}")
            print(f"SSC: {ssc}")
            print(f"TM Score: {tm_score}")
            print(f"RMSD Score: {rmsd_score}")
            # print(f"GDT Score: {gdt}")
            # print(f"Secondary Difference: {secondary_difference}")
            # print(f"Tertiary Difference: {tertiary_difference}")

            one_shot_result = {"pdb_file": pdb_directory, "PDB_ID": filename[:-4],
                               "Recovery": recovery_score, "Pred_Seq": predicted_sequence,
                               "True_Seq": true_sequence, "Perplexity": perplexity,
                               "SSC": ssc, "RMSD": rmsd_score, "TM_Score": tm_score, "GDT_Score": gdt}
                               # "Secondary_Differences": secondary_difference,
                               # "Tertiary_Differences": tertiary_difference}

            results.append(one_shot_result)

    current_time = datetime.now()
    time_string = datetime.now().strftime("%Y%m%d-%H%M%S")
    # data/RNA_Type_Benchmark/CRISPR_RNA_only/
    dataset_name = pdb_directory[24:-1]
    out_dir = Path(output_file_name)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(output_file_name + dataset_name + time_string + ".json")
    with open(output_file_name + "/" + ("CM" if contact_map else "No_CM") + "_" + dataset_name + time_string + ".json", "w") as file:
        for result in results:
            file.write(json.dumps(result) + "\n")
    print("Results updated!")
    print(f"Failed: {failed}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run gRNAde with specified test dataset.")
    parser.add_argument('--data', type=str, required=True, help='Name of the testing dataset to use')
    parser.add_argument('--model', type=str, required=True, help='Name of the model to run')
    parser.add_argument('--model_train', type=str, required=False, default="all",
                        help="Name of the trained model to use")
    parser.add_argument('--num_conformers', type=int, required=False, default=1,
                        help='Name of the testing dataset to use')
    parser.add_argument('--contact_map', type=bool, required=False, default=False)
    parser.add_argument('--get_motifs', type=bool, required=False, default=False)
    parser.add_argument('--pdb_directory', type=str, required=False, default=False)
    args = parser.parse_args()
    model_choice = args.model
    data_choice = args.data
    model_train_choice = args.model_train
    num_conformers = args.num_conformers
    contact_map = args.contact_map
    get_motifs = args.get_motifs
    pdb_file_dir = args.pdb_directory

    if data_choice.lower() == "rna_puzzles":
        pdb_directory = "data/3D_Puzzles/split_pdbs/"
    elif data_choice.lower() == "casp15":
        pdb_directory = "data/CASP15/split_pdbs3/"
    elif data_choice.lower() == "das_split":
        pdb_directory = "data/DAS_Split/"
    elif data_choice.lower() == "structsim_v2_split":
        pdb_directory = "data/structsim_v2_split/"
    elif data_choice.lower() == "rhodesign_test":
        pdb_directory = "RhoDesign/RhoDesign_main/data/test/"
    elif data_choice.lower() == "custom_set":
        if model_choice.lower() == "rhodesign":
            pdb_directory = "data/Custom_Puzzle_Set/Final_Set_PDBs2/"
        else:
            pdb_directory = ""
    elif data_choice.lower() == "cif_test":
        pdb_directory = "data/cif_test/"
    elif data_choice.lower() == "lncrna":
        pdb_directory = "data/RNA_Type_Benchmark/lncRNA_split_pdbs"
    elif data_choice.lower() == "crispr":
        pdb_directory = "data/RNA_Type_Benchmark/CRISPR_RNA_split_pdbs"
    elif data_choice.lower() == "regulatory":
        pdb_directory = "data/RNA_Type_Benchmark/Short_Regulatory_ncRNA_split_pdbs"
    elif data_choice.lower() == "housekeeping":
        pdb_directory = "data/RNA_Type_Benchmark/Housekeeping_ncRNA_split_pdbs"
    elif data_choice.lower() == "coding":
        pdb_directory = "data/RNA_Type_Benchmark/coding_RNA_split_pdbs"
    elif pdb_file_dir is not None:
        pdb_directory = pdb_file_dir

    if model_choice.lower() == "grnade":
        run_grnade(pdb_directory, num_conformers=num_conformers,
                   model_choice=model_train_choice, get_motifs=get_motifs)
    elif model_choice.lower() == "ribodiffusion":
        run_ribodiffusion(pdb_directory, model_choice=model_train_choice, get_motifs=get_motifs)
    elif model_choice.lower() == "rdesign":
        run_rdesign(pdb_directory, model_choice=model_train_choice, get_motifs=get_motifs)
    elif model_choice.lower() == "rhodesign":
        run_rhodesign(pdb_directory, contact_map=contact_map, get_motifs=get_motifs)
