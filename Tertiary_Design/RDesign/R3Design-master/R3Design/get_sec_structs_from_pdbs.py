import json
import numpy as np
import pandas as pd
from src.data.sec_struct_utils import pdb_to_sec_struct
import argparse
import os
import cpdb
from src.constants import (
    RNA_ATOMS, 
    RNA_NUCLEOTIDES, 
    PURINES,
    PYRIMIDINES,
    FILL_VALUE
)

def get_seq_from_pdb_file(filepath):
    df = cpdb.parse(filepath, df=True)

    # create unique residue id
    df["residue_id"] = (
        df["chain_id"]
        + ":"
        + df["residue_name"]
        + ":"
        + df["residue_number"].astype(str)
    )

    # get sequence
    nt_list = [res.split(":")[1] for res in df.residue_id.unique()]
    # replace non-standard nucleotides with placeholder
    nt_list = [nt if nt in RNA_NUCLEOTIDES else "_" for nt in nt_list]
    sequence = "".join(nt_list)
    return sequence

def grnade_get_ss(filepath):
    seq = get_seq_from_pdb_file(filepath)
    return pdb_to_sec_struct(filepath, seq) 


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert pdb files to secondary structures.")
    parser.add_argument('--pdb_path', type=str, required=True, help='Pathway to get PDB files')
    parser.add_argument('--ss_output_path', type=str, required=True, help='Output file and pathway for SSes')
    args = parser.parse_args()
    pdb_path = args.pdb_path
    ss_output_path = args.ss_output_path
    pdb_files = os.listdir(pdb_path)
    pdb_files = [file for file in pdb_files if file[-4:] == ".pdb"]
    sec_structs = []
    sequences = []
    for file in pdb_files:
        pdb_file_pathway = pdb_path+"/"+file
        pdb_file_seq = get_seq_from_pdb_file(pdb_file_pathway)
        try:
            sec_struct = pdb_to_sec_struct(pdb_file_pathway, pdb_file_seq, keep_pseudoknots=True)
        except:
            
            print("failed to find sec_struct")
            sec_struct = None
        sec_structs.append(sec_struct)
        sequences.append(pdb_file_seq)

    with open(ss_output_path+".txt", "w") as file:
        for pdb_file, sec_struct, seq in zip(pdb_files, sec_structs, sequences):
            row = {"pdb_file": pdb_file, "seq": seq, "struct": sec_struct}
            file.write(json.dumps(row)+"\n")