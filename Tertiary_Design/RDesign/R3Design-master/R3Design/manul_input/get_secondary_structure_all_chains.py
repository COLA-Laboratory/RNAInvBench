# get_secondary_structure.py
import os
from moderna import *
from Bio import PDB

def get_chains(pdb_path):
    rna_residue_map = {
    'A': 'A',  # Adenine
    'C': 'C',  # Cytosine
    'G': 'G',  # Guanine
    'U': 'U',  # Uracil
    }
    parser = PDB.PDBParser()
    struct = parser.get_structure("RNA", pdb_path)
    chains = []
    seq = ""
    for model in struct:
        for chain in model:
            chains.append(chain.id)
    return list(set(chains)) # remove duplicates

def get_secondary_structure(pdb_path, chain=False):
    """
    Load a PDB file and return its secondary structure.
    """
    if chain == False:
        chains = get_chains(pdb_path)
    else:
        chains = chain
    secondary_structs = []
    # print(pdb_path+" : "+chains)
    for chain in chains:
        try:
            t = load_template(pdb_path, chain)
            clean_structure(t)
            if get_secstruc(t) == "":
                raise ValueError("Returned Blank")
            secondary_structs.append(get_secstruc(t))
        except Exception as e:
            print("Error processing {}: {}".format(pdb_path, e))
    if len(secondary_structs) > 0:
        return secondary_struct
    else:
        return None

if __name__ == "__main__":
    import sys
    # pdb_path = sys.argv[1]
    # chain_name = sys.argv[2]
    pdb_path = '/home/jack/Tertiary_Models/gRNAde/geometric-rna-design-main/data/3D_Puzzles/4l81.pdb'
    chain_name = 'A'
    ss = get_secondary_structure(pdb_path,chain_name)
    print(ss)
