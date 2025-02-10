from Bio.PDB import PDBParser, PDBIO
import os

# List of standard protein residues (for checking)
protein_residues = set(["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"])

def separate_rna_chains(pdb_file, output_folder):
    parser = PDBParser(QUIET=True)
    io = PDBIO()
    
    # Parse the structure
    if "Puzzle" in pdb_file:
        structure = parser.get_structure(pdb_file[11:16], pdb_file)
        print(pdb_file[12:16])
    else:
        structure = parser.get_structure(pdb_file[2:6], pdb_file)
        print(pdb_file[2:6])
    
    print(pdb_file)
    
    # Loop through chains and process each one
    pdb_chains = structure.get_chains()
    chain_number = 1  # Start a counter for unique chain numbers
    for chain in pdb_chains:
        # Check if the chain contains protein residues
        contains_protein = False
        for residue in chain.get_residues():
            if residue.get_resname() in protein_residues:
                contains_protein = True
                break  # No need to check further once we know it's a protein chain
        
        # If the chain contains protein residues, remove it from RNA processing
        if contains_protein:
            print(f"  -> Chain {chain.get_id()} contains protein atoms, skipping...")
        else:
            # Save the chain as a separate PDB file (RNA-only chain)
            output_filename = f"{output_folder}/{structure.get_id()}_{chain_number}_{chain.get_id()}.pdb"
            io.set_structure(chain)
            io.save(output_filename)
            print(f"  -> Chain {chain.get_id()} saved as RNA chain to {output_filename}")
        
        # Increment chain number for the next chain
        chain_number += 1

pdb_files = []
for root, dirs, files in os.walk("."):
    for file in files:
        if file.endswith(".pdb") and "checkpoint" not in file:
            pdb_files.append(os.path.join(root, file))
for pdb_file in pdb_files:
    separate_rna_chains(pdb_file, "split_pdbs2")
