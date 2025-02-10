from gRNAde import perplexity_from_pdb_file, perplexity_from_directory
import json
import argparse
from gRNAde.geometric_rna_design_main.src.data.featurizer import RNAGraphFeaturizer
import torch.nn.functional as F


featurizer = RNAGraphFeaturizer(
    split = "test",  # set to 'train' to use noise augmentation
    radius = RADIUS,
    top_k = TOP_K,
    num_rbf = NUM_RBF,
    num_posenc = NUM_POSENC,
    max_num_conformers = max_num_conformers,
    noise_scale = NOISE_SCALE
)


def perplexity_from_directory(
        self,
        seq: str,
        directory_filepath: str,
        temperature = 1.0,
        seed = 0
    ):
    """
    Compute perplexity of an RNA sequences for a set of backbones 
    from a directory of PDB files,
    i.e. P (sequence | backbone conformational ensemble)

    Args:
        seq (str): RNA sequence
        directory_filepath (str): filepath to directory of PDB files
        temperature (float): temperature for sampling
        seed (int): random seed for reproducibility

    Returns:
        perplexity (float): perplexity for RNA sequence
    """
    pdb_filelist = []
    for pdb_filepath in os.listdir(directory_filepath):
        if pdb_filepath.endswith(".pdb"):
            pdb_filelist.append(os.path.join(directory_filepath, pdb_filepath))
    featurized_data, raw_data = featurizer.featurize_from_pdb_filelist(pdb_filelist)
    return perplexity(seq, raw_data, featurized_data, temperature, seed)

@torch.no_grad()
def perplexity(
    self, 
    seq: str,
    raw_data: dict, 
    temperature = 1.0,
    seed = 0
):
    """
    Compute perplexity of an RNA sequence conditioned on 
    one or more backbones from raw data.

    Args:
        seq (str): RNA sequence
        raw_data (dict): Raw RNA data dictionary with keys:
            - sequence (str): RNA sequence of length `num_res`.
            - coords_list (Tensor): Backbone coordinates with shape
                `(num_conf, num_res, num_bb_atoms, 3)`.
            - sec_struct_list (List[str]): Secondary structure for each
                conformer in dotbracket notation.
        featurized_data (torch_geometric.data.Data): featurized RNA data
        temperature (float): temperature for sampling
        seed (int): random seed for reproducibility
    
    Returns:
        perplexity (float): perplexity for RNA sequence
    """  
    if raw_data['coords_list'][0].shape[1] == 3:
        # Expected input: num_conf x num_res x num_bb_atoms x 3
        # Backbone atoms: (P, C4', N1 or N9)
        pass
    elif raw_data['coords_list'][0].shape[1] == len(RNA_ATOMS):
        coords_list = []
        for coords in raw_data['coords_list']:
            # Only keep backbone atom coordinates: num_res x num_bb_atoms x 3
            coords = get_backbone_coords(coords, raw_data['sequence'])
            # Do not add structures with missing coordinates for ALL residues
            if not torch.all((coords == FILL_VALUE).sum(axis=(1,2)) > 0):
                coords_list.append(coords)

        if len(coords_list) > 0:
            # Add processed coords_list to self.data_list
            raw_data['coords_list'] = coords_list
    else:
        raise ValueError(f"Invalid number of atoms per nucleotide in input data: {raw_data['coords_list'][0].shape[1]}")

    if featurized_data is None:
        # featurize raw data
        featurized_data = featurizer.featurize(raw_data)

    # transfer data to device
    featurized_data = featurized_data.to(self.device)

    # convert sequence to tensor
    _seq = torch.as_tensor(
        [featurizer.letter_to_num[residue] for residue in seq], 
        device=self.device, 
        dtype=torch.long
    )
    featurized_data.seq = _seq

    # raw logits for perplexity calculation: seq_len x out_dim
    logits = model.forward(featurized_data)

    # compute perplexity
    perplexity = torch.exp(F.cross_entropy(
        logits / temperature, 
        _seq,
        reduction="none"
    ).mean()).cpu().numpy()
    
    return perplexity

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get perplexity from specified RNA file")
    parser.add_argument('--data', type=str, required=True, help='Link to JSON file containing pdb files and sequences')
    args = parser.parse_args()
    data_pathway = args.data
    data = []
    try:
        with open(data_pathway) as file:
            for line in file.readlines():
                line = json.loads(line)
                data.append({"file_pathway": line["pdb_file"], "seq": line["seq"]})
    except:
        print("File not found.")
        
    perplexities = []
    for row in data:
        perplexities.append(calculate_perplexity(row["file_pathway"], row["seq"]))
        print(f"File: {row['file_pathway']}, Perplexity: {perplexities[-1]}")

    with open("perplexities.json", "w") as file:
        for i, row in enumerate(data):
            row["perplexity"] = perplexities[i]
            file.write(json.dumps(row)+"\n")