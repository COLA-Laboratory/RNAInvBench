import torch
from tqdm import tqdm
import numpy as np
import random
from models import *
from utils import *
from diffusion import NoiseScheduleVP
from sampling import get_sampling_fn
from datasets import utils as du
from torch_geometric.data import Batch
import logging
import pickle
import functools
import tree
import copy
import time
import torch.nn.functional as F
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from datetime import datetime
import subprocess
from torchmetrics.functional.classification import binary_matthews_corrcoef



def set_random_seed(config):
    seed = config.seed

    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

    np.random.seed(seed)
    random.seed(seed)

    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def get_optimizer(config, params):
    """Return a flax optimizer object based on `config`."""
    if config.optim.optimizer == 'Adam':
        optimizer = optim.Adam(params, lr=config.optim.lr, betas=(config.optim.beta1, 0.999), eps=config.optim.eps,
                               weight_decay=config.optim.weight_decay)
    elif config.optim.optimizer == 'AdamW':
        optimizer = torch.optim.AdamW(params, lr=config.optim.lr, amsgrad=True, weight_decay=1e-12)
    else:
        raise NotImplementedError(
            f'Optimizer {config.optim.optimizer} not supported yet!'
        )
    return optimizer

"""
def calculate_perplexity(
    model, 
    seq: str,
    raw_data: dict, 
    featurized_data = None, 
    temperature = 1.0,
    seed = 0
):
    """
"""
    Compute perplexity of an RNA sequence conditioned on 
    one or more backbones from raw data.
    % Modified from gRNAde
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
"""
    print(raw_data.keys())
    # ['seq', 'coords', 'node_s', 'node_v', 'edge_s', 'edge_v', 'edge_index', 'mask', 'z_t']
    if raw_data['coords'][0].shape[1] == 3:
        # Expected input: num_conf x num_res x num_bb_atoms x 3
        # Backbone atoms: (P, C4', N1 or N9)
        pass
    elif raw_data['coords'][0].shape[1] == len(RNA_ATOMS):
        coords_list = []
        for coords in raw_data['coords']:
            # Only keep backbone atom coordinates: num_res x num_bb_atoms x 3
            coords = get_backbone_coords(coords, raw_data['seq'])
            # Do not add structures with missing coordinates for ALL residues
            if not torch.all((coords == FILL_VALUE).sum(axis=(1,2)) > 0):
                coords_list.append(coords)

        if len(coords_list) > 0:
            # Add processed coords_list to self.data_list
            raw_data['coords'] = coords_list
    else:
        raise ValueError(f"Invalid number of atoms per nucleotide in input data: {raw_data['coords_list'][0].shape[1]}")

    featurized_data = featurized_data.to(self.device)

    # convert sequence to tensor
    _seq = torch.as_tensor(
        [featurizer.letter_to_num[residue] for residue in seq], 
        device=self.device, 
        dtype=torch.long
    )
    featurized_data.seq = _seq

    log_likelihoods = []
    for t in range(model.num_timesteps):
        noise_level = model.noise_schedule[t]
        # Forward pass through diffusion model
        denoised_logits = model.denoise_at_timestep(featurized_data, t)
        # Compute log-likelihood
        log_likelihood = F.log_softmax(denoised_logits / temperature, dim=-1)
        log_likelihoods.append(log_likelihood.mean())
    
    # Aggregate log-likelihoods
    perplexity = torch.exp(-torch.stack(log_likelihoods).mean()).cpu().numpy()
    
    return perplexity
"""

def calculate_perplexity(model, data, sampling_fn, temperature=1):
    try:
        # Ensure noise_level and time are correctly shaped tensors
        time = torch.tensor([1.0], device=data['node_s'].device)
        noise_level = torch.log(torch.tensor([1.0], device=data['node_s'].device))
        
      
        # Call model with explicit parameters
        logits = model(
            data, 
            time=time, 
            noise_level=noise_level
        )


        
        # log_probs = F.log_softmax(logits / temperature, dim=-1).view(-1, logits.size(-1))
        
        # Cross-entropy using original sequence
        seq_indices = data['seq'].view(-1)
        log_probs_flat = logits.view(-1, logits.size(-1))
        perplexity = torch.exp(F.cross_entropy(
            log_probs_flat / temperature, 
            seq_indices,
            reduction="none"
        ).mean()).cpu().detach().numpy()
        
        return perplexity

    
    except Exception as e:
        print(f"Error in perplexity calculation: {e}")
        raise


def predict_sec_struct(seq):
    eternafold_path = "gRNAde/geometric_rna_design_main/tools/EternaFold/EternaFold-master 2/src/contrafold"
    current_datetime = datetime.now().strftime("%Y%m%d_%H%M%S")
    fasta_file_path = f"fastas_temp/temp_{current_datetime}.fasta"
    cmd = [eternafold_path, "predict", fasta_file_path]
    SeqIO.write(SeqRecord(Seq(seq), id="temp"),fasta_file_path, "fasta")
    output = subprocess.run(cmd, check=True, capture_output=True).stdout.decode("utf-8")
    return [x for x in output.split("\n")[-2]]


def fetch_secondary_struct_consistency(real_seq, pred_seq):
    NUM_TO_LETTER = {0: 'A', 1: 'G', 2: 'C', 3: 'U'}
    real_seq = real_seq.view(-1)
    pred_seq = pred_seq.view(-1)
    real_seq_nucs = ''.join(NUM_TO_LETTER[num] for num in real_seq.cpu().numpy())
    pred_seq_nucs = ''.join(NUM_TO_LETTER[num] for num in pred_seq.cpu().numpy())
    real_sec_struct = predict_sec_struct(real_seq_nucs)
    pred_sec_struct = predict_sec_struct(pred_seq_nucs)
    pred_sec_struct = torch.tensor([1 if c in "()" else 0 for c in pred_sec_struct], dtype=torch.float32)
    real_sec_struct = torch.tensor([1 if c in "()" else 0 for c in real_sec_struct], dtype=torch.float32)
    return binary_matthews_corrcoef(pred_sec_struct,real_sec_struct).float().mean()
    
        
    

def vpsde_inference(config, save_folder,
                    pdb_file='./example/R1107.pdb'):
    """Runs inference for RNA inverse design in a given dir."""
    # Create directory for eval_folder
    os.makedirs(save_folder, exist_ok=True)

    # Initialize model
    model = create_model(config)
    ema = ExponentialMovingAverage(model.parameters(), decay=config.model.ema_decay)
    optimizer = get_optimizer(config, model.parameters())
    state = dict(optimizer=optimizer, model=model, ema=ema, step=0)

    model_size = sum(p.numel() for p in model.parameters()) * 4 / 2 ** 20
    print('model size: {:.1f}MB'.format(model_size))

    # Checkpoint name
    checkpoint_path = 'RiboDiffusion/RiboDiffusion-main/ckpts/exp_inf.pth'

    # Load checkpoint
    state = restore_checkpoint(checkpoint_path, state, device=config.device)
    ema.copy_to(model.parameters())

    # Initialize noise scheduler
    noise_scheduler = NoiseScheduleVP(config.sde.schedule, continuous_beta_0=config.sde.continuous_beta_0,
                                      continuous_beta_1=config.sde.continuous_beta_1)

    # Obtain data scalar and inverse scalar
    inverse_scaler = get_data_inverse_scaler(config)

    # Setup new sampling function for multi-state
    test_sampling_fn = get_sampling_fn(config, noise_scheduler, config.eval.sampling_steps, inverse_scaler)
    pdb2data = functools.partial(du.PDBtoData, num_posenc=config.data.num_posenc,
                                 num_rbf=config.data.num_rbf, knn_num=config.data.knn_num)

    fasta_dir = os.path.join(save_folder, 'fasta')
    os.makedirs(fasta_dir, exist_ok=True)

    # run inference on a single pdb file
    print('Start inference on a single pdb file')
    pdb_id = pdb_file.replace('.pdb', '')
    if '/' in pdb_id:
        pdb_id = pdb_id.split('/')[-1]
    struct_data = pdb2data(pdb_file)
    struct_data = tree.map_structure(lambda x:
                                     x.unsqueeze(0).repeat_interleave(config.eval.n_samples,
                                                                      dim=0).to(config.device),
                                     struct_data)
    samples = test_sampling_fn(model, struct_data)

    # save to fasta dir
    for i in range(len(samples)):
        du.sample_to_fasta(samples[i], pdb_id,
                           os.path.join(fasta_dir, pdb_id + '_' + str(i) + '.fasta'))

    recovery_ = samples.eq(struct_data['seq']).float().mean().item()
    print(f'{pdb_id}, recovery_rate {recovery_:.4f}')
    perplexity = calculate_perplexity(model, struct_data, test_sampling_fn)
    print(f'{pdb_id}, perplexity {perplexity:.4f}')
    ssc = fetch_secondary_struct_consistency(samples, struct_data["seq"])
    print(f'{pdb_id}, secondary_struct_consistency {ssc:.4f}')