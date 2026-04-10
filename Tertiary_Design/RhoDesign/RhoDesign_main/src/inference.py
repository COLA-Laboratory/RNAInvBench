from RhoDesign import RhoDesignModel
from alphabet import Alphabet
import torch
import torch.nn.functional as F
from tqdm import tqdm
import numpy as np
import random
from util import load_structure, extract_coords_from_structure, seq_rec_rate, CoordBatchConverter, load_structure
from alphabet import Alphabet

import os
import argparse
import random
seed = random.randint(1, 10000)
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed_all(seed)

alphabet = Alphabet(['A','G','C','U','X'])

class args_class:  # use the same param as esm-if1, waiting to be adjusted...
    def __init__(self, encoder_embed_dim, decoder_embed_dim, dropout):
        self.local_rank = int(os.getenv("LOCAL_RANK", -1))
        self.device_id = [0, 1, 2, 3, 4, 5, 6, 7]
        self.epochs = 100
        self.lr = 1e-5
        self.batch_size = 1
        self.encoder_embed_dim = encoder_embed_dim
        self.decoder_embed_dim = decoder_embed_dim
        self.dropout = dropout
        self.gvp_top_k_neighbors = 15
        self.gvp_node_hidden_dim_vector = 256
        self.gvp_node_hidden_dim_scalar = 512
        self.gvp_edge_hidden_dim_scalar = 32
        self.gvp_edge_hidden_dim_vector = 1
        self.gvp_num_encoder_layers = 3
        self.gvp_dropout = 0.1
        self.encoder_layers = 3
        self.encoder_attention_heads = 4
        self.attention_dropout = 0.1
        self.encoder_ffn_embed_dim = 512
        self.decoder_layers = 3
        self.decoder_attention_heads = 4
        self.decoder_ffn_embed_dim = 512


def get_sequence_loss(model, coords, _device, seq, ss_ct_map, confidence=None):
    batch_converter = CoordBatchConverter(model.decoder.dictionary)
    batch_coords, confidence, _, tokens, padding_mask, ss_ct_map = batch_converter(
        [(coords, confidence, seq, ss_ct_map)]
    )
    batch_coords = batch_coords.to(_device)
    confidence = confidence.to(_device)
    padding_mask = padding_mask.bool().to(_device)
    ss_ct_map = ss_ct_map.to(_device)

    c = batch_coords[:, :, [0, 1, 2], :]
    adc = batch_coords[:, :, :, :]

    prev_output_tokens = tokens[:, :-1].to(_device)
    target = tokens[:, 1:].to(_device)
    target_padding_mask = (target == alphabet.padding_idx)
    
    logits, _ = model.forward(c, adc, ss_ct_map, padding_mask, confidence, prev_output_tokens)

    # logits = logits.view(-1, logits.size(-1))
    # target = target.view(-1)
    # target_padding_mask = target_padding_mask.view(-1)
    loss = F.cross_entropy(logits, target, reduction='none')
    loss = loss.cpu().detach().numpy()
    mask = (~target_padding_mask).cpu().numpy()
    return loss, mask

def score_sequence(model, coords, _device, seq, ss_ct_map):
    loss, mask = get_sequence_loss(model, coords, _device, seq, ss_ct_map)
    nll = np.sum(loss * mask) / np.sum(mask)
    return nll

def score_backbone(model, coords, seq, _device, ss_ct_map):
    nll = score_sequence(model, coords, _device, seq, ss_ct_map)
    return np.exp(nll)


def eval_ppl(model, pdb_list, ss_ct_map, model_path):
    """
    fpath: path to pdb file
    """
    
    temp=torch.load(model_path) 
    model.load_state_dict(temp)  
    model.eval()

    with torch.no_grad():
        pfile = './../data/test/'
        ssfile = './../data/test_ss/'
        ppl = []
        wrong_ppl = []
        wrong_p = []
        for i, pdb in tqdm(enumerate(pdb_list)):
            fpath = pfile+str(i)+'.pdb'
            coords, seq = extract_coords_from_structure(pdb)
            ppl_v = score_backbone(model,coords,seq,_device, ss_ct_map)
            ppl.append(ppl_v)
    return np.mean(ppl)



def eval(model,pdb_path,ss_path,save_path,_device,temp=1e-5):
    """
    fpath: path to pdb file
    """

    model_path = 'checkpoints/ss_apexp_best.pth'

    name = pdb_path.split('/')[-2]
    
    model_dir=torch.load(model_path) 
    model.load_state_dict(model_dir)  
    model.eval()
    rc = []
    
    ss_ct_map = np.load(ss_path)
    pdb = load_structure(pdb_path)
    coords, seq = extract_coords_from_structure(pdb)

    pred_seq = model.sample(coords,ss_ct_map,_device,temperature=temp)
    rc_value = seq_rec_rate(seq,pred_seq)
    perplexity = eval_ppl(model, [pdb], ss_ct_map, model_path)
    rc.append(rc_value)
    with open(os.path.join(save_path,f'{name}_with2d.fasta'),'w') as f:
        f.write(f'>{name}_with2d'+'\n')
        f.write(pred_seq+'\n')
    
    print('sequence: ' + pred_seq)
    print('recovery rate: ' + str(np.mean(rc)))
    print('perplexity: ' + str(perplexity))
    print('true sequence: ' + seq)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Description of your script')
    parser.add_argument('-pdb', '--pdb_file', type=str,help='path to the pdb file',required=True)
    parser.add_argument('-ss', '--secondary_structure', type=str, help='path to the secondary structure file',required=True)
    parser.add_argument('-save', '--save_path', type=str, help='path to the save directory',required=True)
    parser.add_argument('-device', '--device', default=1, type=int, help='Assign the device to run the model')
    parser.add_argument('-temp', '--temperature', default=1e-5, type=float, help='temperature for sampling')
    args = parser.parse_args()

    pdb = args.pdb_file
    ss = args.secondary_structure
    save_path = args.save_path
    _device = args.device
    temp = args.temperature

    model_args = args_class(512,512,0.1)
    dictionary = Alphabet(['A','G','C','U','X'])
    model = RhoDesignModel(model_args, dictionary).cuda(device=_device)
    eval(model,pdb,ss,save_path,_device,temp)
