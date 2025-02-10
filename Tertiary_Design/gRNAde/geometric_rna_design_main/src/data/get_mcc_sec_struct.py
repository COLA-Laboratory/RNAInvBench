import numpy as np
import pandas as pd
import RNA
import torch
from torchmetrics.functional.classification import binary_matthews_corrcoef


def dotbracket_to_adjacency(
        sec_struct: str,
        keep_pseudoknots: bool = False,
    ) -> np.ndarray:
    """
    Convert secondary structure in dot-bracket notation to 
    adjacency matrix.
    """
    n = len(sec_struct)
    adj = np.zeros((n, n), dtype=np.int8)
        
    if keep_pseudoknots == False:
        stack = []
        for i, db_char in enumerate(sec_struct):
            if db_char == '(':
                stack.append(i)
            elif db_char == ')':
                j = stack.pop()
                adj[i, j] = 1
                adj[j, i] = 1
    else:
        stack={
            '(':[],
            '[':[],
            '<':[],
            '{':[]
        }
        pop={
            ')':'(',
            ']':'[',
            '>':"<",
            '}':'{'
        }
        for i, db_char in enumerate(sec_struct):
            if db_char in stack:
                stack[db_char].append((i, db_char))
            elif db_char in pop:
                forward_bracket = stack[pop[db_char]].pop()
                adj[forward_bracket[0], i] = 1
                adj[i, forward_bracket[0]] = 1    
    return adj



df = pd.read_csv("InverseRNACASP15Results.csv")
casp_seqs = []
with open("CASP15_seqs.txt", "r") as file:
    for line in file.readlines():
        casp_seqs.append(line.strip())

seqs = df["sequences"].to_list()
print(len(seqs))
print(len(casp_seqs))
i = 0
mcc_scores = []
recovery_scores = []
nucleotide_map = {'A': 0, 'G': 1, 'C': 2, 'U': 3, '_': 4}


for seq, casp_seq, casp_struct in zip(seqs, casp_seqs, df["structures"].to_list()):
    seq = seq.split(" ")[0]
    if len(seq) == len(casp_seq):
        mcc_scores.append(
                binary_matthews_corrcoef(
                    torch.tensor(dotbracket_to_adjacency(list(RNA.fold(seq)[0]))),
                    torch.tensor(dotbracket_to_adjacency(list(RNA.fold(casp_seq)[0]))),
                ).float().mean())
        seq_tensor = torch.tensor([nucleotide_map[char] for char in seq], dtype=torch.long)
        casp_seq_tensor = torch.tensor([nucleotide_map[char] for char in casp_seq], dtype=torch.long)

        recovery_scores.append(
            seq_tensor.eq(casp_seq_tensor).float().cpu().numpy().mean()
        )
    # recovery = #calculate recovery
    # ssc = RNA.fold(seq)[0], #calculate ssc
print(np.mean(mcc_scores))
print(np.mean(recovery_scores))