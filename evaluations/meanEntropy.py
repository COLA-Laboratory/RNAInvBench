import os
import pandas as pd
import math
from collections import Counter
import matplotlib.pyplot as plt

def meanEntropy(result_folder='results', plot=False):
    """
    design multiple times for one target structure, use entropy to evaluate the diversity of designed sequences.
    parameters:
    - result_folder (str): folder path of results of one algorithm.
    - plot (str): plot or not. plot gc/mfe/both volin

    return:
    -
    """
    df_meanEntropy = pd.DataFrame(columns=['method', 'mean_entropy'])
    file_list = os.listdir(result_folder)
    # group file by method name
    method_groups = {}
    for file_name in file_list:
        method_name = file_name.split('_')[0]
        if method_name not in method_groups:
            method_groups[method_name] = []
        method_groups[method_name].append(file_name)

    if plot:
        num_violins = len(method_groups)
        fig, axes = plt.subplots(1, num_violins, figsize=(4 * num_violins, 6), sharey=True)

    for i, (method, files) in enumerate(method_groups.items()):
        meanEntropy_list = []
        for file in files:
            if file.endswith('.pkl'):
                file_path = os.path.join(result_folder, file)
                df = pd.read_pickle(file_path)  # read the result pickle
                try:
                    if df['sequence'].any():
                        seq_list = df['sequence'].tolist()
                        try:
                            seqCons = df['seq_constraints'][0]
                        except KeyError:
                            # seq_constraints doesn't exists
                            seqCons = 'N' * len(seq_list[0])
                        mean_etp = cal_entropy(seq_list, seqCons)
                        meanEntropy_list.append(mean_etp)
                except KeyError:
                    continue
        new_row = pd.DataFrame({'method': method, 'mean_entropy': [meanEntropy_list]})
        df_meanEntropy = pd.concat([df_meanEntropy, new_row], ignore_index=True)

        # plot
        if plot:
            ax = axes[i]
            ax.violinplot(meanEntropy_list, showmeans=True, showextrema=True, showmedians=True)
            ax.set_title(method)
            ax.set_xticks([])
            if i == 0:
                ax.set_ylabel("Mean Entropy")
    if plot:
        plt.tight_layout()
        plt.show()

    return df_meanEntropy

def cal_entropy(sequences_list, constraints):
    # calculate the mean entropy of one target structure
    L = len(sequences_list[0])

    # constraint positions
    constraint_positions = set()
    for i, nucleotide in enumerate(constraints):
        if nucleotide != 'N':
            constraint_positions.add(i)
    nc = len(constraint_positions)

    # entropy of every position
    position_entropies = []
    for i in range(L):
        if i not in constraint_positions:
            nucleotides = [seq[i] for seq in sequences_list]
            counts = Counter(nucleotides)
            total = sum(counts.values())
            position_entropy = 0
            for nucleotide, count in counts.items():
                p = count / total
                position_entropy -= p * math.log2(p)
                position_entropies.append(position_entropy)

    # cal the mean entropy
    mean_entropy = sum(position_entropies) / (L - nc)
    return mean_entropy


if __name__ == "__main__":
    df = meanEntropy("../results/eterna100_v2", plot=True)
    pass