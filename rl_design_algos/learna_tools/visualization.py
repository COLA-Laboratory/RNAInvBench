import subprocess
import datetime
import numpy as np
import pandas as pd

from pathlib import Path
from collections import defaultdict

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# from Bio.Seq import Seq
# from Bio import motifs
import logomaker as lm


def display_varna_plot(path):
    image = mpimg.imread(path)
    plt.imshow(image)
    plt.axis("off")   # turns off axes
    plt.axis("tight")  # gets rid of white border
    plt.axis("image")  # square up the image instead of filling the "figure" space
    plt.show()



def plot_with_varna(pairs, sequence, Id, plot_dir='plots', plot_type='radiate', resolution='2.0'):
        Path(plot_dir).mkdir(exist_ok=True, parents=True)
        multiplets = get_multiplets(pairs)
        watson_pairs, wobble_pairs, noncanonical_pairs = type_pairs(pairs, sequence)
        lone_bp = lone_pair(pairs)
        tertiary_bp = multiplets + noncanonical_pairs + lone_bp
        tertiary_bp = [list(x) for x in set(tuple(x) for x in tertiary_bp)]

        str_tertiary = []
        for i,I in enumerate(tertiary_bp):
            if i==0:
                str_tertiary += ('(' + str(I[0]+1) + ',' + str(I[1]+1) + '):color=""#FFFF00""')
            else:
                str_tertiary += (';(' + str(I[0]+1) + ',' + str(I[1]+1) + '):color=""#FFFF00""')

        tertiary_bp = ''.join(str_tertiary)

        ct_path = Path(plot_dir, f"{Id}.ct")

        ct_file_output(pairs, sequence, Id, ct_path)

        varna_path = str(Path(plot_dir, f"{Id}").resolve())
        if plot_type == 'radiate':
            subprocess.Popen(["varna", '-i', str(ct_path.resolve()), '-o', varna_path + '_radiate.png', '-algorithm', 'radiate', '-resolution', resolution, '-bpStyle', 'lw', '-auxBPs', tertiary_bp], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]  # -> to plot structure only: '-basesStyle1', 'fill=#FFFFFF,label=#FFFFFF,number=#FFFFFF', '-applyBasesStyle1on', ','.join([str(i) for i in range(1, len(sequence)+1)]),
        elif plot_type == 'line':
            subprocess.Popen(["varna", '-i', str(ct_path.resolve()), '-o', varna_path + '_line.png', '-algorithm', 'line', '-resolution', resolution, '-bpStyle', 'lw', '-auxBPs', tertiary_bp], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]



def plot_df_with_varna(df, plot_dir='plots', name='', plot_type='radiate', resolution='2.0', show=False):
    if name:
        plot_dir = Path(plot_dir, name)
    else:
        current_time = f'{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}'
        plot_dir = Path(plot_dir, current_time)
    Ids = defaultdict(list)
    for i, row in df.iterrows():
        pairs = db2pairs(row['structure'])
        sequence = row['sequence']
        Ids[str(row['Id'])].append(str(row['Id']))
        Id = str(row['Id']) + '_' + str(len(Ids[str(row['Id'])])-1)
        plot_with_varna(pairs, sequence, Id, plot_dir=plot_dir, plot_type=plot_type, resolution=resolution)
        if show:
            display_varna_plot(Path(plot_dir, f"{Id}" + '_' + f"{plot_type}.png").resolve())


def plot_weblogo(df):
    """
    Seems like the weblogo option of Biopython is not working properly...
    This function simply does not work.
    """
    df['sequence'] = df['sequence'].apply(lambda x: x.replace('U', 'T'))
    instances = df['sequence']
    m = motifs.create(instances)
    m.weblogo('logo.png')

def plot_sequence_logo(df, plotting_dir='plots', name='current_sequence_logo', show=False):
    """
    Create a logo from a dataframe that contains a 'sequence' column.

    All sequences have to be of the same length.
    """

    outdir = Path(plotting_dir, 'logos')
    outdir.mkdir(exist_ok=True, parents=True)
    sequences = df['sequence'].to_list()
    counts_mat = lm.alignment_to_matrix(sequences) / len(sequences)
    counts_mat.head()
    logo = lm.Logo(counts_mat, stack_order='small_on_top', figsize=[30, 6])  # , vsep=0.05, stack_order='small_on_top')
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left','bottom', 'top', 'right'], visible=True, linewidth=2)
    
    # Hide the x-axis labels and ticks
    logo.ax.set_xticks([])
    logo.ax.set_xticklabels([])
    
    # Hide the y-axis labels and ticks
    logo.ax.set_yticks([])
    logo.ax.set_yticklabels([])
    
    save_path = Path(outdir, f"{name}.png")
    plt.savefig(save_path)
    if show:
        plt.show()
    plt.close('all')


def plot_structure_logo(df, plotting_dir='plots', name='current_structure_logo', show=False):
    """
    Create a logo from a dataframe that contains a 'structure' column.

    All structures have to be of the same length.
    """
    outdir = Path(plotting_dir, 'logos')
    outdir.mkdir(exist_ok=True, parents=True)
    sequences = df['structure'].to_list()
    counts_mat = lm.alignment_to_matrix(sequences) / len(sequences)
    counts_mat.head()
    logo = lm.Logo(counts_mat, stack_order='small_on_top')  # , vsep=0.05, stack_order='small_on_top')
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left','bottom'], visible=True, linewidth=2)
    save_path = Path(outdir, f"{name}.png")
    plt.savefig(save_path)
    if show:
        plt.show()
    plt.close('all')


################################################################################
# following code snippets originate from SPOT-RNA repo at https://github.com/jaswindersingh2/SPOT-RNA
# we use them for plotting with varna in the visualization module
################################################################################


# copy-paste from SPOT-RNA2 source code
def lone_pair(pairs):
    lone_pairs = []
    pairs.sort()
    for i, I in enumerate(pairs):
        if ([I[0] - 1, I[1] + 1] not in pairs) and ([I[0] + 1, I[1] - 1] not in pairs):
            lone_pairs.append(I)

    return lone_pairs


# copy-paste from SPOT-RNA2 source code
def type_pairs(pairs, sequence):
    sequence = [i.upper() for i in sequence]
    # seq_pairs = [[sequence[i[0]],sequence[i[1]]] for i in pairs]

    AU_pair = []
    GC_pair = []
    GU_pair = []
    other_pairs = []
    for i in pairs:
        if [sequence[i[0]],sequence[i[1]]] in [["A","U"], ["U","A"]]:
            AU_pair.append(i)
        elif [sequence[i[0]],sequence[i[1]]] in [["G","C"], ["C","G"]]:
            GC_pair.append(i)
        elif [sequence[i[0]],sequence[i[1]]] in [["G","U"], ["U","G"]]:
            GU_pair.append(i)
        else:
            other_pairs.append(i)
    watson_pairs_t = AU_pair + GC_pair
    wobble_pairs_t = GU_pair
    other_pairs_t = other_pairs
        # print(watson_pairs_t, wobble_pairs_t, other_pairs_t)
    return watson_pairs_t, wobble_pairs_t, other_pairs_t



# copy-paste from SPOT-RNA2 source code
def multiplets_pairs(pred_pairs):

    pred_pair = [i[:2] for i in pred_pairs]
    temp_list = flatten(pred_pair)
    temp_list.sort()
    new_list = sorted(set(temp_list))
    dup_list = []
    for i in range(len(new_list)):
        if (temp_list.count(new_list[i]) > 1):
            dup_list.append(new_list[i])

    dub_pairs = []
    for e in pred_pair:
        if e[0] in dup_list:
            dub_pairs.append(e)
        elif e[1] in dup_list:
            dub_pairs.append(e)

    temp3 = []
    for i in dup_list:
        temp4 = []
        for k in dub_pairs:
            if i in k:
                temp4.append(k)
        temp3.append(temp4)

    return temp3


# copy-paste from SPOT-RNA2 source code
def multiplets_free_bp(pred_pairs, y_pred):
    L = len(pred_pairs)
    multiplets_bp = multiplets_pairs(pred_pairs)
    save_multiplets = []
    while len(multiplets_bp) > 0:
        remove_pairs = []
        for i in multiplets_bp:
            save_prob = []
            for j in i:
                save_prob.append(y_pred[j[0], j[1]])
            remove_pairs.append(i[save_prob.index(min(save_prob))])
            save_multiplets.append(i[save_prob.index(min(save_prob))])
        pred_pairs = [k for k in pred_pairs if k not in remove_pairs]
        multiplets_bp = multiplets_pairs(pred_pairs)
    save_multiplets = [list(x) for x in set(tuple(x) for x in save_multiplets)]
    assert L == len(pred_pairs)+len(save_multiplets)
    #print(L, len(pred_pairs), save_multiplets)
    return pred_pairs, save_multiplets


def flatten(x):
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, str):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

def ct_file_output(pairs, seq, id, save_result_path):

    col1 = np.arange(1, len(seq) + 1, 1)
    col2 = np.array([i for i in seq])
    col3 = np.arange(0, len(seq), 1)
    col4 = np.append(np.delete(col1, 0), [0])
    col5 = np.zeros(len(seq), dtype=int)

    for i, I in enumerate(pairs):
        col5[I[0]] = int(I[1]) + 1
        col5[I[1]] = int(I[0]) + 1
    col6 = np.arange(1, len(seq) + 1, 1)
    temp = np.vstack((np.char.mod('%d', col1), col2, np.char.mod('%d', col3), np.char.mod('%d', col4),
                      np.char.mod('%d', col5), np.char.mod('%d', col6))).T
    #os.chdir(save_result_path)
    #print(os.path.join(save_result_path, str(id[0:-1]))+'.spotrna')
    np.savetxt(save_result_path, (temp), delimiter='\t\t', fmt="%s", header=str(len(seq)) + '\t\t' + str(id) + '\t\t' + ' output\n' , comments='')  # changed outpath from original SPOT-RNA code

    return


def get_multiplets(pairs):
    pair_dict = defaultdict(list)
    for p in pairs:
        pair_dict[p[0]].append(p[1])
        pair_dict[p[1]].append(p[0])

    multiplets = []

    for k, v in pair_dict.items():
        if len(v) > 1:
            for p in v:
                multiplets.append(tuple(sorted([k, p])))
    multiplets = list(set(multiplets))

    return multiplets

################################################################################

def db2pairs(structure, start_index=0):
    """
    Converts dot-bracket string into a list of pairs.

    Input:
      structure <string>: A sequence in dot-bracket format.
      start_index <int>: Starting index of first nucleotide (default zero-indexing).

    Returns:
      pairs <list>: A list of tuples of (index1, index2, pk_level).

    """
    level_stacks = defaultdict(list)
    closing_partners = {')': '(', ']': '[', '}': '{', '>': '<'}
    levels = {')': 0, ']': 1, '}': 2, '>': 3}

    pairs = []

    for i, sym in enumerate(structure, start_index):
        if sym == '.':
            continue
        # high order pks are alphabetical characters
        if sym.isalpha():
            if sym.isupper():
                level_stacks[sym].append(i)
            else:
                try:  # in case we have invalid preditions, we continue with next bracket
                    op = level_stacks[sym.upper()].pop()
                    pairs.append((op, i,
                                  ord(sym.upper()) - 61))  # use asci code if letter is used to asign PKs, start with level 4 (A has asci code 65)
                except:
                    continue
        else:
            if sym in closing_partners.values():
                level_stacks[sym].append(i)
            else:
                try:  # in case we have invalid preditions, we continue with next bracket
                    op = level_stacks[closing_partners[sym]].pop()
                    pairs.append([op, i, levels[sym]])
                except:
                    continue
    return sorted(pairs, key=lambda x: x[0])
