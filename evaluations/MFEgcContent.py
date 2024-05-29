import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

import ViennaRNA

def MFE_GCcontent(result_folder='results', plot='both'):
    """
    GC content of designed sequence (stability), use ViennaRNA to predict the MFE structure and compute.
    reconsider if you want to give a target gc content to algorithms as a constraint before you call this evaluation.

    parameters:
    - result_folder (str): folder path of results of one algorithm.
    - plot (str): plot or not. plot gc/mfe/both volin

    return:
    - df_MFEgc (dict): df containing MFE and corresponding gc content.
    """
    df_MFEgc = pd.DataFrame(columns=['method', 'gc_contents', 'predict_structures', 'mfes'])
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
        fig, axes = plt.subplots(1, num_violins, figsize=(4*num_violins, 6), sharey=True)
    if plot == 'both':
        fig2, axes2 = plt.subplots(1, num_violins, figsize=(4*num_violins, 6))

    for i, (method, files) in enumerate(method_groups.items()):
        gc_list = []
        preStruct_list = []
        mfe_list = []
        for file in files:
            if file.endswith('.pkl'):
                file_path = os.path.join(result_folder, file)
                df = pd.read_pickle(file_path)  # read the result pickle
                try:
                    if df['sequence'].any():
                        # compute MFE structure
                        structure_pre, mfe = ViennaRNA.fold(df['sequence'][0])
                        # compute gc content
                        gc = paired_gcContent(df['sequence'][0], structure_pre)
                        gc_list.append(gc)
                        preStruct_list.append(structure_pre)
                        mfe_list.append(mfe)
                except KeyError:
                    continue
        new_row = pd.DataFrame({'method': method, 'gc_contents': [gc_list], 'predict_structures': [preStruct_list], 'mfes': [mfe_list]})
        df_MFEgc = pd.concat([df_MFEgc, new_row], ignore_index=True)

        # plot
        if plot == 'both':
            _plot(i, method, 'gc', axes, gc_list)
            _plot(i, method, 'mfe', axes2, mfe_list)
        else:
            _plot(i, method, plot, axes, eval(plot+'_list'))

    if plot:
        plt.tight_layout()
        plt.show()

    return df_MFEgc


def paired_gcContent(sequence, structure):
    # compute the gc content of paired region
    gc_paired_count = 0
    total_paired_count = 0

    stack = []
    for base, symbol in zip(sequence, structure):
        if symbol == '(':
            # if '(', press into stack
            stack.append(base)
        elif symbol == ')':
            # if ')', pop stack top and judge if it's GC pair
            if stack:
                left_base = stack.pop()
                if (left_base == 'G' and base == 'C') or (left_base == 'C' and base == 'G'):
                    gc_paired_count += 1
                total_paired_count += 1

    if total_paired_count > 0:
        return gc_paired_count / total_paired_count
    else:
        return 0.0

def _plot(i, method, key, axes, list):
    ax = axes[i]
    ax.violinplot(list, showmeans=True, showextrema=True, showmedians=True)
    ax.set_title(method)
    ax.set_xticks([])
    if i == 0:
        if key == 'gc':
            ax.set_ylabel("gc content")
        elif key == 'mfe':
            ax.set_ylabel("MFE")