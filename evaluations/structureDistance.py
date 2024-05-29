import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

def structureDistance(result_folder='results', plot=False):
    """
    structure distance between designed sequences(MFE structure) and target structure

    parameters:
    - result_folder (str): folder path of results of one algorithm.
    - plot (bool): plot or not

    return:
    distance_dict (dict): distances(every structure) of each method
    """
    df_distance = pd.DataFrame(columns=['method', 'distances'])
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

    for i, (method, files) in enumerate(method_groups.items()):
        distances_list = []
        for file in files:
            if file.endswith('.pkl'):
                file_path = os.path.join(result_folder, file)
                df = pd.read_pickle(file_path)  # read the result pickle
                try:
                    if df['distance'].any():
                        distances_list.append(df['distance'][0])
                except KeyError:
                    continue
        new_row = pd.DataFrame({'method': method, 'distances': [distances_list]})
        df_distance = pd.concat([df_distance, new_row], ignore_index=True)

        # plot
        if plot:
            ax = axes[i]
            ax.violinplot(distances_list, showmeans=True, showextrema=True, showmedians=True)
            ax.set_title(method)
            ax.set_xticks([])
            if i == 0:
                ax.set_ylabel("distances")

    if plot:
        plt.tight_layout()
        plt.show()

    return df_distance

if __name__ == "__main__":
    df = structureDistance("../results/eterna100_v2", plot=True)

