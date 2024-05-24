import pandas as pd
import numpy as np
import json
import os
import matplotlib.pyplot as plt
import RNA
import ast
import pickle
import subprocess
import seaborn as sb
from pyfaidx import Fasta

def get_bp_position(structure) :
    pk_pairs,cbrk, bp, nbp= [] , [], [], []

    pairs = {
        'bp' : [],
        'pk' : [],
        'nbp' : []
    }
    for i,elt in enumerate(structure) :
        if elt =='[' :
            pk_pairs.append(i)
        elif elt ==']' :
            pairs['pk'] += [(pk_pairs.pop(),i)]
        elif elt == '(' :
            bp += [i]
        elif elt == ')' :
            pairs['bp'] += [(bp.pop(),i)]
        elif elt == '{' :
            cbrk += [i]
        elif elt == '}' :
            pairs['pk'] += [(cbrk.pop(),i)]
        else :
            pairs['nbp'] += [i]
    return pairs

def bp_distance(target, strc2) :
    pos1 = get_bp_position(target)
    pair1 = set(pos1['bp'] + pos1['pk'])
    pos2 = get_bp_position(strc2)
    pair2 = set(pos2['bp'] + pos2['pk'])

    return len(pair1) + len(pair2) - 2*len(pair1.intersection(pair2))





def main() :

    data_pkbase = Fasta('../../data/PseudoBase++/Pseudobase++_TEST')

    pkbase_peertype = {}
    pkbase_peertype ["H"] = []
    pkbase_peertype ["cH"] = []
    pkbase_peertype ["B"] = []
    pkbase_peertype ["K"] = []

    for item in data_pkbase.items() :
        if 'cH-type' in item[0] :
            print (item)
            seq = list(item[1])[0].__dict__
            pkbase_peertype["cH"].append(seq['seq'])
        elif 'H-type' in item[0] :
            seq = list(item[1])[0].__dict__
            pkbase_peertype["H"].append(seq['seq'])
        elif 'B-type' in item[0] :
            seq = list(item[1])[0].__dict__
            pkbase_peertype["B"].append(seq['seq'])
        elif 'K-type' in item[0] :
            print('===========================')
            print(item[0])
            seq = list(item[1])[0].__dict__
            pkbase_peertype["K"].append(seq['seq'])


    target_type = {}

    targets = pd.read_csv("../../data/PseudoBase++/sorted_pk.csv").values[:,-1]


    for key in targets :
        if key in pkbase_peertype['H'] :
            target_type[key] = 'H-type'
        if key in pkbase_peertype['K'] :
            target_type[key] = 'K-type'
        if key in pkbase_peertype['B'] :
            target_type[key] = 'B-type'
        if key in pkbase_peertype['cH'] :
            target_type[key] = 'cH-type'

    print(len(target_type))



    ###############Plots for the pseudoknots benchmark using ipknots###########################
    ipknot_df = pd.read_csv("../../data/PseudoBase++/benchmark_result_ipknot.csv")

    ###################################Ploting code#############################################
    figure = plt.figure(constrained_layout=True, figsize=(12,6))
    gs = figure.add_gridspec(nrows=2, ncols=2, left=0.05, right=0.48, wspace=0.05)
    ax = figure.add_subplot(gs[0,0])

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    ax.set_xlabel("PK-type", fontsize=12)
    ax.set_ylabel("BP distance", fontsize=12)

    lengths = list(set(ipknot_df ['Length'].values))
    ant_ipknot = ipknot_df[ipknot_df['Tools'] == 'antaRNA']
    levy_ipknot = ipknot_df[ipknot_df['Tools'] == 'aRNAque']

    length_plots = []
    data_ant = []
    data_levy = []

    for i in range(len(lengths)) :
        data_ant += [ant_ipknot[ant_ipknot['Length']==lengths[i]]['BP distance'].values]
        length_plots.append(lengths[i])
        data_levy += [levy_ipknot[levy_ipknot['Length']==lengths[i]]['BP distance'].values]



    plt.title("(A) IPknot", fontsize=15)
    sb.boxplot(ax=ax, y='BP distance', x='PK-type', hue='Tools', data=ipknot_df)
    plt.legend([],[], frameon=False)
    ax2 = figure.add_subplot(gs[1,0])
    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax2.set_xlabel('Length (L)', fontsize=12)
    ax2.set_ylabel('Mean BP distance',  fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.title("(C) IPknot", fontsize=15)
    #ax2.set(yscale="log")
    ax2.scatter(length_plots, [np.mean(dt) for dt in data_ant], label='antaRNA')
    ax2.scatter(length_plots, [np.mean(dt) for dt in data_levy], label='aRNAque')
    #plt.savefig('../../images/PseudoBase++/pk_ipknotaRNAqueVSantaRNA_2.pdf')
    #plt.show()

    ##############################Loading and cleaning data for pkbase benchmark with hotknots##########################

    hotknot_df = pd.read_csv("../../data/PseudoBase++/benchmark_result_hotknot.csv")
    lengths = list(set(hotknot_df ['Length'].values))
    ant_hotknot = hotknot_df[hotknot_df['Tools'] == 'antaRNA']
    levy_hotknot = hotknot_df[hotknot_df['Tools'] == 'aRNAque']

    length_plots = []
    data_ant = []
    data_levy = []

    for i in range(len(lengths)) :
        data_ant += [ant_hotknot[ant_hotknot['Length']==lengths[i]]['BP distance'].values]
        length_plots.append(lengths[i])
        levy_row = levy_hotknot[levy_hotknot['Length']==lengths[i]]['BP distance'].values
        if len(levy_row) == 0 :
            print(lengths[i])
        data_levy += [levy_row]

    #figure = plt.figure(constrained_layout=True, figsize=(9,6))
    #gs = figure.add_gridspec(nrows=2, ncols=1, left=0.05, right=0.48, wspace=0.05)
    ax = figure.add_subplot(gs[0,1])

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_ylabel('BP distance', fontsize=12)
    ax.set_xlabel('PK-type', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.title("(B) HotKnots", fontsize=15)

    #plt.title("Distribution BP distance to the target (20 runs peer target)")
    sb.boxplot(ax=ax, y='BP distance', x='PK-type', hue='Tools', data=hotknot_df)
    ax2 = figure.add_subplot(gs[1,1])
    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax2.set_xlabel(r'Length ($L$)', fontsize=12)
    ax2.set_ylabel('Mean BP distance distribution',fontsize=12)
    print(len(data_ant), len(data_levy))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.title("(D) HotKnots", fontsize=15)
    #ax2.set(yscale="log")
    ax2.scatter(length_plots, [np.mean(dt) for dt in data_ant],  label='antaRNA')
    ax2.scatter(length_plots, [np.mean(dt) for dt in data_levy], label='aRNAque')
    #plt.legend()
    plt.savefig('../../images/PseudoBase++/pk_hotknot_IPknot_aRNAqueVSantaRNA_2.pdf')
    plt.show()


if __name__ == "__main__" :

    main()
