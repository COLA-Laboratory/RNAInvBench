"""
     @author : nonosaha@mis.mpg.de/cyrillecardinale@gmail.com

"""
#import libraries
import numpy as np
import RNA
import ast
import matplotlib.pyplot as plt
import scipy.stats as st
import pandas as pd
import os
import json
import seaborn as sb


def entropy_seq(sequences) :
    p_k = []
    length = len(sequences[0])
    for n in ['A','C','G','U'] :
        p_n = []
        for i in range(length) :

            p = 0.
            for seq in sequences :
                if seq[i] == n :
                    p += 1.
            p_n.append(p/len(sequences))
        p_k += [p_n]
    print(st.entropy(p_k, axis=0))
    return sum(st.entropy(p_k, axis=0))


def main() :


    with open("../../data/Eterna100/V1/eterna_ds.json") as js_file :
        eterna_density = json.load(js_file)
        js_file.close()

    eterna_ds = {
        'Low' : [],
        'High': [],
    }
    eterna_ds_Y = [0., 0.]
    for key in eterna_density.keys() :
        if float(key) <= 0.5 :
            eterna_ds_Y[0] += len(eterna_density[key])
            eterna_ds['Low'] += eterna_density[key]
        else :
            eterna_ds_Y[1] += len(eterna_density[key])
            eterna_ds['High'] += eterna_density[key]

    pk_ds = {}

    pk_df = pd.read_csv("../../data/PseudoBase++/pkbase.csv")
    for tg in pk_df.values[:,-1] :
        pk_ds[tg] = np.round((tg.count('(')*2. + tg.count('[')*2. + tg.count('{')*2.)/len(tg),1)
    pk_ds_X = set(list(pk_ds.values()))

    pkbp_ds = list(pk_ds.values())

    pk_ds_Y = [0., 0.]
    low_length = []
    high_length = []

    for key in pk_ds.keys() :

        if pk_ds[key] <= 0.5 :
            pk_ds_Y[0] += 1
            low_length += [len(key)]
        else :
            pk_ds_Y[-1] += 1
            high_length += [len(key)]

    eterna_1999GC2 = pd.read_csv('../../data/Eterna100/V1/eterna1999_zipf_GC2.csv').values[:,-1].reshape(100,5)
    targets = []

    print("PK low dens. mean length = ", np.median(low_length), len(low_length))
    print("PK high dens. mean length = ", np.median(high_length), len(high_length))
    #plt.boxplot([low_length, high_length], labels=["Low", "High"])
    #plt.show()
    with open("../../data/Eterna100/V1/eterna_1999_zipfGC2.log") as file_ :
        lines = file_.readlines()

        for line in lines :
            if line.startswith("('Solving") :
                targets += [ast.literal_eval(line)[-1]]

    eterna199_levy_ds = {
        'Low' : [],
        'High' : []
    }

    for i,t in enumerate(targets) :
        if t in eterna_ds['Low'] :
            eterna199_levy_ds["Low"] += eterna_1999GC2[i].tolist()
        else:

            eterna199_levy_ds["High"] += eterna_1999GC2[i].tolist()

    plot_df = []
    for key in eterna199_levy_ds.keys() :
        for val in eterna199_levy_ds[key] :
            plot_df += [[key, val+1,"Levy"]]


    eterna_1999GC2_OP = pd.read_csv('../../data/Eterna100/V1/eterna1999op_GC2.csv').values[:,-1].reshape(100,5)

    eterna199_OP_ds = {
    'Low' : [],
    'High' : []
    }

    for i,t in enumerate(targets) :
        if t in eterna_ds['Low'] :
            eterna199_OP_ds["Low"] += eterna_1999GC2_OP[i].tolist()
        else:

            eterna199_OP_ds["High"] += eterna_1999GC2_OP[i].tolist()

    #print("Mean length of low dens. bp: ", np.median([len(tg) for tg in eterna_ds["low"]]))
    #print("Mean length of High dens. bp: ", np.median([len(tg) for tg in eterna_ds["high"]]))
    #plt.boxplot([[len(tg) for tg in eterna_ds["low"]], [len(tg) for tg in eterna_ds["high"]]], labels=["low", "high"])
    #plt.show()

    data_plot3 = []
    for key in eterna_ds.keys():
        for tg in eterna_ds[key] :
            data_plot3 += [[key, len(tg), "Eterna100"]]
    for l in low_length :
        data_plot3 += [["Low", l, "PseudoBase++"]]
    for l in high_length :
        data_plot3 += [["High", l, "PseudoBase++"]]

    data_plot3_df = pd.DataFrame(data_plot3, columns=["Base pair density", "Length", "Dataset"])

    for key in eterna199_OP_ds.keys() :
        for val in eterna199_OP_ds[key] :
            plot_df += [[key, val+1,"OP"]]


    df_plot = pd.DataFrame(plot_df, columns=["Base pair density","Generation ($t$)","Mutation type" ])
    figure = plt.figure(constrained_layout=True, figsize=(11,4))
    gs = figure.add_gridspec(nrows=1, ncols=3, left=0.05, right=0.48, wspace=0.05)
    ax = figure.add_subplot(gs[0,0])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    ax.set_xlabel("Base pair density", fontsize=12)
    ax.set_ylabel(r"Generation ($t$)", fontsize=12)
    ax.set(yscale="log")
    plt.title("(A)", fontsize=15)
    sb_bx = sb.boxplot(data=df_plot, x="Base pair density", y=r"Generation ($t$)", hue="Mutation type",ax=ax, palette={
        "Levy": "deepskyblue",
        "OP": "darkorange"
    })
    handles, _ = sb_bx.get_legend_handles_labels()          # Get the artists.
    sb_bx.legend(handles, ["Lévy mutation", "Local mutation"], loc="upper left") #
    #plt.legend(loc="upper left")
    ax = figure.add_subplot(gs[0,1])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.title("(B)", fontsize=15)
    plt.xlabel('Base pair density', fontsize=12)
    plt.ylabel('# of targets', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.bar(["Low","High"], height=pk_ds_Y, align='center', width=0.3, label='Pseudobase++',color='wheat')
    for i in range(2):
        plt.annotate(str(np.round(pk_ds_Y[i]*100/sum(pk_ds_Y), 2))+"%", xy=(i,pk_ds_Y[i]-((pk_ds_Y[i]*100/sum(pk_ds_Y))-13)), ha='center', va='bottom')
    plt.bar(["Low","High"], height=eterna_ds_Y, align='center', width=0.3, label='Eterna100', color='tan')
    for i in range(2):
        plt.annotate(str(np.round(eterna_ds_Y[i]*100/sum(eterna_ds_Y), 2))+"%", xy=(i,eterna_ds_Y[i]-35), ha='center', va='bottom')

    ax = figure.add_subplot(gs[0,2])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    ax.set_xlabel("Base pair density", fontsize=12)
    ax.set_ylabel(" Target length", fontsize=12)
    plt.title("(C)", fontsize=15)
    sb_bx = sb.boxplot(data=data_plot3_df, x="Base pair density", y="Length", hue="Dataset",ax=ax,palette={
        "Eterna100": "tan",
        "PseudoBase++": "wheat"
    })
    plt.legend(loc="upper left")
    plt.savefig('../../images/levy_analysis.pdf')
    plt.show()

if __name__ == '__main__' :
    main()
