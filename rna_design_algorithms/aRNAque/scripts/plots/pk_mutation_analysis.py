import pandas as pd
import numpy as np
import json
import os
import matplotlib.pyplot as plt
import RNA
import ast
import seaborn as sb



def main() :

    binomial = []
    levy = []
    all_gen_bino = []
    lengths = [25, 32,33,35,42, 52,56,61,63,64,66,72,73,76, 77, 78, 81, 82, 88,92,104,110,120,128,130,180]
    c_starts = []
    mu_starts = []

    pk_data = pd.read_csv("../../data/PseudoBase++/tuning_mutdata.csv")
    tuning_data = {target: [] for target in pk_data.values[:,-1]}
    #print(tuning_data)
    for l in lengths :
        with open("../../data/PseudoBase++/pk_tuning/pk__c"+str(l)+".out") as file_:
            data = []
            targets = []

            for line in file_.readlines() :
                if line.startswith("[(") :
                    data.append(ast.literal_eval(line))

                if line.startswith("('Solving") :
                    targets.append(ast.literal_eval(line)[-1])

            for t in set(targets) :
                tuning_data[t] = []


            print(len(data), len(targets))

            for i in range(len(targets)):
                gens = []
                dt = data[i]
                for j in range(len(dt)):
                    gens += [dt[j][-1]]
                tuning_data[targets[i]] += [gens]

    for key in tuning_data :
        if len(tuning_data[key]) > 0 :
            c = np.arange(1,2.1, step=0.1)
            med = [np.median(d) for d in tuning_data[key]]
            c_starts.append(c[med.index(min(med))])
            if 200>np.min(med)>=10 :
                levy.append([len(key),c[med.index(min(med))]])
                print(key,med)


    print(c_starts)

    tuning_data2 = {target: [] for target in pk_data.values[:,-1]}
    for l in lengths :
        with open("../../data/PseudoBase++/pk_tuning/pk__mu"+str(l)+".out") as file_:
            data = []
            targets = []

            for line in file_.readlines() :
                if line.startswith("[(") :
                    data.append(ast.literal_eval(line))

                if line.startswith("('Solving") :
                    targets.append(ast.literal_eval(line)[-1])

            for t in set(targets) :
                tuning_data2[t] = []


            print(len(data), len(targets))

            for i in range(len(targets)):
                gens = []
                dt = data[i]
                for j in range(len(dt)):
                    gens += [dt[j][-1]]
                tuning_data2[targets[i]] += [gens]
    bins = []
    for key in tuning_data2 :
        if len(tuning_data2[key]) > 0 :
            s=np.round(0.5/len(key), 4)
            mu = np.arange(s,.2+s, step=s)
            med = [np.median(d) for d in tuning_data2[key]]
            mu_starts.append(mu[med.index(min(med))])
            bins.append(mu[med.index(min(med))]*len(key))
            if mu[med.index(min(med))]<0.025 :
                binomial.append([len(key),mu[med.index(min(med))]])
                print(key,med, mu[med.index(min(med))])


    print(binomial)

    binomial = np.array(binomial)
    levy = np.array(levy)
    """
    figure = plt.figure(constrained_layout=True, figsize=(8,4))
    gs = figure.add_gridspec(nrows=1, ncols=2, left=0.05, right=0.48, wspace=0.05)
    ax = figure.add_subplot(gs[0,0])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.title("Binomial mutation")
    plt.xlabel("Target length")
    plt.ylabel(r"$\mu*$")
    plt.plot(binomial[:,0],binomial[:,1], 'o-', color='blue')
    ax = figure.add_subplot(gs[0,1])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.plot(levy[:,0], levy[:,1], 'o-', color='orange')
    plt.title("Levy mutation")
    plt.xlabel("Target length")
    plt.ylabel(r"$c*$")
    plt.savefig("../../images/PseudoBase++/tuning_param2.pdf")
    plt.show()
    plt.boxplot([np.array(c_starts)/sum(c_starts), np.array(mu_starts)/sum(mu_starts)], labels=[r'$c*$',r"$\mu*$"])
    plt.ylabel('Normalized distribution')
    plt.savefig("../../images/PseudoBase++/tuning_param.pdf")
    plt.show()

    figure = plt.figure(constrained_layout=True, figsize=(8,4))
    gs = figure.add_gridspec(nrows=1, ncols=2, left=0.05, right=0.48, wspace=0.05)
    ax = figure.add_subplot(gs[0,0])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.title("(A): Levy mutation")
    plt.xlabel(r"$c*$")
    ax.hist(c_starts, bins=20, label=r"$c*$", alpha=0.8)
    ax = figure.add_subplot(gs[0,1])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.title("(B): Binomial mutation")
    plt.xlabel(r"$\mu*$")
    plt.hist(mu_starts, color="darkorange", bins=10, alpha=0.8)
    plt.savefig("../../images/PseudoBase++/LevyVSBino_mutation_rate.pdf")
    plt.show()
    """

    pk_data_op = []
    with open('../../data/PseudoBase++/pk_ipnot_mu.out') as f :
        for line in f :
            if line.startswith('[(') :
                pk_data_op.append(line)
    ipknot_dataOP = []

    for elt in pk_data_op :
        ipknot_dataOP.append(ast.literal_eval(elt))

    plot_dataOp = []
    #sequence_
    for list_ in ipknot_dataOP :
        l = []

        for elt in list_ :
            l.append(elt[-1])

        if len(l)<20 :
            continue
        plot_dataOp.append(l)
    figure = plt.figure(constrained_layout=True, figsize=(12,4))
    gs = figure.add_gridspec(nrows=1, ncols=2, left=0.05, right=0.48, wspace=0.05)

    ax = figure.add_subplot(gs[0,1])
    ax.set_title(r"(B): Binamial mutation (histogram of $\mu*L$)", fontsize=12)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_ylabel("Number of target structures")
    plt.xlabel(r"$\mu*L$", fontsize=12)
    plt.hist(bins, color="darkorange", bins=50, alpha=0.8)

    """ax2 = plt.axes([0.65, 0.5, 0.3, 0.4])
    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    #plt.title("Tuning Binomial mutation rate ($\mu$)")
    ax2.set_ylabel('Median # of gens',fontsize=10)
    ax2.set_xlabel(r'Mutation rate ($\mu$)',fontsize=10)
    plt.plot([np.round(val,2) for val in np.arange(0.005,0.67, 0.005)],[np.median(dt) for dt in plot_dataOp],'--',color="darkorange", label='# of gens')
    plt.legend(loc='lower right',bbox_to_anchor=(1, 0.5))
    ax3 = ax2.twinx()
    ax3.spines["top"].set_visible(False)
    ax3.set_ylabel(r'Success rate (%)',fontsize=10)
    ax3.plot([np.round(val,2) for val in np.arange(0.005,0.67, 0.005)],[(len(dt)-dt.count(200))/0.5 for dt in plot_dataOp], color="darkorange", label="Success rate")
    plt.legend(loc='upper right', bbox_to_anchor=(1, 0.5))"""
    #plt.savefig("../../images/PseudoBase++/bino_tuning.pdf")
    #plt.show()
    """
    pk_data_levy = []
    with open('../../data/PseudoBase++/pk_ipnot_c.out') as f :
        for line in f :
            if line.startswith('[(') :
                pk_data_levy.append(line)

    ipknot_data = []

    for elt in pk_data_levy :
        ipknot_data.append(ast.literal_eval(elt))

    plot_dataLv = []
    sequence_dt = []
    for list_ in ipknot_data :
        l = []
        s = []

        for elt in list_ :
            l.append(elt[-1])
            s.append(elt[0][0])

        if len(l)<20 :
            continue
        plot_dataLv.append(l)
        sequence_dt.append(s)

    """
    #figure = plt.figure(constrained_layout=True, figsize=(8,4))
    #gs = figure.add_gridspec(nrows=1, ncols=1, left=0.05, right=0.48, wspace=0.05)

    ax = figure.add_subplot(gs[0,0])
    ax.set_title(r"(A): Lévy mutation (histogram of $c*$)", fontsize=12)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    #plt.title("(A): Levy mutation",fontsize=12)
    ax.set_ylabel("Number of target structures")
    plt.xlabel(r"$c*$", fontsize=12)
    ax.hist(c_starts, bins=20, label=r"$c*$", alpha=0.8, color="deepskyblue")

    """ax2 = plt.axes([0.15, 0.5, 0.3, 0.4])
    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)"""
    #plt.title("Tuning Binomial mutation rate ($\mu$)")
    #plt.title("Tuning Zipf's exponent ($s$)")
    """ax2.set_ylabel('Median # of gens',fontsize=10)
    ax2.set_xlabel(r'Zipf exponent ($c$)',fontsize=10)
    plt.plot([np.round(val,1) for val in np.arange(0.2,7.2,0.2)],[np.median(dt) for dt in plot_dataLv],'--', color='deepskyblue', label='# of gens')
    plt.legend(loc='lower right',bbox_to_anchor=(1, 0.5))
    ax3 = ax2.twinx()
    ax3.spines["top"].set_visible(False)
    ax3.set_ylabel(r'Success rate (%)',fontsize=10)
    ax3.plot([np.round(val,1) for val in np.arange(0.2,7.2,0.2)],[len(dt)-dt.count(200) for dt in plot_dataLv], color='deepskyblue', label="Success rate")
    plt.legend(loc='upper right', bbox_to_anchor=(1, 0.5))"""
    plt.savefig("../../images/PseudoBase++/mut_tuning.pdf")
    plt.show()

if __name__ == "__main__" :
    main()
