import pandas as pd
import numpy as np
import json
import os
import matplotlib.pyplot as plt
import RNA
import ast
import seaborn as sb





def main() :
    with open("../../data/Eterna100/V1/mut_analysis/eterna_2004_V2") as file_ :
        data = file_.readlines()
        file_.close()

    all_data = []
    for dt in data :
        all_data += [ast.literal_eval(dt)]

    all_data = np.array(all_data)
    print(all_data.shape)
    zipf_data = all_data.reshape(11,12*100*2)

    clean_data  = []
    for target_data in zipf_data :
        gens = []
        for dt in target_data.reshape(1200,2) :
            gens.append(dt[-1])
        clean_data.append(gens)

    df_data = []
    target = 1
    for gens in clean_data :
        c = 1.
        reshaped_data = np.array(gens).reshape(12,100)
        for gen in reshaped_data :
            for g in gen :
                df_data +=[[g,c,target]]
            c += 0.5
        target +=1


    df_ = pd.DataFrame(df_data, columns=[r'Generation ($t$)', r'Zipf exponent ($s$)', 'Target ID'])

    #figure = plt.figure(constrained_layout=True, figsize=(11,5))
    #gs = figure.add_gridspec(nrows=1, ncols=1, left=0.05, right=0.48, wspace=0.05)
    #ax = figure.add_subplot(gs[0,0])
    #ax.spines["right"].set_visible(False)
    #ax.spines["top"].set_visible(False)
    #sb.boxplot(data=df_, hue=r"Zipf exponent ($s$)", x="Target ID", y=r"Generation ($t$)")
    #plt.show()
    #plt.savefig('../images/zipf_exponent_11.pdf')

    lengths = []
    for dt in zipf_data :
        lengths.append(len(dt[0][0]))
    c_zipf = []
    all_gens_levy = []
    for t in range(1,12) :
        med = df_[df_["Target ID"]==t][r"Generation ($t$)"].values.reshape(12,100)
        d = [np.median(m) for m in med]
        d_start = min(d)
        c = np.arange(1, 7, 0.5)
        if d.count(d_start) > 1 :
            print(d)
        print(np.min(d), d.index(d_start),c[d.index(d_start)])
        all_gens_levy += [med[d.index(d_start)]]
        c_zipf.append(c[d.index(d_start)])

    #plt.yscale("log")
    plt.plot(lengths,c_zipf, 'o-', label=r'$c*$')
    #plt.show()
    binomial = []
    all_gen_bino = []
    for l in lengths :
        with open("../../data/ETerna100/V1/mut_analysis/eterna__mu"+str(l)+".out") as file_:
            data = []
            for line in file_.readlines() :
                if line.startswith("[(") :
                    data.append(ast.literal_eval(line))
            meds = []
            all_gens = []
            for entry in data :
                gens = []
                for elt in entry :
                    if elt != None :
                        #print(elt)
                        gens += [elt[-1]]
                meds += [np.median(gens)]
                all_gens += [gens]

            s =np.around(0.5/l,4)
            mu_range = [np.round(mu, 4) for mu in np.arange(s,0.2+s, step=s )]
            if l == 35 :
                binomial += [mu_range[-2]]
                all_gen_bino += [all_gens[-2]]
                continue
            binomial += [mu_range[meds.index(min(meds))]]
            all_gen_bino += [all_gens[meds.index(min(meds))]]

            print(l, len(data), meds, len(mu_range), meds.index(min(meds)), mu_range[meds.index(min(meds))])

    #plt.yscale("log")
    plt.plot(lengths, binomial, 'o-', label=r'$\mu*$')
    plt.xlabel(r"Lenght $((\mu, c) *L)$")
    plt.legend()
    plt.savefig("../../images/mut_tuning2.pdf")
    plt.show()

    plot_d = []
    for i in range(len(lengths)) :
        for g in all_gen_bino[i] :
            plot_d += [[lengths[i], g, r'$\mu*$']]
        for g in all_gens_levy[i] :
            plot_d += [[lengths[i], g, r'$c*$']]

    df_plot = pd.DataFrame(plot_d, columns=["Length", "Generation", "Best parameter"])
    print(df_plot)
    for elt in all_gen_bino :
        print(len(elt))

    figure = plt.figure(constrained_layout=True, figsize=(11,5))
    gs = figure.add_gridspec(nrows=1, ncols=1, left=0.05, right=0.48, wspace=0.05)
    ax = figure.add_subplot(gs[0,0])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set(yscale="log")
    sb.boxplot(ax=ax, data=df_plot, hue=r"Best parameter", x="Length", y=r"Generation")
    plt.savefig("../../images/mut_tuning.pdf")
    plt.show()


if __name__ == "__main__":
    main()
