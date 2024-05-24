"""
     @author : nonosaha@mis.mpg.de/cyrillecardinale@gmail.com

"""
#import libraries
import numpy as np
import RNA
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb



def main() :

    targets = pd.read_csv('../../data/PseudoBase++/full_sorted_pk.csv')

    """
    gc_data = []
    for gc in [0.25, 0.5, 0.75, 1] :
        df = pd.read_csv('../../data/PseudoBase++/antaRNA/antaRNAGC'+str(gc)+'.csv')
        for line in df.values:
            gc_data +=[[line[1], line[2], line[3], RNA.bp_distance(line[1],line[3]), gc, len(line[1])]]

        print(gc)

    pd.DataFrame(gc_data, columns=["target", "sequence", "structure", "bd_disctance", "GC content", "Length"]).to_csv('../../data/PseudoBase++/antaRNA/ipknot_resultGC_content.csv')
    """

    #data = pd.read_csv('../../data/PseudoBase++/antaRNA/ipknot_resultGC_content.csv')
    data = pd.DataFrame([])
    plot_data = {}
    for gc in [25, 50, 75, 100] :
        plot_data[gc*0.01] = []
        dt = pd.read_csv('../../data/PseudoBase++/antaRNA/antaRNA_ipknotGC'+str(gc)+'.csv')
        data = data.append(dt)
        seqs = dt["sequence"].values.reshape(266, 20)
        gc_emp = []
        for sets in seqs :
            for sq in sets :
                nl = list(sq)
                gc_emp += [(nl.count("G")+nl.count("C"))/len(nl)]
        print(gc, np.mean(gc_emp))
        distances = dt["bp_distance"].values.reshape(266, 20)
        for line in distances :
            nl = list(line)
            plot_data[gc*0.01] +=[nl.count(0)*1.0/len(nl)]

    plt.boxplot(plot_data.values(), labels=plot_data.keys())
    plt.show()



    arnaque_data = {}
    for gc in [25, 50, 75, 100] :
        dt = pd.read_csv('../../data/PseudoBase++/aRNAque/aRNAque_ipknotGC'+str(gc)+'.csv')
        data = data.append(dt)
        arnaque_data [gc*0.01] = []

        distances = dt["bp_distance"].values.reshape(266, 20)
        seqs = dt["sequence"].values.reshape(266, 20)

        gc_emp = []
        for sets in seqs :
            for sq in sets :
                nl = list(sq)
                gc_emp += [(nl.count("G")+nl.count("C"))/len(nl)]
        print(gc, np.mean(gc_emp))
        for line in distances :
            nl = list(line)
            arnaque_data[gc*0.01] +=[nl.count(0)*1.0/len(nl)]

    plt.boxplot(arnaque_data.values(), labels=arnaque_data.keys())
    plt.show()

    print(data)
    data.to_csv("../../data/PseudoBase++/GC_content_data_ipknot.csv")

    figure = plt.figure(constrained_layout=True, figsize=(12,4))
    gs = figure.add_gridspec(nrows=1, ncols=3, left=0.05, right=0.48, wspace=0.05)
    ax = figure.add_subplot(gs[0,0])
    ax.set_title("(A)", fontsize=12)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    sb_bx = sb.boxplot(y='bp_distance', x='GC content', hue='Tool', data=data)
    sb_bx.set(ylabel='BP distance ', xlabel='GC content')
    ax = figure.add_subplot(gs[0,1])
    ax.set_title("(B)", fontsize=12)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    sb_bx = sb.boxplot(y='GC distance', x='GC content', hue='Tool', data=data)
    sb_bx.set(ylabel='GC-content distance', xlabel='GC content')
    plt.legend([],[], frameon=False)
    barplot_data = {
                    "aRNAque": {},
                    "antaRNA": {}
                    }


    for gc in [0.25, 0.5, 0.75, 1] :
        print((len(plot_data[gc])-plot_data[gc].count(0))*100/266.)
        df = data[data["GC content"]==gc]
        df_ar = df[df["Tool"]=="aRNAque"]
        barplot_data["aRNAque"][gc] = len(df_ar[df_ar["bp_distance"] ==0])

        df_ant = df[df["Tool"]=="antaRNA"]
        barplot_data["antaRNA"][gc] = len(df_ant[df_ant["bp_distance"] ==0])

    print(barplot_data)
    ax = figure.add_subplot(gs[0,2])
    ax.set_title("(C)", fontsize=12)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_xlabel("GC content")
    ax.set_ylabel("Number of successes")
    si = 0
    for key in barplot_data.keys():
        X_si = np.array(list(barplot_data[key].keys()))+si
        X = np.array(list(barplot_data[key].keys()))
        Y = list(barplot_data[key].values())
        if key == "aRNAque":
            plt.bar(X_si, Y, label=key, width=.1, alpha=1.,color="peru")
        else :
            plt.bar(X_si, Y, label=key, width=.1, alpha=1.,color="steelblue")


        for i in range(len(X)) :
            plt.annotate(str(Y[i]), xy=(X_si[i],Y[i]), ha='center', va='bottom')
            if key == "aRNAque":
                plt.annotate("("+str(len(arnaque_data[X[i]])-arnaque_data[X[i]].count(0))+")", xy=(X_si[i],Y[i]+220), ha='center', va='bottom')
            else :
                plt.annotate("("+str(len(plot_data[X[i]])-plot_data[X[i]].count(0))+")", xy=(X_si[i],Y[i]+220), ha='center', va='bottom')

        si +=0.1

    plt.savefig("../../images/PseudoBase++/fig6_new.pdf")
    plt.savefig("../../images/PseudoBase++/fig6_new.png")
    plt.show()



if __name__ == '__main__':
    main()
