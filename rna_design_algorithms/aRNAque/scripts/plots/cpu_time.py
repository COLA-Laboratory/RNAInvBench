"""
     @author : nonosaha@mis.mpg.de/cyrillecardinale@gmail.com
"""

#import libraries
import numpy as np
import RNA
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import matplotlib as mpl

# first part with least squares
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression


def f_model(x, n):
    return np.exp(x**n)

def main() :
    """
    antaRNA_data = pd.read_csv('../../data/PseudoBase++/antaRNA/antaRNA_ipknotGC25.csv')

    antaRNA_cputime = pd.read_csv('../../data/PseudoBase++/antaRNA/antaRNA_Hotknots_GC50.csv')

    hotknot_df = pd.read_csv("../../data/PseudoBase++/benchmark_result_hotknot.csv")

    print(antaRNA_cputime)
    #plt.plot(antaRNA_cputime["CPU time"].values[:-1])
    #plt.show()


    aRNAque_df = hotknot_df[hotknot_df['Tools']=="aRNAque"]

    time_data = []
    for target in aRNAque_df["Target"].values :
        time_data += antaRNA_cputime[antaRNA_cputime["target"]==target].values.tolist()

    print(len(aRNAque_df.values), len(time_data))
    print(time_data[0], aRNAque_df.values[0])
    print(len(aRNAque_df.values[0]), len(time_data[0]))

    time_clean_data = []

    for line in time_data :
        print(len(line))
        time_clean_data +=[[line[1], line[2],line[3], line[4], line[5], line[7],line[8], line[-1]]]

    for line in aRNAque_df.values :
        time_clean_data +=[[line[5], line[1],line[2], line[3], 0.5, line[6],line[7], line[-1]]]

    time_df = pd.DataFrame(time_clean_data, columns=["target","sequence", "structure", "BP distance", "GC content", "Length","Tool", "CPU time"])

    cpu_times = []
    for target in set(time_df["target"].values) :
        dt = time_df [time_df["target"]==target]
        dt_arq =  dt[dt["Tool"]=="aRNAque"]
        dt_ant =  dt[dt["Tool"]=="antaRNA"]
        if len(dt_arq) > 0 and len(dt_ant) >0 :
            cpu_times +=[[target, np.median(dt_arq["CPU time"].values), np.median(dt_ant["CPU time"].values), np.median(dt_arq["BP distance"].values), np.median(dt_ant["BP distance"].values), len(target)]]
    cpu_times.sort(key=lambda elt: elt[-1])
    df = pd.DataFrame(cpu_times, columns=["Target", "aRNAque_time", "antaRNA_time", "aRNAque_bpd", "antaRNA_bpd", "Length"])
    df.to_csv("../../data/PseudoBase++/CPU_time_Hotknots.csv")
    """
    #bp_dists = aRNAque_df['BP distance'].values.reshape(253, 20)
    #times = aRNAque_df['Time'].values.reshape(253, 20)
    df = pd.read_csv("../../data/PseudoBase++/CPU_time_Hotknots.csv")
    Y_arnaque = np.array(df["aRNAque_time"].values.tolist(), dtype=float)
    X = np.array(df["Length"].values.tolist(), dtype=float).reshape(-1, 1)

    model = LinearRegression()
    model.fit(X, Y_arnaque)
    print("coefficient of determination", model.score(X, Y_arnaque))
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(X,Y_arnaque, 'yo', X, model.predict(X), '--k')
    plt.show()
    colors = mpl.cm.rainbow(np.linspace(0, 1, len(df)))
    figure = plt.figure(constrained_layout=True, figsize=(8,5))
    gs = figure.add_gridspec(nrows=1, ncols=1, left=0.05, right=0.48, wspace=0.05)

    ax = figure.add_subplot(gs[0,0])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.xlabel("aRNAque average CPU time")
    plt.ylabel('antaRNA average CPU time')
    plt.xscale('log')
    plt.yscale('log')

    dt1 = df[df['antaRNA_bpd']==0]
    fig= plt.scatter(dt1[dt1['aRNAque_bpd']==0]['aRNAque_time'],dt1[dt1['aRNAque_bpd']==0]['antaRNA_time'],c=colors[dt1[dt1['aRNAque_bpd']==0]['aRNAque_time'].index], marker='o', label=r'$BP_{aRNAque} = BP_{antaRNA} = 0$')

    dt2 = df[df['aRNAque_bpd']==0]
    fig2=plt.scatter(dt2[dt2['antaRNA_bpd']>0]['aRNAque_time'],dt2[dt2['antaRNA_bpd']>0]['antaRNA_time'],c=colors[dt2[dt2['antaRNA_bpd']>0]['aRNAque_time'].index], marker='+', label=r'$BP_{aRNAque} = 0; BP_{antaRNA}>0$')

    dt3 = df[df['aRNAque_bpd']>df['antaRNA_bpd']]
    print(len(dt1), len(dt2), len(dt3), len(dt2[dt2['antaRNA_bpd']>0]), len(dt1[dt1['aRNAque_bpd']==0]), len(dt3[dt3["antaRNA_bpd"]==0]))
    plt.scatter(dt3['aRNAque_time'],dt3['antaRNA_time'],c=colors[dt3['aRNAque_time'].index], marker='_', label=r'$BP_{aRNAque}>BP_{antaRNA}$')

    plt.legend(loc='upper left')
    ax2 = plt.axes([0.5, 0.2, 0.4, 0.3])
    #ax2.spines["right"].set_visible(False)
    #ax2.spines["top"].set_visible(False)
    plt.yscale('log')
    ax2.plot(df['Length'].values,df['antaRNA_time'], "o", label='antaRNA')
    ax2.plot(df['Length'].values,df['aRNAque_time'],"o", label='aRNAque')
    plt.ylabel(r'CPU time')
    plt.xlabel('Length')
    plt.legend()
    plt.savefig("../../images/PseudoBase++/CPU_time_Hotknots2.pdf")
    plt.savefig("../../images/PseudoBase++/CPU_time_Hotknots2.png")
    plt.show()



if __name__ == '__main__':
    main()
