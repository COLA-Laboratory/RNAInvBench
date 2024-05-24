import pandas as pd
import numpy
import RNA
import matplotlib.pyplot as plt
import subprocess
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split


def plot_result(eterna_solved, eterna_df, title):
    lengths = eterna_df["Length"].values
    bp_density = eterna_df["bp_density"].values
    all_maker = eterna_df["Eterna_id"].values
    for i,strc in enumerate(all_maker):
        x = lengths[i]
        y = bp_density[i]

        if strc in eterna_solved :
            c = 'green'
        else :
            c = 'red'
        if strc == eterna_solved[1]:
            plt.scatter(x, y, c=c, label ="solved")
        elif i==98:
            plt.scatter(x, y, c=c, label ='unsolved')
        else :
            plt.scatter(x, y, c=c)
        #plt.text(x,y,strc)
    plt.plot( lengths,numpy.ones(len(lengths))*numpy.median(bp_density),'b')
    plt.plot( numpy.ones(len(lengths))*numpy.median(lengths),bp_density, 'b')
    plt.xlabel("Length", fontweight="bold")
    plt.ylabel("bp density ", fontweight="bold")
    plt.title(title + ": " +str(len(eterna_solved))+"%")
    plt.legend()
    plt.grid(True)




def main() :

    df = pd.read_csv("../../data/Eterna100/V1/sorted-full-eterna-dt.csv")
    eterna_strc = []
    for line in df.values:
        eterna_strc.append([line[-1].strip(), line[1]])

    lengths = []
    bp_density = []
    all_maker = []
    for line in df.values :
        lengths.append(len(line[-1]))
        list_ = list(line[-1])
        bp_density.append(list_.count("(")*2/float(len(line[-1])))
        all_maker.append(line[-2])

    eterna_df = pd.DataFrame(data=numpy.array([lengths, bp_density, all_maker]).T, columns=['Length', 'bp_density', 'Eterna_id'])
    rnaevol = []
    result = pd.read_csv('../../data/Eterna100/V1/all_eterna.csv')
    arnaque_df = result[result["Tool"]=="rnainverse"]
    solved = []
    solved_arnaque = []

    df1 = pd.read_csv('../../data/Eterna100/V1/result-eterna-med.csv')

    for line in df1.values :
        if type(line[2]) == type("nono"):
            print(line[1])
            solved_arnaque.append(line[1])
    print(set(arnaque_df["id"].values), len(set(arnaque_df["id"].values)))
    for id in set(arnaque_df["id"].values):
        rst = arnaque_df[arnaque_df["id"]==id]["Hamming Distance"].values

        if 0.0 in rst :
            if 'NON' not in arnaque_df[arnaque_df["id"]==id]["sequences"].values :
                solved += [id]
    if len(solved) > 0 :
        plot_result(solved,eterna_df, "RNAinverse Performance on Eterna100")
        #plt.savefig("../../images/Eterna100/V1/difficulty/RNAinverse.pdf")
        plt.show()






if __name__ =="__main__" :
    main()
