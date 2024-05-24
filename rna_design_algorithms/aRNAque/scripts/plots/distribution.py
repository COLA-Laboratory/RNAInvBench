from matplotlib import pyplot as plt
import numpy as np
import ast

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, SUPPRESS



def nthHarmonic(N,s) :

    harmonic = 1.00
    for i in range(2, N + 1) :
        harmonic += 1 / i**s

    return harmonic

def zipf_rvgen(low, high, N, size, a=6.) :
    choices = np.array(range(low, high+1))
    probs = choices**(-a) / nthHarmonic(N,a)
    return np.random.choice(choices,p=np.array(probs)/sum(probs),size=size)

def bino_gen(mu, length, size) :
    rst = [list(np.random.binomial(1, mu, length)).count(1) for i in range(size)]
    return rst

def parseArguments():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter, argument_default=SUPPRESS)
    parser.add_argument("--length", '-l' ,type=int, help= "Lenght of the Sequence")
    parser.add_argument('-c', type=float, default=1.5, help='Levy of zipf law exponent')
    parser.add_argument('-ps', type=int, default=100, help='Population size')
    return parser.parse_args()

def main() :
    args = parseArguments()
    zipf_dist = list(zipf_rvgen(1,args.length, args.length, args.ps, args.c))

    barplot_data = []

    for i in set(zipf_dist) :
        barplot_data +=[[i,zipf_dist.count(i)]]


    barplot_data = np.array(barplot_data)
    p = barplot_data[:,1]/sum(barplot_data[:,1])
    mean_zipf = sum(p*barplot_data[:,0])
    print(f"The mean of the Zipf's distribution is {mean_zipf}")
    mu = 0.12 #mean_zipf/args.length
    print(f"The Binomial mutation rate is {mu}")


    bino_dist = bino_gen(mu, args.length, args.ps)
    barplot_data_bino = []

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


    for i in set(bino_dist) :
        barplot_data_bino +=[[i,bino_dist.count(i)]]

    barplot_data_bino = np.array(barplot_data_bino)
    figure = plt.figure(constrained_layout=True, figsize=(10,2.5))
    gs = figure.add_gridspec(nrows=1, ncols=3, left=0.05, right=0.48, wspace=0.05)
    ax = figure.add_subplot(gs[0,0])
    ax.set_title("(A)")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    #plt.title("Binomial and Zipf distribution centered at the mean = "+str(np.round(mean_zipf,2)))

    plt.xlabel("Point mutations", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)
    plt.bar(barplot_data_bino[:,0], barplot_data_bino[:,1], width=.8, align='center', log=True, label=r"Binomial with $\mu*="+str(np.round(mu,3))+"$", color="darkorange")
    plt.bar(barplot_data[:,0], barplot_data[:,1], width=.8, align='center', log=True, label=r"Zipf with $c*="+str(args.c)+"$", alpha=0.5, color="deepskyblue")
    plt.legend()
    ax = figure.add_subplot(gs[0,1])
    ax.set_title("(B)")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_ylabel('Median # of gens',fontsize=10)
    ax.set_xlabel(r'$\mu$',fontsize=10)
    success_ratesB = [(len(dt)-dt.count(200))/0.5 for dt in plot_dataOp]
    med_Bino = [np.median(dt) for dt in plot_dataOp]
    mus = [np.round(val,2) for val in np.arange(0.005,0.67, 0.005)]

    ax.plot(mus[med_Bino.index(min(med_Bino))],min(med_Bino),'o', color='darkorange', ms=5)
    print("Best mu = ", mus[med_Bino.index(min(med_Bino))])
    ax.text(mus[med_Bino.index(min(med_Bino))]-0.05,min(med_Bino), r'$\mu*$')
    plt.plot(mus,[np.median(dt) for dt in plot_dataOp],'--',color="darkorange", label='# of gens')
    plt.legend(loc='lower right',bbox_to_anchor=(1, 0.5))
    ax3 = ax.twinx()
    ax3.spines["top"].set_visible(False)
    ax3.set_ylabel(r'Success (%)',fontsize=10)
    ax3.plot([np.round(val,2) for val in np.arange(0.005,0.67, 0.005)],[(len(dt)-dt.count(200))/0.5 for dt in plot_dataOp], color="darkorange", label="Success rate")
    plt.legend(loc='upper right', bbox_to_anchor=(1, 0.5))
    ax = figure.add_subplot(gs[0,2])
    ax.set_title("(C)")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_ylabel('Median # of gens',fontsize=10)
    c_expo = [np.round(val,1) for val in np.arange(0.2,7.2,0.2)]
    med_Levy = [np.median(dt) for dt in plot_dataLv]
    ax.set_xlabel(r'$c$',fontsize=10)
    plt.plot(c_expo,[np.median(dt) for dt in plot_dataLv],'--', color='deepskyblue', label='# of gens')
    success_ratesL = [len(dt)-dt.count(200) for dt in plot_dataLv]
    ax.plot(c_expo[success_ratesL.index(max(success_ratesL))],min(med_Levy),'o', color='deepskyblue', ms=5)
    ax.text(c_expo[success_ratesL.index(max(success_ratesL))]+0.1,min(med_Levy), r'$c*$')
    plt.legend(loc='lower right',bbox_to_anchor=(1, 0.5))
    ax3 = ax.twinx()
    ax3.spines["top"].set_visible(False)
    ax3.set_ylabel(r'Success (%)',fontsize=10)
    ax3.plot(c_expo,success_ratesL, color='deepskyblue', label="Success rate")

    plt.legend(loc='upper right', bbox_to_anchor=(1, 0.5))
    plt.savefig("../../images/illustration.pdf")
    plt.show()




if __name__ =="__main__" :
    main()
