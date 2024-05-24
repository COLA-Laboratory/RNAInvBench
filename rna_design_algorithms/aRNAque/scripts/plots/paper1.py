import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ast
import seaborn as sb
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import matplotlib as mpl
import json
import scipy.stats as stat





def main() :
    #####Plot for the RFAM benchmark data set aRNAque vs. Nupack ############################
    """
    rfam_data = pd.read_csv('../../data/RFAM/benchmark_result_turner2004.csv')
    print(rfam_data)

    rnaevol_rfam = rfam_data[rfam_data['Tool']=='aRNAque']
    nupack_rfam = rfam_data[rfam_data['Tool']=='nupack']

    figure = plt.figure(constrained_layout=True, figsize=(11,7))
    gs = figure.add_gridspec(nrows=2, ncols=2, left=0.05, right=0.48, wspace=0.05)
    ax = figure.add_subplot(gs[0,0])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.xlabel('Target ID')
    plt.ylabel('Median MFE')
    plt.plot([np.median(rnaevol_rfam[rnaevol_rfam['id']==i]['mfe'].values) for i in range(30)], 'o')
    plt.plot([np.median(nupack_rfam[nupack_rfam['id']==i]['mfe'].values) for i in range(30)], 'o')
    #plt.plot([np.median(rnainverse_rfam[rnainverse_rfam['id']==i]['mfe'].values) for i in range(30)], 'o-')


    ax = figure.add_subplot(gs[0,1])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.xlabel('Target ID')
    plt.ylabel('Median ensemble defect')
    #plt.yscale('log')
    plt.plot([np.median(rnaevol_rfam[rnaevol_rfam['id']==i]['ed'].values) for i in range(30)], 'o')
    plt.plot([np.median(nupack_rfam[nupack_rfam['id']==i]['ed'].values) for i in range(30)], 'o')
    #plt.plot([np.median(rnainverse_rfam[rnainverse_rfam['id']==i]['ed'].values) for i in range(30)], 'o-')

    ax = figure.add_subplot(gs[1,0:])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    #sb.set(rc={'figure.figsize': (20., 8.27)})
    #plt.ylabel('Hamming Distance')
    ax.set_ylabel('Hamming Distance')
    sb.boxenplot(ax=ax, y='Hamming Distance', x='id', hue='Tool', data=rfam_data)
    plt.savefig('../../images/RFAM/rfam_NPvsaRNAque.pdf')
    plt.show()


    ### Plot on Eterna Benchmark dataset
    eterna_bench = pd.read_csv('../../data/Eterna100/V1/benchmark_result_turner2004.csv')
    print(eterna_bench)

    rnaevol_eterna = eterna_bench[eterna_bench['Tool']=='aRNAque']
    sentrna_eterna = eterna_bench[eterna_bench['Tool']=='sentrna']
    erd_eterna = eterna_bench[eterna_bench['Tool']=='erd']
    rnainverse_eterna = eterna_bench[eterna_bench['Tool']=='rnainverse']


    dt = [np.median(rnaevol_eterna[rnaevol_eterna['id']==i]['mfe'].values) for i in rnaevol_eterna['id'].values]
    target_ids = list(set(rnaevol_eterna['id'].values))


    figure = plt.figure(constrained_layout=True, figsize=(11,7))
    gs = figure.add_gridspec(nrows=2, ncols=2, left=0.05, right=0.48, wspace=0.05)
    ax = figure.add_subplot(gs[0,0])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.xlabel('Target ID')
    plt.ylabel('Median MFE')
    plt.plot(target_ids,[np.median(rnaevol_eterna[rnaevol_eterna['id']==i]['mfe'].values) for i in set(rnaevol_eterna['id'].values)], 'o')
    plt.plot(target_ids,[np.median(sentrna_eterna[sentrna_eterna['id']==i]['mfe'].values) for i in set(rnaevol_eterna['id'].values)], 'o')
    plt.plot(target_ids,[np.median(erd_eterna[erd_eterna['id']==i]['mfe'].values) for i in set(rnaevol_eterna['id'].values)], 'o')
    plt.plot(target_ids,[np.median(rnainverse_eterna[rnainverse_eterna['id']==i]['mfe'].values) for i in set(rnaevol_eterna['id'].values)], 'o')



    ax = figure.add_subplot(gs[0,1])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.xlabel('Target ID')
    plt.ylabel('Median ensemble defect')
    #plt.yscale('log')
    plt.ylim((-0.5,1))
    plt.plot(target_ids,[np.median(rnaevol_eterna[rnaevol_eterna['id']==i]['ed'].values) for i in set(rnaevol_eterna['id'].values)], 'o')
    plt.plot(target_ids,[np.median(sentrna_eterna[sentrna_eterna['id']==i]['ed'].values) for i in set(rnaevol_eterna['id'].values)], 'o')
    plt.plot(target_ids,[np.median(erd_eterna[erd_eterna['id']==i]['ed'].values) for i in set(rnaevol_eterna['id'].values)], 'o')
    plt.plot(target_ids,[np.median(rnainverse_eterna[rnainverse_eterna['id']==i]['ed'].values) for i in set(rnaevol_eterna['id'].values)], 'o')


    ax = figure.add_subplot(gs[1,0:])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    #sb.set(rc={'figure.figsize': (20., 8.27)})
    #plt.ylabel('Hamming Distance')
    ax.set_ylabel('Hamming Distance')
    sb.boxenplot(ax=ax, y='Hamming Distance', x='id', hue='Tool', data=eterna_bench)
    plt.savefig('../../images/Eterna100/V1/eterna_bench.pdf')
    plt.show()
    """

    ###### Plot on CPU time and quality of the designed sequences #####################

    dt = pd.read_csv('../../data/Eterna100/V1/time_data.csv')
    colors = mpl.cm.rainbow(np.linspace(0, 1, 100))
    figure = plt.figure(constrained_layout=True, figsize=(8,5))
    gs = figure.add_gridspec(nrows=1, ncols=1, left=0.05, right=0.48, wspace=0.05)

    ax = figure.add_subplot(gs[0,0])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.xlabel("aRNAque average CPU time")
    plt.ylabel('RNAinverse average CPU time')
    plt.xscale('log')
    plt.yscale('log')

    dt1 = dt[dt['rnainverse_ham']==0]
    fig= plt.scatter(dt1[dt1['aRNAque_ham']==0]['aRNAque_time'],dt1[dt1['aRNAque_ham']==0]['rnainverse_time'],c=colors[dt1[dt1['aRNAque_ham']==0]['aRNAque_time'].index], marker='o', label=r'$MHD_{aRNAque} = MHD_{RNAinverse} = 0$')

    dt2 = dt[dt['aRNAque_ham']==0]
    fig2=plt.scatter(dt2[dt2['rnainverse_ham']>0]['aRNAque_time'],dt2[dt2['rnainverse_ham']>0]['rnainverse_time'],c=colors[dt2[dt2['rnainverse_ham']>0]['aRNAque_time'].index], marker='+', label=r'$MHD_{aRNAque} = 0; MHD_{RNAinverse}>0$')

    dt3 = dt[dt['aRNAque_ham']>dt['rnainverse_ham']]
    plt.scatter(dt3['aRNAque_time'],dt3['rnainverse_time'],c=colors[dt3['aRNAque_time'].index], marker='_', label=r'$MHD_{aRNAque}<MHD{RNAinverse}$')

    plt.legend(loc='lower right')
    ax2 = plt.axes([0.16, 0.66, 0.4, 0.3])
    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    plt.yscale('log')
    ax2.plot(dt['Length'].values,dt['rnainverse_time'], label='RNAinverse')
    ax2.plot(dt['Length'].values,dt['aRNAque_time'], label='aRNAque')
    plt.ylabel(r'T($s$)')
    plt.xlabel('L')

    plt.legend()

    plt.savefig("../../images/Eterna100/V1/cpuvssuccess.pdf")
    plt.show()







if __name__== "__main__" :
    main()
