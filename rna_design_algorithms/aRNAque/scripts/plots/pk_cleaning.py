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

    nH = len(pkbase_peertype["H"])
    ncH = len(pkbase_peertype["cH"])
    nB = len(pkbase_peertype["B"])
    nK = len(pkbase_peertype["K"])
    print("nb of H = {}, nb of cH = {}, nb of B = {}, nb of K = {}".format(nH, ncH, nB, nK))
    target_type = {}
    #print(pkbase_peertype)
    sorted_targets = pd.read_csv("../../data/PseudoBase++/sorted_pk.csv").values[:,-1]

    target_df = pd.read_csv("../../data/PseudoBase++/pkbase.csv")
    targets = target_df.values[:,-1]
    print(len)
    for key in targets :
        if key in pkbase_peertype['H'] :
            target_type[key] = 'H-type'
        elif key in pkbase_peertype['K'] :
            target_type[key] = 'K-type'
        elif key in pkbase_peertype['B'] :
            target_type[key] = 'B-type'
        elif key in pkbase_peertype['cH'] :
            target_type[key] = 'cH-type'
        else :
            print(key)


    print("H-type",list(target_type.values()).count("H-type"))
    print("K-type",list(target_type.values()).count("K-type"))
    print("B-type",list(target_type.values()).count("B-type"))
    print("cH-type",list(target_type.values( )).count("cH-type"))

    print(len(set(targets)))
    print(len(targets))

    tg_values = target_df.values.tolist()

    tg_values.sort(key=lambda elt: len(elt[-1]))

    print(pd.DataFrame(tg_values))

    pd.DataFrame(tg_values).to_csv("../../data/PseudoBase++/full_sorted_pk.csv", header=None, index=None)




if __name__ == '__main__':
    main()
