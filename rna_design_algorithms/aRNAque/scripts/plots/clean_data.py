"""
     @author : nonosaha@mis.mpg.de/cyrillecardinale@gmail.com

"""
#import libraries
import numpy as np
import RNA
import ast
import matplotlib.pyplot as plt
import scipy.stats as st
from pandas import read_csv
import os
import json


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

    diversities = {
    'levy' : {},
    'op' : {}
    }
    mean_fitness = {
    'levy' : {},
    'op' : {}
    }

    max_fitness = {
    'levy' : {},
    'op' : {}
    }

    distinct_seq = {
    'levy' : {},
    'op' : {}
    }
    distinct_strucs = {
    'levy' : {},
    'op' : {}
    }

    infos = []
    mean_fit2 = []
    for f in range(10) :
        print("folder ", f)
        fit_values  =[]
        root_path = "../../data/Log/Levy/"+str(f)

        diversities['levy'][f] = []
        mean_fitness['levy'][f] = []
        max_fitness['levy'][f] = []

        for i in range(len(os.listdir(root_path))-1) :
            pop = read_csv(root_path+'/gen'+str(i)+'.csv')
            pop1= read_csv(root_path+'/gen'+str(i)+'.csv')
            seqs = pop['RNA_sequence'].values
            diversities['levy'][f] += [entropy_seq(list(set(seqs)))]
            mean_fitness['levy'][f] += [np.mean(pop['Fitness'].values)]
            max_fitness['levy'][f] += [np.max(pop['Fitness'].values)]

        seq = set()
        strucs = set()
        distinct_s = []
        distinct_ss = []
        for i in range(1,len(os.listdir(root_path))-1) :
            pop = read_csv(root_path+'/gen'+str(i-1)+'.csv')
            pop1= read_csv(root_path+'/gen'+str(i)+'.csv')

            seq = set.union(seq, set(pop['RNA_sequence'].values.tolist()+pop1['RNA_sequence'].values.tolist()))
            strucs = set.union(strucs, set(pop['RNA_structure'].values.tolist()+pop1['RNA_structure'].values.tolist()))
            distinct_s += [len(seq)]
            distinct_ss += [len(strucs)]
        distinct_seq['levy'][f] = distinct_s
        distinct_strucs['levy'][f] = distinct_ss

        ##Compute the distinct sequences/ structures over time for the binomial_mutation.
        seq = set()
        strucs = set()
        distinct_s = []
        distinct_ss = []
        root_path = "../../data/Log/OP/"+str(f)
        for i in range(1,len(os.listdir(root_path))-1) :
            pop = read_csv(root_path+'/gen'+str(i-1)+'.csv')
            pop1= read_csv(root_path+'/gen'+str(i)+'.csv')

            seq = set.union(seq, set(pop['RNA_sequence'].values.tolist()+pop1['RNA_sequence'].values.tolist()))
            strucs = set.union(strucs, set(pop['RNA_structure'].values.tolist()+pop1['RNA_structure'].values.tolist()))
            distinct_s += [len(seq)]
            distinct_ss += [len(strucs)]

        distinct_seq['op'][f] = distinct_s
        distinct_strucs['op'][f] = distinct_ss
        diversities['op'][f] = []
        mean_fitness['op'][f] = []
        max_fitness['op'][f] = []

        for i in range(len(os.listdir("../../data/Log/OP/"+str(f)))-1) :
            pop = read_csv("../../data/Log/OP/"+str(f)+'/gen'+str(i)+'.csv')
            seqs = pop['RNA_sequence'].values
            diversities['op'][f] += [entropy_seq(list(set(seqs)))]
            mean_fitness['op'][f] += [np.mean(pop['Fitness'].values)]
            max_fitness['op'][f] += [np.max(pop['Fitness'].values)]



    names = ["distinct_stuctures.json", "distinct_sequences.json", "mean_fitness.json", "max_fitness.json", "diversity.json"]

    with open("../../data/diversity/distinct_stuctures.json", "w") as file_ :
        json.dump(distinct_strucs, file_)
        file_.close()

    with open("../../data/diversity/distinct_sequences.json", "w") as file_ :
        json.dump(distinct_seq, file_)
        file_.close()

    with open("../../data/diversity/entropies.json", "w") as file_ :
        json.dump(diversities, file_)
        file_.close()

    with open("../../data/diversity/mean_fitness.json", "w") as file_ :
        json.dump(mean_fitness, file_)
        file_.close()

    with open("../../data/diversity/max_fitness.json", "w") as file_ :
        json.dump(max_fitness, file_)
        file_.close()


if __name__ == '__main__' :
    main()
