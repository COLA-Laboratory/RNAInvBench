import sys
import subprocess
import ViennaRNA
from utils.sample import sample_incarnation
import pandas as pd

def call_incarnation(target_structure, constraints=None, tries=1, GCcontent=None):
    '''
    Use incarnation sampling as a start sequence of RNAinverse.

    Parameters:
    - target_structure (str): The target RNA secondary structure in dot-bracket notation.
    - constraints (str): The sequence constraints of results
    - tries(int): number of sequences you want to generate (same target structure)
    - GCcontent(float): design sequences' GC content you want to reach

    Returns:
    - sequences (list[str]): list of designed sequences(if tries > 1)
    - distances (list[int]): list of distances betweeen designed sequences' structure and target structure(if tries > 1)
    '''
    if not target_structure:
        print('Error: target structure is empty!')
        sys.exit(1)

    # sample
    init_seqs = sample_incarnation(target_structure, constraints, tries, GCcontent)

    # RNAinverse
    output_Seq = []  # designed sequences list
    output_Dis = []  # distances between designed seqs' structure and target structure
    if tries > 1:
        for i in range(tries):
            seq, dis = ViennaRNA.inverse_fold(init_seqs[i], target_structure)
            output_Seq.append(seq)
            output_Dis.append(dis)
    else:
        seq, dis = ViennaRNA.inverse_fold(init_seqs, target_structure)
        output_Seq = seq
        output_Dis = dis
    
    index = 0
    for x in range(len(output_Dis)):
        if output_Dis[x] < output_Dis[index]:
            index = x
    
    return output_Seq[index], output_Dis[index]

if __name__ == "__main__":
    target_structure = "(((....)))((...))"
    df = pd.read_csv("../Eterna100.txt")
    structures = df["str"].to_list()
    results = []
    for target_structure in structures:
        tries = 5
        output_seq, output_Dis = call_incarnation(target_structure, tries=tries)
    
        results.append([output_seq, output_Dis])
        print("Best Sequence: ", output_seq)
        print("Best Distance: ", output_Dis)

    df = pd.DataFrame(results)
    df.to_csv("incaRNAtion_results.csv")

    # print("----- RNAinverse without start strings --------")
    # for i in range(tries):
    #     start = 'NNNNNNNNNNNNNNNNN'  # without sequence constraints
    #     seq, dis = ViennaRNA.inverse_fold(start, target_structure)
    #     print(seq, dis)



