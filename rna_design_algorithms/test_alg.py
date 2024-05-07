import gzip
import pandas as pd
import os
from utilis.RNA_fold import get_str
from sa import simulated_annealing
from rnaiverse import call_rnainverse
import argparse
def main(args):
    file_name = args.file_name
    algorithm = args.alg_name
    file_path = os.path.join('..', 'data', file_name)
    with gzip.open(file_path, 'rb') as file:
        data = pd.read_pickle(file)

    row = data.iloc[args.task]
    # Converting lists to structures
    dot_bracket_str = get_str(row)
    print("Target Structure:", dot_bracket_str)
    final_sequence,final_structure,final_distances = None,None,None

    if algorithm == 'rnainverse':
        # Example usage:
        target_structure = dot_bracket_str
        try:
            final_sequence,final_structure,final_distances = call_rnainverse(target_structure)
        except:
            print("RNAinverse failed")
    elif algorithm == 'sa':
        target_structure = dot_bracket_str
        best_seq,best_str, best_distances= simulated_annealing(target_structure)
        final_sequence = best_seq[-1]
        final_structure = best_str[-1]
        final_distances= best_distances[-1]
    print(f"{algorithm} Found RNA Structure: {final_structure}")
    print(f"{algorithm} Found RNA Sequence: {final_sequence}")
    print("Final Distance:", final_distances)


if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument("--alg_name", "-a", dest="alg_name", type=str, default='rnainverse', choices=['sa', 'rnainverse'])
    args.add_argument("--file_name", "-f", dest="file_name", type=str, default='inverse_rna_folding_train.plk.gz')
    args.add_argument("--task", "-d", dest="task", type=int,default=0)  #
    args = args.parse_args()
    main(args)
