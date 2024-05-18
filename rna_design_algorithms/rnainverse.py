import ViennaRNA
import time
import pandas as pd

def call_rnainverse(target_structure, constraint=None, tries=1, save_file='results/rnainverse_result'):
    """
    Call RNAinverse from ViennaRNA package to find an RNA sequence that folds into the target structure.
    Parameters:
    - target_structure (str): The target RNA secondary structure in dot-bracket notation.
    - constraints (str): sequence constraints, if we have else None
    - tries (int): try multiple times to get multiple designed sequences for one target structure.
    Returns:
    - str: An RNA sequence that is predicted to fold into the target structure.
    """
    if constraint is None:
        constraint = ''.join(['N'] * len(target_structure))
    sequences_design = []
    distance_design = []
    time_consume = []
    for i in range(tries):
        start_time = time.time()
        seq, dis = ViennaRNA.inverse_fold(constraint, target_structure)
        end_time = time.time()
        sequences_design.append(seq)
        distance_design.append(dis)
        time_consume.append(end_time-start_time)

    # save to file
    data = {'target_structure': target_structure, 'seq_constraints': constraint, 'sequence': sequences_design, 'time': time_consume, 'distance': distance_design}
    df = pd.DataFrame(data)
    save_file = save_file + '.pkl'
    df.to_pickle(save_file)

    return sequences_design

if __name__ == '__main__':
    # Example usage:
    target_structure = "(((....)))((...))"
    rna_sequence = call_rnainverse(target_structure, tries=3)
    print(f"RNA sequence that folds into the target structure: {rna_sequence}")
    df = pd.read_pickle("results/rnainverse_result.pkl")
    pass
