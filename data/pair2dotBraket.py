import pandas as pd
import os
from tqdm import tqdm

directory = "../data/"

def to_dotBracket():
    """
    convert paired representation to dot-bracket representation
    add dot-bracket rep. to the last column of original file
    """
    for file_name in os.listdir(directory):
        if file_name.endswith(".plk.gz"):
            print("processing " + file_name + " ....")
            dotbrackets = []
            df = pd.read_pickle(directory+file_name)
            # every line
            for _, row in tqdm(df.iterrows(), total=len(df), desc='Processing'):
                pairs = row['pairs']
                seq_length = row['length']
                # init dot-bracket with all .
                dotbracket = ['.' for _ in range(seq_length)]
                # traverse all pairs
                for i, j, typ in pairs:
                    # only consider the standard pairs[i, j, 0]
                    if typ == 0:
                        dotbracket[i] = '('
                        dotbracket[j] = ')'
                dotbrackets.append(''.join(dotbracket))

            # add a new column 'dotbracket'
            df['dotbracket'] = dotbrackets
            # save as a new pklfile
            base_name, ext = os.path.splitext(file_name)
            if ext == '.gz':
                base_name, _ = os.path.splitext(base_name)
            new_pkl_file = base_name + '_dotbracket.pkl.gz'
            df.to_pickle(new_pkl_file)

if __name__ == "__main__":
    # # ---------- already executed -----------
    # to_dotBracket()

    # print one file
    file_name = 'inverse_rna_folding_benchmark_dotbracket.pkl.gz'
    df = pd.read_pickle(file_name)
    print(df['pairs'], df['dotbracket'])