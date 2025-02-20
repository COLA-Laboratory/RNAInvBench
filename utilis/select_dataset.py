import pandas as pd

datasets = ['rfam_taneda-test', 'rfam_learn-test', 'anta_strand-test', 'ArchiveII-test', 'eterna100_v2', 'eterna100_v1']

def select_dataset(file, dataset=None, range_len=None, rep=0):
    """
    select data from one dataset

    parameters:
    - file (str): from which pkl file
    - database (str): name of dataset
        - one of ['rfam_taneda-test', 'rfam_learn-test', 'anta_strand-test', 'ArchiveII-test', 'eterna100_v2', 'eterna100_v1']
    - range_len ([min, max]): range of target length
        - default: None, no length limits.
    - rep (int, [0, 1]): representation of structure
        - 0 for dot-bracket (default)
        - 1 for pairs

    return:
    - target_structure
    - sequence constraints
    """
    df = pd.read_pickle(file)

    # database
    if dataset is not None:
        if isinstance(dataset, str):
            dataset = [dataset]

        df = df[df['origin'].isin(dataset)]

    # length range
    if range_len is not None:
        min_len, max_len = range_len
        df = df[(df['length'] >= min_len) & (df['length'] <= max_len)]

    # representation
    if rep == 0:
        target_structure = df['dotbracket'].tolist()
    elif rep == 1:
        target_structure = df['pairs'].tolist()
    else:
        raise ValueError("Invalid value for 'rep'. Must be 0 or 1.")

    # sequence
    sequence_constraints = df['sequence'].tolist()

    return target_structure, sequence_constraints

if __name__ == "__main__":
    # example

    dir = '../data/'
    file = 'inverse_rna_folding_train_dotbracket.pkl.gz'
    path = dir + file

    # dataset = 'eterna100_v1'  # one database
    dataset = ['rfam_taneda-test', 'rfam_learn-test']  # multiple databases
    range_len = [1, 500]
    target_structure, _ = select_dataset(path, dataset, range_len)
    print(target_structure)