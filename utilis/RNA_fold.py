import RNA
import time
def predict_structure(sequence):
    """
    Predict the secondary structure and minimum free energy of an RNA sequence using ViennaRNA's Python API.
    Parameters:
    - sequence: The RNA sequence, in string format.
    Returns:
    - The structure and free energy.

    # test
    sequence = "GCGCUUCGCCGCGCCG"
    predict_structure(sequence)
    """
    # Create an RNA folding object
    fold_compound = RNA.fold_compound(sequence)
    # Calculate the secondary structure and minimum free energy
    start_time = time.time()
    (structure, mfe) = fold_compound.mfe()
    # End timing
    end_time = time.time()
    # Compute total time
    computation_time = end_time - start_time
    # Print the results
    # print(f"Sequence: {sequence}")
    # print(f"Structure: {structure}")
    # print(f"Minimum Free Energy: {mfe} kcal/mol")
    return structure, mfe, computation_time

def str_distance(s1, s2):
    """Calculate the Hamming distance between two strings"""
    dis= sum(el1 != el2 for el1, el2 in zip(s1, s2))
    return dis


def get_str(row):
    # get structure from data
    sequence_length = row['length']
    pairs = row['pairs']
    # Initialising a braces representation string
    dot_bracket = ['.' for _ in range(sequence_length)]
    # 填充配对信息
    for i, j, _ in pairs:
        dot_bracket[i] = '('
        dot_bracket[j] = ')'
    # Converting lists to strings
    dot_bracket_str = ''.join(dot_bracket)
    return dot_bracket_str