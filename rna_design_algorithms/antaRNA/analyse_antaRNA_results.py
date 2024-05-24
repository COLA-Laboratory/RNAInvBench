import numpy as np
import pandas as pd
import ast
import RNA
import subprocess

def add_to_seq_ratio(sequence, the_dict):
    if len(sequence) >= 350:
        the_dict["350+"] += 1
    elif len(sequence) >= 150:
        the_dict["150-349"] += 1
    elif len(sequence) >= 100:
        the_dict["100-149"] += 1
    elif len(sequence) >= 50:
        the_dict["50-99"] += 1
    else:
        the_dict["0-49"] += 1
    return the_dict


def hamming_distance(seq1, seq2):
    """Calculate the Hamming distance between two sequences of equal length."""
    print(seq1, " ", seq2)
    if len(seq1) != len(seq2):
        # raise ValueError("Sequences must be of equal length")
        print("antaRNA failed to produce a valid sequence")
        return 100
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def fold_seq(sequence):
    # Call RNAfold with the sequence as input
    process = subprocess.Popen(['RNAfold'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate(input=sequence + b'\n')

    # Check if there was an error
    if process.returncode != 0:
        print("Error folding RNA sequence:", stderr.decode('utf-8'))
        return None, None

    # Process the output
    output = stdout.decode('utf-8').split('\n')
    if len(output) < 2:
        print("Unexpected output format:", stdout.decode('utf-8'))
        return None, None

    folded_structure = output[1].split(' ')[0]  # The folded structure
    mfe = output[1].split(' ')[-1].strip('()')  # The MFE value

    return folded_structure, float(mfe)


df = pd.read_csv("antaRNA_Results.csv", skiprows=1, header=None)
df2 = pd.read_csv("Eterna100.txt")
e_structs = df2["str"].to_list()
non_nan_rows = df
non_nan_data = []

for index, row in non_nan_rows.iterrows():
    non_nan_values = row.dropna()
    non_nan_data.append((index, non_nan_values))

data = []

# Display the extracted non-NaN data
print("\nNon-NaN data:")
for index, values in non_nan_data:
    print("Row ", index, ": ")
    data.append(values)
seqs = []
for x in range(len(data)):
    seq = ""
    cleaned_list = [y for y in data[x] if y]
    if isinstance(cleaned_list[0], str):
        seq = ast.literal_eval(cleaned_list[0])
    elif isinstance(cleaned_list[1], str):
        seq = ast.literal_eval(cleaned_list[1])
    seqs.append(seq[0][0])

successes = 0
completed_seq_ratios = {"0-49": 0, "50-99": 0, "100-149": 0, "150-349": 0, "350+": 0}
count = 0
for seq in seqs:
    struct, mfe = fold_seq(seq)
    distance = hamming_distance(struct, e_structs[count])
    count += 1
    if distance == 0:
        completed_seq_ratios = add_to_seq_ratio(seq, completed_seq_ratios)
        successes += 1
    print(count, " : ", distance)
print("SUCCESSES: ", successes)
print("Sequences Tested: ", len(seqs))
print("SEQ UNDER 50 SOLVED: ", completed_seq_ratios["0-49"])
print("SEQ UNDER 100 SOLVED: ", completed_seq_ratios["50-99"])
print("SEQ UNDER 150 SOLVED: ", completed_seq_ratios["100-149"])
print("SEQ UNDER 350 SOLVED: ", completed_seq_ratios["150-349"])
print("SEQ ABOVE OR EQUAL TO 350 SOLVED: ", completed_seq_ratios["350+"])
