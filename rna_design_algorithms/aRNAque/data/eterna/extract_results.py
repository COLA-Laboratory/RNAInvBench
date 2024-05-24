import re


def parse_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read()

    data_blocks = re.split(r'={100}\nSolving for the target = ', content)[1:]
    results = []

    for block in data_blocks:
        sequences_info = re.findall(r'([AUGC]+) \| ([()\.]+) \| (-?\d+\.\d+) \| (\d+\.\d+)', block)
        structured_block = [[seq, struct, float(mfe), float(hamming)] for seq, struct, mfe, hamming in sequences_info]
        results.append(structured_block)

    return results


def add_to_seq_ratio(sequence, the_dict):
    print(the_dict)
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



# Usage example:
file_path = 'EternaV2_OP4.out'
parsed_data = parse_file(file_path)
# print(parsed_data)
print("Sequences Tested: ", len(parsed_data))
completed_seq_ratios = {"0-49": 0, "50-99": 0, "100-149": 0, "150-349": 0, "350+": 0}
successes = 0
for x in range(len(parsed_data)):
    if parsed_data[x][0][3] == 1.0 or parsed_data[x][1][3] == 1.0:
            completed_seq_ratios = add_to_seq_ratio(parsed_data[x][0][0], completed_seq_ratios)
            successes += 1

print("SUCCESSES: ", successes)
print("Sequences Tested: ", len(parsed_data))
print("SEQ UNDER 50 SOLVED: ", completed_seq_ratios["0-49"])
print("SEQ UNDER 100 SOLVED: ", completed_seq_ratios["50-99"])
print("SEQ UNDER 150 SOLVED: ", completed_seq_ratios["100-149"])
print("SEQ UNDER 350 SOLVED: ", completed_seq_ratios["150-349"])
print("SEQ ABOVE OR EQUAL TO 350 SOLVED: ", completed_seq_ratios["350+"])