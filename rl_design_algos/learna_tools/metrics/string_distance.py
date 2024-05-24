import numpy as np


def levenshtein_with_n(designed_sequence=None, designed_structure=None, target_sequence=None, target_structure=None):
    """
    Levenshtein distance between target_structure and designed_structure with 'N' for any char
    """
    # Initialisation for the table of calculations
    rows = len(target_structure) + 1
    cols = len(designed_structure) + 1
    table = [[0 for x in range(cols)] for x in range(rows)]

    for i in range(cols):
        table[0][i] = i

    for i in range(rows):
        table[i][0] = i

    # Calculation of the distances in the table
    for col in range(1, cols):
        for row in range(1, rows):
            # if same char cost for substitution is 0, else 1
            if target_structure[row - 1] == designed_structure[col - 1] or target_structure[row - 1] == 'N' or designed_structure[col - 1] == 'N':
                cost = 0
            else:
                cost = 1
            table[row][col] = min(table[row - 1][col] + 1,        # delete
                                  table[row][col - 1] + 1,        # insert
                                  table[row - 1][col - 1] + cost) # substitution

    #bottom right in the table is the levenshtein distance
    return table[row][col] / len(designed_structure)


def number_of_brackets_with_n(designed_sequence=None, designed_structure=None, target_sequence=None, target_structure=None):
    """
    difference of number of brackets between designed_structure and target_structure
    """
    mismatches =  abs(designed_structure.count('(') - target_structure.count('(')) + \
           abs(designed_structure.count(')') - target_structure.count(')')) + \
           abs(designed_structure.count('.') - target_structure.count('.'))

    return (mismatches / 2) - target_structure.count('N')


def hamming_with_n(designed_sequence=None, designed_structure=None, target_sequence=None, target_structure=None):
    l1 = np.asarray(list(designed_structure))
    l2 = np.asarray(list(target_structure))
    distance = np.sum(((l1 != l2) & (l1 != 'N') & (l2 != 'N')).astype(np.int8))
    return distance   #  / len(designed_structure)
