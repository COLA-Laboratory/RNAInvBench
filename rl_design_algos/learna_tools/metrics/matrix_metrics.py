import numpy as np
from collections import defaultdict

def number_of_mismatches(m1, m2):
    """
    counts the number of differences of the two matrices m1 and m2
    """
    missmatches = 0
    for x in range(len(m1)):
        for y in range(len(m1[x])):
            if m1[x][y] == m2[x][y]:
                missmatches += 1

    return missmatches


def db2mat(db):
    """
    Get matrix from dot-bracket notation of structure.
    Does not support partial RNA design with missing structure regions (Ns).

    Input: dot-bracket notation string or list
    Returns: Binary matrix with 1 indicating pair, 0 indicating no pair
    """
    corresponding_brackets = {')': '(', ']': '[', '}': '{', '>': '<'}
    mat = np.zeros((len(db), len(db)))
    pair_stack = defaultdict(list)

    for i, sym in enumerate(db):
        if sym == '.':
            continue
        elif sym in ['(', '[', '{', '<']:
            pair_stack[sym].append(i)
        else:
            open_idx = pair_stack[corresponding_brackets[sym]].pop()
            mat[open_idx, i] = 1
            mat[i, open_idx] = 1
    return mat


def score_matrix(true_mat, pred_mat):
    """
    Score two matrices with common metrics: F1-score, recall, specificity,
    precision, Matthews-Correlation-Coeficient (MCC).
    Does not support partial RNA design with missing structure regions.

    Input: true matrix, predicted matrix.
    Returns: Dict with scores.
    """
    solved = np.all(np.equal(true_mat, pred_mat)).astype(np.int8)
    if solved == 1:
        f1_score = 1
        precision = 1
        recall = 1
        specificity = 1
        mcc = 1
    else:
        # tp = true-positives
        # fp = false-positives
        # tn = true-negatives
        # fn = false-negatives
        tp = np.logical_and(pred_mat, true_mat).sum()
        non_correct = (tp == 0).astype(int)
        tn = np.logical_and(np.logical_not(pred_mat), np.logical_not(true_mat)).sum()
        fp = pred_mat.sum() - tp
        fn = true_mat.sum() - tp

        recall = tp / (tp + fn + 1e-8)
        precision = tp / (tp + fp + 1e-8)
        specificity = tn / (tn + fp + 1e-8)
        f1_score = 2 * tp / (2 * tp + fp + fn)
        mcc = (tp * tn - fp * fn) / np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) + 1e-8)

    mat_metrics = {'f1_score': f1_score, 'solved': solved}
    mat_metrics['precision'] = precision
    mat_metrics['recall'] = recall
    mat_metrics['specificity'] = specificity
    mat_metrics['mcc'] = mcc

    return mat_metrics



if __name__ == '__main__':
    from matplotlib import pyplot as plt

    true_strucs = ['....((((....))))....', '...(([[[))..]]]...', '..(())...']
    pred_strucs = ['...(((((....)))))...', '....(((((...))))).', '..(())...']
    solved = [False, False, True]

    for i, (t, p) in enumerate(zip(true_strucs, pred_strucs)):
        true = db2mat(t)
        pred = db2mat(p)
        scores = score_matrix(true, pred)
        assert scores['solved'] == solved[i]

        fix, axs = plt.subplots(1, 2)
        axs[0].imshow(true)
        axs[1].imshow(pred)
        plt.show()
