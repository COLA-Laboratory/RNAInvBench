import RNA
import random
import math
import time
from utilis.RNA_fold import str_distance, predict_structure
import matplotlib.pyplot as plt


def mutate_sequence(seq):
    """随机变异序列中的一个碱基"""
    nucleotides = ['A', 'C', 'G', 'U']
    pos = random.randint(0, len(seq) - 1)
    original_nucleotide = seq[pos]
    nucleotides.remove(original_nucleotide)
    mutated_nucleotide = random.choice(nucleotides)
    return seq[:pos] + mutated_nucleotide + seq[pos + 1:]


def simulated_annealing(target_structure, initial_sequence, temperature=100.0, cooling_rate=0.95, max_steps=1000):
    current_seq = initial_sequence
    current_structure, _ , _= predict_structure(current_seq)
    current_distance = str_distance(current_structure, target_structure)
    best_distances = [current_distance]
    best_seq = [current_seq]

    for step in range(max_steps):
        temperature *= cooling_rate
        candidate_seq = mutate_sequence(current_seq)
        candidate_structure, _, _ = predict_structure(candidate_seq)
        candidate_distance = str_distance(candidate_structure, target_structure)

        if candidate_distance < current_distance or random.random() < math.exp(
                (current_distance - candidate_distance) / temperature):
            current_seq = candidate_seq
            current_distance = candidate_distance
            current_structure = candidate_structure
        if current_distance < best_distances[-1]:
            best_distances.append(current_distance)
            best_seq.append(current_seq)
        else:
            best_distances.append(best_distances[-1])
            best_seq.append(best_seq[-1])

        if current_distance == 0:
            break

    return current_seq, current_structure, current_distance, best_distances


# define main to get example working
if __name__ == '__main__':
    target_structure = "(((....)))((...))"
    initial_sequence = "GCGCUUCGCCGCGCCG"

    final_sequence, final_structure, final_distance, best_distances = simulated_annealing(target_structure,
                                                                                          initial_sequence)
    print("Final Sequence:", final_sequence)
    print("Final Structure:", final_structure)
    print("Final Distance:", final_distance)

    # 绘制最优解的迭代曲线
    plt.plot(best_distances)
    plt.xlabel('Iteration')
    plt.ylabel('Best Hamming Distance')
    plt.title('Best Solution Iteration Curve')
    plt.show()
