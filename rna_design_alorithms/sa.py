import RNA
import random
import math
import time
from utilis.RNA_fold import str_distance, predict_structure
import matplotlib.pyplot as plt

def generate_initial_sequence(structure_length):
    """
    随机生成与目标结构长度相等的RNA序列。
    """
    nucleotides = ['A', 'U', 'C', 'G']
    return ''.join(random.choice(nucleotides) for _ in range(structure_length))

def mutate_sequence(seq):
    """随机变异序列中的一个碱基"""
    nucleotides = ['A', 'C', 'G', 'U']
    pos = random.randint(0, len(seq) - 1)
    original_nucleotide = seq[pos]
    nucleotides.remove(original_nucleotide)
    mutated_nucleotide = random.choice(nucleotides)
    return seq[:pos] + mutated_nucleotide + seq[pos + 1:]


def simulated_annealing(target_structure,  temperature=100.0, cooling_rate=0.95, max_steps=1000):
    initial_sequence = generate_initial_sequence(len(target_structure))
    current_seq = initial_sequence
    current_structure, _ , _= predict_structure(current_seq)
    current_distance = str_distance(current_structure, target_structure)
    best_distances = [current_distance]
    best_seq = [current_seq]
    best_str = [current_structure]

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
            best_str.append(current_structure)
        else:
            best_distances.append(best_distances[-1])
            best_seq.append(best_seq[-1])
            best_str.append(best_str[-1])

        if current_distance == 0:
            break

    return best_seq,best_str, best_distances


# define main to get example working
if __name__ == '__main__':
    target_structure = "(((....)))((...))"

    best_seq,best_str, best_distances= simulated_annealing(target_structure)
    print("Final Sequence:", best_seq[-1])
    print("Final Structure:", best_str[-1])
    print("Final Distance:", best_distances[-1])

    # 绘制最优解的迭代曲线
    plt.plot(best_distances)
    plt.xlabel('Iteration')
    plt.ylabel('Best Hamming Distance')
    plt.title('Best Solution Iteration Curve')
    plt.show()
