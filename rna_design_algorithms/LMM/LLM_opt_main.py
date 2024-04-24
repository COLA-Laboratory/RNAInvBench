import RNA
import random
import math
import numpy as np
import time
from utilis.RNA_fold import str_distance, predict_structure
import matplotlib.pyplot as plt
from uti import replace_all_u_with_t,replace_all_t_with_u,mask_positions,get_seq,hamming_distance,random_mutate_sequence,count_list_differences,mask_sequence
def generate_initial_sequence(structure_length):
    """
    随机生成与目标结构长度相等的RNA序列。
    """
    nucleotides = ['A', 'U', 'C', 'G']
    return ''.join(random.choice(nucleotides) for _ in range(structure_length))


def compare_structures(target_structure, current_structure):
    # 获取两个结构的长度，它们应该是相同的
    if len(target_structure) != len(current_structure):
        raise ValueError("Structures must be of the same length")
    # 初始化差异位置列表
    diff_positions = []
    # 遍历每个结构的每个位置
    for i in range(len(target_structure)):
        if target_structure[i] != current_structure[i]:
            # 记录不同的位置
            diff_positions.append(i + 1)  # 加1为了使用1-based index
    return diff_positions

def replace_all_u_with_t(sequence):
    """
    Replace all 'U's in a given RNA sequence with 'T's.

    Args:
    sequence (str): RNA sequence to be modified.

    Returns:
    str: Modified RNA sequence with all 'U's replaced by 'T's.
    """
    return sequence.replace('U', 'T')

def replace_all_t_with_u(sequence):
    """
    Replace all 'T's in a given RNA sequence with 'U's.

    Args:
    sequence (str): RNA sequence to be modified.

    Returns:
    str: Modified RNA sequence with all 'T's replaced by 'U's.
    """
    return sequence.replace('T', 'U')


def integrate_sequences(current_seq, positions, target):
    # Replace specific positions in current_seq with <mask>
    modified_seq_parts = [
        "<mask>" if i in positions else char
        for i, char in enumerate(current_seq)
    ]

    # Convert the list back to a string
    modified_seq = ''.join(modified_seq_parts)

    # Concatenate the modified sequence with the target sequence, adding <eos> in between
    final_seq = f"{modified_seq}<eos>{target}"
    return final_seq


def llm_predict_seq(current_seq,position,target):
    import requests
    # Flask服务的URL，根据您的Flask应用的实际运行情况进行修改
    url = 'http://144.173.65.124:5000/predict'
    # 把U改成T
    modified_sequence_u_t = replace_all_u_with_t(current_seq)

    tes_seq = f"{modified_sequence_u_t}<eos>{target}"
    data = {'seq': tes_seq}

    # 发送POST请求
    response = requests.post(url, json=data)
    if response.status_code == 200:
        # 请求成功，打印预测结果
        print('Prediction result:', response.json())
        new_seq_o = response.json()['predicted_sequence'].replace(" ", "")
        # 把T改成U
        new_seq = replace_all_t_with_u(new_seq_o[:len(current_seq)])

    else:
        # 请求失败，打印错误信息
        print('Error:', response.json())


    seq = integrate_sequences(modified_sequence_u_t, position, target)

    # 构造POST请求的数据，这里的键需要与服务端接收的键匹配
    data = {'seq': seq}
    # 发送POST请求
    # 检查响应状态码
    try:
        response = requests.post(url, json=data)
        if response.status_code == 200:
            # 请求成功，打印预测结果
            print('Prediction result:', response.json())
            new_seq_o=response.json()['predicted_sequence'].replace(" ", "")
            # 把T改成U
            new_seq  = replace_all_t_with_u(new_seq_o[:len(current_seq)])

        else:
            # 请求失败，打印错误信息
            print('Error:', response.json())
    except:
        new_seq = current_seq
    return new_seq




def mutate_sequence(seq,current_structure,target_structure):
    """随机变异序列中的一个碱基"""
    differences = compare_structures(target_structure, current_structure)
    num_elements = min(2, len(differences))  # Ensure at least 2 elements are selected
    position= random.sample(differences, num_elements)
    # position= random.sample(range(0, len(current_structure)), 2)
    #  找到低置信度的位置



    new_seq=llm_predict_seq(seq, position, target_structure)

    return new_seq


from uti import random_mutate_sequence, replace_all_u_with_t, get_seq, replace_all_t_with_u


def llm_op(target_structure,  max_steps=100,add_str=True):
    # initial_sequence = generate_initial_sequence(len(target_structure))
    seq= '<mask>'*len(target_structure)
    initial_sequence, _ = get_seq(seq, target_structure, add_str=add_str,position=[0])
    initial_sequence = initial_sequence.replace('.', 'G')
    current_seq = initial_sequence
    current_structure, _ , _= predict_structure(current_seq)
    current_distance = str_distance(current_structure, target_structure)
    best_distances = [current_distance]
    best_seq = [current_seq]
    best_str = [current_structure]

    visited_seqs = set()

    for step in range(max_steps):
        seq = current_seq
        pre_seq, pre_position = get_seq(seq, target_structure, add_str=add_str,position=[0])
        n_candi,candi_seq,candi_str=100,[],[]
        pre_position = pre_position[:10]
        for idx in pre_position:
            seq_mask = mask_sequence(seq, positions=[idx])
            pre_seq= get_seq(seq_mask, target_structure,add_str=add_str)
            pre_seq = [s1 if i!=idx else s2 for i, (s1, s2) in enumerate(zip(seq, pre_seq))]
            candidate_seq=''.join(pre_seq)
            # print(f"Step {step}: Distance {best_distances[-1]} seq: {seq} pre_seq: {pre_seq}")
            visited_seqs.add(candidate_seq)
            # print(len(visited_seqs))
            candi_seq.append(candidate_seq)
            candi_structure, _, _ = predict_structure(candidate_seq)
            candi_str.append(candi_structure)
        # 选择适应值最小的位置作为最后的位置

            candidate_structure, _, _ = predict_structure(candidate_seq)
            candidate_distance = str_distance(candidate_structure, target_structure)

            if candidate_distance < current_distance:
                current_seq = candidate_seq
                current_distance = candidate_distance
                current_structure = candidate_structure
            if current_distance < min(best_distances):
                best_distances.append(current_distance)
                best_seq.append(current_seq)
                best_str.append(current_structure)
            else:
                best_distances.append(best_distances[-1])
                best_seq.append(best_seq[-1])
                best_str.append(best_str[-1])

            if current_distance == 0:
                break
        # print(f"Step {step}: Distance {best_distances[-1]}")
    return best_seq,best_str, best_distances


# define main to get example working
if __name__ == '__main__':
    secondary_structures, sample_solutions = [], []

    file_path = '../samfeo/data/eterna/eterna100_vienna1.txt'

    # Open and read the file
    with open(file_path, 'r') as file:
        # Read the lines from the file
        lines = file.readlines()
        # Skip the header line
        for line in lines[1:]:
            columns = line.split('\t')  # Assuming tab-separated values (TSV format)
            if len(columns) >= 7:  # Check if there are enough columns
                secondary_structures.append(columns[4])  # Secondary Structure column
                sample_solutions.append(columns[6])  # Sample Solution (1) column

    for index, (seq, target) in enumerate(zip(sample_solutions, secondary_structures), start=0):

        target_structure =target#    GUUUGGAAAAACAAAC


        best_seq,best_str, best_distances= llm_op(target_structure)
        print("Final Sequence:", best_seq[-1])
        print("Final Structure:", best_str[-1])
        print("Final Distance:", best_distances[-1])
        print(f"Task {index}: Distance {best_distances[-1]}")

        # 绘制最优解的迭代曲线
        plt.plot(best_distances)
        plt.xlabel('Iteration')
        plt.ylabel('Best Hamming Distance')
        plt.title('Best Solution Iteration Curve')
        plt.show()
