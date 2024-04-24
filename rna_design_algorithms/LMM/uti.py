import RNA
import random
import math,requests
import time
from utilis.RNA_fold import str_distance, predict_structure
import matplotlib.pyplot as plt

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

def is_rna_sequence(seq):
    # Set of RNA nucleotides
    rna_nucleotides = {'A', 'C', 'G', 'U'}
    # Check if every character in the sequence is a valid RNA nucleotide
    return all(char in rna_nucleotides for char in seq)


def get_seq(seq, target,add_str=True,position=[]):
    # Flask服务的URL，根据您的Flask应用的实际运行情况进行修改
    url = 'http://144.173.65.124:5000/predict'
    if add_str:
        tes_seq = f"{seq}<eos>{target}"
    else:
        tes_seq =f"{seq}"
    # 构造POST请求的数据，这里的键需要与服务端接收的键匹配
    data = {'seq': tes_seq}
    # 发送POST请求
    response = requests.post(url, json=data)
    # 检查响应状态码
    if response.status_code == 200:
        # 请求成功，打印预测结果
        # print('Prediction result:', response.json())
        res = response.json()
        try:
            pre_seq = res['sequence'][:len(target)]
            pro = res['probs']
            smallest_n_indices=[]
            # n =len(position)
            # Calculate the minimum probability for each entry
            min_values = [(index, max(prob.values())) for index, prob in enumerate(pro[:len(target)])]
            threshould= 0.1
            below_threshould = [(index, value) for index, value in min_values if value < threshould]
            smallest_n_indices = [index for index, value in below_threshould]

            # Sort by the minimum probability values
            min_values_sorted = sorted(min_values, key=lambda x: x[1])
            # Get the indices of the entries with the smallest 'n' minimum probabilities
            smallest_n_indices = [index for index, _ in min_values_sorted[:len(target)]]
        except:
            1
    else:
        # 请求失败，打印错误信息
        print('Error:', response.json())
        pre_seq = seq
    if len(position)>0:
        return pre_seq,smallest_n_indices
    else:
        return pre_seq



def count_list_differences(list1, list2):
    set1 = set(list1)
    # Count elements in list2 that are also in set1
    return sum(1 for element in list2 if element in set1)




def mask_positions(sequence, num_masks):
    # Convert the string to a list of characters to manipulate it
    sequence_list = list(sequence)
    # Generate 'num_masks' unique random positions to mask
    if num_masks < len(sequence_list):  # Ensure we don't exceed the length of the sequence
        positions = random.sample(range(len(sequence_list)), num_masks)
        # Replace the selected positions with '<mask>'
        for pos in positions:
            sequence_list[pos] = '<mask>'
    return ''.join(sequence_list)  # Convert list back to string

def random_mutate_sequence(seq, num):
    """Randomly mutate specified number of positions in the RNA sequence to a different nucleotide."""
    # Ensure we're working with a mutable version of the sequence
    seq_list = list(seq)
    position = random.sample(range(0,len(seq)), num)
    for i in position:
        # Original nucleotide at the position
        original_nucleotide = seq_list[i]
        # Possible nucleotides excluding the original one
        nucleotides = ['A', 'C', 'G', 'U']
        nucleotides.remove(original_nucleotide)
        # Mutate to a different random nucleotide
        seq_list[i] = random.choice(nucleotides)
    # Join the list back into a string to form the new sequence
    new_seq = ''.join(seq_list)
    return position,new_seq




def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length")
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def mask_sequence(sequence, positions):
    """
    Masks specified positions in a given sequence with '<mask>'.

    Args:
    sequence (str): The original sequence.
    positions (list of int): List of positions to mask.

    Returns:
    str: The sequence with specified positions replaced by '<mask>'.
    """
    # Convert the sequence to a list of characters for mutability
    sequence_list = list(sequence)

    # Loop through each specified position and replace with '<mask>'
    for pos in positions:
        if pos < len(sequence_list):  # Check if the position is within the sequence length
            sequence_list[pos] = '<mask>'
        else:
            raise IndexError("Position out of the sequence range")

    # Join the list back into a string
    masked_sequence = ''.join(sequence_list)
    return masked_sequence






# from transformers import BertTokenizer, TFBertModel
# import tensorflow as tf
#
# # Load pre-trained model tokenizer (vocabulary) and model
# tokenizer = BertTokenizer.from_pretrained('bert-base-uncased')
# model = TFBertModel.from_pretrained('bert-base-uncased')
#
# def get_embeddings(text):
#     """
#     Generate embeddings for the input text using BERT.
#     """
#     # Encode text
#     encoded_input = tokenizer(text, return_tensors='tf', truncation=True, max_length=512)
#     # Compute embeddings
#     output = model(encoded_input)
#     # Use the mean of the last hidden state from BERT as the sentence embedding
#     embeddings = tf.reduce_mean(output.last_hidden_state, axis=1)
#     return embeddings
#
# def cosine_similarity(embeddings1, embeddings2):
#     """
#     Compute the cosine similarity between two embeddings.
#     """
#     # Normalize embeddings
#     normalized_embeddings1 = tf.nn.l2_normalize(embeddings1, axis=1)
#     normalized_embeddings2 = tf.nn.l2_normalize(embeddings2, axis=1)
#     # Compute similarity
#     cosine_similarities = tf.reduce_sum(tf.multiply(normalized_embeddings1, normalized_embeddings2), axis=1)
#     return cosine_similarities.numpy()
#
# # Example texts
# text1 = "This is a sample sentence."
# text2 = "This is a similar sentence."
#
# # Get embeddings
# embeddings1 = get_embeddings(text1)
# embeddings2 = get_embeddings(text2)
#
# # Calculate cosine similarity
# similarity = cosine_similarity(embeddings1, embeddings2)
# print("Cosine Similarity:", similarity)
#
#
