# 加入训练数据进行fine TUNING

import requests,random,json
from uti import is_rna_sequence,replace_all_t_with_u,replace_all_u_with_t,mask_positions,get_seq,hamming_distance,random_mutate_sequence,count_list_differences


def load_multiple_json_objects(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            try:
                data.append(json.loads(line))
            except json.JSONDecodeError as e:
                print(f"Error decoding JSON from line: {line}")
                raise e
    return data

def test_pre_pos(sample_solutions, secondary_structures,num_masks=1,add_str=True):
    all_error,all_pos_error,all_pos_acc=[],[],[]
    for index, (seq, target) in enumerate(zip(sample_solutions, secondary_structures), start=0):
        if is_rna_sequence(seq):#     随机变化n个碱基
            seq_ori=seq
            position,seq=random_mutate_sequence(seq, num=num_masks)
            seq=replace_all_u_with_t(seq)
            #     预测别学了
            pre_seq,smallest_n_indices= get_seq(seq, target,add_str=add_str,position=position)
            pre_seq=replace_all_t_with_u(pre_seq)
            seq=replace_all_t_with_u(seq)
            #  计算序列预测准确性
            error=100*hamming_distance(seq_ori, pre_seq)/num_masks
            # print(f"task : {index}, seq_error: {error}")
            all_error.append(error)

            #  计算位置预测准确性
            pos_acc=100*count_list_differences(position, smallest_n_indices)/len(position)
            print(f"task : {index}, pos_corr_acc: {pos_acc}, position: {position}, pre_position: {smallest_n_indices}")
            all_pos_acc.append(pos_acc)
        else:
            pass
    total_error=sum(all_error) /len(all_error)
    print(f"total_seq_error: {total_error}")

    total_pos_acc=sum(all_pos_acc) /len(all_pos_acc)
    print(f"total_pos_acc: {total_pos_acc}")
    return total_error,total_pos_acc


# 1. read data
task_name='eterna'
secondary_structures ,sample_solutions= [],[]

if task_name=='eterna':
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
elif task_name=='plant':
    file_path = '../samfeo/data/train.json'
    data =load_multiple_json_objects(file_path)

    plant_structures, plant_solutions = [], []
    for entry in data:
        # Assuming 'structure' and 'solution' keys exist
        secondary_structures.append(entry['label'])
        sample_solutions.append(entry['seq'])

total_error,total_pos_error=test_pre_pos(sample_solutions, secondary_structures,num_masks=2,add_str=True)

1








