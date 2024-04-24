# 加入训练数据进行fine TUNING

import requests,random
from uti import is_rna_sequence,replace_all_t_with_u,replace_all_u_with_t,mask_positions,get_seq,hamming_distance




def get_pre_acc(sample_solutions, secondary_structures,num_masks=1,add_str=False):
    all_error=[]
    for index, (seq, target) in enumerate(zip(sample_solutions, secondary_structures), start=0):
        if is_rna_sequence(seq):#     随机变化n个碱基
            seq=replace_all_u_with_t(seq)
            seq_mask = mask_positions(seq, num_masks=num_masks)
            #     预测
            pre_seq= get_seq(seq_mask, target,add_str=add_str)
            pre_seq_E = pre_seq.replace('<mask>', 'E')
            pre_seq_E=replace_all_t_with_u(pre_seq_E)
            seq=replace_all_t_with_u(seq)
            #  计算准确性
            error=100*hamming_distance(seq, pre_seq_E)/num_masks
            print(f"error: {error}")
            all_error.append(error)
        else:
            pass
    total_error=sum(all_error) /len(all_error)
    print(f"total_error: {total_error}")
    return total_error


# 1. read data
file_path = '../samfeo/data/eterna/eterna100_vienna1.txt'
secondary_structures ,sample_solutions= [],[]

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


error=get_pre_acc(sample_solutions, secondary_structures,num_masks=1,add_str=False)










