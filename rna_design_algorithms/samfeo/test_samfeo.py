# 1. read data
path_txt='data/eterna/eterna100.txt'
targets = []
with open(path_txt) as f:
    for line in f:
        targets.append(line.strip())
print(targets)
#2.  input algorithm and conduct it
from samfeo_bk import main
sequence, found_str, target_str, found_dist=main(path=path_txt)
print(found_dist)

#3.  evaluation metrics
alfa=1.e-5
filtered_data = [x for x in found_dist if x < alfa]
count = len(filtered_data)
percentage = (count / len(found_dist)) * 100
print(f"Sucess rate: {percentage:.2f}%")