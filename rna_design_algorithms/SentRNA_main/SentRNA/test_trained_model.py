import pickle
import os

############ PASTE PARAMS HERE ######################
path = '../models/trained_models/matlot/test/MI-54_trial4'
puzzle_name = 'Mat - Lot 2-2 B'
renderer = 'eterna'
moveset = "[GoodPairsMove([['G', 'C'], ['C', 'G']]), MissingPairsMove([['G', 'C'], ['C', 'G']])]"
n_trajs = 300
MI_tolerance = 1e-5
############ PASTE PARAMS HERE ######################

# try:
#     os.system("python run.py --mode test --renderer %s --input_data ../data/test/eterna100.pkl --test_model %s --results_path out.pkl --test_puzzle_name '%s' --MI_tolerance %f"%(renderer, path, puzzle_name, MI_tolerance))
# except:
#     os.system("python run.py --mode test --renderer \"{renderer}\" --input_data \"../data/test/eterna100.pkl\" --test_model \"{path}\" --results_path \"out.pkl\" --test_puzzle_name \"{puzzle_name}\" --MI_tolerance {MI_tolerance}".format(renderer=renderer, path=path, puzzle_name=puzzle_name, MI_tolerance=MI_tolerance))

for i in range(100):
    # os.system("python run.py --mode refine --input_data test_results/out.pkl --results_path refined_%d.pkl --refine_puzzle_name '%s' --move_set \"%s\" --n_trajs %d"%(i, puzzle_name, moveset, n_trajs))
    command = f'python run.py --mode refine --input_data test_results/model_tested.pkl  --results_path refined_{i}.pkl --refine_puzzle_name "{puzzle_name}" --move_set "{moveset}" --n_trajs {n_trajs}'
    os.system(command)
eff = 0
for i in os.listdir('refined'):
    results = pickle.load(open('refined/%s'%(i)))
    if results[-1][-1] == 1:
        eff += 1
print ('Solution efficacy is %f'%(float(eff) / len(os.listdir('refined'))))
