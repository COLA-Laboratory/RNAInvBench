import pickle
import os
import pickle
retrain = False
from sentrna import main, Args

def load_pickle_file(file_path):
    with open(file_path, 'rb') as file:
        data = pickle.load(file)
    return data
# Specify the path to the .pkl file
file_path = r'F:\Python_code\RNA\SentRNA\SentRNA\test_results\model_tested.pkl'
# Load the data from the pickle file
loaded_data = load_pickle_file(file_path)
# Now loaded_data contains the data loaded from the pickle file
for ind in range(len(loaded_data)):
    puzzle_name = loaded_data[ind][0]
    print (puzzle_name)
    renderer = 'eterna'
    moveset = "[GoodPairsMove([['G', 'C'], ['C', 'G']]), MissingPairsMove([['G', 'C'], ['C', 'G']])]"
    MI_tolerance = 1e-5
    model_name = 'refine'
    input_data = 'test_results/model_tested.pkl'
    results_path =f'refined_{ind}.pkl'
    refine_puzzle_name = puzzle_name
    args = Args(
            mode=model_name,
            input_data=input_data,
            results_path=results_path,
            refine_puzzle_name=refine_puzzle_name,
            move_set=moveset,
            MI_tolerance=MI_tolerance,
            renderer=renderer,
            # Set other parameters as required...
        )
    main(args)


eff = 0
for i in os.listdir('refined'):
    with open(os.path.join('refined', i), 'rb') as file:
        results = pickle.load(file)
        print (results[-1][-1])
        if results[-1][-1] >0.999:
            eff += 1
print ('Solution efficacy is %f'%(float(eff) / len(os.listdir('refined'))))
