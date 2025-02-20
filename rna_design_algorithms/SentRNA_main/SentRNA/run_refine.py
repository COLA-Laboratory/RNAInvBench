import numpy as np
import random
import os
import pickle
from util.compute_mfe import *
from util.draw_rna import *
from util.mutinf import *
from util.refine_moves import *
import argparse
from util.featurize_util import *
from util.feedforward import *


def train(results_path, n_layers, hidden_size, nb_epochs, mini_epoch, MI_features_list, nn_training_dataset, val_dataset, test_dataset, MI_tolerance, renderer):
    # Model params
    input_size = len(nn_training_dataset[0][0])
    layer_sizes = [input_size] + [hidden_size] * n_layers + [4]

    import os

    results_path = "model_tested"  # 假设这是你的动态路径部分
    directory = 'results'
    file_path = f'results/MI_features_list.{results_path}.pkl'

    # 如果目录不存在，则创建它
    os.makedirs(directory, exist_ok=True)

    # 现在，你可以安全地保存文件
    with open(file_path, 'wb') as file:  # 注意：使用二进制写模式'wb'
        pickle.dump(MI_features_list, file)

    # pickle.dump(MI_features_list, open('results/MI_features_list.%s.pkl'%(results_path), 'w'))
    # pickle.dump(layer_sizes, open('results/layer_sizes.%s.pkl'%(results_path), 'w'))
    val_accuracies = []
    val_solutions = []
    model = TensorflowClassifierModel(layer_sizes=layer_sizes)
    for i in range(int(nb_epochs / mini_epoch)):
        print ('Mini epoch %d'%(i))
        save_path = '%s_mini-epoch_%d.ckpt'%(results_path, i)
        restore_path = '%s_mini-epoch_%d.ckpt'%(results_path, i - 1)
        val_model = 'test/%s-%d'%(save_path, mini_epoch + 1)
        if '%s-%d.meta'%(restore_path, mini_epoch + 1) in os.listdir('test'):
            checkpoint = 'test/%s-%d'%(restore_path, mini_epoch + 1)
        else:
            checkpoint = None
        model.fit(nn_training_dataset, loss_thresh=1e-9, nb_epochs=mini_epoch, save_path=save_path, checkpoint=checkpoint)
        # Validation
        dot_bracket, seq, fixed_bases = val_dataset
        val_solution, val_input = model.evaluate(dot_bracket, seq, fixed_bases, layer_sizes, MI_features_list, val_model, refine=False, MI_tolerance=MI_tolerance, \
                                                 renderer=renderer)
        val_accuracy = check_answer(val_solution, dot_bracket)
        val_solution2, val_input2 = model.evaluate(dot_bracket, val_solution, fixed_bases, layer_sizes, MI_features_list, val_model, refine=True, \
                                                   MI_tolerance=MI_tolerance, renderer=renderer)
        val_accuracy2 = check_answer(val_solution2, dot_bracket)
        if val_accuracy2 > val_accuracy:
            val_accuracy = val_accuracy2
            val_solution = val_solution2
        val_solutions.append(val_solution)
        val_accuracies.append(val_accuracy)
    val_accuracies = np.array(val_accuracies)
    best_model = np.where(np.array(val_accuracies) == np.max(val_accuracies))[0][0]
    # Testing
    test_model = 'test/%s_mini-epoch_%d.ckpt-%d'%(results_path, best_model, mini_epoch + 1)
    dot_bracket, seq, fixed_bases = test_dataset
    test_solution, test_input = model.evaluate(dot_bracket, seq, fixed_bases, layer_sizes, MI_features_list, test_model, refine=False, MI_tolerance=MI_tolerance, \
                                               renderer=renderer)
    test_accuracy = check_answer(test_solution, dot_bracket)
    test_solution2, test_input2 = model.evaluate(dot_bracket, test_solution, fixed_bases, layer_sizes, MI_features_list, test_model, refine=True, \
                                                 MI_tolerance=MI_tolerance, renderer=renderer)
    test_accuracy2 = check_answer(test_solution2, dot_bracket)
    if test_accuracy2 > test_accuracy:
        test_solution = test_solution2
        test_accuracy = test_accuracy2
    pickle.dump([best_model, val_solutions[best_model], val_accuracies[best_model], test_solution, test_accuracy], open('results/%s.pkl'%(results_path), 'wb'))
    try:
        os.system('mv test/%s_mini-epoch_%d* test/%s'%(results_path, best_model, results_path))
        os.system('rm test/*')
    except:
        import shutil
        import glob
        import os

        # 构造源文件的路径模式
        source_pattern = f'test/{results_path}_mini-epoch_{best_model}*'
        # 目标目录路径
        destination_dir = f'test/{results_path}'

        # 确保目标目录存在
        if not os.path.exists(destination_dir):
            os.makedirs(destination_dir)

        # 移动匹配的文件
        for file_path in glob.glob(source_pattern):
            # 计算目标文件路径
            destination_file_path = os.path.join(destination_dir, os.path.basename(file_path))
            # 移动文件
            shutil.move(file_path, destination_file_path)
        test_dir = 'test/'
        for item in os.listdir(test_dir):
            item_path = os.path.join(test_dir, item)
            if os.path.isfile(item_path):
                os.remove(item_path)  # 删除文件
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)  # 删除目录及其所有内容

    return 0


def test(dataset, model, results_path, puzzle_name, MI_tolerance, renderer):
    test_puzzles = [i[0] for i in dataset]
    model_name = model[model.index('test') + 5:]
    base_dir = model[:model.index('test')][:-1]
    if base_dir == '':
        base_dir = '.'
    for i in os.listdir(model):
        if '.data' in i:
            model_path = i[:i.index('.data')]
            test_model_path = '%s/%s'%(model, model_path)
    try:
        MI_features_list = pickle.load(open('%s/results/MI_features_list.%s.pkl'%(base_dir, model_name)))
        layer_sizes = pickle.load(open('%s/results/layer_sizes.%s.pkl' % (base_dir, model_name)))

    except:
        MI_features_list = pickle.load(open('%s/results/MI_features_list.%s.pkl' % (base_dir, model_name), 'rb'))
        layer_sizes = pickle.load(open('%s/results/layer_sizes.%s.pkl'%(base_dir, model_name), 'rb'))

    model = TensorflowClassifierModel(layer_sizes=layer_sizes)
    output = []
    debug = []
    for puzzle in test_puzzles:
        if puzzle_name != '-1' and puzzle != puzzle_name:
            continue
        dot_bracket, seq, fixed_bases = parse_progression_dataset(dataset, [puzzle], 1, MI_features_list, evaluate=True)
        solution, _ = model.evaluate(dot_bracket, seq, fixed_bases, layer_sizes, MI_features_list, test_model_path, refine=False, MI_tolerance=MI_tolerance, \
                                     renderer=renderer)
        accuracy = check_answer(solution, dot_bracket)
        solution2, _ = model.evaluate(dot_bracket, solution, fixed_bases, layer_sizes, MI_features_list, test_model_path, refine=True, MI_tolerance=MI_tolerance, \
                                      renderer=renderer)
        accuracy2 = check_answer(solution2, dot_bracket)
        if accuracy2 > accuracy:
            accuracy = accuracy2
            solution = solution2
        output.append([puzzle, dot_bracket, solution, accuracy])
        print ([puzzle, dot_bracket, solution, accuracy])
    pickle.dump(output, open('test_results/%s'%(results_path), 'wb'))
    return 0


def refine(dataset, output_path, n_trajs, n_steps, move_set, puzzle_name):
    refined_data = []
    for puzzle in dataset:
        if puzzle_name != '-1' and puzzle[0] != puzzle_name:
            continue
        if puzzle[-1] < 1:
            print ('Trying to refine %s'%(puzzle[0]))
            dot_bracket = puzzle[1]
            input_solution = puzzle[2]
            accuracy, _, _, _ = refine_check_answer(dot_bracket, input_solution)
            for traj in range(n_trajs):
                print ('Trajectory %d'%(traj))
                if accuracy == 1:
                    print ('Found a valid solution')
                    break
                solution = deepcopy(input_solution)
                move_traj = np.random.choice(move_set, n_steps, replace=True)
                for move in move_traj:
                    output, solution = move.apply(dot_bracket, solution)
                    try:
                        output, solution = move.apply(dot_bracket, solution)
                    except:
                        break
                    new_accuracy, _, _ = output
                    if new_accuracy > accuracy:
                        print ('Found better accuracy threshold of %f'%(new_accuracy))
                        print( solution)
                        accuracy = new_accuracy
                        input_solution = solution
                        break
            puzzle_refine = [puzzle[0], dot_bracket, solution, accuracy]
            refined_data.append(puzzle_refine)
        else:
            refined_data.append(puzzle)
    pickle.dump(refined_data, open('refined/%s'%(output_path), 'wb'))
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, argument_default=argparse.SUPPRESS)
    parser.add_argument('--mode', type=str, default='refine', help='Either "train", "test", or "refine"')
    parser.add_argument('--input_data', type=str,default='test_results/out.pkl')
    # parser.add_argument('--input_data', type=str,default='../data/test/eterna100.pkl')
    parser.add_argument('--results_path', type=str,default='refined_0.pkl')
    ######## Train arguments ############
    parser.add_argument('--n_train_puzzles', type=int, default=50, help='How many puzzle to include in the training set.')
    parser.add_argument('--n_solutions_per_puzzle', type=int, default=1, help='How many player solutions per puzzle to use for training.')
    parser.add_argument('--stochastic_fill', type=bool, default=False, help='''When generating the solution trajectory from a player solution, either fill in the puzzle
    sequentially (False) or randomly (True)''')
    parser.add_argument('--renderer', type=str, default='eterna', help='How to render target structures. Choose either "rnaplot" or "eterna"')
    # Long range feature compute arguments
    parser.add_argument('--long_range_input', type=str, default='-1', help='A .pkl file of long-range features to use instead of calculating them.')
    parser.add_argument('--n_long_range_features', type=int, default=20, help='How many long-range features to use as part of the input.')
    # if --long_range_input is supplied, the rest of these are not necessary.
    parser.add_argument('--n_min_solutions', type=int, default=50, help='''A puzzle must have at least this many solutions to be used in the mutual information
    calculation to determine long-range features. Prevents introduction of noise into the calculation when only few player solutions exist.''')
    parser.add_argument('--long_range_output', type=str, default='long_range_features.pkl', help='''What to name the .pkl file of long range features 
    that will be automatically generated during training.''')
    parser.add_argument('--force_add_features', type=bool, default=True, help='''Force the addition of "features_per_puzzle" features to the aggregate list 
    of long-range features for each puzzle. If this is set to False, if the top "features_per_puzzle" features of a puzzle are already part of the aggregate 
    list, that puzzle is skipped and no features are contributed by that puzzle. If set to True, the first non-overlapping "features_per_puzzle" features (ordered by 
    mutual information) will be added.''')
    parser.add_argument('--features_per_puzzle', type=int, default=1, help='''How many long range features each puzzle can contribute.''')
    parser.add_argument('--random_append', type=int, default=0, help='''Before performing the mutual information calculation for a given puzzle, 
    append ("random_append" x number of player solutions) randomly generated sequences to the dataset. This is for instances in which a specific 
    strategy must be used to solve a puzzle, e.g. multiloop coordinated GC pair ordering, that result in identical solutions for the multiloop, 
    resulting in mutual information between all multiloop bases to be 0. By adding randomly generated solutions, we force recognition of correlation 
    between these bases and picking up of the multiloop features.''')
    parser.add_argument('--MI_tolerance', type=float, default=1e-5, help='''When training on a new puzzle using a set of long-range features, any feature found in the puzzle within this cutoff from a feature in the set will be registered as that feature.''')
    # Neural network params
    parser.add_argument('--n_layers', type=int, default=3, help='Number of hidden layers.')
    parser.add_argument('--hidden_size', type=int, default=100, help='Number of nodes in each hidden layer.')
    parser.add_argument('--n_epochs', type=int, default=1000, help='How many epochs to train the model.')
    parser.add_argument('--checkpoint_length', type=int, default=100, help='How many epochs to train before performing a validation.')
    ########## Test arguments ##############
    parser.add_argument('--test_model', type=str, default='../models/trained_models/matlot/test/MI-54_trial4',help='Name of trained model')
    parser.add_argument('--test_puzzle_name', type=str, default='Mat - Lot 2-2 B', help='''Name of the puzzle to use to test the model. If not supplied, 
    the model will be tested on all puzzles in input_data.''')
    ######### Refine arguments ############
    parser.add_argument('--n_trajs', type=int, default=300, help='''Number of refinement trajectories.''')
    parser.add_argument('--n_steps', type=int, default=30, help='''Number of moves per trajectory.''')
    parser.add_argument('--move_set', type=str, default= "[GoodPairsMove([['G', 'C'], ['C', 'G']]), MissingPairsMove([['G', 'C'], ['C', 'G']])]", help='''Moveset used during refinement:
    1) GoodPairsMove(): re-pair two correctly paired bases   
    2) BadPairsMove(): unpair two incorrectly paired bases
    3) MissingPairsMove(): pair two bases that should be paired but are currently unpaired
    4) BoostMove(): boost using G or U-U boosts.''')
    parser.add_argument('--refine_puzzle_name', type=str, default='Mat - Lot 2-2 B', help='''Name of puzzle to refine. If not supplied, will refine 
    all puzzles in input_data.''')
    args = parser.parse_args()
    if '.pkl' in args.input_data:
        try:
            input_data = pickle.load(open(args.input_data))
        except:
            # 使用 with 语句确保文件正确关闭
            with open(args.input_data, 'rb') as file:  # 注意 'rb' 模式
                input_data = pickle.load(file)

    else:
        input_data = load_txt(args.input_data)
    if args.mode == 'train':
        os.system('mkdir results test test/%s'%(args.results_path))
        unique_puzzles = []
        # Create a list of unique puzzle names
        for i in input_data:
            if i[0] not in unique_puzzles:
                unique_puzzles.append(i[0])
        # Create a dictionary of number of solutions per puzzle
        puzzle_solution_count = {}
        for i in input_data:
            if i[0] not in puzzle_solution_count.keys():
                puzzle_solution_count[i[0]] = 1
            else:
                puzzle_solution_count[i[0]] += 1
        train_puzzles = unique_puzzles[:-2]
        val_puzzle = unique_puzzles[-2] 
        test_puzzle = unique_puzzles[-1]
        to_train_on = random.sample(range(len(train_puzzles)), args.n_train_puzzles)
        print("Training Starts Here")
        train_puzzles = [train_puzzles[i] for i in to_train_on]  
        if args.long_range_input == '-1':
            MI_features_master = compute_MI_features(input_data, unique_puzzles, puzzle_solution_count, args.n_min_solutions, args.features_per_puzzle, args.force_add_features, args.random_append, args.renderer)
            try:
                pickle.dump(MI_features_master, open(args.long_range_output, 'w'))
            except:
                pickle.dump(MI_features_master, open(args.long_range_output, 'wb'))

        else:
            MI_features_master = pickle.load(open(args.long_range_input))
        np.random.shuffle(MI_features_master)
        MI_features_list = MI_features_master[:args.n_long_range_features]
        inputs, labels, rewards = parse_progression_dataset(input_data, train_puzzles, args.n_solutions_per_puzzle, MI_features_list, evaluate=False, shuffle=args.stochastic_fill, train_on_solved=True, MI_tolerance=args.MI_tolerance, renderer=args.renderer)
        nn_training_dataset = [inputs, labels, rewards, None]
        val_dataset = parse_progression_dataset(input_data, [val_puzzle], 1, MI_features_list, evaluate=True)
        test_dataset = parse_progression_dataset(input_data, [test_puzzle], 1, MI_features_list, evaluate=True)
        input_size = len(nn_training_dataset[0][0])
        train(args.results_path, args.n_layers, args.hidden_size, args.n_epochs, args.checkpoint_length, MI_features_list, nn_training_dataset, val_dataset, test_dataset, \
              args.MI_tolerance, args.renderer)
    elif args.mode == 'test':
        os.system('mkdir test_results')
        test(input_data, args.test_model, args.results_path, args.test_puzzle_name, args.MI_tolerance, args.renderer)
    elif args.mode == 'refine':
        os.system('mkdir refined')
        move_set = eval(args.move_set)
        refine(input_data, args.results_path, args.n_trajs, args.n_steps, move_set, args.refine_puzzle_name)
    # Removing unnecessary files generated during run
    os.system('rm input rna.ps') 
