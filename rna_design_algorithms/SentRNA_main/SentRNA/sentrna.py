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
import os
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'
import tensorflow as tf


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

class Args:
    def __init__(self, mode='test', input_data='../data/test/eterna100.pkl', results_path='model_tested.pkl',
                 n_train_puzzles=50, n_solutions_per_puzzle=1, stochastic_fill=False, renderer='eterna',
                 long_range_input='-1', n_long_range_features=20, n_min_solutions=50,
                 long_range_output='long_range_features.pkl', force_add_features=True, features_per_puzzle=1,
                 random_append=0, MI_tolerance=1e-5, n_layers=3, hidden_size=100, n_epochs=1000,
                 checkpoint_length=100, test_model='../models/trained_models/matlot/test/MI-54_trial4',
                 test_puzzle_name='-1', n_trajs=300, n_steps=30, move_set='[GoodPairsMove(),BadPairsMove(),MissingPairsMove(),BoostMove()]',
                 refine_puzzle_name='-1'):
        # Initialize all parameters
        self.mode = mode
        self.input_data = input_data
        self.results_path = results_path
        self.n_train_puzzles = n_train_puzzles
        self.n_solutions_per_puzzle = n_solutions_per_puzzle
        self.stochastic_fill = stochastic_fill
        self.renderer = renderer
        self.long_range_input = long_range_input
        self.n_long_range_features = n_long_range_features
        self.n_min_solutions = n_min_solutions
        self.long_range_output = long_range_output
        self.force_add_features = force_add_features
        self.features_per_puzzle = features_per_puzzle
        self.random_append = random_append
        self.MI_tolerance = MI_tolerance
        self.n_layers = n_layers
        self.hidden_size = hidden_size
        self.n_epochs = n_epochs
        self.checkpoint_length = checkpoint_length
        self.test_model = test_model
        self.test_puzzle_name = test_puzzle_name
        self.n_trajs = n_trajs
        self.n_steps = n_steps
        self.move_set = move_set
        self.refine_puzzle_name = refine_puzzle_name


        # Continue initializing other parameters as above...

def main(args):
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
    # os.system('rm input rna.ps')

# args = Args(
#         mode='test',
#         input_data='../data/test/eterna100.pkl',
#         results_path='model_tested.pkl',
#     )
#
# main(args)

# model_name = 'refine'
# input_data = 'test_results/model_tested.pkl'
# results_path = 'model_refined.pkl'
# args = Args(
#         mode=model_name,
#         input_data=input_data,
#         results_path=results_path,
#         # Set other parameters as required...
#     )
# main(args)