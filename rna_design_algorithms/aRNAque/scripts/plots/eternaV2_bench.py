
import numpy as np
import pandas as pd
import json
import ast
import RNA






def main() :

    with open("../../data/Eterna100/V2/eterna100_vienna2.txt") as ft :
        eternaV2 = ft.readlines()
    data_eternaV2 = []
    for line in eternaV2[1:] :
        dt = line.split("\t")
        data_eternaV2.append([dt[1],dt[4], dt[6]])
        print(dt[1], dt[6])

    eterna_V2 = pd.DataFrame(data_eternaV2, columns=["Puzzle Name", "Target", "solvers"])
    print(eterna_V2)

    arnaque_first = []
    targets_first = []
    benchmark_first = []
    with open('../../data/Eterna100/V2/EternaV2_Levy.out') as f :
        for line in f :
            if line.startswith("('Solving for the target") :
                tag = line.split(',')[-1]
                tag = tag.replace(" ", "")
                print(tag[1:-3])
                targets_first.append (tag[1:-3])
            if line.startswith('[(') :
                arnaque_first.append(ast.literal_eval(line))

    print(len(arnaque_first), len(targets_first))

    for i, list_ in enumerate(arnaque_first) :
        try :
            for line in list_ :
                benchmark_first.append([line[0][0], line[0][1],
                                     RNA.hamming_distance(targets_first[i],line[0][1]),
                                     targets_first[i], len(targets_first[i]), 'aRNAque GC2'])
        except KeyError:
            benchmark_first.append([line[0][0], line[0][1],
                                 RNA.hamming_distance(targets_first[i],line[0][1]),
                                 targets_first[i], len(targets_first[i]), 'aRNAque GC2'])

    first_solved = []
    benchmark_first = np.array(benchmark_first)
    print(len(benchmark_first))
    for dt in benchmark_first.reshape((99,20,6)) :
        for d in dt :
            if d[2] == '0' :
                first_solved.append(dt)
                break

    arnaque_second = []
    targets_second = []
    benchmark_second = []
    with open('../../data/Eterna100/V2/EternaV2_LevyGC1.out') as f :
        for line in f :
            if line.startswith("('Solving for the target") :
                tag = line.split(',')[-1]
                tag = tag.replace(" ", "")
                print(tag[1:-3])
                targets_second.append (tag[1:-3])
            if line.startswith('[(') :
                arnaque_second.append(ast.literal_eval(line))

    print(len(arnaque_second), len(targets_second))

    for i, list_ in enumerate(arnaque_second) :
        for line in list_ :
            benchmark_second.append([line[0][0], line[0][1],
                                 RNA.hamming_distance(targets_second[i],line[0][1]),
                                 targets_second[i], len(targets_second[i]), 'aRNAque GC1'])

    benchmark_second = np.array(benchmark_second)
    for dt in benchmark_second.reshape((27,20,6)) :
        for d in dt :
            if d[2] == '0' :
                first_solved.append(dt)
                break


    arnaque_third = []
    targets_third = []
    benchmark_third = []
    with open('../../data/Eterna100/V2/EternaV2_LevyGC1_5000.out') as f :
        for line in f :
            if line.startswith("('Solving for the target") :
                tag = line.split(',')[-1]
                tag = tag.replace(" ", "")
                print(tag[1:-3])
                targets_third.append (tag[1:-3])
            if line.startswith('[(') :
                arnaque_third.append(ast.literal_eval(line))

    print(len(arnaque_third), len(targets_third))

    for i, list_ in enumerate(arnaque_third) :
        for line in list_ :
            benchmark_third.append([line[0][0], line[0][1],
                                 RNA.hamming_distance(targets_third[i],line[0][1]),
                                 targets_third[i], len(targets_third[i]), 'aRNAque GC1'])

    benchmark_third = np.array(benchmark_third)
    for dt in benchmark_third.reshape((23,20,6)) :
        for d in dt :
            if d[2] == '0' :
                first_solved.append(dt)
                break


    arnaque_op = []
    targets_op = []
    benchmark_op = []
    with open('../../data/Eterna100/V2/EternaV2_OP.out') as f :
        for line in f :
            if line.startswith("('Solving for the target") :
                tag = line.split(',')[-1]
                tag = tag.replace(" ", "")
                print(tag[1:-3])
                targets_op.append (tag[1:-3])
            if line.startswith('[(') :
                arnaque_op.append(ast.literal_eval(line))

    print(len(arnaque_op), len(targets_op))

    for i, list_ in enumerate(arnaque_op) :

        for line in list_ :
            if line == None :
                continue
            benchmark_op.append([line[0][0], line[0][1],
                                 RNA.hamming_distance(targets_op[i],line[0][1]),
                                 targets_op[i], len(targets_op[i]), 'OP GC1'])

    benchmark_op = np.array(benchmark_op)
    op_success = []
    for dt in benchmark_op.reshape((27,20,6)) :
        for d in dt :
            if d[2] == '0' :
                op_success += [dt]
                first_solved.append(dt)
                break

    arnaque_op2 = []
    targets_op2 = []
    benchmark_op2 = []
    with open('../../data/Eterna100/V2/EternaV2_OP2.out') as f :
        for line in f :
            if line.startswith("('Solving for the target") :
                tag = line.split(',')[-1]
                tag = tag.replace(" ", "")
                print(tag[1:-3])
                targets_op2.append (tag[1:-3])
            if line.startswith('[(') or line.startswith('[None') :
                arnaque_op2.append(ast.literal_eval(line))

    for i, list_ in enumerate(arnaque_op2) :

        for line in list_ :
            if line == None :
                benchmark_op2.append([None, None,
                                     None,
                                     targets_op2[i], len(targets_op2[i]), 'OP GC1'])
                continue
            benchmark_op2.append([line[0][0], line[0][1],
                                 RNA.hamming_distance(targets_op2[i],line[0][1]),
                                 targets_op2[i], len(targets_op2[i]), 'OP GC1'])

    benchmark_op2 = np.array(benchmark_op2)
    op_success2 = []
    for dt in benchmark_op2.reshape((13,20,6)) :
        #print(dt)
        for d in dt :
            if d[2] == 0 :
                op_success2 += [dt]
                first_solved.append(dt)
                break


    print(op_success2)

    solved_targets = set()
    arnaque_solutions = []

    for lin in first_solved :
        for l in lin :
            if l[2] == '0' or l[2]==0:
                solved_targets.add(l[1])
                arnaque_solutions +=[[eterna_V2[eterna_V2["Target"]==l[1]].values[0][0]]+l.tolist()]



    print("Solved targets: ", len(solved_targets))
    print(arnaque_solutions[0])
    arnaque_solutions.sort(key=lambda elt: len(elt[4]))
    df = pd.DataFrame(arnaque_solutions, columns=["Puzzle Name", "aRNAque Solution","Hamming distance", "RNAfold structure", "Target Structure", "Length", "aRNAque parameter"])
    print (df)
    df.to_csv('../../data/Eterna100/V2/aRNAque_solutions.csv')

    """
    solutions = {}
    clean_solutions = []
    for key in solved_targets :
        solutions[key] = []
        for lin in first_solved :
            for l in lin :
                if key not in l :
                    pass
                elif l[2] == '0' or l[2]==0:
                    solutions[key].append(l[0])


    print(len(solutions) )
    print([len(solutions[key]) for key in solutions.keys()])
    with open("../../data/Eterna100/V2/arnaque_solution.json", 'w') as js_file :
        json.dump(solutions,js_file)
        js_file.close()

    for key in solutions.keys() :
        assert len(key) == len(solutions[key][0]), print(key, len(key), len(solutions[key][0]))
        if len(solutions[key])<5 :
            print(key, solutions[key])
    """


if __name__ == "__main__" :
    main()
