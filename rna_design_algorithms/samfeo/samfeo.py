import os
import sys
import time
import heapq
import argparse

import numpy as np
import pandas as pd

from rna_design_algorithms.samfeo.utils.vienna import position_ed_pd_mfe, position_ed_ned_mfe, mfe
from rna_design_algorithms.samfeo.utils.structure import extract_pairs, struct_dist
from rna_design_algorithms.samfeo.utils.constants import P1, P2, U1, U2

# from multiprocessing import Pool, cpu_count


name2pair = {'cg': ['CG', 'GC'],
             'cggu': ['CG', 'GC', 'GU', 'UG'],
             'cgau': ['CG', 'GC', 'AU', 'UA'],
             'all': ['CG', 'GC', 'AU', 'UA', 'GU', 'UG']}

nuc_others = {'A': 'CGU',
              'C': 'AGU',
              'U': 'ACG',
              'G': 'ACU'}

nuc_pair_others = {'AU': ['UA', 'CG', 'GC', 'UG', 'GU'],
                   'UA': ['AU', 'CG', 'GC', 'UG', 'GU'],
                   'CG': ['AU', 'UA', 'GC', 'UG', 'GU'],
                   'GC': ['AU', 'UA', 'CG', 'UG', 'GU'],
                   'GU': ['AU', 'UA', 'CG', 'GC', 'UG'],
                   'UG': ['AU', 'UA', 'CG', 'GC', 'GU']}

nuc_all = ['A', 'C', 'G', 'U']
nuc_pair_all = ['AU', 'UA', 'CG', 'GC', 'UG', 'GU']

STAY = 2000
STOP = 0.01
MAX_REPEAT = 1000


class RNAStructure:

    def __init__(self, seq, score, v,
                 v_list):  # v_list: positional NED, v: objective value, socore: used for priority queue
        self.seq = seq
        self.score = score
        self.v = v
        self.v_list = v_list

    def __gt__(self, other):
        return self.score > other.score

    def __lt__(self, other):
        return self.score < other.score

    def __eq__(self, other):
        return self.seq == other.seq

    def __ge__(self, other):
        return self.score >= other.score

    def __le__(self, other):
        return self.score <= other.score

    def __str__(self):
        return f"{self.seq}: {self.score: .4f}, {1 - self.score: .4f}"

    def __repr__(self):
        return f"RNAStructure('{self.seq}', {self.score})"

    def __hash__(self):
        return hash(self.seq)


def init_with_pair(t, pos_pairs, pairs_init):
    rna = list("." * len(t))
    assert len(rna) == len(t)
    for i, s in enumerate(t):
        if s == ".":
            rna[i] = 'A'
            if name_pair == 'all':
                rna[i] = np.random.choice(['A', 'C', 'G', 'U'])
        elif s == "(":
            j = pos_pairs[i]
            pair = np.random.choice(pairs_init)
            rna[i] = pair[0]
            rna[j] = pair[1]
        elif s == ")":
            pass
        else:
            raise ValueError(f'the value of structure at position: {i} is not right: {s}!')
    return "".join(rna)


# targeted initilization
def init_k(target, pos_pairs, k):
    # print(f'name_pair: {name_pair}')
    pair_pool = name2pair[name_pair]
    # print(f'pair_pool: {pair_pool}')
    init_0 = init_with_pair(target, pos_pairs, pair_pool)
    p_list = [init_0]
    # if too few pairs then use 'cggu', however this may never happen
    if k > len(pair_pool) ** (len(pos_pairs) / 2) and len(pair_pool) < 4:
        pair_pool = name2pair['cggu']
    # the max number of intial sequences is: len(pair_pool)**(len(pos_pairs)/2)
    while len(p_list) < min(k, len(pair_pool) ** (len(pos_pairs) / 2)):
        init_i = init_with_pair(target, pos_pairs, pair_pool)
        if init_i not in p_list:
            p_list.append(init_i)
    return p_list


def pairs_match(ss):  # find the pairs in a secondary structure, return a dictionary
    assert len(ss) > 5
    pairs = dict()
    stack = []
    for i, s in enumerate(ss):
        if s == ".":
            pass
        elif s == "(":
            stack.append(i)
        elif s == ")":
            j = stack.pop()
            assert j < i
            pairs[j] = i
            pairs[i] = j
        else:
            raise ValueError(f'the value of structure at position: {i} is not right: {s}!')
    return pairs


def mutate_pair(nuc_i, nuc_j, exclude=False):
    pair_ij = nuc_i + nuc_j
    return np.random.choice(nuc_pair_others[pair_ij]) if exclude else np.random.choice(nuc_pair_all)


def mutate_unpair(nuc_i, exclude=False):
    return np.random.choice(list(nuc_others[nuc_i])) if exclude else np.random.choice(nuc_all)


# traditional mutation
def mutate_tradition(seq, pairs, v, v_list, T, pairs_dg=None):
    v_list = [v / T for v in v_list]
    probs = np.exp(v_list) / sum(np.exp(v_list))
    index = np.random.choice(list(range(len(seq))), p=probs)
    seq_next = [nuc for nuc in seq]
    if index in pairs:
        i = min(index, pairs[index])
        j = max(index, pairs[index])
        pair_ij = seq[i] + seq[j]
        pair_new = np.random.choice(nuc_pair_others[pair_ij])
        seq_next[i] = pair_new[0]
        seq_next[j] = pair_new[1]
    else:
        c = np.random.choice(list(nuc_others[seq[index]]))
        assert c != seq[index]
        seq_next[index] = c
    return "".join(seq_next)


# structured mutation
def mutate_structured(seq, pairs, v, v_list, T):
    v_list = [v / T for v in v_list]
    probs = np.exp(v_list) / sum(np.exp(v_list))
    index = np.random.choice(list(range(len(seq))), p=probs)
    pairs_mt = []
    unpairs_mt = []

    if index in pairs:
        i = min(index, pairs[index])
        j = max(index, pairs[index])
        pairs_mt.append((i, j))
        if j - 1 in pairs and pairs[j - 1] == i + 1:
            pairs_mt.append((pairs[j - 1], j - 1))
            if i + 2 not in pairs and j - 2 not in pairs:
                unpairs_mt.append(i + 2)
                unpairs_mt.append(j - 2)
        if i + 1 not in pairs and j - 1 not in pairs:
            unpairs_mt.append(i + 1)
            unpairs_mt.append(j - 1)
    else:
        unpairs_mt.append(index)
        if index - 1 in pairs and pairs[index - 1] > index:
            pairs_mt.append((index - 1, pairs[index - 1]))
            if pairs[index - 1] - 1 not in pairs:
                unpairs_mt.append(pairs[index - 1] - 1)
        elif index + 1 in pairs and pairs[index + 1] < index:
            pairs_mt.append((pairs[index + 1], index + 1))
            if pairs[index + 1] + 1 not in pairs:
                unpairs_mt.append(pairs[index + 1] + 1)

    assert len(pairs_mt) <= 2, pairs_mt
    assert len(unpairs_mt) <= 2, unpairs_mt

    # one pair
    if len(pairs_mt) == 1:
        pairs_selected_index = np.random.choice(range(len(P1)))
        pairs_selected = P1[pairs_selected_index]
    else:  # two pair
        pairs_selected_index = np.random.choice(range(len(P2)))
        pairs_selected = P2[pairs_selected_index]

    # one unpair
    if len(unpairs_mt) == 1:
        unpairs_selected_index = np.random.choice(range(len(U1)))
        unpairs_selected = U1[unpairs_selected_index]
    else:  # two unpair
        unpairs_selected_index = np.random.choice(range(len(U2)))
        unpairs_selected = U2[unpairs_selected_index]

    nuc_list = list(seq)
    for pos_pair, pair in zip(pairs_mt, pairs_selected):
        nuc_list[pos_pair[0]] = pair[0]
        nuc_list[pos_pair[1]] = pair[1]
    for pos_unpair, unpair in zip(unpairs_mt, unpairs_selected):
        nuc_list[pos_unpair] = unpair
    return "".join(nuc_list)


def samfeo(target, f, steps, k, t=1, check_mfe=True, sm=True, freq_print=10):
    global seed_np
    np.random.seed(seed_np)
    # print(f'seed_np: {seed_np}')
    if sm:
        mutate = mutate_structured
    else:
        mutate = mutate_tradition
    # print(f'steps: {steps}, t: {t}, k: {k}, structured mutation: {sm}, ensemble objective: {f.__name__}')

    # targeted initilization
    pairs = pairs_match(target)
    intial_list = init_k(target, pairs, k)
    history = set()
    k_best = []
    log = []
    dist_list = []
    mfe_list = []
    count_umfe = 0
    for p in intial_list:
        v_list, v, ss_list = f(p,
                               target)  # v_list: positional NED, v: objective value, ss_list: (multiple) MFE structures by subopt of ViennaRNA
        rna_struct = RNAStructure(seq=p, score=1 - v, v=v, v_list=v_list)
        rna_struct.dist = min([struct_dist(target, ss_subopt) for ss_subopt in ss_list])  # ss: secondary structure
        rna_struct.subcount = len(ss_list)
        k_best.append(rna_struct)
        history.add(rna_struct.seq)

    # priority queue
    heapq.heapify(k_best)
    for i, rna_struct in enumerate(k_best):
        # print(i, rna_struct)
        log.append(1 - rna_struct.score)
        if rna_struct.dist == 0:  # MFE solution
            mfe_list.append(rna_struct.seq)
        if rna_struct.dist == 0 and rna_struct.subcount == 1:  # UMFE solution
            dist_list.append(-2)
            count_umfe += 1
        else:
            dist_list.append(rna_struct.dist)

    # log of lowest objective value at eachs iterations
    v_min = min(log)
    iter_min = 0
    log_min = [v_min]
    for i in range(steps):
        # sequence selection
        score_list = [rna_struct.score / t * 2 for rna_struct in k_best]  # objective values
        probs_boltzmann_1 = np.exp(score_list) / sum(np.exp(score_list))  # boltzmann distribution
        try:
            p = np.random.choice(k_best, p=probs_boltzmann_1)
        except Exception as e:
            print(e)
            p = np.random.choice(k_best)

        # position sampling and mutation
        seq_next = mutate(p.seq, pairs, p.v, p.v_list, t)
        num_repeat = 0
        while seq_next in history:
            num_repeat += 1
            if num_repeat > len(target) * MAX_REPEAT:
                break
            p = np.random.choice(k_best, p=probs_boltzmann_1)
            seq_next = mutate(p.seq, pairs, p.v, p.v_list, t)
        if num_repeat > len(target) * MAX_REPEAT:
            # print(f'num_repeat: {num_repeat} > {len(target) * MAX_REPEAT}')
            break
        history.add(seq_next)

        # evaluation new sequence
        v_list_next, v_next, ss_list = f(seq_next, target)

        # mfe and umfe solutions as byproducts
        umfe = False
        if check_mfe:
            dist = min([struct_dist(target, ss_subopt) for ss_subopt in ss_list])
            if dist == 0:
                mfe_list.append(seq_next)
                if len(ss_list) == 1:
                    umfe = True
        else:
            dist = len(target)  # set a dummy dist
        if not umfe:
            dist_list.append(dist)
        else:
            dist_list.append(-2)
            count_umfe += 1

        # update priority queue(multi-frontier)
        rna_struct_next = RNAStructure(seq_next, 1 - v_next, v_next, v_list_next)

        if len(k_best) < k:
            heapq.heappush(k_best, rna_struct_next)
        elif rna_struct_next > k_best[0]:
            heapq.heappushpop(k_best, rna_struct_next)
        if v_next <= v_min:
            iter_min = i

        # update log
        v_min = min(v_min, v_next)
        log_min.append(v_min)
        log.append(v_next)
        assert len(dist_list) == len(log)

        # # output information during iteration
        # if (i + 1) % freq_print == 0:
        #     improve = v_min - log_min[-freq_print]
        #     if check_mfe:
        #         print(
        #             f"iter: {i + 1: 5d}\t value: {v_min: .4f}\t mfe count: {len(mfe_list): 5d}\t umfe count: {count_umfe}\t best iter: {iter_min} improve: {improve:.2e}")
        #     else:
        #         print(f"iter: {i + 1: 5d}\t value: {v_min: .4f}\t best iter: {iter_min} improve: {improve:.2e}")

        # stop if convergency condition is satisfied
        if v_min < STOP or (len(log_min) > STAY and v_min - log_min[-STAY] > -1e-6):
            break
    return k_best, log, mfe_list, dist_list


# RNA design in batch
def design(target, func, num_step, k, t, check_mfe, sm):
    # start samfeo
    start_time = time.time()
    k_best, log, mfe_list, dist_list = samfeo(target, func, num_step, k=k, t=t, check_mfe=check_mfe, sm=sm)  # rna and ensemble defect
    end_time = time.time()
    time_consume = end_time - start_time

    rna_best = max(k_best)
    seq = rna_best.seq
    ss_mfe = mfe(seq)[0]
    dist = struct_dist(target, ss_mfe)
    return seq, ss_mfe, target, dist, time_consume


class Args:
    def __init__(self, target, object='pd', k=10, t=1, step=5000, name='', init='cg', repeat=1, start=0, nomfe=False,
                 nosm=False, bp=False, online=False):
        self.target = target
        self.object = object
        self.k = k
        self.t = t
        self.step = step
        self.name = name
        self.init = init
        self.repeat = repeat
        self.start = start
        self.nomfe = nomfe
        self.nosm = nosm
        self.bp = bp
        self.online = online


def call_samfeo(target, object='pd', k=10, t=1, step=5000, init='cg', repeat=1, start=0, nomfe=False, nosm=False, save_file="results/samfeo_result"):

    global name_pair, stop, seed_np
    name_pair = init

    if object == 'ned':  # normalized ensemble defect
        f_obj = position_ed_ned_mfe
    elif object == 'pd':  # probability defect
        f_obj = position_ed_pd_mfe
    else:
        raise ValueError('the objective in not correct!')

    output_Seq = []  # designed sequences
    ouput_Dis = []  # distances
    time_consumes = []
    for i in range(repeat):
        # generate seed
        seed_np = 2020 + (i + start) * 2021
        np.random.seed(seed_np)
        seq, ss_mfe, target, dist, time = design(target, f_obj, step, k=k, t=t, check_mfe=not nomfe, sm=not nosm)
        ouput_Dis.append(dist)
        output_Seq.append(seq)
        time_consumes.append(time)

    # save to file
    data = {'target_structure': target, 'sequence': output_Seq, 'time': time_consumes, 'distance': ouput_Dis}
    df = pd.DataFrame(data)
    save_file = save_file + '.pkl'
    df.to_pickle(save_file)

    return output_Seq


if __name__ == "__main__":
    target_structure = "(((....)))((...))"
    repeat = 10
    output_seq = call_samfeo(target_structure, repeat=repeat)

    for i in range(repeat):
        print(output_seq[i])

    df = pd.read_pickle("results/samfeo_result.pkl")
    pass
