
"""
@author: Nono Saha Cyrille Merleau, email: nonosaha@mis.mpg.de/csaha@aims.edu.gh

"""

from pandas import DataFrame
from numpy import array as nparray, ones
from numpy.random import choice as npchoice
from os import remove, system, getcwd
from multiprocess import Pool, cpu_count
from uuid import uuid4
from utilities.folding_wrapper import *


ROOT_LOG_FOLDER = getcwd()+"/../data/Log"

def getHairpinCoord(sequence, target) :
    cmd_input = str(sequence)+"\n"+str(target)+"\n"
    id_ = str(uuid4())
    system("echo '"+cmd_input+"'| RNAeval -v |grep Hairpin | tr -d A-Z-a-z-'('-')' | cut -d ':' -f 1 >"+ROOT_LOG_FOLDER+"/tmp/hairpin"+id_+".loop")
    hairpins = []
    with open(ROOT_LOG_FOLDER+"/tmp/hairpin"+id_+".loop", "r") as loops :
        data_r = loops.read().split("\n")
        for dt in data_r :
            if dt != '' :
                hairpins.append(tuple(nparray(list(map(int, dt.strip().split(','))))-1))
    loops.close()
    remove(ROOT_LOG_FOLDER+"/tmp/hairpin"+id_+".loop")

    return hairpins


def boost_hairpins(sequence, coord) :
    stems = ["A"]*(coord[1]-coord[0]-1)
    if coord[1]-coord[0]-1 >=4 :
        stems[0] = "G"
        stems.insert(0,"G")
        stems.append("C")
        sequence[coord[0]:coord[1]+1] = stems
    else :
        sequence[coord[0]+1:coord[1]] = stems
    return sequence

def hamming(listOfStructures, landscape) :
    dists = [landscape.fitness(strc) for strc in listOfStructures]
    return dists

def _defect(listOfSequences, landscape) :
    eds = [landscape.ens_defect(s) for s in listOfSequences]
    return eds

def get_bp_position(structure) :
    pk_pairs, bp, nbp= [] , [], []

    pairs = {
        'bp' : [],
        'pk' : [],
        'nbp' : []
    }
    for i,elt in enumerate(structure) :
        if elt =='[' :
            pk_pairs.append(i)
        elif elt ==']' :
            pairs['pk'] += [(pk_pairs.pop(),i)]

        elif elt == '(' :
            bp += [i]
        elif elt == ')' :
            pairs['bp'] += [(bp.pop(),i)]
        else :
            pairs['nbp'] += [i]
    return pairs

#Logging population
def bt_save_population(prev_pop, population,gen, root_path) :
    data = []
    prev_data = []
    for i in range(len(population)) :
        data.append([population[i].RNA_seq, population[i].RNA_structure,population[i].mfe, population[i].fitness])
        prev_data.append([prev_pop[i].RNA_seq, prev_pop[i].RNA_structure,prev_pop[i].mfe, prev_pop[i].fitness])


    dataFrame = DataFrame(data)
    prev_dataFrame = DataFrame(prev_data)
    prev_dataFrame.to_csv(root_path+"/prev_gen"+str(gen)+".csv")
    dataFrame.to_csv(root_path+"/gen"+str(gen)+".csv")


def save_population(population,gen, root_path, file_prefix="") :

    population.to_csv(root_path+"/gen"+file_prefix+str(gen)+".csv")



#Evaluate the energy
def ppeval(listOfSeqs, target) :
    task = uuid4()
    with open(ROOT_LOG_FOLDER+"/tmp/rnaeval_in"+str(task), "w") as file_ :
        for seq in listOfSeqs :
            file_.write(seq.strip()+"\n"+target.strip()+"\n")
        file_.close()

    system("RNAeval -j --infile="+ROOT_LOG_FOLDER+"/tmp/rnaeval_in"+str(task)+" |tr -d A-Z,'(',')'|cut -d ' ' -f 2- >"+ROOT_LOG_FOLDER+"/tmp/result_"+str(task))
    with open(ROOT_LOG_FOLDER+"/tmp/result_"+str(task), "r") as file_ :
        eval_ = file_.read().split()
    remove(ROOT_LOG_FOLDER+"/tmp/result_"+str(task))
    remove(ROOT_LOG_FOLDER+"/tmp/rnaeval_in"+str(task))


    return nparray(eval_, dtype=float).tolist()

def nthHarmonic(N,s) :

    harmonic = 1.00
    for i in range(2, N + 1) :
        harmonic += 1 / i**s

    return harmonic

def zipf_rvgen(low, high, N, size, a=6.) :
    choices = nparray(range(low, high+1))
    probs = choices**(-a) / nthHarmonic(N,a)
    return npchoice(choices,p=nparray(probs)/sum(probs),size=size)

def gen_point_mutation_dist(size,pos,c):

    bp_pos = len(pos['bp']+pos['pk'])
    nbp_pos = len(pos['nbp'])
    dist = {}
    if c !=None :
        if 0<c<7.5:
            if bp_pos > 0 :
                dist['bp'] = zipf_rvgen(1,bp_pos, bp_pos, size, c)
            else :
                dist['bp'] = []
            dist['nbp'] = zipf_rvgen(1,nbp_pos, nbp_pos, size, c)
        else :
            dist['bp'] = modulo_dist(1,bp_pos, size)
            dist['nbp'] = modulo_dist(1,nbp_pos,size)
    else :
        if bp_pos > 0 :
            dist['bp'] = ones(size, dtype=int)
        else :
            dist['bp'] = []
        dist['nbp'] = ones(size, dtype=int)

    return dist


def ppfold(listOfSeqs,tool, nrj_param) :

    if tool =="v" :
        return ppRNAfold(listOfSeqs, nrj_param)
    if tool =="pk" :
        return ppKiss(listOfSeqs)
    if tool =="ip" :
        return ppipknot(listOfSeqs)
    if tool =="hk" :
        return ppHotknots(listOfSeqs)



def checkTarget(target) : 
    # check the number pairs 
    target_list = list(target)

    post = get_bp_position(target)

    if (target_list.count('(') != target_list.count(')')) or (target_list.count('{') != target_list.count('}')) or (target_list.count('[') != target_list.count(']')) : 
        print("ERROR: Number of opened and closed brackets not equal.")
        return False 
    for ch in set(target_list) : 
        if ch not in ['(', ')', '.', '[', ']', '{', '}'] : 
            print ("ERROR: Character ", ch, " is not allowed in the string secondary structure representation.")
            return False
    for p in post["bp"] + post["pk"] : 
        if p[1] - p[0]-1 < 3 : 
            print ("EROOR: Incorrect base-pair at position ", p)
            return False
