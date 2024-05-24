"""

    ######################################################################################################
         A SIMPLE EVOLUTIONARY ALGORITHM GUIDED BY LOCAL AND LEVY MUTATION FOR INVERSE RNA FOLDING
    ######################################################################################################

    @author: Nono Saha Cyrille Merleau, email: nonosaha@mis.mpg.de/csaha@aims.edu.gh
"""
#!/usr/bin/env python3
#Importting necessary python libraries
from pandas import DataFrame
from numpy import where, array, random, round, concatenate
from numpy.random import  choice as npchoice
from os import rmdir, path, mkdir,getcwd
from time import time
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, SUPPRESS
from utilities.Landscape import Landscape
from utilities.utility import (get_bp_position, getHairpinCoord, ppeval, ppfold,
                    gen_point_mutation_dist, boost_hairpins, hamming, _defect,
                    save_population, checkTarget)
from multiprocess import Pool, cpu_count
from json import load
from datetime import datetime
from shutil import rmtree

#GLOBAL Variables
base_paire = ["GC","CG","AU","UA", "GU", "UG"]
nucleotides = ["A", "G","U","C"]
ROOT_LOG_FOLDER = getcwd()+"/../data/Log"


# Initialize a population of RNA.
def init(pop_size, landscape, constraints=None, npseed=None) :

    target = landscape.target
    pos = get_bp_position(target)
    bp_pos = pos['bp']  + pos['bp']
    init_depth =len(target)

    nucleotides = ["A", "U" , "G", "C"]
    base_paire = ["GC","CG","AU","UA", "GU", "UG"]
    pop = []
    i = 0
    if npseed != None :
        random.seed(npseed)
    while i < pop_size :

        if i < 4:
            arn = npchoice(nucleotides[i:i+1],init_depth)
        else :
            arn = npchoice(nucleotides,len(target))

        for bp_cord in bp_pos :
            bp = npchoice(base_paire,1)
            arn[bp_cord[0]] = bp[0][0]
            arn[bp_cord[1]] = bp[0][1]

        if constraints!=None :
                arn[constraints['pos']] = constraints['sequence']
        pop.append(''.join(arn))
        i = len(pop)

    return pop


#Mutation function
def mutate(pop, pos, p_n, p_c,  distribution,constraints=None,loop_boosting=True) :
    bp_pos = pos["bp"] + pos['pk']
    nbp_indexes = pos['nbp']
    mutants = []



    bp_points = distribution['bp']
    nb_points = distribution['nbp']

    for i in range(len(pop)) :
        mutant = array(list(pop[i]))

	    #Mutate a non base pair position
        nbp_indexes_to_mutate = npchoice(list(nbp_indexes), nb_points[i], replace=False)
        mutant[nbp_indexes_to_mutate]= npchoice(nucleotides, nb_points[i],p=p_n)


        #Mutate a base pair position
        if len(bp_pos) > 0 :

            bp_choice = npchoice(base_paire,bp_points[i], p=p_c)
            bp_indexes_to_mutate = npchoice(range(len(bp_pos)), bp_points[i],replace=False)
            bp_to_mutate = array(bp_pos)[bp_indexes_to_mutate]

            for j in range(bp_points[i]) :
                mutant[bp_to_mutate[j][0]] =bp_choice[j][0]
                mutant[bp_to_mutate[j][1]] =bp_choice[j][1]
            if loop_boosting :
                for hp_pos in pos["hairpins"] :
                    mutant = boost_hairpins(mutant, hp_pos)


        if constraints!=None :
            mutant[constraints['pos']] = constraints['sequence']

        mutants.append("".join(mutant))

    return mutants

# Selection function
def select(population, size, method='F') :

    if method == 'F' :
        weights  = array(population["Fitness"], dtype = float)
    elif method == 'NED' :
        evals = array(population["Evals"], dtype = float)
        mfes  = array(population["Mfes"], dtype = float)
        delta_mfe = evals - mfes
        delta_mfe = delta_mfe / sum(delta_mfe)
        weights = (1-delta_mfe)**100
    elif method =='ED' :
        weights = 1./array(population["ED"], dtype = float)
    else :
        print ('The method '+method+ 'does not exist')
        return;

    selected = npchoice(list(population["RNA_sequence"]),size=size,p=weights/sum(weights))
    return list(selected)


#Elitism method that aims to copy the number of fittest sequences to the next generation
def elite(population, size, weight='Fitness') :
    sorted_pop = population.sort_values(by=weight, ascending=False)

    return sorted_pop.values[:size]



#Main Evolutionary algorithm
def simple_EA(param) :
    """
        Add comment here: TODO
    """
    if param['verbose'] :
        print (" Starting of evolution ")

    prev_population = param['init_pop'].copy(deep=True) #Initialize the population of RNA
    population_size =  param['n']

    t = 0
    max_found = len(set(prev_population[prev_population["Fitness"]==str(1.0)]["RNA_sequence"].values))
    b = prev_population.sort_values(by=["Fitness"],ascending=False).values[0]

    if param['log'] :
        root_path = '../log/'+param['process_id']
        if not path.exists(root_path) :
            mkdir(root_path)


    while (t<param["time"]) and (max_found < param['msf']):
        if param['log'] :
            save_population(prev_population, t, root_path)
        if param['verbose'] :
            print("Generation {}\n{}\n{}\t{}\t{}".format(t,b[0],b[1],b[2],b[3]))

        newgeneration = []
        best = elite(prev_population, 10, weight='Fitness')
        selected_ind = select(prev_population,population_size, param['sm'])

        dist = gen_point_mutation_dist(population_size,param['pos'],param['c'])
        mutated = mutate(selected_ind,param['pos'], param['p_n'], param['p_c'],
                            dist,param['constraints'],param['lb'])+ list(best[:,0])

        newgeneration.append(mutated)

        structures, mfes = ppfold(mutated,param['ftool'], param["nrj_param"])
        newgeneration.append(structures)
        newgeneration.append(mfes)
        newgeneration.append(hamming(structures,param['landscape']))
        if param['sm'] == 'NED' :
            evals = ppeval(mutated, param['landscape'].target)
        else :
            evals = [0]*len(mutated)
        newgeneration.append(evals)
        prev_population = DataFrame(array(newgeneration).T, columns=["RNA_sequence", "RNA_structure", "Mfes", "Fitness","Evals"])
        b = prev_population.sort_values(by=["Fitness"],ascending=False).values[0]
        current_fitest = max(array(prev_population["Fitness"], dtype=float))

        t +=1
        max_found = len(set(prev_population[prev_population["Fitness"]==str(1.0)]["RNA_sequence"].values))

    return prev_population, t


def ED_refinement(param) :

    print (" Starting ED refinement ")
    prev_population = param['init_pop'].copy(deep=True)


    population_size =  param['n']
    t = param['EDg']
    landscape = param['landscape']

    while (t > 0) :

        newgeneration = []

        best = elite(prev_population, 10,weight='Fitness')
        max_fitness = max(array(prev_population['Fitness'], dtype=float))
        selected_ind = select(prev_population,population_size,'ED')

        dist = gen_point_mutation_dist(population_size,param['pos'],param['c'])
        mutated = mutate(selected_ind,param['pos'], param['p_n'], param['p_c'],dist,param['constraints'])+list(best[:,0])
        newgeneration.append(mutated)

        structures, mfes = ppfold(mutated,param['ftool'], param['nrj_param'])
        newgeneration.append(structures)
        newgeneration.append(mfes)
        newgeneration.append(hamming(structures,landscape))

        defects = _defect(mutated, landscape)
        newgeneration.append(defects)

        prev_population = DataFrame(array(newgeneration).T, columns=["RNA_sequence", "RNA_structure", "Mfes", "Fitness","ED"])
        seq_found = prev_population[prev_population["Fitness"]==str(max_fitness)]
        if len(seq_found.values) > 0 :
            prev_population = seq_found

        t -=1

    return prev_population



def run(param) :

    random.seed()
    tic = time()
    best_pop, Mgen = simple_EA(param)

    if param['EDg']>0 :
        ens_defects = _defect(list(best_pop["RNA_sequence"]), landscape)
        init_pop = DataFrame(array([list(best_pop["RNA_sequence"]), list(best_pop["RNA_structure"]), list(best_pop["Mfes"]), list(best_pop["Fitness"]), ens_defects]).T, columns=["RNA_sequence", "RNA_structure", "Mfes", "Fitness","ED"])
        param['init_pop'] = init_pop
        best_pop = ED_refinement(param)
        founds = best_pop[best_pop['Fitness']==str(max(array(best_pop['Fitness'], dtype=float)))]

        founds = founds.sort_values(by="ED", ascending=True)
        rst = []
        for i, val in enumerate(founds["RNA_sequence"].unique()) : 
            if i>param["msf"]: 
                break
            rst += [founds[founds["RNA_sequence"]==val].values[0].tolist()]
            if param['verbose'] :
                print (founds[founds["RNA_sequence"]==val].values[0], Mgen)

    else :
        #sorted_pop = population.sort_values(by="Fitness", ascending=False)
        #founds = sorted_pop[:param["msf"]]
        founds = best_pop[best_pop['Fitness']==str(max(array(best_pop['Fitness'], dtype=float)))]
        founds = founds.sort_values(by='Mfes', ascending=False)
        rst = []
        for i, val in enumerate(founds["RNA_sequence"].unique()) : 
            if i>param["msf"]: 
                break
            rst += [founds[founds["RNA_sequence"]==val].values[0].tolist()]
            if param['verbose'] :
                print (founds[founds["RNA_sequence"]==val].values[0], Mgen)
    toc = time()
    if param['verbose']:
        print ("Time = ", toc-tic, "For method : ", param['sm'])
    return (rst, Mgen)


def parseArguments() :

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter, argument_default=SUPPRESS)
    parser.add_argument('--target', '-t', type=str, help='Target RNA secondary structure in dot bracket representation')
    parser.add_argument('--job','-j', type=int,default=1, help="Number of EA runs")
    parser.add_argument('-g', type=int, default=150, help="Number of generation")
    parser.add_argument('-n', type=int,default=100, help="Population Size")
    parser.add_argument('-msf', type=int,default=1, help="maximum sequence found")
    parser.add_argument('-sm', type=str,default ='NED',help="Selection method: the only possible values are {F,NED}")
    parser.add_argument('-bp', type=str, default='GC2',help="Distribution of nucleotide and base pairs. Possible values are {GC,GC1,GC2,GC3,GC4,GC25, GC50, GC75,ALL}, please check the online doc for more details")
    parser.add_argument('--Cs', type=str, default=None,help="sequence constraints: the lenght of the sequence should be the same as the target. Example: target=((....)), C=GNNNANNC")
    parser.add_argument('-EDg', type=int, default=0, help="number of generation for Ensemble defect refinement")
    parser.add_argument('-c', type=float,default=None, help="Exponent of the zipf's distribution" )
    parser.add_argument('--hairpin_boosting', action="store_true", default=False, help="Boost the hairpin loops. When false no hairpins boosting" )
    parser.add_argument('--folding_tool','-ft',type=str,default='v', help="folding tool to be used: v for RNAfold from viennarna, ip for IPknot and hk for Hotknots" )
    parser.add_argument('--log', action="store_true", default=False, help="Store the population for each instance of the inverse folding in a folder" )
    parser.add_argument('--verbose', action="store_true", default=False, help="Print the mean fitness evolution on a standard output " )
    parser.add_argument('--turner1999', action="store_true", default=False, help="Use the old energy parameters" )

    parser.add_argument('-seed', type=int, default=None, help="Seed for the initial population")
    return parser.parse_args()


def main() :

    args = parseArguments()

    if not path.exists(ROOT_LOG_FOLDER+"/tmp/") :
        mkdir(ROOT_LOG_FOLDER+'/tmp')
    
    if args.folding_tool not in {"v", "ip", "hk"} : 
        print("Please make sure the folding tool entered is correct. i.e -ft should be in {v, ip, hk}") 
        exit(1)
    if args.folding_tool in {"ip", "hk"}: 
        if args.EDg > 0 : 
            print("WARNING: No ens. def. refinement for pseudoknotter target")
            args.EDg = 0
        if args.sm == "NED" : 
            print("WARNING: Only hamming distance selection function allowed for pseudoknotted target")
            args.sm = "F"

    if args.bp not in {"GC","GC1","GC2","GC3","GC4","GC25", "GC50", "GC75","ALL"} : 
        print("ERROR: Please, enter a correct base-pair distribution parameter. use --help for help.")
        exit(1)
    if args.sm not in {"NED", "F"} : 
        print("ERROR: Please, enter a correct selection function. e.g. F: for hamming distance and NED: for energy distance")
        exit(1)

    if args.target is not None :
        target = args.target
        if (checkTarget(target)==False): 
           print("Incorrect target secondary structure.")
           exit(1)
        if ('[' in list(target)) or (']' in list(target)) or ('{' in list(target)) or ('' in list(target)): 
            if args.folding_tool == 'v': 
                print("Please, choose the appropriate foolding tool.")
                exit(1)
    main_sequence = args.Cs
    constraints = {}
    if main_sequence !=None :
        if len(main_sequence) == len(target) :
            if set(list(main_sequence)).issubset(set(["A","C","G","U","N"])) :
                constraints["pos"] = where(array(list(main_sequence))!="N")
                constraints["sequence"] = array(list(main_sequence))[constraints["pos"]]
            else :
                print('Constraint sequence error: Please check the help using "python aRNAque.py --help"')

    else :
        constraints = None


    landscape = Landscape(target)
    mut_params = {
        'ALL' : {
            'p_n' : [0.7,0.1,0.1,.1],
            'p_c' : [0.3, 0.2, 0.2, 0.1,0.1,0.1]
        },
        'GC' : {
            'p_n' : [0.25,0.25,0.25,0.25],
            'p_c' : [0.5,0.5,0.0,0.0,0.0,0.0]
        },
        'GC1' : {
            "p_n" : [0.7,0.1,0.1,.1],
            "p_c" : [0.3, 0.2, 0.2, 0.1,0.1,0.1]
        },
        'GC2' : {
            "p_n" : [0.7,0.1,0.1,.1],
            "p_c" : [0.4, 0.5, 0.1, 0.,0.,0.]
        },
        'GC3' : {
            "p_n" : [0.75,0.1,0.1,.05],
            "p_c" : [0.4, 0.5, 0.1, 0.,0.,0.]
        },
        'GC4' : {
            "p_n" : [0.95,0.0,0.05,0.0],
            "p_c" : [0.4, 0.4, 0.2, 0.,0.,0.]
        },
        'GC5' : {
            "p_n" : [0.25,0.65,0.05,0.05],
            "p_c" : [0.4,0.5,0.1,0.0,0.0,0.0]
        },
        'UN' : {
            'p_n' : None,
            'p_c' : None
        },
        'b' :  {
            'p_c': [0.284226552051662, 0.03570581512544118, 0.18006763282289678, 0.284226552051662, 0.03570581512544118, 0.18006763282289678],
            'p_n': [0.6700050529739228, 0.11993631562036493, 0.10168862490269953, 0.10837000650301278]
        },
    }

    best_mutation_params = {}
    try :

        with open(ROOT_LOG_FOLDER+'/../../params/mutation/best_params.json') as js_params :
            best_mutation_params = load(js_params)
            js_params.close()
    except IOError:
        pass

    for key in best_mutation_params.keys() :
        mut_params[key] = {}
        mut_params[key]['p_c'] = best_mutation_params[key][0]
        mut_params[key]['p_n'] = best_mutation_params[key][1]

    pos = get_bp_position(target)
    pos["hairpins"] =  []


    params = []
    log_folder = ''
    if args.log  :
        log_folder = str(datetime.now()).replace(" ", "") + '/'
        mkdir("../log/"+log_folder)

    for i in range(args.job) :
        if args.log :
            pop = ['A'*len(target)]*args.n
        else :
            pop = init(args.n, landscape, constraints=constraints, npseed=args.seed)

        if args.sm == 'NED' :
            evals = ppeval(pop, target)
        else :
            evals = [0]*len(pop)

        strcs, mfes = ppfold(pop,args.folding_tool, args.turner1999)
        fitnesses = hamming(strcs, landscape)
        init_pop = DataFrame(array([pop, strcs, mfes, fitnesses, evals]).T, columns=["RNA_sequence", "RNA_structure", "Mfes", "Fitness","Evals"])

        if len(pos['hairpins']) == 0 :
            pos["hairpins"] =  getHairpinCoord(pop[0],target)



        params += [{
            'init_pop' : init_pop,
            'ftool' : args.folding_tool,
            'p_c' : mut_params[args.bp]['p_c'],
            'p_n' : mut_params[args.bp]['p_n'],
            'lb': args.hairpin_boosting,
            'constraints' : constraints,
            'landscape' : landscape,
            'pos' :pos,
            'c' : args.c,
            'time' : args.g,
            'n': args.n,
            'EDg': args.EDg,
            'sm' : args.sm,
            'msf': args.msf,
            'process_id': log_folder+str(i),
            'log' : args.log,
            'verbose': args.verbose,
            "nrj_param": args.turner1999
        }]



    print ("=====================================================================================================")
    print ("Solving for the target = "+target )
    print ("=====================================================================================================")

    print (len(pos["hairpins"]), " Hairpin loop(s)")
    print (len(pos["bp"])," base pair(s) in total")

    pool = Pool(cpu_count())
    result = pool.map(run,params)
    pool.close()
    #print ("|Sequence | Structure | Fitness | Energy | NED |")
    #result = concatenate(result).tolist()
    #print (result)
    for i, val in enumerate(result):
        print("*"*50)
        print ("Job id:", i, " Number of generations: ", val[-1])
        print("*"*50)

        rst_job = val[0]
        rst_job.sort(key=lambda elt: float(elt[2]))
        for rst in rst_job :
            print (rst[0], "|" ,rst[1] , "|", round(float(rst[2]),2), "|", rst[3])

    try :
        rmtree(ROOT_LOG_FOLDER+'/tmp', ignore_errors=True)
    except OSError as e:
        print("Error while cleaning the tmp directory files", e)




if __name__ == "__main__":
    main()
