import os 
import numpy as np 
import multiprocess as mp 
import argparse
import time as timelib
def ant_exe(params) :
    cmd = "python antarna.py -tGC {}  -t {} MFE -Cstr '{}' -p --pkprogram={} "
        
    target = params['target']
    time = params['time']
    tGC = params['tGC']
    ft = params['ft']
    st = timelib.time()
    p = os.popen(cmd.format(tGC, time,target, ft))
    rst = p.read().split()
    spt = timelib.time()

    print("TEST ", rst)
    
    if len(rst) > 0 : 
        return (rst[-1], spt-st)
    else : 
        return None


def pp_antRNA(target,time, tGC, nb_runs,  ft="pKiss") :

    st = timelib.time()
    pool = mp.Pool(mp.cpu_count())
    result = pool.map(ant_exe, [{"target": target,
                                "time": time, 
                                "tGC": tGC, 
                                "ft": ft, 
                                "taskID": i} for i in range(nb_runs) ])
    spt = timelib.time()

    pool.close() 

    return {target: result, "time": spt-st}


def parse_arguments():
    """Parsing command line
    """
    parser = argparse.ArgumentParser(description="Wrapper for antaRNA program: used for benchmark ")
    parser.add_argument('--number_of_runs', '-nr', default=1, help="Number of sequences to be designed")
    parser.add_argument('--target', '-t', help="Target secondary structure")    
    parser.add_argument('--tGC', '-tGC',type=float, default=0.5, help="GC content")    
    parser.add_argument('--ft', '-ft',type= str ,help="Folding tool ", default="pKiss")
    parser.add_argument('--time', '-it', type=int, default=1200, help="Number of iterations")

    return parser.parse_args()

def main() :
    args = parse_arguments()
    args.target = "...(((...)))..."
    args.ft = "pKiss"
    print (pp_antRNA(args.target,args.time, args.tGC, int(args.number_of_runs), args.ft))
    print("testhehe")



if __name__=="__main__" : 
    main()


