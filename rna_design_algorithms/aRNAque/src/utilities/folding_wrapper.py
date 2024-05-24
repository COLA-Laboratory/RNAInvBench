"""
@author: Nono Saha Cyrille Merleau, email: nonosaha@mis.mpg.de/csaha@aims.edu.gh
"""

import os
from numpy import array
from RNA import fold, read_parameter_file
from uuid import uuid4

ROOT_LOG_FOLDER = os.getcwd()+"/../data/Log"

def ppRNAfold(listOfSeqs, nrj_param) :
    if (nrj_param) :
       read_parameter_file("../params/energy/rna_turner1999.par")

    rst = array ([list(fold(seq)) for seq in listOfSeqs])
    return rst[:,0].tolist(), array(rst[:,1],dtype=float).tolist()

def fold_with_ipknot(seq) :

    out_file = str(uuid4())
    cmd= "echo '>"+out_file+"\n"+seq+"'>"+ROOT_LOG_FOLDER+"/tmp/"+out_file+".fa"
    os.system(cmd)
    ipknot_cmd = os.environ.get("IPKNOT")+"/./ipknot "+ROOT_LOG_FOLDER+"/tmp/"+out_file+".fa"
    p = os.popen(ipknot_cmd)
    rst = p.read().split()
    p.close()
    os.remove(ROOT_LOG_FOLDER+"/tmp/"+out_file+".fa")
    if len(rst) > 0 :
        return rst[-1]
    else :
        print("ERROR during the folding with ipknot")
        return None

def ppipknot(listOfSeqs) :
    result = []
    for s in listOfSeqs :
        result.append(fold_with_ipknot(s))
    result = array(result)

    return result, [0]*len(result)

def pKiss(seq) :

    cmd= "pKiss --strategy 'P' --mode='mfe' {} 2>/dev/null".format(seq)
    p = os.popen(cmd)
    rst = p.read().split()
    if len(rst) > 0 :
        return (rst[-1],rst[-2])
    else :
        print("ERROR during the folding with pKiss")
        return None

def hotknots(seq) :

    cmd= os.environ.get("HOTKNOTS_ROOT")+"./bin/HotKnots -s  {} -m {} 2>/dev/null".format(seq, "CC")
    p = os.popen(cmd)
    rst = p.read().split('\n')
    rst = rst[2].split()
    if len(rst) > 0 :
        return (rst[-2],rst[-1])
    else :
        print("ERROR during the folding with Hotknots")
        return None

def ppHotknots(listOfSeqs) :

    result = []
    for s in listOfSeqs :
        result.append(hotknots(s))
    result = array(result)
    return list(result[:,0]),list(result[:,1])
