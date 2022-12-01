from tool import *
from result_log import *
#from AnalysFlow import *
from AnalysFlow_parallel import *
from AnalysFlow_parallel_cata import *
from tool_reactCount import generateReactions
import json
import pickle
import copy
import sys,os
from math import isnan
import global_var as gvar
import networkx as nx
import time
from inputio import renderSetting


# Print the warning information to NULL
import os # if you have not already done this

if(os.path.exists("./errorAndwarning.log")):
    os.remove("./errorAndwarning.log")

os.close(os.open("./errorAndwarning.log", os.O_CREAT))
fd = os.open('./errorAndwarning.log',os.O_WRONLY)
os.dup2(fd,2)


if __name__ == "__main__":
    """
    CPX-MechGen: 
    Discover reaction mechanism from complex reaction network.

    -Author      : Chen Xin
    -Email       : chenxin199261@gmail.com
    -Create Date : 2021/1/19 
    -Inputs:
        1. Multiple xyz trajectory files. fill their paths in variable 'fname'.

    -Outputs:
        1. blk_evo  : The evolution of monitoring fragments. 
        2. graph.dat: Graph metaData 

    -Data structure:
        tempList1: atoms belonging to which group and group id
                   tempList1 = [atomRec*nstep]
                   atomRec[i] for atom i is [blknum,grp_id]

    """
    p_start0 = time.time()
    Inpfile = "inputfile"
    renderSetting(Inpfile)
    print(time.strftime("StartTime:  %Y-%m-%d %H:%M:%S", time.localtime()))

    # 0.  Initializing 
    debug =  True
    #analy_mode = 1      # 1.raw xyz data to results 
                         # 2.raw xyz data to tempRec
                         # 3.temRecs to results
    

    # get network cleaning lag.
    if(gvar.stableMolTime[1] == gvar.stableMolTime[1]):
        gvar.StableMolLag = \
        int(float(gvar.stableMolTime[0])/float(gvar.timePerFrame[0]))
    elif(gvar.stableMolTime[1] == 'ps' and gvar.timePerFrame[1]=='fs'):
        gvar.StableMolLag = \
        int(float(gvar.stableMolTime[0])/float(gvar.timePerFrame[0])*1000)
    else:
        gvar.StableMolLag = 10
        print("Warning, your trajectory has large time lag")
    print("StableCriLag is ",gvar.StableMolLag)
    

    os.environ["PYTHONHASHSEED"] = "12222"
    gvar.GR = nx.MultiDiGraph()


    bnd_cri = 1.5  # Criterion for bond linkage
    for key1 in gvar.radicri_names:
        for key2 in gvar.radicri_names[key1]:
            gvar.radicri_names[key1][key2] *= bnd_cri

    # Build Mask matrix and catalyst atoms(Maybe) for futher analysis
    if(gvar.trajTyp =="xyz"):
        BuildMaskForXYZ(gvar.trajFiles[0])

    if(gvar.trajTyp =="lmp"):
        BuildMaskForLAMMPS(gvar.trajFiles[0])



#   if(gvar.trajTyp == 'lmp'):  trackBlocksLAMMPS()
#   if(gvar.trajTyp == 'xyz'):  trackBlocksXYZ()
    CountFilesAndLines()
    if(gvar.cata):
        gvar.CataAtom.sort() 
        print("cata")
        parallelPool_cata()
    else:
        parallelPool()
    print(time.strftime("EndTime:  %Y-%m-%d %H:%M:%S", time.localtime()))
    p_end0 = time.time()
    print("Total seconds:",p_end0-p_start0)

# +++++++++++++++++++++++++++++++++++++++++++++++++++
#   if (analy_mode == 2 or analy_mode == 1):
#       tempList1 = trackBlocksXYZ(fnames,trackList)
#   if (analy_mode == 3):
#       tempList1=[]
#       for ftname in ftempnames:
#           print(ftname)
#           with open(ftname, 'rb') as f:
#               tl = pickle.load(f)
#               tempList1.extend(tl)
#   if(analy_mode !=2 ):
#       reactionGen(tempList1)
