from tool import *
from result_log import *
from tool_reactCount import generateReactions
import json
import pickle
import copy
import sys,os
from math import isnan
import global_var as gvar
import networkx as nx
import time

def trackBlocksXYZ():
    Rstep = []
    speRec = []
    tempList=[]    
    printed = []
    istep = 1
    step  = 1
    pbc = gvar.pbcmole
    bnd_cri = 1.5  # Criterion for bond linkage
    fnames = gvar.trajFiles
    bnd_cri = gvar.bnd_cri
    
    # 1. read xyz data in
    for fname in fnames:
        with  open(fname,'r') as f:
            if (istep == 1):
                #gvar.pbcXYZ = PBCRec[istep-1]
                print(gvar.pbcXYZ)
                # First Record:
                p_time = time.time()
                Natom = int(f.readline())
                gvar.atomList = [0]*Natom # BuildUpAtomListi
                temGrpRec= [[0]*2 for i in range(Natom)] 
                f.readline()
                for i in range(Natom):
                    line = f.readline().split()
                    gvar.atomList[i] = \
                    [line[0],[float(line[1]),float(line[2]),float(line[3])],[]]
                   #   ^                 ^                                  ^
                   #  ele.Name,  element coordinate,                  neighbor   
                # Initialize global Mask Matrix
                Element = [row[0] for row in gvar.atomList] 
                radii_array = np.array([*map(gvar.radii_dict.get, Element)])
                radii_array = np.tile(radii_array,(len(Element),1))
                gvar.GlobalMaskMat = (radii_array+radii_array.T)*bnd_cri
                np.fill_diagonal(gvar.GlobalMaskMat,-1)

                # Build atom list.
                grpRec,nMol = buildNeigh_AtomicBased(10) # 20 anstrom into neighbour :::: Can be parallel
                #print(nMol)
                print("--- neighbor %s seconds1 ---" % (time.time() - p_time))
                p_time = time.time()
                gvar.blockList = [[0]*7 for i in range(nMol)] # BuildUpBlockList
                Molid = 1
                # Initial block
                for rec in grpRec:
                    gvar.blockList[Molid-1][0] = Molid
                    gvar.blockList[Molid-1][1] = list(map(int,rec))
                    Molid = Molid+1

                if(pbc):
                    BlockInfoUpdatePBC(bnd_cri)  #:::: Can be parallel
                    printUnknowStruc()           #:::: Can be parallel
                else:
                    # for future develope AIMD-base non-fragment
                    pass
                print(gvar.atomList[12])
                print("--- neighbor %s seconds2 ---" % (time.time() - p_time))
                p_time = time.time()
                P_blockList = copy.deepcopy(gvar.blockList)
            while True:
#               if (debug): print(istep,step,len(blockList))
                line = f.readline()

                # Analyse Finished for one file.  
                if not line:
                    printUnknowStruc()           
                    NodeTranslation("Graph_5000.dot")
                    gvar.DicReactionRec = generateReactions(gvar.GR)
                    PrintMetaDataRaw()
                    break
                Natom =int(line)
                f.readline()
                for i in range(Natom):
                    line = f.readline().split()
                   #gvar.atomList[i] = \
                   #[line[0],[float(line[1]),float(line[2]),float(line[3])],[]]
#
                    gvar.atomList[i][1] =[float(line[1]),float(line[2]),float(line[3])]  
                if(istep%80 == 0):
                    print("number of Mol:",nMol)
                    printUnknowStruc()
                if(istep%1000 == 0):
                    print("--- %s seconds ---" % (time.time() - p_time))
                    p_time = time.time()
                    print(istep)
                if(istep%20 == 0):
                    time_neibor1 = time.time()
                    BlockNeighborUpdate(5)
                    time_neibor2 = time.time()
                    print("Neibor --- %s seconds ---" % (time_neibor2 - time_neibor1))
                if(istep%1 == 0 and istep !=1 ):
                    p_1 = time.time()
                    if(pbc):
                        time_grp1 = time.time()
                        grpRec,nMol = buildDistMart(bnd_cri)
                        time_grp2 = time.time()
                        print("Group --- %s seconds ---" % (time_grp2 - time_grp1))
                    else:
                        pass

                    gvar.blockList = [[0]*7 for i in range(nMol)] # BuildUpBlockList
                    Molid = 1
                    for rec in grpRec:
                        gvar.blockList[Molid-1][0] = Molid
                        gvar.blockList[Molid-1][1] = list(map(int,rec))
                        Molid = Molid+1
                    p_2 = time.time()

                    if(pbc):
                        time_blk1 = time.time()
                        BlockInfoUpdatePBC(bnd_cri)
                        time_blk2 = time.time()
                        print("Block --- %s seconds ---" % (time_blk2 - time_blk1))
                    else:
                        pass
                    p_3 = time.time()

                    # Build Reaction Connection.
                    compare2Step(P_blockList, gvar.blockList, istep-1) 
                    #
                    P_blockList = copy.deepcopy(gvar.blockList)

                istep = 1+istep
    return(tempList)

def trackBlocksLAMMPS():
    Rstep = []
    speRec = []
    istep = 0
    pbc = gvar.pbcmole
    bnd_cri = 1.5  # Criterion for bond linkage
    fnames = gvar.trajFiles
    bnd_cri = gvar.bnd_cri
    COORcount = -1
    Natom = 0
    # 1. read xyz data in
    readFlag = 0 # 1. number of atoms; 2. Box bounds 3. Coord
    for fname in fnames:
        with  open(fname,'r') as f:
            while True:
                line = f.readline()
                # File-read Finished
                if not line:
                    write_dot(gvar.GR,"G1cores.dot")
                    print(gvar.GR.number_of_edges())
                    printUnknowStruc()           
                    NodeTranslation("Graph_5000.dot")
                    gvar.DicReactionRec = generateReactions(gvar.GR)
                    PrintMetaDataRaw()
                    break
                if("NUMBER OF ATOMS" in line): 
                    readFlag = 1
                    istep = istep+1
                    continue

                if("BOUNDS" in line): 
                    readFlag = 2
                    PBCcount = 0
                    pbcXYZList = [[0,0,0],[0,0,0]]
                    continue

                if("x y z" in line): 
                    readFlag = 3
                    lineSList = line.split()
                    lineSList.remove("ATOMS")
                    iatom = lineSList.index('id') - 1
                    ielement = lineSList.index('type') - 1
                    x = lineSList.index('x') - 1
                    y = lineSList.index('y') - 1
                    z = lineSList.index('z') - 1
                    #print(iatom,ielement,x,y,z)
                    COORcount = 0
                    continue
                if(readFlag == 0): continue

                if(readFlag == 1):
                    Natom =int(line)
                    if(istep ==1):gvar.atomList = [0]*Natom
                    readFlag = 0
                    continue

                if(readFlag == 2):
                    pbcXYZList[0][PBCcount] = float(line.split()[0])
                    pbcXYZList[1][PBCcount] = float(line.split()[1])

                    PBCcount = PBCcount +1
                    if(PBCcount==3): 
                        readFlag = 0
                        gvar.pbcXYZ =[tuple(pbcXYZList[0]),tuple(pbcXYZList[1])] 
                    continue
                if(readFlag == 3):
                    CoorRec = line.split()
                    COORcount = COORcount +1
                    element = gvar.lmpAtmDict[int(CoorRec[ielement])]
                    if(istep == 1):
                        gvar.atomList[int(CoorRec[iatom])-1] = \
                        [element,[float(CoorRec[x]),float(CoorRec[y]),float(CoorRec[z])],[]]
                    else:
                        gvar.atomList[int(CoorRec[iatom])-1][1] = \
                        [float(CoorRec[x]),float(CoorRec[y]),float(CoorRec[z])]
                    if(COORcount == Natom):
                        readFlag = 0


                # AtomList RenerFinished. Analysis begin
                if(COORcount == Natom ):
                    if(istep ==1): 
                        # Initialize global Mask Matrix
                        Element = [row[0] for row in gvar.atomList] 
                        radii_array = np.array([*map(gvar.radii_dict.get, Element)])
                        radii_array = np.tile(radii_array,(len(Element),1))
                        gvar.GlobalMaskMat = (radii_array+radii_array.T)*bnd_cri
                        np.fill_diagonal(gvar.GlobalMaskMat,-1)
                        #Build neiborList and calculate first steps of molecules
                        grpRec,nMol = buildNeigh_AtomicBased(10) # 10 anstrom into neighbour :::: Can be parallel

                        #print(nMol)
                        gvar.blockList = [[0]*7 for i in range(nMol)] # BuildUpBlockList
                        Molid = 1
                        # Initial block
                        for rec in grpRec:
                            gvar.blockList[Molid-1][0] = Molid
                            gvar.blockList[Molid-1][1] = list(map(int,rec))
                            Molid = Molid+1

                        if(pbc):
                            BlockInfoUpdatePBC(bnd_cri)  #:::: Can be parallel
                            printUnknowStruc()           #:::: Can be parallel
                        else:
                            # for future develope AIMD-base non-fragment
                            pass
                        print(gvar.atomList[12])
                        P_blockList = copy.deepcopy(gvar.blockList)
                    if(istep%1 == 0 and istep !=1 ):
                        if(pbc):
                            time_grp1 = time.time()
                            grpRec,nMol = buildDistMart(bnd_cri)
                            time_grp2 = time.time()
                        else:
                            pass

                        gvar.blockList = [[0]*7 for i in range(nMol)] # BuildUpBlockList
                        Molid = 1
                        for rec in grpRec:
                            gvar.blockList[Molid-1][0] = Molid
                            gvar.blockList[Molid-1][1] = list(map(int,rec))
                            Molid = Molid+1

                        if(pbc):
                            time_blk1 = time.time()
                            BlockInfoUpdatePBC(bnd_cri)
                            time_blk2 = time.time()
                            print("Block --- %s seconds ---" % (time_blk2 - time_blk1))
                        else:
                            pass

                        # Build Reaction Connection.
                        compare2Step(P_blockList, gvar.blockList, istep-1) 
                        #
                        P_blockList = copy.deepcopy(gvar.blockList)
                    if(istep%80 == 0):
                        print("number of Mol:",nMol)
                        printUnknowStruc()
                    if(istep%1 == 0):
                        time_neibor1 = time.time()
                        BlockNeighborUpdate(5)
                        time_neibor2 = time.time()
                        print("Neibor --- %s seconds ---" % (time_neibor2 - time_neibor1))
                  # print(readFlag,Natom,gvar.pbcXYZ)
                  # print(gvar.atomList[0],istep)
                  # print(gvar.atomList[1],istep)
                #print(CoorRec[ielement],CoorRec[iatom],CoorRec[x],CoorRec[y],CoorRec[z])



    return(tempList)

