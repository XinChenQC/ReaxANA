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
from multiprocessing import Pool,current_process

#==========================================
# Count files and lines in the analysis.  #
#                                         #
#==========================================

def CountFilesAndLines():
    fnames = gvar.trajFiles
    TotLines = 0
    ifile = 0
    linesPerFrame = 0
    FramePerFile = []  
    for fname in fnames:
        TotLinesFile = 0
        with  open(fname,'r') as f:
            #Count lines of 1 traj step:
            if(ifile == 0 and linesPerFrame ==0):
                # Count frameLenth for xyz files 
                if(gvar.trajTyp == 'xyz'):  
                    line = f.readline()
                    linesPerFrame = int(line)+2

                # Count frameLenth for lmp trajectories  files 
                if(gvar.trajTyp == 'lmp'):
                    HookMark = f.readline()
                    linesPerFrame = linesPerFrame+1
                    while True:
                        line = f.readline()
                        if(line == HookMark): break
                        linesPerFrame = linesPerFrame+1

                f.seek(0)
            
            #Count Total Lines and step:
            while True:
                line = f.readline()
                if not line:
                    break
                TotLines = TotLines+1
                TotLinesFile = TotLinesFile+1

        ifile = ifile+1
        FramePerFile.append(TotLinesFile)
        f.close()

    # Count finished
    if(TotLines%linesPerFrame != 0): 
        print("Warning: trajectory file error")
    else:
        TotalFrames = int(TotLines/linesPerFrame)
        print("Your trajectory  have  ",TotalFrames," frames")
        gvar.StepsinTotal = TotalFrames
        for irec in range(len(FramePerFile)):
            FramePerFile[irec] = int(FramePerFile[irec]/linesPerFrame)
    if(gvar.ncores == 1):
        gvar.ProcessArrangement = [(0,0,len(FramePerFile)-1,-1,0)]
        gvar.StepRecPerCore = [(1,TotalFrames)]
        return
    FramesPerCore = int(TotalFrames/float(gvar.ncores))+1
    FirstCoreFrames = int(FramesPerCore*1.2)
    LastCoreFrames = int(FramesPerCore*0.8)
    FramesArrangeRec = []
    for icore in range(gvar.ncores-1):
        FramesArrangeRec.append(FirstCoreFrames+\
        int(icore*(LastCoreFrames-FirstCoreFrames)/(gvar.ncores-1)))
    LastCoreFrames = TotalFrames - sum(FramesArrangeRec)
    FramesArrangeRec.append(LastCoreFrames)
    print(FramesArrangeRec)
    CutPosition = [] #[(fileName,lineNum),(fileName,lineNum)]
    ifile = 0
    icount = 0
    icore = 1
    Acc = FramesArrangeRec[0]
    for fileCount in FramePerFile:
        localCount = 0
        for i in range(fileCount):
            icount = icount+1
            localCount = localCount+1
            if(icount == Acc and icount<TotalFrames):
                print(icount,fnames[ifile],localCount)
                Acc = FramesArrangeRec[icore]+Acc 
                icore = icore + 1
                if(localCount == fileCount):
                    CutPosition.append((ifile,localCount*linesPerFrame,"N"))
                else:
                    CutPosition.append((ifile,localCount*linesPerFrame,"C"))
        ifile = ifile+1

    # Assign trajectories to cores
    StartPoint = 1
    for icore in range(gvar.ncores):
        if (icore == gvar.ncores-1):
            startFrame = StartPoint
            endFrame = TotalFrames 
        else:
            startFrame = StartPoint
            endFrame   = StartPoint+FramesArrangeRec[icore]
            StartPoint = StartPoint+FramesArrangeRec[icore]+1
        gvar.StepRecPerCore.append((startFrame,endFrame))
        if icore == 0:
            # From firstfile firstLine (0,0)
            ProcessArrangement=[(0,0,CutPosition[icore][0],CutPosition[icore][1],icore)]
        elif(icore == gvar.ncores-1):
            # To lastfile lastLine (0,0)
            ProcessArrangement.append((CutPosition[icore-1][0],CutPosition[icore-1][1],\
                                       len(gvar.trajFiles)-1,-1,icore))
        else:
            if(CutPosition[icore-1][2] == "C"):
                ProcessArrangement.append((CutPosition[icore-1][0],CutPosition[icore-1][1],\
                                           CutPosition[icore][0],  CutPosition[icore][1],icore))
            else:
                ProcessArrangement.append((CutPosition[icore][0],0,\
                                           CutPosition[icore][0],CutPosition[icore][1],icore))
    gvar.ProcessArrangement = ProcessArrangement
    print(gvar.ProcessArrangement)
    print(gvar.StepRecPerCore)

# Read files parallel

def parallelReadXYZ(Assign):
    fstart = Assign[0] 
    fend   = Assign[2]
    startLine = Assign[1]
    endLine   = Assign[3]
    BlockRec_first = None
    BlockRec_last  = None
    time.sleep(Assign[4]*90.0)
    P_Files = []
    time_Ini = time.time()
    for ifile in range(fstart,fend+1):
        P_Files.append(gvar.trajFiles[ifile])
    print(current_process().name," processing files: ",\
          P_Files," processing  ",gvar.StepRecPerCore[Assign[4]]," step ")

    iline = 1
    istep = 1
    pbc = gvar.pbcmole
    bnd_cri = gvar.bnd_cri
    COORcount = -1
    Natom = 0
    StepRec = gvar.StepRecPerCore[Assign[4]]
    ifileCount = 0
    for fname in P_Files:
        print(fname,istep,current_process().name)
        fileLineCount = 0
        ifileCount = ifileCount+1
        with  open(fname,'r') as f:
            # Switch to the reading head.
            while True:
                fileLineCount += 1
                if(ifileCount>1): break
                if(ifileCount==1 and (fileLineCount < startLine+1) ):
                    f.readline()
                    continue 
                if(ifileCount==1 and (fileLineCount == startLine+1) ):
                    print(Assign,"start reading",current_process().name)
                    break
                f.readline()
            while True:
                # Jump first line
                line = f.readline()
                # Exit 
                if not line:
                    # ReadFinal
                    if(ifileCount == len(P_Files)):
                        print("file Last typ2")
                        BlockRec_last  = copy.deepcopy(gvar.blockList)
                    else:
                        print("NotLast")
                    print("Last line",fileLineCount,Assign,"fileend")
                    print("Last line:",line,"fileend")
                    break
                Natom = int(line.strip('\n') )
                if(istep==1): gvar.atomList = [0]*Natom # BuildUpAtomListi
                fileLineCount += 1

                line = f.readline()
                fileLineCount += 1

                innerLine = 0
                # Read Atoms line by line
                while innerLine < Natom:
                    line = f.readline().strip().split()
                    fileLineCount += 1
                    # Build current atomlist
                    if(istep ==1):
                        gvar.atomList[innerLine] = \
                        [line[0],[float(line[1]),float(line[2]),float(line[3])],[]]
                    else:
                        gvar.atomList[innerLine][1] =[float(line[1]),float(line[2]),float(line[3])]  
                    innerLine += 1
                if(istep == 1):
                    grpRec,nMol = buildNeigh_AtomicBased(10) # 20 anstrom into neighbour :::: Can be parallel
                else:
                    grpRec,nMol = buildDistMart(bnd_cri)
                print(nMol)
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
                    if(istep == 1): BlockRec_first = copy.deepcopy(gvar.blockList)
                else:
                    # for future develope AIMD-base non-fragment
                    pass

                if(istep >1):
                    compare2Step(P_blockList, gvar.blockList, istep-2+StepRec[0]) 

                # Last Line
                if(ifileCount==len(P_Files) and fileLineCount == endLine+1):
                    BlockRec_last  = copy.deepcopy(gvar.blockList)
                    break

                if(istep%5 == 0):
                    time_neibor1 = time.time()
                    BlockNeighborUpdate(10)
                    time_neibor2 = time.time()
                    #print("Neibor --- %s seconds ---" % (time_neibor2 - time_neibor1))
                if(istep%100 == 0):
                    time_tmp = time.time()
                    print((time_tmp - time_Ini)," seconds pass in calculating recent 100 steps")
                    time_Ini = time_tmp


                # Build Reaction Connection.

                P_blockList = copy.deepcopy(gvar.blockList)
                #

              #   # print(Assign[3],fileLineCount)
              #   # print(istep,len(gvar.SpeciesCount))
                istep += 1
    print(BlockRec_first[0],"first")
    print(BlockRec_last[0],"last")
    NodeTranslation("Graph_1000.dot")
    ReactionCleanSub()
    return([Assign[4], gvar.GR, gvar.DicStuct, gvar.DicMoleInfo, gvar.SpeciesCount,BlockRec_first,BlockRec_last])



def parallelReadLAMMPS(Assign):
    fstart = Assign[0] 
    fend   = Assign[2]
    startLine = Assign[1]
    endLine   = Assign[3]
    time.sleep(Assign[4]*30.0)
    P_Files = []
    for ifile in range(fstart,fend+1):
        P_Files.append(gvar.trajFiles[ifile])
    time_process0 = time.time()
    print(current_process().name," processing files: ",\
          P_Files," processing  ",gvar.StepRecPerCore[Assign[4]]," step ")

    istep = 0
    pbc = gvar.pbcmole
    bnd_cri = 1.5  # Criterion for bond linkage
    bnd_cri = gvar.bnd_cri
    COORcount = -1
    Natom = 0
    Nnode = 0
    ifileCount = 0
    BlockRec_first = None
    BlockRec_last  = None
    # 1. read xyz data in
    readFlag = 0 # 1. number of atoms; 2. Box bounds 3. Coord
    for fname in P_Files:
        print("Processing:",fname)
        fileLineCount = 0
        ifileCount = ifileCount+1
        with  open(fname,'r') as f:
            while True:
                line = f.readline()
                fileLineCount += 1
                #print(ifileCount,fileLineCount,Assign,"T")
                if(ifileCount==1 and (fileLineCount < startLine+1) ):
                    continue 
                #if(ifileCount==1 and (fileLineCount == startLine+1) ):
                    #print(line,Assign)

                # Last Line
                if(ifileCount==len(P_Files) and fileLineCount == endLine+1):
                    BlockRec_last  = gvar.blockList
                    #print(Assign[3],fileLineCount,line)
                    #print(istep,len(gvar.SpeciesCount))
                    break 
                # File-read Finished

                if not line:
                    # Not the last line
                    if(ifileCount==len(P_Files)):
                        BlockRec_last  = gvar.blockList
                        #print(istep,len(gvar.SpeciesCount))
                        
                   #printUnknowStruc()           
                   #NodeTranslation("Graph_5000.dot")
                   #gvar.DicReactionRec = generateReactions(gvar.GR)
                   #PrintMetaDataRaw()
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
                        time0 = time.time()                                      #Time
                        Element = [row[0] for row in gvar.atomList] 
                        radii_array = np.array([*map(gvar.radii_dict.get, Element)])
                        radii_array = np.tile(radii_array,(len(Element),1))
                        gvar.GlobalMaskMat = (radii_array+radii_array.T)*bnd_cri
                        np.fill_diagonal(gvar.GlobalMaskMat,-1)
                        #Build neiborList and calculate first steps of molecules
                        grpRec,nMol = buildNeigh_AtomicBased(10) # 10 anstrom into neighbour :::: Can be parallel
                        time1 = time.time()                                      #Time
                        gvar.Time_DC = gvar.Time_DC + (time1-time0)

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
                        BlockRec_first = gvar.blockList
                        P_blockList = copy.deepcopy(gvar.blockList)
                        time2 = time.time()                                      #Time
                        gvar.Time_MI = gvar.Time_MI + (time2-time1)
                    if(istep%1 == 0 and istep !=1 ):
                        time0 = time.time()                                      #Time
                        if(pbc):
                            time_grp1 = time.time()
                            grpRec,nMol = buildDistMart(bnd_cri)
                            time_grp2 = time.time()
                        else:
                            pass
                        time1 = time.time()                                      #Time
                        gvar.Time_DC =  gvar.Time_DC + (time1-time0)

                        gvar.blockList = [[0]*7 for i in range(nMol)] # BuildUpBlockList
                        Molid = 1
                        for rec in grpRec:
                            gvar.blockList[Molid-1][0] = Molid
                            gvar.blockList[Molid-1][1] = list(map(int,rec))
                            Molid = Molid+1

                        if(pbc):
                            #time_blk1 = time.time()
                            BlockInfoUpdatePBC(bnd_cri)
                            #time_blk2 = time.time()
                            #print("Initial Block time --- %s seconds ---" % (time_blk2 - time_blk1))
                        else:
                            pass
                        time2 = time.time()                                      #Time
                        gvar.Time_MI =  gvar.Time_MI + (time2-time1)

                        # Build Reaction Connection.
                        compare2Step(P_blockList, gvar.blockList, istep-2+gvar.StepRecPerCore[Assign[4]][0]) 
                        time3 = time.time()                                      #Time
                        gvar.Time_BW = gvar.Time_BW + (time3-time2)
                        #
                        P_blockList = copy.deepcopy(gvar.blockList)
                    if(istep%10 == 0):
                        time0 = time.time()                                      #Time
                        printUnknowStruc()
                        time1 = time.time()                                      #Time
                        gvar.Time_TS = gvar.Time_TS + (time1-time0)
                    if(istep%100 == 0):
                        time_100step = time.time()
                        print("Recent 100 steps use:", time_100step-time_process0, " seconds")
                        print("number of Mol:",nMol," at step: ",istep)
                        time_process0 =  time_100step
                      # NodeC = gvar.GR.number_of_nodes()
                      # if(NodeC - Nnode > 1000):
                      #     NodeTranslation("Graph_1000.dot")
                      #     ReactionCleanGreyTrans()
                      # Nnode = NodeC
                    if(istep%1 == 0):
                        time0 = time.time()                                      #Time
                        BlockNeighborUpdate(15)
                        time1 = time.time()                                      #Time
                        gvar.Time_DC =  gvar.Time_DC + (time1-time0)

    #NodeTranslation("Graph_1000.dot")
    #ReactionCleanSub()
    print("========")
    print("Distance calc.:",gvar.Time_DC-gvar.Time_FM)
    print("Fragmentation :",gvar.Time_FM)
    print("Mole Info.    :",gvar.Time_MI)
    print("To smile      :",gvar.Time_TS)
    print("Blockchange   :",gvar.Time_BW)
    return([Assign[4], gvar.GR, gvar.DicStuct,   gvar.DicMoleInfo, gvar.SpeciesCount,BlockRec_first,BlockRec_last])
    #         ^           ^         ^               ^                    ^                    ^            ^
    #     CPU Number    Graph    structureRec     Moleinfomation     NumberOfSpeciesPertime  FirstBlk    LastBlk

# Build Parallel Pool
def parallelPool():
    p = Pool(gvar.ncores)
    G_temp = nx.DiGraph() 
    if(gvar.trajTyp =="xyz"):
        print(gvar.ProcessArrangement)
        AssiT = p.map(parallelReadXYZ,gvar.ProcessArrangement)
        AssiT.sort(key = lambda x: x[0])

    if(gvar.trajTyp =="lmp"):
        AssiT = p.map(parallelReadLAMMPS,gvar.ProcessArrangement)
        AssiT.sort(key = lambda x: x[0])

    # Combine the return structures.  
    for itm in AssiT:
       #print(itm[0],itm[1].number_of_nodes())
       #print(itm[0],gvar.StepRecPerCore[itm[0]])
        # Combine structure dictionary.
        gvar.DicStuct = {**gvar.DicStuct,**itm[2]}
        printUnknowStruc()
      # print("======= Core: ",itm[0],"========")
      # for step in itm[2]:
      #     print(step,itm[2][step][3])

        # Combine moleinfomation dictionary.
        gvar.DicMoleInfo = {**gvar.DicMoleInfo,**itm[3]}

        # Combine species count.
        print(len(gvar.SpeciesCount),len(itm[4]))
        gvar.SpeciesCount = gvar.SpeciesCount + itm[4]

        # Combine reaction event graph.
        gvar.GR.add_nodes_from(list(itm[1].nodes(data=True)) )
        gvar.GR.add_edges_from(list(itm[1].edges(data=True)) )
    print(gvar.GR.number_of_edges())
    for i in range(1,len(AssiT)):
        compare2Step(AssiT[i-1][6],AssiT[i][5],gvar.StepRecPerCore[i-1][1])

   #print(len(gvar.DicStuct),"lineDicStruct")
   #for sml in  gvar.DicStuct:
   #    print(sml,gvar.DicStuct[sml][3])
    NodeTranslation("Graph_1000.dot")
    ReactionClean()


    print(len(gvar.SpeciesCount),gvar.GR.number_of_edges())
    time_1 = time.time()
    write_dot(gvar.GR,"G2cores.dot")
    time_2 = time.time()
    printUnknowStruc()
    time_3 = time.time()
    NodeTranslation("Graph2cores_5000.dot")
    time_4 = time.time()
    gvar.DicReactionRec = generateReactions(gvar.GR)
    time_5 = time.time()
    totalBackup()
    time_6 = time.time()
    PrintMetaDataRaw()
    time_7 = time.time()

    print("savedot:    ", time_2-time_1)
    print("printunknow:", time_3-time_2)
    print("NodeSimpl:  ", time_4-time_3)
    print("ReactionGen ", time_5-time_4)
    print("Backup:     ", time_6-time_5)
    print("PrintRaw:   ", time_7-time_6)
    print("NC         :",gvar.Time_ER)
    print("ER         :",gvar.Time_NC)

    
