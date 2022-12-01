# 2gb3 fnavgl
from __future__ import (
    absolute_import, division, print_function, unicode_literals,
)
# Third-party libraries
import numpy as np
import gc
from math import sqrt
from scipy.spatial.distance import squareform, pdist,cdist
from numpy.linalg import norm
from union_find import *
import global_var as gvar
from tool_contract import *
from tool_removeOscill import *
#import resource

#==============================
def BuildMaskForXYZ(fname):
    bnd_cri = gvar.bnd_cri
    with  open(fname,'r') as f:
        #gvar.pbcXYZ = PBCRec[istep-1]
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


            if(gvar.cataSelectFlag):
                if(line[0] == gvar.cataLabel): gvar.CataAtom.append(i)
        # Initialize global Mask Matrix
        Element = [row[0] for row in gvar.atomList] 
        gvar.GlobalMaskMat = np.array([*map(gvar.radii_dict.get, Element)],dtype=np.float16)
        gvar.GlobalMaskMat = np.tile(gvar.GlobalMaskMat,(len(Element),1))
        gvar.GlobalMaskMat = (gvar.GlobalMaskMat+gvar.GlobalMaskMat.T)*bnd_cri
        np.fill_diagonal(gvar.GlobalMaskMat,0)
        #gvar.GlobalMaskMat = squareform(radii_array)
        #del radii_array
        return

def BuildMaskForLAMMPS(fname):
    bnd_cri = gvar.bnd_cri
    readFlag = 0 # 1. number of atoms; 2. Box bounds 3. Coord
    istep = 0
    with  open(fname,'r') as f:
        while True:
            line = f.readline()
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
                    if(gvar.cataSelectFlag and CoorRec[ielement] == gvar.cataLabel):
                        gvar.CataAtom.append(int(CoorRec[iatom])-1)
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
                    gvar.GlobalMaskMat = np.array([*map(gvar.radii_dict.get, Element)])
                    gvar.GlobalMaskMat = np.tile(gvar.GlobalMaskMat,(len(Element),1))
                    gvar.GlobalMaskMat = (gvar.GlobalMaskMat+gvar.GlobalMaskMat.T)*bnd_cri
                    np.fill_diagonal(gvar.GlobalMaskMat,0)
                    #gvar.GlobalMaskMat = squareform(radii_array)
                    #del radii_array
                    return

#===============================
# Calculate distant between atom 'atm'
# and a list of atoms 'A'
# Input: 
#     atm : atom, [x,y,z]
#     A:  distance list
#            [[X Y Z],
#             [X Y Z],
#             [X Y Z],
#               ...
#             [X Y Z]]
#
#     p: 3D length of PBC box
# Return: distance array
#              
#===============================
def calc_ListAtomPD(atm,A,p=(10,10,10)):
    out = np.empty((3, A.shape[0])) 
    idx = 0

    for o, i in zip(out, A.T):
        out[idx] = cdist([[atm[idx]]],i[:, None], 'cityblock')
        idx += 1
    out[0][out[0] > p[0]/2] -= p[0]
    out[1][out[1] > p[1]/2] -= p[1]
    out[2][out[2] > p[2]/2] -= p[2]
    out16 = np.float16(out)
    del(out)
    return (norm(out16, axis=0))

#===============================
# Pair distance for a list of points
#
# Input: A: distance list
#            [[X Y Z],
#             [X Y Z],
#             [X Y Z],
#               ...
#             [X Y Z]]
#
#        p: 3D length of PBC box
# Return: distanceMatrix
#===============================

def calc_ListPD(A,p=(10,10,10)):
    out = np.empty((3, A.shape[0]*(A.shape[0]-1)//2)) 
    for o, i in zip(out, A.T):
        pdist(i[:, None], 'cityblock', out=o)
    out[0][out[0] > p[0]/2] -= p[0]
    out[1][out[1] > p[1]/2] -= p[1]
    out[2][out[2] > p[2]/2] -= p[2]
    out16 = np.float16(out)
    del(out)
    gc.collect()
    return squareform(norm(out16, axis=0))

def calc_ListPD_nPBC(A):
    dist = pdist(A)
    if(len(dist)>0): 
        MaxD = np.max(dist)
    else:
        MaxD = 2
    return squareform(dist),MaxD
    
#===============================
# Build search neighbor based on Atoms
#
# Input: atomsList, cri: criteria, pbcxyz
# Return: 1. atomList with neighbor
#         2. Groups [[grp1],[grp2],[grp3],....]
#===============================

def buildNeigh_AtomicBased(cri):
    import time
    pbcxyz = gvar.pbcXYZ
    atomList = gvar.atomList
    dim = len(atomList)

    atmC = np.array([row[1] for row in atomList],dtype=np.float16)
    p_1 = time.time()


    gvar.GlobalDistMat= calc_ListPD(atmC,
    (pbcxyz[1][0]-pbcxyz[0][0],
     pbcxyz[1][1]-pbcxyz[0][1],
     pbcxyz[1][2]-pbcxyz[0][2]))
    
    p_2 = time.time()
    # Mask for neibor
    mask = np.full((dim,dim), cri, dtype=np.float16)
    #mask = np.float16(np.ones(int( dim*(dim-1)/2) )*cri)
    resultM = gvar.GlobalDistMat < mask
    del mask

    # Mask for fragment
    LinkMat =   gvar.GlobalDistMat < gvar.GlobalMaskMat

    '''
    for i in range(dim-1):
        I =  (i+2)*(i+1) 
        for j in range(i+1,dim):
            idx = dim*i+j - I//2
            if(resultM[idx]):
                gvar.atomList[i][2].append(j)
    '''
    for i in range(dim):
        LinkMat[i][i] = False
        for j in range(i):
            if(resultM[i,j]):
                gvar.atomList[i][2].append(j)

    uf = groupSplit(LinkMat)
    MolRec = uf.components()
    nMol = len(MolRec)
    del LinkMat
    #del DistRes
    del resultM
    gc.collect()
    p_3 = time.time()
    gvar.Time_FM =gvar.Time_FM + (p_3-p_2)

    return MolRec,nMol

def buildDistMart(cri):
    import time
    pbcxyz = gvar.pbcXYZ
    atomList = gvar.atomList
    p_1 = time.time()
    idx = 0
    #gvar.GlobalDistMat = DistRes
    #GDistMat = gvar.GlobalDistMat 
    for iatm in atomList:
        atmList = iatm[2] # get neighbor of iatm
        atmCT = np.array([atomList[intam][1] for intam in atmList])
        o = calc_ListAtomPD(np.array(iatm[1]),atmCT,
        (pbcxyz[1][0]-pbcxyz[0][0],
         pbcxyz[1][1]-pbcxyz[0][1],
         pbcxyz[1][2]-pbcxyz[0][2]))
        idx2 = 0
        for j in atmList:
            gvar.GlobalDistMat[idx][j] = gvar.GlobalDistMat[j][idx]= o[idx2] 
            idx2 += 1
        idx += 1
    p_2 = time.time()

    LinkMat = gvar.GlobalDistMat < gvar.GlobalMaskMat 
    #LinkMat = gvar.GlobalDistMat < squareform(gvar.GlobalMaskMat)
    uf = groupSplit(LinkMat)
    del LinkMat
    gc.collect()
    MolRec = uf.components()
    nMol = len(MolRec)
    p_3 = time.time()
    gvar.Time_FM =gvar.Time_FM + (p_3-p_2)

    return MolRec,nMol

#===============================
# Update block information list with PBC
#
# Input: BlockList, AtomList,pbcxyz
# Return: 1. atomList with neighbor
#         2. Groups [[grp1],[grp2],[grp3],....]
#===============================

def MolHash(connectMat,element):
    import networkx as nx 
    from collections import Counter
    G = nx.Graph()
    EleCount = Counter(element)
    for i in range(len(connectMat)): 
        G.add_node(i,label=element[i])
        for j in range(i): 
            if connectMat[i][j]: 
                G.add_edge(i,j)
    strR = ""
    if "C" in EleCount.keys():
        if(EleCount["C"]>1):
            strR += "C"+str(EleCount["C"])
        else:
            strR += "C"

    if "H" in EleCount.keys():
        if(EleCount["H"]>1):
            strR += "H"+str(EleCount["H"])
        else:
            strR += "H"

    if "O" in EleCount.keys():
        if(EleCount["O"]>1):
            strR += "O"+str(EleCount["O"])
        else:
            strR += "O"

    if "N" in EleCount.keys():
        if(EleCount["N"]>1):
            strR += "N"+str(EleCount["N"])
        else:
            strR += "N"

    if "F" in EleCount.keys():
        if(EleCount["F"]>1):
            strR += "F"+str(EleCount["F"])
        else:
            strR += "F"
    #print(nx.weisfeiler_lehman_graph_hash(G,node_attr='label'),strR)
    return nx.weisfeiler_lehman_graph_hash(G,node_attr='label'),strR
    

def MolCenter(atomListSub):
    center=[0,0,0]
    nAtom = len(atomListSub)
    for atom in atomListSub:
        center[0]=center[0]+atom[1][0]/nAtom
        center[1]=center[1]+atom[1][1]/nAtom
        center[2]=center[2]+atom[1][2]/nAtom
    return center


def getBlkInfoPBC(atomList, Brec, cri, pbcxyz):
    import time
    pbcxyz = gvar.pbcXYZ
    atmC = np.array([row[1] for row in atomList])
    Element = np.array([row[0] for row in atomList])
    p=(
    pbcxyz[1][0]-pbcxyz[0][0],
    pbcxyz[1][1]-pbcxyz[0][1],
    pbcxyz[1][2]-pbcxyz[0][2]) 
    # Remove PBC
    for i in range(3):
        for j in range(1,len(atmC)):
            m = atmC[j,i]-atmC[0,i]
            dis = abs(m)
            if(dis > p[i]/2 ): atmC[j,i] -= p[i] * np.sign(m)

    distMat,MaxD = calc_ListPD_nPBC(atmC)
    '''
    distMat = np.zeros((len(Brec[1]),len(Brec[1]) ))
    a1 = 0
    for i in Brec[1]:
        a2 = 0
        for j in Brec[1]:
            distMat[a1,a2] = gvar.GlobalDistMat[i,j]
            a2 += 1
        a1 += 1

    MaxD = np.max(distMat)
    '''
    #Build Mask matrix
    radii_list = np.array([*map(gvar.radii_dict.get, Element)])
    radii_list = np.tile(radii_list,(len(Element),1))
    MaskMat = (radii_list+radii_list.T)*cri
    np.fill_diagonal(MaskMat,-1)
    # LinkMat
    LKmat = distMat<MaskMat*1
    # Transform Structure to SMILES
    #[SMILES, title] = xyzfileToSMILE(Element,list(atmC))
    Rhash,formula = MolHash(LKmat,Element)
    Rhash = "H"+Rhash
    SMILES=Rhash
    # Add unknow SMILES to dictionary
    if (Rhash in gvar.DicStuct.keys()):
        pass
    else:
        gvar.DicStuct[Rhash] = [Element,list(atmC),False,"",formula]
    return SMILES,MaxD


def BlockInfoUpdatePBC(cri):
    from collections import Counter
    import hashlib
    pbcxyz = gvar.pbcXYZ
    atomList = gvar.atomList
    BlockList = gvar.blockList
    iblk = 0
    HASH_count = []
    for rec in gvar.blockList:
        atomListSub = []
        rec[1].sort()
        for i in rec[1]:
            atomListSub.append(atomList[i])
        hashD,MaxD = getBlkInfoPBC(atomListSub,rec,cri,pbcxyz)
        #mat = np.array(buildLinkMatSubPBC(atomListSub,1.50,pbcxyz))
        rec[2] = atomListSub 
        rec[3] = hashD

        HASH_count.append(hashD)
        rec[4] = MolCenter(atomListSub)
        rec[5] = MaxD/2
        # Build fragment HASH ID.
        AtmLabelStr ="".join([str(i) for i in rec[1]])
        rec[6] = hashlib.sha1((AtmLabelStr+hashD).encode('utf-8')).hexdigest()[0:20]
        # Atom fragment information to DicMoleInfo
        Label = "S"+rec[6]
        if(Label not in gvar.DicMoleInfo):
            gvar.DicMoleInfo[Label] = rec[1]

    # Count species.
    gvar.SpeciesCount.append(Counter(HASH_count))

def BlockNeighborUpdate(cri):
    import time
    p_1 = time.time()
    pbcxyz = gvar.pbcXYZ
    blockList = gvar.blockList
    atomList = gvar.atomList
    dim = len(blockList)
    MoC = np.array([row[4] for row in blockList])
    # Build MaskMat
    radii_block = \
       np.array([row[5] for row in blockList],dtype=np.float16)
    radii_block = np.tile(radii_block,(dim,1))
    radii_Maskmat = radii_block + radii_block.T + cri
    np.fill_diagonal(radii_Maskmat,-1)

    # calculte distance of center of mass 
    DistRes = calc_ListPD(MoC,
    (pbcxyz[1][0]-pbcxyz[0][0],
     pbcxyz[1][1]-pbcxyz[0][1],
     pbcxyz[1][2]-pbcxyz[0][2]))
    p_2 = time.time()
    resultM = DistRes  < radii_Maskmat

    # Add neighbor to list
    for i in range(dim):
        neiListT = []
        for j in range(i):
            if(resultM[i,j]):
                neiListT = neiListT + blockList[j][1]

        for iatm in blockList[i][1]:
            nei_local = [k for k in blockList[i][1] if k > iatm ]
            atomList[iatm][2] = neiListT + nei_local
    p_3 = time.time()

#+++++++++++++++++++++++++++++
#
# Build reaction links 
#
#+++++++++++++++++++++++++++++

def compare2Step(blk1,blk2,istep):
    import networkx as nx
    from networkx.drawing.nx_pydot import write_dot
    G_t = nx.DiGraph()
    GR = gvar.GR
    hashindex_blk1 = [x[6] for x in blk1]
    hashindex_blk2 = [x[6] for x in blk2]
    set_blk1 = set(hashindex_blk1)
    set_blk2 = set(hashindex_blk2)

    FromBlks = {}
    ToBlks =   {}

    commPart = set_blk1.intersection(set_blk2) 
    FromBlks_idx = set_blk1-commPart
    ToBlks_idx   = set_blk2-commPart

    FromAtomListSet = []
    FromIDX = []
    for ihash in FromBlks_idx:
        idx = hashindex_blk1.index(ihash)
        FromBlks[ihash] = blk1[idx]
        FromAtomListSet.append(set(blk1[idx][1]))
        FromIDX.append(idx)

    ToAtomListSet = []
    ToIDX = []
    for ihash in ToBlks_idx:
        idx = hashindex_blk2.index(ihash)
        ToBlks[ihash]   = blk2[idx]
        ToAtomListSet.append(set(blk2[idx][1]))
        ToIDX.append(idx)
   
    nArrow = 0
    for i in range(len(FromAtomListSet)):
        for j in range(len(ToAtomListSet)):
            overlapSet = FromAtomListSet[i].intersection(ToAtomListSet[j])
            if(len(overlapSet) != 0):
                nArrow += 1
                # For isomerization. ## Need to be optimized. 
                if(len(FromAtomListSet[i]) == len(ToAtomListSet[j])):
                    G_t.add_edge("R"+str(FromIDX[i]), "P"+str(ToIDX[j]))
                # Other Cases.
                else:
                    G_t.add_edge("R"+str(FromIDX[i]), "P"+str(ToIDX[j]))

   # Add curent step reaction to  graph
   #print(nArrow,len(list(nx.weakly_connected_components(G_t))))
    subGrapRec = list(nx.weakly_connected_components(G_t)) 
    for link in subGrapRec:
        RList = []
        PList = []
        NR = 0
        NP = 0
        for itm in link:
            if itm[0] == 'R':
                NR += 1
                RList.append(int(itm[1:]))
            if itm[0] == 'P':
                NP += 1
                PList.append(int(itm[1:]))

        # Multi to Multi
        if(NR >  1 and NP >  1):
            num = 0

            for RC in RList:
                if ("CC_" in str(blk1[RC][6]) or "NC_" in str(blk1[RC][6])):
                    num += int(blk1[RC][6][3:],16)
                else:
                    num += int(blk1[RC][6],16)

            for PD in PList:
                if ("CC_" in str(blk2[PD][6]) or "NC_" in str(blk2[PD][6])):
                    num -= int(blk2[PD][6][3:],16)
                else:
                    num -= int(blk2[PD][6],16)

            InterIndex = "int-"+hex(num)
            InterLabel = "int-"+hex(num)[0:5]
            GR.add_node(InterIndex,label=InterLabel,hashD=InterLabel)

            for RC in RList:
                HashR = "S"+blk1[RC][6]
                GR.add_node(HashR,hashD=blk1[RC][3],
                                  label=gvar.DicStuct[blk1[RC][3]][4])
#               GR.add_node(HashR,label=gvar.DicStuct[blk1[RC][3]][4])
                GR.add_edge(HashR,InterIndex,label = str(istep),color="red")

            for PD in PList:
                HashP = "S"+blk2[PD][6]
                GR.add_node(HashP,hashD=blk2[PD][3],
                                  label=gvar.DicStuct[blk2[PD][3]][4])
#               GR.add_node(HashP,label=gvar.DicStuct[blk2[PD][3]][4])
                GR.add_edge(InterIndex,HashP,label = str(istep),color="blue")

        # Combination reaction
        if(NR >  1 and NP == 1):
            HashP = "S"+blk2[PList[0]][6]
            GR.add_node(HashP,hashD=blk2[PList[0]][3],
                              label=gvar.DicStuct[blk2[PList[0]][3]][4])

            for RC in RList:
                HashR = "S"+blk1[RC][6]
                GR.add_node(HashR,hashD=blk1[RC][3],
                                  label=gvar.DicStuct[blk1[RC][3]][4])
                GR.add_edge(HashR,HashP,label = str(istep),color="red")

        # Split reaction
        if(NR == 1 and NP  > 1):
            HashR = "S"+blk1[RList[0]][6]
            GR.add_node(HashR,hashD=blk1[RList[0]][3],
                              label=gvar.DicStuct[blk1[RList[0]][3]][4])
#           GR.add_node(HashR,label=gvar.DicStuct[blk1[RList[0]][3]][4])

            for PD in PList:
                HashP = "S"+blk2[PD][6]
                GR.add_node(HashP,hashD=blk2[PD][3],
                                  label=gvar.DicStuct[blk2[PD][3]][4])
#               GR.add_node(HashP,label=gvar.DicStuct[blk2[PD][3]][4])
                GR.add_edge(HashR,HashP,label = str(istep),color="blue")

        # Isomerization reaction
        if(NR == 1 and NP == 1):
            HashR = "S"+blk1[RList[0]][6]
            HashP = "S"+blk2[PList[0]][6]
            GR.add_node(HashR,hashD=blk1[RList[0]][3],
                              label=gvar.DicStuct[blk1[RList[0]][3]][4])

            GR.add_node(HashP,hashD=blk2[PList[0]][3],
                              label=gvar.DicStuct[blk2[PList[0]][3]][4])
#           GR.add_node(HashR,label=gvar.DicStuct[blk1[RList[0]][3]][4])
#           GR.add_node(HashP,label=gvar.DicStuct[blk2[PList[0]][3]][4])

            GR.add_edge(HashR,HashP,label = str(istep),color="grey")
#   if(istep == 101):
#       NodeTranslation(GR,Gname = "Graph_1000.dot")

def NodeTranslation(Gname = "Graph_1000.dot"):
    GR = gvar.GR
    for node in GR.nodes(data=True):
        hashD = node[1]['hashD']
        if("int-" not in hashD):
            node[1]["SMILE"] = gvar.DicStuct[hashD][3]
        else:
            node[1]["SMILE"] = "int"
    write_dot(GR,"UnGraph_1000.dot")
    #ReactionClean()
    GR = gvar.GR
    for node in GR.nodes(data=True):
        if('contraction' in node[1]):
            del node[1]['contraction']
    #write_dot(GR,Gname)



def ReactionClean():
    import networkx as nx
    from copy import deepcopy

    G = gvar.GR 
    G_t = deepcopy(G)
    for edge in  G_t.edges(data=True,keys=True):
        for key in G_t[edge[0]][edge[1]]:
            if(isinstance(G_t[edge[0]][edge[1]][key]['label'], str)):
                G_t[edge[0]][edge[1]][key]['label'] = int(G_t[edge[0]][edge[1]][key]['label'].strip('\"'))

    N_nodes = G.number_of_nodes()
    N_edges = G.number_of_edges()
    print("total nodes in original graph: ",N_nodes)
    print("total edges in original graph: ",N_edges)
    N_nodesp = -1
    while(N_nodes - N_nodesp != 0):
        N_nodesp = N_nodes
        p_0 = time.time()
        G_t = remove_useless_trans(G_t,gvar.StableMolLag)
        p_1 = time.time()
        gvar.Time_ER = gvar.Time_ER + (p_1-p_0)
        N_nodes = G_t.number_of_nodes()
        N_edges = G_t.number_of_edges()

#   N_nodesp = -1
#   while(N_nodes - N_nodesp != 0):
        p_0 = time.time()
        G_t = contract_nodes(G_t)
        N_nodesp = N_nodes
        N_nodes = G_t.number_of_nodes()
        N_edges = G_t.number_of_edges()
        G_t = nodeShortRemove(G_t,gvar.StableMolLag)
        p_1 = time.time()
        gvar.Time_NC = gvar.Time_NC + (p_1-p_0)
    #   print("total nodes after GBR: ",N_nodes)
    #   print("total edges after GBR: ",N_edges)
    print("total nodes after removal: ",N_nodes)
    print("total edges after removal: ",N_edges)

    gvar.GR = G_t

def ReactionCleanSub():

    #====================================
    #
    # Clean local reaction network
    #
    #====================================

    import networkx as nx
    from copy import deepcopy

    G = gvar.GR 
    G_t = deepcopy(G)
    for edge in  G_t.edges(data=True,keys=True):
        for key in G_t[edge[0]][edge[1]]:
            if(isinstance(G_t[edge[0]][edge[1]][key]['label'], str)):
                G_t[edge[0]][edge[1]][key]['label'] = int(G_t[edge[0]][edge[1]][key]['label'].strip('\"'))

    N_nodes = G.number_of_nodes()
    N_edges = G.number_of_edges()
    print("total nodes in original graph: ",N_nodes)
    print("total edges in original graph: ",N_edges)
    N_nodesp = -1
    while(N_nodes - N_nodesp != 0):
        N_nodesp = N_nodes
        p_0 = time.time()
        G_t = remove_useless_trans(G_t,int(gvar.StableMolLag/2))
        p_1 = time.time()
        gvar.Time_ER = gvar.Time_ER + (p_1-p_0)
        N_nodes = G_t.number_of_nodes()
        N_edges = G_t.number_of_edges()
        p_0 = time.time()
        G_t = contract_nodes(G_t)
        p_1 = time.time()
        gvar.Time_NC = gvar.Time_NC + (p_1-p_0)
        N_nodesp = N_nodes
        N_nodes = G_t.number_of_nodes()
        N_edges = G_t.number_of_edges()
        G_t = nodeShortRemove(G_t,gvar.StableMolLag)
    print("total nodes after removal: ",N_nodes)
    print("total edges after removal: ",N_edges)


    gvar.GR = G_t

def ReactionCleanGreyTrans():

    #====================================
    #
    # Clean local reaction network
    #
    #====================================

    import networkx as nx
    from copy import deepcopy

    G = gvar.GR 
    G_t = deepcopy(G)
    for edge in  G_t.edges(data=True,keys=True):
        for key in G_t[edge[0]][edge[1]]:
            if(isinstance(G_t[edge[0]][edge[1]][key]['label'], str)):
                G_t[edge[0]][edge[1]][key]['label'] = int(G_t[edge[0]][edge[1]][key]['label'].strip('\"'))

    N_nodes = G.number_of_nodes()
    N_edges = G.number_of_edges()
    print("total nodes in original graph: ",N_nodes)
    print("total edges in original graph: ",N_edges)
    N_nodesp = -1
    while(N_nodes - N_nodesp != 0):
        G_t = contract_nodes_Onlygrey(G_t)
        N_nodesp = N_nodes
        N_nodes = G_t.number_of_nodes()
        N_edges = G_t.number_of_edges()
        #G_t = nodeShortRemove(G_t,gvar.StableMolLag)
        print("total nodes after GreyRemove: ",N_nodes)
        print("total edges after GreyRemove: ",N_edges)

    gvar.GR = G_t

def printUnknowStruc():
    import os
    from trans_smile import xyzfileToSMILE
    if not os.path.exists('specRec'):
            os.umask(0)
            os.makedirs('specRec',mode=0o777)
    count = 0
    for key in gvar.DicStuct:
        if(gvar.DicStuct[key][2]):
            continue
        gvar.DicStuct[key][2] = True
        count += 1
       #FileName = "S_"+ gvar.DicStuct[key][4] +"_"+ str(key)[0:7]+ ".xyz"
       #f = open("./specRec/"+FileName,"w")
        [SMILES, title] = \
        xyzfileToSMILE(gvar.DicStuct[key][0],gvar.DicStuct[key][1])
        if("N([O])[O]" in SMILES and not SMILES.startswith("N([O])[O]") and ".N([O])[O]" not in SMILES ): 
            SMILES = SMILES.replace("N([O])[O]", "N(=O)=O")
        if("[O]N([O])" in SMILES and not SMILES.startswith("[O]N([O])") and ".[O]N([O])" not in SMILES ):
            SMILES = SMILES.replace("[O]N([O])", "O=N(=O)")
        gvar.DicStuct[key][3] = SMILES

   # print(count," structrues are transform to SMILES")


##
#   Read box from LAMMPS files.
#
##
def readPBCfromLAMMPS(fname):
    istep = 0
    PBCrec = []
    with  open(fname,'r') as f:
        while True:
            line = f.readline()
            if not line:break
            if("ITEM: BOX BOUNDS" in line):
                line1 = f.readline().split()
                line1 =[float(line1[0]),float(line1[1])]

                line2 = f.readline().split()
                line2 =[float(line2[0]),float(line2[1])]

                line3 = f.readline().split()
                line3 =[float(line3[0]),float(line3[1])]
                PBCrec.append([(line1[0],line2[0],line3[0]),\
                               (line1[1],line2[1],line3[1])])
            else:
                pass
    return PBCrec


#########
#
# Backup
#
##########
def totalBackup():
    import pickle
    pickle.dump( gvar.DicReactionRec, open( "DicReactRec.p", "wb" ))
    pickle.dump( gvar.DicStuct,       open( "DicStruct.p", "wb" ))
    pickle.dump( gvar.DicMoleInfo,    open( "DicMoleinfo.p", "wb" ))
    pickle.dump( gvar.SpeciesCount,   open( "SpecieCount.p", "wb" ))




if __name__ == "__main__":
   #atm = np.array( [1,1,1] )
   #A   = np.array( [[1,0,1],[9,9,0],[1,9,9]])
   #o = calc_ListAtomPD(atm,A,p=(10,10,10))
   #print(A)
   #print(o)
   #readPBCfromLAMMPS("reax.trj")
   connectMat = [[True,False],[False,True]]
   element = ['C','H']
   r = MolHash(connectMat,element)
   print(r)
   connectMat = [[True,False],[False,True]]
   element = ['H','C']
   r = MolHash(connectMat,element)
   print(r)
   pass

