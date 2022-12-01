from __future__ import (
    absolute_import, division, print_function, unicode_literals,
)
# Third-party libraries
import numpy as np
from math import sqrt
from scipy.spatial.distance import squareform, pdist,cdist
from tool import *
from union_find import * 

# Tools for solving catalysis

def buildNeigh_AtomicBased_cata(cri):
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
    print(p_2-p_1)
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
    print(nMol)
    del LinkMat
    #del DistRes
    del resultM
    gc.collect()

    return MolRec,nMol

def buildDistMart_cata(cri):
    import time
    pbcxyz = gvar.pbcXYZ
    atomList = gvar.atomList
    p_1 = time.time()
    idx = 0
    #gvar.GlobalDistMat = DistRes
    #GDistMat = gvar.GlobalDistMat 
    for iatm in atomList:
        atmList = iatm[2]+gvar.CataAtom # get neighbor of iatm
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
    print(p_2-p_1)


    return MolRec,nMol

def getBlkInfoPBC_cata(atomList, Brec, cri, pbcxyz,cataLabel):
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
    Rhash = cataLabel +"H"+Rhash
    SMILES=Rhash
    # Add unknow SMILES to dictionary
    if (Rhash in gvar.DicStuct.keys()):
        pass
    else:
        gvar.DicStuct[Rhash] = [Element,list(atmC),False,"",formula]
    return SMILES,MaxD

def BlockInfoUpdatePBC_cata(blockList_o,cri):
    from collections import Counter
    import hashlib
    pbcxyz = gvar.pbcXYZ
    atomList = gvar.atomList
    BlockList = blockList_o
    iblk = 0
    HASH_count = []
    TotalBlockTemp=[]
    # First Loop to catch data.
    for rec in BlockList:
        CataAtoms = list(set(gvar.CataAtom)&set(rec[1]))
        # Large block with catalyst.
        if( len(CataAtoms)>0):
            RegularAtoms = list(set(rec[1])-set(CataAtoms))
            if(len(RegularAtoms) == 0): continue
            MaskMatT = gvar.GlobalMaskMat[np.ix_(RegularAtoms,RegularAtoms)]
            DistMatT = gvar.GlobalDistMat[np.ix_(RegularAtoms,RegularAtoms)]
            LinkMatT = DistMatT < MaskMatT
            uf = groupSplit(LinkMatT)
            MolRec = uf.components()
            for irec in MolRec:
                localIdx = []
                for idx in irec:
                    localIdx.append(RegularAtoms[idx])
                TotalBlockTemp.append([iblk,localIdx,"C"])
                iblk += 1
            continue
        TotalBlockTemp.append([iblk,rec[1],"NC"])
        iblk += 1
    nMol = len(TotalBlockTemp)
    gvar.blockList = []    
    for rec in TotalBlockTemp:
        recTemp = [0,0,0,0,0,0,0,0] 
        atomListSub = []
        rec[1].sort()
        recTemp[0] = rec[0]
        recTemp[1] = rec[1]
        recTemp[7] = rec[2]

        for i in rec[1]:
            atomListSub.append(atomList[i])

        if(recTemp[-1] == 'C'):
            catalabel = "CC_"
        else:
            catalabel = "NC_"
            
        hashD,MaxD = getBlkInfoPBC_cata(atomListSub,rec,cri,pbcxyz,catalabel)
        recTemp[2] = atomListSub

        recTemp[3] = hashD
        HASH_count.append(hashD)
        recTemp[4] = MolCenter(atomListSub)
        recTemp[5] = MaxD/2

        # Build fragment HASH ID.
        AtmLabelStr ="".join([str(i) for i in rec[1]])
        recTemp[6] = hashlib.sha1((AtmLabelStr+hashD).encode('utf-8')).hexdigest()
        gvar.blockList.append(recTemp)

        if(recTemp[-1] == 'C'):
            recTemp[6] = "CC_"+recTemp[6]

        else:
            recTemp[6] = "NC_"+recTemp[6]

        # Atom fragment information to DicMoleInfo
        Label = "S"+recTemp[6]
        if(Label not in gvar.DicMoleInfo):
            gvar.DicMoleInfo[Label] = rec[1]
    
    # Count species.
    gvar.SpeciesCount.append(Counter(HASH_count))

def BlockNeighborUpdate_cata(cri):
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


def printUnknowStruc_cata():
    import os
    from trans_smile import xyzfileToSMILE
    if not os.path.exists('specRec'):
            os.umask(0)
            os.makedirs('specRec',mode=0o777)
    count = 0
    for key in gvar.DicStuct:
        print(key,gvar.DicStuct[key][3])
        if(gvar.DicStuct[key][2]):
            continue
        gvar.DicStuct[key][2] = True
        count += 1
       #FileName = "S_"+ gvar.DicStuct[key][4] +"_"+ str(key)[0:7]+ ".xyz"
       #f = open("./specRec/"+FileName,"w")
        [SMILES, title] = \
        xyzfileToSMILE(gvar.DicStuct[key][0],gvar.DicStuct[key][1])
        
        gvar.DicStuct[key][3] =key[0:3]+SMILES
    print(count," structrues are transform to SMILES")

if __name__ == "__main__":
    print("catalyst")
    pass
