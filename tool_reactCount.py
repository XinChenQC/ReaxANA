import networkx as nx

import pydot
from copy import deepcopy
from tool_removeOscill import *
from tool_contract import *
from networkx.drawing.nx_pydot import write_dot

def Sort(sub_li,iplace): 
    # reverse = None (Sorts in Ascending order) 
    # key is set to sort using second element of  
    # sublist lambda has been used 
    sub_li.sort(key = lambda x: x[iplace]) 
    return sub_li 

def GreyReactAbs(G,removeList):    
    # Detect Reactions
    GreyList = []
    
    EdgeList  = list(G.edges(data=True))
    for edge in EdgeList:
        if(edge[2]['color'] == 'grey'):
            NameFrom = G.nodes[edge[0]]['SMILE']
            NameTo = G.nodes[edge[1]]['SMILE']
            removeList.append((edge[0],edge[1],edge[2]['label']))
              
            if(len(NameFrom) == 2 or len(NameTo) == 2 or NameFrom == NameTo):   
                continue
            GreyList.append(["G",(edge[0],edge[1]),(NameFrom,NameTo),edge[2]['label']])

    return GreyList,removeList

def getReactionID(Rlist,Plist,Catalist):
    import hashlib
    Rlist.sort()
    Plist.sort()
    Catalist.sort()
    HashStr = ""
    HashStr = "R " + HashStr+ str(len(Rlist))  #
    for itm in Rlist: 
        HashStr = HashStr+itm+"-"              #              
        
    HashStr = HashStr + " P "+ str(len(Plist)) #
    for itm in Plist:
        HashStr = HashStr+itm+"-"              #
        
    HashStr = HashStr + " cata "
    for itm in Catalist:
        HashStr = HashStr+itm
    return(hashlib.sha1(HashStr.encode('utf-8')).hexdigest()) 

def buildBasicInfo(G):
    RawNodeList = G.nodes(data=True)
    rec = []
    for node in RawNodeList:
        # Nodeinfo=[nodeID,nodelabel,[Inedges],[OutEdges]]
        # Get Label
        nodeInfo = [node[0],node[1]['label']]
        
        # GetInEdges:
        InEdge = []
        edgesRawIn = G.in_edges(node[0],data=True)
        for edge in edgesRawIn:
            InEdge.append((edge[0],edge[2]['color'],int(edge[2]['label'])))
        InEdge = Sort(InEdge,2)
        #print(InEdge)
        nodeInfo.append(InEdge)
        
        # GetOutEdges:
        OutEdge = []
        edgesRawOut = G.out_edges(node[0],data=True)
        for edge in edgesRawOut:
            OutEdge.append((edge[1],edge[2]['color'],int(edge[2]['label'])))
        OutEdge = Sort(OutEdge,2)
        nodeInfo.append(OutEdge)
        # determine the type of each node. 
        #    'center' nodes are the nodes combined by others. Generely transit structure.
        #    'edge' nodes are nodes without splitting and combined from edge nodes.         
        Type='edge'
        for i in range(len(nodeInfo[2]) - 1):
            #print(nodeInfo[0],nodeInfo[2][i][2], nodeInfo[2][i+1][2])
            if(nodeInfo[2][i][2] == nodeInfo[2][i+1][2] and nodeInfo[2][i][1]=='red' and nodeInfo[2][i+1][1]=='red'):
                Type='center'
                break
        for i in range(len(nodeInfo[3]) - 1):
            if(nodeInfo[3][i][2] == nodeInfo[3][i+1][2] and nodeInfo[3][i][1]=='blue' and nodeInfo[3][i+1][1]=='blue'):
                Type='center'
                break
        nodeInfo.append(Type)
        rec.append(nodeInfo)
    return rec


def CombineRedBlue(G,Recs,Tlag):
    # Combine the red and blue together for analysis.
    def SortRecitm(tup):
        tup_t = deepcopy(tup[1:])
        tup_t = (tup[0],) + tuple(sorted(tup_t))
        return tup_t
    
    CenterRecTot = []
    for rec in Recs:
        if(rec[4]=="center"):
            # Combine red
            CenterRec = []
            CenterRec.append(rec[0])
            CombineRed = []
            StepTup = []
            for itm in rec[2]:
                StepTup.append(itm[2],)
            StepTup = list(set(StepTup))
            for stepNum in StepTup:
                Combinetup_temp = (stepNum,)
                for itm in rec[2]:
                    if(itm[2] == stepNum):Combinetup_temp = Combinetup_temp + (itm[0],)
                #print(Combinetup_temp)
                Combinetup_temp = SortRecitm(Combinetup_temp)
                CombineRed.append(Combinetup_temp)
            CombineRed = Sort(CombineRed,0)
            CenterRec.append(CombineRed)
            #print(CombineRed,"Red")
                
            # Combine blue
            CombineBlue = []
            StepTup = []
            for itm in rec[3]:
                StepTup.append(itm[2],)
            StepTup = list(set(StepTup))
            for stepNum in StepTup:
                Combinetup_temp = (stepNum,)
                for itm in rec[3]:
                    if(itm[2] == stepNum):Combinetup_temp = Combinetup_temp + (itm[0],)
                #print(Combinetup_temp)
                Combinetup_temp = SortRecitm(Combinetup_temp)
                CombineBlue.append(Combinetup_temp)
            CombineBlue = Sort(CombineBlue,0)
            CenterRec.append(CombineBlue)
            #print(CombineBlue,"blue")

            nRed  = len(CombineRed)
            nblue = len(CombineBlue)
            if(nRed == nblue and nRed>1): 
                for itag in range(nRed):
                    SplitRec = [rec[0]]
                    SplitRec.append([CombineRed[itag]])
                    SplitRec.append([CombineBlue[itag]])
                    CenterRecTot.append(SplitRec) 
            else:
                CenterRecTot.append(CenterRec)

            #print("=============")
    return(CenterRecTot)

def SimpleNodeshortReactAbs(CombRecs,Tlag,G_t):    
    # Detect Reactions
    removeList = []
    PRreactions = []
    NodermList = []
    for rec in CombRecs:
        #nPair = min(len(rec[1]),len(rec[2]))
        nRed = len(rec[1])
        nblue = len(rec[2])
        if(not (nRed == nblue and nRed ==1)): continue
        ReacL = rec[1][0][1:]
        ProdL = rec[2][0][1:]
        
        if(rec[2][0][0] - rec[1][0][0] < Tlag and rec[2][0][0]-rec[1][0][0] >= 0 and len(ReacL)>1 and len(ProdL)>1):
            # Build specied Name Rec:
            NameReacL=[]
            #print(G_t.nodes[str(ReacL[0])]['label'])
            for node in ReacL:
                NameReacL.append(G_t.nodes[str(node)]['SMILE'].strip("\""))
            NameProdL=[]
            for node in ProdL:
                NameProdL.append(G_t.nodes[str(node)]['SMILE'].strip("\""))
            
            
            PRreactions.append([list(ReacL),list(ProdL),rec[1][0][0],rec[2][0][0],"S",rec[0]])
            NodermList.append(rec[0])
            for Reacrec in ReacL:
                removeList.append((Reacrec,rec[0],rec[1][0][0]))
            for Prodrec in ProdL:
                removeList.append((rec[0],Prodrec,rec[2][0][0]))
                
        elif (len(ReacL)>1 and len(ProdL)>1 and rec[2][0][0] - rec[1][0][0] > Tlag   and rec[2][0][0]-rec[1][0][0] > 0 ):
            # Build specied Name Rec:
            NameReacL=[]
            #print(G_t.nodes[str(ReacL[0])]['label'])
            for node in ReacL:
                NameReacL.append(G_t.nodes[str(node)]['SMILE'].strip("\""))
            
            NameProdL=[]
            for node in ProdL:
                NameProdL.append(G_t.nodes[str(node)]['SMILE'].strip("\""))
            
            
            PRreactions.append([list(ReacL),list(ProdL),rec[1][0][0],rec[2][0][0],"L",rec[0]])
            for Reacrec in ReacL:
                removeList.append((Reacrec,rec[0],rec[1][0][0]))
            for Prodrec in ProdL:
                removeList.append((rec[0],Prodrec,rec[2][0][0]))            
    SReactions = []
    LReactions = []
    for rec in PRreactions:
        sameEle = tuple(set(rec[0]).intersection(set(rec[1])))
        # For short reactions: A+B -> C+D 
        if(rec[-2] =='S'):
            cata = []
            for i_s in sameEle:
                cata.append(i_s)
                rec[0].remove(i_s)
                rec[1].remove(i_s)
                #print(rec[0].index(i_s),rec[1].index(i_s))
            
            # transform to species name. 
            NameReacL=[]
            #print(G_t.nodes[str(ReacL[0])]['label'])
            for node in tuple(rec[0]):
                NameReacL.append(G_t.nodes[str(node)]['SMILE'].strip("\""))
                
            NameProdL=[]
            for node in tuple(rec[1]):
                NameProdL.append(G_t.nodes[str(node)]['SMILE'].strip("\""))  
                
            NameCataL=[]
            for node in tuple(cata):
                NameCataL.append(G_t.nodes[str(node)]['SMILE'].strip("\""))  
                
            #print(["S",tuple(rec[0]),tuple(rec[1]),tuple(cata),tuple(NameReacL),tuple(NameProdL),tuple(NameCataL), rec[2],rec[3]])
            SReactions.append(["S",tuple(rec[0]),tuple(rec[1]),tuple(cata),tuple(NameReacL),tuple(NameProdL),tuple(NameCataL), rec[2],rec[3]])
        # For long reactions: A+B -> AB; AB -> C+D  
        else:
            NameReacL=[]
            #print(G_t.nodes[str(ReacL[0])]['label'])
            for node in tuple(rec[0]):
                NameReacL.append(G_t.nodes[str(node)]['SMILE'].strip("\""))
                
            NameProdL=[]
            for node in tuple(rec[1]):
                NameProdL.append(G_t.nodes[str(node)]['SMILE'].strip("\""))  
            
            LReactions.append(["LC",tuple(rec[0]),tuple([rec[-1],]),tuple(NameReacL),tuple( [G_t.nodes[str(rec[-1])]['SMILE'].strip("\""),] ), rec[2]])
            LReactions.append(["LS",tuple([rec[-1],]),tuple(rec[1]),tuple( [G_t.nodes[str(rec[-1])]['SMILE'].strip("\""),] ),tuple(NameProdL), rec[3]])
            
            
    return(removeList,NodermList,SReactions,LReactions)

def CombAndSplitReactABS(G_t,recs,removeList):
    CombineReact = []
    SplitReact = []

    for itm in recs:
        # Combine Reaction:
        CenterName = G_t.nodes[str(itm[0])]['SMILE'].strip("\"")
        if(len(itm[1]) > 0):
        #   Split Reaction
            for pair in itm[1]:
                if (len(pair) > 2):
                    ReacL = pair[1:]
                    for Reacrec in ReacL:
                        removeList.append((Reacrec,itm[0],pair[0]))
                    # Get Smiles 
                    NameReacL=[]
                    #print(G_t.nodes[str(ReacL[0])]['label'])
                    for node in ReacL:
                        NameReacL.append(G_t.nodes[str(node)]['SMILE'].strip("\""))
                    #print(NameReacL," -> ",CenterName,1)
                    CombineReact.append(["C",tuple(ReacL),tuple([itm[0],]),tuple(NameReacL),tuple( [CenterName,] ), pair[0]])
                    #print(["C",tuple(ReacL),tuple([itm[0],]),tuple(NameReacL),tuple( [CenterName,] ), pair[0]])
                    
         # Combine Reaction
        if(len(itm[2]) > 0):       
            for pair in itm[2]:
                if (len(pair) > 2):
                    ProdL = pair[1:]
                    for Prodrec in ProdL:
                        removeList.append((itm[0],Prodrec,pair[0]))
                     # Get Smiles 
                    NameProdL=[]
                    for node in ProdL:
                        NameProdL.append(G_t.nodes[str(node)]['SMILE'].strip("\""))
                    #print(CenterName," -> ",NameProdL,2)
                    SplitReact.append(["S",tuple([itm[0],]),tuple(ProdL),tuple( [CenterName,] ),tuple(NameProdL), pair[0]])
                    #print(["S",tuple([itm[0],]),tuple(ProdL),tuple( [CenterName,] ),tuple(NameProdL), pair[0]])
    return CombineReact, SplitReact, removeList

def RegulateAndClassify(SR,LR,CR,SpR,GreyR):
    #============================
    #  Print fast reactions 
    #    A+B -> C+D ....
    #============================
    RegulateRec = []
    # ["HashLable",[HashR], [HashP], [HashCata.], [SMILES_R], [SMILES_P], [SMILES_cata], (time1,time2)]
    # A+B -> C+D reactions
    for rec in SR:
        tmpRec = [0,0,0,0,0,0,0,0]
        # hashrec
        tmpRec[1] = list(rec[1])
        tmpRec[2] = list(rec[2])
        tmpRec[3] = list(rec[3])
        # SMILES
        tmpRec[4] = list(rec[4])
        tmpRec[5] = list(rec[5])
        tmpRec[6] = list(rec[6])
        #Time = 
        tmpRec[7] = (rec[7],rec[8])
        tmpRec[0] = getReactionID(tmpRec[4],tmpRec[5],tmpRec[6])
        if(len(tmpRec[4]) > 0 and  len(tmpRec[5])>0 ):
            RegulateRec.append(tmpRec)
        
    # A+B -> C and C-> A+B reactions in Long time
    for rec in LR:
        tmpRec = [0,0,0,0,0,0,0,0]
        # hashrec
        tmpRec[1] = list(rec[1])
        tmpRec[2] = list(rec[2])
        tmpRec[3] = []
        # SMILES
        tmpRec[4] = list(rec[3])
        tmpRec[5] = list(rec[4])
        tmpRec[6] = []
        #Time = 
        tmpRec[7] = (rec[5],-1)
        tmpRec[0] = getReactionID(tmpRec[4],tmpRec[5],tmpRec[6])
        RegulateRec.append(tmpRec)
        
    # A+B -> C
    for rec in CR:
        tmpRec = [0,0,0,0,0,0,0,0]
        # hashrec
        tmpRec[1] = list(rec[1])
        tmpRec[2] = list(rec[2])
        tmpRec[3] = []
        # SMILES
        tmpRec[4] = list(rec[3])
        tmpRec[5] = list(rec[4])
        tmpRec[6] = []
        #Time = 
        tmpRec[7] = (rec[5],-1)
        tmpRec[0] = getReactionID(tmpRec[4],tmpRec[5],tmpRec[6])
        RegulateRec.append(tmpRec)
        
    # C -> A+B    
    for rec in SpR:
        tmpRec = [0,0,0,0,0,0,0,0]
        # hashrec
        tmpRec[1] = list(rec[1])
        tmpRec[2] = list(rec[2])
        tmpRec[3] = []
        # SMILES
        tmpRec[4] = list(rec[3])
        tmpRec[5] = list(rec[4])
        tmpRec[6] = []
        #Time = 
        tmpRec[7] = (rec[5],-1)
        tmpRec[0] = getReactionID(tmpRec[4],tmpRec[5],tmpRec[6])
        RegulateRec.append(tmpRec)

    # A -> A'    
    for rec in GreyR:
        tmpRec = [0,0,0,0,0,0,0,0]
        # hashrec
        tmpRec[1] = [rec[1][0]]
        tmpRec[2] = [rec[1][1]]
        tmpRec[3] = []
        # SMILES
        tmpRec[4] = [rec[2][0].strip('\"')]
        tmpRec[5] = [rec[2][1].strip('\"')]
        tmpRec[6] = []
        #Time = 
        tmpRec[7] = (rec[3],-1)
        tmpRec[0] = getReactionID(tmpRec[4],tmpRec[5],tmpRec[6])
        RegulateRec.append(tmpRec)
        
    ReactRecDict = {}
    
    for itm in RegulateRec:
        if(itm[0] in ReactRecDict):
            ReactRecDict[itm[0]].append(itm)
        else:
            ReactRecDict[itm[0]] = [itm]
            
    ReactRecDictNew =  {}         
    for k in sorted(ReactRecDict, key=lambda k: len(ReactRecDict[k]),reverse=True):
        ReactRecDictNew[k] = ReactRecDict[k]

    ReactRecDict = ReactRecDictNew
    return(ReactRecDict)
    


def generateReactions(G):
    G_t = deepcopy(G)
    Tlag = gvar.StableMolLag 
    #Tlag = 40 
    # Remove Wrong lable
    for edge in  G_t.edges(data=True,keys=True):
        for key in G_t[edge[0]][edge[1]]:
            if(isinstance(G_t[edge[0]][edge[1]][key]['label'], str)):
                G_t[edge[0]][edge[1]][key]['label'] = int(G_t[edge[0]][edge[1]][key]['label'].strip('\"'))


    recs = buildBasicInfo(G_t)
    N_nodes = G_t.number_of_nodes()
    N_edges = G_t.number_of_edges()
    print("total nodes in original graph1: ",N_nodes)
    print("total edges in original graph1: ",N_edges)

    # Count A+B -> C+D and transofrmation reactions
    CombList = CombineRedBlue(G_t,recs,Tlag)
    removeList,NodermList,SR,LR = SimpleNodeshortReactAbs(CombList,Tlag,G_t)
    greyR,removeList = GreyReactAbs(G_t,removeList)

    rmList_UPD = addKeyToRMlist(G_t,removeList)
    rmList_UPD = list(set(rmList_UPD))
    G_t = rmEdgeList(G_t,rmList_UPD)
    rmSinglenode(G_t)

    # Count A->B+C and A+B-> C reactions
    recs = buildBasicInfo(G_t)
    CombList = CombineRedBlue(G_t,recs,Tlag)
    CR, SpR,removeList = CombAndSplitReactABS(G,CombList,[])

    rmList_UPD = addKeyToRMlist(G_t,removeList)
    rmList_UPD = list(set(rmList_UPD))
    G_t = rmEdgeList(G_t,rmList_UPD)
    rmSinglenode(G_t)

    #  Convert to reactionsDictionary
    reactDict = RegulateAndClassify(SR,LR,CR,SpR,greyR)

    N_nodes = G_t.number_of_nodes()
    N_edges = G_t.number_of_edges()
    print("total nodes in original graph1: ",N_nodes)
    print("total edges in original graph1: ",N_edges)
    return reactDict


if __name__ == '__main__':
    G = nx.drawing.nx_pydot.read_dot("Graph2cores_5000.dot")
    G.remove_node('\\n')
    DicreacR =  generateReactions(G)
    for k in DicreacR:
        # print reaction event 1 by 1:
        ievent = 0
        for event in DicreacR[k]:
            EventStr = "   "+str(ievent)+":  "
            ievent = 1+ievent

            for Rhash in event[1]:
                 if ("int" in Rhash): print(event)
            for Phash in event[2]:
                 if ("int" in Phash): print(event)

