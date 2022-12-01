import networkx as nx
import pydot
from copy import deepcopy
from networkx.drawing.nx_pydot import write_dot
from collections import Counter
import time

def Sort(sub_li,iplace): 
    sub_li.sort(key = lambda x: x[iplace]) 
    return sub_li 


def buildBasicInfo(G):
    """Build basic information for nodes.

    Args:
        G (graph) : Reaction graph.
    Returns:
        rec    : List of records or each node:[nodeID,nodelabel,[Inedges],[OutEdges]]
    """
    RawNodeList = G.nodes(data=True)
    rec = []
    for node in RawNodeList:
        # Nodeinfo=[nodeID,nodelabel,[Inedges],[OutEdges]]
        # Get Label
        nodeInfo = [node[0],node[1]['label']]
        
        # GetInEdges:
        InEdge = []
#       StepCount = []
        edgesRawIn = G.in_edges(node[0],data=True)
        for edge in edgesRawIn:
#           StepCount.append(edge[2]['label'])
            InEdge.append((edge[0],edge[2]['color'],edge[2]['label']))
        InEdge = Sort(InEdge,2)
#       InStepCount = Counter(StepCount)
        #print(InEdge)
        nodeInfo.append(InEdge)
        
        # GetOutEdges:
        OutEdge = []
#       StepCount = []
        edgesRawOut = G.out_edges(node[0],data=True)
        for edge in edgesRawOut:
#           StepCount.append(edge[2]['label'])
            OutEdge.append((edge[1],edge[2]['color'],edge[2]['label']))
        OutEdge = Sort(OutEdge,2)
#       OutStepCount = Counter(StepCount)
        nodeInfo.append(OutEdge)

        # determine the type of each node. 
        #    'center' nodes are the nodes combined by others. Generely transit structure.
        #    'edge' nodes are nodes without splitting and combined from edge nodes.         
        Type='edge'
        for i in range(len(nodeInfo[2]) - 1):
            Type='edge'
            if(nodeInfo[2][i][2] == nodeInfo[2][i+1][2] and nodeInfo[2][i][1]=='red' and nodeInfo[2][i+1][1]=='red'):
                Type='center'
                break
            if(nodeInfo[2][i][2] == nodeInfo[2][i+1][2] and nodeInfo[2][i][1]=='blue' and nodeInfo[2][i+1][1]=='blue'):
                Type='center'
                break
        for i in range(len(nodeInfo[3]) - 1):
            Type='edge'
            if(nodeInfo[3][i][2] == nodeInfo[3][i+1][2] and nodeInfo[3][i][1]=='red' and nodeInfo[3][i+1][1]=='red'):
                Type='center'
                break
            if(nodeInfo[3][i][2] == nodeInfo[3][i+1][2] and nodeInfo[3][i][1]=='blue' and nodeInfo[3][i+1][1]=='blue'):
                Type='center'
                break
        nodeInfo.append(Type)
#       nodeInfo.append([InStepCount,OutStepCount])
        rec.append(nodeInfo)
    return rec
    
    
def getConjInfo(G,Totrec):
    """Append conjunction information to basic information list.

    Args:
        G (graph) : Reaction graph.
        Totrac    : information list for nodes.
    Returns:
        Totrac    : information list for nodes with conjunction information.
        conjDict  : dictionary for conjunction information.
    """
    Conj = []
    conjDict={}
    for rec in Totrec:
        redTo = []
        for outnode in rec[3]:
            if(outnode[1] == 'red'):
                redTo.append((outnode[0],outnode[2]))
        #Build conj. List
        conjList = []
        for arrow in redTo:
            edgesRawIn = G.in_edges(arrow[0],data=True)
            for nodeIn in edgesRawIn:
                if( int(nodeIn[2]['label']) == arrow[1] 
                    and nodeIn[2]['color'] == 'red'                   
                    and nodeIn[0] != rec[0]):
                    conjList.append( (nodeIn[0],arrow[1]) )
        conjDict[rec[0]]=conjList
        rec.append(conjList)
    return Totrec,conjDict

def removeFakegrey(Totrec):
    """Append conjunction information to basic information list.

    Args:
        G (graph) : Reaction graph.
        Totrac    : information list for nodes.
    Returns:
        Totrac    : information list for nodes with conjunction information.
        conjDict  : dictionary for conjunction information.
    """
    rmList = []
    for rec in Totrec:
        # Fix the inherit bug of 3 molecules collision
        if(len(rec[2]) > 2):
            # record of step #. 
            li=[]
            for itm in rec[2]:
                li.append(itm[2])
                largeNum = 0
            for number in li:
                if li.count(number) > 2: 
                    largeNum = number
                    break
            # Remove the fake grey transformation.
            for itm in rec[2]:
                if(itm[1] == 'grey' and itm[2] == largeNum):
                    rmList.append((itm[0],rec[0],itm[2]))

    return(rmList)
                
                

def rmSinglenode(G):                                    
    """Remove isolated node

    Args:
        G (graph) : Reaction graph.
    Returns:
        G (graph) : Updated graph
    """
    RawNodeList = G.nodes()
    rmNodeList = []
    for node in RawNodeList:                            
        edgesRawIn = G.in_edges(node,data=True)      
        edgesRawOut = G.out_edges(node,data=True)
        if(len(edgesRawIn)==0 and len(edgesRawOut)==0): 
            rmNodeList.append(node)
    for node in rmNodeList:
        G.remove_node(node)
            
            
def rmEdgeList(G,rmList):
    """Remove isolated node

    Args:
        G (graph) : Reaction graph.
        rmList    : remove edges list
    Returns:
        G (graph) : Updated graph.
    """
    count = 0
    for itm in rmList:
        count = count +1
        G.remove_edge(itm[0],itm[1],key=itm[3])
    return G
            
    
def addKeyToRMlist(G,rmList):
    """add key information to remove list

    Args:
        G (graph) : Reaction graph 
        rmList    : remove edges list
    Returns:
        rmListUPD : Updated remove list 
    """
    rmListUPD = []
    for itm in rmList:
        dicts = G.get_edge_data(itm[0], itm[1])
        if(dicts is None): print(itm)
        for key in dicts:
            if(dicts[key]['label'] == itm[2]):
                rmListUPD.append((itm[0],itm[1],itm[2],key))
    return(rmListUPD)

def addKeyToRMlist_grey(G,rmList):
    """add key information to remove list
       (For fake gray specially)
    Args:
        G (graph) : Reaction graph 
        rmList    : remove edges list
    Returns:
        rmListUPD : Updated remove list 
    """
    rmListUPD = []
    for itm in rmList:
        dicts = G.get_edge_data(itm[0], itm[1])
        if(dicts is None): print(itm)
        for key in dicts:
            if(dicts[key]['label'] == itm[2] and dicts[key]['color'] == 'grey'  ):
                rmListUPD.append((itm[0],itm[1],itm[2],key))
    return(rmListUPD)

def reaction_filter_grey_forward_direct(G,totrec,step_cri=30):
    """Remove grey edge in forward direct way

     Out-edge start remove Grey edges. Remove Edge in straight way:
      2    4    6
      |    |    |
    ->1    3    5
    Args:
        G (graph) : Reaction graph.
        totrec    : Node Information Records.
        step_cri  : Remove critia for steps. default: 30 steps
    Returns:
        rmList    : List of edges to remove.
    """
    rmList = []
    grayPass = []
    for rec in totrec:
        if (rec[0] in grayPass):
            continue
        recTemp = deepcopy(rec)
        # remove red and blue from record, Only grey transformation kept. 
        for trans in rec[2]:
            if(trans[1] == 'red' or trans[1] == 'blue'):
                recTemp[2].remove(trans)
        for trans in rec[3]:
            if(trans[1] == 'red' or trans[1] == 'blue'):
                recTemp[3].remove(trans)
            
        if( len(recTemp[2]) == 0 or len(recTemp[3])==0):
            continue
        rmList_t = []
        for i in range(min( len(recTemp[2]),len(recTemp[3]) )):
            if (recTemp[2][i][2] -  recTemp[3][i][2] < step_cri and
                recTemp[2][i][2] >  recTemp[3][i][2] and
                recTemp[3][i][0] == recTemp[2][i][0]):
                rmList_t.append( (recTemp[2][i][0],recTemp[0],recTemp[2][i][2],"straight") )
                rmList_t.append( (recTemp[0],recTemp[3][i][0],recTemp[3][i][2],"straight") )
                grayPass.append(  recTemp[2][i][0])
        for itm in rmList_t:
            rmList.append((itm[0],itm[1],itm[2]))
    return rmList

def reaction_filter_grey_forward_zigzag(G,totrec,step_cri=30):
    """Remove grey edge in forward zigzag way

    out-edge start remove Grey edges. Remove Edge in zigzag way:
    ->1   3   5
        /   /      
      2   4   6
    Args:
        G (graph) : Reaction graph.
        totrec    : Node Information Records.
        step_cri  : Remove critia for steps. default: 30 steps
    Returns:
        rmList    : List of edges to remove.
    """
    rmList = []
    grayPass = []
    for rec in totrec:
        if (rec[0] in grayPass):
            continue
        recTemp = deepcopy(rec)
        # remove red and blue from record, Only grey transformation kept. 
        for trans in rec[2]:
            if(trans[1] == 'red' or trans[1] == 'blue'):
                recTemp[2].remove(trans)
        for trans in rec[3]:
            if(trans[1] == 'red' or trans[1] == 'blue'):
                recTemp[3].remove(trans)
            
        if( len(recTemp[2]) == 0 or len(recTemp[3])==0):
            continue

        rmList_t = []
        for i in range(min( len(recTemp[2])-1,len(recTemp[3]) )):
            if (recTemp[2][i+1][2] -  recTemp[3][i][2] < step_cri and
                recTemp[2][i][2]   <  recTemp[3][i][2] and
                recTemp[3][i][0] == recTemp[2][i][0]):
                
                rmList_t.append( (recTemp[0],recTemp[3][i][0],recTemp[3][i][2],"zigzag" ) )
                rmList_t.append( (recTemp[2][i+1][0],recTemp[0],recTemp[2][i+1][2],"zigzag" ) )
                grayPass.append(  recTemp[3][i][0])

               # print("===========")
               # print((recTemp[0],recTemp[3][i][0],recTemp[3][i][2],"zigzag" ))
               # print((recTemp[2][i+1][0],recTemp[0],recTemp[2][i+1][2],"zigzag" ))
               # print("===========")
        for itm in rmList_t:
            rmList.append((itm[0],itm[1],itm[2]))
    return rmList

def reaction_filter_grey_reversed_direct(G,totrec,step_cri=30):
    """Remove grey edge in reversed direct way

    out-edge start remove Grey edges. Remove Edge in zigzag way:
       -> 1    3    5
          |    |    |
          2    4    6
    Args:
        G (graph) : Reaction graph.
        totrec    : Node Information Records.
        step_cri  : Remove critia for steps. default: 30 steps
    Returns:
        rmList    : List of edges to remove.
    """
    rmList = []
    grayPass = []
    for rec in totrec:
        if (rec[0] in grayPass):
            continue
        recTemp = deepcopy(rec)
        # remove red and blue from record, Only grey transformation kept. 
        for trans in rec[2]:
            if(trans[1] == 'red' or trans[1] == 'blue'):
                recTemp[2].remove(trans)
        for trans in rec[3]:
            if(trans[1] == 'red' or trans[1] == 'blue'):
                recTemp[3].remove(trans)
            
        if( len(recTemp[2]) == 0 or len(recTemp[3])==0):
            continue
        rmList_t = []
        for i in range(min( len(recTemp[2]),len(recTemp[3]) )):
            if (recTemp[3][i][2] -  recTemp[2][i][2] < step_cri and
                recTemp[3][i][2] >  recTemp[2][i][2] and
                recTemp[3][i][0] == recTemp[2][i][0]):
                rmList_t.append( (recTemp[2][i][0],recTemp[0],recTemp[2][i][2],"straight") )
                rmList_t.append( (recTemp[0],recTemp[3][i][0],recTemp[3][i][2],"straight") )
                grayPass.append(  recTemp[2][i][0])
        for itm in rmList_t:
            rmList.append((itm[0],itm[1],itm[2]))
    return rmList


def reaction_filter_grey_reversed_zigzag(G,totrec,step_cri=30):
    """Remove grey edge in reversed zigzag way

    Remove Edge in zigzag way:
      ->  1    3    5
             /    /      
          2    4    6
    Args:
        G (graph) : Reaction graph.
        totrec    : Node Information Records.
        step_cri  : Remove critia for steps. default: 30 steps
    Returns:
        rmList    : List of edges to remove.
    """
    rmList = []
    grayPass = []
    for rec in totrec:
        if (rec[0] in grayPass):
            continue
        recTemp = deepcopy(rec)
        # remove red and blue from record, Only grey transformation kept. 
        for trans in rec[2]:
            if(trans[1] == 'red' or trans[1] == 'blue'):
                recTemp[2].remove(trans)
        for trans in rec[3]:
            if(trans[1] == 'red' or trans[1] == 'blue'):
                recTemp[3].remove(trans)
            
        if( len(recTemp[2]) == 0 or len(recTemp[3])==0):
            continue

        rmList_t = []
        for i in range(min( len(recTemp[2]),len(recTemp[3])-1 )):
            if (recTemp[3][i+1][2] -  recTemp[2][i][2] < step_cri and
                recTemp[2][i][2]   >  recTemp[3][i][2] and
                recTemp[3][i+1][0] == recTemp[2][i][0]):
                
                rmList_t.append( (recTemp[0],recTemp[3][i+1][0],recTemp[3][i+1][2],"zigzag" ) )
                rmList_t.append( (recTemp[2][i][0],recTemp[0],recTemp[2][i][2],"zigzag" ) )
                grayPass.append(  recTemp[3][i][0])

        for itm in rmList_t:
            rmList.append((itm[0],itm[1],itm[2]))
    return rmList


def reaction_filter_br_forward_direct(G,totrec,step_cri=30):
    """Remove blue and red edges in forward direct way.

    Remove Edge in straight way:
          2b    4b    6b
          |     |     |
        ->1r    3r    5r
    Args:
        G (graph) : Reaction graph.
        totrec    : Node Information Records.
        step_cri  : Remove critia for steps. default: 30 steps
    Returns:
        rmList    : List of edges to remove.
    """
    rmList = []
    rbPass = []

    for rec in totrec:
        if (rec[0] in rbPass or rec[4] =='center'):
            continue
        recTemp = deepcopy(rec)
        # remove grey from record, Only red/blue transformation kept. 
        for trans in rec[2]:
            if(trans[1] == 'grey'):
                recTemp[2].remove(trans)
                
        for trans in rec[3]:
            if(trans[1] == 'grey'):
                recTemp[3].remove(trans)

        if( len(recTemp[2]) == 0 or len(recTemp[3])==0):
            continue
        rmList_t = []
        for i in range(min( len(rec[2]),len(rec[3]) )):
            if (rec[2][i][2]-rec[3][i][2] < step_cri and
                rec[2][i][2]>rec[3][i][2] and 
                rec[2][i][1] == 'blue' and 
                rec[3][i][1] == 'red' and
                rec[3][i][0] == rec[2][i][0]):                

                # Detect the removal actions.
                conj_tag1 = None
                for rec_conj in rec[5]:
                    # same step out-edge
                    if(rec_conj[1] == rec[3][i][2]):
                        conj_tag1 = rec_conj[0]
                        break
                if(conj_tag1 is None):
                    inList = [] 
                else:
                    inList = G.in_edges(conj_tag1,data=True)
                inList = G.in_edges(conj_tag1,data=True)
                remove = False
                for itm in inList:
                    if(itm[0] == rec[2][i][0] and itm[1] == conj_tag1 and int(itm[2]['label']) == rec[2][i][2]):
                        remove = True
                        break
                if (remove):
                    rmList_t.append( (rec[2][i][0],rec[0],rec[2][i][2],"straight") )
                    rmList_t.append( (rec[0],rec[3][i][0],rec[3][i][2],"straight") )
                    rmList_t.append( (rec[2][i][0],conj_tag1,rec[2][i][2],"straight") )
                    rmList_t.append( (conj_tag1,rec[3][i][0],rec[3][i][2],"straight") )
                    
        for itm in rmList_t:
            rmList.append((itm[0],itm[1],itm[2]))
    return rmList


def reaction_filter_br_forward_zigzag(G,totrec,step_cri=30):
    """Remove blue and red edges in forward zigzag way.

     Remove Edge in Zigzag way:
    ->1b  3b  5b
         /  /    
      2r  4r  6r
    Args:
        G (graph) : Reaction graph.
        totrec    : Node Information Records.
        step_cri  : Remove critia for steps. default: 30 steps
    Returns:
        rmList    : List of edges to remove.
    """
    rmList = []
    rbPass = []
    for rec in totrec:
        if (rec[0] in rbPass or rec[4] =='center'):
            continue
        recTemp = deepcopy(rec)
        # remove grey from record, Only red/blue transformation kept. 
        for trans in rec[2]:
            if(trans[1] == 'grey'):
                recTemp[2].remove(trans)
                
        for trans in rec[3]:
            if(trans[1] == 'grey'):
                recTemp[3].remove(trans)

        if( len(recTemp[2]) == 0 or len(recTemp[3])==0):
            continue
        rmList_t = []
        for i in range(min( len(rec[2])-1,len(rec[3]) )):
            if (rec[2][i+1][2] - rec[3][i][2] < step_cri and
                rec[2][i][2] < rec[3][i][2] and
                rec[3][i][1] == 'red' and 
                rec[2][i+1][1] == 'blue' and
                rec[3][i][0] == rec[2][i+1][0]):

                # Detect the removal actions. 
                conj_tag1 = None
                for rec_conj in rec[5]:
                    # same step out-edge
                    if(rec_conj[1] == rec[3][i][2]):
                        conj_tag1 = rec_conj[0]
                        break
                if(conj_tag1 is None):
                    inList = [] 
                else:
                    inList = G.in_edges(conj_tag1,data=True)
                remove = False
                
                for itm in inList:
                    if(itm[0] == rec[2][i+1][0] and itm[1] == conj_tag1 and int(itm[2]['label']) == rec[2][i+1][2]):
                        remove = True
                        break

                if (remove):
                    rmList_t.append( (rec[0],rec[3][i][0],rec[3][i][2],"zigzag" ) )
                    rmList_t.append( (rec[2][i+1][0],rec[0],rec[2][i+1][2],"zigzag" ) )
                    rmList_t.append( (conj_tag1,rec[3][i][0],rec[3][i][2],"zigzag" ) )
                    rmList_t.append( (rec[2][i+1][0],conj_tag1,rec[2][i+1][2],"zigzag" ) )
                    rbPass.append(conj_tag1)                
                    
        for itm in rmList_t:
            rmList.append((itm[0],itm[1],itm[2]))
    return rmList


def reaction_filter_br_reversed_direct(G,totrec,step_cri=30):
    """Remove blue and red edges in reversed direct way.

     Remove Edges in straight way:
     ->1b  3b  5b
       |   |   |  
       2r  4r  6r
    Args:
        G (graph) : Reaction graph.
        totrec    : Node Information Records.
        step_cri  : Remove critia for steps. default: 30 steps
    Returns:
        rmList    : List of edges to remove.
    """
    rmList = []
    rbPass = []
    for rec in totrec:
        if (rec[0] in rbPass or rec[4] =='center'):
            continue
        recTemp = deepcopy(rec)
        # remove grey from record, Only red/blue transformation kept. 
        for trans in rec[2]:
            if(trans[1] == 'grey'):
                recTemp[2].remove(trans)
                
        for trans in rec[3]:
            if(trans[1] == 'grey'):
                recTemp[3].remove(trans)

        if( len(recTemp[2]) == 0 or len(recTemp[3])==0):
            continue
        rmList_t = []
        for i in range(min( len(rec[2]),len(rec[3]) )):
            if (rec[3][i][2]-rec[2][i][2] < step_cri and
                rec[3][i][2]>rec[2][i][2] and 
                rec[2][i][1] == 'blue' and 
                rec[3][i][1] == 'red' and
                rec[3][i][0] == rec[2][i][0]):           

                # Detect the removal actions.
                conj_tag1 = None
                for rec_conj in rec[5]:
                    # same step out-edge
                    if(rec_conj[1] == rec[3][i][2]):
                        conj_tag1 = rec_conj[0]
                        break
                if(conj_tag1 is None):
                    inList = [] 
                else:
                    inList = G.in_edges(conj_tag1,data=True)
                remove = False
                for itm in inList:
                    if(itm[0] == rec[2][i][0] and itm[1] == conj_tag1 and int(itm[2]['label']) == rec[2][i][2]):
                        remove = True
                        break
                if (remove):
                    rmList_t.append( (rec[2][i][0],rec[0],rec[2][i][2],"straight") )
                    rmList_t.append( (rec[0],rec[3][i][0],rec[3][i][2],"straight") )
                    rmList_t.append( (rec[2][i][0],conj_tag1,rec[2][i][2],"straight") )
                    rmList_t.append( (conj_tag1,rec[3][i][0],rec[3][i][2],"straight") )
                    
        for itm in rmList_t:
            rmList.append((itm[0],itm[1],itm[2]))
    return rmList        


def reaction_filter_br_reversed_zigzag(G,totrec,step_cri=30):
    """Remove blue and red edges in reversed zigzag way.

    Out-edge start remove Grey edges. Remove Edge in zigzag way:
      2b   4b   6b
         \    \        
    ->1r   3r   5r
    Args:
        G (graph) : Reaction graph.
        totrec    : Node Information Records.
        step_cri  : Remove critia for steps. default: 30 steps
    Returns:
        rmList    : List of edges to remove.
    """
    rmList = []
    rbPass = []
    for rec in totrec:
        if (rec[0] in rbPass or rec[4] =='center'):
            continue
        recTemp = deepcopy(rec)
        # remove grey from record, Only red/blue transformation kept. 
        for trans in rec[2]:
            if(trans[1] == 'grey'):
                recTemp[2].remove(trans)
                
        for trans in rec[3]:
            if(trans[1] == 'grey'):
                recTemp[3].remove(trans)

        if( len(recTemp[2]) == 0 or len(recTemp[3])==0):
            continue
        rmList_t = []
        for i in range(min( len(rec[2]),len(rec[3])-1 )):
            if (rec[3][i+1][2] - rec[2][i][2] < step_cri and
                rec[2][i][2] > rec[3][i][2] and
                rec[3][i+1][1] == 'red' and 
                rec[2][i][1] == 'blue' and
                rec[3][i+1][0] == rec[2][i][0]):

                # Detect the removal actions. 
                conj_tag1 = None
                for rec_conj in rec[5]:
                    # same step out-edge
                    if(rec_conj[1] == rec[3][i+1][2]):
                        conj_tag1 = rec_conj[0]
                        break
                if(conj_tag1 is None):
                    inList = [] 
                else:
                    inList = G.in_edges(conj_tag1,data=True)

                remove = False
                
                for itm in inList:
                    if(itm[0] == rec[2][i][0] and itm[1] == conj_tag1 and int(itm[2]['label']) == rec[2][i][2]):
                        remove = True
                        break

                if (remove):
                    rmList_t.append( (rec[0],rec[3][i+1][0],rec[3][i+1][2],"zigzag" ) )
                    rmList_t.append( (rec[2][i][0],rec[0],rec[2][i][2],"zigzag" ) )
                    rmList_t.append( (conj_tag1,rec[3][i+1][0],rec[3][i+1][2],"zigzag" ) )
                    rmList_t.append( (rec[2][i][0],conj_tag1,rec[2][i][2],"zigzag" ) )
                    rbPass.append(conj_tag1)                
                    
        for itm in rmList_t:
            rmList.append((itm[0],itm[1],itm[2]))
    return rmList

def remove_useless_trans(G,Tlag=10):
    G_t = deepcopy(G)
    N_nodes0 = G.number_of_nodes()
    N_edges0 = G.number_of_edges()
    # Remove fake grey edges:
    rec = buildBasicInfo(G_t)
    Totrec,conjDic = getConjInfo(G_t,rec)
    rmList = removeFakegrey(Totrec)
    rmList_UPD = addKeyToRMlist_grey(G_t,rmList)
    rmList_UPD = list(set(rmList_UPD))
    G_t = rmEdgeList(G_t,rmList_UPD)
    Tinfo = 0
    Trm   = 0
    conv = False
    count=0
    while not conv:
        count=count+1
        T00 = time.time()
        rec = buildBasicInfo(G_t)
        Totrec,conjDic = getConjInfo(G_t,rec)
        rmList = reaction_filter_grey_forward_direct(G_t,Totrec,Tlag)
        T01 = time.time()
        rmList_UPD = addKeyToRMlist(G_t,rmList)
        rmList_UPD = list(set(rmList_UPD))
        G_t = rmEdgeList(G_t,rmList_UPD)
        T02 = time.time()
        Tinfo = Tinfo+(T01-T00)
        Trm = Trm+(T02-T01)
        N_edges = G_t.number_of_edges()
#       print(N_edges,N_edges0)
        if(N_edges-N_edges0 == 0): 
            conv = True
        else:
            conv = False
            N_edges0 = N_edges
        
#   print(count,"  ===1===")
    
    
    conv = False
    count=0
    while not conv:
        count=count+1
        T00 = time.time()
        rec = buildBasicInfo(G_t)
        Totrec,conjDic = getConjInfo(G_t,rec)
        rmList = reaction_filter_grey_reversed_direct(G_t,Totrec,Tlag)
        T01 = time.time()
        rmList_UPD = addKeyToRMlist(G_t,rmList)
        rmList_UPD = list(set(rmList_UPD))
        G_t = rmEdgeList(G_t,rmList_UPD)
        T02 = time.time()
        Tinfo = Tinfo+(T01-T00)
        Trm = Trm+(T02-T01)
        N_edges = G_t.number_of_edges()
#       print(N_edges,N_edges0)
        if(N_edges-N_edges0 == 0): 
            conv = True
        else:
            conv = False
            N_edges0 = N_edges
        
#   print(count,"  ===2===")    
    conv = False
    count=0
    while not conv:
        count=count+1
        T00 = time.time()
        rec = buildBasicInfo(G_t)
        Totrec,conjDic = getConjInfo(G_t,rec)
        rmList = reaction_filter_grey_forward_zigzag(G_t,Totrec,Tlag)
        T01 = time.time()
        rmList_UPD = addKeyToRMlist(G_t,rmList)
        rmList_UPD = list(set(rmList_UPD))
        G_t = rmEdgeList(G_t,rmList_UPD)
        T02 = time.time()
        Tinfo = Tinfo+(T01-T00)
        Trm = Trm+(T02-T01)
        N_edges = G_t.number_of_edges()
#       print(N_edges,N_edges0)
        if(N_edges-N_edges0 == 0): 
            conv = True
        else:
            conv = False
            N_edges0 = N_edges
        
#   print(count,"  ===3===")    


    conv = False
    count=0
    while not conv:
        count=count+1
        T00 = time.time()
        rec = buildBasicInfo(G_t)
        Totrec,conjDic = getConjInfo(G_t,rec)
        rmList = reaction_filter_grey_reversed_zigzag(G_t,Totrec,Tlag)
        T01 = time.time()
        rmList_UPD = addKeyToRMlist(G_t,rmList)
        rmList_UPD = list(set(rmList_UPD))
        G_t = rmEdgeList(G_t,rmList_UPD)
        T02 = time.time()
        Tinfo = Tinfo+(T01-T00)
        Trm = Trm+(T02-T01)
        N_edges = G_t.number_of_edges()
#       print(N_edges,N_edges0)
        if(N_edges-N_edges0 == 0): 
            conv = True
        else:
            conv = False
            N_edges0 = N_edges
#   print(count,"  ===4===")    


    conv = False
    count=0
    while not conv:
        count=count+1
        T00 = time.time()
        rec = buildBasicInfo(G_t)
        Totrec,conjDic = getConjInfo(G_t,rec)
        rmList = reaction_filter_br_forward_direct(G_t,Totrec,Tlag)
        T01 = time.time()
        rmList_UPD = addKeyToRMlist(G_t,rmList)
        rmList_UPD = list(set(rmList_UPD))
        G_t = rmEdgeList(G_t,rmList_UPD)
        T02 = time.time()
        Tinfo = Tinfo+(T01-T00)
        Trm = Trm+(T02-T01)
        N_edges = G_t.number_of_edges()
#       print(N_edges,N_edges0)
        if(N_edges-N_edges0 == 0): 
            conv = True
        else:
            conv = False
            N_edges0 = N_edges
        
#   print(count,"  ===1===")

    conv = False
    count=0
    while not conv:
        count=count+1
        T00 = time.time()
        rec = buildBasicInfo(G_t)
        Totrec,conjDic = getConjInfo(G_t,rec)
        rmList = reaction_filter_br_reversed_direct(G_t,Totrec,Tlag)
        T01 = time.time()
        rmList_UPD = addKeyToRMlist(G_t,rmList)
        rmList_UPD = list(set(rmList_UPD))
        G_t = rmEdgeList(G_t,rmList_UPD)
        T02 = time.time()
        Tinfo = Tinfo+(T01-T00)
        Trm = Trm+(T02-T01)
        N_edges = G_t.number_of_edges()
#       print(N_edges,N_edges0)
        if(N_edges-N_edges0 == 0): 
            conv = True
        else:
            conv = False
            N_edges0 = N_edges
        
#   print(count,"  ===2===") 
    conv = False
    count=0
    while not conv:
        count=count+1
        T00 = time.time()
        rec = buildBasicInfo(G_t)
        Totrec,conjDic = getConjInfo(G_t,rec)
        rmList = reaction_filter_br_forward_zigzag(G_t,Totrec,Tlag)
        T01 = time.time()
        rmList_UPD = addKeyToRMlist(G_t,rmList)
        rmList_UPD = list(set(rmList_UPD))
        G_t = rmEdgeList(G_t,rmList_UPD)
        T02 = time.time()
        Tinfo = Tinfo+(T01-T00)
        Trm = Trm+(T02-T01)
        N_edges = G_t.number_of_edges()
#       print(N_edges,N_edges0)
        if(N_edges-N_edges0 == 0): 
            conv = True
        else:
            conv = False
            N_edges0 = N_edges
        
#   print(count,"  ===3===")    
   

    conv = False
    count=0
    while not conv:
        count=count+1
        T00 = time.time()
        rec = buildBasicInfo(G_t)
        Totrec,conjDic = getConjInfo(G_t,rec)
        rmList = reaction_filter_br_reversed_zigzag(G_t,Totrec,Tlag)
        T01 = time.time()
        rmList_UPD = addKeyToRMlist(G_t,rmList)
        rmList_UPD = list(set(rmList_UPD))
        G_t = rmEdgeList(G_t,rmList_UPD)
        T02 = time.time()
        Tinfo = Tinfo+(T01-T00)
        Trm = Trm+(T02-T01)
        N_edges = G_t.number_of_edges()
#       print(N_edges,N_edges0)
        if(N_edges-N_edges0 == 0): 
            conv = True
        else:
            conv = False
            N_edges0 = N_edges
#   print(count,"  ===4===")
    rmSinglenode(G_t)
    return G_t



def nodeShortRemove(G,Tlag):
    # Combine the red and blue together for analysis.
    Recs = buildBasicInfo(G)

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
                if(len(Combinetup_temp)>2):
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
                if(len(Combinetup_temp)>2):
                    Combinetup_temp = SortRecitm(Combinetup_temp)
                    CombineBlue.append(Combinetup_temp)
            CombineBlue = Sort(CombineBlue,0)
            CenterRec.append(CombineBlue)


            #print(CombineBlue,"blue")
            if(len(CenterRec[1])>0 and len(CenterRec[2])>0):
                CenterRecTot.append(CenterRec)
            
    # Decide Remove or not?
    removeList = [] 
    for rec in CenterRecTot:
        nLoop = min(len(rec[1]),len(rec[2]))
        for i in range(nLoop):
            Flag = False
            if(abs(rec[2][i][0] - rec[1][i][0]) < Tlag and len(rec[2][i]) == len(rec[1][i]) and rec[2][i][0] > rec[1][i][0]):
                Flag = True
                for j in range(len(rec[2][i])-1  ):
                    if(rec[2][i][j+1] != rec[1][i][j+1] ):
                        Flag = False
                        break
            if(Flag):
                for itm in range(len(rec[1][i])-1):
                    removeList.append((rec[1][i][itm+1],rec[0],rec[1][i][0]))
                for itm in range(len(rec[2][i])-1):
                    removeList.append((rec[0],rec[2][i][itm+1],rec[2][i][0]))
                    

    rmList_UPD = addKeyToRMlist(G,removeList)
    rmList_UPD = list(set(rmList_UPD))
    G_t = rmEdgeList(G,rmList_UPD)
    rmSinglenode(G_t)
    return(G_t) 

if __name__ == '__main__':
    G = nx.drawing.nx_pydot.read_dot("/home/xchen/Develop/CPX-MechGen_wasted/reactionGraphsmall56-G.data")
    dic = nx.get_node_attributes(G,'label') 
    N_nodes0 = G.number_of_nodes()
    N_edges0 = G.number_of_edges()
    print("total nodes in original graph: ",N_nodes0)
    print("total edges in original graph: ",N_edges0)
    G_t = remove_useless_trans(G)
    N_nodes = G_t.number_of_nodes()
    N_edges = G_t.number_of_edges()
    print("total nodes in original graph: ",N_nodes)
    print("total edges in original graph: ",N_edges)
    write_dot(G_t, "reduce_small56.dot") 

