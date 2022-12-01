import networkx as nx
import pydot
from copy import deepcopy
from networkx.drawing.nx_pydot import write_dot
import time
import global_var as gvar

def contract_grey(G,passPair,lag_cri=30):
    """Contract short lifetime grey contraction.

    Args:
        G (graph) : Reaction graph.
        lag_cri : Lifetime criteria. (default=30)
        passPair: Transformation pairs need to be neglected
    Returns:
        contrac_result: Grey transformation contaction Suggestion,[(a,b),(e,f), .....]
        currentpassPair:Pair should be passed.
    """
    
    #1. Get all node pairs contain grey transformations.
    greyList =[]
    # [(node1,node2),(node_n,node_n-1),....]
    # 
    RawEdgeList = G.edges(data=True)
    for rec in RawEdgeList:
        if(rec[2]['color'] == 'grey' and \
          (rec[0],rec[1]) not in greyList and \
          (rec[1],rec[0]) not in greyList and \
          (rec[0],rec[1]) not in passPair and \
          (rec[1],rec[0]) not in passPair and
          rec[1] != rec[0]):        
            greyList.append((rec[0],rec[1]))
    
    #2. Build node information dictionary.
    builtRec={}
    for pair in greyList:
        for inode in pair:
            
            if(inode in builtRec):
                continue
            
            rec_in = G.in_edges(inode,data=True)
            rec_out = G.out_edges(inode,data=True)
            #print(rec_in,rec_out)
            
            GreyDicRec=[]
            # Build In record
            for rec1 in rec_in:
                if(rec1[2]['color']=='grey'):
                    Instep = int(rec1[2]['label'])
                    mini_sub = 200000
                    iloc = 0
                    rec_out_l = list(rec_out)
                    flag = False
                    
                    for i in range(len(rec_out_l)):
                        Outstep = int(rec_out_l[i][2]['label'])
                        if(Outstep - Instep >0 and Outstep - Instep < mini_sub):
                            iloc = i
                            mini_sub = Outstep - Instep
                            flag = True 
                    if(flag):    
                        GreyDicRec.append([rec_out_l[iloc][1], "in",Instep ,int(rec_out_l[iloc][2]['label'])])
                    else:
                        if(len(rec_out_l) == 0): 
                            GreyDicRec.append([-1, "in",Instep ,None])
                        else:
                            GreyDicRec.append([rec_out_l[iloc][1], "in",Instep ,None])
            # Build Out record
            for rec1 in rec_out:
                if(rec1[2]['color']=='grey'):
                    Outstep = int(rec1[2]['label'])
                    mini_sub = 200000
                    iloc = 0
                    rec_in_l = list(rec_in)
                    flag = False
                    
                    for i in range(len(rec_in_l)):
                        Instep = int(rec_in_l[i][2]['label'])
                        if(Outstep - Instep >0 and Outstep - Instep < mini_sub):
                            iloc = i
                            mini_sub = Outstep - Instep
                            flag = True
                    if(flag):    
                        GreyDicRec.append([rec_in_l[iloc][0],"Out",Outstep ,int(rec_in_l[iloc][2]['label']) ])
                    else:
                        if(len(rec_in_l) == 0): 
                            GreyDicRec.append([-1, "Out",Outstep ,None])
                        else:
                            GreyDicRec.append([rec_in_l[iloc][0],"Out",Outstep ,None])         
            builtRec[inode] = GreyDicRec
    #for key in builtRec:
        #print(key,builtRec[key])
    #3. Pair Opinions
    # pair1:(a-b;a-b;b-a)
    # pair2:(c-d;c-d;none)
    #print("start pair deciding")
    passRec = []
    pairOpinionList = []
    for pair in greyList:
        # 3.1 Add search pair ot pass pair.
        if (pair[0] in passRec or pair[1] in passRec):
            continue
        passRec.append(pair[0])
        passRec.append(pair[1])
        
        # 3.2 start getting opinion.
        edge1_lifeTime = 10000
        edge2_lifeTime = 10000
        
        opinion_list = []
        
        # 0->1
        #print(lag_cri,"lag_cri")
        edges_1 = G.get_edge_data(pair[0],pair[1])
        if(edges_1 == None): edges_1=[]
        for edge_key in edges_1:
            grey_step = int(edges_1[edge_key]['label'])
            # Life time of 0
            opinion =None
            for rec_t in builtRec[pair[0]]:
                #print (rec_t[1],rec_t[2])
                if(rec_t[3] is None and rec_t[2] == grey_step):
                    edge1_lifeTime = lag_cri+10
                    break
                
                if(rec_t[1] == 'Out' and rec_t[2] == grey_step):
                    edge1_lifeTime = rec_t[2] - rec_t[3]
            # Life time of 1
        
            for rec_t in builtRec[pair[1]]:
                if(rec_t[2] == grey_step and rec_t[3] is None):
                    edge2_lifeTime = lag_cri+10
                    break
                
                if(rec_t[1] == 'in' and rec_t[2] == grey_step):
                    edge2_lifeTime = rec_t[3] - rec_t[2]
                    #print(edge2_lifeTime)

            if(opinion =="no"): 
                opinion ="no"
                opinion_list.append(opinion)
                break 
            if(opinion is None and edge2_lifeTime > edge1_lifeTime and edge1_lifeTime < lag_cri): 
                opinion ="a->b"
                opinion_list.append(opinion)
            if(opinion is None and edge1_lifeTime >  edge2_lifeTime and edge2_lifeTime < lag_cri): 
                opinion ="b->a"
                opinion_list.append(opinion)
            if(opinion is None and edge1_lifeTime == edge2_lifeTime and edge2_lifeTime < lag_cri \
                 and len(G.nodes[pair[0]]["SMILE"]) < len(G.nodes[pair[1]]["SMILE"])  ): 
                opinion ="b->a"
                opinion_list.append(opinion)
            if(opinion is None and edge1_lifeTime == edge2_lifeTime and edge2_lifeTime < lag_cri \
                 and len(G.nodes[pair[0]]["SMILE"]) > len(G.nodes[pair[1]]["SMILE"])  ): 
                opinion ="a->b"
                opinion_list.append(opinion)
        if('no' in opinion_list):
            pairOpinionList.append( [pair,opinion_list])
            break
        #print("0->1 opinion: ", opinion_list)
                    
        # 1->0
        edges_2 = G.get_edge_data(pair[1],pair[0])
        if(edges_2 == None): edges_2=[]
        for edge_key in edges_2:
            grey_step = int(edges_2[edge_key]['label'])
            # Life time of 0
            opinion =None
            for rec_t in builtRec[pair[0]]:
                #print (rec_t[1],rec_t[2])
                if(rec_t[3] is None and rec_t[2] == grey_step):
                    edge1_lifeTime = lag_cri+10
                    break
                
                if(rec_t[1] == 'in' and rec_t[2] == grey_step):
                    edge1_lifeTime = rec_t[3] - rec_t[2]
            # Life time of 1

            for rec_t in builtRec[pair[1]]:
                
                if(rec_t[3] is None and rec_t[2] == grey_step):
                    edge2_lifeTime = lag_cri+10
                    break
                
                if(rec_t[1] == 'Out' and rec_t[2] == grey_step):
                    edge2_lifeTime = rec_t[2] - rec_t[3]
                    
            if(opinion =="no"): 
                opinion ="no"
                opinion_list.append(opinion)
                break 
            if(opinion is None and edge2_lifeTime >  edge1_lifeTime and edge1_lifeTime < lag_cri): 
                opinion ="a->b"
                opinion_list.append(opinion)
            if(opinion is None and edge1_lifeTime >  edge2_lifeTime and edge2_lifeTime < lag_cri): 
                opinion ="b->a"
                opinion_list.append(opinion)

            if(opinion is None and edge1_lifeTime == edge2_lifeTime and edge2_lifeTime < lag_cri \
                 and len(G.nodes[pair[0]]["SMILE"]) < len(G.nodes[pair[1]]["SMILE"])  ): 
                opinion ="b->a"
                opinion_list.append(opinion)

            if(opinion is None and edge1_lifeTime == edge2_lifeTime and edge2_lifeTime < lag_cri \
                 and len(G.nodes[pair[0]]["SMILE"]) > len(G.nodes[pair[1]]["SMILE"])  ): 
                opinion ="a->b"
                opinion_list.append(opinion)

        if('no' in opinion_list):
            pairOpinionList.append( pair,opinion_list)
            break
            
        pairOpinionList.append( [pair,opinion_list])
    # 4. Final Opinion. 
    contrac_result =  []
    currentpassPair = []
    for opinion in pairOpinionList:
        if ( len(opinion[1]) ==0 ):
            currentpassPair.append((opinion[0][0],opinion[0][1]))
            continue 
        if ('no' in opinion[1] ):
            currentpassPair.append((opinion[0][0],opinion[0][1]))
            continue 
        result = dict((i, opinion[1].count(i)) for i in opinion[1])
        if(len(result)>1):
            if(result['a->b'] >= result['b->a']):
                contrac_result.append((opinion[0][0],opinion[0][1]))
            else:
                contrac_result.append((opinion[0][1],opinion[0][0]))
        else:
            if('a->b'in result):
                contrac_result.append((opinion[0][0],opinion[0][1]))
            else:
                contrac_result.append((opinion[0][1],opinion[0][0]))
    return contrac_result,currentpassPair


def contract_blue(G,lag_cri=5):
    """Contract blue transformation.

    Args:
        G (graph) : Reaction graph.
        lag_cri : Lifetime criteria. (default=30)
    Returns:
        cont_pair: Blue transformation contaction Suggestion,[(a,b),(e,f), .....]
    """
    #1. Get all node pairs contain blue transformations.
    blueList =[] 
    passNode = []
    # [(node1,node2),(node_n,node_n-1),....]
    RawEdgeList = G.edges(data=True)
    for rec in RawEdgeList:
        if(rec[2]['color'] == 'blue' and \
          (rec[0],rec[1])  not in blueList and \
           len(G.out_edges(rec[1])) == 2 and \
           len(G.in_edges(rec[1]))  == 1 and \
           G.number_of_edges(rec[0],rec[1]) == 1):
            if(isinstance(rec[2]['label'], int)):
                blueList.append([rec[0], rec[1], rec[2]['label'] ])
            else:
                blueList.append([rec[0], rec[1], int(rec[2]['label'] ) ])
            
    # 1, Check To_node information.
    cont_pair = []
    for edge in blueList:
        rec_out = G.out_edges(edge[1],data=True)
        if(rec_out is None or edge[0] in passNode or edge[1] in passNode): continue
        if(len(rec_out) == 2):
            rec_out = list(rec_out)
            if (rec_out[0][2]['color'] == 'blue' and \
                rec_out[1][2]['color'] == 'blue' and \
                rec_out[0][2]['label'] ==  rec_out[1][2]['label'] and \
                rec_out[0][2]['label'] - edge[2] <lag_cri):
                edge.append( [rec_out[0][1],rec_out[1][1],rec_out[0][2]['label']])
                passNode.append(rec_out[0][1])
                passNode.append(rec_out[0][1])
                passNode.append(edge[1])
                passNode.append(edge[0])
                cont_pair.append(edge)
    return(cont_pair)


def contract_red(G,lag_cri=5):
    """Contract combination transformation.

    Args:
        G (graph) : Reaction graph.
        lag_cri : Lifetime criteria. (default=30)
    Returns:
        cont_pair: Red transformation contaction Suggestion,[(a,b),(e,f), .....]
    """
    #1. Get all node pairs contain blue transformations.
    redList =[] 
    passNode = []
    # [(node1,node2),(node_n,node_n-1),....]
    # 
    RawEdgeList = G.edges(data=True)
    for rec in RawEdgeList:
        if(rec[2]['color'] == 'red' and \
          (rec[0],rec[1])  not in redList and \
           len(G.in_edges(rec[0])) == 2 and \
           len(G.out_edges(rec[0]))== 1 and \
           G.number_of_edges(rec[0],rec[1]) == 1):
            
            if(isinstance(rec[2]['label'], int)):
                redList.append([rec[0], rec[1], rec[2]['label'] ])
            else:
                redList.append([rec[0], rec[1], int(rec[2]['label']) ])
            
    # 1, Check To_node information.
    cont_pair = []
    for edge in redList:
        rec_in = G.in_edges(edge[0],data=True)
        if(rec_in is None or edge[1] in passNode or edge[0] in passNode): continue
        if(len(rec_in) == 2):
            rec_in = list(rec_in)
            if (rec_in[0][2]['color'] == 'red' and \
                rec_in[1][2]['color'] == 'red' and \
                rec_in[0][2]['label'] ==  rec_in[1][2]['label'] and \
                edge[2] - rec_in[0][2]['label'] < lag_cri):
                edge.append( [rec_in[0][0],rec_in[1][0],rec_in[0][2]['label']])
                passNode.append(rec_in[0][0])
                passNode.append(rec_in[1][0])
                passNode.append(edge[1])
                passNode.append(edge[0])
                cont_pair.append(edge)
    return(cont_pair)

def contract_nodes(G):
    """Nodes contraction driver.

    Args:
        G (graph) : Reaction graph.
    Returns:
        G_t (graph) : Reaction Graph after contraction.
    """
    G_t = deepcopy(G)
    # Contract grey.
    N_nodes = len(G_t.nodes(data=True))
    N_nodesp = -1
    passPair = []
    T0 = time.time()
    while (N_nodes-N_nodesp != 0):
        N_nodesp = N_nodes
        contrR,currPass = contract_grey(G_t,passPair,gvar.StableMolLag)
        for cont in contrR:
            G_t = nx.contracted_nodes(G_t,cont[1], cont[0],self_loops=False,copy=False)
        for itm in currPass:
            passPair.append(itm)
        N_nodes = len(G_t.nodes(data=True))   
#       print(N_nodesp,N_nodes,"number of nodes after grey")


    N_nodes = len(G_t.nodes(data=True))
    N_nodesp = -1
    T1 = time.time()

    while (N_nodes-N_nodesp != 0):
        #Contract Blue edges.        
        N_nodesp = N_nodes
        cont_pairs = contract_blue(G_t,gvar.StableMolLag)
        for pair in cont_pairs:
            G_t = nx.contracted_nodes(G_t,pair[0], pair[1],self_loops=False,copy=False)
            for key in G_t[pair[0]][pair[3][0]]:
                if( G_t[pair[0]][pair[3][0]][key]['label'] == pair[3][2]):
                    G_t[pair[0]][pair[3][0]][key]['label'] = pair[2]

            for key in G_t[pair[0]][pair[3][1]]:
                if( G_t[pair[0]][pair[3][1]][key]['label'] == pair[3][2]):
                    G_t[pair[0]][pair[3][1]][key]['label'] = pair[2]

        cont_pairs = contract_red(G_t,gvar.StableMolLag)
        #Contract Red edges.        
        for pair in cont_pairs:
            G_t = nx.contracted_nodes(G_t,pair[1], pair[0],self_loops=False)
            for key in G_t[pair[3][0]][pair[1]]:
                if( G_t[pair[3][0]][pair[1]][key]['label'] == pair[3][2]):
                    G_t[pair[3][0]][pair[1]][key]['label'] = pair[2]

            for key in  G_t[pair[3][1]][pair[1]]:
                if( G_t[pair[3][1]][pair[1]][key]['label'] == pair[3][2]):
                    G_t[pair[3][1]][pair[1]][key]['label'] = pair[2]

        N_nodes = len(G_t.nodes(data=True))   
    T2 = time.time()
    return G_t


def contract_nodes_Onlygrey(G):
    """Nodes contraction driver.

    Args:
        G (graph) : Reaction graph.
    Returns:
        G_t (graph) : Reaction Graph after contraction.
    """
    G_t = deepcopy(G)
    # Contract grey.
    N_nodes = len(G_t.nodes(data=True))
    N_nodesp = -1
    passPair = []
    T0 = time.time()
    while (N_nodes-N_nodesp != 0):
        T00 = time.time()
        N_nodesp = N_nodes
        contrR,currPass = contract_grey(G_t,passPair,gvar.StableMolLag)
        T01 = time.time()
        for cont in contrR:
            G_t = nx.contracted_nodes(G_t,cont[1], cont[0],self_loops=False)
        for itm in currPass:
            passPair.append(itm)
        T02 = time.time()
        
        N_nodes = len(G_t.nodes(data=True))   


    N_nodes = len(G_t.nodes(data=True))
    N_nodesp = -1
    T1 = time.time()

    return G_t


# For test
if __name__ == '__main__':
    #G = nx.drawing.nx_pydot.read_dot("/home/xchen/Develop/jupyter-lab/ChenXin-2021/ReactionNet-analysis/34813.dot")
    G = nx.drawing.nx_pydot.read_dot("/home/xchen/Develop/CPX-MechGen/CPX-MechGen/example/354781.dot")
    G_t = deepcopy(G)
    for edge in  G_t.edges(data=True,keys=True):
        for key in G_t[edge[0]][edge[1]]:
            if(isinstance(G_t[edge[0]][edge[1]][key]['label'], str)):
                G_t[edge[0]][edge[1]][key]['label'] = int(G_t[edge[0]][edge[1]][key]['label'].strip('\"'))
    dic = nx.get_node_attributes(G,'label') 
    N_nodes0 = G.number_of_nodes()
    N_edges0 = G.number_of_edges()
    print("total nodes in original graph: ",N_nodes0)
    print("total edges in original graph: ",N_edges0)
    G_t = contract_nodes(G_t)
    G_t = contract_nodes(G_t)
    N_nodes = G_t.number_of_nodes()
    N_edges = G_t.number_of_edges()
    print("total nodes in contracted graph: ",N_nodes)
    print("total edges in contracted graph: ",N_edges)
    write_dot(G_t, "exp_cont.dot") 




