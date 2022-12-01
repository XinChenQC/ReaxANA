from tool_contract import *
from tool_removeOscill import *

from timeit import default_timer as timer
from datetime import timedelta


def prtNodeInfo(G,inode):
    rec_in = G.in_edges(inode,data=True)
    rec_out = G.out_edges(inode,data=True)
    print(inode)
    print(" inRec:")
    for rec1 in rec_in:
        print(rec1)
    print("outRec:")
    for rec2 in rec_out:
        print(rec2)
    print("================")

if __name__ == '__main__':
    t0 = timer()  # Timer 0
    G = nx.drawing.nx_pydot.read_dot("reduce_small56.dot")
#   G = nx.drawing.nx_pydot.read_dot("./reduce_small56_uc.dot")
    dic = nx.get_node_attributes(G,'label') 
    G_t = deepcopy(G)
    for edge in  G_t.edges(data=True,keys=True):
        for key in G_t[edge[0]][edge[1]]:
            if(isinstance(G_t[edge[0]][edge[1]][key]['label'], str)):
                G_t[edge[0]][edge[1]][key]['label'] = int(G_t[edge[0]][edge[1]][key]['label'].strip('\"'))
    t1 = timer()  # Timer 1
    #prtNodeInfo(G_t,'7402809223411185418')

    print("Loadin time:  "+str(timedelta(seconds=t1-t0)))
    N_nodes = G.number_of_nodes()
    N_edges = G.number_of_edges()
    print("total nodes in original graph: ",N_nodes)
    print("total edges in original graph: ",N_edges)
    N_nodesp = -1
    s1 = timer()  # Timer 1
    while(N_nodes - N_nodesp != 0):
        N_nodesp = N_nodes
        G_t = remove_useless_trans(G_t)
        N_nodes = G_t.number_of_nodes()
        N_edges = G_t.number_of_edges()
        print("total nodes after removal: ",N_nodes)
        print("total edges after removal: ",N_edges)
        #prtNodeInfo(G_t,'7402809223411185418')
    s2 = timer()  # Timer 1
    print("Clean time:  "+str(timedelta(seconds=s2-s1)))
    write_dot(G_t, "reduce_small56_r.dot") 


    N_nodesp = -1
    while(N_nodes - N_nodesp != 0):
        G_t = contract_nodes(G_t)
        N_nodesp = N_nodes
        N_nodes = G_t.number_of_nodes()
        N_edges = G_t.number_of_edges()
        G_t = nodeShortRemove(G_t,100)
        print("total nodes after GBR: ",N_nodes)
        print("total edges after GBR: ",N_edges)
        #prtNodeInfo(G_t,'7402809223411185418')
    s3 = timer()  # Timer 1
    print("Contract time:  "+str(timedelta(seconds=s3-s2)))
    write_dot(G_t, "reduce_small56_cf.dot")
#   G_t = contract_nodes(G_t)
#   N_nodes = G_t.number_of_nodes()
#   N_edges = G_t.number_of_edges()
#   print("total nodes in original graph: ",N_nodes)
#   print("total edges in original graph: ",N_edges)
#   write_dot(G_t, "reduce_small56.dot") 




