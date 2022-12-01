import global_var as gvar
from result_log import *


def recoverData():
    import networkx as nx
    import pickle
    gvar.DicReactionRec= pickle.load(open( "DicReactRec1.p", "rb" ))
    gvar.DicStuct      = pickle.load(open( "DicStruct1.p",   "rb" ))
    gvar.DicMoleInfo   = pickle.load(open( "DicMoleinfo1.p", "rb" ))
    gvar.SpeciesCount  = pickle.load(open( "SpecieCount1.p", "rb" ))
    gvar.GR = nx.drawing.nx_pydot.read_dot("G_5000.dot")

    EmptyList = []
    for node in gvar.GR.nodes(data=True):
        if(len(node[1])==0): EmptyList.append(node[0])
    for node in EmptyList:
        gvar.GR.remove_node(node)

    for node in gvar.GR.nodes(data=True):
        print(node[1])
        hashD = node[1]['hashD']
    nx.drawing.nx_pydot.write_dot(gvar.GR,"G_5000_t.dot")



if __name__ == "__main__":
    recoverData()
    PrintMetaDataRaw()

