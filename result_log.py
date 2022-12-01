import networkx as nx
import pydot
from copy import deepcopy
import global_var as gvar
import sys,os
import json


def print_XYZandSVG(key):
    from rdkit import Chem
    from rdkit.Chem import Draw
    DicStuct = gvar.DicStuct
    SMILES = DicStuct[key][3]
    if(len(SMILES)==0): SMILES="FAILTOPARSE" 
    # Print XYZ File
    FileName = "S_"+ DicStuct[key][4] +"_"+ str(key)[0:7]+ ".xyz"
    f = open("./specRec/"+FileName,"w")
    f.write( str(len(gvar.DicStuct[key][0])) + "\n")
    f.write( SMILES+","+gvar.DicStuct[key][4] + "\n")
    for i in range(len(gvar.DicStuct[key][0])):
        f.write(str(gvar.DicStuct[key][0][i])+"     "+\
                str(gvar.DicStuct[key][1][i][0])+"  "+\
                str(gvar.DicStuct[key][1][i][1])+"  "+\
                str(gvar.DicStuct[key][1][i][2])+"\n")
    f.close()
    # Print SVG File
    if(len(SMILES)==0 or SMILES=="FAILTOPARSE"): return
    if("CC_" == SMILES[0:3] or "NC_"== SMILES[0:3]):  # For catalysis
        mol = Chem.MolFromSmiles(SMILES[3:])
    else:
        mol = Chem.MolFromSmiles(SMILES)
    FileName = "S_"+ DicStuct[key][4] +"_"+ str(key)[0:7]+ ".svg"
    try:
        Draw.MolToFile(mol,"./specRec/"+FileName,size=(200, 200),fitImage=True)
    except:
         print("Save failed: ",SMILES)

def PrintMetaDataRaw():
    fileResult = open("result.log",'w') 
    G=gvar.GR
    DicreacR = gvar.DicReactionRec
    specDict = []
    TimeSeriedDict = {}
    MaxSmileLength = 0
    for node in G.nodes(data=True):
        hashD = node[1]['hashD']
        if("int-" not in hashD and len(node[1]['SMILE'])>0):
            specDict.append(hashD)
        # Print Failed to solve structrues.
        if(len(node[1]['SMILE'])==0 ): print_XYZandSVG(hashD)
    specDict = list(set(specDict))

    for spec in specDict:
        specRev = []
        SMILES_label = gvar.DicStuct[spec][3]
            
        for step in gvar.SpeciesCount:
            if(spec in step):
                specRev.append(step[spec])
            else:
                specRev.append(0)
        TimeSeriedDict[spec] = [gvar.DicStuct[spec][3],specRev]
    # Record species revolution
    smileRev={} 
    for key in TimeSeriedDict:
        if (TimeSeriedDict[key][0] in smileRev):
            smileRev[TimeSeriedDict[key][0]] = \
            [ smileRev[TimeSeriedDict[key][0]][i] + TimeSeriedDict[key][1][i] \
            for i in range(len(smileRev[TimeSeriedDict[key][0]]))]
        else:
            smileRev[TimeSeriedDict[key][0]] = TimeSeriedDict[key][1]

    # Collect species infomation according to SMILES
    SMILES_INFO = {}
    for spec in specDict:
        SMILES = gvar.DicStuct[spec][3]
        if (len(SMILES) > MaxSmileLength):
            MaxSmileLength = len(SMILES)
        if(SMILES not in SMILES_INFO):
            abandnt  = max(smileRev[SMILES]) # Abundant of species
            fName    = "S_"+ gvar.DicStuct[spec][4] +"_"+ str(spec)[0:7]
            formul   = gvar.DicStuct[spec][4]
            SMILES_INFO[SMILES] = [abandnt,fName,formul,"1"]
            print_XYZandSVG(spec)

    SMILES_INFO = dict(sorted(SMILES_INFO.items(), key=lambda item: item[1][0],reverse=True))
    # Print species information
    fileResult.write("[species]\n")
    fileResult.write('#----------------------------------------------------------------------------------------\n')
    fileResult.write(f'{"#SMILES":{MaxSmileLength+2}} {"filename":24} {"Formula":12} {"Abundance":7}\n')
    fileResult.write('#----------------------------------------------------------------------------------------\n')
    for SMILES in SMILES_INFO:
        fileResult.write(f'{SMILES:{MaxSmileLength+2}}{SMILES_INFO[SMILES][1]:24}{SMILES_INFO[SMILES][2]:12}{SMILES_INFO[SMILES][0]:7}\n')
        #                    ^SMILES ^filename ^Formula ^Maximum Abundance                          
    fileResult.write("[species_end]\n\n")
    # Print specie Revolution
    fileResult.write("[speciesrev]\n")
    fileResult.write(f'{len(smileRev):5}{len(smileRev[SMILES]):10}\n')
    Totstep = len(smileRev[SMILES])
    for key in SMILES_INFO:
        #  print Rec 10 itm per line
        i = 1
        SMILES_INFO[key][3] = smileRev[key]
        fileResult.write(f'{key:{MaxSmileLength+2}}\n')
        for itm in smileRev[key]:
            if(i%17 == 0):
                fileResult.write(f'{itm:5}\n')
            else:
                fileResult.write(f'{itm:5}')
            i += 1
        fileResult.write('\n')
    fileResult.write("[speciesrev_end]\n")


    # Print formular Revolution
    formula_INFO = {} # Formula: [SMILES1,SMILES2,... ]
    for key in SMILES_INFO:
        if(SMILES_INFO[key][2] not in formula_INFO ):
            formula_INFO[SMILES_INFO[key][2]] = [key]
        else:
            formula_INFO[SMILES_INFO[key][2]].append(key)

    formula_REC = {} # Formula: [number ]
    for key in formula_INFO:
        formula_REC[key] = [ 0 for i in range(Totstep)]
        for SMILES in formula_INFO[key]:
            formula_REC[key] = \
            [smileRev[SMILES][i]+formula_REC[key][i] for i in range(Totstep)]

    fileResult.write("[formularev]\n")
    for key in formula_REC:
        #  print Rec 10 itm per line
        i = 1
        fileResult.write(f'{key:{MaxSmileLength+2}}\n')
        for itm in formula_REC[key]:
            if(i%17 == 0):
                fileResult.write(f'{itm:5}\n')
            else:
                fileResult.write(f'{itm:5}')
            i += 1
        fileResult.write('\n')

    fileResult.write("[formularev_end]\n")
    fileResult.write("[Reactionrev]\n")
    iReac = 0
    MOLE_rec ={}
    for k in DicreacR:
        # print reaction title
        iReac += 1
        RectionStr = " #"+str(iReac)
        for SML_reatant in DicreacR[k][0][4]:
            RectionStr =  RectionStr +" " +SML_reatant+" +"
        RectionStr = RectionStr.rstrip("+")

        RectionStr = RectionStr +" --> "

        for SML_prod in DicreacR[k][0][5]:
            RectionStr =  RectionStr +" " +SML_prod+" +"
        RectionStr = RectionStr.rstrip("+")

        if(len(DicreacR[k][0][6]) != 0 ):
            RectionStr = RectionStr + "   cata. "+DicreacR[k][0][6][0]
        RectionStr = RectionStr +" # of events: "+str(len(DicreacR[k]))+"\n"
        fileResult.write(RectionStr)
        if('int' in RectionStr): continue
        # print reaction event 1 by 1:
        ievent = 0
        for event in DicreacR[k]:
            EventStr = "   "+str(ievent)+":  "
            ievent = 1+ievent

            for Rhash in event[1]:
                EventStr = EventStr +"  " +str(gvar.DicMoleInfo[Rhash])+ "  +"
                MOLE_rec[Rhash] =gvar.DicMoleInfo[Rhash]
            EventStr = EventStr.rstrip("+")
            EventStr = EventStr + " --> "

            for Phash in event[2]:
                EventStr = EventStr +"  " +str(gvar.DicMoleInfo[Phash])+ "  +"
                MOLE_rec[Phash] =gvar.DicMoleInfo[Phash]
            EventStr = EventStr.rstrip("+")

            if(len(event[3])!=0):
                EventStr = EventStr +"  cata. "+str(gvar.DicMoleInfo[event[3][0]])
                MOLE_rec[event[3][0]] = gvar.DicMoleInfo[event[3][0]]
            fileResult.write(EventStr+"\n")

    fileResult.write("[Reactionrev_end]\n")

    # Render JSON files.
    FileJSON = open("data.json",'w',encoding='utf-8')
    JSONDic={"SMILES"  : SMILES_INFO,
             "FORMULA" : formula_REC,
             "REACTION": DicreacR,
             "MOLECULE": MOLE_rec} # SMILES, formula,reactions
    json.dump(JSONDic,FileJSON,ensure_ascii=False,indent=None)

    print(len(specDict))
    fileResult.close()
            

if __name__ == '__main__':
    G = nx.drawing.nx_pydot.read_dot("./Graph_5000.dot")
    get_importantSpec(G)
