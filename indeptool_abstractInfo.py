import copy
import sys
filename  = sys.argv[1]

spctag    = "[species]"
spcRevtag = "[speciesrev]"
fmlRevtag = "[formularev]"
TrajLen = 0
N_smiles = 0

R_Flag =""

spcDic = {}
fmlDic = {}

def prtspcDicOrigin():
    time = [ i*0.2 for i in range(TrajLen)]
    # print species record:
    for rec in spcDic:
        fileName = "S_"+rec+"_spec.data"
        fileData = open(fileName,"w")
        for i in range(TrajLen):
            pass
            fileData.write(f'  {time[i]:.3f}{spcDic[rec][2][i]:10d}\n')
        fileData.close()
    

    # print formula record:
    for rec in fmlDic:
        fileName = "F_"+rec+"_formula.data"
        fileData = open(fileName,"w")
        for i in range(TrajLen):
            fileData.write(f'  {time[i]:.3f}{fmlDic[rec][i]:10d}\n' )
        fileData.close()
# Read result file
with open(filename, "r") as f:
    while True:
        line = f.readline()
        if not line:
            break
        line = line.strip()

        if(len(line) == 0 ): continue
        if(line[0]   =="#"): continue

        if("[species]"    in line): 
            R_Flag = "SPCINFO"
            continue

        if("[speciesrev]" in line):
            R_Flag = "SPCREV"
            continue

        if("[formularev]" in line):
            R_Flag = "FMLREV"
            continue

        if("_end]"        in line): 
            R_Flag = ""
            continue

        # Read species information in
        # SMILES:[Formula,abd,[rev]]
        if(R_Flag == "SPCINFO"):
            Sline = line.split()
            smiles    = Sline[0].strip()
            formula   = Sline[2]
            Abundance = Sline[3]
            spcDic[smiles] = [formula,Abundance]

        # Read species revolution information in
        if(R_Flag == "SPCREV"):
            Sline = line.split()
            N_smiles = int(Sline[0])
            TrajLen  = int(Sline[1])
            R_Flag = "SMILEreading_SPCREV"
            continue

        if(R_Flag == "SMILEreading_SPCREV"):
            SMILES = line
            R_Flag = "Valuereading_SPCREV"
            icount = 0
            value = []
            continue
        # Read value
        if(R_Flag == "Valuereading_SPCREV"):
            Sline = line.split()
            for num in (Sline):
                icount += 1
                value.append(int(num))
            if(icount == TrajLen ):
                icount = 0
                R_Flag = "SMILEreading_SPCREV"
                spcDic[SMILES].append(value)
                continue

        # Read formula revolution information in
        if(R_Flag == "FMLREV"):
            Sline = line.split()
            N_formula = int(Sline[0])
            TrajLen  = int(Sline[1])
            R_Flag = "FORMULAreading_FMLREV"
            continue

        if(R_Flag == "FORMULAreading_FMLREV"):
            formula = line
            R_Flag = "Valuereading_FMLREV"
            icount = 0
            value = []
            continue

        # Read value
        if(R_Flag == "Valuereading_FMLREV"):
            Sline = line.split()
            for num in (Sline):
                icount += 1
                value.append(int(num))
            if(icount == TrajLen ):
                icount = 0
                R_Flag = "FORMULAreading_FMLREV"
                fmlDic[formula] = value
                continue

# print for origin.

prtspcDicOrigin()



