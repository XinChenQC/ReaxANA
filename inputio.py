import global_var as gvar
import sys

'''
 Get setting from import file
 ncores:           number of cores. (integer)
 trajtype:         trajectory type. (string, xyz or lammps)
 molelifetime:     lifetime criteria for a stable molecule (float)
 pbc:              is the molecule cut by pbc boundary (boolean, True, False)
 pbcbox:           Define the size of box.  (float, x0 x1    y0 y1    z0 z1)
 time:             time per frame (float string, 1 fs)
'''

def renderSetting(Inpfile):
    inpFile = open(Inpfile,'r')
    lines = inpFile.readlines()
    
    # Read settings by tags
    Flagncore    = False
    Flagpbc      = False
    Flagtrjtyp   = False
    Flaglftime   = False
    catalystAtom = []
    for line in lines:
        lineR = line.strip()
        if(len(lineR) == 0): continue
        if(lineR[0] == "#"): continue
        optionLine = lineR.split()

        if(optionLine[0].lower()=='ncores'):     gvar.ncores   = int(optionLine[1]);Flagncore=True
        if(optionLine[0].lower()=='trajtype'):   gvar.trajTyp  = optionLine[1];Flagtrjtyp=True
        if(optionLine[0].lower()=='trajfiles'):  gvar.trajFiles= optionLine[1:] 
        if(optionLine[0].lower()=='molelifttime'):  
            gvar.stableMolTime[0] = float(optionLine[1])
            gvar.stableMolTime[1] = optionLine[2].lower()
            Flaglftime   = True

        if(optionLine[0].lower()=='timeperframe'):
            gvar.timePerFrame = [0,0]
            gvar.timePerFrame[0] = float(optionLine[1])
            gvar.timePerFrame[1] = optionLine[2].lower()

        if(optionLine[0].lower()=='pbcmole' and "t" in optionLine[1].lower()): 
            gvar.pbcmole   = True 
            Flagpbc   = True

        if(optionLine[0].lower()=='pbcmole' and "f" in optionLine[1].lower()): 
            gvar.pbcmole   = False
            Flagpbc   = True

        if(optionLine[0].lower()=='pbcbox'): 
           #[(x0, y0, z0),(x1, y1, z1)]
            gvar.pbcXYZ = [(float(optionLine[1]), float(optionLine[3]), float(optionLine[5]) ),\
                           (float(optionLine[2]), float(optionLine[4]), float(optionLine[6]) )]
        if(optionLine[0].lower()=='lmpelement'):  
            for i in range(len(optionLine)-1):
                gvar.lmpAtmDict[i+1] = optionLine[i+1]

        # For analysing systems with catalysts:
        if(optionLine[0].lower()=='catalysis' and "t" in optionLine[1].lower()): 
            gvar.cata   = True

        if(optionLine[0].lower()=='selectype'): 
            cataSelectTyp   = int(optionLine[1])

        if(optionLine[0].lower()=='cataselection'): 
            cataAtomselection = optionLine[1]
    
    # Check the consistency of the setting
    if(len(gvar.trajFiles) == 0):
        sys.exit("Error: you have to provide filename for the trajectories")

    if(gvar.timePerFrame is None):
        sys.exit("Error: you have to provide time per frame")

    if(gvar.trajTyp == 'xyz' and gvar.pbcXYZ is None):
        sys.exit("Error: you have to provide pbcbox boundary for xyz-format file")
    # have catalysis
    if(gvar.cata):
        if(cataSelectTyp == 1):
            gvar.cataSelectFlag = True
            gvar.cataLabel = str(cataAtomselection)
        if(cataSelectTyp == 2):
            if('-' in cataAtomselection):
                a1 = int(cataAtomselection.split('-')[0])-1
                a2 = int(cataAtomselection.split('-')[1])
                gvar.CataAtom = list(range(a1,a2))
            if(',' in cataAtomselection):
                lableList = cataAtomselection.split(',')
                gvar.CataAtom = [int(x)-1 for x in lableList]


    # Default setting warning for User
    if(not Flagncore):   print("Warning: The CPU cores number is set to 1")
    if(not Flagpbc):     print("Warning: Molecules can be break by PBC boundary")
    if(not Flaglftime):  print("Warning: Molecular lifetime criterion is set 1ps")
    if(not Flagtrjtyp):  print("Warning: Trajectry format is set to LAMMPS .trj format")
