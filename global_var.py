#++++++++++++++++
#
#  Global variables
#
#++++++++++++++++

atomList = None
# [line[0],[float(line[1]),float(line[2]),float(line[3])],[]]
#   ^                      ^                               ^
#  ele.Name,     element coordinate,                   neighbor

blockList = None
# [line[0],  [list],      []]
#    ^         ^          
#  MolID      AtomList


SpeciesCount = []


#+++++++++++++++++++++
#  Reaction event network graph
#+++++++++++++++++++++

GR = None

#++++++++++++++++++++++
#
#  Dictionaries
#
#++++++++++++++++++++++
ELEMENTS_dict = {
        'GHOST':  0, 'H'    :  1, 'He'   :  2, 'Li'   :  3, 'Be'   :  4,
        'B'    :  5, 'C'    :  6, 'N'    :  7, 'O'    :  8, 'F'    :  9,
        'Ne'   : 10, 'Na'   : 11, 'Mg'   : 12, 'Al'   : 13, 'Si'   : 14,
        'P'    : 15, 'S'    : 16, 'Cl'   : 17, 'Ar'   : 18, 'K'    : 19,
        'Ca'   : 20, 'Sc'   : 21, 'Ti'   : 22, 'V'    : 23, 'Cr'   : 24,
        'Mn'   : 25, 'Fe'   : 26, 'Co'   : 27, 'Ni'   : 28, 'Cu'   : 29,
        'Zn'   : 30, 'Ga'   : 31, 'Ge'   : 32, 'As'   : 33, 'Se'   : 34,
        'Br'   : 35, 'Kr'   : 36, 'Rb'   : 37, 'Sr'   : 38, 'Y'    : 39,
        'Zr'   : 40, 'Nb'   : 41, 'Mo'   : 42, 'Tc'   : 43, 'Ru'   : 44,
        'Rh'   : 45, 'Pd'   : 46, 'Ag'   : 47, 'Cd'   : 48, 'In'   : 49,
        'Sn'   : 50, 'Sb'   : 51, 'Te'   : 52, 'I'    : 53, 'Xe'   : 54,
        'Cs'   : 55, 'Ba'   : 56, 'La'   : 57, 'Ce'   : 58, 'Pr'   : 59,
        'Nd'   : 60, 'Pm'   : 61, 'Sm'   : 62, 'Eu'   : 63, 'Gd'   : 64,
        'Tb'   : 65, 'Dy'   : 66, 'Ho'   : 67, 'Er'   : 68, 'Tm'   : 69,
        'Yb'   : 70, 'Lu'   : 71, 'Hf'   : 72, 'Ta'   : 73, 'W'    : 74,
        'Re'   : 75, 'Os'   : 76, 'Ir'   : 77, 'Pt'   : 78, 'Au'   : 79,
        'Hg'   : 80, 'Tl'   : 81, 'Pb'   : 82, 'Bi'   : 83, 'Po'   : 84,
        'At'   : 85, 'Rn'   : 86, 'Fr'   : 87, 'Ra'   : 88, 'Ac'   : 89,
        'Th'   : 90, 'Pa'   : 91, 'U'    : 92, 'Np'   : 93, 'Pu'   : 94,
        'Am'   : 95, 'Cm'   : 96, 'Bk'   : 97, 'Cf'   : 98, 'Es'   : 99,
        'Fm'   :100, 'Md'   :101, 'No'   :102, 'Lr'   :103, 'Rf'   :104,
        'Db'   :105, 'Sg'   :106, 'Bh'   :107, 'Hs'   :108, 'Mt'   :109,
        'E110' :110, 'E111' :111, 'E112' :112, 'E113' :113, 'E114' :114,
        'E115' :115, 'E116' :116, 'E117' :117, 'E118' :118 }

radii_dict = { 
        'H': 0.31, 'He': 0.28,
        'Li': 1.21, 'Be': 0.96, 'B': 0.84, 'C': 0.69, 'N': 0.71, 
        'O': 0.66, 'F': 0.64, 'Ne': 0.58,
        'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07,
        'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
        'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 'V': 1.53,
        'Cr': 1.39, 'Mn': 1.39, 'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24,
        'Cu': 1.32, 'Zn': 1.22                                                                 
 #       'H':0.4, 'C':0.67, 'N':0.65, 'O':0.57 
}                                                                          
#mass from https://www.periodic-table.org/Radium-atomic-mass/
mass_dict = {  
        'H': 1.008, 'He': 4.003,
        'Li': 6.941, 'Be': 9.012, 'B': 10.811, 'C': 12.011, 'N': 14.007, 
        'O': 15.999, 'F': 18.998, 'Ne': 20.180,
        'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.086, 'P': 30.974,
        'S': 32.065, 'Cl': 35.453, 'Ar': 39.948,
        'K': 39.098, 'Ca': 40.078, 'Sc': 44.956, 'Ti': 47.867, 'V': 50.942,
        'Cr': 51.996, 'Mn': 54.938, 'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693,
        'Cu': 63.546, 'Zn': 65.409
}  

 


radicri_names = {
    'C': {

        'C': 1.34,'H': 1.07,
        'O': 1.27,'Ni': 1.91 
    },

    'H': {
        'C': 1.07,'H': 0.80,
        'O': 0.97,'Ni': 1.64
    },

    'O': {
        'C': 1.27,'H': 0.97,
        'O': 1.14,'Ni': 1.81
    }, 
   'Ni': {
        'C': 1.27,'H': 0.97,
        'O': 1.14,'Ni': 2.48
    }, 
}

#+++++++++
#
# User defined variables
#
#+++++++++

# number of cores
ncores = 4 

trajTyp   = 'lmp' # trajectory type.
trajFiles = []


# LAMMPS elementName to number
# example = {'1':"C", '2':"H", "3":"O"   }
#
#
lmpAtmDict = {} 

# PBC setting
# [(0, 0, 0),(210, 210, 210)]
#
pbcXYZ  = None      # PBC
pbcmole = False     # molecule cut by pbc ?
timePerFrame =  None # default 1fs
stableMolTime = [1,'ps'] #stable molecule lifetime criteria
StableMolLag = 10        # frames to determine lag 

bnd_cri = 1.5


#++++++++++
#  intermediate variable
#
#++++++++++
GlobalMaskMat=None
GlobalDistMat=None

DicStuct = {}
# {HASH: [element],  [coor],    printed?,   SMILE,  formul}
#           ^          ^           ^          ^
#     Element list   Coordinate  True/Fale



DicMoleInfo = {}
#{MolHash: [atomLable]   }
#            ^
#        atomLable-list

DicReactionRec = {}
#{ReactionHash: ["HashLable":
#1.  [HashR,HashP,HashCata.,SMILES_R,SMILES_P,SMILES_cata,(time1,time2)] 
#2.  [HashR,HashP,HashCata.,SMILES_R,SMILES_P,SMILES_cata,(time1,time2)]}
#


#+++++++++++++++++++++
#  Parallel Arrangement
#+++++++++++++++++++++
StepsinTotal = 0
ProcessArrangement = []
#====
#        int          int        int         int
# [
#  (fileNumberStart,startLine,fileNumberEnd,endLine),  Processor 1
#  (fileNumberStart,startLine,fileNumberEnd,endLine)   Processor 2
# ]
#====


StepRecPerCore = []
# [(start_frame,end_frame),(start_frame,end_frame),(start_frame,end_frame)]
#
#


#======
#  Catalysis part
#=====
cata     = False
CataAtom = []
cataSelectFlag = False
cataLabel = 0

#=====
# Time profiling
# 
#=====

Time_DC = 0
Time_MI = 0
Time_TS = 0
Time_FM = 0
Time_BW = 0
Time_ER = 0
Time_NC = 0
