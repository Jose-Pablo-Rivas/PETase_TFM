from itertools import groupby
import sys,os,re,yasara, time
import pandas as pd
import numpy as np
import subprocess

#### Functions ####

def yasaraInitialization():
    yasara.info.mode = 'txt'
    yasara.Console("Off")
    yasara.Processors(6, 0)


def getDistances(prot, LigandC, LigandO, LigandD):
    SER = 'Name OG Res 160'
    MET = 'Name N Res 161'
    TYR = 'Name N Res 87'
    distances_SER = []
    prov_SER = []
    Met = []
    distances_MET = []
    Tyr = []
    distances_TYR = []
    energy = []
    #Load the Yob file
    yasara.LoadYOb(prot)
    oxygen = []
    proteins = []
    Dihedral = []
    conformation = []
        
    for n in range(len(LigandC)):
        
        #Get distance1 - Ser160 (Catalytic) 
        distance_SER = str(yasara.Distance(SER, LigandC[n]))
        distance_SER = distance_SER.replace('[','')
        distance_SER = distance_SER.replace(']','')
        distance_SER = float(distance_SER)
        distance_SER = float(distance_SER)
        proteins.append(prot)
        prov_SER.append(round(distance_SER, 3))
        
    #Get Energy:
    line = yasara.run("BFactorRes unk")[0]
    line = float(-line)
    line = round(line, 3)
    energy.append(line)
        
    #Get distance2 - Met161 (Oxyanion Hole)
    oxygen.append(n)
    distance_MET = str(yasara.Distance(MET, LigandO[n]))
    distance_MET = distance_MET.replace('[','')
    distance_MET = distance_MET.replace(']','')
    distances_MET.append(round(float(distance_MET), 3))
        
    #Get distance3 - Tyr87 (Oxyanion Hole)
    distance_TYR = str(yasara.Distance(TYR, LigandO[n]))
    distance_TYR = distance_TYR.replace('[','')
    distance_TYR = distance_TYR.replace(']','')
    distances_TYR.append(round(float(distance_TYR), 3))

    # Get dihedrals
    dihedral = get_dihedrals(LigandD)
    Dihedral.append(dihedral)
    conformation.append(get_conformation(dihedral))
    # Get hydrogen bonds
    line = 'ListHBoAtom  %s, Name N Res 87, results = 6' % LigandO[prov_SER.index(min(prov_SER))]
    if yasara.run(line) != []:
        # ph_tyr.append(yasara.run(line)[:6])
        Tyr.append('YES')
    else:
        Tyr.append('NO')
    line = 'ListHBoAtom  %s, Name N Res 161, results = 6' % LigandO[prov_SER.index(min(prov_SER))]
    if yasara.run(line) != []:
        # ph_met.append(yasara.run(line)[:6])
        Met.append('YES')
    else:
         Met.append('NO')

    yasara.DelObj('All')
    return distances_SER, Met, distances_MET, Tyr, distances_TYR, energy, proteins, Dihedral, conformation


def get_dihedrals(LigandD):
    Dihedrals=[]
    #inialize counters
    line = 'Dihedral {0},{1},{2},{3}'.format(LigandD[0], LigandD[1], LigandD[2], LigandD[3])
    dihedral = yasara.run(line)[0]
    Dihedrals.append(abs(round(float(dihedral), 3)))
    return Dihedrals


def get_conformation(dihedral):
    if 50 <= dihedral[0] <= 90:
        return 'GAUCHE'
    elif 150 <= dihedral[0] <= 180:
        return 'TRANS'
    else:
        return 'UNDEFINED'

def get_conformation(dihedral):
    if 50 <= dihedral <= 90:
        return 'gauche'
    elif 150 <= dihedral <= 180:
        return 'trans'
    else:
        return 'undefined'

def get_variant(prot):
    return prot.split('_')[1]

def get_temperature(prot):
    return prot.split('_')[2]

# Get the command line arguments:


#Specify here the atoms for which you want to calculate the distance:
LigandC = ['Name C2 Res UNK','Name C7 Res UNK','Name C12 Res UNK','Name C17 Res UNK']
LigandO = ['Name O3 Res UNK', 'Name O4 Res UNK', 'Name O7 Res UNK', 'Name O8 Res UNK']

#Specify here the atoms for which you want to calculate the dihedral angles:
LigandD = ['Name O2 Res UNK', 'Name C8 Res UNK', 'Name C11 Res UNK', 'Name O5 Res UNK']


#### MAIN ####


start = time.time()
yasaraInitialization()
proteins = sys.argv[1:]
DIST_SER = []
MET = []
DIST_MET = []
TYR = []
DIST_TYR = []
ENE = []
STRUCT = []
N = []
DIHEDRAL = []
CONF = []
VAR = []
TEMP = []

##Main functions

for prot in proteins:
    Serine = "160"
    Metionine = '161'
    Tyrosine = '87'
    #Tyr_Met_oxyanion = hydrogen_bond_oxyanion(prot, LigandO)
    #if Tyr_Met_oxyanion != []:
    distances_SER, Met, distances_MET, Tyr, distances_TYR, energy, proteins, Dihedral, conformation = getDistances(
        prot, LigandC, LigandO, LigandD)
    DIHEDRAL = DIHEDRAL + Dihedral[0]
    CONF = CONF + conformation
    DIST_SER = DIST_SER + distances_SER
    MET = MET + Met
    DIST_MET = DIST_MET + distances_MET
    TYR = TYR + Tyr
    DIST_TYR = DIST_TYR + distances_TYR
    ENE = ENE + energy
    #Get the variant
    variant=get_variant(prot)
    VAR = VAR + variant
    #Get the temperature
    temperature=get_temperature(prot)
    TEMP = TEMP + temperature
    PROT = PROT + prot
df = pd.DataFrame(list(zip(PROT, DIST_SER, ENE, MET, DIST_MET, TYR, DIST_TYR, DIHEDRAL, CONF, VAR, TEMP)), columns =['STRUCTURE', 'DISTANCE_SER','ENERGY', 'MET', 'DISTANCE_MET', 'TYR', 'DISTANCE_TYR', 'DIHEDRAL', 'CONFORMATION', 'VARIANT', 'TEMPERATURE'])
df = df.sort_values(by=['DISTANCE_SER'])
df.to_csv('DOCKING_METRICS_ANALYSIS.csv', index=False)
