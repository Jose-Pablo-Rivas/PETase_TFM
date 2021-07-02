import yasara, sys, os, time
import numpy as np
import pandas as pd
from itertools import groupby
import re
import subprocess


#### FUNCTIONS ####

def yasaraInitialization():
    yasara.info.mode = "txt"
    yasara.Console("Off")
    yasara.Processors(6, 0)


def distances(protein):
    SER_HB = 'OG res 160'
    ASP_1 = 'OD1 res 206'
    ASP_2 = 'OD2 res 206'
    HIS_1 = 'NE2 res 237'
    HIS_2 = 'ND1 res 237'
    HB_SER_HIS = []
    DIST_SER_HIS = []
    HB_HIS_ASP = []
    DIST_HIS_ASP = []
    FRAME = []
    yasara.LoadYOb(protein)
    FRAME.append(protein)
    # Check if there is a hydrogen bond between the SER and the HIS
    line = "ListHBoAtom  {0}, {1}".format(SER_HB, HIS_1)
    if yasara.run(line) != []:
        HB_SER_HIS.append('YES')
    else:
        HB_SER_HIS.append('NO')
    # Append distance between SER and HIS
    distance_1 = str(yasara.Distance('Name HG res 160', 'Name NE2 res 237'))
    distance_1 = distance_1.replace('[', '')
    distance_1 = distance_1.replace(']', '')
    distance_1 = float(distance_1)
    DIST_SER_HIS.append(round(distance_1, 3))

    # Check if there is a hydrogen bond between the HIS and the ASP
    line_1 = "ListHBoAtom  {0}, {1}".format(HIS_2, ASP_1)
    line_2 = "ListHBoAtom {0}, {1}".format(HIS_2, ASP_2)
    if yasara.run(line_1) != [] or yasara.run(line_2) != []:
        HB_HIS_ASP.append('YES')
    else:
        HB_HIS_ASP.append('NO')
    # Append distance between HIS and ASP
    distance_2 = str(yasara.Distance('Name HD1 res 237', 'Name OD1 res 206'))
    distance_2 = distance_2.replace('[', '')
    distance_2 = distance_2.replace(']', '')
    distance_2 = float(distance_2)
    distance_3 = str(yasara.Distance('Name HD1 res 237', 'Name OD2 res 206'))
    distance_3 = distance_3.replace('[', '')
    distance_3 = distance_3.replace(']', '')
    distance_3 = float(distance_3)
    if distance_2 < distance_3:
        DIST_HIS_ASP.append(round(distance_2, 3))
    elif distance_2 > distance_3:
        DIST_HIS_ASP.append(round(distance_3, 3))
    yasara.DelObj('All')
    return HB_SER_HIS, DIST_SER_HIS, HB_HIS_ASP, DIST_HIS_ASP, FRAME


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
    ##yasara.info.mode = 'txt'
    ##yasara.Console("Off")
    yasara.LoadYOb(prot)
    oxygen = []
    proteins = []
    Dihedral = []
    conformation = []

    for n in range(len(LigandC)):

        # Get distance1 - Ser160 (Catalytic)
        distance_SER = str(yasara.Distance(SER, LigandC[n]))
        distance_SER = distance_SER.replace('[', '')
        distance_SER = distance_SER.replace(']', '')
        distance_SER = float(distance_SER)
        proteins.append(prot)
        prov_SER.append(round(distance_SER, 3))

    distances_SER.append(min(prov_SER))

    # Get Energy:
    line = yasara.run("BFactorRes unk")[0]
    line = float(-line)
    line = round(line, 3)
    energy.append(line)

    # Get distance2 - Met161 (Oxyanion Hole)
    distance_MET = str(yasara.Distance(MET, LigandO[prov_SER.index(min(prov_SER))]))
    distance_MET = distance_MET.replace('[', '')
    distance_MET = distance_MET.replace(']', '')
    distances_MET.append(round(float(distance_MET), 3))

    # Get distance3 - Tyr87 (Oxyanion Hole)
    distance_TYR = str(yasara.Distance(TYR, LigandO[prov_SER.index(min(prov_SER))]))
    distance_TYR = distance_TYR.replace('[', '')
    distance_TYR = distance_TYR.replace(']', '')
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

def dihedrosCalc(protein, residuo):
    # Specify here the atoms for which you want to calculate dihedral angles:
    Atoms_Dihedrals = ['Name CA Res {}'.format(residuo), 'Name CB Res {}'.format(residuo),
                       'Name CG Res {}'.format(residuo), 'Name CD2 Res {}'.format(residuo)]

    Dihedrals = []
    # Load the yob
    yasara.LoadYOb(protein)

    line = 'Dihedral {0}, {1}, {2}, {3}'.format(Atoms_Dihedrals[0], Atoms_Dihedrals[1], Atoms_Dihedrals[2],
                                                Atoms_Dihedrals[3])
    dihedral = yasara.run(line)[0]
    Dihedrals.append(round(float(dihedral), 3))

    ## remove yasara objects

    yasara.DelObj('All')

    return Dihedrals

#Specify here the atoms for which you want to calculate the distance:
LigandC = ['Name C2 Res UNK','Name C7 Res UNK','Name C12 Res UNK','Name C17 Res UNK']
LigandO = ['Name O3 Res UNK', 'Name O4 Res UNK', 'Name O7 Res UNK', 'Name O8 Res UNK']

#Specify here the atoms for which you want to calculate the dihedral angles:
LigandD = ['Name O2 Res UNK', 'Name C8 Res UNK', 'Name C11 Res UNK', 'Name O5 Res UNK']


#### MAIN ####


start = time.time()
yasaraInitialization()
proteins = sys.argv[1:]
HSH = []
DIS_1 = []
HHA = []
DIS_2 = []
F = []
n=0
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
DIHEDRAL_LIST = ['159', '185', '237']
DIHEDRAL_159 = []
DIHEDRAL_185 = []
DIHEDRAL_237 = []

for protein in proteins:
    if n < 10001:
        HB_SER_HIS, DIST_SER_HIS, HB_HIS_ASP, DIST_HIS_ASP, FRAME = distances(protein)
        HSH = HSH + HB_SER_HIS
        DIS_1 = DIS_1 + DIST_SER_HIS
        HHA = HHA + HB_HIS_ASP
        DIS_2 = DIS_2 + DIST_HIS_ASP
        F = F + FRAME
        Serine = "160"
        Metionine = '161'
        Tyrosine = '87'
        # Tyr_Met_oxyanion = hydrogen_bond_oxyanion(prot, LigandO)
        # if Tyr_Met_oxyanion != []:
        distances_SER, Met, distances_MET, Tyr, distances_TYR, energy, proteins, Dihedral, conformation = getDistances(
            protein, LigandC, LigandO, LigandD)
        DIHEDRAL = DIHEDRAL + Dihedral[0]
        CONF = CONF + conformation
        DIST_SER = DIST_SER + distances_SER
        MET = MET + Met
        DIST_MET = DIST_MET + distances_MET
        TYR = TYR + Tyr
        DIST_TYR = DIST_TYR + distances_TYR
        ENE = ENE + energy
        N.append(n)
        for i in range(len(DIHEDRAL_LIST)):
            if DIHEDRAL_LIST[i] == '159':
                DIHEDRAL_159.append(dihedrosCalc(protein, DIHEDRAL_LIST[i])[0])
            elif DIHEDRAL_LIST[i] == '185':
                DIHEDRAL_185.append(dihedrosCalc(protein, DIHEDRAL_LIST[i])[0])
            elif DIHEDRAL_LIST[i] == '237':
                DIHEDRAL_237.append(dihedrosCalc(protein, DIHEDRAL_LIST[i])[0])
        n+=1
    else:
        break

print(len(HSH))
print(len(HHA))
print(len(F))
df = pd.DataFrame(list(zip(N, HSH, DIS_1, HHA, DIS_2, DIST_SER, ENE, MET, DIST_MET, TYR, DIST_TYR, DIHEDRAL, CONF, DIHEDRAL_159, DIHEDRAL_185, DIHEDRAL_237, F)),
                  columns =['FRAME', 'HB_SER_HIS', 'DISTANCE_SER_HIS', 'HB_HIS_ASP', 'DISTANCE_HIS_ASP', 'DISTANCE_SER', 'ENERGY', 'OX_HB_MET', 'DISTANCE_MET', 'OX_HB_TYR', 'DISTANCE_TYR', 'DIHEDRAL', 'CONFORMATION', 'DIHEDRAL_159', 'DIHEDRAL_185', 'DIHEDRAL_237', 'STRUCTURE'])
df = df.sort_values(by=['FRAME'])
df.to_csv('PETase_PET_303.csv', index=False)
print("time (min): ", (time.time() - start)/60.)
yasara.Exit()
