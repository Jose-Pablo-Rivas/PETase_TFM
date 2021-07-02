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

def catalytic(protein):
    SER_HB = 'OG res 160'
    ASP_1 = 'OD1 res 206'
    ASP_2 = 'OD2 res 206'
    HIS_1 = 'NE2 res 237'
    HIS_2 = 'ND1 res 237'
    CATALYTIC=[]
    yasara.LoadPDB(protein)
    # Check if there is a hydrogen bond between the SER and the HIS
    line_1 = "ListHBoAtom  {0}, {1}".format(SER_HB, HIS_1)
    line_2 = "ListHBoAtom  {0}, {1}".format(HIS_2, ASP_1)
    line_3 = "ListHBoAtom {0}, {1}".format(HIS_2, ASP_2)
    if yasara.run(line_1) != [] and (yasara.run(line_2) != [] or yasara.run(line_3) != []):
        CATALYTIC.append('YES')
    else:
        CATALYTIC.append('NO')
    yasara.DelObj('All')
    return CATALYTIC


##########
## MAIN ##
##########

directory = sys.argv[1]

cluster=[]
representatives=[]
n_members=[]
catalytic_structure=[]
system=[]
Temperature=[]

for filename in os.listdir(directory):
    if filename.endswith('.pdb'):
        CAT = catalytic(filename)
        catalytic_structure.append(CAT)
        n = [float(n) for n in re.findall(r'-?\d+\.?d*', filename)]
        cluster.append(int(n[0]))
        representatives.append(int(n[1]))
        n_members.append(int(n[2]))
        system.append(directory.split('/')[2])
        Temperature.append(directory.split('/')[1].split('_')[1])
df = pd.DataFrame(list(zip(cluster,representatives,n_members, catalytic_structure, system, Temperature)), columns=['Cluster', 'Structure_Representative', 'Number_Members', 'Is_Catalytic', 'Variant', 'Temperature'])
df = df.sort_values(by=['FRAME'], ascending=False)
df.to_csv('Cluster_{0}_{1}.csv'.format(system[0], Temperature[0]), index=False)
print("time (min): ", (time.time() - start)/60.)
yasara.Exit()