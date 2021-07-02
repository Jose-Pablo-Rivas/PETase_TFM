import sys
import os
from Bio.PDB import *
import ast

# Vamos a cragar los dos archivos pdb. El primero que debemos poner como argumento es el WT.pd y segundo el DuraPETase.pdb

WT = sys.argv[1]
DP = sys.argv[2]

# Vamos a sacar la información relevante de nuestros archivos (los b-factors)
def bfactors_dict(WT, DP):
    bfact_WT = {}
    bfact_DP = {}
    p = PDBParser()
    q = PDBParser()
    struc_WT = p.get_structure('WT', WT)
    struc_DP = q.get_structure('DP', DP)
    c1 = 32
    c2 = 32
    # itero sobre la estructura del WT
    for model in struc_WT:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    # obtengo los bfactors de esos átomos
                    # a.get_bfactor()
                    if residue.get_resname()+str(c1) not in bfact_WT:
                        bfact_WT[residue.get_resname()+str(c1)] = [[atom.get_bfactor()]]
                    else:
                        bfact_WT[residue.get_resname()+str(c1)].append([atom.get_bfactor()])
                c1 += 1
    # itero sobre la estructura del DP
    for model in struc_DP:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    # obtengo los bfactors de esos átomos
                    # a.get_bfactor()
                    if residue.get_resname()+str(c2) not in bfact_DP:
                        bfact_DP[residue.get_resname()+str(c2)] = [[atom.get_bfactor()]]
                    else:
                        bfact_DP[residue.get_resname()+str(c2)].append([atom.get_bfactor()])
                c2 += 1
    #Hemos creado los diccionarios con key=Res** y como valores los bfactors de los átomos de ese residuo.

    return bfact_WT, bfact_DP

def bfact_dif_common_res(dict_WT, dict_DP):

    #Esta función va a sacar un diccionario común que contenga los residuos compartidos por ambas proteínas
    #y la diferencia en sus bfactors. Los no compartidos llevarán el nombre del residuo del WT y valores de bfactor=0

    dict_common={}

    #En lugar de seguir la misma nomenclatura que antes para las keys, para a utilizar solo el númeor del residuo

    count = 32

    #Vamos a iterar a través a las keys de los diccionarios para ver keys comunes y poder añadirlas al diccionario común

    for key in dict_WT:
        dict_common[count] = []
        if key in dict_DP:
            for i in range(len(dict_WT[key])):
                dict_common[count].append([round(dict_DP[key][i][0] - dict_WT[key][i][0], 2)])
        else:
            for _ in range(len(dict_WT[key])):
                dict_common[count].append([0.00])
        count+=1

    #vamos a decirle a función que nos devuelva el diccionario común
    return dict_common


def pdbBfactor(pdb, data_dict):
    pdb_io = PDBIO()
    pdb_parser = PDBParser()
    pdb_file= pdb
    structure = pdb_parser.get_structure(" ", pdb_file)
    count=32
    for model in structure:
        for chain in model:
            for residue in chain:
                for i, atom in enumerate(residue.get_atoms()):
                    atom.bfactor = data_dict[count][i][0]
                count+=1

    pdb_io.set_structure(structure)
    pdb_io.save(pdb_file + "-new.pdb")
'''
        if line[0:6] == "ATOM   " or line[0:7] == "HETATM ":
            resnum = line[23:26]strip()
            if resnum in list(data_dict.keys()):
                out.append("%s%6s%s.2F%s" % (line[:60]))
'''

# Esta es el cuerpo principal del script

bfact_wt, bfact_dp = bfactors_dict(WT, DP)
bf_common = bfact_dif_common_res(bfact_wt, bfact_dp)
print(bfact_wt)
pdbBfactor(WT, bf_common)

# vamos a crear el fichero txt con los nuevos valores de bfactor
'''
f = open('bfactors.txt', 'wt')
f.writelines(''.join(str(j) + '\n' for j in bfactor_list))
f.close()
'''