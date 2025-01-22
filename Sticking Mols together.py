# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 12:03:44 2024

@author: Saul Reid-Guest
Sticking molecules together
"""

import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


m = Chem.MolFromSmiles('PCNCS')
molsmiles = Chem.MolToSmiles(m)
print(molsmiles)
Start = rdkit.Chem.Draw.ShowMol(Chem.MolFromSmiles(molsmiles))
print(Start)

ATOM_ID  = input("Target atom: ")
ATOM_Pos = []
for i in range(len(molsmiles)):
    if molsmiles[i] == ATOM_ID:
        ATOM_Pos.append(i)
        
if len(ATOM_Pos) >= 2:
    print(ATOM_Pos)
    ATOM_Pos_out = int(input("Which position would you like to attach at? "))

else:
    ATOM_Pos_out = ATOM_Pos[0]

print(ATOM_Pos_out)

def splitting (mol_in, mol_add, position):
    if position == 0:
        mol_out = mol_add + mol_in
    elif position == (len(mol_in)-1):
        mol_out = mol_in + mol_add
    else:
        mol_out = mol_in[0:(position+1)] + "(" + mol_add + ")" + mol_in[(position+1):]
    return (mol_out)

product = splitting(molsmiles, input("What would you like to add? "), ATOM_Pos_out)
print(product)
End = rdkit.Chem.Draw.ShowMol(Chem.MolFromSmiles(product))