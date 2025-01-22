# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:51:19 2024

@author: Saul Reid-Guest
"""
from rdkit import Chem
m = Chem.MolFromSmiles('C=CC#N')
with Chem.SDWriter('test1.sdf') as w:
    ID = "name1"
    mol = m
    mol.SetProp("_Name",ID)
    w.write(mol)

suppl = Chem.SDMolSupplier("test.sdf")
m2 = Chem.MolToSmiles(suppl[0])
print (m2)