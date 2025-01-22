# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 11:09:14 2024

@author: Saul Reid-Guest
"""
# import numpy as np
# import pandas as pd
import matplotlib.pyplot as plt
# import matplotlib.image as mpimg
# import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
Xyn = input("Do you want to use a smiles sting that is \n (1) Hardcoded in \n (2) User input \n (3) From SDF \n : ")
m = Chem.MolFromSmiles("C=CC#N")
if Xyn == "2":
    m = input("Please insert smiles string ")
elif Xyn == "3":
    file = input("Please insert file address - excluding quotation marks and replacing / with \ :")
    print(file)
    suppl = Chem.SDMolSupplier(file)
    suppl2 = Chem.SDMolSupplier("test.sdf")
    m = suppl[0]

print((Chem.MolToSmiles(m)))
img = Draw.MolToImage(m)
imgplot = plt.imshow(img)

