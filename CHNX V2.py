# -*- coding: utf-8 -*-
"""
CHNX V2
Created on Fri Oct 11 15:51:01 2024

@author: Saul Reid-Guest
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

#define CHN
ID_C = "[C]"
ID_H = "[H]"
ID_N = "[N]"

def variable_inputs(): # collect data
    # collect data from user
    Xyn = input("Do you want to calculate X? (Y/N): ") # create options for additional inputs/outputs
    MWyn = input("Do you have the MW? (Y/N): ")
    string_C = input("Percentage C: ")
    string_H = input("Percentage H: ")
    string_N = input("Percentage N: ")
    outlast = [string_C , string_H , string_N]
    if Xyn[0] == "Y" or Xyn[0] == "y": # optional input for additional element to analyse
        string_X = input("Percentage X: ")
        string_XID = input("Identity of X: ")
        XID = "[" + string_XID + "]" # Add Square brackets so that it is read as a lone element
        outlast.append(string_X)
        outlast.append(XID)
    if MWyn[0] == "Y" or MWyn[0] == "y": # optional input
        string_MW = input("Molecular weight: ")
        outlast.append(float(string_MW))
    
    return (outlast)
# no of returned items is 3 , 4, 5 , 6 define for each
test = variable_inputs()

MW = 0 # create variable set at a null value
if len(test) == 4 or len(test) == 6: # if MW= True set its value
    MW = test[-1]
percentages = [float(test[0]) , float(test[1]), float(test[2])] # create list containing only the percentages
Atom_ID = [ID_C , ID_H , ID_N]
if len(test) == 5 or len(test) == 6: # only if X is included
    percentages.append(float((test[3]))) # add % x to list
    Atom_ID.append(test[4]) # add ID of X to list
if sum(percentages) > 2: # if the sum of the percentages is greater than 2 it converts to decimal form
    tpercent = [] # open temporary list of percentages
    for x in range(len(percentages)): # for the whole list of percentages
        converttodecimal = percentages[x]/100 # change from percentage to decimal form
        tpercent.append(converttodecimal) # add the decimal form to the temporary list
    percentages = tpercent # overwrite the original list with the decimal forms
   
empirical = True # true if empirical formula is all that can be outputted
if MW != 0: # empirical is not the only output
    empirical = False
emp = []
for x in range(len(percentages)):
    tn = percentages[x] / Descriptors.MolWt(Chem.MolFromSmiles(Atom_ID[x])) # calculate %/AW before working out the ratio
    emp.append(tn)

emp_sorted = sorted(emp)
# list_pos = 0
emp_min_l = [] # find the minima so that it can be set to 1
for x in range(len(emp_sorted)):
    if emp_sorted[x] != 0:
        emp_min_l.append(emp_sorted[x]) # remove all values that = 0 to avoid maths errors

emp_min = emp_min_l[0]        #set the minima from a list of non-zero values

emp_ratio = [] # create list for ratio of empirical values
for i in range(len(emp)): # calculate the ratio for all items in the list
    tratio = emp[i] / emp_min 
    emp_ratio.append(tratio) # add empirical ratio to list for outputs

MW_adjusted = []
if empirical == False: # if the molecular weight is defined calculate molecular formula
    for i in range(len(percentages)):
        tn_MW = round((percentages[i] * MW) / Descriptors.MolWt(Chem.MolFromSmiles(Atom_ID[i]))) # ( % x MW ) / atomic weight
        MW_adjusted.append(tn_MW) # add identity to list

#calculate double bond equivalence
# RBD = C - H+X/2  + N+P/2 +1 - see: https://crawfordscientific.com/chromatography-blog/post/rdb-rule-rings-and-double-bond-equivalents
def DBE (MWs, Atoms): # Input must be in the order CHNX for this to work
    Valid = True # assume the input for x is valid - checked later
    nC = MWs[0] # assign the number of carbons
    nH = MWs[1] # assign number of nitrogens
    nN = MWs[2] # assign number of Hydrogens
    if Atoms[-1] == "[P]":
        nN = MWs[2] + MWs[-1] # add Phosphorous atoms to the nitrogen term
    elif Atoms[-1] == "[Br]" or Atoms[-1] == "[Cl]" or Atoms[-1] == "[F]" or Atoms[-1] == "[Li]":
        nH = MWs[1] + MWs[-1] # add monovalent species (e.g. halides) to the hydrogen term
    elif Atoms[-1] != "[N]" and Atoms[-1] != "[O]":
        Valid = False # mark that the atom x is specified and does not have its identity accounted for in this calculation
    DBE = nC - ((nH)/2) + (nN/2) +1
    return (DBE , Valid)


duplicate = True
if empirical == False:
    for x in range(len(emp_ratio)):
        if round(emp_ratio[x]) != round(MW_adjusted[x]):
            duplicate = False
            print("ALERT", emp_ratio, MW_adjusted)
    print("\n", "The molecular formula is:")
    for x in range(len(MW_adjusted)):
        print(Atom_ID[x] , round(MW_adjusted[x]))
    if duplicate == False:
        print("Empirical formula is")
        for x in range(len(emp_ratio)):
            print(Atom_ID[x] , round(emp_ratio[x]))
    DBE, DBE_validity = DBE(MW_adjusted , Atom_ID)
    print("DBE = ", DBE)
    if DBE_validity == False:
        print("CAUTION the DBE does not account for the x atom")

if empirical == True:
    print("\n", "The empirical formula is:")
    for x in range(len(emp_ratio)):
        print(Atom_ID[x] , round(emp_ratio[x]))
