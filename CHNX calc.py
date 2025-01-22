# -*- coding: utf-8 -*-
"""
Created on Sun Dec 24 16:45:06 2023

@author: Saul Reid-Guest
"""

# aims: CHNX calculator

# Set atomic weights
C = 12.011
H = 1.00784
N = 14.0067



def CHN(CP, HP , NP , XP, XW): # intake percentages mass of substituent elements and atomic mass of x
    if CP == 0: # if the mass of x is unspecified, ratio is set to 0
        CR = 0
    else:   
        CR = CP/C        # calculate the ratio of carbon in the sample
    if HP == 0: # if the mass of x is unspecified, ratio is set to 0
        HR = 0
    else:
        HR = HP/H # ^ but hydrogen
    if NP == 0: # if the mass of x is unspecified, ratio is set to 0
        NR = 0
    else:
        NR = NP/N # ^ but N
    if XW == 0 or XP == 0: # if the mass of x is unspecified, ratio is set to 0
        XR = 0
    else:
        XR = XP/XW #  calculate the ratio of X in the sample
    
    return (CR , HR , NR , XR) # output the 4 ratios

def variable_inputs(): # collect data
    # collect data from user
    string_MW = input("Molecular weight ")
    string_C = input("Percentage C ")
    string_H = input("Percentage H ")
    string_N = input("Percentage N ")
    string_X = input("Percentage X ")
    string_XW = input("Atomic weight of X ")
    
    MW = float(string_MW) # make the numbers from strings to floats
    CI = float(string_C)
    HI = float(string_H)
    NI = float(string_N)
    XI = float(string_X)
    XW = float(string_XW)
    
    total = CI + HI + NI + XI # calculate the total percentage
    if total > 1.5: # Convert percentage to decimal if needed
        CI = CI/100
        HI = HI/100
        NI = NI/100
        XI = XI/100
        
    return (CI , HI , NI , XI , XW, MW) # output the data collected

result = variable_inputs()
print (result)
MW = result[5]

E_list = ["C", "H", "N", "X" ]

if MW >= 1:
    out1 = CHN((result[0]*MW) , (result[1]*MW), (result[2]*MW), (result[3]*MW), result[4])
    print("MW specified - predicted CHNX in compound is: ")
    for x in range (0, len(out1)):
        print(E_list[x] , round(out1[x], 0))


else:
    out1 = CHN(result[0] , result[1], result[2], result[3], result[4])
    
    low = 100
    for x in range(0, len(out1)):
        if out1[x] < low and out1[x] != 0:
            low = out1[x]
       
    out2 = [] 
    for x in range (0 , len(out1)):
        if out1[x] > 0:
            getratioed = out1[x]/low
        else:
            getratioed = 0
        out2.append(getratioed)
        
    
    print("\n The empirical formula is: ")
    for x in range (0 , 4):
        print(E_list[x] , round(out2[x], 2))
        