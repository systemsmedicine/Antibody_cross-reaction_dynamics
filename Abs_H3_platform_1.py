# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 15:23:19 2019

@author: Gustavo
"""
from funct_shift import *
import matplotlib.pyplot as plt
import math   
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import pandas as pd
from collections import OrderedDict, defaultdict
import random
from matplotlib.font_manager import FontProperties
import csv

indica= '1'

prop_cycle = plt.rcParams['axes.prop_cycle'] #Colors
colors = prop_cycle.by_key()['color']
Abs_color=colors[2]
s3 = 90
colors.extend(colors)

multiplier = 5
Abs_half_life = []
Abs_list_Vic = []
Abs_list_ind = []
Abs_list_Phi = []
Abs_list_Mas = []
Abs_list_Txs = []
Abs_list_HK = []
days_abs = 90
stop_abs = (24/6)*days_abs
half_abs = 240  #  240 h (10 d)
step_size = 6 

Plas_total = pd.read_csv('Plas_final_{}.csv'.format(indica))
#Plasma_out = np.empty([1, 1])
P_Vic = np.array([Plas_total.iloc[:,1]]).T
P_ind = np.array([Plas_total.iloc[:,2]]).T
P_Phi = np.array([Plas_total.iloc[:,3]]).T
P_Mas = np.array([Plas_total.iloc[:,4]]).T
P_Txs = np.array([Plas_total.iloc[:,5]]).T
P_Hk8 = np.array([Plas_total.iloc[:,6]]).T


for i in range(len(P_Phi)):  
#    if i>=315:
#        print(i)
    Abs_hr = P_Phi[i,0] * multiplier
    life = np.arange(Abs_hr)
    Abs_list = py_discrete_decay(life, half_abs, step_size, stop_abs) 
    Abs_half_life.append(Abs_list)  
    for i in range(len(Abs_half_life)):
        aux1 = Abs_half_life[i]
        if aux1[0] == 0 or (len(Abs_half_life) > len(Abs_half_life[i])):
            continue
        aux = Abs_half_life[i][len(Abs_half_life) - (i+1)]
        if aux.size != 0:
            Abs_hr = Abs_hr + aux
    Abs_list_Phi.append(Abs_hr)
    

Abs_half_life = []
for i in range(len(P_Vic)):  
    Abs_hr = P_Vic[i,0] * multiplier
    life = np.arange(Abs_hr)
    Abs_list = py_discrete_decay(life, half_abs, step_size, stop_abs) 
    Abs_half_life.append(Abs_list)  
    for i in range(len(Abs_half_life)):
        aux1 = Abs_half_life[i]
        if aux1[0] == 0 or (len(Abs_half_life) > len(Abs_half_life[i])):
            continue
        aux = Abs_half_life[i][len(Abs_half_life) - (i+1)]
        if aux.size != 0:
            Abs_hr = Abs_hr + aux
    Abs_list_Vic.append(Abs_hr)

    
Abs_half_life = []
for i in range(len(P_ind)):  
    Abs_hr = P_ind[i,0] * multiplier
    life = np.arange(Abs_hr)
    Abs_list = py_discrete_decay(life, half_abs, step_size, stop_abs) 
    Abs_half_life.append(Abs_list)  
    for i in range(len(Abs_half_life)):
        aux1 = Abs_half_life[i]
        if aux1[0] == 0 or (len(Abs_half_life) > len(Abs_half_life[i])):
            continue
        aux = Abs_half_life[i][len(Abs_half_life) - (i+1)]
        if aux.size != 0:
            Abs_hr = Abs_hr + aux
    Abs_list_ind.append(Abs_hr)


Abs_half_life = []        
for i in range(len(P_Mas)):  
    Abs_hr = P_Mas[i,0] * multiplier
    life = np.arange(Abs_hr)
    Abs_list = py_discrete_decay(life, half_abs, step_size, stop_abs) 
    Abs_half_life.append(Abs_list)  
    for i in range(len(Abs_half_life)):
        aux1 = Abs_half_life[i]
        if aux1[0] == 0 or (len(Abs_half_life) > len(Abs_half_life[i])):
            continue
        aux = Abs_half_life[i][len(Abs_half_life) - (i+1)]
        if aux.size != 0:
            Abs_hr = Abs_hr + aux
    Abs_list_Mas.append(Abs_hr)

Abs_half_life = []        
for i in range(len(P_Txs)):  
    Abs_hr = P_Txs[i,0] * multiplier
    life = np.arange(Abs_hr)
    Abs_list = py_discrete_decay(life, half_abs, step_size, stop_abs) 
    Abs_half_life.append(Abs_list)  
    for i in range(len(Abs_half_life)):
        aux1 = Abs_half_life[i]
        if aux1[0] == 0 or (len(Abs_half_life) > len(Abs_half_life[i])):
            continue
        aux = Abs_half_life[i][len(Abs_half_life) - (i+1)]
        if aux.size != 0:
            Abs_hr = Abs_hr + aux
    Abs_list_Txs.append(Abs_hr)


Abs_half_life = []
for i in range(len(P_Hk8)):  
    Abs_hr = P_Hk8[i,0] * multiplier
    life = np.arange(Abs_hr)
    Abs_list = py_discrete_decay(life, half_abs, step_size, stop_abs) 
    Abs_half_life.append(Abs_list)  
    for i in range(len(Abs_half_life)):
        aux1 = Abs_half_life[i]
        if aux1[0] == 0 or (len(Abs_half_life) > len(Abs_half_life[i])):
            continue
        aux = Abs_half_life[i][len(Abs_half_life) - (i+1)]
        if aux.size != 0:
            Abs_hr = Abs_hr + aux
    Abs_list_HK.append(Abs_hr)


time = np.arange(len(P_Hk8))
ax5=plt.figure(5, figsize=(11,7), facecolor='w', edgecolor='k')
plt.semilogy((time)*6/24, Abs_list_Phi, '--', color=colors[0], 
                                                      label='Philippines-1982')
plt.semilogy((time)*6/24, Abs_list_Vic, '-*', color=colors[3], 
                                                         label='Victoria-2011')
plt.semilogy((time)*6/24, Abs_list_ind, '-.', color=colors[1], 
                                                          label='Indiana-2011')
plt.semilogy((time)*6/24, Abs_list_Mas, '-',color=colors[2], 
                                                    label='Massachusetts-1982')
plt.semilogy((time)*6/24, Abs_list_Txs, marker="3", markersize=8, 
                                           color=colors[4], label='Texas-2004')
plt.semilogy((time)*6/24, Abs_list_HK, '-+', color=colors[5], 
                                                                 label='HK-68')
plt.ylim(bottom=-1)
plt.xlabel("Time (days)", fontsize=16)
plt.ylabel("Abs cells (counts)", fontsize=16)
plt.title('Antibodies per H3 strain, platform', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(loc='lower right',fontsize=12)
ax5.savefig('Abs_platform_{}.pdf'.format(indica), format='pdf', dpi=1400)
plt.show()




Abs_Phi = np.array(Abs_list_Phi)
Abs_Vic = np.array(Abs_list_Vic)
Abs_ind = np.array(Abs_list_ind)
Abs_Mas = np.array(Abs_list_Mas)
Abs_Txs = np.array(Abs_list_Txs)
Abs_Hk8 = np.array(Abs_list_HK)
Abs_total = [Abs_Vic, Abs_ind, Abs_Phi, Abs_Mas, Abs_Txs, Abs_Hk8]


with open('Abs_total_platform_{}.csv'.format(indica), 'w') as writeFile:
    writer = csv.writer(writeFile)
    writer.writerows(Abs_total)
writeFile.close()

index = 168
first = [Abs_Vic[index], Abs_ind[index], Abs_Phi[index], Abs_Mas[index],
        Abs_Txs[index], Abs_Hk8[index]]

file = open('Abs_first_platform_{}.txt'.format(indica),'w') 
file.write(str(first)) 
file.close() 


index = len(Abs_Phi)-1
last = [Abs_Vic[index], Abs_ind[index], Abs_Phi[index], Abs_Mas[index],
        Abs_Txs[index], Abs_Hk8[index]]

file = open('Abs_last_platform_{}.txt'.format(indica),'w') 
file.write(str(last)) 
file.close() 








   