# -*- coding: utf-8 -*-
"""
Created on 05.07.2019 

Platform for testing different infection frameworks

@author: Gustavo Hernandez-Mejia
"""

from funct_shift import *
import matplotlib.pyplot as plt
import math   
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import pandas as pd
import random
from matplotlib.font_manager import FontProperties
import csv

#______________________________________________________________________________
Indicator = '1'
NOAb=10000                   # Number of B cells 
"""           
                           
                           Variable 
                          Prelocation
"""
nbits=14           #   NUMBER OF BITS PER STRING
dsi = 168 #        day of second infection (42 days)
std = 168 + 168 #  stop day 42 d after second infection


vic_hr=np.delete(vic_hr, 0,0) 
indi_hr=np.delete(indi_hr, 0,0) 
Philipp_hr=np.delete(Philipp_hr, 0,0) 
Massa_hr=np.delete(Massa_hr, 0,0) 
Texas_hr=np.delete(Texas_hr, 0,0) 
HK68_hr=np.delete(HK68_hr, 0,0) 

Memory_cell=np.delete(Memory_cell, 0,0) 
fitness_dynamic = np.delete(fitness_dynamic, 0,0) 

p_vic_hr=np.delete(p_vic_hr, 0,0) 
p_indi_hr=np.delete(p_indi_hr, 0,0) 
p_Philipp_hr=np.delete(p_Philipp_hr, 0,0) 
p_Massa_hr=np.delete(p_Massa_hr, 0,0) 
p_Texas_hr=np.delete(p_Texas_hr, 0,0) 
p_HK68_hr=np.delete(p_HK68_hr, 0,0) 
#                          Antibody
Ab_vic_hr=np.delete(Ab_vic_hr, 0,0) 
Ab_indi_hr=np.delete(Ab_indi_hr, 0,0) 
Ab_Philipp_hr=np.delete(Ab_Philipp_hr, 0,0) 
Ab_Massa_hr=np.delete(Ab_Massa_hr, 0,0) 
Ab_Texas_hr=np.delete(Ab_Texas_hr, 0,0) 
Ab_HK68_hr=np.delete(Ab_HK68_hr, 0,0) 

Abs_total = np.delete(Abs_total, 0,0) 
Bcel_total = np.delete(Bcel_total, 0,0) 
Plas_total = np.delete(Plas_total, 0,0) 


#______________________________________________________________________________

"""                         Virus Generation
    Uniformly generate random binary strings for Philippines-82 strains
  
"""
#___________________   Params for virus generation Phill 82
NOV=10                                #   Number of virus
x_range=[shiftx +133,shiftx +141]       # range of x-axis covering
y_range=[shift +-181,shift +-173]     # range of y-axis covering
centro=[shiftx +137,shift +-177]       # Philippines 82 centroid
centro_bin = np.array([[dec2bin_sig(centro[0], nbits), 
                                               dec2bin_sig(centro[1], nbits)]])

#___________________           First Virus generation
for i in range(0,NOV):
    Vic_rand=np.concatenate((Vic_rand, np.array([[random.randint(x_range[0], 
                x_range[1]), random.randint(y_range[0],y_range[1])]])), axis=0)
Vic_rand=np.delete(Vic_rand, 0,0) #delets the very first position


#___________________          First Virus binary strings
for i in range(0,NOV):
    Vic_rand_bin=np.concatenate((Vic_rand_bin, np.array([[dec2bin_sig
         (int(Vic_rand[i,0]), nbits), dec2bin_sig(int(Vic_rand[i,1]), nbits)]]
                                                                    )), axis=0)
Vic_rand_bin=np.delete(Vic_rand_bin, 0,0) #delets the very first position

#______________________________________________________________________________

"""                         Second Virus Generation
    Uniformly generate random binary strings for Victoria-2011 strain
  
"""
##___________________   Params for virus generation Phill 82
NOV=10                                #   Number of virus
centro_bin_2 = np.array([[dec2bin_sig(centro_Victo[0], nbits), 
                                         dec2bin_sig(centro_Victo[1], nbits)]])

#___________________       Virus generation
for i in range(0,NOV):
    Vic_rand_second=np.concatenate((Vic_rand_second, np.array([[random.randint
                        (x_range_vic[0],x_range_vic[1]), random.randint
                                   (y_range_vic[0],y_range_vic[1])]])), axis=0)
Vic_rand_second=np.delete(Vic_rand_second, 0,0) #delets the very first position


#___________________       Virus binary strings
for i in range(0,NOV):
    Vic_rand_bin_sec=np.concatenate((Vic_rand_bin_sec, np.array([[dec2bin_sig
         (int(Vic_rand_second[i,0]), nbits), dec2bin_sig(int(Vic_rand_second[i,1]), nbits)]]
                                                                    )), axis=0)
Vic_rand_bin_sec=np.delete(Vic_rand_bin_sec, 0,0) #delets the very first position

#______________________________________________________________________________

"""                            GERMINAL CENTER
    
    Initial B-cells Generation
    Randomly, uniformly generated binary strings to represent the initial 
    B-cells population
  
"""
#___________________            B-cells generation

Bcell_rand=np.empty([1, 2])
for i in range(0,NOAb):
    Bcell_rand=np.concatenate((Bcell_rand, np.array([[str(np.random.randint(2, 
         size=(nbits,))),str(np.random.randint(2, size=(nbits,)))]]) ), axis=0)
Bcell_rand=np.delete(Bcell_rand, 0,0) 

#        Delete array components from genetation of binary strings
for i in range(0,NOAb):
    for j in range(0,2):
        aux = Bcell_rand[i,j]
        aux = aux.replace("[", "")               #delete [ 
        aux = aux.replace("]", "")               #delete  ]
        Bcell_rand[i,j] = aux.replace(" ", "")    #delete " "
 
Bcell_valid=np.empty([1, 3])       
for i in range(len(Bcell_rand)):
    L1 = list(Bcell_rand[i,0])
    L2 = list(Bcell_rand[i,1])
    L1[0]= '0'
    L2[0]= '0'
    Bcell_valid=np.concatenate((Bcell_valid, np.array([[''.join(L1), 
                                                ''.join(L2), i]])), axis=0)
Bcell_valid=np.delete(Bcell_valid, 0,0)  
 
Bcell_rand = Bcell_valid[:,0:2]


"""                        GERMINAL CENTER DYNAMICS
    
                  every 6 hours a generation is represented

  
"""
#    Threshold for Principal and secondary matching areas (alpha & gamma)
alpha = 3 # 3    refers to BETA in code
gamma = 3 # 3    refers to changes in secondary area ALPHA

hr = 1                  #  COUNTER OF HOURS OF SIMULATION (hr*6)

nbit_secon = 7 #  7   number of bits in secondary area + 1
nbit_princi = 7 #   7   #   number of bits in principal area without sign

param = 4 # number of matching indicators, Principal and secondary, for X and Y

plasma_list_step = []
plasma_half_life = []
phil_half_life = []

Abs_list_step = []
Abs_half_life = []

#                                PRINCIPAL LOOP
while hr <= dsi + 168:

 if hr == dsi:   
     ##########################################################################
     
     xticks=[3570, 3580, 3590, 3600, 3610, 3620, 3630, 3640]
     xticks_y=[3260, 3270, 3280, 3290, 3300, 3310]
    
     xticks1=['3570\n00 1101 1111 0010', '3580\n00 1101 1111 1100',
              '3590\n00 1110 0000 0110',
             '3600\n00 1110 0001 0000', '3610\n00 1110 0001 1010', 
             '3620\n00 1110 0010 0100', 
             '3630\n00 1110 0010 1110', '3640\n00 1110 0011 1000']
     xticksy=['3260            \n00 1100 1011 1100',
              '3270            \n00 1100 1100 0110', 
             '3280            \n00 1100 1101 0000',
             '3290            \n00 1100 1101 1010',
             '3300            \n00 1100 1110 0100', 
             '3310            \n00 1100 1110 1110']
     s2=1900
     ax=plt.figure(figsize=(11,9), facecolor='w', edgecolor='k')
     plt.scatter(bin2dec_twos(int(Bcell_rand[0,0],2), nbits), bin2dec_twos(int
          (Bcell_rand[0,1],2), nbits), marker="1", s=80, 
                                               color=Abs_color,label='B-cells')
     plt.scatter(centro[0], centro[1], s=s1, color=colors[0], marker='o',alpha=alpha1)
     plt.scatter(centro[0], centro[1], s=s2, color=colors[0], marker='o',alpha=alpha2)
     plt.scatter(centro[0], centro[1], marker="2", s=s_s1, color=colors[0], 
                 label='Philippines-1982')
     plt.scatter(centro_Victo[0], centro_Victo[1], s=s1, color=colors[3],
                                                            marker='o',alpha=alpha1)
     plt.scatter(centro_Victo[0], centro_Victo[1], s=s2, color=colors[3],
                                                            marker='o',alpha=alpha2)
     plt.scatter(centro_Victo[0], centro_Victo[1], marker="+", s=s_s1, 
                                             color=colors[3], label='Victoria-2011')
     plt.scatter(centro_indi[0], centro_indi[1], s=s1, color=colors[1],
                                                            marker='o',alpha=alpha1)
     plt.scatter(centro_indi[0], centro_indi[1], s=s2, color=colors[1],
                                                            marker='o',alpha=alpha2)
     plt.scatter(centro_indi[0], centro_indi[1], marker="4", s=s_s1, 
                                              color=colors[1], label='Indiana-2011')
     plt.scatter(centro_Massa[0], centro_Massa[1], s=s1, color=colors[2], 
                                                            marker='o',alpha=alpha1)
     plt.scatter(centro_Massa[0], centro_Massa[1], s=s2, color=colors[2], 
                                                            marker='o',alpha=alpha2)
     plt.scatter(centro_Massa[0], centro_Massa[1], marker="3", s=s_s1, 
                                        color=colors[2], label='Massachusetts-1982') 
     plt.scatter(centro_Texas[0], centro_Texas[1], s=s1, color=colors[4], 
                                                            marker='o',alpha=alpha1)
     plt.scatter(centro_Texas[0], centro_Texas[1], s=s2, color=colors[4],
                                                             marker='o',alpha=alpha2)
     plt.scatter(centro_Texas[0], centro_Texas[1], marker="4", s=s_s1, 
                                                color=colors[4], label='Texas-2004')
     plt.scatter(centro_HK68[0], centro_HK68[1], s=s1, color=colors[5],
                                                            marker='o',alpha=alpha1)
     plt.scatter(centro_HK68[0], centro_HK68[1], s=s2, color=colors[5], 
                                                            marker='o',alpha=alpha2)
     plt.scatter(centro_HK68[0], centro_HK68[1], marker="x", s=s_s1, 
                                            color=colors[5], label='Hong Kong-1968')
     for i in range(len(Vic_rand)):
         plt.scatter(Vic_rand[i,0], Vic_rand[i,1], marker="2",s=s_s1,color=colors[0],
                     alpha=0.8)
#     for i in range(len(Vic_rand_second)):
#         plt.scatter(Vic_rand_second[i,0], Vic_rand_second[i,1], marker="+", s=s_s1, 
#                                   color=colors[3],alpha=0.8)
     for i in range(len(Bcell_rand)):
         plt.scatter(bin2dec_twos(int(Bcell_rand[i,0],2), nbits), bin2dec_twos(int
          (Bcell_rand[i,1],2), nbits), marker="1", s=80, color=Abs_color,alpha=0.6)
     #for i in range(len(last_Abs)):
     #    plt.scatter(bin2dec_twos(int(last_Abs[i,0],2), nbits), bin2dec_twos(int
     #     (last_Abs[i,1],2), nbits), marker="3", s=s_s2, color=colors[8],alpha=0.6)
     plt.xlabel("AA difference", fontsize=16)
     plt.xlim(3570, 3650)
     plt.ylabel("AA difference", fontsize=16)
     plt.ylim(3240, 3325) 
     plt.title('Cross-reactome H3 bit-map, platform', fontsize=20)
     plt.legend(loc='lower right',fontsize=12)
     ax.savefig('Map_frt_platform_{}.pdf'.format(Indicator), format='pdf', dpi=1400)
     #ax.savefig('Map_10mil_{}.svg'.format(Indicator), format='svg', dpi=1400)
     plt.show()
     
###############################################################################
     
#    This module takes high affinity B-cells of the first infection to make 
#    them have a larger half-life and not participare in GCs cycles
     
     Vic_rand_bin = Vic_rand_bin_sec
     B_cells_Philipp2=np.empty([1, 2])
     Bcell_valid2 = np.empty([1, 2])

     for i in range(len(Bcell_rand)):   
         x_Bcell=bin2dec_twos(int(Bcell_rand[i,0],2), nbits)       
         y_Bcell=bin2dec_twos(int(Bcell_rand[i,1],2), nbits)       

         if (x_Bcell >= x_range_Philipp[0] and x_Bcell <= x_range_Philipp[1] and 
          y_Bcell >= y_range_Philipp[0] and y_Bcell <= y_range_Philipp[1]):
            B_cells_Philipp2=np.concatenate((B_cells_Philipp2, 
                    np.array([[Bcell_rand[i,0], Bcell_rand[i,1]]])), axis=0)
            Bcell_Philipp_num += 1
         else:
             Bcell_valid2=np.concatenate((Bcell_valid2, 
                    np.array([[Bcell_rand[i,0], Bcell_rand[i,1]]])), axis=0)
             
     B_cells_Philipp2=np.delete(B_cells_Philipp2, 0,0) 
     Bcell_valid2=np.delete(Bcell_valid2, 0,0) 
#     Bcell_rand = Bcell_valid2

     Bcell_rand=np.empty([1, 2])
     Bcell_rand2=np.empty([1, 2])
     for i in range(len(B_cells_Philipp2)):
         Bcell_rand2=np.concatenate((Bcell_rand2, np.array([[str(np.random.randint(2, 
              size=(nbits,))),str(np.random.randint(2, size=(nbits,)))]]) ), axis=0)
     Bcell_rand2=np.delete(Bcell_rand2, 0,0) 
     
     for i in range(len(B_cells_Philipp2)):
#     for i in range(0,NOAb):
         for j in range(0,2):
             aux = Bcell_rand2[i,j]
             aux = aux.replace("[", "")               #delete [ 
             aux = aux.replace("]", "")               #delete  ]
             Bcell_rand2[i,j] = aux.replace(" ", "")    #delete " "
      
     Bcell_valid=np.empty([1, 3])       
     for i in range(len(Bcell_rand2)):
         L1 = list(Bcell_rand2[i,0])
         L2 = list(Bcell_rand2[i,1])
         L1[0]= '0'
         L2[0]= '0'
         Bcell_valid=np.concatenate((Bcell_valid, np.array([[''.join(L1), 
                                                   ''.join(L2), i]])), axis=0)
     Bcell_valid=np.delete(Bcell_valid, 0,0)    
     Bcell_rand2 = Bcell_valid[:,0:2]
     
     Bcell_rand = np.concatenate((Bcell_valid2, Bcell_rand2), axis=0)
     Bcell_rand=np.delete(Bcell_rand, 0,0) 
###############################################################################     
     
#     Change of affinity threshold for second infection in ALPHA
     alpha = 3 # 3   refers to BETA in code
     gamma = 4 # 3   refers to changes in secondary area ALPHA

###############################################################################
   
    
 Bcell_valid=np.empty([1, 3]) # Bcell_valid
 Antb_x = Bcell_rand[:,0]
 Antb_y = Bcell_rand[:,1]

 
 for i in range(len(Bcell_rand)):
     Bcell_valid=np.concatenate((Bcell_valid, np.array([[Bcell_rand[i,0], 
                                                Bcell_rand[i,1], i]])), axis=0)
 Bcell_valid=np.delete(Bcell_valid, 0,0)  
 
# _____________     Principal matching area & FITNESS
  
 Match_score = np.ones((len(Bcell_valid), len(Vic_rand_bin)))
 Fitness = np.ones((len(Bcell_valid), 1))   
 
 if hr == 167 or hr == 333:
     print(hr)
 
 for i in range(len(Bcell_valid)):  
    x_Anb = Bcell_valid[i,0]
    y_Anb = Bcell_valid[i,1]
    
    for j in range(len(Vic_rand_bin)):  
        x_Virion = Vic_rand_bin[j,0]
        y_Virion = Vic_rand_bin[j,1]
        princi_mat_x = principal_area_comp(x_Virion, x_Anb, alpha, nbit_secon)
        princi_mat_y = principal_area_comp(y_Virion, y_Anb, alpha, nbit_secon)
        second_mat_x = second_area_comp(x_Virion, x_Anb, gamma, nbit_princi)
        second_mat_y = second_area_comp(y_Virion, y_Anb, gamma, nbit_princi)
        #   Matching score (MAx 58)
        Match_score[i,j] = (princi_mat_x + second_mat_x + # match_principal(29)
                            princi_mat_y + second_mat_y)  # match_secondary(29)

    Fitness[i,0] = math.exp(np.average(Match_score[i,:])) # mean de matching

 Mean_fitness = np.average(Fitness)
 
 fitness_dynamic = np.concatenate((fitness_dynamic, np.array([[Mean_fitness]]))
                                                                      , axis=0)
# 
 Bcell_valid =np.hstack((Bcell_valid,Fitness))  

# _____________                  REPLICATION
 High_affin=np.empty([1, 2])
 for i in range(len(Bcell_valid)):  
     if float(Bcell_valid[i,3]) >= Mean_fitness:
         High_affin=np.concatenate((High_affin, np.array([[Bcell_valid[i,0], 
                                                  Bcell_valid[i,1]]])), axis=0)
 High_affin=np.delete(High_affin, 0,0) #delets the very first position
 
 #  NEW B-CELLS FROM THE HIGHEST AFFINITY POPULATION
 New_Bcell=np.empty([1, 2]) # New_Bcell
 
# _____________                  MUTATION 
 for i in range(len(High_affin)):  
     New_Bcell=np.concatenate((New_Bcell, np.array([mutation(High_affin[i,0],
                                           High_affin[i,1], nbits)]) ), axis=0)
 
 New_Bcell=np.delete(New_Bcell, 0,0)  
 
 # Double b_cell production after activation (high affinity)
 High_affin_family=np.concatenate((High_affin, High_affin, New_Bcell), axis=0)

#______________________________________________________________________________

#            Probability of becoming a Memory or Plasma cell 
 
 Plasma_hr = np.empty([1, 2])
 
# if hr > 153:
#     print('Check')
 for i in range(len(High_affin)):
     if random.random() <= 0.5: 
#         Plasma_cell = np.concatenate((Plasma_cell, np.array([[High_affin[i,0],
#                                                High_affin[i,1], hr]])), axis=0)
         Plasma_hr = np.concatenate((Plasma_hr, np.array([[High_affin[i,0],
                                                   High_affin[i,1]]])), axis=0)
#     else:
#         Memory_cell = np.concatenate((Memory_cell, np.array([[High_affin[i,0],
#                                                High_affin[i,1], hr]])), axis=0)

 Plasma_hr=np.delete(Plasma_hr, 0,0)
#______________________________________________________________________________
 
#                       Half-life of B cells Phil (3d)  
 
 if hr == dsi:   
     days1 = 80
     time_stop = (24/6)*days1
     Phil_Plasma = 1200  #  1200 h (50d)
     step_size = 6     # step of 6h
     Phil_list = half_life_decay(B_cells_Philipp2, Phil_Plasma, time_stop, 
                                                                    step_size) 
 
#                       Half-life of Plasma cells (3d)   
 
 days = 70
 time_stop = (24/6)*days
 half_Plasma = 72  #  72 h (3d)
 step_size = 6     # step of 6h
 
 Plasma_list = half_life_decay(Plasma_hr, half_Plasma, time_stop, step_size) 
 plasma_half_life.append(Plasma_list)  
 
 long = 20 # origin 20
 if hr != 1:
     if len(plasma_half_life) <= long:
         cuantity = len(plasma_half_life)
     else:
         cuantity = long
        
     for i in range(cuantity):
         aux1 = plasma_half_life[i]
         if len(aux1[0]) == 0:
             continue
         if (len(plasma_half_life[len(plasma_half_life) - (i+1)])) <= i:
             continue
         aux = plasma_half_life[len(plasma_half_life) - (i+1)][i]
         if aux.size != 0:
             Plasma_hr = np.concatenate((Plasma_hr, aux), axis=0)


#____________   REGISTERING PLASMA-CELLS COUNT PER GENERATION (6h)  ___________
     
     Plas_vic_num=0
     Plas_indi_num=0
     Plas_Philipp_num=0
     Plas_Massa_num=0
     Plas_Texas_num=0
     Plas_HK68_num=0
     
     P_Abs_vic = np.empty([1, 2])
     P_Abs_indi = np.empty([1, 2])
     P_Abs_Philipp = np.empty([1, 2])
     P_Abs_Massa = np.empty([1, 2])
     P_Abs_HK68 = np.empty([1, 2])
     P_Abs_Texas = np.empty([1, 2])
     Valid_Abs_hr = np.empty([1, 2])
     
#     if hr >= 168:
#         Plasma_hr = np.concatenate((Plasma_hr, B_long_Philipp), axis=0)
#         
    
     for i in range(len(Plasma_hr)):   
        x_Pcell=bin2dec_twos(int(Plasma_hr[i,0],2), nbits)       
        y_Pcell=bin2dec_twos(int(Plasma_hr[i,1],2), nbits)       
    
        if (x_Pcell >= x_range_vic[0] and x_Pcell <= x_range_vic[1] and 
            y_Pcell >= y_range_vic[0] and y_Pcell <= y_range_vic[1]):
              P_Abs_vic = np.concatenate((P_Abs_vic, np.array([[Plasma_hr[i,0], 
                                                   Plasma_hr[i,1] ]])), axis=0)
              Plas_vic_num += 1
                
        if (x_Pcell >= x_range_indi[0] and x_Pcell <= x_range_indi[1] and 
            y_Pcell >= y_range_indi[0] and y_Pcell <= y_range_indi[1]):
            P_Abs_indi = np.concatenate((P_Abs_indi, np.array([[Plasma_hr[i,0], 
                                                   Plasma_hr[i,1] ]])), axis=0)
            Plas_indi_num += 1
        if (x_Pcell >= x_range_Philipp[0] and x_Pcell <= x_range_Philipp[1] and 
            y_Pcell >= y_range_Philipp[0] and y_Pcell <= y_range_Philipp[1]):
              P_Abs_Philipp = np.concatenate((P_Abs_Philipp, np.array([[Plasma_hr[i,0], 
                                                   Plasma_hr[i,1] ]])), axis=0)
              Plas_Philipp_num += 1
                
#                if hr >= 168:
#                    P_Abs_Philipp = np.concatenate((P_Abs_Philipp, Phil_list[hr -168]), axis=0)
#                    Plas_Philipp_num = Plas_Philipp_num + len(Phil_list[hr -168])   
                
        if (x_Pcell >= x_range_Massa[0] and x_Pcell <= x_range_Massa[1] and 
            y_Pcell >= y_range_Massa[0] and y_Pcell <= y_range_Massa[1]):
                P_Abs_Massa = np.concatenate((P_Abs_Massa, np.array([[Plasma_hr[i,0], 
                                                       Plasma_hr[i,1] ]])), axis=0)
                Plas_Massa_num += 1
        if (x_Pcell >= x_range_Texas[0] and x_Pcell <= x_range_Texas[1] and 
            y_Pcell >= y_range_Texas[0] and y_Pcell <= y_range_Texas[1]):
                P_Abs_Texas = np.concatenate((P_Abs_Texas, np.array([[Plasma_hr[i,0], 
                                                       Plasma_hr[i,1] ]])), axis=0)
                Plas_Texas_num += 1
        if (x_Pcell >= x_range_HK68[0] and x_Pcell <= x_range_HK68[1] and 
            y_Pcell >= y_range_HK68[0] and y_Pcell <= y_range_HK68[1]):
                P_Abs_HK68 = np.concatenate((P_Abs_HK68, np.array([[Plasma_hr[i,0], 
                                                       Plasma_hr[i,1] ]])), axis=0)
                Plas_HK68_num += 1
     
     
     P_Abs_vic=np.delete(P_Abs_vic, 0,0) 
     P_Abs_indi=np.delete(P_Abs_indi, 0,0) 
     P_Abs_Philipp=np.delete(P_Abs_Philipp, 0,0) 
     P_Abs_Massa=np.delete(P_Abs_Massa, 0,0) 
     P_Abs_HK68=np.delete(P_Abs_HK68, 0,0) 
     P_Abs_Texas=np.delete(P_Abs_Texas, 0,0)    
     
     if hr == dsi:   
         days1 = 80
         time_stop = (24/6)*days1
         Phil_Plasma = 480  #  480 h (20d)
         step_size = 6     # step of 6h
         Phil_list1 = half_life_decay(P_Abs_Philipp, Phil_Plasma, time_stop, step_size) 
     
     if hr >= dsi:
#         P_Abs_Philipp = np.concatenate((P_Abs_Philipp, 
#                                                  Phil_list1[hr -168] ), axis=0)
#         Plas_Philipp_num = Plas_Philipp_num + len(Phil_list1[hr -168]) 
         Plas_Philipp_num = len(Phil_list1[hr -dsi]) 
     
     p_vic_hr=np.concatenate((p_vic_hr, np.array([[hr, Plas_vic_num]])), axis=0)
     p_indi_hr=np.concatenate((p_indi_hr, np.array([[hr, Plas_indi_num]])), axis=0)
     p_Philipp_hr=np.concatenate((p_Philipp_hr, np.array([[hr, Plas_Philipp_num]])),axis=0)
     p_Massa_hr=np.concatenate((p_Massa_hr, np.array([[hr, Plas_Massa_num]])), axis=0)
     p_Texas_hr=np.concatenate((p_Texas_hr, np.array([[hr, Plas_Texas_num]])), axis=0)
     p_HK68_hr=np.concatenate((p_HK68_hr, np.array([[hr, Plas_HK68_num]])), axis=0)
     
     Plas_total = np.concatenate((Plas_total, np.array([[hr, Plas_vic_num, 
        Plas_indi_num, Plas_Philipp_num, Plas_Massa_num, Plas_Texas_num, 
                                                     Plas_HK68_num]])), axis=0)


#______________________________________________________________________________
#            

# _________________   Uniformly at Random a B-cell is selected to die 
#                     In its place, the new B-cell is added

 for i in range(len(High_affin_family)): 
     point = np.random.randint(len(Bcell_rand), size=(1,))
     point = point[0]  
     Bcell_rand[point] = High_affin_family[i]

# 
# _____________                  MUTATION 
 
 point = np.random.randint(len(Bcell_rand), size=(1,))
 point = point[0]
 mutated = mutation(Bcell_rand[point,0], Bcell_rand[point,1], nbits)
 Bcell_rand[point,0] = mutated[0]
 Bcell_rand[point,1] = mutated[1] 
#______________________________________________________________________________
 
 
#____________   REGISTERING B-CELLS COUNT PER GENERATION (6h) _________________

 #Bcells_counting(Bcell_rand, nbits)
 
 Bcell_vic_num=0
 Bcell_indi_num=0
 Bcell_Philipp_num=0
 Bcell_Massa_num=0
 Bcell_Texas_num=0
 Bcell_HK68_num=0
 
 B_long_Philipp = np.empty([1, 2])

 for i in range(len(Bcell_rand)):   
    x_Bcell=bin2dec_twos(int(Bcell_rand[i,0],2), nbits)       
    y_Bcell=bin2dec_twos(int(Bcell_rand[i,1],2), nbits)       

    if (x_Bcell >= x_range_vic[0] and x_Bcell <= x_range_vic[1] and 
        y_Bcell >= y_range_vic[0] and y_Bcell <= y_range_vic[1]):
#            B_cells_vic=np.concatenate((B_cells_vic, np.array([[bin2dec_twos
#          (int(Bcell_rand[i,0],2), nbits), bin2dec_twos(int(Bcell_rand[i,1],2),
#                                                         nbits),hr]])), axis=0)
            Bcell_vic_num += 1
            
    if (x_Bcell >= x_range_indi[0] and x_Bcell <= x_range_indi[1] and 
        y_Bcell >= y_range_indi[0] and y_Bcell <= y_range_indi[1]):
#            B_cells_indi=np.concatenate((B_cells_indi, np.array([[bin2dec_twos
#          (int(Bcell_rand[i,0],2), nbits), bin2dec_twos(int(Bcell_rand[i,1],2),
#                                                        nbits),hr ]])), axis=0)
            Bcell_indi_num += 1
    if (x_Bcell >= x_range_Philipp[0] and x_Bcell <= x_range_Philipp[1] and 
        y_Bcell >= y_range_Philipp[0] and y_Bcell <= y_range_Philipp[1]):
            B_cells_Philipp=np.concatenate((B_cells_Philipp, np.array([[bin2dec_twos
            (int(Bcell_rand[i,0],2), nbits), bin2dec_twos(int(Bcell_rand[i,1],2),
                                                        nbits) ,hr]])), axis=0)
            B_long_Philipp = np.concatenate((B_long_Philipp, np.array([[Bcell_rand[i,0], 
                                                       Bcell_rand[i,1] ]])), axis=0)
            Bcell_Philipp_num += 1
    if (x_Bcell >= x_range_Massa[0] and x_Bcell <= x_range_Massa[1] and 
        y_Bcell >= y_range_Massa[0] and y_Bcell <= y_range_Massa[1]):
#            B_cells_Massa=np.concatenate((B_cells_Massa, np.array([[bin2dec_twos
#            (int(Bcell_rand[i,0],2), nbits), bin2dec_twos(int(Bcell_rand[i,1],2),
#                                                        nbits) ,hr]])), axis=0)
            Bcell_Massa_num += 1
    if (x_Bcell >= x_range_Texas[0] and x_Bcell <= x_range_Texas[1] and 
        y_Bcell >= y_range_Texas[0] and y_Bcell <= y_range_Texas[1]):
#            B_cells_Texas=np.concatenate((B_cells_Texas, np.array([[bin2dec_twos
#            (int(Bcell_rand[i,0],2), nbits), bin2dec_twos(int(Bcell_rand[i,1],2),
#                                                        nbits) ,hr]])), axis=0)
            Bcell_Texas_num += 1
    if (x_Bcell >= x_range_HK68[0] and x_Bcell <= x_range_HK68[1] and 
        y_Bcell >= y_range_HK68[0] and y_Bcell <= y_range_HK68[1]):
#            B_cells_HK68=np.concatenate((B_cells_HK68, np.array([[bin2dec_twos
#            (int(Bcell_rand[i,0],2), nbits), bin2dec_twos(int(Bcell_rand[i,1],2),
#                                                        nbits),hr ]])), axis=0)
            Bcell_HK68_num += 1
 
 B_long_Philipp=np.delete(B_long_Philipp, 0,0) 
 if hr >= dsi:
     B_long_Philipp = np.concatenate((B_long_Philipp, 
                                          Phil_list[hr -dsi] ), axis=0)
     Bcell_Philipp_num = Bcell_Philipp_num + len(Phil_list[hr -dsi])    
    
    
    
    
 vic_hr=np.concatenate((vic_hr, np.array([[hr, Bcell_vic_num]])), axis=0)
 indi_hr=np.concatenate((indi_hr, np.array([[hr, Bcell_indi_num]])), axis=0)
 Philipp_hr=np.concatenate((Philipp_hr, np.array([[hr, Bcell_Philipp_num]])),axis=0)
 Massa_hr=np.concatenate((Massa_hr, np.array([[hr, Bcell_Massa_num]])), axis=0)
 Texas_hr=np.concatenate((Texas_hr, np.array([[hr, Bcell_Texas_num]])), axis=0)
 HK68_hr=np.concatenate((HK68_hr, np.array([[hr, Bcell_HK68_num]])), axis=0)
 
 Bcel_total = np.concatenate((Bcel_total, np.array([[hr, Bcell_vic_num, 
            Bcell_indi_num, Bcell_Philipp_num, Bcell_Massa_num, 
                                   Bcell_Texas_num, Bcell_HK68_num]])), axis=0)

#______________________________________________________________________________
#  SAVE DATA
 
 with open('Bcel_fina_{}.csv'.format(Indicator), 'w') as writeFile:
     writer = csv.writer(writeFile)
     writer.writerows(Bcel_total)
 writeFile.close()
 
 with open('Plas_final_{}.csv'.format(Indicator), 'w') as writeFile:
     writer = csv.writer(writeFile)
     writer.writerows(Plas_total)
 writeFile.close()
 
 Bcell_valid=np.empty([1, 2]) 
 for i in range(len(Bcell_rand)):
     Bcell_valid=np.concatenate((Bcell_valid, np.array([[bin2dec_twos(int(
     Bcell_rand[i,0],2), nbits), bin2dec_twos(int(Bcell_rand[i,1],2), nbits) 
                                                                  ]])), axis=0)

 Bcell_valid=np.delete(Bcell_valid, 0,0)  
 with open('Bcell_{}.csv'.format(Indicator), 'w') as writeFile:
     writer = csv.writer(writeFile)
     writer.writerows(Bcell_valid)
 writeFile.close()

 print('hr',hr)
 hr += 1   
#                        CLOSING PRINCIPAL LOOP
 
###############################################################################


#                              PLOTTING
s3 = 80
colors.extend(colors)
colors.extend(colors)
# _____________________________________________________________________________


"""                               PLOTTING
"""

xticks=[3570, 3580, 3590, 3600, 3610, 3620, 3630, 3640]
xticks_y=[3260, 3270, 3280, 3290, 3300, 3310]

xticks1=['3570\n00 1101 1111 0010', '3580\n00 1101 1111 1100',
         '3590\n00 1110 0000 0110',
         '3600\n00 1110 0001 0000', '3610\n00 1110 0001 1010',
         '3620\n00 1110 0010 0100', 
         '3630\n00 1110 0010 1110', '3640\n00 1110 0011 1000']
xticksy=['3260            \n00 1100 1011 1100',
         '3270            \n00 1100 1100 0110', 
         '3280            \n00 1100 1101 0000',
         '3290            \n00 1100 1101 1010',
         '3300            \n00 1100 1110 0100', 
         '3310            \n00 1100 1110 1110']
s2=1900
ax=plt.figure(figsize=(11,9), facecolor='w', edgecolor='k')
plt.scatter(bin2dec_twos(int(Bcell_rand[0,0],2), nbits), bin2dec_twos(int
     (Bcell_rand[0,1],2), nbits), marker="1", s=80, color=Abs_color,label='B-cells')
plt.scatter(centro[0], centro[1], s=s1, color=colors[0], marker='o',alpha=alpha1)
plt.scatter(centro[0], centro[1], s=s2, color=colors[0], marker='o',alpha=alpha2)
plt.scatter(centro[0], centro[1], marker="2", s=s_s1, color=colors[0], 
            label='Philippines-1982')
plt.scatter(centro_Victo[0], centro_Victo[1], s=s1, color=colors[3],
                                                       marker='o',alpha=alpha1)
plt.scatter(centro_Victo[0], centro_Victo[1], s=s2, color=colors[3],
                                                       marker='o',alpha=alpha2)
plt.scatter(centro_Victo[0], centro_Victo[1], marker="+", s=s_s1, 
                                        color=colors[3], label='Victoria-2011')
plt.scatter(centro_indi[0], centro_indi[1], s=s1, color=colors[1],
                                                       marker='o',alpha=alpha1)
plt.scatter(centro_indi[0], centro_indi[1], s=s2, color=colors[1],
                                                       marker='o',alpha=alpha2)
plt.scatter(centro_indi[0], centro_indi[1], marker="4", s=s_s1, 
                                         color=colors[1], label='Indiana-2011')
plt.scatter(centro_Massa[0], centro_Massa[1], s=s1, color=colors[2], 
                                                       marker='o',alpha=alpha1)
plt.scatter(centro_Massa[0], centro_Massa[1], s=s2, color=colors[2], 
                                                       marker='o',alpha=alpha2)
plt.scatter(centro_Massa[0], centro_Massa[1], marker="3", s=s_s1, 
                                   color=colors[2], label='Massachusetts-1982') 
plt.scatter(centro_Texas[0], centro_Texas[1], s=s1, color=colors[4], 
                                                       marker='o',alpha=alpha1)
plt.scatter(centro_Texas[0], centro_Texas[1], s=s2, color=colors[4],
                                                       marker='o',alpha=alpha2)
plt.scatter(centro_Texas[0], centro_Texas[1], marker="4", s=s_s1, 
                                           color=colors[4], label='Texas-2004')
plt.scatter(centro_HK68[0], centro_HK68[1], s=s1, color=colors[5],
                                                       marker='o',alpha=alpha1)
plt.scatter(centro_HK68[0], centro_HK68[1], s=s2, color=colors[5], 
                                                       marker='o',alpha=alpha2)
plt.scatter(centro_HK68[0], centro_HK68[1], marker="x", s=s_s1, 
                                       color=colors[5], label='Hong Kong-1968')
for i in range(len(Vic_rand)):
    plt.scatter(Vic_rand[i,0], Vic_rand[i,1], marker="2",s=s_s1,color=colors[0],
                alpha=0.8)
for i in range(len(Vic_rand_second)):
    plt.scatter(Vic_rand_second[i,0], Vic_rand_second[i,1], marker="+", s=s_s1, 
                              color=colors[3],alpha=0.8)
for i in range(len(Bcell_rand)):
    plt.scatter(bin2dec_twos(int(Bcell_rand[i,0],2), nbits), bin2dec_twos(int
     (Bcell_rand[i,1],2), nbits), marker="1", s=80, color=Abs_color,alpha=0.6)
#for i in range(len(last_Abs)):
#    plt.scatter(bin2dec_twos(int(last_Abs[i,0],2), nbits), bin2dec_twos(int
#     (last_Abs[i,1],2), nbits), marker="3", s=s_s2, color=colors[8],alpha=0.6)
plt.xlabel("AA difference", fontsize=16)
plt.xlim(3570, 3650)
#plt.xlim(3450, 3900)
#plt.xticks(xticks, xticks1, rotation='vertical', fontsize=14)
#plt.margins(0.2)
plt.ylabel("AA difference", fontsize=16)
#plt.ylim(3100, 3480) 
plt.ylim(3240, 3325) 
#plt.yticks(xticks_y, xticksy, rotation='horizontal', fontsize=14)
plt.title('Cross-reactome H3 bit-map, platform', fontsize=20)
plt.legend(loc='lower right',fontsize=12)
ax.savefig('Map_seco_platform_{}.pdf'.format(Indicator), format='pdf', dpi=1400)
#ax.savefig('Map_10mil_{}.svg'.format(Indicator), format='svg', dpi=1400)
plt.show()


#
#______________________________________________________________________________
#
ax1=plt.figure(2, figsize=(11,7), facecolor='w', edgecolor='k')
plt.semilogy((Philipp_hr[:,0])*6/24, Philipp_hr[:,1], '--', color=colors[0], 
                                                      label='Philippines-1982')
plt.semilogy((vic_hr[:,0])*6/24, vic_hr[:,1], '-*', color=colors[3], 
                                                         label='Victoria-2011')
plt.semilogy((indi_hr[:,0])*6/24, indi_hr[:,1], '-.', color=colors[1], 
                                                          label='Indiana-2011')
plt.semilogy((Massa_hr[:,0])*6/24, Massa_hr[:,1], '-',color=colors[2], 
                                                    label='Massachusetts-1982')
plt.semilogy((Texas_hr[:,0])*6/24, Texas_hr[:,1], marker="3", markersize=8, 
                                           color=colors[4], label='Texas-2004')
plt.semilogy((HK68_hr[:,0])*6/24, HK68_hr[:,1], '-+', color=colors[5], 
                                                                 label='HK-68')
plt.xlabel("Time (days)", fontsize=16)
plt.ylabel("B cells (counts)", fontsize=16)
plt.title('Matching B-cells per H3N2 strain, platform', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(loc='upper left',fontsize=12)
ax1.savefig('B_platform_{}.pdf'.format(Indicator), format='pdf', dpi=1400)
plt.show()
#

ax4=plt.figure(4, figsize=(11,7), facecolor='w', edgecolor='k')
plt.semilogy((p_Philipp_hr[:,0])*6/24, p_Philipp_hr[:,1], '--', color=colors[0], 
                                                      label='Philippines-1982')
plt.semilogy((p_vic_hr[:,0])*6/24, p_vic_hr[:,1], '-*', color=colors[3], 
                                                         label='Victoria-2011')
plt.semilogy((p_indi_hr[:,0])*6/24, p_indi_hr[:,1], '-.', color=colors[1], 
                                                          label='Indiana-2011')
plt.semilogy((p_Massa_hr[:,0])*6/24, p_Massa_hr[:,1], '-',color=colors[2], 
                                                    label='Massachusetts-1982')
plt.semilogy((p_Texas_hr[:,0])*6/24, p_Texas_hr[:,1], marker="3", markersize=8, 
                                           color=colors[4], label='Texas-2004')
plt.semilogy((p_HK68_hr[:,0])*6/24, p_HK68_hr[:,1], '-+', color=colors[5], 
                                                                 label='HK-68')
plt.xlabel("Time (days)", fontsize=16)
plt.ylabel("Plasma cells (counts)", fontsize=16)
plt.title('Plasma cells per H3 strain, platform', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(loc='upper left',fontsize=12)
ax4.savefig('Pl_platform_{}.pdf'.format(Indicator), format='pdf', dpi=1400)
plt.show()


print('Done Gus 10mil {}'.format(Indicator))


#_____________________       EOF      _________________________________________ 
