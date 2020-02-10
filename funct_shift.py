# -*- coding: utf-8 -*-
"""
Created on Mar  29 2019


@author: Gustavo Hernandez-Mejia
"""

#import matplotlib.pyplot as plt
#import math   
#from mpl_toolkits.mplot3d import axes3d
import numpy as np
import matplotlib.pyplot as plt
import math   
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import pandas as pd
import random

#import pandas as pd
#import random
#______________________________________________________________________________

#_______________      Hemagglutinin (HAs)              ________________________
#shift = 3500-25
#shift = 3500-35 # 20  ·33 ok in 2nd  35 ok for both 1-2nd
#shiftx = 3500-13 # 17  13 ok with 2nd 

shift = 3500-32 # 20  ·33 ok in 2nd  35 ok for both 1-2nd
shiftx = 3500-24 # 17  13 ok with 2nd 

#27-17 : Overgrow in 1st, 2nd decays
#17-17 : Overgrow in 1st, 2nd decays
H1 = [[shift +66,	shiftx +157, 'H1',	'H1-California-2009'],
      [shift +106,	shiftx +134, 'H1',	'H1-Puerto Rico-1934'],
      [shift +116,	shiftx +153, 'H1',	'H1-New Caledonia-1999'],
      [shift +99,	shiftx +148, 'H1',	'H1-Fort Monmouth-1947'],
      [shift +82,	shiftx +126, 'H1',	'H1-South Carolina-1918'],
      [shift +101,	shiftx +181, 'H1',	'H1-Jiangsu-2011']]


H3 = [[shift +-173,	shiftx +159, 'H3',	'H3-Victoria-2011'],
      [shift +-177,	shiftx +137, 'H3',	'H3-Philippines-1982'],
      [shift +-189,	shiftx +147, 'H3',	'H3-Indiana-2011'],
      [shift +-178,	shiftx +110, 'H3',	'H3-Massachusetts-2011'],
      [shift +-205,	shiftx +106, 'H3',	'H3-Texas-2004'],
      [shift +-174,	shiftx +121, 'H3',	'H3-Hong Kong-1968']]


B = [shift +-44,	shiftx +-393,	'B',	'Yamagata-Unknown']

HAs = [[ shift +177,shiftx +79, 	'H2',	'H2-Unknown-1999'],
      [shift + 180,	shiftx +88,	'H2',	'H2-Japan-1957'],
      [shift +-192, shiftx +19,	'H4',	'H4-Czechoslovakia-1956'],
      [shift + 152,	shiftx +34,	'H5',	'H5-Vietnam-2004'],
      [shift + 160,	shiftx +  32,    'H5',	'H5-Indonesia-2005'],
      [shift +  64,	shiftx +   4,	'H6',	'H6-Unknown-2002'],
      [shift +-210, shiftx +-120,	'H7',	'H7-Jalisco-2012'],
      [shift +-172,	shiftx +-149,	'H7',	'H7-Shanghai-2013'],
      [shift +  46,	shiftx +-177,	'H8',	'H8-Sweden-2002'],
      [shift +  22,	shiftx +-120,	'H9',	'H9-Hong Kong-1999'],
      [shift +-150,	shiftx + -99,	'H10',	'H10-Unknown-2010'],
      [shift + 113,	shiftx + -90,	'H11',	'H11-Netherlands-1999'],
      [shift +   5,	shiftx +-191,	'H12',	'H12-Unknown-2007'],
      [shift + 135,	shiftx +-185,	'H13',	'H13-Sweden-1999'],
      [shift + 135,	shiftx +-185,	'H13',	'H13-Sweden-1999'],
      [shift +-226,	shiftx +  19,	'H14',	'H14-Astrakhan-1982'],
      [shift +-197,	shiftx +-148,	'H15',	'H15-shearwater-1983'],
      [shift + 162,	shiftx +-170,	'H16',	'H16-Sweden-1999'],
      [shift + 239,	shiftx + -82, 	'H17',	'H17-Guatemala-2010'],      
      [shift + 258,	shiftx + -32,	'H18',	'H18-Peru-2010']]

#______________________________________________________________________________

#_________________________      MARKERS      __________________________________

mark = ['2','3','v','8','s','o','+','*','x',',','p','.','<','>',
        '+','*','x',',','p','.','<','>','2','3','v','8','s','o',
        'o','+','*','x',',']


#______________________________________________________________________________

#__________________         ANTIGEN/VIRUS 
# variables for plotting
s_s1=150
alpha1=0.11
alpha2=0.17
s1=5990
s2=3100
s_s2=70

prop_cycle = plt.rcParams['axes.prop_cycle'] #Colors
colors = prop_cycle.by_key()['color']
Vic11=colors[0]
Vic11_2=colors[1]
Abs_color=colors[2]

# variables for virus
Vic_rand=np.empty([1, 2])
Vic_rand_bin=np.empty([1, 2])
Vic_rand_second=np.empty([1, 2])
Vic_rand_bin_sec=np.empty([1, 2])

# variables for viral strains
B_cells_vic=np.empty([1, 3])
B_cells_indi=np.empty([1, 3])
B_cells_Philipp=np.empty([1, 3])
B_cells_Massa=np.empty([1, 3])
B_cells_Texas=np.empty([1, 3])
B_cells_HK68=np.empty([1, 3])

# variables for initial viral strains (prior to infection)
B_cells_vic_init=np.empty([1, 2])
B_cells_indi_init=np.empty([1, 2])
B_cells_Philipp_init=np.empty([1, 2])
B_cells_Massa_init=np.empty([1, 2])
B_cells_Texas_init=np.empty([1, 2])
B_cells_HK68_init=np.empty([1, 2])
# variables for final viral strains (after infection)
B_cells_vic_final=np.empty([1, 2])
B_cells_indi_final=np.empty([1, 2])
B_cells_Philipp_final=np.empty([1, 2])
B_cells_Massa_final=np.empty([1, 2])
B_cells_Texas_final=np.empty([1, 2])
B_cells_HK68_final=np.empty([1, 2])



#________________________     RANGES AND CENTROIDS CLUSTER
#                                 Victoria-2011
x_range_vic=[shiftx +154,shiftx +164]
y_range_vic=[shift +-178,shift +-168]                      
centro_Victo=[shiftx +159,shift +-173] 
#                                 Indiana-2011
x_range_indi=[shiftx +142,shiftx +152]  
y_range_indi=[shift +-194,shift +-184]                           
centro_indi=[shiftx +147,shift +-189]
#                               Philippines 1982
x_range_Philipp=[shiftx +132,shiftx +142]  
y_range_Philipp=[shift +-182,shift +-172]                          
centro_Philipp=[shiftx +137,shift +-177]
#                              Massachusetts-2011
x_range_Massa=[shiftx +105,shiftx +115]  
y_range_Massa=[shift +-183,shift +-173]                          
centro_Massa=[shiftx +110,shift +-178] 
#                                 Texas-2004
x_range_Texas=[shiftx +101,shiftx +111]   
y_range_Texas=[shift +-210,shift +-200]                        
centro_Texas=[shiftx +106,shift +-205] 
#                               Hong Kong-1968
x_range_HK68=[shiftx +116,shiftx +126]   
y_range_HK68=[shift +-179,shift +-169]                        
centro_HK68=[shiftx +121,shift +-174] 


#__________________            B CELLS - ABS
# variables for B cell acumulation per strain, per generation time (6h)
vic_hr = np.empty([1, 2])
indi_hr = np.empty([1, 2])
Philipp_hr = np.empty([1, 2])
Massa_hr = np.empty([1, 2])
Texas_hr = np.empty([1, 2])
HK68_hr = np.empty([1, 2])
#________    B CELLS - ABS amount
Plasma_cell = np.empty([1, 3])
Memory_cell = np.empty([1, 3])
Abs_serum = np.empty([1, 3])
# auxiliar for deleting random positions from arrey
point_clear = np.empty([1, 1]) 
fitness_dynamic = np.empty([1, 1])


#__________________            PLASMA CELLS
#P_cells_vic=np.empty([1, 3])
#P_cells_indi=np.empty([1, 3])
#P_cells_Philipp=np.empty([1, 3])
#P_cells_Massa=np.empty([1, 3])
#P_cells_Texas=np.empty([1, 3])
#P_cells_HK68=np.empty([1, 3])

p_vic_hr = np.empty([1, 2])
p_indi_hr = np.empty([1, 2])
p_Philipp_hr = np.empty([1, 2])
p_Massa_hr = np.empty([1, 2])
p_Texas_hr = np.empty([1, 2])
p_HK68_hr = np.empty([1, 2])

#__________________            ANTIBODIES
Abs_cells_vic=np.empty([1, 3])
Abs_cells_indi=np.empty([1, 3])
Abs_cells_Philipp=np.empty([1, 3])
Abs_cells_Massa=np.empty([1, 3])
Abs_cells_Texas=np.empty([1, 3])
Abs_cells_HK68=np.empty([1, 3])

Ab_vic_hr = np.empty([1, 2])
Ab_indi_hr = np.empty([1, 2])
Ab_Philipp_hr = np.empty([1, 2])
Ab_Massa_hr = np.empty([1, 2])
Ab_Texas_hr = np.empty([1, 2])
Ab_HK68_hr = np.empty([1, 2])

Abs_total = np.empty([1, 7])
Abs_total = np.concatenate((Abs_total, np.array([['hr', 'Vic', 
            'Indi', 'Phil', 'Mass', 'Tx', 'HK68']])), axis=0)
Bcel_total = np.empty([1, 7])
Bcel_total = np.concatenate((Bcel_total, np.array([['hr', 'Vic', 
            'Indi', 'Phil', 'Mass', 'Tx', 'HK68']])), axis=0)
Plas_total = np.empty([1, 7])
Plas_total = np.concatenate((Plas_total, np.array([['hr', 'Vic', 
            'Indi', 'Phil', 'Mass', 'Tx', 'HK68']])), axis=0)


#______________________________________________________________________________

""" Function to convert a signed integer base 10 into a two's complement signed
    number of size N (nbits)
    N=10, the actual size of virus coordinate strings. 
    The complete representation of virus is 2N long.
    
    input: s_dec: decimal signed number in the range (-(2^N)/2, ((2^N)/2)-1)
           nbits: N
    output: signed two's complement number of size N

"""
def dec2bin_sig(s_dec, nbits):
    s_bin = bin(s_dec & int("1"*nbits, 2))[2:]
    return ("{0:0>%s}" % (nbits)).format(s_bin)

#______________________________________________________________________________

""" Convert a two's complement signed binary integer of size N (nbits)
    to a decimal (base 10) signed integer
    The complete representation of virus is 2N long.
    Numbers should be in the range (-(2^N)/2, ((2^N)/2)-1)
    
    input. two_bin: signed two's complement number of size N
           nbits: N 
    output. decimal signed number

"""
def bin2dec_twos(two_bin, nbits): 
    """compute the 2's complement of int value val"""
    if (two_bin & (1 << (nbits - 1))) != 0: # sign bit 
        two_bin = two_bin - (1 << nbits)    # compute negative value
    return two_bin                          # return positive value 

#______________________________________________________________________________

"""                        Antigen matching score
    The matching score is computed for the antigen coordinates x and y separa-
    telly.
    
    The philosophy of the match score is a hybrid of binary sequences with
    matching domains and the longest common substrings with a certain threshold. 
    SEE --> Khailaile et al._2014_JImmunol
        --> Luo_2015_PNAS
        --> Robert_2018_CurrOpBiotechno
        --> Smith_1999_PNAS

    The matching domains are the string sign (SI),
    the principal matching area (P1, P2, ..., P5),
    and the secondary (S1, S2, ... S4) matching area.
    
     -----------------------------------------------------------------------
     |  Sn  |  P5  |  P4  |  P3  |  P2  |  P1  |  S4  |  S3  |  S2  |  S1  |  
     -----------------------------------------------------------------------
    
    SEE THE MAIN TEXT FOR FURTHER DETAILS
    
    
    
1. The antibody population is independent in bit-map space of the 
    virus-Vic11-cluster covering area 
2. 
3. 
  
   


             Matching length of the principal area
   For a virion-antibody pair, the function returns the match length of the 
   principal area, that is 5 bits of the 10-long bit string. The 5 bits are, 
   from left to right the following after the sign bit (MSB).
   
   Internally, the matching length is weighted by beta (beta1), this parameter
   is equivalent to the bit position of the most left bit of the common 
   substring of the principal matching area.
   Note that due to the threshold (alpha) with value 3, the shortest common 
   substring, at the rightest position, will take the beta value of 3.
     ------------------------------------------------------------------ 
     |     |  5  |  4  |  3  |  P2  |  P1  |     |      |      |      |  
     ------------------------------------------------------------------ 

   input. virion,  antib, alpha (threshold for common substring), 
          nbit_secon: length of the secondary matching area
   output. matching: matching length weighted by beta
"""
#def principal_area_comp(virion, antib, alpha, nbit_secon):
#
#    space =  alpha 
#    beta1 = 0
#    if virion[1:(len(virion)-(nbit_secon-1))] == antib[1:(len(virion)-
#                                                              (nbit_secon-1))]:
#        return (len(virion) - nbit_secon)**2 # Equal Substring 
#    for i in range(space+1):
#        if (virion[(i+1):(i+alpha+1)] == antib[(i+1):(i+alpha+1)]):
#            space = space + 1
#            beta1 = 1
#            
#    space = space - 1    
#    matching = space * beta1
#    return matching


def principal_area_comp(virion, antib, alpha, nbit_secon):

    space =  alpha 
    beta1 = 0
    if virion[1:(len(virion)-(nbit_secon-1))] == antib[1:(len(virion)-
                                                              (nbit_secon-1))]:
        return (len(virion) - nbit_secon) # Equal Substring 
    for i in range(space+1):
        if (virion[(i+1):(i+1+alpha)] == antib[(i+1):(i+1+alpha)]):
            space = space + 1
            beta1 = 1
            
    space = space - 1    
    matching = space * beta1
    return matching



"""              Matching length of the secondary area
   For a virion-antibody pair, the function returns the match length of the 
   secondary area.
   --------------------------------------------------------------------------
   |      |      |      |      |      |     |   4   |   3   |   2   |   1   |  
   -------------------------------------------------------------------------- 

   input. virion,  antib, gamm (threshold for common substring), 
          nbit_secon: length of the principal matching area
   output. matching: matching length
"""
def second_area_comp(virion, antib, gamm, nbit_secon):

    #space = len(virion) - gamm - nbit_secon 
    space = gamm
    area = gamm
    if virion[nbit_secon+1:len(virion)] == antib[nbit_secon+1:len(virion)]:
        return (len(virion) - nbit_secon-1) # Equal Substring         # Equal Substring 
    for i in range(space+1):
        if (virion[(i+nbit_secon+1):(i+nbit_secon+1+gamm)] == 
                                    antib[(i+nbit_secon+1):(i+nbit_secon+1+gamm)]):
            area = area + 1
    area = area - 1
    if area == gamm:
        area = 0
#    if area == 2 and gamm ==3:
#        area = 0
#    if area == 3 and gamm ==4:
#        area = 0
    return area


#______________________________________________________________________________

"""                           MUTATION
   For a given antibody of 2*N long string (x_Anb, y_Anb), a random position
   in the 2*N-long range is selected and changed to its complement value with 
   the "not" logic function.
   N = 10
   
   input. 2*N-long antibody (x_Anb, y_Anb)
   output. mutated 2*N-long antibody (x_Anb, y_Anb)
  
"""

def mutation(x_Anb, y_Anb, nbit):
    point = np.random.randint(2*nbit, size=(1,))  # 20 = 2*N
    point = point[0]
    #print("point: ",point)
    if point < nbit:
        ff = int (not int(x_Anb[point]))
        lis = list(x_Anb)
        lis[point] = str(ff)
        lis[0] = str(0)
        x_Anb = "".join(lis)
    else:
        point = point - nbit
        ff = int (not int(y_Anb[point]))
        lis = list(y_Anb)
        lis[point] = str(ff)
        lis[0] = str(0)
        y_Anb = "".join(lis)
        
    return x_Anb, y_Anb
    
#______________________________________________________________________________


"""                           Half-life
  
"""
def half_life_decay(popul, half_l, stop, step):
    
    event_list_step = []
    event_list_step.append(popul)
    point_clear = np.empty([1, 1])
    point_clear=np.delete(point_clear, 0,0)
    L = np.log(2)/half_l
    prob_dying = 1-np.exp(-L*step)
    
    while step <= stop and popul.size != 0:
        dead_count = 0
        for i in range(len(popul)):
            if random.random() <= prob_dying:
                point_clear = np.concatenate((point_clear, np.array([[i]])), 
                                                                   axis=None)
                dead_count += 1
        if point_clear.size != 0:
            popul = np.delete(popul, point_clear, 0) 
            
        point_clear = np.empty([1, 1])      
        point_clear=np.delete(point_clear, 0,0)
        event_list_step.append(popul)  #  Keeps record of individuals
        
    
        step += 1  
    return event_list_step


"""                           Half-life 2 Only numbers
  
"""
def py_discrete_decay(N_pop, hal_lif, step, time):
    
    event_list_step = []
    event_list_step.append(N_pop)
    point_clear = np.empty([1, 1])
    point_clear = np.delete(point_clear, 0,0)
    quantity = np.empty([1, 1])
    quantity[0,0] = len(N_pop)
    L = np.log(2)/hal_lif
    prob_dying = 1-np.exp(-L*step)
    
    while step <= time and N_pop.size != 0:
        dead_count = 0
        for i in range(len(N_pop)):
            if random.random() <= prob_dying:
                point_clear = np.concatenate((point_clear, np.array([[i]])), 
                                                                   axis=None)
                dead_count += 1
        if point_clear.size != 0:
            N_pop = np.delete(N_pop, point_clear, 0) 
            
        point_clear = np.empty([1, 1])      
        point_clear=np.delete(point_clear, 0,0)
        event_list_step.append(N_pop)  #  Keeps record of individuals
        quantity = np.concatenate((quantity, np.array([[len(N_pop)]])), 
                                                                   axis=None)
        step += 1  
    
    return quantity



#_____________________       EOF      _________________________________________