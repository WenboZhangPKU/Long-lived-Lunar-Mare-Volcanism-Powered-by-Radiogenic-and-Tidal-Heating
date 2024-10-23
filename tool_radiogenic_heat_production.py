#------------------------------------------------------------#
#Different schemes of internal heating for Bulk Silicate Moon
#-----------------------------------------------------------#
import numpy as NP
from scipy.constants import Julian_year

def radio_heat_production(time, model):
    time_billion_year = 4.0 - time/Julian_year/1e9
    lamda =[9.8458e-10,1.5507e-10,5.5452e-10,4.9511e-11]
    H =[5.69e-4,9.46e-5,2.92e-5,2.64e-5] #heat production rate, W/kg   
    HPR_time = 0 # W/kg
    if model=='Zhang2013':
        C1=[25.7*0.0071e-9,25.7*0.9928e-9,0,102.8e-9]
        HPR_time=C1[0]*H[0]*NP.exp(time_billion_year*lamda[0]*1e9)+C1[1]*H[1]*NP.exp(time_billion_year*lamda[1]*1e9) \
                    +C1[2]*H[2]*NP.exp(time_billion_year*lamda[2]*1e9)+C1[3]*H[3]*NP.exp(time_billion_year*lamda[3]*1e9)
        
    elif model=='Laneuville2018':
        C2=[79.5e-9/3.7*0.0071,79.5e-9/3.7*0.9928,79.5e-9/3.7*4500*1.19e-4,79.5e-9]
        HPR_time=C2[0]*H[0]*NP.exp(time_billion_year*lamda[0]*1e9)+C2[1]*H[1]*NP.exp(time_billion_year*lamda[1]*1e9) \
                        +C2[2]*H[2]*NP.exp(time_billion_year*lamda[2]*1e9)+C2[3]*H[3]*NP.exp(time_billion_year*lamda[3]*1e9)
    
    elif model=='Zhang2017':
        C3=[25.7*0.0071e-9,25.7*0.9928e-9,2500*25.7e-9*1.19e-4,102.8e-9]
        HPR_time=C3[0]*H[0]*NP.exp(time_billion_year*lamda[0]*1e9)+C3[1]*H[1]*NP.exp(time_billion_year*lamda[1]*1e9) \
                        +C3[2]*H[2]*NP.exp(time_billion_year*lamda[2]*1e9)+C3[3]*H[3]*NP.exp(time_billion_year*lamda[3]*1e9)
    
    elif model=='TS2002':
        C4=[31.0e-9*0.0071,31.0e-9*0.9928,310e-6*1.19e-4,124e-9]
        HPR_time=(C4[0]*H[0]*NP.exp(time_billion_year*lamda[0]*1e9)+C4[1]*H[1]*NP.exp(time_billion_year*lamda[1]*1e9) \
                +C4[2]*H[2]*NP.exp(time_billion_year*lamda[2]*1e9)+C4[3]*H[3]*NP.exp(time_billion_year*lamda[3]*1e9))

    elif model=='Li2019':
        C5=[25.7*0.0071e-9,25.7*0.9928e-9,2500*102.8e-9*1.19e-4,102.8e-9]
        HPR_time=C5[0]*H[0]*NP.exp(time_billion_year*lamda[0]*1e9)+C5[1]*H[1]*NP.exp(time_billion_year*lamda[1]*1e9) \
                        +C5[2]*H[2]*NP.exp(time_billion_year*lamda[2]*1e9)+C5[3]*H[3]*NP.exp(time_billion_year*lamda[3]*1e9)

    elif model=='MS1995':
        C6=[20.3e-9*0.0071,20.3e-9*0.9928,20.3e-9*11800*1.19e-4,79.5e-9]
        HPR_time=(C6[0]*H[0]*NP.exp(time_billion_year*lamda[0]*1e9)+C6[1]*H[1]*NP.exp(time_billion_year*lamda[1]*1e9) \
                    +C6[2]*H[2]*NP.exp(time_billion_year*lamda[2]*1e9)+C6[3]*H[3]*NP.exp(time_billion_year*lamda[3]*1e9))
        
    elif model=='TS2002noK':
        C7=[31.0e-9*0.0071,31.0e-9*0.9928,0*1.19e-4,124e-9]
        HPR_time=(C7[0]*H[0]*NP.exp(time_billion_year*lamda[0]*1e9)+C7[1]*H[1]*NP.exp(time_billion_year*lamda[1]*1e9) \
                +C7[2]*H[2]*NP.exp(time_billion_year*lamda[2]*1e9)+C7[3]*H[3]*NP.exp(time_billion_year*lamda[3]*1e9))

    elif model=='MS1995noK':
        C8=[20.3e-9*0.0071,20.3e-9*0.9928,0.0,79.5e-9]
        HPR_time=(C8[0]*H[0]*NP.exp(time_billion_year*lamda[0]*1e9)+C8[1]*H[1]*NP.exp(time_billion_year*lamda[1]*1e9) \
                    +C8[2]*H[2]*NP.exp(time_billion_year*lamda[2]*1e9)+C8[3]*H[3]*NP.exp(time_billion_year*lamda[3]*1e9)) 

    elif model=='Taylor2006':
        C9=[27.5e-9*0.0071,27.5e-9*0.9928,27.5e-9*2000*1.19e-4,110e-9]
        HPR_time=(C9[0]*H[0]*NP.exp(time_billion_year*lamda[0]*1e9)+C9[1]*H[1]*NP.exp(time_billion_year*lamda[1]*1e9) \
                    +C9[2]*H[2]*NP.exp(time_billion_year*lamda[2]*1e9)+C9[3]*H[3]*NP.exp(time_billion_year*lamda[3]*1e9)) 
    
    elif model=='TW2014':
        C10=[20.3e-9*0.0071,20.3e-9*0.9928,36.9e-6*1.19e-4,79.5e-9]
        HPR_time=(C10[0]*H[0]*NP.exp(time_billion_year*lamda[0]*1e9)+C10[1]*H[1]*NP.exp(time_billion_year*lamda[1]*1e9) \
                    +C10[2]*H[2]*NP.exp(time_billion_year*lamda[2]*1e9)+C10[3]*H[3]*NP.exp(time_billion_year*lamda[3]*1e9)) 
    
    elif model=='none':
        HPR_time = 0.0
    return HPR_time