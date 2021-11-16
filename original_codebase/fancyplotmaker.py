#IMPORTED CODES AND SETTINGS
import numpy as np
import corner
import matplotlib.pyplot as plt
import csv
import time
import warnings
import xlsxwriter
import json
from numba import jit
import os
import pymultinest
from pymultinest.solve import solve
from scipy.special import erfcinv
from scipy.special import lambertw as W
from scipy import stats as stat
try: os.mkdir('chains')
except OSError: pass
warnings.simplefilter('ignore', np.RankWarning)
warnings.simplefilter('ignore', RuntimeWarning)
warnings.simplefilter('ignore', FutureWarning)
#PERSONAL FAVOURITE GRAPH SETTINGS
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top'] = 'True'
plt.rcParams['ytick.right'] = 'True'
plt.rcParams['xtick.labelsize']= 12
plt.rcParams['ytick.labelsize']= 12
plt.rcParams['text.usetex'] = 'True'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
#Uploading the model stellar catalogue for modelling
stellar = open('StellarCompositionsSortFE.csv')
stellarlist =  [row for row in csv.reader(stellar)]
stellararray = np.asarray(stellarlist)
stellarfloat = stellararray.astype(np.float)
#Uploading WD system names
Names = ['GD61']
N=0
Al_obs=0
Al_err=0
Ti_obs=0
Ti_err=0
Ca_obs=-7.90
Ca_err=0.0634
Ni_obs=0
Ni_err=0
Fe_obs=-7.60
Fe_err=0.0667
Cr_obs=0
Cr_err=0
Mg_obs=-6.69
Mg_err=0.0467
Si_obs=-6.82
Si_err=0.0367
Na_obs=0
Na_err=0
O_obs=-5.95
O_err=0.0434
C_obs=0
C_err=0
Nz_obs=0
Nz_err=0
t_Al=173500																																																					
t_Ti=78200
t_Ca=78200
t_Ni=85500
t_Fe=85500
t_Cr=85500
t_Mg=180800
t_Si=143800
t_Na=180800
t_O=170600
t_C=173000
t_Nz=173000

#UPLOADING THE OBSERVATIONS AND STELLAR CATALOGUE
#Uploading the Observational Compositions in the form Al,Alerror,Ti,Tierror log10 to H or He
observations = np.array([Al_obs,Al_err,Ti_obs,Ti_err,Ca_obs,Ca_err,Ni_obs,Ni_err,Fe_obs,Fe_err,Cr_obs,Cr_err,Mg_obs,Mg_err,Si_obs,Si_err,Na_obs,Na_err,O_obs,O_err,C_obs,C_err,Nz_obs,Nz_err])
#Uploading the Observational Sinking Timescales in years in the form Al,Ti,Ca,Ni,Fe,Cr,Mg,Si,Na,O,C,N
timescales = np.array([t_Al,t_Ti,t_Ca,t_Ni,t_Fe,t_Cr,t_Mg,t_Si,t_Na,t_O,t_C,t_Nz])
#Uploading the Observational Errors for the plot
AlMg_err = ((Al_err)**2 + (Mg_err)**2)**(0.5)
TiMg_err = ((Ti_err)**2 + (Mg_err)**2)**(0.5)
CaMg_err = ((Ca_err)**2 + (Mg_err)**2)**(0.5)
NiMg_err = ((Ni_err)**2 + (Mg_err)**2)**(0.5)
FeMg_err = ((Fe_err)**2 + (Mg_err)**2)**(0.5)
CrMg_err = ((Cr_err)**2 + (Mg_err)**2)**(0.5)
SiMg_err = ((Si_err)**2 + (Mg_err)**2)**(0.5)
NaMg_err = ((Na_err)**2 + (Mg_err)**2)**(0.5)
OMg_err = ((O_err)**2 + (Mg_err)**2)**(0.5)
CMg_err = ((C_err)**2 + (Mg_err)**2)**(0.5)
NzMg_err = ((Nz_err)**2 + (Mg_err)**2)**(0.5)
errors = np.array([AlMg_err,TiMg_err,CaMg_err,NiMg_err,FeMg_err,CrMg_err,SiMg_err,NaMg_err,OMg_err,CMg_err,NzMg_err])
vars()['whitedwarftimescales'+str(0)] = timescales[:]
        
   
vars()['whitedwarf'+str(0)] = observations[:]
        
    
vars()['whitedwarfabundances'+str(0)] = np.transpose([vars()['whitedwarf'+str(0)][0],vars()['whitedwarf'+str(0)][2],vars()['whitedwarf'+str(0)][4],vars()['whitedwarf'+str(0)][6],vars()['whitedwarf'+str(0)][8],vars()['whitedwarf'+str(0)][10],vars()['whitedwarf'+str(0)][12],vars()['whitedwarf'+str(0)][14],vars()['whitedwarf'+str(0)][16],vars()['whitedwarf'+str(0)][18],vars()['whitedwarf'+str(0)][20],vars()['whitedwarf'+str(0)][22]])
vars()['whitedwarferrors'+str(0)] = np.transpose([vars()['whitedwarf'+str(0)][1],vars()['whitedwarf'+str(0)][3],vars()['whitedwarf'+str(0)][5],vars()['whitedwarf'+str(0)][7],vars()['whitedwarf'+str(0)][9],vars()['whitedwarf'+str(0)][11],vars()['whitedwarf'+str(0)][13],vars()['whitedwarf'+str(0)][15],vars()['whitedwarf'+str(0)][17],vars()['whitedwarf'+str(0)][19],vars()['whitedwarf'+str(0)][21],vars()['whitedwarf'+str(0)][23]])
    

vars()['errors'+str(0)] = errors[:]
#COMPRESSING THE OBSERVATIONS
#Compress the observations to only the observed elements in each system
if vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][1],vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][3],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][5],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][8]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][1],vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][3],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][5],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][8]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][1],vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][3],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][5],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][8]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][1],vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][3],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][5],vars()['whitedwarferrors'+str(N)][6]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][1],vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][3],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][5],vars()['whitedwarfabundances'+str(N)][6]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][1],vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][3],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][5],vars()['whitedwarftimescales'+str(N)][6]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][3],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][5],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][8]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][3],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][5],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][8]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][3],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][5],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][8]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][1],vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][5],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][8]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][1],vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][5],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][8]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][1],vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][5],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][8]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][1],vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][3],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][8]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][1],vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][3],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][8]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][1],vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][3],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][8]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][5],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][8]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][5],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][8]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][5],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][8]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][1],vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][8]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][1],vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][8]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][1],vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][8]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][3],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][5],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][8]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][3],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][8]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][3],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][8]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][1],vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][3],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][6]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][1],vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][3],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][6]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][1],vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][3],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][6]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][3],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][5],vars()['whitedwarferrors'+str(N)][6]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][3],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][5],vars()['whitedwarfabundances'+str(N)][6]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][3],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][5],vars()['whitedwarftimescales'+str(N)][6]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][1],vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][5],vars()['whitedwarferrors'+str(N)][6]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][1],vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][5],vars()['whitedwarfabundances'+str(N)][6]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][1],vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][5],vars()['whitedwarftimescales'+str(N)][6]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][3],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][6]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][3],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][6]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][3],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][6]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][8]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][8]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][8]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][5],vars()['whitedwarferrors'+str(N)][6]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][5],vars()['whitedwarfabundances'+str(N)][6]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][5],vars()['whitedwarftimescales'+str(N)][6]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][1],vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][6]])
        whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][1],vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][6]])
        whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][1],vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][6]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][6]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][6]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][6]])
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][0],vars()['whitedwarferrors'+str(N)][1],vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][3],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][5],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][7],vars()['whitedwarferrors'+str(N)][8],vars()['whitedwarferrors'+str(N)][9]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][0],vars()['whitedwarfabundances'+str(N)][1],vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][3],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][5],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][7],vars()['whitedwarfabundances'+str(N)][8],vars()['whitedwarfabundances'+str(N)][9]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][0],vars()['whitedwarftimescales'+str(N)][1],vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][3],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][5],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][7],vars()['whitedwarftimescales'+str(N)][8],vars()['whitedwarftimescales'+str(N)][9]])
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][0],vars()['whitedwarferrors'+str(N)][1],vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][3],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][5],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][7],vars()['whitedwarferrors'+str(N)][9]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][0],vars()['whitedwarfabundances'+str(N)][1],vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][3],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][5],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][7],vars()['whitedwarfabundances'+str(N)][9]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][0],vars()['whitedwarftimescales'+str(N)][1],vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][3],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][5],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][7],vars()['whitedwarftimescales'+str(N)][9]])
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][0],vars()['whitedwarferrors'+str(N)][1],vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][3],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][5],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][7],vars()['whitedwarferrors'+str(N)][8]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][0],vars()['whitedwarfabundances'+str(N)][1],vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][3],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][5],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][7],vars()['whitedwarfabundances'+str(N)][8]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][0],vars()['whitedwarftimescales'+str(N)][1],vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][3],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][5],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][7],vars()['whitedwarftimescales'+str(N)][8]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][1],vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][5],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][7],vars()['whitedwarferrors'+str(N)][8],vars()['whitedwarferrors'+str(N)][9]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][1],vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][5],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][7],vars()['whitedwarfabundances'+str(N)][8],vars()['whitedwarfabundances'+str(N)][9]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][1],vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][5],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][7],vars()['whitedwarftimescales'+str(N)][8],vars()['whitedwarftimescales'+str(N)][9]])
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][0],vars()['whitedwarferrors'+str(N)][1],vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][5],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][7],vars()['whitedwarferrors'+str(N)][9]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][0],vars()['whitedwarfabundances'+str(N)][1],vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][5],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][7],vars()['whitedwarfabundances'+str(N)][9]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][0],vars()['whitedwarftimescales'+str(N)][1],vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][5],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][7],vars()['whitedwarftimescales'+str(N)][9]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][1],vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][3],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][5],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][7],vars()['whitedwarferrors'+str(N)][9]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][1],vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][3],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][5],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][7],vars()['whitedwarfabundances'+str(N)][9]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][1],vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][3],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][5],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][7],vars()['whitedwarftimescales'+str(N)][9]])
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][0],vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][3],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][5],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][7],vars()['whitedwarferrors'+str(N)][9]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][0],vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][3],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][5],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][7],vars()['whitedwarfabundances'+str(N)][9]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][0],vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][3],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][5],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][7],vars()['whitedwarftimescales'+str(N)][9]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][1],vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][3],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][5],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][7]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][1],vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][3],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][5],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][7]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][1],vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][3],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][5],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][7]])
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] == 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][0],vars()['whitedwarferrors'+str(N)][3],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][5],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][7],vars()['whitedwarferrors'+str(N)][9]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][0],vars()['whitedwarfabundances'+str(N)][3],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][5],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][7],vars()['whitedwarfabundances'+str(N)][9]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][0],vars()['whitedwarftimescales'+str(N)][3],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][5],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][7],vars()['whitedwarftimescales'+str(N)][9]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][1],vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][5],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][7],vars()['whitedwarferrors'+str(N)][9]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][1],vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][5],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][7],vars()['whitedwarfabundances'+str(N)][9]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][1],vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][5],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][7],vars()['whitedwarftimescales'+str(N)][9]])
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][0],vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][3],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][7],vars()['whitedwarferrors'+str(N)][9]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][0],vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][3],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][7],vars()['whitedwarfabundances'+str(N)][9]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][0],vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][3],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][7],vars()['whitedwarftimescales'+str(N)][9]])
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][0],vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][7],vars()['whitedwarferrors'+str(N)][9]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][0],vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][7],vars()['whitedwarfabundances'+str(N)][9]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][0],vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][7],vars()['whitedwarftimescales'+str(N)][9]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][7],vars()['whitedwarferrors'+str(N)][9]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][7],vars()['whitedwarfabundances'+str(N)][9]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][7],vars()['whitedwarftimescales'+str(N)][9]])
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] != 0 and vars()['whitedwarfabundances'+str(N)][11] != 0:
          whitedwarferrors = np.transpose([vars()['whitedwarferrors'+str(N)][2],vars()['whitedwarferrors'+str(N)][3],vars()['whitedwarferrors'+str(N)][4],vars()['whitedwarferrors'+str(N)][6],vars()['whitedwarferrors'+str(N)][7],vars()['whitedwarferrors'+str(N)][9],vars()['whitedwarferrors'+str(N)][10],vars()['whitedwarferrors'+str(N)][11]])
          whitedwarfabundances = np.transpose([vars()['whitedwarfabundances'+str(N)][2],vars()['whitedwarfabundances'+str(N)][3],vars()['whitedwarfabundances'+str(N)][4],vars()['whitedwarfabundances'+str(N)][6],vars()['whitedwarfabundances'+str(N)][7],vars()['whitedwarfabundances'+str(N)][9],vars()['whitedwarfabundances'+str(N)][10],vars()['whitedwarfabundances'+str(N)][11]])
          whitedwarftimescales=np.transpose([vars()['whitedwarftimescales'+str(N)][2],vars()['whitedwarftimescales'+str(N)][3],vars()['whitedwarftimescales'+str(N)][4],vars()['whitedwarftimescales'+str(N)][6],vars()['whitedwarftimescales'+str(N)][7],vars()['whitedwarftimescales'+str(N)][9],vars()['whitedwarftimescales'+str(N)][10],vars()['whitedwarftimescales'+str(N)][11]])              
else:
          print("Error: Observations cannot be analysed due to insufficent elemental abundances")
#SETTING UP THE MODEL
#Set up the model to only the observed elements in each system
if vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        planetesimal_abundance0 = np.zeros(7)
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        planetesimal_abundance0 = np.zeros(6)
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        planetesimal_abundance0 = np.zeros(6)
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        planetesimal_abundance0 = np.zeros(6)
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        planetesimal_abundance0 = np.zeros(6)
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        planetesimal_abundance0 = np.zeros(5)
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
       planetesimal_abundance0 = np.zeros(5)
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        planetesimal_abundance0 = np.zeros(5)
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        planetesimal_abundance0 = np.zeros(5)
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        planetesimal_abundance0 = np.zeros(5)
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        planetesimal_abundance0 = np.zeros(5)
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        planetesimal_abundance0 = np.zeros(4)
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        planetesimal_abundance0 = np.zeros(4)
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        planetesimal_abundance0 = np.zeros(4)
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        planetesimal_abundance0 = np.zeros(4)
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        planetesimal_abundance0 = np.zeros(3)
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          planetesimal_abundance0 = np.zeros(10)
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          planetesimal_abundance0 = np.zeros(9)
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          planetesimal_abundance0 = np.zeros(9)
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          planetesimal_abundance0 = np.zeros(8)
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          planetesimal_abundance0 = np.zeros(8)          
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          planetesimal_abundance0 = np.zeros(8)          
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          planetesimal_abundance0 = np.zeros(8)          
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          planetesimal_abundance0 = np.zeros(7)          
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] == 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          planetesimal_abundance0 = np.zeros(7)           
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          planetesimal_abundance0 = np.zeros(7)          
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          planetesimal_abundance0 = np.zeros(7)           
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          planetesimal_abundance0 = np.zeros(6)           
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
          planetesimal_abundance0 = np.zeros(5)           
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] != 0 and vars()['whitedwarfabundances'+str(N)][11] != 0:
          planetesimal_abundance0 = np.zeros(8)       
else:
        print("Error: Stellar Catalogue cannot be reformed")    
if vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ti','Ca','Ni','Fe','Cr','Mg','Na']
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ti','Ca','Ni','Fe','Cr','Mg']
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ca','Ni','Fe','Cr','Mg','Na']
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ti','Ca','Fe','Cr','Mg','Na']
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ti','Ca','Ni','Fe','Mg','Na']
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ca','Fe','Cr','Mg','Na']
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ti','Ca','Fe','Mg','Na']
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ca','Ni','Fe','Mg','Na']
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ti','Ca','Ni','Fe','Mg']
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ca','Ni','Fe','Cr','Mg']
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ti','Ca','Fe','Cr','Mg']
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ca','Ni','Fe','Mg']
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ca','Fe','Mg','Na']
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ca','Fe','Cr','Mg']
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ti','Ca','Fe','Mg']
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ca','Fe','Mg']
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Al','Ti','Ca','Ni','Fe','Cr','Mg','Si','Na','O']
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Al','Ti','Ca','Ni','Fe','Cr','Mg','Si','O']
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Al','Ti','Ca','Ni','Fe','Cr','Mg','Si','Na']
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ti','Ca','Fe','Cr','Mg','Si','Na','O']
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Al','Ti','Ca','Fe','Cr','Mg','Si','O']      
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ti','Ca','Ni','Fe','Cr','Mg','Si','O']        
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Al','Ca','Ni','Fe','Cr','Mg','Si','O']         
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ti','Ca','Ni','Fe','Cr','Mg','Si']       
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] == 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Al','Ni','Fe','Cr','Mg','Si','O']           
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ti','Ca','Fe','Cr','Mg','Si','O']       
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Al','Ca','Ni','Fe','Mg','Si','O']           
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Al','Ca','Fe','Mg','Si','O']          
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        x = ['Ca','Fe','Mg','Si','O']       
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] != 0 and vars()['whitedwarfabundances'+str(N)][11] != 0:
        x = ['Ca','Ni','Fe','Mg','Si','O','C','N']      
else:
        print("Error: Elements cannot be assigned names")  
    #Defining the functions used in the model


@jit(nopython=True)    
def T_disc(d_formation, t_formation):
                #define the constants
                Mstar = 2.34
                s0 = (33*1.496*10**11) 
                sigma = (5.67*10**(-8)) 
                k = (1.38*10**(-23))
                mH = (1.67*10**(-27))
                alpha = 0.01
                gamma = 1.7
                mu = 2.4
                k0 = 0.3
                sun = (2.0e30)
                G = (6.674*10**(-11))
                #stellar evolution tracks from Siess 2000
                Tstar = 259.5822261222*(Mstar**3)-1424.0943845798*(Mstar**2)+2808.5873546861*(Mstar)+2644.9291158227
                Rstar = (6.955*10**8)*(0.1543971235*( Mstar**8) - 1.2349012901*(Mstar**7)+2.7854995*(Mstar**6)+3.0573592544*(Mstar**5)-25.2247969873*(Mstar**4)+47.9947080829*(Mstar**3)-42.5877850958*(Mstar**2)+18.4521562096*(Mstar)-0.3995121438)
                #disc model constants from Chambers 2009
                Te = 1380
                Tvis = (((27*k0)/(64*sigma))**(1/3))*(((alpha*gamma*k)/(mu*mH))**(1/3))*(((7*sun)/(100*3.14*(s0**2)))**(2/3))*(((G*sun)/(s0**3))**(1/6))*(Mstar**(5/6))
                Trad = ((4/7)**(1/4))*(((((Tstar)**8)*((Rstar)**(4))*k)/(G*sun*mu*mH*(s0**(3))))**(1/7))*(Mstar**(-1/7))
                Svis = (7*sun)/(100*3.14*(s0**2))*(Mstar) 
                Se = Svis*((Tvis/Te)**(14/19)) 
                Srad= Svis*(Tvis/Trad)
                tvis = (1/(16*3.14))*((mu*mH)/(alpha*gamma*k))*(sun/10)*(((G*sun)/((s0)**3))**(0.5))*(1/(Svis))*(1/(Tvis))*(Mstar**(3/2)) 
                invtvis = 1/tvis
                #Chambers Disc Model
                if 0 < d_formation < (((s0)*((Se/Svis)**(95/63))*((1+(t_formation*(3.1536*10**13)*(invtvis)))**(-19/36)))/(1.496*10**11)):
                   T = (Tvis**(5/19))*(Te**(14/19))*((((1.496*10**11)*d_formation)/s0)**(-9/38))*((1+(t_formation*(3.1536*10**13)*(invtvis)))**(-1/8)) 
                elif (((s0)*((Se/Svis)**(95/63))*((1+(t_formation*(3.1536*10**13)*(invtvis)))**(-19/36)))/(1.496*10**11)) < d_formation < (((s0)*((Srad/Svis)**(70/33))*((1+(t_formation*(3.1536*10**13)*(invtvis)))**(-133/132)))/(1.496*10**11)):
                   T = Tvis*((((1.496*10**11)*d_formation)/s0)**(-9/10))*((1+((t_formation*(3.1536*10**13))/(tvis)))**(-19/40))
                elif d_formation > (((s0)*((Srad/Svis)**(70/33))*((1+(t_formation*(3.1536*10**13)*(invtvis)))**(-133/132)))/(1.496*10**11)):
                   T = Trad*(((((1.496*10**11)*d_formation)/s0)**(-3/7))) 
                else:
                   T = 5000
                return T
            
@jit(nopython=True)
def S_disc(t_formation):
                #define the constants
                Mstar = 2.34
                s0 = (33*1.496*10**11) 
                sigma = (5.67*10**(-8)) 
                k = (1.38*10**(-23))
                mH = (1.67*10**(-27))
                alpha = 0.01
                gamma = 1.7
                mu = 2.4
                k0 = 0.3
                sun = (2.0e30)
                G = (6.674*10**(-11))
                #stellar evolution tracks from Siess 2000
                Tstar = 259.5822261222*(Mstar**3)-1424.0943845798*(Mstar**2)+2808.5873546861*(Mstar)+2644.9291158227
                Rstar = (6.955*10**8)*(0.1543971235*( Mstar**8) - 1.2349012901*(Mstar**7)+2.7854995*(Mstar**6)+3.0573592544*(Mstar**5)-25.2247969873*(Mstar**4)+47.9947080829*(Mstar**3)-42.5877850958*(Mstar**2)+18.4521562096*(Mstar)-0.3995121438)
                #disc model constants from Chambers 2009
                Tvis = (((27*k0)/(64*sigma))**(1/3))*(((alpha*gamma*k)/(mu*mH))**(1/3))*(((7*sun)/(100*3.14*(s0**2)))**(2/3))*(((G*sun)/(s0**3))**(1/6))*(Mstar**(5/6))
                Trad = ((4/7)**(1/4))*(((((Tstar)**8)*((Rstar)**(4))*k)/(G*sun*mu*mH*(s0**(3))))**(1/7))*(Mstar**(-1/7))
                Svis = (7*sun)/(100*3.14*(s0**2))*(Mstar) 
                Srad= Svis*(Tvis/Trad)
                tvis = (1/(16*3.14))*((mu*mH)/(alpha*gamma*k))*(sun/10)*(((G*sun)/((s0)**3))**(0.5))*(1/(Svis))*(1/(Tvis))*(Mstar**(3/2)) 
                invtvis = 1/tvis        
                t1 = tvis*(((Tvis/Trad)**(112/73))-1)
                r1 = s0*((Srad/Svis)**(70/33))*((1+t1*invtvis)**(-133/132))
                T1 = Trad*(r1/s0)**(-3/7)
                S1 = Srad*((r1/s0)**(-15/14))*((1+t1*invtvis)**(-19/16));
                Mdot1 = ((0.3*sun*Mstar*invtvis)/16)*(T1/Tvis)*(S1/Svis)*((r1/s0)**1.5)
                trad = (0.7*Mstar*sun*((Trad/Tvis)**(21/73)))/(13*Mdot1)
                invtrad = 1/trad;
                #Chambers Disc Model
                if 0 <= t_formation*(3.1536*10**13) <= t1:
                   s = (s0*(1+(t_formation*(3.1536*10**13)*invtvis))**(6/16))            
                elif t_formation*(3.1536*10**13) > t1:
                   s = (s0*((Tvis/Trad)**(42/73))*((1+((t_formation*(3.1536*10**13)-t1)*invtrad))**(14/13)))
                else:
                   s = s0
                return s/(1.496*10**11)   
    
@jit(nopython=True)
def frange(x, y, step):
            while x < y:
                yield x
                x += step
      
@jit(nopython=True)
def erfcc(x):
            z = abs(x)
            t = 1. / (1. + 0.5*z)
            r = t * np.exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+
                t*(.09678418+t*(-.18628806+t*(.27886807+
                t*(-1.13520398+t*(1.48851587+t*(-.82215223+
                t*.17087277)))))))))
            if (x >= 0.):
                 return r
            else:
                 return 2. - r
    
@jit(nopython=True)    
def normcdf(x, mu, sigma):
            t = x-mu;
            y = 0.5*erfcc(-t/(sigma*((2.0)**0.5)));
            if y>1.0:
               y = 1.0;
            return y    
    
@jit(nopython=True)
def fznd(d_formation,z_formation):
            x = np.zeros(601)
            for i in frange(d_formation-3,d_formation+0.01,0.01):
                x[int(round((i-(d_formation-3))/(0.01)))] = normcdf(i+0.005, d_formation, z_formation) - normcdf(i-0.005, d_formation, z_formation)
            for i in frange(d_formation+0.01,d_formation+3.01,0.01):
                x[int(round((i-(d_formation-3))/(0.01)))] = x[int(round(((d_formation+3)-i)/(0.01)))]
            return x
    
    
@jit(nopython=True)
def Al_c(T,t_formation):
                 def Al_cond_a(x):
                     Al_cond_a = (3.10807552511606477612696196001603328486749024117674622402773820795118808746337890625E-14)*x**(6)+(-2.98880529944537063480852551587588043779941671118649537675082683563232421875E-10)*x**(5)+(0.000001195980623269447509876800377071059955369491945020854473114013671875)*x**(4)+(-0.0025493162128629974494398169326814240776002407073974609375)*x**(3)+(3.05321509324545647956483662710525095462799072265625)*x**(2)+(-1948.202926753468318565865047276020050048828125)*x+(517456.53613930637948215007781982421875)
                     return Al_cond_a
                 if T - (16.04*np.log(t_formation+0.01058)+72.96259) <= 1500:
                    Al_cond = 1
                 elif 1500 < T - (16.04*np.log(t_formation+0.01058)+72.96259) < 1707:    
                    Al_cond = Al_cond_a(T - (16.04*np.log(t_formation+0.01058)+72.96259))
                 else:
                    Al_cond = 0
                 return Al_cond
    
@jit(nopython=True)
def Al(d_formation,z_formation,t_formation):
            if z_formation <= 0:
               Al_final = Al_c(T_disc(d_formation, t_formation),t_formation)
            else:
              Al_f = np.zeros(601)
              x = fznd(d_formation,z_formation)
              for j in frange (d_formation-3,d_formation+3.01,0.01):
                  Al_f[int(round((j-(d_formation-3))/(0.01)))] = x[int(round((j-(d_formation-3))/(0.01)))]*Al_c(T_disc(j, t_formation),t_formation)
              Al_final = 0
              for i in range(0,len(Al_f)): 
                  Al_final = Al_final + Al_f[i]
            return Al_final
        
@jit(nopython=True)
def Ti_c(T,t_formation):
                 def Ti_cond_a(x):
                     Ti_cond_a = ((-3.3156670111210073797065924276213958720032300478940123777960025108768604695796966552734375E-17)*x**8)+((4.203993324422311826043359377802558009035345520931770124661852605640888214111328125E-13)*x**7)+((-2.3308242069385726646509548153074486975810941657982766628265380859375E-9)*x**6)+((0.00000738067440514030478622246878028789751624572090804576873779296875)*x**5)+((-0.01459951816407163531497115371848849463276565074920654296875)*x**4)+((18.472934446428066479484186857007443904876708984375)*x**3)+((-14601.158369747523465775884687900543212890625)*x**2)+((6591344.180495082400739192962646484375)*x)-1301102300.1227576732635498046875 
                     return Ti_cond_a
                 if T - (16.04*np.log(t_formation+0.0087753)+75.96259) <= 1468:
                    Ti_cond = 1
                 elif 1468 < T - (16.04*np.log(t_formation+0.0087753)+75.96259) < 1674:    
                    Ti_cond = Ti_cond_a(T - (16.04*np.log(t_formation+0.0087753)+75.96259))
                 else:
                    Ti_cond = 0
                 return Ti_cond
          
@jit(nopython=True)  
def Ti(d_formation, z_formation, t_formation): 
            if z_formation <= 0:
               Ti_final = Ti_c(T_disc(d_formation, t_formation),t_formation)
            else:
                Ti_f = np.zeros(601)
                x = fznd(d_formation,z_formation)
                for j in frange (d_formation-3,d_formation+3.01,0.01):
                    Ti_f[int(round((j-(d_formation-3))/(0.01)))] = x[int(round((j-(d_formation-3))/(0.01)))]*Ti_c(T_disc(j, t_formation),t_formation)
                Ti_final = 0
                for i in range(0,len(Ti_f)):
                    Ti_final = Ti_final + Ti_f[i]
            return Ti_final
                
@jit(nopython=True)        
def Ca_c(T,t_formation):
                 def Ca_cond_a(x):
                     Ca_cond_a = ((1.15409860888684715070558147524663013805356300913418380133063043806756872836490644901671900212081144326652093479168570411275140941143035888671875E-41)*x**16) + ((-7.171665870324334403620738008342045875966057968309169581629460765338530039275550912576033579768383478814097742315425421111285686492919921875E-38)*x**15) + ((1.2486097871082906305004603787250304179946426439046147644930143684236979496464802796533845145876551185892822104506194591522216796875E-34)*x**14) + ((4.09425585589608066200986512139754006976727890677916899640251763168390083934020397238816302287744974819361232221126556396484375E-32)*x**13) + ((-1.97313282299698333438597679824231542449623624115544106697797917372750001077479076183607276107068173587322235107421875E-28)*x**12) + ((-1.8485424392426202673313455796478418847601333433741379248512214217140849081832953970661037601530551910400390625E-25)*x**11) + ((2.172565761056839443273346761743906029390306724271352062374741860584226316177591797895729541778564453125E-22)*x**10) + ((5.5671791492154865113878770421949079644284899772400877944467101343661852297373116016387939453125E-19)*x**9) + ((1.43556093043120870651237403236267293930612860722133283797319336372311227023601531982421875E-16)*x**8) + ((-9.665669613106708433014416674686200754186560235581282540806569159030914306640625E-13)*x**7) + ((-1.2936014489262972379678462677885485143658428341950639151036739349365234375E-9)*x**6) + ((8.200728257700856034889355140882205574826002703048288822174072265625E-7)*x**5) + ((0.0034863053420420317234096341252325146342627704143524169921875)*x**4) + ((-0.0092027180751187119545075887572238571010529994964599609375)*x**3) + ((-8525.7967774152275524102151393890380859375)*x**2) + ((8250207.196432891301810741424560546875)*x) - 2402808821.69645404815673828125
                     return Ca_cond_a
                 def Ca_cond_b(x):
                     Ca_cond_b = ((1.0819736872713684141293723526945425982857134942566375098431792673831453966839578038452529727585357786490878888481514952246698157978244125843048095703125E-43)*x**16) + ((-1.00121002204770569644344717507797641907136690805543250751134176009899429725778395239087485103271537027931625818411021100473590195178985595703125E-39)*x**15) + ((3.143752730835549102437107619052008930266790777750619857771379440303695646789305264916830249220802695475640575750730931758880615234375E-36)*x**14) + ((-2.23858730393389934223833753257197491496828266028773194063879543211711744845308910415913015989897161261978908441960811614990234375E-33)*x**13) + ((-6.256171036799533489090208823588447925165592004448965596623878159027800405740725279979397299712218227796256542205810546875E-30)*x**12) + ((4.82504655593846055647991338708181686010254373981655842588378357329765319487513419716151474858634173870086669921875E-27)*x**11) + ((1.794041001768184297981856119470577692687516066053918472476805712374769985473221822758205235004425048828125E-23)*x**10) + ((-1.69323485681020684065239217360233171851616412602120902665598679848191210339791723527014255523681640625E-21)*x**9) + ((-4.919640238019071843730290295297970651334654583236756508757281380894710309803485870361328125E-17)*x**8) + ((-3.1548596414976621099022750727316664358574650750544066113434382714331150054931640625E-14)*x**7) + ((1.1457545307366237122341474325839063598542200139718261198140680789947509765625E-10)*x**6) + ((1.5593163186288970291031963193162379610612333635799586772918701171875E-7)*x**5) + ((-0.00027866287790703596349839443746532197110354900360107421875)*x**4) + ((-0.46333233646740412670084197088726796209812164306640625)*x**3) + ((1219.143638702181988264783285558223724365234375)*x**2) + ((-932219.0221518310718238353729248046875)*x) + 250622657.8589092195034027099609375
                     return Ca_cond_b
                 if T - (16.95*np.log(t_formation+0.0166)+36.43) <= 1325:
                    Ca_cond = 1
                 elif 1325 < T - (16.95*np.log(t_formation+0.0166)+36.43) <= 1459:    
                    Ca_cond = Ca_cond_a(T - (16.95*np.log(t_formation+0.0166)+36.43))
                 elif 1459 < T - (16.95*np.log(t_formation+0.0166)+36.43) < 1742:
                    Ca_cond = Ca_cond_b(T - (16.95*np.log(t_formation+0.0166)+36.43))
                 else:
                    Ca_cond = 0
                 return Ca_cond
    
@jit(nopython=True)
def Ca(d_formation, z_formation, t_formation):        
            if z_formation <= 0:
               Ca_final = Ca_c(T_disc(d_formation, t_formation),t_formation)
            else:
                Ca_f = np.zeros(601)
                x = fznd(d_formation,z_formation)
                for j in frange (d_formation-3,d_formation+3.01,0.01):
                    Ca_f[int(round((j-(d_formation-3))/(0.01)))] = x[int(round((j-(d_formation-3))/(0.01)))]*Ca_c(T_disc(j, t_formation),t_formation)
                Ca_final = 0
                for i in range(0,len(Ca_f)):
                    Ca_final = Ca_final + Ca_f[i]
            return Ca_final
    
@jit(nopython=True)
def Ni_c(T,t_formation):
                 def Ni_cond_a(x):
                     Ni_cond_a = ((-0.0001823061527836233983342062447974285532836802303791046142578125)*x**2)+((0.46120119247048896315988031346932984888553619384765625)*x) - 290.6519019506703216393361799418926239013671875
                     return Ni_cond_a
                 if T - (12.0530*np.log(t_formation+0.0115)-0.1378177) <= 1280:
                    Ni_cond = 1
                 elif 1280 < T - (12.0530*np.log(t_formation+0.0115)-0.1378177) < 1340.2:    
                    Ni_cond = Ni_cond_a(T - (12.0530*np.log(t_formation+0.0115)-0.1378177))
                 else:
                    Ni_cond = 0
                 return Ni_cond
    
@jit(nopython=True)
def Ni(d_formation, z_formation, t_formation):        
            if z_formation <= 0:
               Ni_final = Ni_c(T_disc(d_formation, t_formation),t_formation)
            else:
                Ni_f = np.zeros(601)
                x = fznd(d_formation,z_formation)
                for j in frange (d_formation-3,d_formation+3.01,0.01):
                    Ni_f[int(round((j-(d_formation-3))/(0.01)))] = x[int(round((j-(d_formation-3))/(0.01)))]*Ni_c(T_disc(j, t_formation),t_formation)
                Ni_final = 0
                for i in range(0,len(Ni_f)):
                    Ni_final = Ni_final + Ni_f[i]
            return Ni_final
     
@jit(nopython=True)    
def Fe_c(T,t_formation):
                 def Fe_cond_a(x):
                     Fe_cond_a = ((5.24668541798360132296085575751383166015658077174776963147451169788837432861328125E-13)*x**6) + ((-3.806244438155727007635671112727171472300824461854062974452972412109375E-9)*x**5) + ((0.00001149918229030552719390108340480338711131480522453784942626953125)*x**4) + ((-0.0185188563629463016912968242877468583174049854278564453125)*x**3) + ((16.767539428254625732961358153261244297027587890625)*x**2) + ((-8093.1022260413728872663341462612152099609375)*x) + 1626848.51622126321308314800262451171875
                     return Fe_cond_a
                 if T - (14.73911*np.log(t_formation+0.018)+31.4985195) <= 1132:
                    Fe_cond = 1
                 elif 1132 < T - (14.73911*np.log(t_formation+0.018)+31.4985195) < 1302.8:    
                    Fe_cond = Fe_cond_a(T - (14.73911*np.log(t_formation+0.018)+31.4985195))
                 else:
                    Fe_cond = 0
                 return Fe_cond
        
@jit(nopython=True)    
def Fe(d_formation, z_formation, t_formation):        
            if z_formation <= 0:
               Fe_final = Fe_c(T_disc(d_formation, t_formation),t_formation)
            else:
                Fe_f = np.zeros(601)
                x = fznd(d_formation,z_formation)
                for j in frange (d_formation-3,d_formation+3.01,0.01):
                    Fe_f[int(round((j-(d_formation-3))/(0.01)))] = x[int(round((j-(d_formation-3))/(0.01)))]*Fe_c(T_disc(j, t_formation),t_formation)
                Fe_final = 0
                for i in range(0,len(Fe_f)):
                    Fe_final = Fe_final + Fe_f[i]
            return Fe_final
    
    
@jit(nopython=True)   
def Cr_c(T,t_formation):
                 def Cr_cond_a(x):
                     Cr_cond_a = ((-7.1438142274896419837267507710285574721426581079608553181969909928739070892333984375E-14)*x**6) + ((5.174629256034451019490644086275092650151208317765849642455577850341796875E-10)*x**5) + ((-0.000001558660606587738834392410647300408754745149053633213043212890625)*x**4) + ((0.0024990965912296153085547789629572434932924807071685791015625)*x**3) + ((-2.249715867410567060602488709264434874057769775390625)*x**2) + ((1078.186945705621610613889060914516448974609375)*x) - 214933.2085798117332160472869873046875
                     return Cr_cond_a
                 if T - (14.23915*np.log(t_formation+0.012)+31.1731504) <= 1100:
                    Cr_cond = 1
                 elif 1100 < T - (14.23915*np.log(t_formation+0.012)+31.1731504) < 1302.8:    
                    Cr_cond = Cr_cond_a(T - (14.23915*np.log(t_formation+0.012)+31.1731504))
                 else:
                    Cr_cond = 0
                 return Cr_cond
       
@jit(nopython=True)
def Cr(d_formation, z_formation, t_formation):        
            if z_formation <= 0:
               Cr_final = Cr_c(T_disc(d_formation, t_formation),t_formation)
            else:
                Cr_f = np.zeros(601)
                x = fznd(d_formation,z_formation)
                for j in frange (d_formation-3,d_formation+3.01,0.01):
                    Cr_f[int(round((j-(d_formation-3))/(0.01)))] = x[int(round((j-(d_formation-3))/(0.01)))]*Cr_c(T_disc(j, t_formation),t_formation)
                Cr_final = 0
                for i in range(0,len(Cr_f)):
                    Cr_final = Cr_final + Cr_f[i]
            return Cr_final
    
@jit(nopython=True)
def Na_c(T,t_formation):
                 def Na_cond_a(x):
                     Na_cond_a = ((-1.2021554780914999564779079206945156923559171746607177055921056307852268218994140625E-15)*x**6) + ((7.187263925563124665321723856900947946350910466861705572227947413921356201171875E-12)*x**5) + ((-1.7970122907779646145371015854806662215281676253653131425380706787109375E-8)*x**4) + ((0.00002405558262328599038811123567160876746129360981285572052001953125)*x**3) + ((-0.01817489681265127554610216975561343133449554443359375)*x**2) + ((7.337489774058976621518013416789472103118896484375)*x) - 1232.983219250172396641573868691921234130859375
                     return Na_cond_a
                 if T - (9.89878457*np.log(t_formation+0.0135)+21.539284458524076) <= 775:
                    Na_cond = 1
                 elif 775 < T - (9.89878457*np.log(t_formation+0.0135)+21.539284458524076) < 1154:    
                    Na_cond = Na_cond_a(T-(9.89878457*np.log(t_formation+0.0135)+21.539284458524076))
                 else:
                    Na_cond = 0
                 return Na_cond
        
@jit(nopython=True)    
def Na(d_formation, z_formation, t_formation):        
            if z_formation <= 0:
               Na_final = Na_c(T_disc(d_formation, t_formation),t_formation)
            else:
                Na_f = np.zeros(601)
                x = fznd(d_formation,z_formation)
                for j in frange (d_formation-3,d_formation+3.01,0.01):
                    Na_f[int(round((j-(d_formation-3))/(0.01)))] = x[int(round((j-(d_formation-3))/(0.01)))]*Na_c(T_disc(j, t_formation),t_formation)
                Na_final = 0
                for i in range(0,len(Na_f)):
                    Na_final = Na_final + Na_f[i]
            return Na_final
    
    
@jit(nopython=True)
def Mg_c(T,t_formation):
                 def Mg_cond_a(x):
                     Mg_cond_a = ((-5.52168686722831612357559355459842953366134121750974372844211757183074951171875E-13)*x**6) + ((4.31852848409746917802742304805761752728443525484181009232997894287109375E-9)*x**5) + ((-0.00001403729180829655235674714719440459020916023291647434234619140625)*x**4) + ((0.0242775491483312004514782955766349914483726024627685546875)*x**3) +((-23.5665589121942247174956719391047954559326171875)*x**2) + ((12175.820469518892423366196453571319580078125)*x) - 2616129.30626156367361545562744140625
                     return Mg_cond_a
                 def Mg_cond_b(x):
                     Mg_cond_b = ((4.405059759212852885170124206340963349498411616433912740831146948039531707763671875E-13)*x**6) + ((-3.66542360387390544713185922616645318061756597671774215996265411376953125E-9)*x**5) + ((0.00001270797581448418760337315536190772036206908524036407470703125)*x**4) + ((-0.0234973480078012504634887847032587160356342792510986328125)*x**3) +((24.43857328915293436466527055017650127410888671875)*x**2) + ((-13555.7660498977129464037716388702392578125)*x) + 3132954.10722694732248783111572265625
                     return Mg_cond_b
                 def Mg_cond_c(x):
                     Mg_cond_c = ((4.0166182802625228425363953713435463628911967931024717959331837846548296511173248291015625E-17)*x**6) + ((-3.677178849584043798790532709340448760970167241257655632580281235277652740478515625E-13)*x**5) + ((1.402579102607708028199615669083001601169513605782412923872470855712890625E-9)*x**4) + ((-0.00000285305722156571719519174874790667928436960210092365741729736328125)*x**3) +((0.0032643007740747824983740255078146219602786004543304443359375)*x**2) + ((-1.9917968667421785955440327597898431122303009033203125)*x) + 506.3697455135483096455573104321956634521484375
                     return Mg_cond_c
                 def Mg_cond_d(x):
                     Mg_cond_d = ((-4.322355443821983153320469836534448750769198878764987769261762817762928534648381173610687255859375E-20)*x**6) + ((4.222125961795620550880260805717440128507087697233568501786749038728885352611541748046875E-16)*x**5) + ((-1.7137119394006084939594905812348525999712356426840642598108388483524322509765625E-12)*x**4) + ((3.69832045126722481988669745888408113554390865829191170632839202880859375E-9)*x**3) +((-0.00000447383939882000267423793615773064402674208395183086395263671875)*x**2) + ((0.0028749438505247700004208155633023125119507312774658203125)*x) - 0.766266938315708134865644751698710024356842041015625 
                     return Mg_cond_d                            
                 if T - (11.47587238*np.log(t_formation+0.016)+24.72092363391385) <= 1192:
                    Mg_cond = 1
                 elif 1192 < T - (11.47587238*np.log(t_formation+0.016)+24.72092363391385) <= 1328:    
                    Mg_cond = Mg_cond_a(T-(11.47587238*np.log(t_formation+0.016)+24.72092363391385))
                 elif 1328 < T - (11.47587238*np.log(t_formation+0.016)+24.72092363391385) <= 1418:    
                    Mg_cond = Mg_cond_b(T-(11.47587238*np.log(t_formation+0.016)+24.72092363391385))
                 elif 1418 < T - (11.47587238*np.log(t_formation+0.016)+24.72092363391385) < 1589:    
                    Mg_cond = Mg_cond_c(T-(11.47587238*np.log(t_formation+0.016)+24.72092363391385))   
                 elif 1589 < T - (11.47587238*np.log(t_formation+0.016)+24.72092363391385) <= 1742:    
                    Mg_cond = Mg_cond_d(T-(11.47587238*np.log(t_formation+0.016)+24.72092363391385))
                 else:
                    Mg_cond = 0.000000001
                 return Mg_cond
        
@jit(nopython=True)     
def Mg(d_formation, z_formation, t_formation):        
            if z_formation <= 0:
               Mg_final = Mg_c(T_disc(d_formation, t_formation),t_formation)
            else:
                Mg_f = np.zeros(601)
                x = fznd(d_formation,z_formation)
                for j in frange (d_formation-3,d_formation+3.01,0.01):
                    Mg_f[int(round((j-(d_formation-3))/(0.01)))] = x[int(round((j-(d_formation-3))/(0.01)))]*Mg_c(T_disc(j, t_formation),t_formation)
                Mg_final = 0
                for i in range(0,len(Mg_f)):
                    Mg_final = Mg_final + Mg_f[i]
            return Mg_final
           
    
@jit(nopython=True)
def Si_c(T,t_formation):
                 def Si_cond_a(x):
                     Si_cond_a = ((3.75405933875938719382681288727743669357805723768392589756084430572185967598698719260852385559701360762119293212890625E-29)*x**12) + ((-3.0073884486408280684201364634866873033452455776110822153585053870977240875195235503269941546022891998291015625E-25)*x**11) + ((9.0048634903801527390104656965941521071313378107324833083649658072022958776869927532970905303955078125E-22)*x**10) + ((-1.0055507157322091416017799916894941717250070124803286615577935236842677113600075244903564453125E-18)*x**9) + ((-4.68256797790428830541734175072154763543180191164017056593138477182947099208831787109375E-16)*x**8) + ((1.733165623947980018272931314359921518321316913358032252290286123752593994140625E-12)*x**7) + ((3.1183996747164098176358373737805489145369364223370212130248546600341796875E-10)*x**6) + ((-0.0000029352518505059235843558142209364092423129477538168430328369140625)*x**5) + ((0.00034179871225393197720077154144746600650250911712646484375)*x**4) + ((5.08636370263309256500861010863445699214935302734375)*x**3) + ((-6179.9044497130735180689953267574310302734375)*x**2) + ((3055402.305034776218235492706298828125)*x) - 576895389.895605564117431640625  
                     return Si_cond_a
                 def Si_cond_b(x):
                     Si_cond_b = ((4.5563314260167819619266633277951851567804970729336905629751608949059594943164185065764826276790699921548366546630859375E-30)*x**12) + ((-3.25127374190906226743013716982109707247917689589838406869117995757491844377451428726999438367784023284912109375E-26)*x**11) + ((7.476737884724708011238013267493669396758286267733162504061338038408024431191734038293361663818359375E-23)*x**10) + ((-2.10773696020689016862573072502591138012545328834638333937064047507448094620485790073871612548828125E-20)*x**9) + ((-1.24373125490839847554958172631850555723606515452946130739775298934546299278736114501953125E-16)*x**8) + ((2.07715768792233028885973054561482808302928619592631065415844204835593700408935546875E-14)*x**7) + ((2.534289823619733228017312305808961259145917210844345390796661376953125E-10)*x**6) + ((5.730895024525212758937972698687624006907981311087496578693389892578125E-8)*x**5) + ((-0.000524995631139943581822127072200601105578243732452392578125)*x**4) + ((-0.2125910157539518952063417600584216415882110595703125)*x**3) + ((1358.807431324278468309785239398479461669921875)*x**2) + ((-1194830.141528136096894741058349609375)*x) + 340880684.0765764713287353515625
                     return Si_cond_b
                 def Si_cond_c(x):
                     Si_cond_c = ((-7.87547923713578326849266449590958455509684883768230522434858118992253209215014151818568498286676771158454357646405696868896484375E-34)*x**12) + ((4.44894114469649605757060633203539537727562066224864568057749877616261289568311103337361345211320440284907817840576171875E-30)*x**11) + ((-5.8916600511480599732198835114156212439868172205058274690367025431289409176383031763180042617022991180419921875E-27)*x**10) + ((-6.3353218927850779094225301869976700985626828307228857232357693688752642291461825152509845793247222900390625E-24)*x**9) + ((8.82176972658880308528122015502828427699107877519698532424452519506985481712035834789276123046875E-21)*x**8) + ((1.99891589419553941703297795942724007838759681088790232106333633055328391492366790771484375E-17)*x**7) + ((-7.441720401059376579603569265131483356666005041062106339921911057899706065654754638671875E-16)*x**6) + ((-4.74392713046148624602930991294336311725476917899868567474186420440673828125E-11)*x**5) + ((-4.526420380549866291779000028634716290554251827416010200977325439453125E-8)*x**4) + ((0.00008453614336463759009430984914246209882549010217189788818359375)*x**3) + ((0.1692234652417285156733584017274552024900913238525390625)*x**2) + ((-310.4804291826909548035473562777042388916015625)*x) + 130553.891433640732429921627044677734375 
                     return Si_cond_c
                 def Si_cond_d(x):
                     Si_cond_d = ((5.9421496515397558612364974277208047379685395956736740987300488448164952901917816887798651117984227720825174401397816836833953857421875E-36)*x**12) + ((-3.4010726649437779111754225605654327248814506563438613654728528349375717146583703097244277824273694932344369590282440185546875E-32)*x**11) + ((4.4535528205847168260778797709914591341048848370527784935661703784025787942134864960674889289293787442147731781005859375E-29)*x**10) + ((5.3813161468261490299537846963942239916043596734736858222388237996320660416760262023672112263739109039306640625E-26)*x**9) + ((-6.8667253731837614112368643141385226348901829403198956068172793647252927229374108719639480113983154296875E-23)*x**8) + ((-1.7517514700538795571877282752353191642183546780363084978404675240426513482816517353057861328125E-19)*x**7) + ((-8.031445720278585992155535549636169995640615367447912031906298579997383058071136474609375E-18)*x**6) + ((4.261063070826374874595249039270690268239351328194430834628292359411716461181640625E-13)*x**5) + ((4.56624329609124617133834951308878936926394231932135880924761295318603515625E-10)*x**4) + ((-7.787798610968512971235101945144396751175008830614387989044189453125E-7)*x**3) + ((-0.001761966259236581212876959767754669883288443088531494140625)*x**2) + ((3.184030127967345524808706613839603960514068603515625)*x) - 1352.292046121405519443214870989322662353515625
                     return Si_cond_d 
                 def Si_cond_e(x):
                     Si_cond_e = ((3.4988671676416799638211296425745044932671150254426578747878778853216885382137769909995363308888825792930532543323352001607418060302734375E-37)*x**12) + ((-2.94809050133692169787083747785070161610149730412277859145594715744668020993842638219502826080198332192594534717500209808349609375E-33)*x**11) + ((8.0225610820918431954301653174759149887550729948478202984966663012636892010902463223676051029542577452957630157470703125E-30)*x**10) + ((-2.77994962557293169052993397868117179035169724293219257854241343257858097050283507911672131740488111972808837890625E-27)*x**9) + ((-1.83865513199147077725657600868732585845843121120492214907237251732274874171935152844525873661041259765625E-23)*x**8) + ((3.960842381981730639304219211609898182470268478819183980386975374443636610521934926509857177734375E-21)*x**7) + ((5.203102948297463771313907760808190830847746143061893775438875309191644191741943359375E-17)*x**6) + ((1.29110718171406867211619895738518398224593981826746613705836352892220020294189453125E-14)*x**5) + ((-1.497994383385505629543742641879614484190508250094353570602834224700927734375E-10)*x**4) + ((-6.82886903037696894757348155431675440496519513544626533985137939453125E-8)*x**3) + ((0.00053456709545030294157552663847354779136367142200469970703125)*x**2) + ((-0.5566119858692413213674399230512790381908416748046875)*x) + 187.67455451644627828500233590602874755859375
                     return Si_cond_e
                 if T - (10.50780681*np.log(t_formation+0.016)+22.6365505802529) <= 1141:
                    Si_cond = 1
                 elif 1141 < T - (10.50780681*np.log(t_formation+0.016)+22.6365505802529) <= 1367:    
                    Si_cond = Si_cond_a(T-(10.50780681*np.log(t_formation+0.016)+22.6365505802529))
                 elif 1367 < T - (10.50780681*np.log(t_formation+0.016)+22.6365505802529) <= 1483:    
                    Si_cond = Si_cond_b(T-(10.50780681*np.log(t_formation+0.016)+22.6365505802529))
                 elif 1483 < T - (10.50780681*np.log(t_formation+0.016)+22.6365505802529) < 1549:    
                    Si_cond = Si_cond_c(T-(10.50780681*np.log(t_formation+0.016)+22.6365505802529))   
                 elif 1549 < T - (10.50780681*np.log(t_formation+0.016)+22.6365505802529) <= 1615:    
                    Si_cond = Si_cond_d(T-(10.50780681*np.log(t_formation+0.016)+22.6365505802529))
                 elif 1615 < T - (10.50780681*np.log(t_formation+0.016)+22.6365505802529) <= 1733:
                    Si_cond = Si_cond_e(T-(10.50780681*np.log(t_formation+0.016)+22.6365505802529))
                 else:
                    Si_cond = 0
                 return Si_cond
         
@jit(nopython=True)    
def Si(d_formation, z_formation, t_formation):        
            if z_formation <= 0:
               Si_final = Si_c(T_disc(d_formation, t_formation),t_formation)
            else:
                Si_f = np.zeros(601)
                x = fznd(d_formation,z_formation)
                for j in frange (d_formation-3,d_formation+3.01,0.01):
                    Si_f[int(round((j-(d_formation-3))/(0.01)))] = x[int(round((j-(d_formation-3))/(0.01)))]*Si_c(T_disc(j, t_formation),t_formation)
                Si_final = 0
                for i in range(0,len(Si_f)):
                    Si_final = Si_final + Si_f[i]
            return Si_final 
    
@jit(nopython=True)
def C_c(T,t_formation):
                 def C_cond_a(x):
                     C_cond_a = ((8.84376730791660702116031235434001332739763212575923034819425083696842193603515625E-13)*x**10) + ((-7.38542591609163222308431502420578151912877729046158492565155029296875E-10)*x**9) + ((2.324153491478158013963553445158094490352596039883792400360107421875E-7)*x**8) + ((-0.0000279885896856371226763642023893652321930858306586742401123046875)*x**7) + ((-0.001211050864123853175458034314715405344031751155853271484375)*x**6) + ((0.58599055646336861702394571693730540573596954345703125)*x**5) + ((-9.63255533408176489729157765395939350128173828125)*x**4) + ((-11422.089305932286151801235973834991455078125)*x**3) + ((1594435.72884632088243961334228515625)*x**2) + ((-90438929.49676106870174407958984375)*x) + 1960950344.6000106334686279296875
                     return C_cond_a            
                 if T + 24.9 - (0.95851056*np.log(t_formation+0.008)+4.627989704197322) <= 131.5:
                    C_cond = 1
                 elif 131.5 < T + 24.9 - (0.95851056*np.log(t_formation+0.008)+4.627989704197322) < 144.2 :    
                    C_cond = C_cond_a(T + 24.9 - (0.95851056*np.log(t_formation+0.008)+4.627989704197322))
                 else:
                    C_cond = 0
                 return C_cond
        
@jit(nopython=True)
def C(d_formation, z_formation, t_formation):
            if z_formation <= 0:
               C_final = C_c(T_disc(d_formation, t_formation),t_formation)
            else:
                C_f = np.zeros(601)
                x = fznd(d_formation,z_formation)
                for j in frange (d_formation-3,d_formation+3.01,0.01):
                    C_f[int(round((j-(d_formation-3))/(0.01)))] = x[int(round((j-(d_formation-3))/(0.01)))]*C_c(T_disc(j, t_formation),t_formation)
                C_final = 0
                for i in range(0,len(C_f)):
                    C_final = C_final + C_f[i]
            return C_final    
        
    
@jit(nopython=True)
def Nz_c(T,t_formation):             
            def Nz_cond_a(x):
                Nz_cond_a = ((8.84376730791660702116031235434001332739763212575923034819425083696842193603515625E-13)*x**10) + ((-7.38542591609163222308431502420578151912877729046158492565155029296875E-10)*x**9) + ((2.324153491478158013963553445158094490352596039883792400360107421875E-7)*x**8) + ((-0.0000279885896856371226763642023893652321930858306586742401123046875)*x**7) + ((-0.001211050864123853175458034314715405344031751155853271484375)*x**6) + ((0.58599055646336861702394571693730540573596954345703125)*x**5) + ((-9.63255533408176489729157765395939350128173828125)*x**4) + ((-11422.089305932286151801235973834991455078125)*x**3) + ((1594435.72884632088243961334228515625)*x**2) + ((-90438929.49676106870174407958984375)*x) + 1960950344.6000106334686279296875
                return Nz_cond_a
            if T + 60.4 - ((-0.2)*t_formation**2)+((-1.1)*t_formation) <= 131.5:
                    N_cond = 1
            elif 131.5 < T + 60.4 - ((-0.2)*t_formation**2)+((-1.1)*t_formation) < 144.2 :    
                    N_cond = Nz_cond_a(T + 60.4 - ((-0.2)*t_formation**2)+((-1.1)*t_formation))
            else:
                    N_cond = 0
            return N_cond
        
@jit(nopython=True)    
def Nz(d_formation, z_formation, t_formation):        
            if z_formation <= 0:
               N_final = Nz_c(T_disc(d_formation, t_formation),t_formation)
            else:
                N_f = np.zeros(601)
                x = fznd(d_formation,z_formation)
                for j in frange (d_formation-3,d_formation+3.01,0.01):
                    N_f[int(round((j-(d_formation-3))/(0.01)))] = x[int(round((j-(d_formation-3))/(0.01)))]*Nz_c(T_disc(j, t_formation),t_formation)
                N_final = 0
                for i in range(0,len(N_f)):
                    N_final = N_final + N_f[i]
            return N_final     
        
@jit(nopython=True)        
def O(T,t_formation,d_formation, fe_star, z_formation, Al, Ti, Ca, Ni, Fe, Cr, Mg, Si, Na):
            def O_cond_a(x):
                O_cond_a = ((3.0130447356937058120666950537328292354655179252631569397635757923126220703125E-11)*x**9) + ((-2.807084698605554265497685913731407136850748429424129426479339599609375E-8)*x**8) + ((0.00001063495621534568083728367060558639423106797039508819580078125)*x**7) + ((-0.001953600072816120969410036423141718842089176177978515625)*x**6) + ((0.118782323977384363100640030097565613687038421630859375)*x**5) + ((20.507024807688964074259274639189243316650390625)*x**4) + ((-5057.821515620658828993327915668487548828125)*x**3) + ((473930.6469112207996658980846405029296875)*x**2) + ((-21956923.4987984411418437957763671875)*x) + 416374115.999850690364837646484375
                return O_cond_a
            def O_cond_b(x):
                O_cond_b = ((-8.01816773465877753609994077705340824082210189017581964342727685046696706194779835641384124755859375E-21)*x**12) + ((1.0070898330734577065393492579489358624709802108123692099272972200196818448603153228759765625E-17)*x**11) + ((-4.70477569404908565433640843044907841217160647528938710593138239346444606781005859375E-15)*x**10) + ((8.019035866412573751495925971206519557686609456226278780377469956874847412109375E-13)*x**9) + ((7.052309561630735782856521997878271844351072417111936374567449092864990234375E-11)*x**8) + ((-3.59520457904541848210907433770755492474791026324965059757232666015625E-8)*x**7) + ((-0.000001434755643641682662121623718920471191040633129887282848358154296875)*x**6) + ((0.0015651353277274580692857153252361968043260276317596435546875)*x**5) + ((-0.0157179314971404616996242253890159190632402896881103515625)*x**4) + ((-69.6575055560066260795792913995683193206787109375)*x**3) + ((13076.950970530795530066825449466705322265625)*x**2) + ((-1011028.2469208813272416591644287109375)*x) + 29949276.2144540511071681976318359375
                return O_cond_b
            def O_cond_c(x):
                O_cond_c = ((-1.935858687654345781991909840331329427273363148123906995386688320272453210581405091961215703122434206306934356689453125E-28)*x**14) + ((4.60031890387486682816812327205108474327971749643549747409456597969749924725846312867361120879650115966796875E-25)*x**13) + ((-4.04574720873363699965459756008350530694121082013269309076396174074119471697486005723476409912109375E-22)*x**12) + ((1.279783157740830553011536130703199131090575326909135702758979480364587288931943476200103759765625E-19)*x**11) + ((2.1961741772431609635338623975292609306017658151745884642647155260419822297990322113037109375E-17)*x**10) + ((-1.892855896132653454907952313138831504740170093292750408409119700081646442413330078125E-14)*x**9) + ((-2.40397966795433912086858782724643778840045715838869000435806810855865478515625E-12)*x**8) + ((2.65923217268813523813786637325744244275682603984023444354534149169921875E-9)*x**7) + ((3.52574456327930609681755268203229292112155235372483730316162109375E-7)*x**6) + ((-0.0003972032161178485503739976625325880377204157412052154541015625)*x**5) + ((-0.0236000284548833205722218053779215551912784576416015625)*x**4) + ((65.7878558825435533208292326889932155609130859375)*x**3) + ((-19576.4395998625623178668320178985595703125)*x**2) + ((2509944.56871254742145538330078125)*x) - 124785614.3718341290950775146484375
                return O_cond_c        
            if T - 1.17393656*np.log(t_formation+0.005)+6.219888463073661 <= 128:
                    O_cond = 1
            elif 128 < T - 1.17393656*np.log(t_formation+0.005)+6.219888463073661 <= 142.6 :    
                    O_cond = ((1.5*(Al/Mg)*stellarfloat[int(round(fe_star))][0])+(2*(Ti/Mg)*(stellarfloat[int(round(fe_star))][1]))+((Ca/Mg)*(stellarfloat[int(round(fe_star))][2]))+(1)+(2*(Si/Mg)*(stellarfloat[int(round(fe_star))][6]))+(1.5*(Cr/Mg)*(stellarfloat[int(round(fe_star))][5]))+((1.38)*(Fe/Mg)*(stellarfloat[int(round(fe_star))][4]))+(0.5*(Na/Mg)*(stellarfloat[int(round(fe_star))][7])))/((1/Mg)*(stellarfloat[int(round(fe_star))][8])) + (O_cond_a(T - 1.17393656*np.log(t_formation+0.005)+6.219888463073661)*(1-((1.5*(Al/Mg)*stellarfloat[int(round(fe_star))][0])+(2*(Ti/Mg)*(stellarfloat[int(round(fe_star))][1]))+((Ca/Mg)*(stellarfloat[int(round(fe_star))][2]))+(1)+(2*(Si/Mg)*(stellarfloat[int(round(fe_star))][6]))+(1.5*(Cr/Mg)*(stellarfloat[int(round(fe_star))][5]))+((1.38)*(Fe/Mg)*(stellarfloat[int(round(fe_star))][4]))+(0.5*(Na/Mg)*(stellarfloat[int(round(fe_star))][7])))/((1/Mg)*(stellarfloat[int(round(fe_star))][8]))))
            elif 142.6 + 1.17393656*np.log(t_formation+0.005)-6.219888463073661 < T <= 180 :
                    O_cond =  ((1.5*(Al/Mg)*stellarfloat[int(round(fe_star))][0])+(2*(Ti/Mg)*(stellarfloat[int(round(fe_star))][1]))+((Ca/Mg)*(stellarfloat[int(round(fe_star))][2]))+(1)+(2*(Si/Mg)*(stellarfloat[int(round(fe_star))][6]))+(1.5*(Cr/Mg)*(stellarfloat[int(round(fe_star))][5]))+((1.38)*(Fe/Mg)*(stellarfloat[int(round(fe_star))][4]))+(0.5*(Na/Mg)*(stellarfloat[int(round(fe_star))][7])))/((1/Mg)*(stellarfloat[int(round(fe_star))][8])) 
            elif 180 < T <= 220 : 
                    O_cond =  ((1.5*(Al/Mg)*stellarfloat[int(round(fe_star))][0])+(2*(Ti/Mg)*(stellarfloat[int(round(fe_star))][1]))+((Ca/Mg)*(stellarfloat[int(round(fe_star))][2]))+(1)+(2*(Si/Mg)*(stellarfloat[int(round(fe_star))][6]))+(1.5*(Cr/Mg)*(stellarfloat[int(round(fe_star))][5]))+((O_cond_b(T))*(Fe/Mg)*(stellarfloat[int(round(fe_star))][4]))+(0.5*(Na/Mg)*(stellarfloat[int(round(fe_star))][7])))/((1/Mg)*(stellarfloat[int(round(fe_star))][8]))       
            elif 220 < T <= 325 : 
                    O_cond =  ((1.5*(Al/Mg)*stellarfloat[int(round(fe_star))][0])+(2*(Ti/Mg)*(stellarfloat[int(round(fe_star))][1]))+((Ca/Mg)*(stellarfloat[int(round(fe_star))][2]))+(1)+(2*(Si/Mg)*(stellarfloat[int(round(fe_star))][6]))+(1.5*(Cr/Mg)*(stellarfloat[int(round(fe_star))][5]))+((0.651)*(Fe/Mg)*(stellarfloat[int(round(fe_star))][4]))+(0.5*(Na/Mg)*(stellarfloat[int(round(fe_star))][7])))/((1/Mg)*(stellarfloat[int(round(fe_star))][8]))
            elif 325 < T <= 409 :    
                    O_cond = ((1.5*(Al/Mg)*stellarfloat[int(round(fe_star))][0])+(2*(Ti/Mg)*(stellarfloat[int(round(fe_star))][1]))+((Ca/Mg)*(stellarfloat[int(round(fe_star))][2]))+(1)+(2*(Si/Mg)*(stellarfloat[int(round(fe_star))][6]))+(1.5*(Cr/Mg)*(stellarfloat[int(round(fe_star))][5]))+((O_cond_c(T))*(Fe/Mg)*(stellarfloat[int(round(fe_star))][4]))+(0.5*(Na/Mg)*(stellarfloat[int(round(fe_star))][7])))/((1/Mg)*(stellarfloat[int(round(fe_star))][8]))
            else:
                    O_cond = ((1.5*(Al/Mg)*stellarfloat[int(round(fe_star))][0])+(2*(Ti/Mg)*(stellarfloat[int(round(fe_star))][1]))+((Ca/Mg)*(stellarfloat[int(round(fe_star))][2]))+(1)+(2*(Si/Mg)*(stellarfloat[int(round(fe_star))][6]))+(1.5*(Cr/Mg)*(stellarfloat[int(round(fe_star))][5]))+(0.5*(Na/Mg)*(stellarfloat[int(round(fe_star))][7])))/((1/Mg)*(stellarfloat[int(round(fe_star))][8]))
            return O_cond    
      
    
def E_Al(N_c,N_o,f_c,f_o):
            if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
               f_m = 1 - f_c - f_o
               N_m = 1 - N_c - N_o
               x_m = (0.0153-(0.0641)*(N_o)-(0)*(N_c))
               if x_m >= 0:
                   x_m = x_m
                   x_o = 0.0641
                   x_c = 0
               else:
                   x_m = 0
                   x_o = 0.0153/N_o
                   x_c = 0
               E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
               E_2 = (0.001*(0.0641)+0.829*((0.0153-(0.0641)*(0.001)-(0)*(0.17))/(0.829))+0.17*(0))
               E_final = (E_1)/(E_2)
            else:
                E_final = 0
            return E_final
    
def E_Ti(N_c,N_o,f_c,f_o):
            if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
               f_m = 1 - f_c - f_o
               N_m = 1 - N_c - N_o
               x_m = (0.00044-(0.0042)*(N_o)-(0)*(N_c))
               if x_m >= 0:
                   x_m = x_m
                   x_o = 0.0042
                   x_c = 0
               else:
                   x_m = 0
                   x_o = 0.00044/N_o
                   x_c = 0
               E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
               E_2 = (0.001*(0.0042)+0.829*((0.00044-(0.0042)*(0.001)-(0)*(0.17))/(0.829))+0.17*(0))
               E_final = (E_1)/(E_2)
            else:
                E_final = 0
            return E_final
    
def E_Ca(N_c,N_o,f_c,f_o):
            if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
               f_m = 1 - f_c - f_o
               N_m = 1 - N_c - N_o
               x_m = (0.011099-(0.04452)*(N_o)-(0)*(N_c))
               if x_m >= 0:
                   x_m = x_m
                   x_o = 0.04452
                   x_c = 0
               else:
                   x_m = 0
                   x_o = 0.011099/N_o
                   x_c = 0
               E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
               E_2 = (0.001*(0.04452)+0.829*((0.011099-(0.04452)*(0.001)-(0)*(0.17))/(0.829))+0.17*(0))
               E_final = (E_1)/(E_2)
            else:
                E_final = 0
            return E_final
    
def E_Na(N_c,N_o,f_c,f_o):
            if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >=N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
               f_m = 1 - f_c - f_o
               N_m = 1 - N_c - N_o
               x_m = (0.002037-(0.01771)*(N_o)-(0)*(N_c))
               if x_m >= 0:
                   x_m = x_m
                   x_o = 0.01771
                   x_c = 0
               else:
                   x_m = 0
                   x_o = 0.002037/N_o
                   x_c = 0
               E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
               E_2 = (0.001*(0.01771)+0.829*((0.002037-(0.01771)*(0.001)-(0)*(0.17))/(0.829))+0.17*(0))
               E_final = (E_1)/(E_2)
            else:
                E_final = 0
            return E_final
    
def E_O(N_c,N_o,f_c,f_o):
            if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
               f_m = 1 - f_c - f_o
               N_m = 1 - N_c - N_o
               x_m = (0.482879-(0.6011)*(N_o)-(0)*(N_c))
               if x_m >= 0:
                   x_m = x_m
                   x_o = 0.6011
                   x_c = 0
               else:
                   x_m = 0
                   x_o = 0.482879/N_o
                   x_c = 0
               E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
               E_2 = (0.001*(0.6011)+0.829*((0.482879-(0.6011)*(0.001)-(0)*(0.17))/(0.829))+0.17*(0))
               E_final = (E_1)/(E_2)
            else:
                E_final = 0
            return E_final
    
    
def E_Mg(N_c,N_o,f_c,f_o):
            if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
               f_m = 1 - f_c - f_o
               N_m = 1 - N_c - N_o
               x_m = (0.16482-(0.04167)*(N_o)-(0)*(N_c))
               if x_m >= 0:
                   x_m = x_m
                   x_o = 0.04167
                   x_c = 0
               else:
                   x_m = 0
                   x_o = 0.16482/N_o
                   x_c = 0
               E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
               E_2 = (0.001*(0.04167)+0.829*((0.16482-(0.04167)*(0.001)-(0)*(0.17))/(0.829))+0.17*(0))
               E_final = (E_1)/(E_2)
            else:
                E_final = 0.0001
            return E_final
    
def E_Fe(N_c,N_o,f_c,f_o):
            if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
               f_m = 1 - f_c - f_o
               N_m = 1 - N_c - N_o
               x_m = (0.149057-(0.0314)*(N_o)-(0.7676)*(N_c))
               if x_m >= 0:
                   x_m = x_m
                   x_o = 0.0314
                   x_c = 0.7676
               else:
                   x_m = 0
                   x_o = (0.149057-0.7676*(N_c))/N_o
                   x_c = 0.7676
               E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
               E_2 = (0.001*(0.0314)+0.829*((0.149057-(0.0314)*(0.001)-(0.7676)*(0.17))/(0.829))+0.17*(0.7676))
               E_final = (E_1)/(E_2)
            else:
                E_final = 0
            return E_final
    
def E_Ni(N_c,N_o,f_c,f_o):
            if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
               f_m = 1 - f_c - f_o
               N_m = 1 - N_c - N_o
               x_m = (0.008066-(0.0000371)*(N_o)-(0.0444)*(N_c))
               if x_m >= 0:
                   x_m = x_m
                   x_o = 0.0000371
                   x_c = 0.0444
               else:
                   x_m = 0
                   x_o = 0.0000371
                   x_c = (0.008066-(0.0000371)*N_o)/(N_c)
               E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
               E_2 = (0.001*(0.0000371)+0.829*((0.008066-(0.0000371)*(0.001)-(0.0444)*(0.17))/(0.829))+0.17*(0.0444))
               E_final = (E_1)/(E_2)
            else:
                E_final = 0
            return E_final
    
def E_Cr(N_c,N_o,f_c,f_o):
            if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
               f_m = 1 - f_c - f_o
               N_m = 1 - N_c - N_o
               x_m = (0.002351-(0.000139)*(N_o)-(0.00868)*(N_c))
               if x_m >= 0:
                   x_m = x_m
                   x_o = 0.000139
                   x_c = 0.00868
               else:
                   x_m = 0
                   x_o = (0.002351-0.00868*(N_c))/(N_o)
                   x_c = 0.00868
               E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
               E_2 = (0.001*(0.000139)+0.829*((0.002351-(0.000139)*(0.001)-(0.00868)*(0.17))/(0.829))+0.17*(0.00868))
               E_final = (E_1)/(E_2)
            else:
                E_final = 0
            return E_final
    
def E_Si(N_c,N_o,f_c,f_o):
            if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
               f_m = 1 - f_c - f_o
               N_m = 1 - N_c - N_o
               x_m = (0.149118-(0.181509)*(N_o)-(0.1071)*(N_c))
               if x_m >= 0:
                   x_m = x_m
                   x_o = 0.181509
                   x_c = 0.1071
               else:
                   x_m = 0
                   x_o = (0.149118-0.1071*(N_c))/(N_o)
                   x_c = 0.1071
               E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
               E_2 = (0.001*(0.181509)+0.829*((0.149118-(0.181509)*(0.001)-(0.1071)*(0.17))/(0.829))+0.17*(0.1071))
               E_final = (E_1)/(E_2)
            else:
                E_final = 0
            return E_final
    
def E_C(N_c,N_o,f_c,f_o):
            if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
               f_m = 1 - f_c - f_o
               N_m = 1 - N_c - N_o
               x_m = (0.001581-(0)*(N_o)-(0.0083482)*(N_c))
               if x_m >= 0:
                   x_m = x_m
                   x_o = 0
                   x_c = 0.0083482
               else:
                   x_m = 0
                   x_o = 0
                   x_c = 0.001581/N_c
               E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
               E_2 = (0.001*(0)+0.829*((0.001581-(0)*(0.001)-(0.0083482)*(0.17))/(0.829))+0.17*(0.0083482))
               E_final = (E_1)/(E_2)
            else:
                E_final = 0
            return E_final
    
def E_Nz(N_c,N_o,f_c,f_o):
            if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
               f_m = 1 - f_c - f_o
               N_m = 1 - N_c - N_o
               x_m = (0.000046429-(0)*(N_o)-(0.0002684)*(N_c))
               if x_m >= 0:
                   x_m = x_m
                   x_o = 0
                   x_c = 0.0002684
               else:
                   x_m = 0
                   x_o = 0
                   x_c = 0.000046429/N_c
               E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
               E_2 = (0.001*(0)+0.829*((0.000046429-(0)*(0.001)-(0.0002684)*(0.17))/(0.829))+0.17*(0.0002684))
               E_final = (E_1)/(E_2)
            else:
                E_final = 0
            return E_final
      
#Defining Analysis Functions     
def confidence_intervals(sample_draws, array, length):
            
        prob = np.empty(sample_draws)
        prob.fill(1.0/sample_draws)   # Each sample assumed to be equally likely
            
        sig_1 = 0.5 + 0.6826/2.0  # Cummulative propabiltiy of +1 sigma level
        sig_2 = 0.5 + 0.954/2.0   # Cummulative propabiltiy of +2 sigma level
        sig_3 = 0.5 + 0.997/2.0   # Cummulative propabiltiy of +3 sigma level
            
        if (length > 0):  # If array is a vector (e.g. a line / curve where we want contours for each x value)
                
            arr_low3 = np.zeros(shape=(length))
            arr_low2 = np.zeros(shape=(length))
            arr_low1 = np.zeros(shape=(length))
            arr_median = np.zeros(shape=(length))
            arr_high1 = np.zeros(shape=(length))
            arr_high2 = np.zeros(shape=(length))
            arr_high3 = np.zeros(shape=(length))
                
            for i in range(length):
            
                arr_ordered = list(zip(prob[:], array[:, i]))
                arr_ordered.sort(key=lambda x: x[1])
                arr_ordered = np.array(arr_ordered)
        
                arr_ordered[:,0] = arr_ordered[:,0].cumsum()
        
                arr_ordered_interp = lambda x: np.interp(x, arr_ordered[:,0], arr_ordered[:,1],
                                                         left=arr_ordered[0,1], right=arr_ordered[-1,1])
        
                arr_low3[i] = arr_ordered_interp(1-sig_3)
                arr_low2[i] = arr_ordered_interp(1-sig_2)
                arr_low1[i] = arr_ordered_interp(1-sig_1)
                arr_median[i] = arr_ordered_interp(0.5)
                arr_high1[i] = arr_ordered_interp(sig_1)
                arr_high2[i] = arr_ordered_interp(sig_2) 
                arr_high3[i] = arr_ordered_interp(sig_3) 
                
            return arr_low3, arr_low2, arr_low1, arr_median, arr_high1, arr_high2, arr_high3
                
        if (length == 0):  # If quantity is just a float
                
            arr_ordered = list(zip(prob[:], array[:]))
            arr_ordered.sort(key=lambda x: x[1])
            arr_ordered = np.array(arr_ordered)
        
            arr_ordered[:,0] = arr_ordered[:,0].cumsum()
        
            arr_ordered_interp = lambda x: np.interp(x, arr_ordered[:,0], arr_ordered[:,1],
                                                     left=arr_ordered[0,1], right=arr_ordered[-1,1])
    
            arr_low3 = arr_ordered_interp(1-sig_3)
            arr_low2 = arr_ordered_interp(1-sig_2)
            arr_low1 = arr_ordered_interp(1-sig_1)
            arr_median = arr_ordered_interp(0.5)
            arr_high1 = arr_ordered_interp(sig_1)
            arr_high2 = arr_ordered_interp(sig_2) 
            arr_high3 = arr_ordered_interp(sig_3) 
                    
            return arr_low3, arr_low2, arr_low1, arr_median, arr_high1, arr_high2, arr_high3   
        
        
def Z_to_sigma(ln_Z1, ln_Z2):  
        np.set_printoptions(precision=50)
        B = np.exp(ln_Z1 - ln_Z2) 
        p = np.real(np.exp(W((-1.0/(B*np.exp(1))),-1)))
        sigma = np.sqrt(2)*erfcinv(p)
        print("p-value = ", p)
        print("n_sigma = ", sigma)
        return B, sigma  
    
def Detection_significance(additional_parameter, err_data, N_data, n_params, n_removed, basename_M1, basename_M2):
    
        norm_log = (-0.5*np.log(2.0*np.pi*err_data*err_data)).sum()
        
        retrieved_M1 = pymultinest.Analyzer(n_params = n_params, outputfiles_basename = basename_M1)
        retrieved_M2 = pymultinest.Analyzer(n_params = (n_params-n_removed), outputfiles_basename = basename_M2)
        s_M1 = retrieved_M1.get_stats()
        s_M2 = retrieved_M2.get_stats()
    
        best_fit_M1 = retrieved_M1.get_best_fit()                       # Best-fitting parameter combination from model 1
        max_likelihood_M1 = best_fit_M1['log_likelihood']               # Maximum log-likelihood for model 1
        best_chi_square_M1 = -2.0 * (max_likelihood_M1 - norm_log)      # Best chi-square for model 1

    
        ln_Z_M1 = s_M1['global evidence']     # Natural log of Bayesian evidence from model 1 
        
        print("Evidence with " + additional_parameter + " = " + str(ln_Z_M1))
                    
        best_fit_M2 = retrieved_M2.get_best_fit()                                      # Best-fitting parameter combination from model 2
        max_likelihood_M2 = best_fit_M2['log_likelihood']                              # Maximum log-likelihood for model 2
        best_chi_square_M2 = -2.0 * (max_likelihood_M2 - norm_log)                     # Best chi-square for model 2

        
        ln_Z_M2 = s_M2['global evidence']     # Natural log of Bayesian evidence from model 1 
        
        print("Evidence without " + additional_parameter + " = " + str(ln_Z_M2))
                    
        Bayes_factor, n_sigma = Z_to_sigma(ln_Z_M1, ln_Z_M2)
        
        print("Bayes_factor of " + additional_parameter + " = " + str(Bayes_factor))
        print("Sigma significance of " + additional_parameter + " = " + str(n_sigma))
        
        return ln_Z_M1, ln_Z_M2, Bayes_factor, n_sigma, best_chi_square_M1, best_chi_square_M2   



parameters =["Stellar metallicity indices","time","Distance","core","log(Pollution Fraction)","event"]
n_params = len(parameters)
prefix = "chains/"+str(N)+"model27"

#RUNNING MULTINEST SAMPLING
#result = solve(LogLikelihood=loglike, Prior=prior, n_dims=n_params, outputfiles_basename=prefix, verbose=True, n_live_points=1000)
#print('evidence: %(logZ).7f +- %(logZerr).7f' % result)
#print("MultiNest Model Completed") 
#with open('%sparams.json' % prefix, 'w') as f:
#    json.dump(parameters, f, indent=2)

#CREATE CORNER PLOT
#a = pymultinest.Analyzer(n_params = n_params, outputfiles_basename = prefix)
#stats = a.get_stats()
#data = a.get_equal_weighted_posterior()[:,0:6]
#corner.corner(data, labels=parameters, label_kwargs=dict(fontsize=12), levels = (0.39346934,0.86466472,0.988891), smooth=True)



#CREATE COMPOSITIONAL PLOT
plt.figure(figsize=(10,5))
x_axis = ['Al/Mg','Ti/Mg','Ca/Mg','Ni/Mg','Fe/Mg','Cr/Mg','Si/Mg','Na/Mg','O/Mg','C/Mg','N/Mg']
values = pymultinest.Analyzer(n_params = n_params, outputfiles_basename = prefix).get_equal_weighted_posterior()
N_values = len(values[:,0])
sample_draws = min(10000, N_values)
sample = np.random.choice(len(values), sample_draws, replace=False)
YOUR_MODEL_OUTPUT_stored = np.zeros(shape=(sample_draws, len(x_axis))) 
for i in range(sample_draws):
    fe_star = values[sample[i],0]
    t_sinceaccretion = values[sample[i],1]
    d_formation = values[sample[i],2]
    t_formation = 1.5
    z_formation = 0.05
    N_c = 0.17
    N_o = 0.01
    f_c = values[sample[i],3]
    f_o = 0
    pollutionfraction = values[sample[i],4]
    t_disc = 10**(values[sample[i],5])
    Al_hsc = Al(10**(d_formation),z_formation,t_formation)
    Ti_hsc = Ti(10**(d_formation),z_formation,t_formation)
    Ca_hsc = Ca(10**(d_formation),z_formation,t_formation)
    Ni_hsc = Ni(10**(d_formation),z_formation,t_formation)
    Fe_hsc = Fe(10**(d_formation),z_formation,t_formation)
    Cr_hsc = Cr(10**(d_formation),z_formation,t_formation)
    Mg_hsc = Mg(10**(d_formation),z_formation,t_formation)
    Si_hsc = Si(10**(d_formation),z_formation,t_formation)
    Na_hsc = Na(10**(d_formation),z_formation,t_formation)
    O_hsc = O(T_disc(10**(d_formation),t_formation),t_formation,10**(d_formation),fe_star,z_formation,Al_hsc,Ti_hsc,Ca_hsc,Ni_hsc,Fe_hsc,Cr_hsc,Mg_hsc,Si_hsc,Na_hsc)
    C_hsc = C(10**(d_formation),z_formation,t_formation)  
    Nz_hsc = Nz(10**(d_formation),z_formation,t_formation)
    Al_dm = E_Al(N_c,N_o,f_c,f_o)
    Ti_dm = E_Ti(N_c,N_o,f_c,f_o)
    Ca_dm = E_Ca(N_c,N_o,f_c,f_o)
    Ni_dm = E_Ni(N_c,N_o,f_c,f_o)
    Fe_dm = E_Fe(N_c,N_o,f_c,f_o)
    Cr_dm = E_Cr(N_c,N_o,f_c,f_o)
    Mg_dm = E_Mg(N_c,N_o,f_c,f_o)
    Si_dm = E_Si(N_c,N_o,f_c,f_o)
    Na_dm = E_Na(N_c,N_o,f_c,f_o) 
    O_dm = E_O(N_c,N_o,f_c,f_o)
    C_dm = E_C(N_c,N_o,f_c,f_o)
    Nz_dm = E_Nz(N_c,N_o,f_c,f_o)
                                        
    if t_sinceaccretion*1000000 == 0:
                                         if np.log10(stellarfloat[int(round(fe_star))][9]*(C_dm/Mg_dm)*(C_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][10])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6])))- 0.8507 > -2 and np.log10(stellarfloat[int(round(fe_star))][10]*(Nz_dm/Mg_dm)*(Nz_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][11])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 0.3113 > -2:
                                           YOUR_MODEL_OUTPUT_stored[i,:] = [np.log10(stellarfloat[int(round(fe_star))][0]*(Al_dm/Mg_dm)*(Al_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][0])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.145, np.log10(stellarfloat[int(round(fe_star))][1]*(Ti_dm/Mg_dm)*(Ti_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][1])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 2.587, np.log10(stellarfloat[int(round(fe_star))][2]*(Ca_dm/Mg_dm)*(Ca_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][2])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.169, np.log10(stellarfloat[int(round(fe_star))][3]*(Ni_dm/Mg_dm)*(Ni_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][3])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.289, np.log10(stellarfloat[int(round(fe_star))][4]*(Fe_dm/Mg_dm)*(Fe_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][4])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.047, np.log10(stellarfloat[int(round(fe_star))][5]*(Cr_dm/Mg_dm)*(Cr_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][5])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.868, np.log10(stellarfloat[int(round(fe_star))][6]*(Si_dm/Mg_dm)*(Si_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][7])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.011, np.log10(stellarfloat[int(round(fe_star))][7]*(Na_dm/Mg_dm)*(Na_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][8])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.346, np.log10(stellarfloat[int(round(fe_star))][8]*(O_dm/Mg_dm)*(O_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][9])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 1.2284, np.log10(stellarfloat[int(round(fe_star))][9]*(C_dm/Mg_dm)*(C_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][10])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6])))- 0.8507, np.log10(stellarfloat[int(round(fe_star))][10]*(Nz_dm/Mg_dm)*(Nz_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][11])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 0.3113]    
                                         elif np.log10(stellarfloat[int(round(fe_star))][9]*(C_dm/Mg_dm)*(C_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][10])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6])))- 0.8507 > -2 and np.log10(stellarfloat[int(round(fe_star))][10]*(Nz_dm/Mg_dm)*(Nz_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][11])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 0.3113 <= -2:
                                           YOUR_MODEL_OUTPUT_stored[i,:] = [np.log10(stellarfloat[int(round(fe_star))][0]*(Al_dm/Mg_dm)*(Al_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][0])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.145, np.log10(stellarfloat[int(round(fe_star))][1]*(Ti_dm/Mg_dm)*(Ti_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][1])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 2.587, np.log10(stellarfloat[int(round(fe_star))][2]*(Ca_dm/Mg_dm)*(Ca_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][2])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.169, np.log10(stellarfloat[int(round(fe_star))][3]*(Ni_dm/Mg_dm)*(Ni_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][3])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.289, np.log10(stellarfloat[int(round(fe_star))][4]*(Fe_dm/Mg_dm)*(Fe_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][4])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.047, np.log10(stellarfloat[int(round(fe_star))][5]*(Cr_dm/Mg_dm)*(Cr_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][5])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.868, np.log10(stellarfloat[int(round(fe_star))][6]*(Si_dm/Mg_dm)*(Si_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][7])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.011, np.log10(stellarfloat[int(round(fe_star))][7]*(Na_dm/Mg_dm)*(Na_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][8])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.346, np.log10(stellarfloat[int(round(fe_star))][8]*(O_dm/Mg_dm)*(O_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][9])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 1.2284, np.log10(stellarfloat[int(round(fe_star))][9]*(C_dm/Mg_dm)*(C_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][10])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6])))- 0.8507, -2] 
                                         elif np.log10(stellarfloat[int(round(fe_star))][9]*(C_dm/Mg_dm)*(C_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][10])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6])))- 0.8507 <= -2 and np.log10(stellarfloat[int(round(fe_star))][10]*(Nz_dm/Mg_dm)*(Nz_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][11])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 0.3113 > -2:
                                           YOUR_MODEL_OUTPUT_stored[i,:] = [np.log10(stellarfloat[int(round(fe_star))][0]*(Al_dm/Mg_dm)*(Al_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][0])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.145, np.log10(stellarfloat[int(round(fe_star))][1]*(Ti_dm/Mg_dm)*(Ti_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][1])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 2.587, np.log10(stellarfloat[int(round(fe_star))][2]*(Ca_dm/Mg_dm)*(Ca_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][2])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.169, np.log10(stellarfloat[int(round(fe_star))][3]*(Ni_dm/Mg_dm)*(Ni_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][3])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.289, np.log10(stellarfloat[int(round(fe_star))][4]*(Fe_dm/Mg_dm)*(Fe_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][4])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.047, np.log10(stellarfloat[int(round(fe_star))][5]*(Cr_dm/Mg_dm)*(Cr_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][5])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.868, np.log10(stellarfloat[int(round(fe_star))][6]*(Si_dm/Mg_dm)*(Si_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][7])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.011, np.log10(stellarfloat[int(round(fe_star))][7]*(Na_dm/Mg_dm)*(Na_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][8])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.346, np.log10(stellarfloat[int(round(fe_star))][8]*(O_dm/Mg_dm)*(O_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][9])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 1.2284, -2, np.log10(stellarfloat[int(round(fe_star))][10]*(Nz_dm/Mg_dm)*(Nz_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][11])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 0.3113]    
                                         else:
                                           YOUR_MODEL_OUTPUT_stored[i,:] = [np.log10(stellarfloat[int(round(fe_star))][0]*(Al_dm/Mg_dm)*(Al_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][0])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.145, np.log10(stellarfloat[int(round(fe_star))][1]*(Ti_dm/Mg_dm)*(Ti_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][1])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 2.587, np.log10(stellarfloat[int(round(fe_star))][2]*(Ca_dm/Mg_dm)*(Ca_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][2])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.169, np.log10(stellarfloat[int(round(fe_star))][3]*(Ni_dm/Mg_dm)*(Ni_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][3])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.289, np.log10(stellarfloat[int(round(fe_star))][4]*(Fe_dm/Mg_dm)*(Fe_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][4])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.047, np.log10(stellarfloat[int(round(fe_star))][5]*(Cr_dm/Mg_dm)*(Cr_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][5])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.868, np.log10(stellarfloat[int(round(fe_star))][6]*(Si_dm/Mg_dm)*(Si_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][7])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.011, np.log10(stellarfloat[int(round(fe_star))][7]*(Na_dm/Mg_dm)*(Na_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][8])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.346, np.log10(stellarfloat[int(round(fe_star))][8]*(O_dm/Mg_dm)*(O_hsc/Mg_hsc)*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][9])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 1.2284, -2, -2] 
    elif 0 < t_sinceaccretion*1000000 < t_disc:
                                         if np.log10(stellarfloat[int(round(fe_star))][9]*(C_dm/Mg_dm)*(C_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][10])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][10])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][10])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6])))- 0.8507 > -2 and np.log10(stellarfloat[int(round(fe_star))][10]*(Nz_dm/Mg_dm)*(Nz_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][11])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][11])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][11])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 0.3113 > -2:
                                           YOUR_MODEL_OUTPUT_stored[i,:] = [np.log10(stellarfloat[int(round(fe_star))][0]*(Al_dm/Mg_dm)*(Al_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][0])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][0])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][0])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.145, np.log10(stellarfloat[int(round(fe_star))][1]*(Ti_dm/Mg_dm)*(Ti_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][1])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][1])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][1])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 2.587, np.log10(stellarfloat[int(round(fe_star))][2]*(Ca_dm/Mg_dm)*(Ca_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][2])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][2])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][2])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.169, np.log10(stellarfloat[int(round(fe_star))][3]*(Ni_dm/Mg_dm)*(Ni_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][3])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][3])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][3])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.289, np.log10(stellarfloat[int(round(fe_star))][4]*(Fe_dm/Mg_dm)*(Fe_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][4])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][4])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][4])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.047, np.log10(stellarfloat[int(round(fe_star))][5]*(Cr_dm/Mg_dm)*(Cr_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][5])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][5])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][5])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.868, np.log10(stellarfloat[int(round(fe_star))][6]*(Si_dm/Mg_dm)*(Si_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][7])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][7])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][7])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.011, np.log10(stellarfloat[int(round(fe_star))][7]*(Na_dm/Mg_dm)*(Na_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][8])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][8])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][8])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.346, np.log10(stellarfloat[int(round(fe_star))][8]*(O_dm/Mg_dm)*(O_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][9])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][9])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][9])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 1.2284, np.log10(stellarfloat[int(round(fe_star))][9]*(C_dm/Mg_dm)*(C_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][10])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][10])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][10])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6])))- 0.8507, np.log10(stellarfloat[int(round(fe_star))][10]*(Nz_dm/Mg_dm)*(Nz_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][11])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][11])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][11])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 0.3113]    
                                         elif np.log10(stellarfloat[int(round(fe_star))][9]*(C_dm/Mg_dm)*(C_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][10])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][10])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][10])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6])))- 0.8507 > -2 and np.log10(stellarfloat[int(round(fe_star))][10]*(Nz_dm/Mg_dm)*(Nz_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][11])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][11])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][11])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 0.3113 <= -2:
                                           YOUR_MODEL_OUTPUT_stored[i,:] = [np.log10(stellarfloat[int(round(fe_star))][0]*(Al_dm/Mg_dm)*(Al_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][0])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][0])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][0])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.145, np.log10(stellarfloat[int(round(fe_star))][1]*(Ti_dm/Mg_dm)*(Ti_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][1])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][1])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][1])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 2.587, np.log10(stellarfloat[int(round(fe_star))][2]*(Ca_dm/Mg_dm)*(Ca_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][2])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][2])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][2])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.169, np.log10(stellarfloat[int(round(fe_star))][3]*(Ni_dm/Mg_dm)*(Ni_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][3])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][3])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][3])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.289, np.log10(stellarfloat[int(round(fe_star))][4]*(Fe_dm/Mg_dm)*(Fe_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][4])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][4])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][4])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.047, np.log10(stellarfloat[int(round(fe_star))][5]*(Cr_dm/Mg_dm)*(Cr_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][5])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][5])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][5])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.868, np.log10(stellarfloat[int(round(fe_star))][6]*(Si_dm/Mg_dm)*(Si_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][7])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][7])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][7])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.011, np.log10(stellarfloat[int(round(fe_star))][7]*(Na_dm/Mg_dm)*(Na_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][8])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][8])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][8])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.346, np.log10(stellarfloat[int(round(fe_star))][8]*(O_dm/Mg_dm)*(O_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][9])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][9])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][9])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 1.2284, np.log10(stellarfloat[int(round(fe_star))][9]*(C_dm/Mg_dm)*(C_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][10])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][10])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][10])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6])))- 0.8507, -2] 
                                         elif np.log10(stellarfloat[int(round(fe_star))][9]*(C_dm/Mg_dm)*(C_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][10])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][10])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][10])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6])))- 0.8507 <= -2 and np.log10(stellarfloat[int(round(fe_star))][10]*(Nz_dm/Mg_dm)*(Nz_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][11])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][11])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][11])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 0.3113 > -2:
                                           YOUR_MODEL_OUTPUT_stored[i,:] = [np.log10(stellarfloat[int(round(fe_star))][0]*(Al_dm/Mg_dm)*(Al_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][0])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][0])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][0])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.145, np.log10(stellarfloat[int(round(fe_star))][1]*(Ti_dm/Mg_dm)*(Ti_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][1])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][1])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][1])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 2.587, np.log10(stellarfloat[int(round(fe_star))][2]*(Ca_dm/Mg_dm)*(Ca_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][2])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][2])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][2])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.169, np.log10(stellarfloat[int(round(fe_star))][3]*(Ni_dm/Mg_dm)*(Ni_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][3])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][3])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][3])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.289, np.log10(stellarfloat[int(round(fe_star))][4]*(Fe_dm/Mg_dm)*(Fe_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][4])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][4])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][4])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.047, np.log10(stellarfloat[int(round(fe_star))][5]*(Cr_dm/Mg_dm)*(Cr_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][5])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][5])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][5])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.868, np.log10(stellarfloat[int(round(fe_star))][6]*(Si_dm/Mg_dm)*(Si_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][7])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][7])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][7])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.011, np.log10(stellarfloat[int(round(fe_star))][7]*(Na_dm/Mg_dm)*(Na_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][8])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][8])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][8])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.346, np.log10(stellarfloat[int(round(fe_star))][8]*(O_dm/Mg_dm)*(O_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][9])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][9])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][9])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 1.2284, -2, np.log10(stellarfloat[int(round(fe_star))][10]*(Nz_dm/Mg_dm)*(Nz_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][11])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][11])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][11])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 0.3113]    
                                         else:
                                           YOUR_MODEL_OUTPUT_stored[i,:] = [np.log10(stellarfloat[int(round(fe_star))][0]*(Al_dm/Mg_dm)*(Al_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][0])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][0])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][0])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.145, np.log10(stellarfloat[int(round(fe_star))][1]*(Ti_dm/Mg_dm)*(Ti_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][1])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][1])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][1])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 2.587, np.log10(stellarfloat[int(round(fe_star))][2]*(Ca_dm/Mg_dm)*(Ca_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][2])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][2])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][2])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.169, np.log10(stellarfloat[int(round(fe_star))][3]*(Ni_dm/Mg_dm)*(Ni_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][3])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][3])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][3])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.289, np.log10(stellarfloat[int(round(fe_star))][4]*(Fe_dm/Mg_dm)*(Fe_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][4])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][4])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][4])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.047, np.log10(stellarfloat[int(round(fe_star))][5]*(Cr_dm/Mg_dm)*(Cr_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][5])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][5])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][5])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.868, np.log10(stellarfloat[int(round(fe_star))][6]*(Si_dm/Mg_dm)*(Si_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][7])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][7])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][7])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.011, np.log10(stellarfloat[int(round(fe_star))][7]*(Na_dm/Mg_dm)*(Na_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][8])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][8])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][8])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.346, np.log10(stellarfloat[int(round(fe_star))][8]*(O_dm/Mg_dm)*(O_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][9])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][9])))/(1-np.exp(-(1000000*t_sinceaccretion)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-((0*1000000)/vars()['whitedwarftimescales'+str(N)][9])+((0*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 1.2284, -2, -2] 
    else:
                                         if np.log10(stellarfloat[int(round(fe_star))][9]*(C_dm/Mg_dm)*(C_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][10])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][10])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][10])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6])))- 0.8507 > -2 and np.log10(stellarfloat[int(round(fe_star))][10]*(Nz_dm/Mg_dm)*(Nz_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][11])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][11])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][11])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 0.3113 > -2:
                                           YOUR_MODEL_OUTPUT_stored[i,:] = [np.log10(stellarfloat[int(round(fe_star))][0]*(Al_dm/Mg_dm)*(Al_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][0])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][0])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][0])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.145, np.log10(stellarfloat[int(round(fe_star))][1]*(Ti_dm/Mg_dm)*(Ti_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][1])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][1])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][1])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 2.587, np.log10(stellarfloat[int(round(fe_star))][2]*(Ca_dm/Mg_dm)*(Ca_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][2])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][2])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][2])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.169, np.log10(stellarfloat[int(round(fe_star))][3]*(Ni_dm/Mg_dm)*(Ni_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][3])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][3])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][3])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.289, np.log10(stellarfloat[int(round(fe_star))][4]*(Fe_dm/Mg_dm)*(Fe_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][4])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][4])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][4])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.047, np.log10(stellarfloat[int(round(fe_star))][5]*(Cr_dm/Mg_dm)*(Cr_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][5])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][5])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][5])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.868, np.log10(stellarfloat[int(round(fe_star))][6]*(Si_dm/Mg_dm)*(Si_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][7])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][7])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][7])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.011, np.log10(stellarfloat[int(round(fe_star))][7]*(Na_dm/Mg_dm)*(Na_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][8])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][8])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][8])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.346, np.log10(stellarfloat[int(round(fe_star))][8]*(O_dm/Mg_dm)*(O_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][9])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][9])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][9])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 1.2284, np.log10(stellarfloat[int(round(fe_star))][9]*(C_dm/Mg_dm)*(C_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][10])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][10])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][10])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6])))- 0.8507, np.log10(stellarfloat[int(round(fe_star))][10]*(Nz_dm/Mg_dm)*(Nz_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][11])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][11])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][11])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 0.3113]    
                                         elif np.log10(stellarfloat[int(round(fe_star))][9]*(C_dm/Mg_dm)*(C_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][10])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][10])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][10])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6])))- 0.8507 > -2 and np.log10(stellarfloat[int(round(fe_star))][10]*(Nz_dm/Mg_dm)*(Nz_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][11])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][11])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][11])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 0.3113 <= -2:
                                           YOUR_MODEL_OUTPUT_stored[i,:] = [np.log10(stellarfloat[int(round(fe_star))][0]*(Al_dm/Mg_dm)*(Al_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][0])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][0])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][0])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.145, np.log10(stellarfloat[int(round(fe_star))][1]*(Ti_dm/Mg_dm)*(Ti_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][1])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][1])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][1])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 2.587, np.log10(stellarfloat[int(round(fe_star))][2]*(Ca_dm/Mg_dm)*(Ca_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][2])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][2])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][2])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.169, np.log10(stellarfloat[int(round(fe_star))][3]*(Ni_dm/Mg_dm)*(Ni_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][3])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][3])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][3])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.289, np.log10(stellarfloat[int(round(fe_star))][4]*(Fe_dm/Mg_dm)*(Fe_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][4])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][4])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][4])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.047, np.log10(stellarfloat[int(round(fe_star))][5]*(Cr_dm/Mg_dm)*(Cr_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][5])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][5])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][5])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.868, np.log10(stellarfloat[int(round(fe_star))][6]*(Si_dm/Mg_dm)*(Si_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][7])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][7])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][7])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.011, np.log10(stellarfloat[int(round(fe_star))][7]*(Na_dm/Mg_dm)*(Na_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][8])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][8])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][8])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.346, np.log10(stellarfloat[int(round(fe_star))][8]*(O_dm/Mg_dm)*(O_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][9])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][9])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][9])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 1.2284, np.log10(stellarfloat[int(round(fe_star))][9]*(C_dm/Mg_dm)*(C_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][10])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][10])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][10])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6])))- 0.8507, -2] 
                                         elif np.log10(stellarfloat[int(round(fe_star))][9]*(C_dm/Mg_dm)*(C_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][10])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][10])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][10])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6])))- 0.8507 <= -2 and np.log10(stellarfloat[int(round(fe_star))][10]*(Nz_dm/Mg_dm)*(Nz_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][11])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][11])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][11])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 0.3113 > -2:
                                           YOUR_MODEL_OUTPUT_stored[i,:] = [np.log10(stellarfloat[int(round(fe_star))][0]*(Al_dm/Mg_dm)*(Al_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][0])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][0])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][0])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.145, np.log10(stellarfloat[int(round(fe_star))][1]*(Ti_dm/Mg_dm)*(Ti_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][1])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][1])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][1])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 2.587, np.log10(stellarfloat[int(round(fe_star))][2]*(Ca_dm/Mg_dm)*(Ca_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][2])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][2])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][2])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.169, np.log10(stellarfloat[int(round(fe_star))][3]*(Ni_dm/Mg_dm)*(Ni_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][3])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][3])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][3])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.289, np.log10(stellarfloat[int(round(fe_star))][4]*(Fe_dm/Mg_dm)*(Fe_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][4])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][4])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][4])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.047, np.log10(stellarfloat[int(round(fe_star))][5]*(Cr_dm/Mg_dm)*(Cr_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][5])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][5])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][5])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.868, np.log10(stellarfloat[int(round(fe_star))][6]*(Si_dm/Mg_dm)*(Si_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][7])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][7])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][7])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.011, np.log10(stellarfloat[int(round(fe_star))][7]*(Na_dm/Mg_dm)*(Na_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][8])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][8])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][8])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.346, np.log10(stellarfloat[int(round(fe_star))][8]*(O_dm/Mg_dm)*(O_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][9])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][9])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][9])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 1.2284, -2, np.log10(stellarfloat[int(round(fe_star))][10]*(Nz_dm/Mg_dm)*(Nz_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][11])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][11])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][11])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 0.3113]    
                                         else:
                                           YOUR_MODEL_OUTPUT_stored[i,:] = [np.log10(stellarfloat[int(round(fe_star))][0]*(Al_dm/Mg_dm)*(Al_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][0])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][0])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][0])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.145, np.log10(stellarfloat[int(round(fe_star))][1]*(Ti_dm/Mg_dm)*(Ti_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][1])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][1])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][1])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 2.587, np.log10(stellarfloat[int(round(fe_star))][2]*(Ca_dm/Mg_dm)*(Ca_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][2])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][2])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][2])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.169, np.log10(stellarfloat[int(round(fe_star))][3]*(Ni_dm/Mg_dm)*(Ni_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][3])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][3])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][3])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.289, np.log10(stellarfloat[int(round(fe_star))][4]*(Fe_dm/Mg_dm)*(Fe_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][4])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][4])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][4])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.047, np.log10(stellarfloat[int(round(fe_star))][5]*(Cr_dm/Mg_dm)*(Cr_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][5])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][5])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][5])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.868, np.log10(stellarfloat[int(round(fe_star))][6]*(Si_dm/Mg_dm)*(Si_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][7])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][7])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][7])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 0.011, np.log10(stellarfloat[int(round(fe_star))][7]*(Na_dm/Mg_dm)*(Na_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][8])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][8])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][8])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) + 1.346, np.log10(stellarfloat[int(round(fe_star))][8]*(O_dm/Mg_dm)*(O_hsc/Mg_hsc)*((vars()['whitedwarftimescales'+str(N)][9])/(vars()['whitedwarftimescales'+str(N)][6]))*((1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][9])))/(1-np.exp(-(t_disc)/(vars()['whitedwarftimescales'+str(N)][6]))))*np.exp(-(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][9])+(((t_sinceaccretion-(t_disc/1000000))*1000000)/vars()['whitedwarftimescales'+str(N)][6]))) - 1.2284, -2, -2] 
    
model_low3, model_low2, model_low1, model_median, model_high1, model_high2, model_high3 = confidence_intervals(sample_draws, YOUR_MODEL_OUTPUT_stored, len(x_axis))
plt.plot([], label=str(Names[N])+' Observational Data', color='k', marker='o', linestyle='none')
plt.plot(x_axis, model_median, color='k', lw=2, label='Median Model')
plt.fill_between(x_axis,model_high1,model_low1, edgecolor='none', color='green', alpha = 0.5, label='1 Sigma Confidence Interval',zorder=2)
if vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, whitedwarfabundances[0]-whitedwarfabundances[5]+2.587, whitedwarfabundances[1]-whitedwarfabundances[5]+1.169, whitedwarfabundances[2]-whitedwarfabundances[5]+1.289, whitedwarfabundances[3]-whitedwarfabundances[5]+0.047, whitedwarfabundances[4]-whitedwarfabundances[5]+1.868, np.nan, whitedwarfabundances[6]-whitedwarfabundances[5]+1.346, np.nan, np.nan, np.nan]
        whitedwarferrorsplot=[0, vars()['errors'+str(N)][1], vars()['errors'+str(N)][2], vars()['errors'+str(N)][3], vars()['errors'+str(N)][4], vars()['errors'+str(N)][5], 0, vars()['errors'+str(N)][7], 0, 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, whitedwarfabundances[0]-whitedwarfabundances[5]+2.587, whitedwarfabundances[1]-whitedwarfabundances[5]+1.169, whitedwarfabundances[2]-whitedwarfabundances[5]+1.289, whitedwarfabundances[3]-whitedwarfabundances[5]+0.047, whitedwarfabundances[4]-whitedwarfabundances[5]+1.868, np.nan, np.nan, np.nan, np.nan, np.nan]
        whitedwarferrorsplot=[0, vars()['errors'+str(N)][1], vars()['errors'+str(N)][2], vars()['errors'+str(N)][3], vars()['errors'+str(N)][4], vars()['errors'+str(N)][5], 0, 0, 0, 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, np.nan, whitedwarfabundances[0]-whitedwarfabundances[4]+1.169, whitedwarfabundances[1]-whitedwarfabundances[4]+1.289, whitedwarfabundances[2]-whitedwarfabundances[4]+0.047, whitedwarfabundances[3]-whitedwarfabundances[4]+1.868, np.nan, whitedwarfabundances[5]-whitedwarfabundances[4]+1.346, np.nan, np.nan, np.nan]
        whitedwarferrorsplot=[0, 0, vars()['errors'+str(N)][2], vars()['errors'+str(N)][3], vars()['errors'+str(N)][4], vars()['errors'+str(N)][5], 0, vars()['errors'+str(N)][7], 0, 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, whitedwarfabundances[0]-whitedwarfabundances[4]+2.587, whitedwarfabundances[1]-whitedwarfabundances[4]+1.169, np.nan, whitedwarfabundances[2]-whitedwarfabundances[4]+0.047, whitedwarfabundances[3]-whitedwarfabundances[4]+1.868, np.nan, whitedwarfabundances[5]-whitedwarfabundances[4]+1.346, np.nan, np.nan, np.nan]
        whitedwarferrorsplot=[0, vars()['errors'+str(N)][1], vars()['errors'+str(N)][2], 0, vars()['errors'+str(N)][4], vars()['errors'+str(N)][5], 0, vars()['errors'+str(N)][7], 0, 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, whitedwarfabundances[0]-whitedwarfabundances[4]+2.587, whitedwarfabundances[1]-whitedwarfabundances[4]+1.169, whitedwarfabundances[2]-whitedwarfabundances[4]+1.289, whitedwarfabundances[3]-whitedwarfabundances[4]+0.047, np.nan, np.nan, whitedwarfabundances[5]-whitedwarfabundances[4]+1.346, np.nan, np.nan, np.nan]
        whitedwarferrorsplot=[0, vars()['errors'+str(N)][1], vars()['errors'+str(N)][2], vars()['errors'+str(N)][3], vars()['errors'+str(N)][4], 0, 0,vars()['errors'+str(N)][7], 0, 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, np.nan, whitedwarfabundances[0]-whitedwarfabundances[3]+1.169, np.nan, whitedwarfabundances[1]-whitedwarfabundances[3]+0.047, whitedwarfabundances[2]-whitedwarfabundances[3]+1.868, np.nan, whitedwarfabundances[4]-whitedwarfabundances[3]+1.346, np.nan, np.nan, np.nan]
        whitedwarferrorsplot=[0, 0, vars()['errors'+str(N)][2], 0, vars()['errors'+str(N)][4], vars()['errors'+str(N)][5], 0, vars()['errors'+str(N)][7], 0, 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, whitedwarfabundances[0]-whitedwarfabundances[3]+2.587, whitedwarfabundances[1]-whitedwarfabundances[3]+1.169, np.nan, whitedwarfabundances[2]-whitedwarfabundances[3]+0.047, np.nan, np.nan, whitedwarfabundances[4]-whitedwarfabundances[3]+1.346, np.nan, np.nan, np.nan]
        whitedwarferrorsplot=[0, vars()['errors'+str(N)][1], vars()['errors'+str(N)][2], 0, vars()['errors'+str(N)][4], 0, 0, vars()['errors'+str(N)][7], 0, 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, np.nan, whitedwarfabundances[0]-whitedwarfabundances[3]+1.169, whitedwarfabundances[1]-whitedwarfabundances[3]+1.289, whitedwarfabundances[2]-whitedwarfabundances[3]+0.047, np.nan, np.nan, whitedwarfabundances[4]-whitedwarfabundances[3]+1.346, np.nan, np.nan, np.nan]
        whitedwarferrorsplot=[0, 0, vars()['errors'+str(N)][2], vars()['errors'+str(N)][3], vars()['errors'+str(N)][4], 0, 0, vars()['errors'+str(N)][7], 0, 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, whitedwarfabundances[0]-whitedwarfabundances[4]+2.587, whitedwarfabundances[1]-whitedwarfabundances[4]+1.169, whitedwarfabundances[2]-whitedwarfabundances[4]+1.289, whitedwarfabundances[3]-whitedwarfabundances[4]+0.047, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
        whitedwarferrorsplot=[0, vars()['errors'+str(N)][1], vars()['errors'+str(N)][2], vars()['errors'+str(N)][3], vars()['errors'+str(N)][4], 0, 0, 0, 0, 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, np.nan, whitedwarfabundances[0]-whitedwarfabundances[4]+1.169, whitedwarfabundances[1]-whitedwarfabundances[4]+1.289, whitedwarfabundances[2]-whitedwarfabundances[4]+0.047, whitedwarfabundances[3]-whitedwarfabundances[4]+1.868, np.nan, np.nan, np.nan, np.nan, np.nan]
        whitedwarferrorsplot=[0, 0, vars()['errors'+str(N)][2], vars()['errors'+str(N)][3], vars()['errors'+str(N)][4], vars()['errors'+str(N)][5], 0, 0, 0, 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, whitedwarfabundances[0]-whitedwarfabundances[4]+2.587, whitedwarfabundances[1]-whitedwarfabundances[4]+1.169, np.nan, whitedwarfabundances[2]-whitedwarfabundances[4]+0.047, whitedwarfabundances[3]-whitedwarfabundances[4]+1.868, np.nan, np.nan, np.nan, np.nan, np.nan]
        whitedwarferrorsplot=[0, vars()['errors'+str(N)][1], vars()['errors'+str(N)][2],0, vars()['errors'+str(N)][4], vars()['errors'+str(N)][5], 0, 0, 0, 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
            whitedwarfabundancesplot=[np.nan, np.nan, whitedwarfabundances[0]-whitedwarfabundances[3]+1.169, whitedwarfabundances[1]-whitedwarfabundances[3]+1.289, whitedwarfabundances[2]-whitedwarfabundances[3]+0.047,np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
            whitedwarferrorsplot=[0, 0, vars()['errors'+str(N)][2], vars()['errors'+str(N)][3], vars()['errors'+str(N)][4], 0, 0, 0, 0, 0, 0]
            plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, np.nan, whitedwarfabundances[0]-whitedwarfabundances[2]+1.169, np.nan, whitedwarfabundances[1]-whitedwarfabundances[2]+0.047, np.nan, np.nan, whitedwarfabundances[3]-whitedwarfabundances[2]+1.346, np.nan, np.nan, np.nan]
        whitedwarferrorsplot=[0, 0, vars()['errors'+str(N)][2], 0, vars()['errors'+str(N)][4], 0, 0, vars()['errors'+str(N)][7], 0, 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, np.nan, whitedwarfabundances[0]-whitedwarfabundances[3]+1.169, np.nan, whitedwarfabundances[1]-whitedwarfabundances[3]+0.047, whitedwarfabundances[2]-whitedwarfabundances[3]+1.868, np.nan, np.nan, np.nan, np.nan, np.nan]
        whitedwarferrorsplot=[0, 0, vars()['errors'+str(N)][2], 0, vars()['errors'+str(N)][4], vars()['errors'+str(N)][5], 0, 0, 0, 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, whitedwarfabundances[0]-whitedwarfabundances[3]+2.587, whitedwarfabundances[1]-whitedwarfabundances[3]+1.169, np.nan, whitedwarfabundances[2]-whitedwarfabundances[3]+0.047, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
        whitedwarferrorsplot=[0, vars()['errors'+str(N)][1], vars()['errors'+str(N)][2], 0, vars()['errors'+str(N)][4], 0, 0, 0, 0, 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] == 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, np.nan, whitedwarfabundances[0]-whitedwarfabundances[2]+1.169, np.nan, whitedwarfabundances[1]-whitedwarfabundances[2]+0.047, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
        whitedwarferrorsplot=[0, 0, vars()['errors'+str(N)][2], 0, vars()['errors'+str(N)][4], 0, 0, 0, 0, 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[whitedwarfabundances[0]-whitedwarfabundances[6]+1.145, whitedwarfabundances[1]-whitedwarfabundances[6]+2.587, whitedwarfabundances[2]-whitedwarfabundances[6]+1.169, whitedwarfabundances[3]-whitedwarfabundances[6]+1.289, whitedwarfabundances[4]-whitedwarfabundances[6]+0.047, whitedwarfabundances[5]-whitedwarfabundances[6]+1.868, whitedwarfabundances[7]-whitedwarfabundances[6]+0.011, whitedwarfabundances[8]-whitedwarfabundances[6]+1.346, whitedwarfabundances[9]-whitedwarfabundances[6]-1.2284, np.nan, np.nan]
        whitedwarferrorsplot=[vars()['errors'+str(N)][0], vars()['errors'+str(N)][1], vars()['errors'+str(N)][2], vars()['errors'+str(N)][3], vars()['errors'+str(N)][4], vars()['errors'+str(N)][5], vars()['errors'+str(N)][6], vars()['errors'+str(N)][7], vars()['errors'+str(N)][8], 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[whitedwarfabundances[0]-whitedwarfabundances[6]+1.145, whitedwarfabundances[1]-whitedwarfabundances[6]+2.587, whitedwarfabundances[2]-whitedwarfabundances[6]+1.169, whitedwarfabundances[3]-whitedwarfabundances[6]+1.289, whitedwarfabundances[4]-whitedwarfabundances[6]+0.047, whitedwarfabundances[5]-whitedwarfabundances[6]+1.868, whitedwarfabundances[7]-whitedwarfabundances[6]+0.011, np.nan, whitedwarfabundances[8]-whitedwarfabundances[6]-1.2284, np.nan, np.nan]
        whitedwarferrorsplot=[vars()['errors'+str(N)][0], vars()['errors'+str(N)][1], vars()['errors'+str(N)][2], vars()['errors'+str(N)][3], vars()['errors'+str(N)][4], vars()['errors'+str(N)][5], vars()['errors'+str(N)][6], 0, vars()['errors'+str(N)][8], 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')       
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[whitedwarfabundances[0]-whitedwarfabundances[6]+1.145, whitedwarfabundances[1]-whitedwarfabundances[6]+2.587, whitedwarfabundances[2]-whitedwarfabundances[6]+1.169, whitedwarfabundances[3]-whitedwarfabundances[6]+1.289, whitedwarfabundances[4]-whitedwarfabundances[6]+0.047, whitedwarfabundances[5]-whitedwarfabundances[6]+1.868, whitedwarfabundances[7]-whitedwarfabundances[6]+0.011, whitedwarfabundances[8]-whitedwarfabundances[6]+1.346, np.nan, np.nan, np.nan]
        whitedwarferrorsplot=[vars()['errors'+str(N)][0], vars()['errors'+str(N)][1], vars()['errors'+str(N)][2], vars()['errors'+str(N)][3], vars()['errors'+str(N)][4], vars()['errors'+str(N)][5], vars()['errors'+str(N)][6], vars()['errors'+str(N)][7], 0, 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] != 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, whitedwarfabundances[0]-whitedwarfabundances[4]+2.587, whitedwarfabundances[1]-whitedwarfabundances[4]+1.169, np.nan, whitedwarfabundances[2]-whitedwarfabundances[4]+0.047, whitedwarfabundances[3]-whitedwarfabundances[4]+1.868, whitedwarfabundances[5]-whitedwarfabundances[4]+0.011, whitedwarfabundances[6]-whitedwarfabundances[4]+1.346, whitedwarfabundances[7]-whitedwarfabundances[4]-1.2284, np.nan, np.nan]
        whitedwarferrorsplot=[0, vars()['errors'+str(N)][1], vars()['errors'+str(N)][2], 0, vars()['errors'+str(N)][4], vars()['errors'+str(N)][5], vars()['errors'+str(N)][6], vars()['errors'+str(N)][7], vars()['errors'+str(N)][8], 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[whitedwarfabundances[0]-whitedwarfabundances[5]+1.145, whitedwarfabundances[1]-whitedwarfabundances[5]+2.587, whitedwarfabundances[2]-whitedwarfabundances[5]+1.169, np.nan, whitedwarfabundances[3]-whitedwarfabundances[5]+0.047, whitedwarfabundances[4]-whitedwarfabundances[5]+1.868, whitedwarfabundances[6]-whitedwarfabundances[5]+0.011, np.nan, whitedwarfabundances[7]-whitedwarfabundances[5]-1.2284, np.nan, np.nan]
        whitedwarferrorsplot=[vars()['errors'+str(N)][0], vars()['errors'+str(N)][1], vars()['errors'+str(N)][2], 0, vars()['errors'+str(N)][4], vars()['errors'+str(N)][5], vars()['errors'+str(N)][6], 0, vars()['errors'+str(N)][8], 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, whitedwarfabundances[0]-whitedwarfabundances[5]+2.587, whitedwarfabundances[1]-whitedwarfabundances[5]+1.169, whitedwarfabundances[2]-whitedwarfabundances[5]+1.289, whitedwarfabundances[3]-whitedwarfabundances[5]+0.047, whitedwarfabundances[4]-whitedwarfabundances[5]+1.868, whitedwarfabundances[6]-whitedwarfabundances[5]+0.011, np.nan, whitedwarfabundances[7]-whitedwarfabundances[5]-1.2284, np.nan, np.nan]
        whitedwarferrorsplot=[0, vars()['errors'+str(N)][1], vars()['errors'+str(N)][2], vars()['errors'+str(N)][3], vars()['errors'+str(N)][4], vars()['errors'+str(N)][5], vars()['errors'+str(N)][6], 0, vars()['errors'+str(N)][8], 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[whitedwarfabundances[0]-whitedwarfabundances[5]+1.145, np.nan, whitedwarfabundances[1]-whitedwarfabundances[5]+1.169, whitedwarfabundances[2]-whitedwarfabundances[5]+1.289, whitedwarfabundances[3]-whitedwarfabundances[5]+0.047, whitedwarfabundances[4]-whitedwarfabundances[5]+1.868, whitedwarfabundances[6]-whitedwarfabundances[5]+0.011, np.nan, whitedwarfabundances[7]-whitedwarfabundances[5]-1.2284, np.nan, np.nan]
        whitedwarferrorsplot=[vars()['errors'+str(N)][0], 0, vars()['errors'+str(N)][2], vars()['errors'+str(N)][3], vars()['errors'+str(N)][4], vars()['errors'+str(N)][5], vars()['errors'+str(N)][6], 0, vars()['errors'+str(N)][8], 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] == 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, whitedwarfabundances[0]-whitedwarfabundances[5]+2.587, whitedwarfabundances[1]-whitedwarfabundances[5]+1.169, whitedwarfabundances[2]-whitedwarfabundances[5]+1.289, whitedwarfabundances[3]-whitedwarfabundances[5]+0.047, whitedwarfabundances[4]-whitedwarfabundances[5]+1.868, whitedwarfabundances[6]-whitedwarfabundances[5]+0.011, np.nan, np.nan, np.nan, np.nan]
        whitedwarferrorsplot=[0, vars()['errors'+str(N)][1], vars()['errors'+str(N)][2], vars()['errors'+str(N)][3], vars()['errors'+str(N)][4], vars()['errors'+str(N)][5], vars()['errors'+str(N)][6], 0, 0, 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] == 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[whitedwarfabundances[0]-whitedwarfabundances[4]+1.145, np.nan, np.nan, whitedwarfabundances[1]-whitedwarfabundances[4]+1.289, whitedwarfabundances[2]-whitedwarfabundances[4]+0.047, whitedwarfabundances[3]-whitedwarfabundances[4]+1.868, whitedwarfabundances[5]-whitedwarfabundances[4]+0.011, np.nan, whitedwarfabundances[6]-whitedwarfabundances[4]-1.2284, np.nan, np.nan]
        whitedwarferrorsplot=[vars()['errors'+str(N)][0], 0, 0, vars()['errors'+str(N)][3], vars()['errors'+str(N)][4], vars()['errors'+str(N)][5], vars()['errors'+str(N)][6],0, vars()['errors'+str(N)][8], 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] != 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] != 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, whitedwarfabundances[0]-whitedwarfabundances[4]+2.587, whitedwarfabundances[1]-whitedwarfabundances[4]+1.169, np.nan, whitedwarfabundances[2]-whitedwarfabundances[4]+0.047, whitedwarfabundances[3]-whitedwarfabundances[4]+1.868, whitedwarfabundances[5]-whitedwarfabundances[4]+0.011, np.nan, whitedwarfabundances[6]-whitedwarfabundances[4]-1.2284, np.nan, np.nan]
        whitedwarferrorsplot=[0, vars()['errors'+str(N)][1], vars()['errors'+str(N)][2], 0, vars()['errors'+str(N)][4], vars()['errors'+str(N)][5], vars()['errors'+str(N)][6], 0, vars()['errors'+str(N)][8], 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[whitedwarfabundances[0]-whitedwarfabundances[4]+1.145, np.nan, whitedwarfabundances[1]-whitedwarfabundances[4]+1.169, whitedwarfabundances[2]-whitedwarfabundances[4]+1.289, whitedwarfabundances[3]-whitedwarfabundances[4]+0.047, np.nan, whitedwarfabundances[5]-whitedwarfabundances[4]+0.011, np.nan, whitedwarfabundances[6]-whitedwarfabundances[4]-1.2284, np.nan, np.nan]
        whitedwarferrorsplot=[vars()['errors'+str(N)][0], 0, vars()['errors'+str(N)][2], vars()['errors'+str(N)][3], vars()['errors'+str(N)][4], 0, vars()['errors'+str(N)][6], 0, vars()['errors'+str(N)][8], 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] != 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[whitedwarfabundances[0]-whitedwarfabundances[3]+1.145, np.nan, whitedwarfabundances[1]-whitedwarfabundances[3]+1.169, np.nan, whitedwarfabundances[2]-whitedwarfabundances[3]+0.047, np.nan, whitedwarfabundances[4]-whitedwarfabundances[3]+0.011, np.nan, whitedwarfabundances[5]-whitedwarfabundances[3]-1.2284, np.nan, np.nan]
        whitedwarferrorsplot=[vars()['errors'+str(N)][0], 0, vars()['errors'+str(N)][2], 0, vars()['errors'+str(N)][4], 0, vars()['errors'+str(N)][6], 0, vars()['errors'+str(N)][8], 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] == 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] == 0 and vars()['whitedwarfabundances'+str(N)][11] == 0:
        whitedwarfabundancesplot=[np.nan, np.nan, whitedwarfabundances[0]-whitedwarfabundances[2]+1.169, np.nan, whitedwarfabundances[1]-whitedwarfabundances[2]+0.047, np.nan, whitedwarfabundances[3]-whitedwarfabundances[2]+0.011, np.nan, whitedwarfabundances[4]-whitedwarfabundances[2]-1.2284, np.nan, np.nan]
        whitedwarferrorsplot=[0, 0, vars()['errors'+str(N)][2], 0, vars()['errors'+str(N)][4], 0, vars()['errors'+str(N)][6], 0, vars()['errors'+str(N)][8], 0, 0]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')
elif vars()['whitedwarfabundances'+str(N)][0] == 0 and vars()['whitedwarfabundances'+str(N)][1] == 0 and vars()['whitedwarfabundances'+str(N)][2] != 0 and vars()['whitedwarfabundances'+str(N)][3] != 0 and vars()['whitedwarfabundances'+str(N)][4] != 0 and vars()['whitedwarfabundances'+str(N)][5] == 0 and vars()['whitedwarfabundances'+str(N)][6] != 0 and vars()['whitedwarfabundances'+str(N)][7] != 0 and vars()['whitedwarfabundances'+str(N)][8] == 0 and vars()['whitedwarfabundances'+str(N)][9] != 0 and vars()['whitedwarfabundances'+str(N)][10] != 0 and vars()['whitedwarfabundances'+str(N)][11] != 0:
        whitedwarfabundancesplot=[np.nan, np.nan, whitedwarfabundances[0]-whitedwarfabundances[3]+1.169, whitedwarfabundances[1]-whitedwarfabundances[3]+1.289, whitedwarfabundances[2]-whitedwarfabundances[3]+0.047, np.nan, whitedwarfabundances[4]-whitedwarfabundances[3]+0.011, np.nan, whitedwarfabundances[5]-whitedwarfabundances[3]-1.2284, whitedwarfabundances[6]-whitedwarfabundances[3]-0.8507, whitedwarfabundances[7]-whitedwarfabundances[3]-0.3113]
        whitedwarferrorsplot=[0, 0, vars()['errors'+str(N)][2], vars()['errors'+str(N)][3], vars()['errors'+str(N)][4], 0, vars()['errors'+str(N)][6], 0, vars()['errors'+str(N)][8], vars()['errors'+str(N)][9], vars()['errors'+str(N)][10]]
        plt.errorbar(x_axis, whitedwarfabundancesplot, yerr=whitedwarferrorsplot, fmt=".k", capsize=2, capthick=1, marker='o')                                
else:
        print("Error: Failed to plot white dwarf observational data")
upper = [0.114145511, 0.131417374, 0.152552378, 0.138930603, 0.151802017, 0.132419952, 0.095546878, 0.230676051, 0.281268694, 0.172593095, 0.318840642]
lower = [-0.300126301, -0.062910874, -0.13177587, -0.161069397, -0.252526231, -0.251908295, -0.12878137, -0.203652197, -0.228731306, -0.234382055, -0.321159357]
solar = [-0.015854489, -0.042910874, -0.05177587, -0.011069397, -0.032526231, -0.021908295, -0.00878137, -0.013652197, -0.098731306, 0.008264847, -0.061159358]
plt.fill_between(x_axis,upper,lower, color='silver', label='Stellar 98\% Composition Range',zorder=1)
plt.xlabel('Lithophiles \, \, \, \, \, \, \, \, \, \, \, \, Siderophiles \, \, \, \, \, \, \, Volatile Lithophiles \, \, \, \, Atmophiles', fontsize=14, horizontalalignment='left', x=0.08)
plt.ylabel('$\mathrm{log((X/Mg)/(X/Mg)}_{\mathrm{mean \, stellar}})$', fontsize=16)
plt.yticks(np.arange(-2, 2.5, step=0.5))
plt.plot(x_axis, solar, label='Solar Composition', color='k', linestyle='--')
plt.errorbar(x_axis, [0.09,0.72,np.nan,-0.77,np.nan,0.61,np.nan,np.nan,np.nan,np.nan,-1.57], yerr=0.15, uplims=True, linestyle='none', marker='o', color='k')
#plt.errorbar(x_axis, [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan], yerr=0.2, uplims=True, linestyle='none', marker='o', color='k')
#\ \ \ \ \ \ \ \Delta \textrm{Crust} = -0.12^{+0.02}_{-0.02}
#    'r'$\textrm{Steady State Accretion required to}\ 2.15\sigma$'  + '\n' +
#+ '\n' + r'$\textrm{Crustal Differentiation required to}\ 2.18\sigma$'+ '\n' 
# + '\n' + r'$\textrm{Steady State Accretion required to}\ 1.16\sigma$'+ '\n' +
plt.text(5,-1.78,r'$\Delta t = -0.01^{+0.86}_{-0.04} \textrm{Myrs} \ \ \ \ \ \ \ T = 157^{+26}_{-19} \textrm{K}\ \ \ \ \ \ \ \Delta \textrm{Core} = -0.16^{+0.01}_{-0.01}\ \ \ \ \ \ \ \textrm{log}(M) = 18.49^{+0.07}_{-0.29} \textrm{kg}   $', bbox=dict(edgecolor='k',facecolor = 'white', linewidth=1.5 ), ha='center')
plt.text(-0.25,1.82, r'$\textrm{Water Ice required to}\ 2.43\sigma$' + '\n' + r'$\textrm{Heating required to}\ 4.91\sigma$' + '\n' + r'$\textrm{Core Differentiation required to}\ 8.44\sigma $' + '\n' + r'$\textrm{Average} \,\chi^{2}\,\textrm{per data point} = 0.05$' , bbox=dict(edgecolor='k',facecolor = 'white', linewidth=1.5 ), va = 'top')
plt.ylim(-2,2)
plt.legend(frameon = False)
plt.savefig(prefix + 'PWDComposition.pdf')

#norm_log = (-0.5*np.log(2.0*np.pi*whitedwarferrors*whitedwarferrors)).sum()
#retrieved_M1 = pymultinest.Analyzer(n_params = n_params, outputfiles_basename = prefix)
#s_M1 = retrieved_M1.get_stats()
#best_fit_M1 = retrieved_M1.get_best_fit()
#max_likelihood_M1 = best_fit_M1['log_likelihood']                                
#best_chi_square_M1 = -2.0 * (max_likelihood_M1 - norm_log)
#print(best_chi_square_M1)
