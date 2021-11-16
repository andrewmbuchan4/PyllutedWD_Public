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
import xlrd
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
  
def discrete_cmap(N, base_cmap=None):
          base = plt.cm.get_cmap(base_cmap)

          color_list = base(np.linspace(0, 1, N))

          cmap_name = base.name + str(N)

          return base.from_list(cmap_name, color_list, N)

for i in range(0,208):
    location = 'chains/' + str(i) +'PWDOutputs.xlsx'
    workbook = xlrd.open_workbook(location)
    worksheet = workbook.sheet_by_index(0)
    vars()['bins'+str(i)]  = [0]*121
    for j in range(0,121):
       vars()['bins'+str(i)][j] = worksheet.cell(12+j,4).value
for i in range(0,208):
    vars()['bins'+str(i)] = np.asarray(vars()['bins'+str(i)])
totalbins = [0]*121
for i in range(0,208):
    totalbins = totalbins + vars()['bins'+str(i)]
totalbins = totalbins/208 
xbar = [0,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600,625,650,675,700,725,750,775,800,825,850,875,900,925,950,975,1000,1025,1050,1075,1100,1125,1150,1175,1200,1225,1250,1275,1300,1325,1350,1375,1400,1425,1450,1475,1500,1525,1550,1575,1600,1625,1650,1675,1700,1725,1750,1775,1800,1825,1850,1875,1900,1925,1950,1975,2000,2025,2050,2075,2100,2125,2150,2175,2200,2225,2250,2275,2300,2325,2350,2375,2400,2425,2450,2475,2500,2525,2550,2575,2600,2625,2650,2675,2700,2725,2750,2775,2800,2825,2850,2875,2900,2925,2950,2975,3000]
plt.figure(figsize=(6,4))
plt.bar(xbar,totalbins,width=25,color='k',edgecolor='k',label='Hollands et al. 2017 data')
plt.axvline(x=220,linestyle='dotted',color='grey',lw=2)   
plt.axvline(x=1250,linestyle='dotted',color='grey',lw=2) 
plt.text(0,(0.25*1.1),'Icy', ha='center')  
plt.text(700,(0.25*1.1),'Dry', ha='center') 
plt.text(2175,(0.25*1.1),'Extreme Heating', ha='center')
plt.ylim(0,(0.3))
plt.xlabel('Formation Temperature/K', fontsize=14)
plt.ylabel('Probability', fontsize=14)
plt.legend(frameon = False, loc=4)
plt.savefig('PopulationFormationTemperatureplot.pdf')


for i in range(0,208):
    location = 'chains/' + str(i) +'PWDOutputs.xlsx'
    workbook = xlrd.open_workbook(location)
    worksheet = workbook.sheet_by_index(0)
    vars()['bins'+str(i)]  = [0]*201
    for j in range(0,201):
       vars()['bins'+str(i)][j] = worksheet.cell(12+j,6).value
for i in range(0,208):
    vars()['bins'+str(i)] = np.asarray(vars()['bins'+str(i)])
totalbins = [0]*201
for i in range(0,208):
    totalbins = totalbins + vars()['bins'+str(i)]
totalbins = totalbins/208 
xbar = [-100,-99,-98,-97,-96,-95,-94,-93,-92,-91,-90,-89,-88,-87,-86,-85,-84,-83,-82,-81,-80,-79,-78,-77,-76,-75,-74,-73,-72,-71,-70,-69,-68,-67,-66,-65,-64,-63,-62,-61,-60,-59,-58,-57,-56,-55,-54,-53,-52,-51,-50,-49,-48,-47,-46,-45,-44,-43,-42,-41,-40,-39,-38,-37,-36,-35,-34,-33,-32,-31,-30,-29,-28,-27,-26,-25,-24,-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100]
plt.figure(figsize=(6,4))
plt.bar(xbar,totalbins,width=1,color='k',edgecolor='k',label='Hollands et al. 2017 data')
plt.xlabel('\% \ \ \ \ \ Core \ \ \ \ \ Lost \ \ \ \ \ \ \ \ \ \ \ \ \ \ \% Mantle+Crust Lost', fontsize=14)
plt.ylabel('Probability', fontsize=14)
plt.ylim(0,(0.035))
plt.text(-56,(0.032),'Core Depleted, Mantle Enhanced', ha='center') 
plt.text(56,(0.032),'Core Enhanced, Mantle Depleted', ha='center')                                            
plt.legend(frameon = False, fontsize=9, loc=7)
plt.axvline(x=0,linestyle='dotted',color='grey',lw=2)
plt.savefig('PopulationCoreDifferentiationplot.pdf')


for i in range(0,208):
    location = 'chains/' + str(i) +'PWDOutputs.xlsx'
    workbook = xlrd.open_workbook(location)
    worksheet = workbook.sheet_by_index(0)
    vars()['bins'+str(i)]  = [0]*201
    for j in range(0,201):
       vars()['bins'+str(i)][j] = worksheet.cell(12+j,8).value
for i in range(0,208):
    vars()['bins'+str(i)] = np.asarray(vars()['bins'+str(i)])
totalbins = [0]*201
for i in range(0,208):
    totalbins = totalbins + vars()['bins'+str(i)]
totalbins = totalbins/208 
xbar = [-100,-99,-98,-97,-96,-95,-94,-93,-92,-91,-90,-89,-88,-87,-86,-85,-84,-83,-82,-81,-80,-79,-78,-77,-76,-75,-74,-73,-72,-71,-70,-69,-68,-67,-66,-65,-64,-63,-62,-61,-60,-59,-58,-57,-56,-55,-54,-53,-52,-51,-50,-49,-48,-47,-46,-45,-44,-43,-42,-41,-40,-39,-38,-37,-36,-35,-34,-33,-32,-31,-30,-29,-28,-27,-26,-25,-24,-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100]
plt.figure(figsize=(6,4))
plt.bar(xbar,totalbins,width=1,color='k',edgecolor='k',label='Hollands et al. 2017 data')
plt.xlabel(' \% \ \ \ \ \ Crust \ \ \ \ \ Lost  \ \ \ \ \ \ \ \ \ \ \ \ \% Mantle+Core Lost', fontsize=14)
plt.ylabel('Probability', fontsize=14)
plt.ylim(0,(0.05*1.25))
plt.text(-50,(0.05*1.05),'Crust Depleted', ha='center') 
plt.text(50,(0.05*1.05),'Crust Enhanced', ha='center')                              
plt.axvline(x=0,linestyle='dotted',color='grey',lw=2)
plt.legend(frameon = False, fontsize=9, loc=7)
plt.savefig('PopulationCrustDifferentiationplot.pdf')


for i in range(0,208):
    location = 'chains/' + str(i) +'PWDOutputs.xlsx'
    workbook = xlrd.open_workbook(location)
    worksheet = workbook.sheet_by_index(0)
    vars()['bins'+str(i)]  = [0]*171
    for j in range(0,171):
       vars()['bins'+str(i)][j] = worksheet.cell(12+j,10).value
for i in range(0,208):
    vars()['bins'+str(i)] = np.asarray(vars()['bins'+str(i)])
totalbins = [0]*171
for i in range(0,208):
    totalbins = totalbins + vars()['bins'+str(i)]
totalbins = totalbins/208
xbar = [8,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8,8.9,9,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,10,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11,11.1,11.2,11.3,11.4,11.5,11.6,11.7,11.8,11.9,12,12.1,12.2,12.3,12.4,12.5,12.6,12.7,12.8,12.9,13,13.1,13.2,13.3,13.4,13.5,13.6,13.7,13.8,13.9,14,14.1,14.2,14.3,14.4,14.5,14.6,14.7,14.8,14.9,15,15.1,15.2,15.3,15.4,15.5,15.6,15.7,15.8,15.9,16,16.1,16.2,16.3,16.4,16.5,16.6,16.7,16.8,16.9,17,17.1,17.2,17.3,17.4,17.5,17.6,17.7,17.8,17.9,18,18.1,18.2,18.3,18.4,18.5,18.6,18.7,18.8,18.9,19,19.1,19.2,19.3,19.4,19.5,19.6,19.7,19.8,19.9,20,20.1,20.2,20.3,20.4,20.5,20.6,20.7,20.8,20.9,21,21.1,21.2,21.3,21.4,21.5,21.6,21.7,21.8,21.9,22,22.1,22.2,22.3,22.4,22.5,22.6,22.7,22.8,22.9,23,23.1,23.2,23.3,23.4,23.5,23.6,23.7,23.8,23.9,24,24.1,24.2,24.3,24.4,24.5,24.6,24.7,24.8,24.9,25]
plt.figure(figsize=(6,4))
plt.bar(xbar,totalbins,width=0.10,color='k',edgecolor='k',label='Hollands et al. 2017 data')
plt.xlim(8,25) 
plt.ylim(0,(0.075*1.25))
plt.xticks(np.arange(8, 26, step=2))
plt.xlabel('Log(Mass of Pollutant/kg)', fontsize=14)
plt.ylabel('Probability', fontsize=14)
h = [0.075,0]
plt.text(24.78,(np.max(h)*1.2),'Earth', va='top', ha='right' ,rotation='vertical')
plt.text(23.81,(np.max(h)*1.2),'Mars', va='top', rotation='vertical')
plt.text(22.87,(np.max(h)*1.2),'Moon', va='top', rotation='vertical')
plt.text(22.11,(np.max(h)*1.2),'Pluto', va='top', rotation='vertical')
plt.text(20.97,(np.max(h)*1.2),'Ceres', va='top', rotation='vertical')
plt.text(20.41,(np.max(h)*1.2),'Vesta', va='top', rotation='vertical')
plt.text(19.94,(np.max(h)*1.2),'Hygiea', va='top', rotation='vertical')
plt.text(19.38,(np.max(h)*1.2),'Psyche', va='top', rotation='vertical')
plt.text(18.93,(np.max(h)*1.2),'Flora', va='top', rotation='vertical')
plt.text(17.72,(np.max(h)*1.2),'Epimetheus', va='top', rotation='vertical')
plt.text(16.72,(np.max(h)*1.2),'Ophelia', va='top', rotation='vertical')
plt.text(16.02,(np.max(h)*1.2),'Phobos', va='top', rotation='vertical') 
plt.text(15.17,(np.max(h)*1.2),'Deimos', va='top', rotation='vertical')
plt.text(14.34,(np.max(h)*1.2),'Comet Halley', va='top', rotation='vertical')
plt.text(13.70,(np.max(h)*1.2),'Toutatis', va='top', rotation='vertical')   
plt.text(13.00,(np.max(h)*1.2),'Comet 67P', va='top', rotation='vertical') 
plt.text(12.18,(np.max(h)*1.2),'Comet SL9', va='top', rotation='vertical')                              
plt.text(11.65,(np.max(h)*1.2),'Ryugu', va='top', rotation='vertical')
plt.text(11.15,(np.max(h)*1.2),'Bennu', va='top', rotation='vertical') 
plt.text(10.55,(np.max(h)*1.2),'Itokawa', va='top', rotation='vertical')   
plt.text(9.46,(np.max(h)*1.2),'1994 WR12', va='top', rotation='vertical')                              
plt.legend(loc='lower left', frameon = False, handletextpad=0.5, fontsize=9)
plt.savefig('PopulationMassplot.pdf')




location = 'chains/0PWDOutputs.xlsx'
workbook = xlrd.open_workbook(location)
worksheet = workbook.sheet_by_index(0)
vars()['bins']  = [0]*1600
for j in range(0,1600):
       vars()['bins'][j] = worksheet.cell(12+j,1).value
totalbins1 = vars()['bins']



for i in range(0,208):
    location = 'chains/' + str(i) +'PWDOutputs.xlsx'
    workbook = xlrd.open_workbook(location)
    worksheet = workbook.sheet_by_index(0)
    vars()['bins'+str(i)]  = [0]*1600
    for j in range(0,1600):
       vars()['bins'+str(i)][j] = worksheet.cell(12+j,2).value
for i in range(0,208):
    vars()['bins'+str(i)] = np.asarray(vars()['bins'+str(i)])
totalbins2 = [0]*1600
for i in range(0,208):
    totalbins2 = totalbins2 + vars()['bins'+str(i)]
totalbins2 = totalbins2/208 


x = totalbins1
x = np.asarray(x)
y = [0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.2,4.4,4.6,4.8,5,5.2,5.4,5.6,5.8,6,6.2,6.4,6.6,6.8,7,7.2,7.4,7.6,7.8]
y = y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y + y
y = np.asarray(y)
z = totalbins2
z = np.asarray(z)
fig = plt.figure(figsize=(7,7)) 
plt.gca().set_aspect('equal', adjustable='box')
a = np.array([0,1,2,3,4,5,6,7,8])
b = a
c = [7.28,7.28]
d = [8,7.28]
plt.plot(a,b,color='grey',linestyle='--')
plt.plot(c,d,color='grey',linestyle='--')
plt.plot([], label='Hollands et al. 2017 data', color='k', marker='s', linestyle='none', markersize=5)
plt.scatter(x,y,65,z,cmap=discrete_cmap(8, 'Greys'),marker='s')
cbar = plt.colorbar(format='%.3f', fraction = 0.045)
cbar.set_label('Probability', fontsize=14, rotation = 270, labelpad=18)
#cbar.set_ticks([0.00,0.025,0.05,0.1,0.15,0.2,0.25])
cbar.set_ticklabels([0.00,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.011,0.012,0.013,0.014,0.015])
plt.xlim(0,7.8)
plt.ylim(0,7.8)
plt.text(4,(0.5),'Declining Phase', ha='center') 
t_Mg = 10**(6.58)
if 10**(6.58) < 0.3:
    plt.text(((7.5+np.log10(5*10**(6.58)))/2),(7.5),'Steady State', ha='center')  
    plt.text(((7.5+np.log10(5*10**(6.58)))/2),(7.1),'Phase', ha='center')     
elif 0.3 <= 10**(6.58) < 2.5: 
                               plt.text(((7.5+np.log10(5*t_Mg))/2),(7.5),'Steady State', ha='center')  
                               plt.text(((7.5+np.log10(5*t_Mg))/2),(7.1),'Phase', ha='center')  
                               plt.text((1.7),(7.5),'Build-Up', ha='center')
                               plt.text(1.7,(7.1),'Phase', ha='center') 
                               plt.arrow(1.2,7.4,-0.9,0,head_width=0.12,linewidth=1.25, fc='k')
elif 2.5 <= 10**(6.58) <=200000:
                               plt.text((np.log10(5*t_Mg)/2),(7.5),'Build-Up', ha='center')
                               plt.text((np.log10(5*t_Mg)/2),(7.1),'Phase', ha='center')
                               plt.text(((7.5+np.log10(5*t_Mg))/2),(7.5),'Steady State', ha='center')  
                               plt.text(((7.5+np.log10(5*t_Mg))/2),(7.1),'Phase', ha='center') 
elif 200000 < 10**(6.58) <= 1500000:
                               plt.text((np.log10(5*t_Mg)/2),(7.5),'Build-Up', ha='center')
                               plt.text((np.log10(5*t_Mg)/2),(7.1),'Phase', ha='center') 
                               plt.text((7),(5.8),'Steady State', ha='center')
                               plt.text(7,(5.4),'Phase', ha='center') 
                               plt.arrow(7,6.2,0,1,head_width=0.12,linewidth=1.25, fc='k')   
else:
                               plt.text((np.log10(5*t_Mg)/2),(7.5),'Build-Up', ha='center')
                               plt.text((np.log10(5*t_Mg)/2),(7.1),'Phase', ha='center') 
                               plt.text((7),(5.8),'Steady State', ha='center')
                               plt.text(7,(5.4),'Phase', ha='center') 
                               plt.arrow(7.55,6.0,0,1.5,head_width=0.12,linewidth=1.25, fc='k',zorder=3)                                                    
plt.legend(loc='lower right', frameon = False, handletextpad=0.5, fontsize=9)
plt.xlabel('log(Time since Accretion Started/Yrs)', fontsize=14)
plt.ylabel('log(Accretion Event Lifetime/Yrs)', fontsize=14) 

plt.savefig('PopulationTimesinceplot.pdf')


for i in range(0,40):
 vars()['t_acc'+str(i)]  = np.sum(z[(0+(40*i)):(40+(40*i))])
totalbins = [t_acc0,t_acc1,t_acc2,t_acc3,t_acc4,t_acc5,t_acc6,t_acc7,t_acc8,t_acc9,t_acc10,t_acc11,t_acc12,t_acc13,t_acc14,t_acc15,t_acc16,t_acc17,t_acc18,t_acc19,t_acc20,t_acc21,t_acc22,t_acc23,t_acc24,t_acc25,t_acc26,t_acc27,t_acc28,t_acc29,t_acc30,t_acc31,t_acc32,t_acc33,t_acc34,t_acc35,t_acc36,t_acc37,t_acc38,t_acc39]
xbar = [0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.2,4.4,4.6,4.8,5,5.2,5.4,5.6,5.8,6,6.2,6.4,6.6,6.8,7,7.2,7.4,7.6,7.8]
plt.figure(figsize=(6,4))
plt.bar(xbar,totalbins,width=0.2,color='k',edgecolor='k',label='Hollands et al. 2017 data')
plt.axvline(x=5.484,linestyle='dotted',color='grey',lw=2) 
plt.axvline(x=6.144,linestyle='dotted',color='grey',lw=2)   
plt.axvline(x=6.794,linestyle='dotted',color='grey',lw=2) 
plt.ylim(0,(0.20))
plt.xlim(0,8)
plt.text(2.75,(0.17),'log(Accretion Event Lifetime/Yrs) $ = 6.14 \pm 0.65$', ha='center') 
plt.xlabel('log(Accretion Event Lifetime/Yrs)', fontsize=14)
plt.ylabel('Probability', fontsize=14)
plt.legend(frameon = False, loc=0)
plt.savefig('PopulationDiscLifeplot.pdf')


