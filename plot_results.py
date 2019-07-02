# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 14:02:08 2018

@author: MM
"""
import os
os.chdir(r'D:\MMautner\UrbanGW')

import flopy
import numpy as np
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.dates as mdates
import pandas as pd
import calendar
import pickle
from gwscripts.dataprocessing import gengriddata as gen
from gwscripts.optimization import measureobjectives as mo
import seaborn as sns
import datetime

sns.set(style="white", palette="muted", color_codes=True)
plt.rcParams['legend.fontsize'] = 20
plt.rcParams['axes.titlesize'] = 22
plt.rcParams['axes.labelsize'] = 22
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['figure.titlesize'] = 24
plt.rcParams.update({'font.size': 20})

#%%

xll = 455000
yll = 2107000
xur = 539000
yur = 2175000
cellsize = 500

STRT_YEAR = 1984
END_YEAR = 2014

ncol = int((xur-xll)/cellsize) # Number of rows
nrow = int((yur-yll)/cellsize) # Number of columns

# Load datasets
ACTIVE_LYR1 = gen.openASC('data_processed\ACTIVE_VM_LYR1.asc')
ACTIVE_LYR2 = gen.openASC('data_processed\ACTIVE_VM_LYR2.asc')
TH1 = gen.openASC('data_processed\THICK1_VM.asc')
DEM = gen.openASC('data_processed\DEM_VM.asc')
GEO = gen.openASC('data_processed\GEO_VM.asc')

# Plotting defaults
l = [4,2,2,2]
c = ['k','goldenrod','blue','darkgreen']
mark = ['-','-','-','--']

scenario_list = ['Historical','WWTP','Leak','Basin']
#%% Head Dictionary
S_heads = {}
for s_name in scenario_list:
    S_heads[s_name] = bf.HeadFile('model_output\VM_'+s_name+'.hds')
    
#%% Heads Contour
fig, axes = plt.subplots(2, 2, figsize=(7,6.3))
a=0
mapTitle = ['Historical','Increased WW Reuse','Repair Leaks','Recharge Basins']
plt.set_cmap('rainbow_r')
axes = axes.flat
cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])

hds = S_heads['Historical']
hist_i = hds.get_data(mflay=1,kstpkper=(30,0))
hist_e = hds.get_data(mflay=1,kstpkper=(30,359))

new_hist_i = np.ones(hist_i.shape)*np.nan#min(h[h>0])
new_hist_i[GEO==2] = hist_i[GEO==2]
#new_hist_i[GEO==3] = hist_i[GEO==3]
new_hist_i[GEO==4] = hist_i[GEO==4]
#new_hist_i[GEO==5] = hist_i[GEO==5]
new_hist_i[ACTIVE_LYR2!=1] = np.nan
new_hist_e = np.ones(hist_e.shape)*np.nan#min(h[h>0])
new_hist_e[GEO==2] = hist_e[GEO==2]
#new_hist_e[GEO==3] = hist_e[GEO==3]
new_hist_e[GEO==4] = hist_e[GEO==4]
#new_hist_e[GEO==5] = hist_e[GEO==5]
new_hist_e[ACTIVE_LYR2!=1] = np.nan

hist_change = new_hist_e-new_hist_i

for i,s_name in enumerate(scenario_list):
    hds = S_heads[s_name]
    h_i = hds.get_data(mflay=1,kstpkper=(30,0))
    h_e = hds.get_data(mflay=1,kstpkper=(30,359))
    
    new_h_i = np.ones(h_i.shape)*np.nan#min(h[h>0])
    new_h_i[GEO==2] = h_i[GEO==2]
#    new_h_i[GEO==3] = h_i[GEO==3]
    new_h_i[GEO==4] = h_i[GEO==4]
#    new_h_i[GEO==5] = h_i[GEO==5]
    new_h_i[ACTIVE_LYR2!=1] = np.nan
    new_h_e = np.ones(h_e.shape)*np.nan#min(h[h>0])
    new_h_e[GEO==2] = h_e[GEO==2]
#    new_h_e[GEO==3] = h_e[GEO==3]
    new_h_e[GEO==4] = h_e[GEO==4]
#    new_h_e[GEO==5] = h_e[GEO==5]
    new_h_e[ACTIVE_LYR2!=1] = np.nan
    
    h_change = new_h_e-new_h_i
    hist_compare = h_change - hist_change
    
    im = axes[a].imshow(hist_compare,vmin=0,vmax=20)
    CS = axes[a].contour(ACTIVE_LYR1, colors='k', linewidths=2)
    axes[a].xaxis.set_visible(False)
    axes[a].yaxis.set_visible(False)
    axes[a].set_title(mapTitle[a].format(i+1))
    
    a+=1
    
fig.subplots_adjust(right=0.8)
fig.colorbar(im, cax=cbar_ax, label='Change in Groundwater Head (m)')
plt.savefig('model_files\output\plots\Hist_change_all.svg')
plt.savefig('model_files\output\plots\Hist_change_all.png', dpi=600)
plt.show()

#%% Multi-Heads Countour

print('Plotting head changes over model period...')
plt.set_cmap('viridis_r')
mapTitles = ['Historical','Increased WW Reuse','Repair Leaks','Recharge Basins']

# Get head values at the nth and mth time steps
first_heads = S_heads['Historical']

def calc_head_change(hds, GEO, ACTIVE, n, m, g_units, lyr):
    # Get head values at the nth and mth time steps
    h_i = hds.get_data(mflay=lyr,kstpkper=n)
    h_e = hds.get_data(mflay=lyr,kstpkper=m)
    
    # Create array of heads only for geologic units g_units
    new_h_i = np.ones(h_i.shape)*np.nan
    for g in g_units:
        new_h_i[GEO==g] = h_i[GEO==g]
    new_h_i[ACTIVE!=1] = np.nan # Set cells outside model area to NaN
    
    new_h_e = np.ones(h_e.shape)*np.nan
    for g in g_units:
        new_h_e[GEO==g] = h_e[GEO==g]
    new_h_e[ACTIVE!=1] = np.nan # Set cells outside model area to NaN
    
    # Find difference between nth and mth time steps
    head_change = new_h_e - new_h_i
    return head_change
 
first_change = calc_head_change(first_heads, GEO, ACTIVE_LYR2,  n = (30,0), m = (30,359), g_units = [2,3,4,5], lyr = 1)

for i, s_name in enumerate(scenario_list):
    scen_heads = S_heads[s_name]
    
    scen_change = calc_head_change(scen_heads, GEO, ACTIVE_LYR2,  n = (30,0), m = (30,359), g_units = [2,3,4,5], lyr = 1)
    
    change_compare = scen_change - first_change
    
    fig, axes = plt.subplots()
    cbar_ax = fig.add_axes()
    
    im = axes.imshow(change_compare,vmin=0,vmax=20)
    CS = axes.contour(ACTIVE_LYR2, colors='k', linewidths=2)
    
    axes.xaxis.set_visible(False)
    axes.yaxis.set_visible(False)
    axes.set_title(mapTitles[i].format(i+1),fontweight="bold", size=20)
    
    fig.colorbar(im, cax=cbar_ax, label='Change in Groundwater Head (m)')
#        plt.savefig('model_output\plots\head-change_'+s_name+'-'+scenario_list[0]+'.svg')
    plt.savefig('model_files\output\plots\head-change_'+s_name+'-'+scenario_list[0]+'.png', dpi=600)
    

#%% Budget
df_1Bdget = {}
df_extra = {}

for s_name in scenario_list:
    mf_list = flopy.utils.MfListBudget('model_output\VM_'+s_name+".list")
    incremental, cumulative = mf_list.get_budget()

    df_1Bdget[s_name], df_extra[s_name] = mf_list.get_dataframes(start_datetime="04-30-1984")
    
#    mthly_Bdget = df_1Bdget[s_name].drop(['CONSTANT_HEAD_IN','TOTAL_IN','CONSTANT_HEAD_OUT','RECHARGE_OUT','TOTAL_OUT','IN-OUT','PERCENT_DISCREPANCY'], axis=1)
#    
#    mthly_Bdget['STORAGE_OUT'] = mthly_Bdget['STORAGE_OUT'].apply(lambda x: x*-1)
#    mthly_Bdget['WELLS_OUT'] = mthly_Bdget['WELLS_OUT'].apply(lambda x: x*-1)
#    mthly_Bdget = mthly_Bdget.multiply(30/1000000)
#    cols = mthly_Bdget.columns.tolist()
#    # reorder columns
#    cols = [cols[1]] + [cols[2]] + [cols[0]] + [cols[4]] + [cols[3]] 
#    # "commit" the reordering
#    mthly_Bdget = mthly_Bdget[cols]
    
#    ax = mthly_Bdget['01-31-1985':'12-31-2013'].plot.area(stacked=True,figsize=(8,9),color=['blue','xkcd:azure','lightblue','red','lightpink'])
#    plt.ylabel(r'Volume ($hm^3$)')
#    plt.title('Groundwater Budget')
#    plt.legend(['Leaks','Precipitation','Storage: In','Pumping','Storage: Out'],loc=4)
#    
#    plt.savefig('model_output\plots\WB_'+modelname+'.png')
#    plt.savefig('model_output\plots\WB_'+modelname+'.svg')
#    plt.close()


#%% Cumulative overdraft
plt.figure(figsize=(12,7.2))

i = 0
for s_name in scenario_list:

    df_extra[s_name]['IN'] = df_extra[s_name]['RECHARGE_IN'].divide(1000000) + df_extra[s_name]['WELLS_IN'].divide(1000000)
    df_extra[s_name]['OUT'] = df_extra[s_name]['WELLS_OUT'].divide(1000000)
    df_extra[s_name]['INOUTCUMSUM'] = df_extra[s_name]['IN'] - df_extra[s_name]['OUT']
    
    df_extra[s_name].INOUTCUMSUM['01-31-1985':'12-31-2013'].plot(linewidth=l[i],color=c[i],style=mark[i])
    i+=1


plt.ylabel(r'Volume (million m$^3$)')
#plt.xlabel(r'Year')
#plt.ylim(-8,10)
#plt.title('Cumulative In - Out')
plt.legend(['Historical','Increase WW Reuse','Repair Leaks','Recharge Basins'],loc='upper right')
plt.gcf().subplots_adjust(left=0.15,right=.95,bottom=0.1,top=.95)

plt.savefig('model_files\output\plots\INOUTCumSum_0607.png', dpi=600)
plt.savefig('model_files\output\plots\INOUTCumSum_0607.svg')
plt.show()

#%% Time series by location
t = pd.DatetimeIndex(freq='M',start='01/30/1985',end='12/31/2013')
coords = [52,61] # subs, [29,77] # pump, [90,25] # mtn, [59,54] #

hTimeS = np.zeros((348,4))

z = [0,1,1,1]

for i,s_name in enumerate(scenario_list):
    j = 0
    for y in range(1,30):
        for m in range(0,12):
            h = S_heads[s_name].get_data(kstpkper=((calendar.monthrange(1984+y,m+1)[1]-1),y*12+m),mflay=1)#S_heads[s_name].get_data(kstpkper=((calendar.monthrange(1984+y,m+1)[1]-1),y*12+m),mflay=1)
            hTimeS[j,i] = h[coords[0],coords[1]]
            j+=1
    
    plt.plot(t, hTimeS[:,i],linewidth=l[i],color=c[i])#,mark[i],zorder=z[i],markersize=5)
    
plt.xlabel('Year')
plt.ylabel('Head Elevation')
plt.title('Head elevation over time at pumping location')
plt.legend(['Historical','Increased WW Reuse','Repair Leaks','Recharge Basins'])
plt.savefig('model_files\output\plots\HDS_Indicator_Pump.svg')
plt.close()

#%%
fig = plt.figure()
ax = fig.gca(projection='3d')

#for s_name in ['Historical','WWTP','Leak','Basin']:
hds = bf.HeadFile('model_output\VM_Test.hds')
h = hds.get_data(mflay=[0,1],kstpkper=(30,359))

HNEW = np.ones(h.shape)*2100#min(h[h>0])
for i, head in enumerate(h):
    HNEW[i][head>-900] = head[head>-900]

# plot a 3D wireframe like in the example mplot3d/wire3d_demo
Z = HNEW[1]
x = np.arange(0,ncol*cellsize,cellsize)
y = np.arange(0,nrow*cellsize,cellsize)
X, Y = np.meshgrid(x, y)

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm)

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5, label='Elevation (masl)')
plt.xlabel('Easting')
plt.ylabel('Northing')
plt.show()

#%% Calculate Objectives
WEL_INFO = {}

n_scenario = 4
n_obj = 3
energy_array = np.zeros((n_scenario,2))
subs_array = np.zeros((n_scenario,2))
mound_array = np.zeros((n_scenario,3))

for s, s_name in enumerate(scenario_list):
    
    with open('model_files\output\objective_data\WEL_INFO_'+s_name+'.pickle', 'rb') as handle:
        WEL_INFO = pickle.load(handle)
    with open('model_files\output\objective_data\LU_'+s_name+'.pickle', 'rb') as handle:
        LU = pickle.load(handle)
    heads = S_heads[s_name]
    
    energy_array[s,:] = mo.measureEnergy(heads,WEL_INFO,DEM)
    subs_array[s,:] = mo.measureSubidence(heads,DEM,ACTIVE_LYR1,TH1)
    mound_array[s,:] = mo.measureMound(heads,DEM,ACTIVE_LYR1,LU,[132,252])

#%%
energy = energy_array[:,0]
subs = subs_array[:,0]/subs_array[:,1]
mound = mound_array[:,0]/min(mound_array[:,0])
cells = (360-14)*np.sum(np.sum(ACTIVE_LYR2)) # Number of cells in Model Area over periods 15-360

barWidth = 0.5
r = np.arange(n_scenario)*0.5 # bar position
y_label = ['Cumulative Pumping Energy (kWh)','Average Depth to GW in Clay Layer (m)','Normalized GW Mounding in Urban Areas']
obj_title = ['Energy Use','Subsidence Avoidance','Urban Mounding']
#ylims = [[4.5E9,6E9],[36,41],[0.9,1.25]]

normalized_o = np.zeros((n_scenario,n_obj))

fig, axes = plt.subplots(nrows=1, ncols=3,figsize=(12,7.2))

for o, obj in enumerate([energy,subs,mound]):
    normalized_o[:,o] = obj / (obj.max(axis=0) - obj.min(axis=0))
    plt.subplot(1,n_obj,o+1)
    plt.bar(r, obj, width=barWidth, edgecolor='white', color=c)
    plt.xticks(r, ['Historical','WWTP','Leak','Basin'],rotation=35,ha='right')
    plt.ylabel(y_label[o])
#    plt.ylim(ylims[o])
    plt.title(obj_title[o],fontweight='bold', y=1.04)

#fig.tight_layout()
# Flip subsidence measure to be minimizing
#normalized_o[:,1] = 1 - normalized_o[:,1]
plt.gcf().subplots_adjust(wspace=0.45,left=0.09,right=.97,bottom=0.15,top=.9)
plt.savefig('model_files\output\plots\Objectives_190607.svg')
plt.savefig('model_files\output\plots\Objectives_190607.png', dpi=600)
plt.show()

##%%
#modelname = 'model_output\VM_WWTP'
#heads = bf.HeadFile(modelname+'.hds')
#
#h = heads.get_data(mflay=[0,1],kstpkper=(30,24))
#plt.set_cmap('coolwarm_r')
#
#HNEW = np.ones(h.shape)*np.nan#min(h[h>0])
#for i, head in enumerate(h):
#    HNEW[i][head>-900] = head[head>-900]
#
#for i, hdslayer in enumerate(HNEW):
#    im = plt.imshow(hdslayer, vmin=2100)
#    ctr = plt.contour(hdslayer, colors='k', linewidths=0.5)
#
#plt.colorbar(im, label='Groundwater Head (m)')

#%% Head Obsevations
obsformation = pd.read_csv('data_raw\obs_formation.csv')
df = pd.read_fwf('model_files\output\VM_Historical-IH2750.hob.out')
df.columns = ['simulated','observed','obs_name','nan']
df = df.drop('nan',axis=1)
df['time_series'] = np.nan
df['obs_id'] = np.nan
df['dem'] = np.nan
df['geo'] = np.nan
df['ddsimulated'] = np.nan
df['ddobserved'] = np.nan

obsinfo = hobs.obs_data
for i in obsinfo:
    oname = i.obsname
    
    t = np.ones(len(df[df['obs_name'].str.contains(oname)].index))*np.nan
    t[:i.nobs] = i.time_series_data.astype(int).view(int)

    df.loc[df['obs_name'].str.contains(oname),'time_series'] = t
    df.loc[df['obs_name'].str.contains(oname),'obs_id'] = oname
    
    df.loc[df['obs_name'].str.contains(oname),'ddsimulated'] = df[df['obs_name'].str.contains(oname)]['simulated'] - df[df['obs_name'].str.contains(oname)]['simulated'].values[0]
    df.loc[df['obs_name'].str.contains(oname),'ddobserved'] = df[df['obs_name'].str.contains(oname)]['observed'] - df[df['obs_name'].str.contains(oname)]['observed'].values[0]

geology = ['Lacustrine','Alluvial','Basalt','Volcaniclastic','Andesite']
for i, r in obsformation.iterrows():
    df.loc[df['obs_id']==r['IDPOZO'],'dem'] = r['DEM_c3_INE']
    df.loc[df['obs_id']==r['IDPOZO'],'geo'] = geology[r['GEOLOGY_ZO']-1]

#%%
ax = sns.scatterplot(x='simulated', y='observed', hue='geo',data=df)
ax.set_ylim(2100,2700)
ax.set_xlim(2100,2700)

#%%
sns.set_style("whitegrid")
g = sns.lmplot(x='simulated', y='observed', hue='geo',data=df,legend=False,palette=dict(Alluvial=(1,0.867,0), Basalt=(0,0.788,0.498) ,Volcaniclastic=(0.9, 0, 0.455) ),size=5,aspect=1.2,scatter_kws={'edgecolor':None,'s':10})
plt.plot(np.linspace(2000,2800,1000), np.linspace(2000,2800,1000), 'k',linestyle=':')
plt.legend(['Tarango (Volcaniclastic)','Alluvial','Fractured Basalt','One-to-one'],loc='best')
plt.xlim(2100,2650)
plt.ylim(2100,2650)
plt.ylabel('Observed Head (masl)')
plt.xlabel('Simulated Head (masl)')
plt.savefig('C:\Users\MM\Google Drive\Davis\Research\Papers\Multi-Objective_Spatial_ValleyMexico\Figures\Sim_vs_Obs-IH2750.jpeg', dpi=600)
plt.close()
#%%
sns.set_style("whitegrid")
g = sns.lmplot(x='ddsimulated', y='ddobserved', hue='geo',data=df,legend=False,palette=dict(Alluvial=(1,0.867,0), Basalt=(0,0.788,0.498) ,Volcaniclastic=(0.9, 0, 0.455) ),size=5,aspect=1.2,scatter_kws={'edgecolor':None,'s':10})
plt.plot(np.linspace(-100,100,1000), np.linspace(-100,100,1000), 'k',linestyle=':')
plt.legend(['Tarango (Volcaniclastic)','Alluvial','Fractured Basalt','One-to-one'],loc='best')
plt.xlim(-60,60)
plt.ylim(-60,60)
plt.ylabel('Observed Head (masl)')
plt.xlabel('Simulated Head (masl)')
plt.savefig('C:\Users\MM\Google Drive\Davis\Research\Papers\Multi-Objective_Spatial_ValleyMexico\Figures\Sim_vs_Obs-ddn-IH2750.jpeg', dpi=600)
plt.close()
#%%

for o in obsformation['IDPOZO'].values:
    ddn_data = df[df['obs_id']==o].copy()
    ddn_data['time_series'] = pd.to_timedelta(ddn_data['time_series'],'d') + pd.to_datetime('1984-01-01')
    ddn_data = ddn_data.set_index('time_series')
    fig, axes = plt.subplots(2, 1, figsize=(5,12))
    ddn_data['simulated'].plot(ax=axes[0])
    ddn_data['observed'].plot(ax=axes[0])
    axes[0].legend(['Simulated','Observed'])
    axes[0].xaxis.label.set_visible(False)
    
    ddn_data['ddsimulated'].plot(ax=axes[1])
    ddn_data['ddobserved'].plot(ax=axes[1])
    axes[1].set_ylim([-35, 10])
    axes[1].legend(['Simulated','Observed'])
    axes[1].xaxis.label.set_visible(False)
    
    plt.savefig('model_files\output\plots\observations\IH_2750\drawdown_'+o+'.jpeg', dpi=300)
    plt.close()
    