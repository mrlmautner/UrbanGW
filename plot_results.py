# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 14:02:08 2018

@author: MM
"""
import flopy
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import pickle
import flopy.utils.binaryfile as bf

sns.set(style="white", palette="muted", color_codes=True)

xll = 455000
yll = 2107000
xur = 539000
yur = 2175000
cellsize = 500

STRT_YEAR = 1984
END_YEAR = 2014

ncol = int((xur-xll)/cellsize) # Number of rows
nrow = int((yur-yll)/cellsize) # Number of columns

#%% Head Dictionary
def get_heads(scenario_list):
    S_heads = {}
    for s_name in scenario_list:
        S_heads[s_name] = bf.HeadFile('model_output\VM_'+s_name+'.hds')
        
    return S_heads
    
#%% Heads Contour
def plt_head_change(s_heads, GEO, ACTIVE_LYR1, ACTIVE_LYR2, n=(30,0), m=(30,359), g_units = [2,3,4], lyr = 1):
    '''
    Calculates the raster of the change in head from the nth time step to the
    mth time step for each model in S_heads. Then plots the difference between
    the change in head over time in the first model in S_heads and all the
    other models in S_heads. Plots the heads for the geologic units g_units and
    in the active area .
    '''
    fig, axes = plt.subplots(2, 2, figsize=(7,6.3))
    mapTitle = ['Historical','Increased WW Reuse','Repair Leaks','Recharge Basins']
    plt.set_cmap('rainbow_r')
    axes = axes.flat
    cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
    
    # Get head values at the nth and mth time steps
    hds = s_heads[s_heads.keys()[0]]
    hist_i = hds.get_data(mflay=lyr,kstpkper=n)
    hist_e = hds.get_data(mflay=lyr,kstpkper=m)
    
    # Create array of heads only for geologic units g_units
    new_hist_i = np.ones(hist_i.shape)*np.nan
    for g in g_units:
        new_hist_i[GEO==g] = hist_i[GEO==g]
    new_hist_i[ACTIVE_LYR2!=lyr] = np.nan # Set cells outside model area to NaN
    
    new_hist_e = np.ones(hist_e.shape)*np.nan
    for g in g_units:
        new_hist_e[GEO==g] = hist_e[GEO==g]
    new_hist_e[ACTIVE_LYR2!=lyr] = np.nan # Set cells outside model area to NaN
    
    # Find difference between nth and mth time steps
    hist_change = new_hist_e-new_hist_i
    
    for i,s_name in enumerate(scenario_list):
        hds = s_heads[s_name]
        h_i = hds.get_data(mflay=1,kstpkper=n)
        h_e = hds.get_data(mflay=1,kstpkper=m)
        
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
        
        im = axes[i].imshow(hist_compare,vmin=0,vmax=20)
        CS = axes[i].contour(ACTIVE_LYR1, colors='k', linewidths=2)
        axes[i].xaxis.set_visible(False)
        axes[i].yaxis.set_visible(False)
        axes[i].set_title(mapTitle[i].format(i+1))
        
        fig.subplots_adjust(right=0.8)
        fig.colorbar(im, cax=cbar_ax, label='Change in Groundwater Head (m)')
        plt.savefig('model_output\plots\Hist_change_all.svg')
        plt.savefig('model_output\plots\Hist_change_all.png', dpi=600)
        plt.show()

##%% Budget
#df_1Bdget = {}
#df_extra = {}
#
#for s_name in scenario_list:
#    mf_list = flopy.utils.MfListBudget('model_output\VM_'+s_name+".list")
#    incremental, cumulative = mf_list.get_budget()
#
#    df_1Bdget[s_name], df_extra[s_name] = mf_list.get_dataframes(start_datetime="04-30-1984")
#    
##    mthly_Bdget = df_1Bdget[s_name].drop(['CONSTANT_HEAD_IN','TOTAL_IN','CONSTANT_HEAD_OUT','RECHARGE_OUT','TOTAL_OUT','IN-OUT','PERCENT_DISCREPANCY'], axis=1)
##    
##    mthly_Bdget['STORAGE_OUT'] = mthly_Bdget['STORAGE_OUT'].apply(lambda x: x*-1)
##    mthly_Bdget['WELLS_OUT'] = mthly_Bdget['WELLS_OUT'].apply(lambda x: x*-1)
##    mthly_Bdget = mthly_Bdget.multiply(30/1000000)
##    cols = mthly_Bdget.columns.tolist()
##    # reorder columns
##    cols = [cols[1]] + [cols[2]] + [cols[0]] + [cols[4]] + [cols[3]] 
##    # "commit" the reordering
##    mthly_Bdget = mthly_Bdget[cols]
#    
##    ax = mthly_Bdget['01-31-1985':'12-31-2013'].plot.area(stacked=True,figsize=(8,9),color=['blue','xkcd:azure','lightblue','red','lightpink'])
##    plt.ylabel(r'Volume ($hm^3$)')
##    plt.title('Groundwater Budget')
##    plt.legend(['Leaks','Precipitation','Storage: In','Pumping','Storage: Out'],loc=4)
##    
##    plt.savefig('model_output\plots\WB_'+modelname+'.png')
##    plt.savefig('model_output\plots\WB_'+modelname+'.svg')
##    plt.close()
#
#
##%% Cumulative overdraft
#i = 0
#
#for s_name in scenario_list:
#
#    df_extra[s_name]['IN'] = df_extra[s_name]['RECHARGE_IN'].divide(1000000) + df_extra[s_name]['WELLS_IN'].divide(1000000)
#    df_extra[s_name]['OUT'] = df_extra[s_name]['WELLS_OUT'].divide(1000000)
#    df_extra[s_name]['INOUTCUMSUM'] = df_extra[s_name]['IN'] - df_extra[s_name]['OUT']
#    
#    df_extra[s_name].INOUTCUMSUM['01-31-1985':'12-31-2013'].plot(linewidth=l[i],color=c[i],style=mark[i])
#    i+=1
#    
#plt.ylabel(r'Volume ($hm^3$)')
##plt.ylim(-8,10)
##plt.title('Cumulative In - Out')
#plt.legend(['Historical','Increase WW Reuse','Repair Leaks','Recharge Basins'])
#
#plt.savefig('model_output\plots\INOUTCumSum.png', dpi=600)
#plt.savefig('model_output\plots\INOUTCumSum.svg')
#plt.show()
#
##%% Time series by location
#t = pd.DatetimeIndex(freq='M',start='01/30/1985',end='12/31/2013')
#coords = [52,61] # subs, [29,77] # pump, [90,25] # mtn, [59,54] #
#
#hTimeS = np.zeros((348,4))
#
#z = [0,1,1,1]
#
#for i,s_name in enumerate(scenario_list):
#    j = 0
#    for y in range(1,30):
#        for m in range(0,12):
#            h = S_heads[s_name].get_data(kstpkper=((calendar.monthrange(1984+y,m+1)[1]-1),y*12+m),mflay=1)#S_heads[s_name].get_data(kstpkper=((calendar.monthrange(1984+y,m+1)[1]-1),y*12+m),mflay=1)
#            hTimeS[j,i] = h[coords[0],coords[1]]
#            j+=1
#    
#    plt.plot(t, hTimeS[:,i],linewidth=l[i],color=c[i])#,mark[i],zorder=z[i],markersize=5)
#    
#plt.xlabel('Year')
#plt.ylabel('Head Elevation')
#plt.title('Head elevation over time at pumping location')
#plt.legend(['Historical','Increased WW Reuse','Repair Leaks','Recharge Basins'])
#plt.savefig('model_output\plots\HDS_Indicator_Pump.svg')
#plt.close()
#
##%%
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#
##for s_name in ['Historical','WWTP','Leak','Basin']:
#hds = bf.HeadFile('model_output\VM_Test.hds')
#h = hds.get_data(mflay=[0,1],kstpkper=(30,359))
#
#HNEW = np.ones(h.shape)*2100#min(h[h>0])
#for i, head in enumerate(h):
#    HNEW[i][head>-900] = head[head>-900]
#
## plot a 3D wireframe like in the example mplot3d/wire3d_demo
#Z = HNEW[1]
#x = np.arange(0,ncol*cellsize,cellsize)
#y = np.arange(0,nrow*cellsize,cellsize)
#X, Y = np.meshgrid(x, y)
#
## Plot the surface.
#surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm)
#
## Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.5, aspect=5, label='Elevation (masl)')
#plt.xlabel('Easting')
#plt.ylabel('Northing')
#plt.show()
#
##%% Calculate Objectives
#WEL_INFO = {}
#
#n_scenario = 4
#n_obj = 3
#energy_array = np.zeros((n_scenario,2))
#subs_array = np.zeros((n_scenario,2))
#mound_array = np.zeros((n_scenario,3))
#
#for s, s_name in enumerate(scenario_list):
#    
#    with open('model_output\objective_data\WEL_INFO_'+s_name+'.pickle', 'rb') as handle:
#        WEL_INFO = pickle.load(handle)
#    with open('model_output\objective_data\LU_'+s_name+'.pickle', 'rb') as handle:
#        LU = pickle.load(handle)
#    heads = S_heads[s_name]
#    
#    energy_array[s,:] = mo.measureEnergy(heads,WEL_INFO,DEM)
#    subs_array[s,:] = mo.measureSubidence(heads,DEM,ACTIVE_LYR1,TH1)
#    mound_array[s,:] = mo.measureMound(heads,DEM,ACTIVE_LYR1,LU,[132,252])
#
##%%
#    
#energy = energy_array[:,0]
#subs = subs_array[:,0]/subs_array[:,1]
#mound = mound_array[:,0]/min(mound_array[:,0])
#cells = (360-14)*np.sum(np.sum(ACTIVE_LYR2)) # Number of cells in Model Area over periods 15-360
#

# Plotting defaults
l = [4,2,2,2]
c = ['k','goldenrod','blue','darkgreen']
mark = ['-','-','-','--']

scenario_list = ['Historical','WWTP','Leak','Basin']

#barWidth = 0.5
#r = np.arange(n_scenario)*0.5 # bar position
#y_label = ['Kilowatt Hours','Average Head Below top of Clay Layer','Ratio to Historical Mounding']
#obj_title = ['Energy Use','Subsidence Avoidance','Mounding']
#ylims = [[4.5E9,6E9],[36,41],[0.9,1.25]]
#
#normalized_o = np.zeros((n_scenario,n_obj))
#
#fig, axes = plt.subplots(nrows=1, ncols=3,figsize=(10,6))
#
#for o, obj in enumerate([energy,subs,mound]):
#    normalized_o[:,o] = obj / (obj.max(axis=0) - obj.min(axis=0))
#    plt.subplot(1,n_obj,o+1)
#    plt.bar(r, obj, width=barWidth, edgecolor='white', color=c)
#    plt.xticks(r, ['Historical','WWTP','Leak','Basin'],rotation='vertical')
#    plt.ylabel(y_label[o])
#    plt.ylim(ylims[o])
#    plt.title(obj_title[o])
#
#fig.tight_layout()
## Flip subsidence measure to be minimizing
##normalized_o[:,1] = 1 - normalized_o[:,1]
#plt.savefig('model_output\plots\Objectives.svg')
#plt.savefig('model_output\plots\Objectives.png', dpi=600)
#plt.show()

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

def parallel_axis(nondom_results,obj_labels,filename):
    # Plots a normalized parallel axis
    plt.figure()
    for ppoint in nondom_results:
        ppoint = (ppoint - nondom_results.min(axis=0)) / (nondom_results.max(axis=0) - nondom_results.min(axis=0))
        plt.plot(range(len(obj_labels)), ppoint, 'steelblue')
    
    plt.gca().set_xticks(range(len(obj_labels)))
    plt.gca().set_xticklabels(obj_labels)
    plt.show()
    plt.savefig(filename)
    