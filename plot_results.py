# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 14:02:08 2018

@author: MM
"""
import flopy
import numpy as np
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf
from matplotlib import cm
import pandas as pd
import calendar
import seaborn as sns
import pickle
from scipy import stats

sns.set(style="white", palette="muted", color_codes=True)
plt.rcParams['legend.fontsize'] = 20
plt.rcParams['axes.titlesize'] = 22
plt.rcParams['axes.labelsize'] = 22
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['figure.titlesize'] = 24
plt.rcParams.update({'font.size': 20})

def get_heads(scenario_list):
    S_heads = {}
    for s_name in scenario_list:
        S_heads[s_name] = bf.HeadFile('model_files\modflow\VM_'+s_name+'.hds')

    return S_heads

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

def plt_head_change(scenario_list, mapTitles, s_heads, GEO, ACTIVE, n = (8,23), m = (8,359), g_units = [2,3,4,5], lyr = 1):
    '''
    Calculates the raster of the change in head from the nth time step to the
    mth time step for each model in S_heads. Then plots the difference between
    the change in head over time in the first model in S_heads and all the
    other models in S_heads. Plots the heads for the geologic units g_units and
    in the active area .
    '''
    print('Plotting head changes over model period...')
    plt.set_cmap('viridis_r')

    # Get head values at the nth and mth time steps
    first_heads = s_heads[list(s_heads.keys())[0]]

    first_change = calc_head_change(first_heads, GEO, ACTIVE, n, m, g_units, lyr)

    for i, s_name in enumerate(scenario_list):
        scen_heads = s_heads[s_name]

        scen_change = calc_head_change(scen_heads, GEO, ACTIVE, n, m, g_units, lyr)

        change_compare = scen_change - first_change

        fig, axes = plt.subplots()
        cbar_ax = fig.add_axes()

        im = axes.imshow(change_compare,vmin=0,vmax=20)
        CS = axes.contour(ACTIVE, colors='k', linewidths=2)

        axes.xaxis.set_visible(False)
        axes.yaxis.set_visible(False)
        axes.set_title(mapTitles[i].format(i+1))

        fig.colorbar(im, cax=cbar_ax, label='Change in Groundwater Head (m)')

        plt.savefig('model_files\output\plots\head-change_'+s_name+'-'+scenario_list[0]+'.png', dpi=600)
        plt.close()

def plt_scen_objectives(scenario_names, num_scen, objectives):
    '''
    objectives is a list with an array of length number of scenarios for each objective
    '''
    print('Plotting scenario performance under objectives...')
    plt.rcParams['legend.fontsize'] = 20
    plt.rcParams['axes.titlesize'] = 22
    plt.rcParams['axes.labelsize'] = 22
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 18
    plt.rcParams['figure.titlesize'] = 24
    plt.rcParams.update({'font.size': 20})
    
    c = ['k','goldenrod','blue','darkgreen']
    barWidth = 0.1
    r = np.arange(num_scen)*0.1 # bar position
    y_label = ['Pumping Energy (kWh)','Depth to Groundwater in Clay (m)','Percent of Urban Cells Flooded']
    obj_title = ['Energy Use','Subsidence Avoidance','Urban Mounding']

    normalized_o = np.zeros((num_scen, len(objectives)))

    fig, axes = plt.subplots(nrows=1, ncols=3,figsize=(12,7.2))

    for o, obj in enumerate(objectives):
        normalized_o[:,o] = obj / (obj.max(axis=0) - obj.min(axis=0))
        plt.subplot(1,len(objectives),o+1)
        plt.bar(r, obj, width=barWidth, edgecolor='white', color=c)
        plt.xticks(r, scenario_names, rotation=35, ha='right')
        plt.ylabel(y_label[o])
        plt.title(obj_title[o], fontweight='bold', y=1.04)

    # Flip subsidence measure to be minimizing
    plt.gcf().subplots_adjust(wspace=0.45,left=0.09,right=.97,bottom=0.15,top=.9)
    plt.savefig('model_files\output\plots\Objectives.png', dpi=600)
    plt.close()

def parallel_axis(nondom_results, obj_labels, opt_run):
    # Plots a normalized parallel axis
    plt.figure()
    for ppoint in nondom_results:
        ppoint = (ppoint - nondom_results.min(axis=0)) / (nondom_results.max(axis=0) - nondom_results.min(axis=0))
        plt.plot(range(len(obj_labels)), ppoint, 'steelblue')

    plt.gca().set_xticks(range(len(obj_labels)))
    plt.gca().set_xticklabels(obj_labels)
    plt.savefig('parallelaxis_' + opt_run + '.png')
    plt.close()

def get_budgets(scenario_list, mapTitles, s_heads):
    # Budget
    df_Bdget = {}
    df_CumSum = {}
    
    for s_name in scenario_list:
        mf_list = flopy.utils.MfListBudget('model_files\modflow\VM_'+s_name+".list")
        incremental, cumulative = mf_list.get_budget()
    
        df_Bdget[s_name], df_CumSum[s_name] = mf_list.get_dataframes(start_datetime="04-30-1984")
        
        mthly_Bdget = df_Bdget[s_name].drop(['CONSTANT_HEAD_IN', 'TOTAL_IN', 'CONSTANT_HEAD_OUT', 'RECHARGE_OUT', 'TOTAL_OUT', 'IN-OUT', 'PERCENT_DISCREPANCY'], axis=1)
        
        mthly_Bdget['STORAGE_OUT'] = mthly_Bdget['STORAGE_OUT'].apply(lambda x: x*-1)
        mthly_Bdget['WELLS_OUT'] = mthly_Bdget['WELLS_OUT'].apply(lambda x: x*-1)
#        mthly_Bdget['DRAINS_OUT'] = mthly_Bdget['DRAINS_OUT'].apply(lambda x: x*-1)
        mthly_Bdget = mthly_Bdget.multiply(30/1000000)
        cols = mthly_Bdget.columns.tolist()
        # reorder columns
        cols = [cols[1]] + [cols[2]] + [cols[0]] + [cols[4]] + [cols[3]] 
        # "commit" the reordering
        mthly_Bdget = mthly_Bdget[cols]
        
    for s_name in scenario_list:
        df_CumSum[s_name]['IN'] = df_CumSum[s_name]['RECHARGE_IN'].divide(1000000) + df_CumSum[s_name]['WELLS_IN'].divide(1000000)
        df_CumSum[s_name]['OUT'] = df_CumSum[s_name]['WELLS_OUT'].divide(1000000) #+ df_CumSum[s_name]['DRAINS_OUT'].divide(1000000)
        df_CumSum[s_name]['INOUTCUMSUM'] = df_CumSum[s_name]['IN'] - df_CumSum[s_name]['OUT']
    
    return df_Bdget, mthly_Bdget, df_CumSum
    
def plt_cum_sum(filename, scenario_list, mapTitles, df_CumSum, start='01-31-1985', end='12-31-2013'):
    # Plotting defaults
    l = [4,2,2,2]
    c = ['k','goldenrod','blue','darkgreen']
    mark = ['-','-','-','--']
    
    for i,s_name in enumerate(scenario_list):
        df_CumSum[s_name].INOUTCUMSUM[start:end].plot(linewidth=l[i],color=c[i],style=mark[i])
        
    plt.ylabel(r'Volume (million m$^3$)')
    #plt.xlabel(r'Year')
    #plt.ylim(-8,10)
    #plt.title('Cumulative In - Out')
    plt.legend(mapTitles,loc='lower left')
    plt.gcf().subplots_adjust(left=0.15,right=.95,bottom=0.1,top=.95)
    
    plt.savefig(filename, dpi=600)
    plt.show()
    
##from gwscripts.dataprocessing import gengriddata as gen
#from gwscripts.optimization import measureobjectives as mo
#

#scenario_list = ['Historical','WWTP','Leak','Basin']
##%% Head Dictionary
#S_heads = {}
#for s_name in scenario_list:
#    S_heads[s_name] = bf.HeadFile('model_files\modflow\VM_'+s_name+'.hds')
#
##%% Multi-Heads Countour
#
#print('Plotting head changes over model period...')
#plt.set_cmap('viridis_r')
#mapTitles = ['Historical','Increased WW Reuse','Repair Leaks','Recharge Basins']
#
## Get head values at the nth and mth time steps
#first_heads = S_heads['Historical']
#
#def calc_head_change(hds, GEO, ACTIVE, n, m, g_units, lyr):
#    # Get head values at the nth and mth time steps
#    h_i = hds.get_data(mflay=lyr,kstpkper=n)
#    h_e = hds.get_data(mflay=lyr,kstpkper=m)
#    
#    # Create array of heads only for geologic units g_units
#    new_h_i = np.ones(h_i.shape)*np.nan
#    for g in g_units:
#        new_h_i[GEO==g] = h_i[GEO==g]
#    new_h_i[ACTIVE!=1] = np.nan # Set cells outside model area to NaN
#    
#    new_h_e = np.ones(h_e.shape)*np.nan
#    for g in g_units:
#        new_h_e[GEO==g] = h_e[GEO==g]
#    new_h_e[ACTIVE!=1] = np.nan # Set cells outside model area to NaN
#    
#    # Find difference between nth and mth time steps
#    head_change = new_h_e - new_h_i
#    return head_change
# 
#first_change = calc_head_change(first_heads, GEO, ACTIVE_LYR2,  n = (9,14), m = (9,359), g_units = [2,3,4,5], lyr = 1)
#
#for i, s_name in enumerate(scenario_list):
#    scen_heads = S_heads[s_name]
#    
#    scen_change = calc_head_change(scen_heads, GEO, ACTIVE_LYR2,  n = (9,14), m = (9,359), g_units = [2,3,4,5], lyr = 1)
#    
#    change_compare = scen_change - first_change
#    
#    fig, axes = plt.subplots()
#    cbar_ax = fig.add_axes()
#    
#    im = axes.imshow(change_compare,vmin=0,vmax=20)
#    CS = axes.contour(ACTIVE_LYR2, colors='k', linewidths=2)
#    
#    axes.xaxis.set_visible(False)
#    axes.yaxis.set_visible(False)
#    axes.set_title(mapTitles[i].format(i+1),fontweight="bold", size=20)
#    
#    fig.colorbar(im, cax=cbar_ax, label='Change in Groundwater Head (m)')
##        plt.savefig('model_files\output\plots\head-change_'+s_name+'-'+scenario_list[0]+'.svg')
#    plt.savefig('model_files\output\plots\head-change_'+s_name+'-'+scenario_list[0]+'.png', dpi=600)
#    
#
##%% Budget
#df_1Bdget = {}
#df_extra = {}
#
#for s_name in scenario_list:
#    mf_list = flopy.utils.MfListBudget('model_files\modflow\VM_'+s_name+".list")
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
##    plt.savefig('model_files\output\plots\WB_'+modelname+'.png')
##    plt.savefig('model_files\output\plots\WB_'+modelname+'.svg')
##    plt.close()
#
#
##%% Cumulative overdraft
#plt.figure(figsize=(12,7.2))
#
#i = 0
#for s_name in scenario_list:
#
#    df_extra[s_name]['IN'] = df_extra[s_name]['RECHARGE_IN'].divide(1000000) + df_extra[s_name]['WELLS_IN'].divide(1000000)
#    df_extra[s_name]['OUT'] = df_extra[s_name]['WELLS_OUT'].divide(1000000)
#    df_extra[s_name]['INOUTCUMSUM'] = df_extra[s_name]['IN'] - df_extra[s_name]['OUT']
#    
#    df_extra[s_name].INOUTCUMSUM['01-31-1985':'12-31-2013'].plot(linewidth=l[i],color=c[i],style=mark[i])
#    i+=1
#
#
#plt.ylabel(r'Volume (million m$^3$)')
##plt.xlabel(r'Year')
##plt.ylim(-8,10)
##plt.title('Cumulative In - Out')
#plt.legend(['Historical','Increase WW Reuse','Repair Leaks','Recharge Basins'],loc='upper right')
#plt.gcf().subplots_adjust(left=0.15,right=.95,bottom=0.1,top=.95)
#
#plt.savefig('model_files\output\plots\INOUTCumSum_0607.png', dpi=600)
#plt.savefig('model_files\output\plots\INOUTCumSum_0607.svg')
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
#plt.savefig('model_files\output\plots\HDS_Indicator_Pump.svg')
#plt.close()
#
##%%
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#
##for s_name in ['Historical','WWTP','Leak','Basin']:
#hds = bf.HeadFile('model_files\modflow\VM_Test.hds')
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
#    with open('model_files\output\objective_data\WEL_INFO_'+s_name+'.pickle', 'rb') as handle:
#        WEL_INFO = pickle.load(handle)
#    with open('model_files\output\objective_data\LU_'+s_name+'.pickle', 'rb') as handle:
#        LU = pickle.load(handle)
#    heads = S_heads[s_name]
#    
#    energy_array[s,:] = mo.measureEnergy(heads,WEL_INFO,DEM)
#    subs_array[s,:] = mo.measureSubidence(heads,DEM,ACTIVE_LYR1,TH1)
#    mound_array[s,:] = mo.measureMound(heads,DEM,ACTIVE_LYR1,LU,[132,252])
#
##%%
#energy = energy_array[:,0]
#subs = subs_array[:,0]/subs_array[:,1]
#mound = mound_array[:,0]/min(mound_array[:,0])
#cells = (360-14)*np.sum(np.sum(ACTIVE_LYR2)) # Number of cells in Model Area over periods 15-360
#
#barWidth = 0.5
#r = np.arange(n_scenario)*0.5 # bar position
#y_label = ['Cumulative Pumping Energy (kWh)','Average Depth to GW in Clay Layer (m)','Normalized GW Mounding in Urban Areas']
#obj_title = ['Energy Use','Subsidence Avoidance','Urban Mounding']
##ylims = [[4.5E9,6E9],[36,41],[0.9,1.25]]
#
#normalized_o = np.zeros((n_scenario,n_obj))
#
#fig, axes = plt.subplots(nrows=1, ncols=3,figsize=(12,7.2))
#
#for o, obj in enumerate([energy,subs,mound]):
#    normalized_o[:,o] = obj / (obj.max(axis=0) - obj.min(axis=0))
#    plt.subplot(1,n_obj,o+1)
#    plt.bar(r, obj, width=barWidth, edgecolor='white', color=c)
#    plt.xticks(r, ['Historical','WWTP','Leak','Basin'],rotation=35,ha='right')
#    plt.ylabel(y_label[o])
##    plt.ylim(ylims[o])
#    plt.title(obj_title[o],fontweight='bold', y=1.04)
#
##fig.tight_layout()
## Flip subsidence measure to be minimizing
##normalized_o[:,1] = 1 - normalized_o[:,1]
#plt.gcf().subplots_adjust(wspace=0.45,left=0.09,right=.97,bottom=0.15,top=.9)
#plt.savefig('model_files\output\plots\Objectives_190702.svg')
#plt.savefig('model_files\output\plots\Objectives_190702.png', dpi=600)
#plt.show()
#
def process_hobs(s_name):
    '''
    Imports head observations from .hob.out file which gives simulated and observed head values
    '''
    obsformation = pd.read_csv('data_raw\obs_formation.csv')
    df = pd.read_fwf('model_files\modflow\VM_'+s_name+'.hob.out',widths=[22,19,22])
    df.columns = ['simulated','observed','obs_name']
    df['time_series'] = np.nan
    df['obs_id'] = np.nan
    df['dem'] = np.nan
    df['geo'] = np.nan
    df['abssimulated'] = np.nan
    df['absobserved'] = np.nan
    df['abssimulated1'] = np.nan
    df['absobserved1'] = np.nan
    
    with open(r'model_files\modflow\VM_OBS.pickle', 'rb') as handle:
        obsinfo = pickle.load(handle)
    
    for i in obsinfo:
        oname = i.obsname
        
        t = np.ones(len(df[df['obs_name'].str.contains(oname)].index))*np.nan
        t[:i.nobs] = i.time_series_data.astype(int).view(int)
    
        df.loc[df['obs_name'].str.contains(oname),'time_series'] = t
        df.loc[df['obs_name'].str.contains(oname),'obs_id'] = oname
        
        df.loc[df['obs_name'].str.contains(oname),'abssimulated'] = df[df['obs_name'].str.contains(oname)][1:]['simulated'] + df[df['obs_name'].str.contains(oname)]['simulated'].values[0]
        df.loc[df['obs_name']==(oname),'abssimulated'] = df[df['obs_name']==oname]['simulated']
        df.loc[df['obs_name'].str.contains(oname),'absobserved'] = df[df['obs_name'].str.contains(oname)][1:]['observed'] + df[df['obs_name'].str.contains(oname)]['observed'].values[0]
        df.loc[df['obs_name']==(oname),'absobserved'] = df[df['obs_name']==oname]['observed']
        
        df.loc[df['obs_name']==(oname+'_1'),'abssimulated1'] = df[df['obs_name']==(oname+'_1')]['simulated']
        df.loc[df['obs_name']==(oname+'_1'),'absobserved1'] = df[df['obs_name']==(oname+'_1')]['observed']
    
    geology = ['Lacustrine','Alluvial','Basalt','Volcaniclastic','Andesite']
    for i, r in obsformation.iterrows():
        df.loc[df['obs_id']==r['IDPOZO'],'dem'] = r['DEM_c3_INE']
        df.loc[df['obs_id']==r['IDPOZO'],'geo'] = geology[r['GEOLOGY_ZO']-1]
    return df, obsformation

def plt_simvsobs(s_name, filename, df=0, obsformation=0):
    
    if not isinstance(df, pd.DataFrame):
        df, obsformation = process_hobs(s_name)
        
    # get coeffs of linear fit
    x = df['absobserved1'].values
    x = x[~np.isnan(x)]
    y = df['abssimulated1'].values
    y = y[~np.isnan(y)]
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    
    sns.set_style("whitegrid")
    g = sns.lmplot(x='absobserved1', y='abssimulated1', hue='geo',data=df,legend=False,palette=dict(Alluvial=(1,0.867,0), Basalt=(0,0.788,0.498) ,Volcaniclastic=(0.9, 0, 0.455) ),size=5,aspect=1.2,scatter_kws={'edgecolor':"none",'s':10, 'alpha':0.3})
    plt.plot(np.linspace(2000,2800,1000), np.linspace(2000,2800,1000), 'k',linestyle=':')
    plt.plot(np.linspace(2000,2800,1000), intercept + slope*np.linspace(2000,2800,1000), 'grey', linewidth=2)
    plt.legend(['Tarango (Volcaniclastic)','Alluvial','Fractured Basalt','One-to-one',"y = {0:.3f}x + {1:.1f}".format(slope,intercept)], loc='upper left')
    plt.xlim(2100,2650)
    plt.ylim(2100,2650)
    plt.xlabel('Observed Head (masl)')
    plt.ylabel('Simulated Head (masl)')
    plt.savefig(filename, dpi=600)
    #plt.close()
    
    return slope, intercept, r_value, p_value, std_err 

def plt_simvsobsddn(s_name, filename, df=0, obsformation=0):
    
    if not isinstance(df, pd.DataFrame):
        df, obsformation = process_hobs(s_name)

    df.loc[df['simulated']>1000,'simulated'] = 0
    df.loc[df['observed']>1000, 'observed'] = 0
    
    # get coeffs of linear fit
    x = df['observed'].values
    y = df['simulated'].values
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    
    sns.set_style("whitegrid")
    g = sns.lmplot(x='observed', y='simulated', hue='geo',data=df,legend=False,palette=dict(Alluvial=(1,0.867,0), Basalt=(0,0.788,0.498) ,Volcaniclastic=(0.9, 0, 0.455) ),size=5,aspect=1.2,scatter_kws={'edgecolor':"none",'s':10, 'alpha':0.3})
    plt.plot(np.linspace(-60,60,1000), np.linspace(-60,60,1000), 'k',linestyle=':')
    plt.plot(np.linspace(-60,60,1000), intercept + slope*np.linspace(-60,60,1000), 'grey', linewidth=2)
    plt.legend(['Tarango (Volcaniclastic)','Alluvial','Fractured Basalt','One-to-one',"y = {0:.3f}x + {1:.1f}".format(slope,intercept)],loc='upper left')
    plt.xlim(-60,60)
    plt.ylim(-60,60)
    plt.xlabel('Observed Drawdown (m)')
    plt.ylabel('Simulated Drawdown (m)')
    plt.savefig(filename, dpi=600)
    
    return slope, intercept, r_value, p_value, std_err 
    
    #plt.close()
## use line_kws to set line label for legend
#ax = sns.regplot(x='abssimulated', y='absobserved', data=df, color='b', 
# line_kws={'label':"y={0:.1f}x+{1:.1f}".format(slope,intercept)})
#ax.legend()
#
#plt.show()
##sns.scatterplot(x='abssimulated', y='absobserved', hue='geo',data=df)
##ax.set_ylim(2100,2700)
##ax.set_xlim(2100,2700)
#
#



def plt_wellhydrographs(s_name, filelocation, df=0, obsformation=0):
    
    if not isinstance(df, pd.DataFrame):
        df, obsformation = process_hobs(s_name)
        
    for o in obsformation['IDPOZO']:
        ddn_data = df[df['obs_id']==o].copy()
        ddn_data['time_series'] = pd.to_timedelta(ddn_data['time_series'],'d') + pd.to_datetime('1984-01-01')
        ddn_data = ddn_data.set_index('time_series')
        ddn_data['simulated'][0] = 0
        ddn_data['observed'][0] = 0
        
        fig, axes = plt.subplots(2, 1, figsize=(5,12))
        ddn_data['abssimulated'].plot(ax=axes[0])
        ddn_data['absobserved'].plot(ax=axes[0])
        axes[0].legend(['Simulated','Observed'])
        axes[0].xaxis.label.set_visible(False)
        
        ddn_data['simulated'].plot(ax=axes[1])
        ddn_data['observed'].plot(ax=axes[1])
        axes[1].set_ylim([-50, 20])
        axes[1].legend(['Simulated','Observed'])
        axes[1].xaxis.label.set_visible(False)
        
        plt.tight_layout()
        plt.savefig(filelocation+o+'.png', dpi=300)
        plt.close()
        