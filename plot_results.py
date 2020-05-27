# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 14:02:08 2018

@author: MM
"""
import flopy
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import pickle
import flopy.utils.binaryfile as bf
from scipy import stats
from pathlib import Path
from mpl_toolkits.mplot3d import Axes3D

sns.set(style="white", palette="muted", color_codes=True)
plt.rcParams['legend.fontsize'] = 20
plt.rcParams['axes.titlesize'] = 22
plt.rcParams['axes.labelsize'] = 22
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['figure.titlesize'] = 24
plt.rcParams.update({'font.size': 20})

xll = 455000
yll = 2107000
xur = 539000
yur = 2175000
cellsize = 500

STRT_YEAR = 1984
END_YEAR = 2014

ncol = int((xur-xll)/cellsize) # Number of rows
nrow = int((yur-yll)/cellsize) # Number of columns

## Head Dictionary
def get_heads(alt_list,safolder=-1):
    S_heads = {}
    if safolder >= 0: 
        headpath = Path.cwd().joinpath('model_files').joinpath('modflow').joinpath(str(safolder))
    else:
        headpath = Path.cwd().joinpath('model_files').joinpath('modflow')
    for name in alt_list:
        S_heads[name] = bf.HeadFile(headpath.joinpath(name+'.hds'))

    return S_heads

def calc_head_change(hds, GEO, ACTIVE, n, m, g_units, lyr):
    # Get head values at the nth and mth time steps
    h_i = hds.get_data(mflay=lyr,kstpkper=n)
    h_e = hds.get_data(mflay=lyr,kstpkper=m)
    
## Heads Contour
def plt_head_change(s_heads, GEO, ACTIVE_LYR1, ACTIVE_LYR2, n=(30,0), m=(30,359), g_units = [2,3,4], lyr = 1):
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

def plt_head_change(alt_list, mapTitles, s_heads, GEO, ACTIVE, n = (8,23), m = (8,359), g_units = [2,3,4,5], lyr = 1):
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
    
    for i,name in enumerate(alt_list):
        hds = s_heads[name]
        h_i = hds.get_data(mflay=1,kstpkper=n)
        h_e = hds.get_data(mflay=1,kstpkper=m)

    print('Plotting head changes over model period...')
    plt.set_cmap('viridis_r')

    # Get head values at the nth and mth time steps
    first_heads = s_heads[list(s_heads.keys())[0]]

    first_change = calc_head_change(first_heads, GEO, ACTIVE, n, m, g_units, lyr)

    for i, name in enumerate(alt_list):
        alt_heads = s_heads[name]

        alt_change = calc_head_change(alt_heads, GEO, ACTIVE, n, m, g_units, lyr)

        change_compare = alt_change - first_change

        fig, axes = plt.subplots()
        cbar_ax = fig.add_axes()

        im = axes.imshow(change_compare,vmin=0,vmax=20)
        CS = axes.contour(ACTIVE, colors='k', linewidths=2)

        axes.xaxis.set_visible(False)
        axes.yaxis.set_visible(False)
        axes.set_title(mapTitles[i].format(i+1))

        fig.colorbar(im, cax=cbar_ax, label='Change in Groundwater Head (m)')

        plt.savefig(Path.cwd().joinpath('model_files').joinpath('output').joinpath('plots').joinpath('head-change_'+name+'-'+alt_list[0]+'.eps'), dpi=600)
        plt.savefig(Path.cwd().joinpath('model_files').joinpath('output').joinpath('plots').joinpath('head-change_'+name+'-'+alt_list[0]+'.png'), dpi=600)
        plt.close()

def plt_alt_objectives(alt_names, num_alt, objectives):
    '''
    objectives is a list with an array of length number of alternatives for each objective
    '''
    print('Plotting alternative performance under objectives...')
    plt.rcParams['legend.fontsize'] = 20
    plt.rcParams['axes.titlesize'] = 22
    plt.rcParams['axes.labelsize'] = 22
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 18
    plt.rcParams['figure.titlesize'] = 24
    plt.rcParams.update({'font.size': 20})
    
    c = ['k','goldenrod','blue','darkgreen']
    barWidth = 0.1
    r = np.arange(num_alt)*0.1 # bar position
    y_label = ['Pumping Energy (kWh)','Depth to Groundwater in Clay (m)','Percent of Urban Cells Flooded']
    obj_title = ['Energy Use','Subsidence Avoidance','Urban Flooding']

    normalized_o = np.zeros((num_alt, len(objectives)))

    fig, axes = plt.subplots(nrows=1, ncols=3,figsize=(12,7.2))

    for o, obj in enumerate(objectives):
        normalized_o[:,o] = obj / (obj.max(axis=0) - obj.min(axis=0))
        plt.subplot(1,len(objectives),o+1)
        plt.bar(r, obj, width=barWidth, edgecolor='white', color=c)
        plt.xticks(r, alt_names, rotation=35, ha='right')
        plt.ylabel(y_label[o])
        plt.title(obj_title[o], fontweight='bold', y=1.04)

    # Flip subsidence measure to be minimizing
    plt.gcf().subplots_adjust(wspace=0.45,left=0.09,right=.97,bottom=0.15,top=.9)
    plt.savefig(Path.cwd() / 'model_files' / 'output' / 'plots' / 'Objectives.eps', dpi=600)
    plt.savefig(Path.cwd() / 'model_files' / 'output' / 'plots' / 'Objectives.png', dpi=600)
    plt.close()

def parallel_axis(nondom_results, obj_labels, opt_run):
    # Plots a normalized parallel axis
    plt.figure()
    for ppoint in nondom_results:
        ppoint = (ppoint - nondom_results.min(axis=0)) / (nondom_results.max(axis=0) - nondom_results.min(axis=0))
        plt.plot(range(len(obj_labels)), ppoint, 'steelblue')

    plt.gca().set_xticks(range(len(obj_labels)))
    plt.gca().set_xticklabels(obj_labels)
    plt.savefig(Path.cwd().joinpath('parallelaxis_' + opt_run + '.png'))
    plt.close()

def get_budgets(alt_list, mapTitles, s_heads):
    # Budget
    df_Bdget = {}
    df_CumSum = {}
    
    for name in alt_list:
        mf_list = flopy.utils.MfListBudget(Path.cwd().joinpath('model_files').joinpath('modflow').joinpath(name+".list"))
        incremental, cumulative = mf_list.get_budget()
    
        df_Bdget[name], df_CumSum[name] = mf_list.get_dataframes(start_datetime="04-30-1984")
        
        mthly_Bdget = df_Bdget[name].drop(['CONSTANT_HEAD_IN', 'TOTAL_IN', 'CONSTANT_HEAD_OUT', 'RECHARGE_OUT', 'TOTAL_OUT', 'IN-OUT', 'PERCENT_DISCREPANCY'], axis=1)
        
        mthly_Bdget['STORAGE_OUT'] = mthly_Bdget['STORAGE_OUT'].apply(lambda x: x*-1)
        mthly_Bdget['WELLS_OUT'] = mthly_Bdget['WELLS_OUT'].apply(lambda x: x*-1)
#        mthly_Bdget['DRAINS_OUT'] = mthly_Bdget['DRAINS_OUT'].apply(lambda x: x*-1)
        mthly_Bdget = mthly_Bdget.multiply(30/1000000)
        cols = mthly_Bdget.columns.tolist()
        # reorder columns
        cols = [cols[1]] + [cols[2]] + [cols[0]] + [cols[4]] + [cols[3]] 
        # "commit" the reordering
        mthly_Bdget = mthly_Bdget[cols]
        
        im = axes[i].imshow(hist_compare,vmin=0,vmax=20)
        CS = axes[i].contour(ACTIVE_LYR1, colors='k', linewidths=2)
        axes[i].xaxis.set_visible(False)
        axes[i].yaxis.set_visible(False)
        axes[i].set_title(mapTitle[i].format(i+1))
        
        fig.subplots_adjust(right=0.8)
        fig.colorbar(im, cax=cbar_ax, label='Change in Groundwater Head (m)')
        plt.savefig(Path.cwd() / 'model_output' / 'plots' / 'Hist_change_all.svg')
        plt.savefig(Path.cwd() / 'model_output' / 'plots' / 'Hist_change_all.png', dpi=600)
        plt.show()

    for name in alt_list:
        df_CumSum[name]['IN'] = df_CumSum[name]['RECHARGE_IN'].divide(1000000) + df_CumSum[name]['WELLS_IN'].divide(1000000)
        df_CumSum[name]['OUT'] = df_CumSum[name]['WELLS_OUT'].divide(1000000) #+ df_CumSum[name]['DRAINS_OUT'].divide(1000000)
        df_CumSum[name]['INOUTCUMSUM'] = df_CumSum[name]['IN'] - df_CumSum[name]['OUT']
    
    return df_Bdget, mthly_Bdget, df_CumSum
    
def plt_cum_sum(filename, alt_list, mapTitles, df_CumSum, start='01-31-1985', end='12-31-2013'):
    # Plotting defaults
    l = [4,2,2,2]
    c = ['k','goldenrod','blue','darkgreen']
    mark = ['-','-','-','--']
    
    for i,name in enumerate(alt_list):
        df_CumSum[name].INOUTCUMSUM[start:end].plot(linewidth=l[i],color=c[i],style=mark[i])
        
    plt.ylabel(r'Volume (million m$^3$)')
    #plt.xlabel(r'Year')
    #plt.ylim(-8,10)
    #plt.title('Cumulative In - Out')
    plt.legend(mapTitles,loc='lower left')
    plt.gcf().subplots_adjust(left=0.15,right=.95,bottom=0.1,top=.95)
    
    plt.savefig(filename+'.png', dpi=600)
    plt.savefig(filename+'.eps', dpi=600)
    
    plt.show()

def process_hobs(folder, name, legend=['Lacustrine','Alluvial','Basalt','Volcaniclastic','Andesite'], obsinfo_loaded=True):
    '''
    Imports head observations from .hob.out file which gives simulated and observed head values
    '''
    obsstats = pd.read_csv(Path.cwd() / 'model_files' / 'modflow' / 'OBS_stats.csv')
    obsformation = pd.read_csv(Path.cwd() / 'data_raw' / 'obs_formation.csv')
    df = pd.read_fwf(Path.cwd().joinpath('model_files').joinpath('output').joinpath('sa').joinpath('hob').joinpath(folder).joinpath(name+'.hob_out'),widths=[22,19,22])
    df.columns = ['simulated','observed','obs_name']
    df['LAT'] = np.nan
    df['LON'] = np.nan
    df['time_series'] = np.nan
    df['obs_id'] = np.nan
    df['stat'] = obsstats
    df['formation'] = np.nan
    df['abssimulated'] = np.nan
    df['absobserved'] = np.nan
    df['abssimulated1'] = np.nan
    df['absobserved1'] = np.nan
    
    if obsinfo_loaded:
        print('Opening observation input file...')
        with open(str(Path.cwd() / 'model_files' / 'modflow' / 'OBS.pickle'), 'rb') as handle:
            obsinfo = pickle.load(handle)
    else:
        print('Processing observation input file...')
        mf = flopy.modflow.Modflow.load(str(Path.cwd().joinpath('model_files').joinpath('modflow').joinpath('Historical.nam')))
        hob = flopy.modflow.ModflowHob.load(str(Path.cwd() / 'model_files' / 'modflow' / 'OBS.ob_hob'), mf)
        winfofile = str(Path.cwd() / 'model_files' / 'modflow' / 'OBS.pickle')
        with open(winfofile, 'wb') as handle:
            pickle.dump(hob.obs_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
        obsinfo = hob.obs_data
    
    for i in obsinfo:
        oname = i.obsname
        
        t = np.ones(len(df[df['obs_name'].str.contains(oname)].index))*np.nan
        t[:i.nobs] = [x[0] for x in i.time_series_data]
    
        df.loc[df['obs_name'].str.contains(oname),'time_series'] = t
        df.loc[df['obs_name'].str.contains(oname),'obs_id'] = oname
        
        df.loc[df['obs_name'].str.contains(oname),'abssimulated'] = df[df['obs_name'].str.contains(oname)][1:]['simulated'] + df[df['obs_name'].str.contains(oname)]['simulated'].values[0]
        df.loc[df['obs_name']==(oname),'abssimulated'] = df[df['obs_name']==oname]['simulated']
        df.loc[df['obs_name'].str.contains(oname),'absobserved'] = df[df['obs_name'].str.contains(oname)][1:]['observed'] + df[df['obs_name'].str.contains(oname)]['observed'].values[0]
        df.loc[df['obs_name']==(oname),'absobserved'] = df[df['obs_name']==oname]['observed']
        
        df.loc[df['obs_name']==(oname+'_1'),'abssimulated'] = df[df['obs_name']==(oname+'_1')]['simulated']
        df.loc[df['obs_name']==(oname+'_1'),'absobserved'] = df[df['obs_name']==(oname+'_1')]['observed']
        df.loc[df['obs_name']==(oname+'_1'),'abssimulated1'] = df[df['obs_name']==(oname+'_1')]['simulated']
        df.loc[df['obs_name']==(oname+'_1'),'absobserved1'] = df[df['obs_name']==(oname+'_1')]['observed']
    
    for i, r in obsformation.iterrows():
        df.loc[df['obs_id']==r['IDPOZO'],'formation'] = legend[r['ZONE']-1]
        df.loc[df['obs_id']==r['IDPOZO'],'LON'] = r['X']
        df.loc[df['obs_id']==r['IDPOZO'],'LAT'] = r['Y']
    
    return df, obsinfo, obsstats, obsformation

def plt_cluster(df, km_cluster, filename):
    fig = plt.figure(figsize=(16, 12))
    ax = fig.add_subplot(111, projection='3d')
    
    xs = df['LON']
    ys = df['LAT']
    zs = df['Z']
    
#    if time_marker:
#        markers = ['o','^','s']
#        for i in range(3):
#            ax.scatter(xs[df['t_group'].values==i], ys[df['t_group'].values==i], zs[df['t_group'].values==i], s=50, alpha=0.6, edgecolors=None, c=km_cluster[df['t_group'].values==i], cmap='Set1', marker=markers[i])
#    else:
#        ax.scatter(xs, ys, zs, s=50, alpha=0.6, c=km_cluster, cmap='Set1', edgecolors=None)
    
    ax.scatter(xs, ys, zs, s=50, alpha=0.6, c=km_cluster, cmap='Set1', edgecolors=None)
    
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('Elevation (masl)')
    
    ax.set_xlim(xs.min(), xs.max())
    ax.set_ylim(ys.min(), ys.max())
    ax.set_zlim(zs.min(), zs.max())
    
    fig.tight_layout()

    fig.savefig(filename)
    plt.close()

def plt_wellhydrographs(folder, name, filelocation, df=0, obsformation=0, obsinfo_loaded=True, timestep='d', startdate='1984-01-01', ddn_lim=[-50, 20], legend=['Lacustrine','Alluvial','Basalt','Volcaniclastic','Andesite']):
    
    if not isinstance(df, pd.DataFrame):
        df, obsinfo, obsstats, obsformation = process_hobs(folder, name, legend=legend, obsinfo_loaded=obsinfo_loaded)
    
    filelocation = Path.cwd().joinpath('model_files').joinpath('output').joinpath('plots').joinpath('observations').joinpath(filelocation)
    filelocation.mkdir(exist_ok=True)
    
    for l in legend:
        filelocation.joinpath(l).mkdir(exist_ok=True)
    
    for o in obsformation['IDPOZO']:
        print(o)
        ddn_data = df[df['obs_id']==o].copy()
        ddn_data['time_series'] = pd.to_timedelta(ddn_data['time_series'], timestep) + pd.to_datetime(startdate)
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
        axes[1].set_ylim(ddn_lim)
        axes[1].legend(['Simulated','Observed'])
        axes[1].xaxis.label.set_visible(False)
        
        plt.tight_layout()
        filename = filelocation.joinpath(ddn_data['formation'][0]).joinpath(o+'.png')
        print(filename)
        plt.savefig(str(filename), dpi=600)
        plt.close()
        
def plt_simvsobs(name, filename, df=0, obsformation=0, obsinfo_loaded=True):
    
    filename = str(Path.cwd() / 'model_files' / 'output' / 'plots' / 'calibration' / filename)
    
    if not isinstance(df, pd.DataFrame):
        df, obsinfo, obsstats, obsformation = process_hobs(name, obsinfo_loaded=obsinfo_loaded)
        
    # get coeffs of linear fit
    x = df['absobserved1'].values
    x = x[~np.isnan(x)]
    y = df['abssimulated1'].values
    y = y[~np.isnan(y)]
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    
    sns.set_style("whitegrid")
    g = sns.lmplot(x='absobserved1', y='abssimulated1', hue='formation',data=df,legend=False,palette=dict(Alluvial=(1,0.867,0), Basalt=(0,0.788,0.498), Volcaniclastic=(0.9, 0, 0.455) ), size=8, aspect=1.2, scatter_kws={'edgecolor':"none",'s':10, 'alpha':0.3})
    plt.plot(np.linspace(2000,2800,1000), np.linspace(2000,2800,1000), 'k',linestyle=':')
    plt.plot(np.linspace(2000,2800,1000), intercept + slope*np.linspace(2000,2800,1000), 'grey', linewidth=2)
    plt.legend(['Tarango (Volcaniclastic)','Alluvial','Fractured Basalt','One-to-one',"y = {0:.3f}x + {1:.1f}, R = {2:.3f} ".format(slope, intercept, r_value)], loc='upper left')
    plt.xlim(2100,2650)
    plt.ylim(2100,2650)
    plt.xlabel('Observed Head (masl)')
    plt.ylabel('Simulated Head (masl)')
    plt.savefig(filename+'.png', dpi=600)
    plt.savefig(filename+'.eps', dpi=600)
    #plt.close()
    
    return slope, intercept, r_value, p_value, std_err 

def plt_simvsobsddn(name, filename, df=0, obsformation=0, obsinfo_loaded=True):
    
    filename = str(Path.cwd() / 'model_files' / 'output' / 'plots' / 'calibration' / filename)
    
    if not isinstance(df, pd.DataFrame):
        df, obsinfo, obsstats, obsformation = process_hobs(name, obsinfo_loaded=obsinfo_loaded)

    df.loc[df['simulated']>1000,'simulated'] = 0
    df.loc[df['observed']>1000, 'observed'] = 0
    
    # get coeffs of linear fit
    x = df['observed'].values
    y = df['simulated'].values
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    
    sns.set_style("whitegrid")
    g = sns.lmplot(x='observed', y='simulated', hue='formation',data=df,legend=False,palette=dict(Alluvial=(1,0.867,0), Basalt=(0,0.788,0.498), Volcaniclastic=(0.9, 0, 0.455)), size=8, aspect=1.2, scatter_kws={'edgecolor':"none",'s':10, 'alpha':0.3})
    plt.plot(np.linspace(-60,60,1000), np.linspace(-60,60,1000), 'k',linestyle=':')
    plt.plot(np.linspace(-60,60,1000), intercept + slope*np.linspace(-60,60,1000), 'grey', linewidth=2)
    plt.legend(['Tarango (Volcaniclastic)','Alluvial','Fractured Basalt','One-to-one',"y = {0:.3f}x + {1:.1f}, R = {2:.3f} ".format(slope, intercept, r_value)],loc='upper left')
    plt.xlim(-60,60)
    plt.ylim(-60,60)
    plt.xlabel('Observed Drawdown (m)')
    plt.ylabel('Simulated Drawdown (m)')
    plt.savefig(filename+'.png', dpi=600)
    plt.savefig(filename+'.eps', dpi=600)
    
    return slope, intercept, r_value, p_value, std_err 
#
#def plt_sensitivity(sensitivity):
#    