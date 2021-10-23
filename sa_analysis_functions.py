# -*- coding: utf-8 -*-
"""
Created on Mon May 11 13:44:44 2020

@author: MM
"""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pathlib import Path
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
        
sns.set_style("darkgrid", {"axes.facecolor": ".93"})

plt.rcParams['legend.fontsize'] = 20
plt.rcParams['axes.titlesize'] = 22
plt.rcParams['axes.labelsize'] = 22
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['figure.titlesize'] = 24
plt.rcParams.update({'font.size': 20})

class dataset():

    # Initializer / Instance attributes
    def __init__(self, filename, n_samples=100000, clusters=[[1],[2],[3],[4],[5],[1,2,3,4,5]], altnames=['Hist','WWTP','Basin','Leak'], altnameslong=['Historical','Wastewater Reuse','Infiltration Basins','Repair Leaks'], objsymbol=['E','W','F'], objnames=['Pumping Energy','Water Quality','Urban Flooding'], parnames=['HK_1','HK_2','HK_3','HK_4','HK_5','SS_1','SS_2','SS_3','SS_4','SS_5','SY_1','SY_2','SY_3','SY_4','SY_5','VANI_1','VANI_2','VANI_3','VANI_4','VANI_5','Q_1','Q_2','Q_3','LK_1','LK_2','LK_3','TWU_1','TWU_2','TWU_3','RCH_1','RCH_2','RCH_3','IN_1'], parnameslong=['Lacustrine Hydraulic Conductivity','Alluvial Hydraulic Conductivity','Basaltic Hydraulic Conductivity','Volcaniclastic Hydraulic Conductivity','Andesitic Hydraulic Conductivity','Lacustrine Specific Storage','Alluvial Specific Storage','Basaltic Specific Storage','Volcaniclastic Specific Storage','Andesitic Specific Storage','Lacustrine Specific Yield','Alluvial Specific Yield','Basaltic Specific Yield','Volcaniclastic Specific Yield','Andesitic Specific Yield','Vertical Anisotropy of\nLacustrine Hydraulic Conductivity','Vertical Anisotropy of\nAlluvial Hydraulic Conductivity','Vertical Anisotropy of\nBasaltic Hydraulic Conductivity','Vertical Anisotropy of\nVolcaniclastic Hydraulic Conductivity','Vertical Anisotropy of\nAndesitic Hydraulic Conductivity','Urban to Periurban\nPumping Multiplier 1990','Urban to Periurban\nPumping Multiplier 2000','Urban to Periurban\nPumping Multiplier 2010','Leak Multiplier 1990','Leak Multiplier 2000','Leak Multiplier 2010','Total Water Use\nMultiplier 1990','Total Water Use\nMultiplier 2000','Total Water Use\nMultiplier 2010','Urban Land Use\nRecharge Multiplier','Natural Land Use\nRecharge Multiplier','Wetland Land Use\nRecharge Multiplier','Leak Infiltration Rate']):
        self.filename = filename # Assign name
        self.n_samples = n_samples
        self.objsymbol = objsymbol
        self.objnames = objnames
        self.altnames = altnames
        self.altnameslong = altnameslong
        self.obj_alt = []
        for o in objsymbol:
            for a in self.altnames:
                self.obj_alt.extend([o+'-'+a])
        self.n_obj_alt = len(self.obj_alt)
        self.clusters = clusters
        for i in range(len(self.clusters)):
            temp = 'C-'
            for j in range(5-len(self.clusters[i])):
                temp += '0'
            for j in self.clusters[i]:
                temp += str(j)
            self.clusters[i] = temp
        self.n_clusters = len(self.clusters)
        self.parnames = parnames
        self.n_params = len(self.parnames)
        self.parnameslong = parnameslong

        self.allnames = ['sample']
        for i in [self.parnames, self.obj_alt, self.clusters]:
            self.allnames.extend(i)
        
        self.alldata = pd.read_csv(Path.cwd().joinpath('SA_data').joinpath(self.filename+'-samples_evald.csv'), names=self.allnames)
        self.ran_cluster_kde = False
        self.ran_norm_param = False

    def plt_cluster(self):
        # Visualize spatial clusters in 3 dimensions
        sns.set_style("whitegrid")
        colors = ['xkcd:maroon','xkcd:orange','xkcd:lightblue','xkcd:olive','xkcd:violet']
        df = pd.read_csv(Path.cwd().joinpath('SA_data').joinpath(self.filename+'-cluster_detailed.csv'))
        fig = plt.figure(figsize=(16, 12))
        ax = fig.add_subplot(111, projection='3d')
        
        df['Colors'] = df['Cluster'].apply(lambda x: colors[x])
        
        xs = df['LON']
        ys = df['LAT']
        zs = df['Z']

        ax.scatter(xs, ys, zs, s=50, alpha=0.6, c=df['Colors'], edgecolors=None)
        
        ax.set_xlabel('Longitude', labelpad=20)
        ax.set_ylabel('Latitude', labelpad=20)
        ax.set_zlabel('Elevation (masl)', labelpad=20)
        
        ax.set_xlim(xs.min(), xs.max())
        ax.set_ylim(ys.min(), ys.max())
        ax.set_zlim(zs.min(), zs.max())
        
        fig.savefig(Path.cwd().joinpath('images').joinpath('sa').joinpath(self.filename+'_Cluster.png'))
        fig.savefig(Path.cwd().joinpath('images').joinpath('sa').joinpath(self.filename+'_Cluster.svg'))
        plt.close()
    
    def normalize_paramrange(self):
        # Parameter KDEs by error percent threshold for select cluster
        # Normalize parameter ranges by min/max of each parameter range
        self.norm_params = self.alldata[self.parnames].copy()
        for i in self.parnames:
            temp_array = self.norm_params[i].values
            self.norm_params[i] = (temp_array - np.amin(temp_array)) / (np.amax(temp_array) - np.amin(temp_array))
        self.norm_params = pd.concat([self.norm_params, self.alldata[self.allnames[self.n_params+1:]]],axis=1)
    
    def cluster_kde(self, params=-1, threshold=-1, n_samples=100000, plot=True):
        # Set the threshold to all clean samples if no threshold provided
        if threshold < 0:
            self.threshold = self.alldata.shape[0]/self.n_samples
        else:
            self.threshold = threshold
        
        if not self.ran_norm_param:
            self.normalize_paramrange()
        
        cols = self.clusters.copy()
        cols.append('Threshold')
        
        # Load the indices of the threshold percent of samples with the least error in each cluster
        self.thresh_samples = pd.read_csv(Path.cwd().joinpath('SA_data').joinpath(self.filename+'-threshold_samples.csv'), names=cols)
        
        if type(params) is int:
            params = self.parnames
            
        self.param_by_cluster = pd.DataFrame()
        for p in params:
            for i in self.clusters:
                # Create dataframe with parameter values for each cluster
                params = self.norm_params[p].values[self.thresh_samples[i][self.thresh_samples['Threshold'] == threshold].values.astype('int')]
                cluster_p_dist = pd.DataFrame(params,columns=['Value'])
                cluster_p_dist['Cluster'] = i
                cluster_p_dist['Param'] = p
                self.param_by_cluster = pd.concat([self.param_by_cluster,cluster_p_dist])
        
        if plot:
            # Subplot for each parameter, KDE for each cluster
            g = sns.FacetGrid(self.param_by_cluster, col='Param', hue='Cluster', palette='viridis_r', col_wrap=5)
            g = (g.map(sns.kdeplot, 'Value'))
            g.set(xlim=(0, 1), xticks=[0,1], xticklabels=['min','max'], ylim=(0,4.5), ylabel='Density', xlabel='', )
            ## Iterate thorugh each axis
            for ax in g.axes.flat:
                # Make title more human-readable and larger
                if ax.get_title():
                    ax.set_title(ax.get_title().split('=')[1], fontsize='xx-large')
                if ax.get_ylabel():
                    ax.set_ylabel(ax.get_ylabel(), fontsize='xx-large')
                if ax.get_xticklabels():
                    ax.set_xticklabels(ax.get_xticklabels(), fontsize='x-large')
                if ax.get_yticklabels():
                    ax.set_yticklabels(ax.get_yticklabels(), fontsize='x-large')
            
            plt.subplots_adjust(left=0.1)
            plt.legend(fontsize = 'x-large', title_fontsize = 'xx-large', loc='lower right', facecolor='white', ncol=2, bbox_to_anchor=(1.05, 1))
            
            plt.savefig(Path.cwd().joinpath('images').joinpath('sa').joinpath('KDE-by-param_th-'+str(self.threshold)+'.png'))
            plt.savefig(Path.cwd().joinpath('images').joinpath('sa').joinpath('KDE-by-param_th-'+str(self.threshold)+'.svg'))
            plt.close()
        
        self.ran_cluster_kde = True
    
    def error_scatter(self, params=-1, threshold=-1, n_samples=100000, show_all=False):
        # Scatter plot showing model sum of squared weight residuals showing all simulation runs in gray and filtered best performing parameter sets in orange
        if not self.ran_cluster_kde:
            self.cluster_kde(params=params, threshold=threshold, n_samples=n_samples, plot=False)
        
        # Create list of labels for figures
        if type(params) is int:
            params = self.parnames
            paramlabels = self.parnameslong
        else:
            paramlabels = []
            for p in params:
                paramlabels.extend([self.parnameslong[self.parnames.index(p)]])
            
        color = 'orange'
        ob_loc = Path.cwd() / 'Error'
        ob_loc.mkdir(exist_ok=True)
            
        for p, param in enumerate(params):
                            
            fig, axes = plt.subplots(nrows=1, ncols=len(self.clusters), sharey=True, figsize=(20, 6))
            fig.suptitle('Sum of Squared Weighted Residuals versus ' + paramlabels[p])
            fig.subplots_adjust(left=0.12, top=0.7)
            
            # Create dataframe with parameter values for each cluster
            for i, c in enumerate(self.clusters):
                df = pd.DataFrame(self.param_by_cluster.loc[self.param_by_cluster['Param'].values==param].values, columns=['Value','Cluster','Param'])
                df = pd.DataFrame(df.loc[df['Cluster'].values==c].values, columns=['Value','Cluster','Param'])
                df['Value'] = self.norm_params[param].values[self.thresh_samples[c][self.thresh_samples['Threshold']==self.threshold].values.astype('int')].copy()
                df[c] = self.alldata[c].values[self.thresh_samples[c][self.thresh_samples['Threshold']==self.threshold].values.astype('int')].copy()
                
                if show_all:
                    self.norm_params.plot(kind='scatter', x=param, y=c, c='silver', s=5, ax=axes[i], colorbar=False)
                df.plot(kind='scatter', x='Value', y=c, c=color, s=5, ax=axes[i], colorbar=False, title=c)
                
                axes[i].set(yscale='log')
#                        axes[a,i].set_ylim(bottom=yllim[o], top=yulim[o])
                axes[i].set_ylabel('')
                
                axes[i].set_xticks([0,1])
                axes[i].set_xticklabels(['min','max'])
                axes[i].set_xlabel(param)
                
                if i == 0:
                    axes[i].set_ylabel('Error')
            
            plt.savefig(Path.cwd().joinpath('images').joinpath('sa').joinpath('error').joinpath(param+'_Error.png'), dpi=600)
            plt.savefig(Path.cwd().joinpath('images').joinpath('sa').joinpath('error').joinpath(param+'_Error.svg'))
            plt.close()
    
    def objective_scatter(self, params=-1, threshold=-1, n_samples=100000, show_all=False):
        # Scatter plot showing objective values showing all simulation runs in gray and filtered best performing parameter sets according to the error metric above in green (), red (), blue ()
        if not self.ran_cluster_kde:
            self.cluster_kde(params=params, threshold=threshold, n_samples=n_samples, plot=False)
        
        # Create list of labels for figures
        if type(params) is int:
            params = self.parnames
            paramlabels = self.parnameslong
        else:
            paramlabels = []
            for p in params:
                paramlabels.extend([self.parnameslong[self.parnames.index(p)]])
            
        color = ['green','red','blue']
        altlabels=['Historical\n','Wastewater\nReuse','Infiltration\nBasins','Repair\nLeaks']
        # Scatter Plots - Objectives
        yllim = [1e8,1e-2,1e-2]
        yulim = [1e11,1e2,1e2] 
        for o, ob in enumerate(self.objsymbol):  
            ob_loc = Path.cwd() / 'images' / ob
            ob_loc.mkdir(exist_ok=True)
            
            for p, param in enumerate(params):
                                
                fig, axes = plt.subplots(nrows=4, ncols=len(self.clusters), sharey=True, figsize=(20, 22))
                fig.suptitle('Objective ' + self.objnames[o] + ' versus ' + paramlabels[p])
                
                for a, alt in enumerate(self.altnames):
                    # Create dataframe with parameter values for each threshold group by cluster
                    for i, c in enumerate(self.clusters):
                        df = pd.DataFrame(self.param_by_cluster.loc[self.param_by_cluster['Param'].values==param].values, columns=['Value','Cluster','Param'])
                        df = pd.DataFrame(df.loc[df['Cluster'].values==c].values, columns=['Value','Cluster','Param'])
                        df['Value'] = self.norm_params[param].values[self.thresh_samples[c][self.thresh_samples['Threshold']==self.threshold].values.astype('int')].copy()
                        df[ob] = self.alldata[ob+'-'+alt].values[self.thresh_samples[c][self.thresh_samples['Threshold']==self.threshold].values.astype('int')].copy()
                        
                        if a == 0:
                            if show_all:
                                self.norm_params.plot(kind='scatter', x=param, y=ob+'-'+alt, c='silver', s=5, ax=axes[a,i], colorbar=False)
                            df.plot(kind='scatter', x='Value', y=ob, c=color[o], s=5, ax=axes[a,i], colorbar=False, title=c)
                        else:
                            if show_all:
                                self.norm_params.plot(kind='scatter', x=param, y=ob+'-'+alt, c='silver', s=5, ax=axes[a,i], colorbar=False)
                            df.plot(kind='scatter', x='Value', y=ob, c=color[o], s=5, ax=axes[a,i], colorbar=False)
        
                        axes[a,i].set(yscale='log')
                        axes[a,i].set_ylim(bottom=yllim[o], top=yulim[o])
                        axes[a,i].set_ylabel('')
                        if a != 3:
                            axes[a,i].set_xticks([])
                            axes[a,i].set_xlabel('')
                        else:
                            axes[a,i].set_xticks([0,1])
                            axes[a,i].set_xticklabels(['min','max'])
                            axes[a,i].set_xlabel(param)
                        
                        if i == 0:
                            axes[a,i].set_ylabel(altlabels[a])
                
                fig.text(0.03, 0.5, self.objnames[o], ha='center', va='center', rotation='vertical')
                fig.subplots_adjust(left=0.12, top=0.85)
                plt.savefig(Path.cwd().joinpath('images').joinpath('sa').joinpath('scatter').joinpath(ob).joinpath(param+'_'+ob+'.png'), dpi=600)
                plt.savefig(Path.cwd().joinpath('images').joinpath('sa').joinpath('scatter').joinpath(ob).joinpath(param+'_'+ob+'.svg'))
                plt.close()

    def historical_scatter(self, params=-1, threshold=-1, n_samples=100000, show_all=False):
        if not self.ran_cluster_kde:
            self.cluster_kde(params=params, threshold=threshold, n_samples=n_samples, plot=False)
        
        # Create list of labels for figures
        if type(params) is int:
            params = self.parnames
            paramlabels = self.parnameslong
        else:
            paramlabels = []
            for p in params:
                paramlabels.extend([self.parnameslong[self.parnames.index(p)]])
                
        color = ['orange','green','red','blue']
        objLong = ['Model Error\n','Pumping Energy\n(kWh)','Water Quality Risk\n(% cells at risk)','Urban Flooding\n(% cells at risk)']
        objShort = ['Error','Energy','Water','Flood']
        
        # Scatter Plots - Objectives
        yllim = [2e3,1e8,1e-2,1e-2]
        yulim = [3e14,7e10,1e2,1e2]
        
        ob_loc = Path.cwd() / 'images' / 'sa' / 'scatter' / 'historical'
        ob_loc.mkdir(exist_ok=True)
            
        for p, param in enumerate(params):
            
            fig, axes = plt.subplots(nrows=4, ncols=len(self.clusters), sharex=True, sharey=False, figsize=(19, 21))
            fig.suptitle('Historical Alternative Error and Objectives vs\n' + self.parnameslong[p])
            
            for o, obj in enumerate(objShort):
                # Create dataframe with parameter values for each threshold group by cluster
                for i, c in enumerate(self.clusters):
                    df = pd.DataFrame(self.param_by_cluster.loc[self.param_by_cluster['Param'].values==param].values, columns=['Value','Cluster','Param'])
                    df = pd.DataFrame(df.loc[df['Cluster'].values==c].values, columns=['Value','Cluster','Param'])
                    df['Value'] = self.norm_params[param].values[self.thresh_samples[c][self.thresh_samples['Threshold']==self.threshold].values.astype('int')].copy()
                    
                    if o==0:
                        df[obj] = self.alldata[c].values[self.thresh_samples[c][self.thresh_samples['Threshold']==self.threshold].values.astype('int')].copy()
                        if show_all:
                                self.norm_params.plot(kind='scatter', x=param, y=c, c='silver', s=5, ax=axes[o,i], colorbar=False)
                        df.plot(kind='scatter', x='Value', y=obj, c=color[o], s=5, ax=axes[o,i], colorbar=False, title=c)
                    else:
                        df[obj] = self.alldata[self.objsymbol[o-1]+'-Hist'].values[self.thresh_samples[c][self.thresh_samples['Threshold']==self.threshold].values.astype('int')].copy()
                        if show_all:
                            self.norm_params.plot(kind='scatter', x=param, y=self.objsymbol[o-1]+'-Hist', c='silver', s=5, ax=axes[o,i], colorbar=False)
                        df.plot(kind='scatter', x='Value', y=obj, c=color[o], s=5, ax=axes[o,i], colorbar=False)
    
                    axes[o,i].set(yscale='log')
                    axes[o,i].set_ylim(bottom=yllim[o], top=yulim[o])
                    axes[o,i].set_ylabel('')
                    if o != 3:
                        axes[o,i].set_xticks([])
                        axes[o,i].set_xlabel('')
                    else:
                        axes[o,i].set_xticks([0,1])
                        axes[o,i].set_xticklabels(['min','max'])
                    
                    if i == 0:
                        axes[o,i].set_ylabel(objLong[o])
                        axes[o,i].set_ylabel(objLong[o])
                    else:
                        axes[o,i].set_yticks([])
                        axes[o,i].set_ylabel('')
            
            plt.savefig(ob_loc.joinpath(param+'.png'), dpi=600)
            plt.savefig(ob_loc.joinpath(param+'.svg'))
            plt.close()
    
    def dominant_alternative_bar(self, threshold):
        objnames=['\nPumping Energy','\nWater Quality','\nUrban Flooding']
        altnameslong=['Historical','Wastewater\nReuse','Infiltration\nBasins','Repair\nLeaks']
        
        data_clust_obj = pd.DataFrame()
                
        for c, cluster in enumerate(self.clusters):
            data_objmax = self.alldata[cluster].copy()
            
            for o, ob in enumerate(objnames):
                objscen = [self.objsymbol[o] + '-' + s for s in self.altnames]
                minscen = self.alldata[objscen]
                minscen.columns = altnameslong
                minscen = minscen.idxmin(axis="columns")
                minscen = minscen.to_frame(ob)
                data_objmax = pd.concat([data_objmax,minscen], axis=1)
        
            df = data_objmax.nsmallest(int(self.n_samples*threshold), cluster)
            df = pd.melt(df, id_vars=cluster, var_name='Objective', value_name='Alternative')
                
            df = df.drop(columns=[cluster])
            df['Cluster'] = '\n' + cluster
            
            data_clust_obj = pd.concat([data_clust_obj, df], axis=0)
        
        x = data_clust_obj.pivot_table(index=['Cluster','Alternative','Objective'], aggfunc='size')
        data = pd.DataFrame(np.array(list(map(list, x.index))), columns=['Cluster','Alternative','Objective'])
        data['Count'] = x.values
        g = sns.catplot(kind='bar', x='Cluster', y='Count', data=data, hue='Alternative', row='Objective', height=3, aspect=3, sharex=True, sharey=True, margin_titles=True, legend_out=True)
        
        g.tight_layout()
        g.savefig(Path.cwd().joinpath('images').joinpath('sa').joinpath('Bar_Dominant-Alternative.png'))
        plt.close()

    def param_sensitivity(self, sens_type='delta', cluster_index=[0,1,2,3,4,30], threshold='0.1', constant='Alternative', def_clust=5, def_alt=0, def_obj=0, def_params=0):
        
        if constant=='Historical':
            names = self.obj_alt.copy()
            names.extend(['A-Error'])
            objnames = ['\nEnergy','\nWater Quality','\nUrban Flooding','\nError Metric']
            altnameslong = ['Historical','Wastewater\nReuse','Infiltration\nBasins','Repair\nLeaks','Historical']
            objsymbol = ['E','W','F','A']
            altnames = ['Hist','WWTP','Basin','Leak','Error']
            if def_params==0:
                def_params = self.parnames
        else:
            names = self.obj_alt.copy()
            objsymbol = self.objsymbol
            altnames = self.altnames
            objnames = self.objnames
            altnameslong = self.altnameslong
        
        sensitivity_df = pd.DataFrame()
        for c, ci in enumerate(cluster_index):
            if constant=='Historical':
                df = pd.read_csv(Path.cwd().joinpath('SA_data').joinpath('sens').joinpath(self.filename).joinpath(sens_type+'-'+threshold+'-'+'{:02d}'.format(ci)+'.csv'),names=names)
            else:
                df = pd.read_csv(Path.cwd().joinpath('SA_data').joinpath('sens').joinpath(self.filename+'_original').joinpath(sens_type+'-'+threshold+'-'+'{:02d}'.format(ci)+'.csv'),names=names)
            df['Param'] = self.parnames
            x = df.rank()
            x = pd.melt(x, id_vars='Param', var_name='Obj-Alt', value_name='Rank')
            df = pd.melt(df, id_vars='Param', var_name='Obj-Alt', value_name='Sensitivity')
            df['Objective'] = df['Obj-Alt'].apply(lambda x: objnames[objsymbol.index(x[0])])
            df['Alternative'] = df['Obj-Alt'].apply(lambda x: altnameslong[altnames.index(x[2:])])
            df['Cluster'] = self.clusters[c]
            
            sensitivity_df = pd.concat([sensitivity_df,df], axis=0)
            
        if constant == 'Cluster':
            df = pd.DataFrame(sensitivity_df.loc[sensitivity_df['Cluster'].values==self.clusters[def_clust]].values, columns=['Param','Obj-Alt','Sensitivity','Objective','Alternative','Cluster'])
            g = sns.catplot(kind='bar', x='Param', y='Sensitivity', data=df, hue='Alternative', row='Objective', height=3, aspect=5, sharex=True, sharey=True, margin_titles=True, legend_out=True)
            plt.xticks(rotation=70)
            g.set(ylim=(0,0.3))
            
            g.tight_layout()
            g.fig.subplots_adjust(top=0.9)
            g.fig.suptitle('Parameter Sensitivity for Low Error Sample '+threshold+': Cluster '+self.clusters[def_clust])
            g.savefig(Path.cwd().joinpath('images').joinpath('sa').joinpath('bar-param-sensitivity').joinpath('cluster').joinpath('Bar_Param-Sensitivity_Cluster-'+self.clusters[def_clust]+'.png'))
            g.savefig(Path.cwd().joinpath('images').joinpath('sa').joinpath('bar-param-sensitivity').joinpath('cluster').joinpath('Bar_Param-Sensitivity_Cluster-'+self.clusters[def_clust]+'.svg'))
            plt.close()
            
        elif constant == 'Alternative':
            df = pd.DataFrame(sensitivity_df.loc[sensitivity_df['Alternative'].values==self.altnameslong[def_alt]].values, columns=['Param','Obj-Alt','Sensitivity','Objective','Alternative','Cluster'])
            g = sns.catplot(kind='bar', x='Param', y='Sensitivity', data=df, hue='Cluster', row='Objective', height=3, aspect=5, sharex=True, sharey=True, margin_titles=True, legend_out=True)
            plt.xticks(rotation=70)
            g.set(ylim=(0,0.3))
            
            g.tight_layout()
            g.fig.subplots_adjust(top=0.9)
            g.fig.suptitle('Parameter Sensitivity for Low Error Sample '+threshold+': '+self.altnameslong[def_alt]+' Alternative')
            g.savefig(Path.cwd().joinpath('images').joinpath('sa').joinpath('bar-param-sensitivity').joinpath('alternative').joinpath('Bar_Param-Sensitivity_Alternative-'+self.altnames[def_alt]+'.png'))
            g.savefig(Path.cwd().joinpath('images').joinpath('sa').joinpath('bar-param-sensitivity').joinpath('alternative').joinpath('Bar_Param-Sensitivity_Alternative-'+self.altnames[def_alt]+'.svg'))
            plt.close()
        
        elif constant == 'Objective':
            df = pd.DataFrame(sensitivity_df.loc[sensitivity_df['Objective'].values==objnames[def_obj]].values, columns=['Param','Obj-Alt','Sensitivity','Objective','Alternative','Cluster'])
            df['Cluster'] = '\n' + df['Cluster']
            
            g = sns.catplot(kind='bar', x='Param', y='Sensitivity', data=df, row='Cluster', hue='Alternative', height=3, aspect=5, sharex=True, sharey=True, margin_titles=True, legend_out=True, palette=sns.xkcd_palette(['white', 'light grey', 'medium grey', 'black']), edgecolor='k')
            plt.xticks(rotation=70)
            g.set(ylim=(0,0.3))
            
            g.tight_layout()
            g.fig.subplots_adjust(top=0.9)
            g.fig.suptitle('Parameter Sensitivity for Low Error Sample '+threshold+': '+objnames[def_obj]+' Objective')
            g.savefig(Path.cwd().joinpath('images').joinpath('sa').joinpath('bar-param-sensitivity').joinpath('objective').joinpath('Bar_Param-Sensitivity_Objective-'+self.objsymbol[def_obj]+'.png'))
            g.savefig(Path.cwd().joinpath('images').joinpath('sa').joinpath('bar-param-sensitivity').joinpath('objective').joinpath('Bar_Param-Sensitivity_Objective-'+self.objsymbol[def_obj]+'.svg'))
            plt.close()

        elif constant == 'Historical':
            
            for i in self.parnames:
                if not (i in def_params):
                    sensitivity_df = sensitivity_df[sensitivity_df['Param'] != i]
            df = pd.DataFrame(sensitivity_df.loc[sensitivity_df['Alternative'].values==altnameslong[def_alt]].values, columns=['Param','Obj-Alt','Sensitivity','Objective','Alternative','Cluster'])
            g = sns.catplot(kind='bar', x='Cluster', y='Sensitivity', data=df, col='Param', row='Objective', height=3, aspect=0.8, sharex=True, sharey=True, margin_titles=True, legend=True, row_order=['\nError Metric','\nEnergy','\nWater Quality','\nUrban Flooding'])
            plt.xticks(rotation=70)
            g.set(ylim=(0,0.31))
            g.set(xticklabels=[])
            g.set(xlabel=None)
            
            # Iterate thorugh each axis
            for ax in g.axes.flat:
            
                # Make title more human-readable and larger
                if ax.get_title():
                    ax.set_title(ax.get_title().split('=')[1])
            
                # Make right ylabel more human-readable and larger
                # Only the 2nd and 4th axes have something in ax.texts
                if ax.texts:
                    # This contains the right ylabel text
                    txt = ax.texts[0]
                    ax.text(txt.get_unitless_position()[0], txt.get_unitless_position()[1], txt.get_text().split('=')[1], transform=ax.transAxes, rotation=270, va='center')
                    # Remove the original text
                    ax.texts[0].remove()
            
            g.tight_layout()
            g.savefig(Path.cwd().joinpath('images').joinpath('sa').joinpath('bar-param-sensitivity').joinpath('Bar_Param-Sensitivity_Historical-Error-Object.png'))
            g.savefig(Path.cwd().joinpath('images').joinpath('sa').joinpath('bar-param-sensitivity').joinpath('Bar_Param-Sensitivity_Historical-Error-Object.svg'))
            plt.close()
            
    def obj_spread(self, threshold=0.1, max_val=1, step_val=0.05, spread=True, ranking=False):
        if not self.ran_cluster_kde:
            self.cluster_kde(threshold=threshold, plot=False)
            
        param = 'HK_4'
        o_dict = {}
        for o, ob in enumerate(self.objsymbol):
            o_dict[ob] = {}
            for j in range(0,3):
                o_dict[ob]['AbsDiff-'+str(j+1)] = pd.DataFrame()
                o_dict[ob]['PerDiff-'+str(j+1)] = pd.DataFrame()
            
            for c in self.clusters:
                df = pd.DataFrame(self.param_by_cluster.loc[self.param_by_cluster['Param'].values==param].values, columns=['Value','Cluster','Param'])
                df = pd.DataFrame(df.loc[df['Cluster'].values==c].values, columns=['Value','Cluster','Param'])
                for a in self.altnames:
                    df[a] = self.alldata[ob+'-'+a].values[self.thresh_samples[c][self.thresh_samples['Threshold']==self.threshold].values.astype('int')].copy()
                
                x = df[self.altnames].rank(axis=1,method='first')
                ranked_array = np.array(x)
                
                for j in range(0,3):
                    df['AbsDiff-'+str(j+1)] = np.nan
                    df['PerDiff-'+str(j+1)] = np.nan
                    
                    for row in df.itertuples():
                        df['AbsDiff-'+str(j+1)][row.Index] = row[int(np.arange(0,4)[ranked_array[row.Index]==1+(j+1)]+4)] - row[int(np.arange(0,4)[ranked_array[row.Index]==1]+4)]
                        df['PerDiff-'+str(j+1)][row.Index] = df['AbsDiff-'+str(j+1)][row.Index] / row[int(np.arange(0,4)[ranked_array[row.Index]==1]+4)]
                    
                    o_dict[ob]['AbsDiff-'+str(j+1)][c] = df['AbsDiff-'+str(j+1)].values
                    o_dict[ob]['PerDiff-'+str(j+1)][c] = df['PerDiff-'+str(j+1)].values * 100
                
                # Incorporate ranked information into pandas dataframe with identifying information
                x['Objective'] = ob
                x['Cluster'] = c
                
                self.ranked_df = pd.concat([self.ranked_df,x])
             
            if spread:    
                for j in ['1','2','3']:
                    plt.figure(figsize=(8,5))
                    sns.histplot(data=o_dict[ob]['PerDiff-'+j], element="poly", palette="bright", bins=list(np.arange(0,max_val,step_val)), alpha=0.07)
                    plt.title(ob+' PerDiff-'+j)
                    plt.savefig(Path.cwd().joinpath('images').joinpath('percent_diff').joinpath('max-'+str(int(max_val*100))+'_step-'+str(int(step_val*100))).joinpath(ob+'_PerDiff-'+j+'.png'))
                    plt.savefig(Path.cwd().joinpath('images').joinpath('percent_diff').joinpath('max-'+str(int(max_val*100))+'_step-'+str(int(step_val*100))).joinpath(ob+'_PerDiff-'+j+'.svg'))
                    plt.close()
        
        if ranking:
            z = self.ranked_df.melt(id_vars=['Objective','Cluster'], var_name='Alternative', value_name='Rank')
            z['Rank'] = z['Rank'].astype(int).astype(str)
            z = z.sort_values(by=['Rank'])
            
            sns.set_palette(palette="magma")
            
            g = sns.displot(z, x="Cluster", y="Rank", col="Alternative",  row="Objective", palette=sns.color_palette("viridis", as_cmap=True))
            g.set_titles('')
            g.set_xlabels('')
            g.set_xticklabels('')
            g.tight_layout()

    def heatmap_by_param(self,  params=-1, threshold=-1, n_samples=100000, c='C-12345'):
        # Heatmap - Params
        objnames = ['Pumping\nEnergy','Water\nQuality Risk','Urban\nFlooding']
        altnameslong = ['Historical','WW Reuse','Basins','Repair Leaks']
        cluster_loc = Path.cwd().joinpath('images').joinpath('heatmaps').joinpath('Ranking_Heatmap-params_' + c)
        cluster_loc.mkdir(exist_ok=True)
        
        if not self.ran_cluster_kde:
            self.cluster_kde(params=params, threshold=threshold, n_samples=n_samples, plot=False)
        
        # Create list of labels for figures
        if type(params) is int:
            params = self.parnames
            paramlabels = self.parnameslong
        else:
            paramlabels = []
            for p in params:
                paramlabels.extend([self.parnameslong[self.parnames.index(p)]])
        
        for p, param in enumerate(params):
            print(param)
            
            o_dict = {}
            ranked_df = pd.DataFrame()
            
            for o, ob in enumerate(self.objsymbol):
                o_dict[ob] = {}
                df = pd.DataFrame(self.param_by_cluster.loc[self.param_by_cluster['Param'].values==param].values, columns=['Value','Cluster','Param'])
                df = pd.DataFrame(df.loc[df['Cluster'].values==c].values, columns=['Value','Cluster','Param'])
                
                for a, alt in enumerate(self.altnames):
                    df[altnameslong[a]] = self.alldata[ob+'-'+alt].values[self.thresh_samples[c][self.thresh_samples['Threshold']==threshold].values.astype('int')].copy()
                
                # Rank each alternative within each parameter set 
                x = df[altnameslong].rank(axis=1,method='first')
                
                # Incorporate ranked information into pandas dataframe with identifying information, assign each datapoint to a decimal bin within the parameter range
                x['Objective'] = objnames[o]
                temp = (df['Value']*10).astype(float)
                x[param] = (np.floor(temp)/10).astype(str)
                
                ranked_df = pd.concat([ranked_df,x])
            
            # Reshape dataset to allow for histogram
            ranked_df_melt = ranked_df.melt(id_vars=['Objective',param], var_name='Alternative', value_name='Rank')
            ranked_df_melt['Rank'] = ranked_df_melt['Rank'].astype(int).astype(str)
            ranked_df_melt = ranked_df_melt.sort_values(by=['Rank'])
            
            # Pivot table to create histogram
            ranked_df_table = pd.pivot_table(ranked_df_melt, values='Objective', index=['Objective', 'Alternative', 'Rank'], columns=[param], aggfunc=len)
            
            # Fill any missing regions in the parameter range with NAN
            for column in np.arange(0,1,0.1):
                try:
                    temp = ranked_df_table[str(np.round(column,1))]
                except:
                    ranked_df_table[str(np.round(column,1))] = np.nan
            
            # Convert from count to percentage of bin total
            for o in objnames:
                temp = ranked_df_table.loc[o]
                for column in temp:
                    temp[column] = temp[column]/(temp[column].sum()/4)
            
            # Fill any missing values in alternative rank space with NAN
            for o in objnames:
                for a in altnameslong:
                    for r in range(1,5):
                        try:
                            temp = ranked_df_table.loc[o,a,str(r)]
                        except:
                            ranked_df_table.loc[o,a,str(r)] = np.nan
            
            # Reshape dataset to use seaborn facetgrid
            ranked_df_normalized = ranked_df_table.reset_index().melt(id_vars=['Objective', 'Alternative', 'Rank'], value_name='Percent Sample')
            
            # Define heatmap function that can be used with seaborn facetgrid
            def draw_heatmap(*args, **kwargs):
                data = kwargs.pop('data')
                d = data.pivot(index=args[1], columns=args[0], values=args[2])
                sns.heatmap(d, **kwargs)
            
            g = sns.FacetGrid(ranked_df_normalized, col='Alternative', row='Objective', height=3, aspect=1.2, margin_titles=True, col_order=altnameslong, row_order=objnames)
            g.map_dataframe(draw_heatmap, param, 'Rank', 'Percent Sample', cbar=False, vmin=0, vmax=1, cmap=sns.color_palette("viridis", as_cmap=True), square=False)
            g.set_titles(col_template="{col_name}", row_template="{row_name}", size='large')
            g.set_xticklabels(['min','','','','','','','','','max'], rotation=0, size='medium')
            g.tight_layout()
            
            plt.savefig(cluster_loc.joinpath(param+'_ranking-heatmap_'+c+'.png'), dpi=600)
            plt.savefig(cluster_loc.joinpath(param+'_ranking-heatmap_'+c+'.svg'))
            plt.close()

#            return ranked_df_melt, ranked_df_table, ranked_df_normalized
            
## Heatmap - Clusters
## Not yet implemented as a function. The following code was used to make the cluster-based heatmaps.
#param = 'HK_2'
#ranked_df = pd.DataFrame()
#objnameslong = ['Pumping\nEnergy','Water\nQuality Risk','Urban\nFlooding']
#altnameslong = ['Historical','WW Reuse','Basins','Repair Leaks']
#
#for o, ob in enumerate(objsymbol):
#    print(ob)
#    
#    for c in clusters:
#        print(c)
#        df = pd.DataFrame(param_by_cluster.loc[param_by_cluster['Param'].values==param].values, columns=['Value','Cluster','Param'])
#        df = pd.DataFrame(df.loc[df['Cluster'].values==c].values, columns=['Value','Cluster','Param'])
#        for a in altnames:
#            df[altnameslong[altnames.index(a)]] = alldata[ob+'-'+a].values[thresh_samples[c][thresh_samples['Threshold']==threshold].values.astype('int')].copy()
#        
#        # Convert 
#        x = df[altnameslong].rank(axis=1,method='first')
#        ranked_array = np.array(x)
#        
#        # Incorporate ranked information into pandas dataframe with identifying information
#        x['Objective'] = objnameslong[objsymbol.index(ob)]
#        x['\nCluster'] = c
#        
#        ranked_df = pd.concat([ranked_df,x])
#
#z = ranked_df.melt(id_vars=['Objective','\nCluster'], var_name='Alternative', value_name='Rank')
#z['Rank'] = z['Rank'].astype(int).astype(str)
#z = z.sort_values(by=['Rank'])
#
#g = sns.displot(z, x='\nCluster', y='Rank', col='Alternative', row='Objective', height=3, aspect=1, cmap=sns.color_palette("viridis", as_cmap=True), cbar=False, vmin=0, vmax=5000, facet_kws=dict(row_order=objnameslong, col_order=altnameslong, margin_titles=True))
#g.set_titles(col_template="{col_name}", row_template="{row_name}", size='large')
#g.set_xticklabels(clusters,rotation=90,size='medium')
#
#g.tight_layout()

##%% Box Plots
## This plotting was explored but not implemented for publication
## Objective
#names = [allnames[0]]
#names.extend(objnames[0:4])
#obj_data = all_data[names].melt('Index', var_name='cols', value_name='vals')
#f, ax = plt.subplots(figsize=(7, 6))
#ax.set_yscale("log")
#g = sns.boxplot(x="cols", y="vals", data=obj_data.sort_values("cols"))
##
### Objective with groups by full error
##names = [allnames[0]]
##names.extend(objnames[0:4])
##obj_data = all_data[names].melt('Index', var_name='cols', value_name='vals')
##obj_data['Low_Error'] = np.tile(np.where(all_data['Low_Error'] == True, True, False),len(objnames[0:4]))
##f, ax = plt.subplots(figsize=(7, 6))
##ax.set_yscale("log")
##g = sns.boxplot(x="cols", y="vals", hue='Low_Error', data=obj_data.sort_values("cols"))
#
## Error by cluster
#names = [allnames[0]]
#names.extend(self.clusters)
#err_data = all_data[names].melt('Index', var_name='cols', value_name='vals')
#f, ax = plt.subplots(figsize=(7, 6))
#ax.set_yscale("log")
#g = sns.boxplot(x="cols", y="vals", data=err_data.sort_values("cols"))
#plt.setp(ax.get_xticklabels(), rotation=90)
#plt.title('Error')
##
### Error by cluster with groups by full error
##names = [allnames[0]]
##names.extend(self.clusters)
##err_data = all_data[names].melt('Index', var_name='cols', value_name='vals')
##err_data['Low_Error'] = np.tile(np.where(all_data['Low_Error'] == True, True, False),len(self.clusters))
##f, ax = plt.subplots(figsize=(7, 6))
##ax.set_yscale("log")
##g = sns.boxplot(x="cols", y="vals", hue='Low_Error', data=err_data.sort_values("cols"))
#