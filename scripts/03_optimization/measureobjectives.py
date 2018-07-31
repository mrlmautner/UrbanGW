# -*- coding: utf-8 -*-
"""
Created on Thu May 24 18:59:47 2018

@author: MM
"""

import os
import flopy
import numpy as np
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf
import time
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import calendar
import pandas as pd

xll = 455000
yll = 2107000
xur = 539000
yur = 2175000
cellsize = 500

ncol = int((xur-xll)/cellsize) # Number of rows
nrow = int((yur-yll)/cellsize) # Number of columns

os.chdir(r'C:\Users\MM\Google Drive\Tlalpan\Simulation_Wrapper\Results')
# Assign name and create modflow model object
modelname = 'VM_Basin'
#%% Head Dictionary
S_heads = {}
for s_name in ['Historical','WWTP','Leak','Basin']:
    S_heads[s_name] = bf.HeadFile('VM_'+s_name+'.hds')
    
#%% Drawdown dictionary
S_drawdown = {}
for s_name in ['Historical','WWTP','Leak','Basin']:
    S_drawdown[s_name] = bf.HeadFile('VM_'+s_name+'.ddn', text='drawdown')
#hds = bf.HeadFile(modelname+'.hds')

#%% Drawdown Contour
fig, axes = plt.subplots(2, 2, figsize=(7,6.3))
a=0
mapTitle = ['Historical','Increased WW Reuse','Repair Leaks','Recharge Basins']
plt.set_cmap('coolwarm')
axes = axes.flat

for s_name in ['Historical','WWTP','Leak','Basin']:
    ddn = S_drawdown[s_name]
    d = ddn.get_data(mflay=[0,1],kstpkper=(30,359))
    
    DNEW = np.ones(d.shape)*np.nan#min(h[h>0])
    for i, ddown in enumerate(d):
        DNEW[i][ddown>-900] = ddown[ddown>-900]
    
    for i, ddnlayer in enumerate(DNEW):
        im = axes[a].imshow(ddnlayer, vmin=-10, vmax=10)
        axes[a].set_title(mapTitle[a].format(i+1))
        ctr = axes[a].contour(ddnlayer, colors='k', linewidths=0.5)
    
    a+=1
    
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
fig.colorbar(im, cax=cbar_ax, label='Groundwater Drawdown (m)')
plt.savefig('DDown_Contour_All.svg')
plt.savefig('DDown_Contour_All.png')
plt.show()

#%% Heads Contour
fig, axes = plt.subplots(2, 2, figsize=(7,6.3))
a=0
mapTitle = ['Historical','Increased WW Reuse','Repair Leaks','Recharge Basins']
plt.set_cmap('coolwarm')
axes = axes.flat

for s_name in ['Historical','WWTP','Leak','Basin']:
    hds = S_heads[s_name]
    h = hds.get_data(mflay=[0,1],kstpkper=(30,359))
    
    HNEW = np.ones(h.shape)*np.nan#min(h[h>0])
    for i, head in enumerate(h):
        HNEW[i][head>-900] = head[head>-900]
    
    for i, hdslayer in enumerate(HNEW):
        im = axes[a].imshow(hdslayer, vmin=2100, vmax=2400)
        axes[a].set_title(mapTitle[a].format(i+1))
        ctr = axes[a].contour(hdslayer, colors='k', linewidths=0.5)
    
    a+=1
    
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
fig.colorbar(im, cax=cbar_ax, label='Groundwater Drawdown (m)')
plt.savefig('HDS_Contour_All.svg')
plt.savefig('HDS_Contour_All.png')
plt.show()

#%% Budget
modelname = 'Leak'
mf_list = flopy.utils.MfListBudget('VM_'+modelname+".list")
incremental, cumulative = mf_list.get_budget()

#%%
df_1Bdget, df_extra = mf_list.get_dataframes(start_datetime="04-30-1984")

mthly_Bdget = df_1Bdget.drop(['CONSTANT_HEAD_IN','TOTAL_IN','CONSTANT_HEAD_OUT','RECHARGE_OUT','TOTAL_OUT','IN-OUT','PERCENT_DISCREPANCY'], axis=1)

mthly_Bdget['STORAGE_OUT'] = mthly_Bdget['STORAGE_OUT'].apply(lambda x: x*-1)
mthly_Bdget['WELLS_OUT'] = mthly_Bdget['WELLS_OUT'].apply(lambda x: x*-1)
mthly_Bdget = mthly_Bdget.multiply(30/1000000)
cols = mthly_Bdget.columns.tolist()
# reorder columns
cols = [cols[1]] + [cols[2]] + [cols[0]] + [cols[4]] + [cols[3]] 
# "commit" the reordering
mthly_Bdget = mthly_Bdget[cols]

ax = mthly_Bdget['01-31-1985':'12-31-2013'].plot.area(stacked=True,figsize=(8,9),color=['blue','xkcd:azure','lightblue','red','lightpink'])
plt.ylabel(r'Volume ($hm^3$)')
#plt.ylim(-8,10)
plt.title('Groundwater Budget')
plt.legend(['Leaks','Precipitation','Storage: In','Pumping','Storage: Out'],loc=4)

plt.savefig('WB_'+modelname+'.png')
plt.savefig('WB_'+modelname+'.svg')
plt.show()


#%% Cumulative overdraft
l = [4,1,1,1]
c = ['lightblue','r','g','k']
i = 0

#df_Bdget = {}
#df_extra = {}
#incremental = {}
#cumulative = {}

for s_name in ['Historical','WWTP','Leak','Basin']:
    
#    mf_list= flopy.utils.MfListBudget('VM_'+s_name+".list")
#    incremental[s_name], cumulative[s_name] = mf_list.get_budget()
#    
#    df_Bdget[s_name], df_extra[s_name] = mf_list.get_dataframes(start_datetime="01-31-1984")
    df_extra[s_name]['IN'] = df_extra[s_name]['RECHARGE_IN'].divide(1000000) + df_extra[s_name]['WELLS_IN'].divide(1000000)
    df_extra[s_name]['OUT'] = df_extra[s_name]['WELLS_OUT'].divide(1000000)
    df_extra[s_name]['INOUTCUMSUM'] = df_extra[s_name]['IN'] - df_extra[s_name]['OUT']
    
    df_extra[s_name].INOUTCUMSUM['01-31-1985':'12-31-2013'].plot(linewidth=l[i],color=c[i])
    i+=1
    
plt.ylabel(r'Volume ($hm^3$)')
#plt.ylim(-8,10)
#plt.title('Cumulative In - Out')
plt.legend(['Historical','Increase WW Reuse','Repair Leaks','Recharge Basins'])

plt.savefig('INOUTCumSum.png')
plt.savefig('INOUTCumSum.svg')
plt.show()
#%% Time series by location
t = pd.DatetimeIndex(freq='M',start='04/30/1984',end='12/31/2013')
coords = [90,25] # mtn, [67,97] # pump3, [52,61] # subs, [84,65] # pump2, [59,54] # pump1, 

i = 0
dTimeS = np.zeros((357,4))

l = [4,1,1,1]
mark = ['-','-^','-*','-s']
c = ['lightblue','r','g','k']
z = [0,1,1,1]

for s_name in ['Historical','WWTP','Leak','Basin']: #]:
    j = 0
    for m in range(3,12):
        d = S_drawdown[s_name].get_data(kstpkper=((calendar.monthrange(1984,m+1)[1]-1),m),mflay=1)
        dTimeS[j,i] = d[coords[0],coords[1]]
        j+=1
    for y in range(1,30):
        for m in range(0,12):
            d = S_drawdown[s_name].get_data(kstpkper=((calendar.monthrange(1984+y,m+1)[1]-1),y*12+m),mflay=1)
            dTimeS[j,i] = d[coords[0],coords[1]]
            j+=1
    
    plt.plot(t, dTimeS[:,i],linewidth=l[i],color=c[i])#,mark[i],zorder=z[i],markersize=5)
    i+=1
    
plt.xlabel('Year')
plt.ylabel('Drawdown (m)')
plt.ylim(plt.ylim()[::-1])
plt.title('Drawdown at mountain region location')
plt.legend(['Historical','Increased WW Reuse','Repair Leaks','Recharge Basins'])
plt.savefig('DDN_Indicator_Mtn.svg')
plt.show()

#%%
t = ['January','February','March','April','May','June','July','August','September','October','November','December']
i = 0
width = 0.4
AM = np.zeros((12,2))
for s in ['9022-EL GUARDA-DF-Precip.csv','15167-EL TEJOCOTE-MEX-Precip.csv']:
    fileName = 'F:/Tlalpan Model/Data/Precipitation/Raw/' + s
    df = pd.read_csv(fileName, delimiter=',', header=None, names=['Month', 'Day', 'Year', 'Precip'], skiprows=8)
    df.Precip = pd.to_numeric(df.Precip,'coerce')
    df = df.set_index(pd.DatetimeIndex(pd.to_datetime((df.Year*10000+df.Month*100+df.Day).apply(str),format='%Y%m%d')))
    temp = df.Precip.resample('M').sum()
    
    # Remove any months with more than 10 days total of missing data
    nanCount = df.Precip.fillna(0.0).resample('M').count()
    noNanCount = df.Precip.resample('M').count()
    x = nanCount - noNanCount
    y = x.index[x > 10]
    temp[y] = np.nan
    
    # Remove any years with more than 3 days of missing data in the wet season
    y = x.index[x > 3]
    y = y[y.month.isin([6,7,8,9])]
    temp[y] = np.nan
    
    avgMthly = np.zeros(12)
    for m in range(0,12):
        mthly = np.ones(30)*np.nan
        for y in range(1984,2014):
            mthly[(y-1984)] = temp[str(y)+'-'+format(m+1, '02d')+'-'+str(calendar.monthrange(y,m+1)[1])]
        
        avgMthly[m] = np.nanmean(mthly)
        
    AM[:,i] = avgMthly
    i+=1

fig = plt.figure(figsize=(6,7))
ax = fig.add_subplot(111)
ax.bar(np.arange(0,12)-width,AM[:,0],width,color='navy',label='-Ymin',align='edge')
ax.bar(np.arange(0,12),AM[:,1],width,color='orange',label='Ymax',align='edge')
plt.xticks(np.arange(0,12), t,rotation=45)
ax.get_yaxis().set_tick_params(which='both', direction='in')
ax.grid(linestyle='--', linewidth=0.5)
plt.legend(['Station 9022', 'Station 15167'])
plt.title('Precipitation at Meteorological Stations')
plt.savefig('Precip.png')
plt.close()

#%% Time series by location
#t = pd.DatetimeIndex(freq='M',start='04/30/1984',end='12/31/2013')
#coords = [90,25] # mtn, [29,77] # pump, [52,61] # subs, [59,54]
#
#i = 0
#hTimeS = np.zeros((357,4))
#
#l = [4,1,1,1]
#mark = ['-','-^','-*','-s']
#c = ['lightblue','r','g','k']
#z = [0,1,1,1]
#
#for s_name in ['Historical','WWTP','Leak','Basin']: #]:
#    j = 0
#    for m in range(3,12):
#        h = S_heads[s_name].get_data(kstpkper=((calendar.monthrange(1984,m+1)[1]-1),m),mflay=1)
#        hTimeS[j,i] = h[coords[0],coords[1]]
#        j+=1
#    for y in range(1,30):
#        for m in range(0,12):
#            h = S_heads[s_name].get_data(kstpkper=((calendar.monthrange(1984+y,m+1)[1]-1),y*12+m),mflay=1)
#            hTimeS[j,i] = d[coords[0],coords[1]]
#            j+=1
#    
#    plt.plot(t, hTimeS[:,i],linewidth=l[i],color=c[i])#,mark[i],zorder=z[i],markersize=5)
#    i+=1
#    
#plt.xlabel('Year')
#plt.ylabel('Head Elevation')
#plt.title('Head elevation over time at mountainous location')
#plt.legend(['Historical','Increased WW Reuse','Repair Leaks','Recharge Basins'])
#plt.savefig('HDS_Indicator_Mountain.svg')
#plt.close()

#%% Time series by location
#subsTimeS = np.zeros(30)
#t = range(1984,2014)
#for year in range(0,30):
#    per = year*12+2
#    h = hds.get_data(kstpkper=(30,per),mflay=1)
#    subsTimeS[year] = h[59,54]
#
#plt.plot(t, subsTimeS)
#
#plt.xlabel('Year')
#plt.ylabel('Head Elevation')
#plt.title('Head elevation over time at high subsidence location')
#plt.savefig('Indicator_sub_'+modelname+'.png')
#plt.close()
#
#pumpTimeS = np.zeros(30)
#t = range(1984,2014)
#for year in range(0,30):
#    per = year*12+7
#    h = hds.get_data(kstpkper=(30,per),mflay=1)
#    pumpTimeS[year] = h[29,77]
#
#plt.plot(t, pumpTimeS)
#
#plt.xlabel('Year')
#plt.ylabel('Head Elevation')
#plt.title('Head elevation over time at high pumping location')
#plt.savefig('Indicator_pump_'+modelname+'.png')
#plt.close()
#
#
#
##%% Drawdown in mountains
#t = pd.DatetimeIndex(freq='M',start='02/29/1984',end='12/31/2013')
##t = range(1984,2014)
#i = 0
#mtnTimeS = np.zeros((359,4))
#
#l = [4,0.5,0.5,0.5]
#mark = ['-','-^','-*','-s']
#c = ['lightblue','k','k','k']
#z = [0,1,1,1]
#
#for s_name in ['Historical']:#,'WWTP','Leak','Basin']:
#    j = 0
#    for m in range(3,12):
#        d = S_drawdown[s_name].get_data(kstpkper=((calendar.monthrange(1984,m+1)[1]-1),m),mflay=1)
#        mtnTimeS[j,i] = d[90,27]
#        j+=1
#    for y in range(1,30):
#        for m in range(0,12):
#            d = S_drawdown[s_name].get_data(kstpkper=((calendar.monthrange(1984+y,m+1)[1]-1),y*12+m),mflay=1)
#            mtnTimeS[j,i] = d[90,27]
#            j+=1
#    
#    plt.plot(t, mtnTimeS[:,i],mark[i],color=c[i],zorder=z[i],linewidth=l[i],markersize=5)
#    i+=1
#    
#plt.xlabel('Year')
#plt.ylabel('Head Elevation')
#plt.legend(['Historical'])#,'Increased WW Reuse','Repair Leaks','Recharge Basins'])
#plt.title('Drawdown over time at mountainous location')
#plt.savefig('DDN_Mountainous_Historical.png')
#plt.close()

#%% Head Elevation Contour
#for s_name in ['TR','WWTP','Leak','Basin']:
#    hds = S_heads[s_name]
#    h = hds.get_data(mflay=[0,1],kstpkper=(30,246))
#    
#    fig, axes = plt.subplots(2, 1, figsize=(8, 8))
#    plt.set_cmap('coolwarm_r')
#    axes = axes.flat
#    HNEW = np.ones(h.shape)*np.nan#min(h[h>0])
#    for i, heads in enumerate(h):
#        HNEW[i][heads>0] = heads[heads>0]
#    for i, hdslayer in enumerate(HNEW):
#        im = axes[i].imshow(hdslayer, vmin=2100, vmax=2300)
#        axes[i].set_title('Layer {}'.format(i+1))
#        ctr = axes[i].contour(hdslayer, colors='k', linewidths=0.5)
#    #
#    #fig.delaxes(axes[-1])
#    fig.subplots_adjust(right=0.8)
#    cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
#    fig.colorbar(im, cax=cbar_ax, label='Head (masl)')
#    
#    plt.savefig('Heads_'+s_name+'.svg')
#    plt.close()

#%%
fig = plt.figure()
ax = fig.gca(projection='3d')

#for s_name in ['Historical','WWTP','Leak','Basin']:
hds = S_heads['Historical']
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

#%%
mf = flopy.modflow.Modflow.load('VM_Historical.nam')
# Create the headfile and budget file objects
headobj = bf.HeadFile('VM_Historical.hds')
times = headobj.get_times()
cbb = bf.CellBudgetFile('VM_Historical.cbc')

#%%
time = times[701]
# Make the plots
head = headobj.get_data(totim=time)
#Print statistics
print('Head statistics')
print('  min: ', head.min())
print('  max: ', head.max())
print('  std: ', head.std())

# Extract flow right face and flow front face
frf = cbb.get_data(text='FLOW RIGHT FACE', totim=time)[0]
fff = cbb.get_data(text='FLOW FRONT FACE', totim=time)[0]
flf = cbb.get_data(text='FLOW LOWER FACE', totim=time)[0]

#Create the plot
f = plt.figure()
ax = f.gca(projection='3d')
#
#modelmap = flopy.plot.ModelMap(model=mf, layer=1)
#qm = modelmap.plot_ibound()
#cs = modelmap.contour_array(head)
#norm = plt.colors.Normalize()
#norm.autoscale(Z)
#sm = cm.ScalarMappable(cmap=cm.coolwarm, norm=norm)
#sm.set_array([])
#ax.quiver(X,Y,Z,frf, fff, flf,length=100,normalize=True,color=cm(norm(o)))
q = ax.quiver(X,Y,Z,frf[1], fff[1], flf[1],length=0.05, cmap=cm.Blues)
q.set_array(np.sort(np.reshape(flf[1],136*168)))

plt.show()

#%%
sat_thck = np.ones((3,136,168))
sat_thck[0,:,:] = np.ones((136,168))*50
sat_thck[1,:,:] = np.ones((136,168))*350
sat_thck[2,:,:] = np.ones((136,168))*1500
qx,qy,qz = flopy.plot.plotutil.centered_specific_discharge(frf, fff, flf, np.ones((168))*cellsize, np.ones((136))*cellsize, sat_thck)

f = plt.figure()
ax = f.gca(projection='3d')

q = ax.quiver(X,Y,Z,qx[1], qy[1], qz[1],length=10000, cmap=cm.coolwarm)
q.set_array(np.sort(np.reshape(qx[1],136*168)))

plt.show()