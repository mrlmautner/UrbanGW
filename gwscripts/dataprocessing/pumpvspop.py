# -*- coding: utf-8 -*-
"""
Created on Sat Mar 17 19:37:14 2018

@author: Marina Mautner
"""
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

df = pd.read_csv('C:/Users/MM/Google Drive/ValleMexico/Data/Pumping/PUMP_POP_INEGI_REPDA_20180320.csv')
df = df.set_index(pd.to_datetime(df.Year.apply(lambda x: str(x)+'/01/01')))

PSeries = df.PopZMVM_Census.dropna()
#PSeries = PSeries.append(df.PopZMVM_Proj.dropna())
t = df.Year[PSeries.index].values - 1895
P = PSeries.values

t1 = np.array([1960,1970,1980,1990,1995,2000,2005,2010])
p1 = df.PopZMVM_Census.dropna()['1960/01/01':]
t2 = np.array([2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,
               2023,2024,2025,2026,2027,2028,2029,2030])
p2 = df.PopZMVM_Proj.dropna()

# Fit population curve
# Census
g1 = np.polyfit(t1, p1, 2)
f1 = np.poly1d(g1)

# CONAPO
g2 = np.polyfit(t2, p2, 2)
f2 = np.poly1d(g2)

x = np.arange(1960,2030,1)
Pt1 = np.zeros(len(x))
Pt2 = np.zeros(len(x))
for i in range(0,len(x)):
    Pt1[i] = f1(x[i])
    Pt2[i] = f2(x[i])
    
plt.plot(np.array(df.Year),np.array(df.PopZMVM_Proj),'.',color='black')
plt.plot(t1,p1,'s',color='grey',zorder=3)
plt.plot(x,Pt1,'--',color='red')
plt.plot(x,Pt2,'--',color='blue')
plt.grid(True)
plt.xlabel('Year')
plt.ylabel('Population')
plt.legend(['CONAPO 2014 Projections','1895 to 2010 Census Data', 'Census Model','CONAPO Projection Model'])

#%%
PM = pd.DataFrame(np.transpose(np.array([x,Pt1])),columns=['Year','P'])
PM = PM.set_index(pd.to_datetime(PM.Year.apply(lambda x: str(int(x))+'/01/01')))
PM.to_csv('C:/Users/MM/Google Drive/ValleMexico/Data/Pumping/POP_CM.csv')

Pop = PM.P['2005-01-01':'2016-01-01']
Pop['2011-01-01'] = np.nan
Pop = Pop.dropna().values
Q = df.REPDA_m3d.dropna().values

WEL = pd.DataFrame(np.transpose(np.array([Q,Pop])/10**6),index=df.REPDA_m3d.dropna().index,columns=['Q','P'])
#
##Population
#g = np.polyfit(Pop, Q, 1)
#f = np.poly1d(g)

#Time
year = np.arange(2005,2017)
year = np.delete(year,6)
g = np.polyfit(year, Q, 1)
f = np.poly1d(g)

#residuals = f(Pop) - Q
residuals = f(year) - Q
ss_res = np.sum(residuals**2)
ss_tot = np.sum((Q-np.mean(Q))**2)
r_squared = 1 - (ss_res / ss_tot)
#rmse = np.sqrt(ss_res/len(Pop))
rmse = np.sqrt(ss_res/len(year))

plt.grid(True)
plt.xlim(1980,2020)
#sns.regplot(x=Pop, y='Q', data=WEL,marker='.',order=1)
sns.regplot(x=year, y='Q', data=WEL,marker='.',order=1)
#plt.xlabel('Population (millions)')
plt.xlabel('Year')
plt.ylabel('Groundwater Concessions (' r'$hm^3/d$' ')')
print(rmse)
print(r_squared)
print(f)

