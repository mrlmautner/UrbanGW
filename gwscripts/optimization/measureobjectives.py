# -*- coding: utf-8 -*-
"""
Created on Thu May 24 18:59:47 2018

@author: MM
"""
import numpy as np
import pandas as pd
import math
import calendar

def measureEnergy(heads,supply_dict,dem):
    E = 0
    
    efficiency = 0.7
    MJc = 9.81 # MegaJoules to lift 1 MegaLiter 1 meter
    kWhc = 3.6 # MegaJoules per kWh
    
    pd.DatetimeIndex(freq='M',start='01/31/1986',end='12/31/2013')
    coords = np.zeros(2)
    
    for i, p in supply_dict.items():
        # Start in year 2
        if i > 23:
            d = calendar.monthrange(1984+math.floor(i/12),(i%12)+1)[1] # number of days in stress period
            h = heads.get_data(kstpkper=(8,i),mflay=1) # heads binary file for last time step in stress period
            
            # loop through all well values
            for n in p:
                # Only measure energy use for pumping wells (negative flow)
                if n[3] < 0:
                    coords[0] = n[1]
                    coords[1] = n[2]
                    wellHead = h[int(coords[0]),int(coords[1])]
                    surface = dem[int(coords[0]),int(coords[1])]
                    depthtowater = max(surface-wellHead,0)
                    pumping = n[3] * (-1) * 0.001 # m3/d --> ML/d
                    E += (efficiency * depthtowater * pumping * MJc / kWhc)*d # kWh per month (stress period) of pumping
    
    return E

def measureSubidence(heads,dem,active,bottom):
    sub = np.zeros(360) # total difference between dem and h
    h = list(sub.copy()) # head map for each stress period
    cells = np.sum(np.sum(active[0])) # Number of cells in Clay layer
    
    # loop through all stress periods
    for t in range(0,360):
        # check during the reference month of the last year (March)
        h1 = heads.get_data(kstpkper=(8,t),mflay=0) # heads binary file for dth day of stress period
        h2 = heads.get_data(kstpkper=(8,t),mflay=1) # heads binary file for dth day of stress period
        
        h[t] = np.multiply(np.less(np.multiply(active[0],h2),0),bottom[1]) + np.multiply(np.less(np.multiply(active[0],h1),0).astype(float) - np.less(np.multiply(active[0],h2),0).astype(float),h2) + np.multiply(np.greater(np.multiply(active[0],h1),0),h1)
        
        sub_temp = np.multiply(active[0],dem) - h[t] # Only add subsidence measure if the hydraulic head is less than the ground surface
        sub_temp[sub_temp < 0] = 0
        
        sub[t] = np.sum(np.sum(sub_temp))

    return np.average(sub[23:]),cells,sub

def measureMound(heads,dem,active,LU,PhasePer):
    mound = 0 # cumulative head above DEM during model period
    urban_cells = 0
    mound_cells = 0
    
    for t in range(23,360):
#        d = calendar.monthrange(1984+math.floor(t/12),(t%12)+1)[1] # number of days in stress period
        h1 = heads.get_data(kstpkper=(8,t),mflay=0) # heads binary file for last time step of stress period
        h2 = heads.get_data(kstpkper=(8,t),mflay=1) # heads binary file for last time step of stress period
      
        if t < PhasePer[0]:
            LUset = '1990'
        elif t < PhasePer[1]:
            LUset = '2000'
        else:
            LUset = '2010'
        
        lu_temp = LU[LUset]['ARRAY']['URBAN']
        h = np.multiply(active[0],h1) + np.multiply((active[1]-active[0]),h2)
        above_dem = h > dem
        lu_active = np.sum(np.sum(np.multiply(lu_temp,active[1]))) # Total urban area of active cells
        lu_above = np.sum(np.sum(lu_temp[above_dem])) # Total area of cells with mounding above dem
        
        urban_cells += lu_active
        mound_cells += lu_above
    
    return mound,mound_cells,urban_cells

def get_objectives(heads, wellinfo, landuse, dem, active, bottom):

    energy = measureEnergy(heads, wellinfo, dem)
    sub, sub_cells, sub_monthly = measureSubidence(heads, dem, active, bottom)
    mound, mound_cells, urban_cells = measureMound(heads, dem, active, landuse, [132,252])
    
    return energy, sub/sub_cells, mound_cells/urban_cells
