# -*- coding: utf-8 -*-
"""
Created on Thu May 24 18:59:47 2018

@author: MM
"""
import numpy as np
import pandas as pd
import math
import calendar

def measureEnergy(heads,supply_dict,dem,bottom):
    E = 0
    
    efficiency = 0.7
    MJc = 9.81 # MegaJoules to lift 1 MegaLiter 1 meter
    kWhc = 3.6 # MegaJoules per kWh
    
    pd.DatetimeIndex(freq='M',start='01/31/1986',end='12/31/2013')
    coords = np.zeros(2)
    
    for i, p in supply_dict.items():
        # Start in year 2
        if i > 347:
            d = calendar.monthrange(1984+math.floor(i/12),(i%12)+1)[1] # number of days in stress period
            h = heads.get_data(kstpkper=(8,i),mflay=1) # heads binary file for last time step in stress period
            
            # loop through all well values
            for n in p:
                # Only measure energy use for pumping wells (negative flow)
                if n[3] < 0:
                    coords[0] = n[1]
                    coords[1] = n[2]
                    wellHead = max(h[int(coords[0]),int(coords[1])], bottom[1][int(coords[0]),int(coords[1])])
                    surface = dem[int(coords[0]),int(coords[1])]
                    depthtowater = max(surface-wellHead,1)
                    pumping = n[3] * (-1) * 0.001 # m3/d --> ML/d
                    E += (efficiency * depthtowater * pumping * MJc / kWhc)*d # kWh per month (stress period) of pumping
    
    return E

def measureSubidence(heads,dem,active,bottom):
    sub = np.zeros(360) # total difference between dem and h
    h = list(sub.copy()) # head map for each stress period
    cells_all = np.sum(np.sum(active[0])) # Number of cells in Clay layer
    cells_below = np.zeros(360) # Number of cells with head below ground surface
    
    # loop through all stress periods
    for t in range(0,360):
        # check during the last time step of each stress period
        h1 = heads.get_data(kstpkper=(8,t),mflay=0) # heads binary file layer 1
        h2 = heads.get_data(kstpkper=(8,t),mflay=1) # heads binary file layer 2
        
        h[t] = np.multiply(np.less(np.multiply(active[0],h2),0),bottom[1]) + np.multiply(np.less(np.multiply(active[0],h1),0).astype(float) - np.less(np.multiply(active[0],h2),0).astype(float),h2) + np.multiply(np.greater(np.multiply(active[0],h1),0),h1) # Matrix that has bottom of layer 2 if cell is dry, hydraulic head in layer 2 if layer 1 is dry less dry cells in layer 2, and hydraulic head in layer 1 if higher than bottom of layer 1
        
        sub_temp = np.multiply(active[0],dem) - h[t] # Subtract the hydraulic head from the ground surface

        sub_temp[sub_temp < 0] = 0 # Only include cells that have a hydraulic head less than the ground surface
        cells_below[t] = np.sum(np.sum(np.greater(sub_temp,0))) # Count all the cells that have a depth to groundwater greater than 0
        
        sub[t] = np.sum(np.sum(sub_temp)) # Sum the total depth to groundwater over all cells with head below the ground surface

    return sub,cells_below

def measureWaterQuality(heads,dem,active,bottom):
    wqual = np.zeros(360) # total difference between bottom of clay layer and head in layer 2
    h = list(wqual.copy()) # head map for each stress period
    cells_below = np.zeros(360) # Number of cells with head below bottom of clay layer
    
    # loop through all stress periods
    for t in range(348,360):
        # check during the last time step of each stress period
        h2 = np.multiply(active[0],heads.get_data(kstpkper=(8,t),mflay=1)) # heads binary file layer 2
        
        b0 = np.multiply(active[0],bottom[0])
        b1 = np.multiply(active[0],bottom[1])
        
        h[t] = np.multiply(np.less(h2,b1),b1) + np.multiply(np.multiply(np.less(h2,b0), np.greater(h2,b1)),h2) + np.multiply(np.greater(h2,b0),b0) # Matrix containing bottom of layer 2 if layer 2 is dry, plus the head in layer 2 if less than the bottom of layer 1, plus the bottom of layer 1 if the head is greater than layer 1
        
        cells_below[t] = np.sum(np.sum(np.less(h[t],b0)))

        wq_temp = b0 - h[t] # Subtract the head matrix from the bottom of the active clay layer
        
        wqual[t] = np.sum(np.sum(wq_temp)) # Sum the total depth to groundwater below the bottom of layer 1 in layer 2 over all cells

    return wqual,cells_below,h

def measureMound(heads,dem,active,LU,PhasePer):
    mound = 0 # cumulative head above DEM during model period
    urban_cells = 0
    mound_cells = 0
    
    for t in range(348,360):
        h1 = heads.get_data(kstpkper=(8,t),mflay=0) # heads binary file for last time step of stress period
        h2 = heads.get_data(kstpkper=(8,t),mflay=1) # heads binary file for last time step of stress period
      
        if t < PhasePer[0]:
            LUset = '1990'
        elif t < PhasePer[1]:
            LUset = '2000'
        else:
            LUset = '2010'
        
        lu_temp = np.multiply(LU[LUset]['ARRAY']['URBAN'],active[1])
        
        lu_active = np.sum(np.sum(lu_temp)) # Total urban area of active cells
        urban_cells += lu_active
                
        h = np.multiply(active[0],h1) + np.multiply((active[1]-active[0]),h2)
        above_dem = np.greater(h,np.multiply(dem,active[1]))
        lu_above = np.sum(np.sum(np.multiply(lu_temp,above_dem))) # Total area of cells with mounding above dem
        mound_cells += lu_above
    
    return mound,mound_cells,urban_cells

def get_objectives(heads, wellinfo, landuse, dem, active, bottom):
    
    cells_clay = np.sum(np.sum(active[0])) # Number of cells in Clay layer
#    energy = measureEnergy(heads, wellinfo, dem)
#    wq, wq_cells, hwq = measureWaterQuality(heads, dem, active, bottom)
#    mound, mound_cells, urban_cells = measureMound(heads, dem, active, landuse, [132,252])
    try:
        energy = measureEnergy(heads, wellinfo, dem, bottom)
    except:
        energy = np.nan
    try:
        sub, sub_cells = measureSubidence(heads, dem, active, bottom)
    except:
        sub, sub_cells = np.nan, np.nan
    try:
        wq, wq_cells, hwq = measureWaterQuality(heads, dem, active, bottom)
    except:
        wq, wq_cells, hwq = np.nan, np.nan, np.nan
    try:
        mound, mound_cells, urban_cells = measureMound(heads, dem, active, landuse, [132,252])
    except:
        mound, mound_cells, urban_cells = np.nan, np.nan, np.nan

    return energy, np.sum(wq_cells)/(cells_clay*12), mound_cells/urban_cells

def calculate_SOSWR(heads, stats):
    '''
    Calculates the sum-of-squared, weighted residual given weights (stats) and array (heads) with three columns: simulated head, observed head, observation ID
    '''
    soswr = 0
    maxerror = 0
    for i in range(len(stats)):
        currenterror = (heads[i+1][1] - heads[i+1][0])
        maxerror = max([maxerror, np.abs(currenterror)])
        soswr += ((1/stats[i]) * currenterror)**2
    
    return soswr, maxerror