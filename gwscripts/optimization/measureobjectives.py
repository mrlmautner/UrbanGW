# -*- coding: utf-8 -*-
"""
Created on Thu May 24 18:59:47 2018

@author: MM
"""
import numpy as np
import pandas as pd
import math
import calendar

def measureEnergy(heads,SUPPLY_DICT,DEM,bottom):
    E = 0
    
    efficiency = 0.7
    MJc = 9.81 # MegaJoules to lift 1 MegaLiter 1 meter
    kWhc = 3.6 # MegaJoules per kWh
    
    pd.DatetimeIndex(freq='M',start='01/31/1986',end='12/31/2013')
    coords = np.zeros(2)
    
    for i, p in SUPPLY_DICT.items():
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
                    b = bottom[int(coords[0]),int(coords[1])]
                    if wellHead < b:
                        wellHead = b
                    surface = DEM[int(coords[0]),int(coords[1])]
                    depthtowater = max(surface-wellHead,0)
                    pumping = n[3] * (-1) * 0.001 # m3/d --> ML/d
                    E += (efficiency * depthtowater * pumping * MJc / kWhc)*d # kWh per month (stress period) of pumping
    
    return E

def measureSubidence(heads,DEM,ACTIVE_LYR1,THICK1,bottom):
    sub = 0 # total head
    cells = np.sum(np.sum(ACTIVE_LYR1)) # Number of cells in Clay layer
    
    # check during the reference month of the last year (March)
    h = heads.get_data(kstpkper=(8,359),mflay=1) # heads binary file for dth day of stress period
    # loop through all cells
    for i in range(int(ACTIVE_LYR1.shape[0])):    
        for j in range(int(ACTIVE_LYR1.shape[1])):
            # Only evaluate cells under the clay layer (where Layer 1 is active)
            if ACTIVE_LYR1[i,j]:
                # If layer 1 is dry (h < 0), then use bottom of layer 2
                if h[i,j]<0:
                    sub += DEM[i,j] - bottom[i,j]
                else:
                    # Only add subsidence measure if the hydraulic head is less than the ground surface
                    if h[i,j] < DEM[i,j]:
                        sub += DEM[i,j] - h[i,j]
    
    return sub,cells

def measureMound(heads,DEM,ACTIVE_LYR1,LU,PhasePer):
    mound = 0 # cumulative head above DEM during model period
    maxmound = 0
    cells = 0
    urban_cells = 0
    
    for t in range(23,360):
        d = calendar.monthrange(1984+math.floor(t/12),(t%12)+1)[1] # number of days in stress period
        h = heads.get_data(kstpkper=(8,t),mflay=0) # heads binary file for dth day of stress period
        
        if t < PhasePer[0]:
            LUset = '1990'
        elif t < PhasePer[1]:
            LUset = '2000'
        else:
            LUset = '2010'
        
        # loop through all cells
        for i in range(int(DEM.shape[0])):    
            for j in range(int(DEM.shape[1])):
                # Only evaluate cells where the head is above ground level and urban cell
                if LU[LUset]['ARRAY']['URBAN'][i,j] > 0:
                    urban_cells += LU[LUset]['ARRAY']['URBAN'][i,j]
                    if DEM[i,j] < h[i,j]:
                        cells += LU[LUset]['ARRAY']['URBAN'][i,j]
                        mound += h[i,j] - DEM[i,j]
                        if (h[i,j] - DEM[i,j]) > maxmound:
                            maxmound = h[i,j] - DEM[i,j]
    
    return mound,cells,urban_cells,maxmound

def get_objectives(heads, wellinfo, landuse, dem, active, thickness):

    try:
        energy = measureEnergy(heads, wellinfo, dem, bottom)
    except:
        energy = np.nan
        print('Energy error')
        
    try:
        sub, sub_cells = measureSubidence(heads, dem, active, thickness, bottom)
        sub_obj = sub/sub_cells
    except:
        sub_obj = np.nan
        print('Subsidence error')
        
    try:
        mound, mound_cells, urban_cells, maxmound = measureMound(heads, dem, active, landuse, [132,252])
        mound_obj = mound_cells/urban_cells
    except:
        mound_obj = np.nan
        print('Mound error')
        
    return energy, sub_obj, mound_obj

def calculate_SOSWR(heads, stats):
    '''
    Calculates the sum-of-squared, weighted residual given weights (stats) and array (heads) with three columns: simulated head, observed head, observation ID
    '''
    soswr = 0
    maxerror = 0
    for i in range(len(stats)):
        currenterror = heads[i+1][1] - heads[i+1][0]
        maxerror = max([maxerror, np.abs(currenterror)])
        soswr += ((1/stats[i]) * currenterror)**2
    
    return soswr, maxerror