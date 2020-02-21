# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 16:38:07 2018

@author: MM
"""

import numpy as np
from osgeo import gdal
import pandas as pd

# Generate .asc header
def getHeader(ncols,nrows,xll,yll,cellsize,noData):
    header = "ncols        %s\n" % ncols
    header += "nrows        %s\n" % nrows
    header += "xllcorner    %s\n" % xll
    header += "yllcorner    %s\n" % yll
    header += "cellsize     %s\n" % cellsize
    header += "NODATA_value     %s" % noData
    return header

def openASC(filename,band=1):
    ds = gdal.Open(filename)
    band = ds.GetRasterBand(band)
    dsAsArray = band.ReadAsArray()
    return dsAsArray

def sampleGD(filename,newfile,header,xll,yll,xur,yur,cellsize,band=1):
    ds = gdal.Warp(newfile,filename,outputBounds=[xll, yll, xur, yur],xRes=cellsize,yRes=cellsize,resampleAlg=gdal.GRA_Mode)
    band = ds.GetRasterBand(band)
    dsAsArray = band.ReadAsArray()
    
    np.savetxt(newfile, dsAsArray, header=header, fmt="%d",comments='')
    return dsAsArray

def averageGD(filename,newfile,header,xll,yll,xur,yur,cellsize,band=1):
    ds = gdal.Warp(newfile,filename,outputBounds=[xll, yll, xur, yur], xRes=cellsize, yRes=cellsize, resampleAlg=gdal.GRA_Average)
    band = ds.GetRasterBand(band)
    dsAsArray = band.ReadAsArray()
    dsAsArray = np.clip(dsAsArray, 0, np.inf)
    
    np.savetxt(newfile, dsAsArray, header=header, fmt="%1.2f",comments='')
    return dsAsArray

def generateWellData(filename, pandasWD, year):
    mths = ["ENE","FEB","MAR","ABR","MAY","JUN","JUL","AGO","SEP","OCT","NOV","DIC"]
    days_m = [31,28,31,30,31,30,31,31,30,31,30,31]
    wellFeatures = ['X','Y','CVE_ENT','CVE_MUN','NOM_MUN','CODE']
    
    newWDset = pd.read_csv(filename,skiprows=6,usecols=range(2,21),index_col=2,na_values='#VALUE!')
    newWDset[mths] = newWDset[mths].apply(pd.to_numeric, errors='coerce')
    newWDset = newWDset.drop(labels='nan')
    
    temp = wellFeatures.copy()
    temp.append('Flow_m3')
    stackedWDset = pd.DataFrame(columns=temp)
    
    for i,m in enumerate(mths):
        temp = wellFeatures.copy()
        temp.append(m)
        tempPD = newWDset[temp].rename(columns={m:'Flow_m3'})
        tempPD.Flow_m3 *= -1/days_m[i]
        tempPD['S_PER'] = pd.Series(np.ones(newWDset.shape[0])*(year-1984)*12+1+i,index=newWDset.index)
        tempPD['E_PER'] = pd.Series(np.ones(newWDset.shape[0])*(year-1984)*12+1+i+1,index=newWDset.index)
        tempPD['S_DAY'] = pd.Series(np.ones(newWDset.shape[0])*(year-1984)*365+sum(days_m[:i]),index=newWDset.index)
        tempPD['E_DAY'] = pd.Series(np.ones(newWDset.shape[0])*(year-1984)*365+sum(days_m[:(i+1)]),index=newWDset.index)
        tempPD['FLOW_M3-S'] = pd.Series(np.ones(newWDset.shape[0]),index=newWDset.index)
        tempPD['START_YR'] = pd.Series(np.ones(newWDset.shape[0])*year,index=newWDset.index)
        tempPD['END_YR'] = pd.Series(np.ones(newWDset.shape[0])*(year+1),index=newWDset.index)
        tempPD['DATASET'] = pd.Series(['SACM']*newWDset.shape[0],index=newWDset.index)
        stackedWDset = stackedWDset.append(tempPD)
    
    stackedWDset = stackedWDset.reindex(index=stackedWDset.index.astype(int))
    stackedWDset = stackedWDset[stackedWDset.Flow_m3 < 0]
        
    pandasWD = pandasWD.append(stackedWDset)
    
    return pandasWD