# Script to generate points, lines, and polygons demarcating Kueppers-Powell forest inventory plot boundaries

#%%
# Load libraries
import pandas as pd
import geopandas as gpd
import math
import numpy as np
import os
import shapely
import re

#%%
# Set in directory
indir = '/Volumes/GoogleDrive/My Drive/Research/RMBL/RMBL-East River Watershed Forest Data/Data/Inventory Plots/Inventory_Plots_GPS_Data/2021'

sensordirs = [
    'KATZJ060810B',
    'KATZJ062410A',
    'KATZJ062510A',
    'KATZJ071909A',
    'KATZJ072311B',
    'KATZJ072812A',
    'KATZ080512A',
    'KATZ080517A',
    'WORSHAMM081308A',
    'WORSHAMM081812B',
    'WORSHAMM082008A',]

pd.set_option('display.max_rows', 200)

#%%
sensorfiles = []
for root, dirs, files in os.walk(indir):
    for d in dirs:
        if d in sensordirs:
            sensfile = os.path.join(root, d, 'Student__Project.shp')
            sensorfiles.append(sensfile)
#%%
allshots=[]
for f in sensorfiles:
    shots = gpd.read_file(f)
    allshots.append(shots)

allshots_gdf = gpd.GeoDataFrame(pd.concat(allshots, ignore_index=True))

#%%
allshots_gdf['Other'] = allshots_gdf['Other'].str.replace('GTH1', 'GT1')
allshots_gdf['Other'] = allshots_gdf['Other'].str.replace('ER-PVG1', 'SR-PVG1')

allshots_gdf['Site'] = allshots_gdf['Other'].str[:7].str.strip().str.replace(' ', '-')

allshots_gdf['Site'] = allshots_gdf['Site'].str.upper()

#%%
allshots_gdf = allshots_gdf[~allshots_gdf['Site'].str.contains('5')]
allshots_gdf = allshots_gdf[~allshots_gdf['Site'].str.contains('6')]

#%%
# Add inverse horizontal precision field for mean center weighting
allshots_gdf['Inv_Horz_Prec'] = 1/allshots_gdf['Horz_Prec']

#%%
# Calculate mean center of point clusters, grouping by corner direction
def cleansites(gpdf):
    
    meanpoints = []

    sites = list(gpdf['Site'].unique())
    
    def geomeancenter(pnt):
        lon = np.sum(pnt.geometry.x * pnt.Inv_Horz_Prec)/np.sum(pnt.Inv_Horz_Prec)
        lat = np.sum(pnt.geometry.y * pnt.Inv_Horz_Prec)/np.sum(pnt.Inv_Horz_Prec)
        geomcenter = shapely.geometry.Point([lon,lat])
        return geomcenter
    
    def meanhorzprec(pnt):
        mhp = np.mean(pnt['Horz_Prec'])
        return mhp

    def meanelv(pnt):
        mel = np.mean(pnt['GNSS_Heigh'])
        return mel
    
    def meanvertprec(pnt):
        mvp = np.mean(pnt['Vert_Prec'])
        return mvp
    
    for s in sites:
        ploti = gpdf[gpdf['Site'] == s]
        cgroups = geomeancenter(ploti)
        mhp = meanhorzprec(ploti)
        mel = meanelv(ploti)
        mvp = meanvertprec(ploti)

        gs = {'SITE_ID':s, 'HORZ_PREC':mhp, 'geometry':cgroups, 'ELEVATION':mel, 'VERT_PREC':mvp}
        meanpoints.append(gs)
    
    cgroups_gs = gpd.GeoDataFrame(meanpoints, crs='EPSG:32613').reset_index(drop=True)
    cgroups_gs['X_COORD'] = cgroups_gs.geometry.x
    cgroups_gs['Y_COORD'] = cgroups_gs.geometry.y
        
    return cgroups_gs[['SITE_ID', 'X_COORD', 'Y_COORD', 'HORZ_PREC', 'ELEVATION', 'VERT_PREC', 'geometry']]

#%%
cleaned_sensites = cleansites(allshots_gdf)

#%%
# Write out all sensor sites
outdir = '/Users/hmworsham/Desktop/'
filename = os.path.join(outdir, 'Kueppers_EastRiver_SensorSites_2021_WGS84UTM13N.shp')
cleaned_sensites.to_file(filename)

#%%
grouped = cleaned_sensites.groupby('SITE_ID')
gdfs=[]
for key, obj in grouped:
    gdfs.append(obj)

#%%
for key, obj in grouped:
    outdir='/Users/hmworsham/Desktop/'
    plotid = obj['SITE_ID'].values[0]
    filename = os.path.join(outdir, plotid+'_SensorSite_WGS84UTM13N')
    obj.to_file(filename)
