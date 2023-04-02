# Script to generate points, lines, and polygons demarcating Kueppers-Powell forest inventory plot boundaries

#%%
# Load libraries
import pandas as pd
import geopandas as gpd
import math
import numpy as np
import os
import re

#%%
# Set arcpy workspace as output
indir = '/Volumes/GoogleDrive/My Drive/Research/RMBL/RMBL-East River Watershed Forest Data/Data/Inventory Plots/Inventory_Plots_GPS_Data/2021'

cornerdirs = [
    'KATZJ061410A',
    'KATZJ062808A',
    'KATZJ063010B',
    'KATZJ070909A',
    'KATZJ071908A',
    'KATZJ072110A',
    'KATZJ072608B',
    'KATZJ072609A',
    'KATZJ072610A',
    'KATZJ072610C',
    'KATZJ072908A',
    'KATZJ072910A',
    'KATZ080208A',
    'KATZ080410A',
    'KATZ080709A_X']

pd.set_option('display.max_rows', 200)

#%%
cornerfiles = []
for root, dirs, files in os.walk(indir):
    for d in dirs:
        if d in cornerdirs:
            cornerfile = os.path.join(root, d, 'Student__Project.shp')
            cornerfiles.append(cornerfile)

allshots=[]
for f in cornerfiles:
    shots = gpd.read_file(f)
    allshots.append(shots)

allshots_gdf = gpd.GeoDataFrame(pd.concat(allshots, ignore_index=True))

allshots_gdf['Other'] = allshots_gdf['Other'].str.replace('SL', 'SR')
allshots_gdf['Other'] = allshots_gdf['Other'].str.replace('ER-PVG1', 'SR-PVG1')

allshots_gdf['Plot'] = allshots_gdf['Other'].str[:7].str.replace(' ', '-')

allshots_gdf['CornersTmp'] = None
allshots_gdf['CornersTmp'] = allshots_gdf['Other'].str[7:]
allshots_gdf['CornersTmp'] = allshots_gdf['CornersTmp'].str.upper()

allshots_gdf = allshots_gdf[~allshots_gdf['CornersTmp'].str.contains('50')]

#%%
# Add canonical 'Corners' field populate with cardinal value for easy sorting
allshots_gdf['Corners'] = None

def corner(field):
    if " NORTH " in field:
        return "1N"
    if " EAST " in field:
        return "2E"
    if " SOUTH " in field:
        return "3S"
    if " WEST " in field:
        return "4W"
    if " NORTHWEST " in field:
        return "1NW"
    if " NORTHEAST " in field:
        return "2NE"
    if " SOUTHEAST " in field:
        return "3SE"
    if " SOUTHWEST " in field:
        return "4SW"
    if " SE " in field:
        return "3SE"
    if " SW " in field:
        return "4SW"
    if " NW " in field:
        return "1NW"
    if " NE " in field:
        return "2NE"
    else:
        return "Error"

#%%
cancorners = list(map(corner, allshots_gdf['CornersTmp']))
cancorners
allshots_gdf['Corners'] = cancorners

# Add inverse horizontal precision field for mean center weighting
allshots_gdf['Inv_Horz_Prec'] = 1/allshots_gdf['Horz_Prec']

#%%
# Calculate mean center of point clusters, grouping by corner direction
def cleancorners(gpdf):
    
    meanpoints = []

    plots = list(gpdf['Plot'].unique())
    
    def geomeancenter(corner):
        lon = np.sum(corner.geometry.x * corner.Inv_Horz_Prec)/np.sum(corner.Inv_Horz_Prec)
        lat = np.sum(corner.geometry.y * corner.Inv_Horz_Prec)/np.sum(corner.Inv_Horz_Prec)
        geomcenter = shapely.geometry.Point([lon,lat])
        return geomcenter
    
    def meanhorzprec(corner):
        mhp = np.mean(corner['Horz_Prec'])
        return mhp

    for p in plots:
        ploti = gpdf[gpdf['Plot'] == p]
        cgroups = ploti.groupby(ploti['Corners']).apply(geomeancenter)
        mhp = ploti.groupby(ploti['Corners']).apply(meanhorzprec)
        gs = {'PLOT_ID':p, 'HORZ_PREC':mhp, 'geometry':cgroups}
        cgroups_gs = gpd.GeoDataFrame(gs, crs='EPSG:32613').reset_index()
        cgroups_gs['X_COORD'] = cgroups_gs.geometry.x
        cgroups_gs['Y_COORD'] = cgroups_gs.geometry.y
        cgroups_gs.rename(columns={'Corners':'CORNER'}, inplace=True)
        meanpoints.append(cgroups_gs[['PLOT_ID', 'CORNER', 'X_COORD', 'Y_COORD', 'HORZ_PREC', 'geometry']])
        
    return meanpoints

#%%
cleaned_corners = cleancorners(allshots_gdf)

cleaned_corners_all = gpd.GeoDataFrame(pd.concat(cleaned_corners)).reset_index(drop=True)
cleaned_corners[-1]

#%%
# Write out plot corners per plot

for i in cleaned_corners:
    outdir='/Users/hmworsham/Desktop/Corner_Points/'
    plotid = i['PLOT_ID'][0]
    filename = os.path.join(outdir, plotid+'_Bound_pts_WGS84UTM13N')
    i.to_file(filename)

#%%
# Create polygons
polygons=[]
for i in cleaned_corners:
    plotid = i['PLOT_ID'][0]
    ch = i.unary_union.convex_hull
    gs = {'PLOT_ID':plotid, 'geometry':ch}
    polygons.append(gs)

#%%
all_poly = gpd.GeoDataFrame(polygons, crs=cleaned_corners_all.crs).reset_index(drop=True)

#%%
all_poly['AREA_GEO'] = all_poly.area
all_poly['PERIM_GEO'] = all_poly.length
all_poly['CENTROID_X'] = all_poly.centroid.x
all_poly['CENTROID_Y'] = all_poly.centroid.y
all_poly['EXT_MIN_X'] = all_poly.bounds['minx']
all_poly['EXT_MIN_Y'] = all_poly.bounds['miny']
all_poly['EXT_MAX_X'] = all_poly.bounds['maxx']
all_poly['EXT_MAX_Y'] = all_poly.bounds['maxy']
all_poly['EXT_MAX_Y'] = all_poly.bounds['maxy']
all_poly['EXT_MAX_Y'] = all_poly.bounds['maxy']
all_poly['GEOMCTR_X'] = all_poly.to_crs('EPSG:4326').centroid.x
all_poly['GEOMCTR_Y'] = all_poly.to_crs('EPSG:4326').centroid.y
all_poly['PLOT_NAME'] = [
    'Point Lookout north slope 2', 
    'Carbon Creek valley 1',
    'Snodgrass northeast slope 3',
    'Carbon Creek valley 3',
    'Poverty Gulch 1',
    'Snodgrass northeast slope 1', 
    'Coal Valley north 2', 
    'Mount Emmons 1', 
    'Coal Valley south 1',
    'Coal Valley north 1'
]
all_poly = all_poly[list(allplots_format.columns)]

all_poly['INSTALL_YR'] = [
    2019,
    2019,
    2021,
    2021,
    2021,
    2021,
    2021,
    2021,
    2021,
    2021
]

#%%
all_plots = allplots_format[~allplots_format['PLOT_ID'].isin(all_poly['PLOT_ID'])]

all_plots['INSTALL_YR'] = [
    2019,
    2019,
    2019,
    2019,
    2020,
    2018,
    2020,
    2019,
    2020,
    2018,
    2020,
    2018
]

all_plots21 = all_plots.append(all_poly).copy().reset_index(drop=True)

#%%
outdir = '/Users/hmworsham/Desktop/AllPlots/'
filename = os.path.join(outdir, 'Kueppers_EastRiver_AllPlots_2021_WGS84UTM13N.shp')
all_plots21.to_file(filename)

#%%
grouped = all_poly.groupby('PLOT_ID')
gdfs=[]
for key, obj in grouped:
    gdfs.append(obj)

#%%
for key, obj in grouped:
    outdir='/Users/hmworsham/Desktop/Polygons/'
    plotid = obj['PLOT_ID'].values[0]
    filename = os.path.join(outdir, plotid+'_Bound_poly_WGS84UTM13N')
    obj.to_file(filename)
#%%
# Split into linestrings
all_lines = all_poly.copy()
all_lines.geometry = all_poly.boundary

linegroups = all_lines.groupby('PLOT_ID')

for key, obj in linegroups:
    outdir='/Users/hmworsham/Desktop/Lines/'
    plotid = obj['PLOT_ID'].values[0]
    filename = os.path.join(outdir, plotid+'_Bound_lines_WGS84UTM13N')
    obj.to_file(filename)

#%%
all_centers = all_poly.copy()
all_centers.geometry = all_poly.centroid
all_centers = all_centers[['PLOT_ID', 'CENTROID_X', 'CENTROID_Y', 'geometry']]
all_centers.rename(columns={'CENTROID_X':'X_COORD', 'CENTROID_Y':'Y_COORD'},inplace=True)
centgroups = all_centers.groupby('PLOT_ID')

#%%
for key, obj in centgroups:
    outdir='/Users/hmworsham/Desktop/Center_Points/'
    plotid = obj['PLOT_ID'].values[0]
    filename = os.path.join(outdir, plotid+'_PlotCenter_WGS84UTM13N')
    obj.to_file(filename)

#%%
# NEED TO MATCH FORMATS OF EXISTING FILES
allplots_format = gpd.read_file('/Volumes/GoogleDrive/My Drive/Research/RMBL/RMBL-East River Watershed Forest Data/Data/Geospatial/Kueppers_EastRiver_Plot_Shapefiles_2020_WGS84UTM13N/AllPlots/Kueppers_EastRiver_AllPlots_2020_WGS84UTM13N.shp')

cornerpts_format = gpd.read_file('/Volumes/GoogleDrive/My Drive/Research/RMBL/RMBL-East River Watershed Forest Data/Data/Geospatial/Kueppers_EastRiver_Plot_Shapefiles_2020_WGS84UTM13N/Corner_Points/CC-UC1_Bound_pts_WGS84UTM13N.shp')

lines_format = gpd.read_file('/Volumes/GoogleDrive/My Drive/Research/RMBL/RMBL-East River Watershed Forest Data/Data/Geospatial/Kueppers_EastRiver_Plot_Shapefiles_2020_WGS84UTM13N/Lines/CC-UC1_Bound_lines_WGS84UTM13N.shp')

centers_format = gpd.read_file('/Volumes/GoogleDrive/My Drive/Research/RMBL/RMBL-East River Watershed Forest Data/Data/Geospatial/Kueppers_EastRiver_Plot_Shapefiles_2020_WGS84UTM13N/Center_Points/CC-UC1_PlotCenter_WGS84UTM13N.shp')

##################################
#%%


# Generate lines from mean center points
pt2lnin = ap.env.workspace + "\\" + i[0:-12] + 'MeanCtr.shp'
pt2lnout = pt2lnin[0:-4] + '_line.shp'
ap.PointsToLine_management(pt2lnin, pt2lnout, '', 'Corners', 'CLOSE')

# Generate polygon from lines
ln2plin = pt2lnout
ln2plout = pt2lnout[0:-9] + '_poly.shp'
ap.FeatureToPolygon_management(ln2plin, ln2plout, '', 'ATTRIBUTES')

# Calculate geometries
ap.AddGeometryAttributes_management(ap.env.workspace + "\\PL1_Corners_MeanCtr_poly.shp", 'AREA; AREA_GEODESIC; CENTROID; EXTENT; PERIMETER_LENGTH; PERIMETER_LENGTH_GEODESIC', 'METERS', 'SQUARE_METERS')

# Gather the points, lines, and polys into separate groups
allpoints = ap.ListFeatureClasses(wild_card = '*MeanCtr.shp') + ap.ListFeatureClasses(wild_card = '*pts*')
alllines = ap.ListFeatureClasses(wild_card = '*line*')
allpolys = ap.ListFeatureClasses(wild_card = '*poly*')

# Reproject shapefiles into WGS 1984 UTM Zone 13N spatial reference (used in RMBL SDP)
outdir = 'Y:\\Desktop\\RMBL\\Projects\\Watershed_Spatial_Dataset\\Output\\Kueppers_Plot_Bnd_2020_WGS84UTM13N\\'

for sf in allpoints:
    ap.Project_management(sf, outdir + sf.split('_')[0] + '_Bound_pts_WGS84UTM13N.shp', ap.SpatialReference('WGS 1984 UTM Zone 13N'), 'NAD_1983_to_WGS_1984_5')

for sf in alllines:
    ap.Project_management(sf, outdir + sf.split('_')[0] + '_Bound_lines_WGS84UTM13N.shp', ap.SpatialReference('WGS 1984 UTM Zone 13N'), 'NAD_1983_to_WGS_1984_5')

for sf in allpolys:
    ap.Project_management(sf, outdir + sf.split('_')[0] + '_Bound_poly_WGS84UTM13N.shp', ap.SpatialReference('WGS 1984 UTM Zone 13N'), 'NAD_1983_to_WGS_1984_5')

# add fields to hold geometric center lat/long in decimal degrees
ap.env.workspace = outdir
wgsutm_allpolys = ap.ListFeatureClasses(wild_card = '*poly*')
wgsutm_allpolys

for sf in wgsutm_allpolys:
    ap.DeleteField_management(sf, 'GEOMCTR_X')
    ap.DeleteField_management(sf, 'GEOMCTR_Y')
    ap.AddField_management(sf, 'GEOMCTR_X', 'LONG', 12, 7, field_is_nullable = 'NULLABLE', field_is_required = 'NON_REQUIRED')
    ap.AddField_management(sf, 'GEOMCTR_Y', 'LONG', 12, 7, field_is_nullable = 'NULLABLE', field_is_required = 'NON_REQUIRED')

# calculate geometric center of plot polygons
#wgsutm_allpolys[0].split('_')[0]
for sf in wgsutm_allpolys:
    ap.MeanCenter_stats(sf, outdir + "\\" + sf.split('_')[0] + '_PlotCenter_WGS84UTM13N.shp')

# Merge polygons into one 2020 plots shapefile
ap.env.workspace = outdir
wgsutm_allpolys = ap.ListFeatureClasses(wild_card = '*poly*')
wgsutm_allpolys

allplotsin = wgsutm_allpolys
allplotsout = outdir + 'Kueppers_AllPlots_2020_WGS84UTM13N.shp'
ap.Merge_management(allplotsin, allplotsout)
#ap.Project_management(allplotsout, allplotsout[:-15] + 'WGS84UTM13N.shp', ap.SpatialReference('WGS 1984 UTM Zone 13N'))
# %%
