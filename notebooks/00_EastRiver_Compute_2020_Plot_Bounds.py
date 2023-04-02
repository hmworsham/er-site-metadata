# Script to generate points, lines, and polygons demarcating Kueppers-Powell forest inventory plot boundaries


# Load libraries
import arcpy as ap 
import pandas as pd
import geopandas as gpd
import math
import os
import re

# Set arcpy workspace as output
ap.env.workspace = r"Y:\Desktop\RMBL\Projects\Watershed_Spatial_Dataset\Scratch\Plot_Bnd_2020"

allshots = ap.ListFeatureClasses(wild_card='*AllShots*')

for sf in allshots:
    print(allshots.index(sf), sf)
allshots
i = ap.ListFeatureClasses(wild_card='*PL1_Corners*')[0]
i

# Add temporary field to clean plot corner descriptions from a bunch of different fields 
ap.AddField_management(i, 'Comment', 'TEXT', '', '', 50, '', 'NULLABLE', 'NON_REQUIRED')
ap.AddField_management(i, 'Other', 'TEXT', '', '', 50, '', 'NULLABLE', 'NON_REQUIRED')
ap.AddField_management(i, 'Other2', 'TEXT', '', '', 50, '', 'NULLABLE', 'NON_REQUIRED')
ap.AddField_management(i, 'CornersTmp', 'TEXT', '', '', 50, '', 'NULLABLE', 'NON_REQUIRED')

# Function to populate the new CornersTemp field with values from possible fields
def cnrpop(shp):
    fields = ('Other', 'Other2', 'Comment', 'CornersTmp')
    cursor = ap.da.UpdateCursor(shp, fields)
    for row in cursor:
        if row[0] != ' ':
            print row[0]
            row[3] = row[0]
        elif row[1] != ' ':
            print row[1]
            row[3] = row[1]
        elif row[2] != '':
            print row[2]
            row[3] = row[2]
        else:
            row[3] = 'NA'
        cursor.updateRow(row)

cnrpop(i)

# Add canonical 'Corners' field populate with cardinal value for easy sorting
ap.AddField_management(i, 'Corners', 'TEXT', '', '', 10, '', 'NULLABLE', 'NON_REQUIRED')

grepexpr = """
def corner(field):
    if "north" in field:
        return "1N"
    if "east" in field:
        return "2E"
    if "south" in field:
        return "3S"
    if "west" in field:
        return "4W"
    if "n corner" in field:
        return "1N"
    if "\se corner" in field:
        return "2E"
    if "s corner" in field:
        return "3S"
    if "\sw corner" in field:
        return "4W"
    if "nw" in field:
        return "1NW"
    if "\sne\s" in field:
        return "2NE"
    if "ne\d" in field:
        return "2NE"
    if "se" in field:
        return "3SE"
    if "sw" in field:
        return "4SW"
    if "NW" in field:
        return "1NW"
    if "NE" in field:
        return "2NE"
    if "SE" in field:
        return "3SE"
    if "SW" in field:
        return "4SW"
    else:
        return "2NE" """
        
ap.CalculateField_management(i, "Corners", "corner(!CornersTmp!)", 'PYTHON', grepexpr)

# Add inverse horizontal precision field for mean center weighting
ap.AddField_management(i, "Inv_Horz_P", "FLOAT", "", "", "", "", "NULLABLE", "REQUIRED")
ap.CalculateField_management(i, "Inv_Horz_P", '1 / !Horz_Prec!', 'PYTHON')

# Calculate mean center of point clusters, grouping by corner direction
ap.MeanCenter_stats(i, ap.env.workspace + "\\" + i[0:-12] + 'MeanCtr', 'Inv_Horz_P', 'Corners', 'Horz_Prec')

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