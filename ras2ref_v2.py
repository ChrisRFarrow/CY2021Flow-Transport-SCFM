import geopandas as gpd
import numpy as np
#import matplotlib.pyplot as plt
#import re
import os
import sys
import pandas as pd
import rasterio
from shapely.geometry import Point
import georaster as gr  # 1.26.1

'''
- xxx


'''


def create_new_dir(directory):
    # directory = os.path.dirname(file_path)
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
        print(f'Created a new directory {directory}\n')


def create_output_folders():
    create_new_dir('output')
    # create_new_dir('output/png')
    # create_new_dir('output/shp')
    # create_new_dir(f'output/png/conc_{var}')
    # create_new_dir(f'output/shp/conc_{var}')


def extract_raster(raster_path, xy):
    # raster_path : path to raster
    # xy: list or array of tuples of x,y i.e [(x1,y1),(x2,y2)...(xn,yn)]
    # returns list of length xy or sampled values
    raster = rasterio.open(raster_path)
    nodata = raster.nodata

    values = []
    for i, xyi in enumerate(xy):
        try:
            value = list(raster.sample([xyi]))[0][0]
        except:
            value = np.nan

        if value == nodata:  # if value is nodata from rasdter, set to nan
            value = np.nan

        values.append(value)

    # values = [item[0] for item in values] # list comprehension to get the value
    return values


if __name__ == "__main__":

    # [1] Load input files ----------------------------------------------------
    raster_file = sys.argv[1]
    #raster_file = "input/plume_cy2020/No3_20201.tif"

    # Specify an output file --------------------------------------------------
    ofile = sys.argv[2]
    #ofile = "output/No3_20201_check060922.dat"

    # reading input point csv file --------------------------------------------
    xy_file = f'input/grid_xy.csv'
    print(f'\nReading input cell coordinate csv file: {xy_file}\n')
    df_coor = pd.read_csv(xy_file)
    nrows, ncols = df_coor['I'].max(), df_coor['J'].max()
    print(f'nrow = {nrows}, ncols={ncols}\n')
    xy = [(df_coor.X, df_coor.Y)
          for i, df_coor in df_coor.iterrows()]

    # Extract concentration values from a raster file
    print(
        f'Reading input raster file: {raster_file} and extracting concentrations at grid cells\n')
    res = extract_raster(raster_file, xy)

    # Export to shapefile for checking ----------------------------------------
    # # polygon: Cols = ['I', 'J', 'CELLACTIVE', 'geometry']
    shp_file_poly = f'input/grid.shp'
    print(f'Reading input polygon grid file: {shp_file_poly}\n')
    gdf_poly = gpd.read_file(shp_file_poly)
    #gdf_poly = gdf_poly[['I', 'J', 'CELLACTIVE', 'geometry']]
    # gdf_poly.to_file('input/grid.shp')

    gdf_poly['Conc'] = res
    gdf_poly['Conc'].loc[gdf_poly['Conc'] < 0] = 0
    gdf_poly['Conc'].loc[gdf_poly['Conc'] > 1e12] = 0

    # Saving cocentration to polygon grid file for checking -------------------
    ofile1 = ofile.split('.')[0] + '.shp'
    print(
        f'Saving cocentration to polygon grid file for checking: {ofile1}\n')
    gdf_poly.to_file(ofile1)

    # save to csv file for checking -------------------------------------------
    ofile2 = ofile.split('.')[0] + '.csv'
    gdf_poly['X'] = df_coor['X']
    gdf_poly['Y'] = df_coor['Y']
    dfout2 = gdf_poly[['I', 'J', 'CELLACTIVE', 'X', 'Y', 'Conc']]
    print(f'Saving to output .csv file: {ofile2}\n')
    dfout2.to_csv(ofile2)

    # Reshape
    #arr = np.reshape(res, (663, 1164))
    arr = np.reshape(gdf_poly['Conc'].to_numpy(), (nrows, ncols))
    #arr[arr > 1e12] = 0
    #arr[arr < -1e12] = 0

    # arr[270:663, 640:1164] = 0  # remove plumes from other areas (for NO3 only? )
    # arr to ref file
    print(f'Saving to output .dat file: {ofile}\n')
    np.savetxt(ofile, arr, delimiter=" ", fmt='%.6f')


# Reference/Source
# https://github.com/rosskush/spatialpy/blob/main/spatialpy/utils/extraction.py
