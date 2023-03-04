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


def extract_raster2gdf(raster_path, gdf):
    # raster_path : path to raster
    # gdf: Geo Pandas GeoDataFrame with coordinates stored in the 'geometry' attribute
    # returns list of length gdf or sampled values
    xy = [(dfrow['geometry'].x, dfrow['geometry'].y)
          for i, dfrow in gdf.iterrows()]
    values = extract_raster(raster_path, xy)

    return values


def raster2pts(rasobj, column='Value'):
    crs = rasobj.crs
    array = rasobj.read(1)

    nrow, ncol = rasobj.height, rasobj.width

    data = []
    i = 0
    for r in range(nrow):
        for c in range(ncol):

            x, y = rasobj.xy(r, c)
            z = array[r, c]
            data.append([z, Point(x, y)])

            i += 1

    gdf = gpd.GeoDataFrame(data, columns=[column, 'geometry'], crs=crs)

    return gdf


def read_ref(ifile, nr, nc):  # Return an array
    # Generate an empty arr
    arr = np.empty([nr, nc])
    arr[:, :] = np.nan

    count = 0
    val = []
    i = 0
    with open(ifile) as fp:
        for line in fp:
            line_spit = line.split()
            line_spit_conv = [float(line_spit[k])
                              for k in range(len(line_spit))]
            # print(f'i = {i}, count = {count}, {line.split()[0]}, {len(val)}\n')
            val = val + line_spit_conv
            count += 1
            if len(val) == nc:
                arr[i, :] = val
                val = []
                i += 1
    return arr


def writearray(array, ncols, fname, dtype):
    """
    Write a 2D array to an output file in 10E12.4 format or 10I10 format
    Input parameters
    ----------------
    array is a 2d array with shape(nrow,ncol)
    Integer or double
    fname is the name of the output file
    dtype is a string
    Can be 'int' or 'double'
    """
    assert(len(array.shape) == 2)
    assert(dtype == 'int' or dtype == 'double')
    nrow, ncol = array.shape

    # print 'writing array to ' + fname

    """
	lookupdict={}
	with open('C:/Projects/Google Drive/PhD/Inversion_of_CategoricalFields/Task1/hylookuptable.dat','r') as fin:
		for line in fin.readlines():
			key,val = line.split()
			# print key,val
			lookupdict.update({int(key):np.double(val)})
	fin.close()
	# print lookupdict
	newarray = np.zeros((ny,nx),dtype='double')
	for key, val in lookupdict.iteritems(): newarray[array==key] = np.double(val)
	"""
    with open(fname, 'w') as fout:
        for i in range(0, nrow):
            jprev = 0
            jnew = 0
            while(jnew < ncol):
                if(jprev+ncols) > ncol:
                    jnew = ncol
                else:
                    jnew = jprev+ncols
                line = ''
                # print jnew,jprev
                for k in range(jprev, jnew):
                    if(dtype == 'int'):
                        line = line+'{:3d}'.format(array[i][k])
                    elif(dtype == 'double'):
                        # line = line + '{:ncols.4e}'.format(array[i][k])
                        # line = line + f'{array[i][k]:12.4e}'
                        # line = line + f'{array[i][k]:15.7e}'  # bot*.ref
                        line = line + f'{array[i][k]:14.6e}'  # thk*.ref
                        # f'{number:9.4f}'

                jprev = jnew
                fout.write(f'{line}\n')

    fout.close()
    print(f'Saved {fname}\n')


if __name__ == "__main__":

    # [1] Load input file -----------------------------------------------------
    #raster_file = f'input/plume_cy2020/C14_20201.tif'
    #ofile = f'output/C14_20201.dat'

    shp_file_point = f'input/grid_point.shp'  # points
    shp_file_poly = f'input/grid2.shp'        # polygon

    #raster_file = sys.argv[1]
    raster_file = "input/plume_cy2020/No3_20201.tif"

    #ofile = sys.argv[2]
    ofile = "output/No3_20201_check060922.dat"

    # read input file ---------------------------------------------------------
    # dfin = pd.read_csv(raster_file)
    xy_file = f'input/grid_xy.csv'
    df_coor = pd.read_csv(xy_file)

    #gdf = gdf[gdf['K'] == 1]
    gdf = gpd.read_file(shp_file_point)
    gdf = gdf[['I', 'J', 'CELLACTIVE', 'geometry']]

    # returns list of length gdf or sampled values
    xy = [(gdf['geometry'].x, gdf['geometry'].y)
          for i, gdf in gdf.iterrows()]
    xy = [(df_coor.X, df_coor.Y)
          for i, df_coor in df_coor.iterrows()]
    # Extract values
    res = extract_raster(raster_file, xy)

    # Export to shapefile for checking
    gdf_poly = gpd.read_file(shp_file_poly)
    gdf = gdf[['I', 'J', 'CELLACTIVE', 'geometry']]
    gdf_poly['Conc'] = res
    gdf_poly['Conc'].loc[gdf_poly['Conc'] < -1e12] = 0
    gdf_poly['Conc'].loc[gdf_poly['Conc'] > 1e12] = 0
    gdf_poly.to_file(f'output/shp/NO3.shp')

    # my_image = gr.SingleBandRaster(raster_file, load_data=extent, latlon=True)
    #my_image = gr.SingleBandRaster(raster_file)
    # print(my_image.extent)
    #print(my_image.nx, my_image.ny)
    #xres, yres = my_image.xres, my_image.yres
    #v = my_image.r
    #v = v.astype('float')

    # Use rasterio
    #my_image2 = rasterio.open(raster_file)

    # Testing a point

    #xt, yt, = 576942.5, 153787.5
    # xt, yt, = 574842.5, 154567.5  # cell 5, 745, no raster pixel
    # xt, yt, = 574692.5, 154755  # cell 1,735, no raster pixel

    # xt, yt, = 574827.5, 154507.5  # cell 9, 744, val=15,677.7

    # cell 110, 954, val=0 (QGIS), val=94.8 (hp's script)
    #xt, yt, = 577972.5, 152992.5

    #val = my_image.value_at_coords(xt, yt, window=1)

    #print(f'v={val} at the test point x={xt}, y={yt}')

    # val = []
    # for i in range(df_coor.shape[0]):
    #     val_tmp = my_image.value_at_coords(
    #         df_coor['X'][i], df_coor['Y'][i], window=1)
    #     # print(val_tmp.size)
    #     if val_tmp is None:
    #         val_tmp = 0
    #     #
    #     val.append(val_tmp)

    # # Assign concentration value to column val
    # df_coor['val'] = val
    # df_coor['val'][df_coor['val'] > 1e12] = 0
    # df_coor['val'][df_coor['val'] < -1e12] = 0

    # # save
    ofile2 = ofile.split('.')[0] + '.csv'
    dfout2 = gdf[['I', 'J', 'CELLACTIVE', 'Conc']]
    dfout2.to_csv(ofile2)

    # Reshape
    #arr = np.reshape(df_coor.val.to_numpy(), (663, 1164))
    arr = np.reshape(res, (663, 1164))
    arr[arr > 1e12] = 0
    arr[arr < -1e12] = 0

    # arr[270:663, 640:1164] = 0  # remove plumes from other areas (for NO3 only? )
    # arr to ref file
    #fname = f'output/Hexavalent_Chromium_CY2020.dat'
    np.savetxt(ofile, arr, delimiter=" ", fmt='%.6f')

# Reference/Source
# https://github.com/rosskush/spatialpy/blob/main/spatialpy/utils/extraction.py
