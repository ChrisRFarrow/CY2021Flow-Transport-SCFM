import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import re
import os
import sys
import pandas as pd
import rasterio
from shapely.geometry import Point
import georaster as gr

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


def get_gdf_extent(gdf):
    bounds = gdf.bounds
    ulc = (bounds['minx'].min(), bounds['maxy'].max())
    lrc = (bounds['maxx'].max(), bounds['miny'].min())
    return ulc, lrc


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

    # shp_file = f'input/grid_point.shp'

    raster_file = sys.argv[1]
    ofile = sys.argv[2]

    # read input file ---------------------------------------------------------
    # dfin = pd.read_csv(raster_file)
    xy_file = f'input/grid_xy.csv'
    df_coor = pd.read_csv(xy_file)

    #gdf = gpd.read_file(shp_file)
    # returns list of length gdf or sampled values
    # xy = [(gdf['geometry'].x, gdf['geometry'].y)
    #      for i, gdf in gdf.iterrows()]
    xy = [(df_coor.X, df_coor.Y)
          for i, df_coor in df_coor.iterrows()]
    # Extract values
    res = extract_raster(raster_file, xy)

    # my_image = gr.SingleBandRaster(raster_file, load_data=extent, latlon=True)
    my_image = gr.SingleBandRaster(raster_file)
    print(my_image.extent)
    print(my_image.nx, my_image.ny)
    xres, yres = my_image.xres, my_image.yres
    #v = my_image.r
    #v = v.astype('float')

    # Testing a point

    #xt, yt, = 576942.5, 153787.5
    xt, yt, = 574692.5, 154755  # cell 1,735, no raster pixel

    xt, yt, = 574842.5, 154567.5  # cell 5, 745, no raster pixel
    val = my_image.value_at_coords(xt, yt, window=1)
    print(f'v={val} at the test point x={xt}, y={yt}')
    val = []
    for i in range(df_coor.shape[0]):
        val_tmp = my_image.value_at_coords(
            df_coor['X'][i], df_coor['Y'][i], window=1)
        # print(val_tmp.size)
        if val_tmp is None:
            val_tmp = 0
        #
        val.append(val_tmp)

    # Assign concentration value to column val
    df_coor['val'] = val
    df_coor['val'][df_coor['val'] > 1e12] = 0
    df_coor['val'][df_coor['val'] < -1e12] = 0
    # save

    # df_coor.to_csv(f'output/Hexavalent_Chromium.csv')

    # Reshape
    arr = np.reshape(df_coor.val.to_numpy(), (663, 1164))
    arr[270:663, 640:1164] = 0  # remove plumes from other areas
    # arr to ref file
    #fname = f'output/Hexavalent_Chromium_CY2020.dat'
    np.savetxt(ofile, arr, delimiter=" ", fmt='%.6f')

# Reference/Source
# https://github.com/rosskush/spatialpy/blob/main/spatialpy/utils/extraction.py
