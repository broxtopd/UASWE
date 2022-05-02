# Downscales 4 km UA Snowpack data to create 800 m maps of SWE and Snow Depth 
#
# USAGE: python Downscale.py InFiles OutDir overwrite 
	
# INPUTS: InFiles - a list of input UA SWE files
#         OutDir - directory where downscaled files will be written
#		  Overwrite - Flag whether or not to overwrite existing files
# OUTPUT: Creates a series of georeferenced rasters (tif) in the OutDir
#
# EXAMPLES:	python Downscale.py 4km/4km_SWE_Depth_WY2016_v01.nc 800m 1
# 			Converts the contents of 4km_SWE_Depth_WY2016_v01.nc to 800m SWE 
#			maps, with the overwrite flag set to 1
#			python Downscale.py 4km/4km_SWE_Depth_WY2015_v01.nc 800m,4km/4km_SWE_Depth_WY2016_v01.nc 800m 0
#			Converts the contents of 4km_SWE_Depth_WY2015_v01.nc and 
#			4km_SWE_Depth_WY2016_v01.nc to 800 m SWE maps (no overwrite)
# 
# Note: Download 4 km UA Snowpack NetCDF files from https://nsidc.org/data/nsidc-0719

###############################################################################
# Copyright (c) 2022, Patrick Broxton
# 
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
###############################################################################

import sys
import os
import netCDF4 as nc
from osgeo import gdal
import scipy.io as sio
import numpy as np
from datetime import datetime, timedelta
from PIL import Image

# Function to interpolate missing values 
def interpolate_missing_pixels(image: np.ndarray, mask: np.ndarray, method: str = 'nearest', fill_value: int = 0
):
    from scipy import interpolate

    h, w = image.shape[:2]
    xx, yy = np.meshgrid(np.arange(w), np.arange(h))

    known_x = xx[~mask]
    known_y = yy[~mask]
    known_v = image[~mask]
    missing_x = xx[mask]
    missing_y = yy[mask]

    interp_values = interpolate.griddata(
        (known_x, known_y), known_v, (missing_x, missing_y),
        method=method, fill_value=fill_value
    )

    interp_image = image.copy()
    interp_image[missing_y, missing_x] = interp_values

    return interp_image

if __name__ == "__main__":

	np.seterr(divide='ignore', invalid='ignore')

	InFiles_list = sys.argv[1]
	InFiles = InFiles_list.split(',')
	OutDir = sys.argv[2]
	overwrite = int(sys.argv[3])


	# Load Downscaling Factors
	md = sio.loadmat('Downscale_Data/Multipliers_800m.mat')
	ds = gdal.Open('Downscale_Data/mask_800m.tif')
	band = ds.GetRasterBand(1)
	mask = band.ReadAsArray()
	geotransform = ds.GetGeoTransform()
	projection = ds.GetProjection()
	ds = None

	# Define Grids
	resolution = 0.00833333333
	ulx = -125.0207999999999942
	uly = 49.9375000000000000
	lrx = -66.4791891046392891
	lry = 24.0625246505202881

	# New Grid
	prism_sz = [round((uly-lry) / resolution), round((lrx-ulx) / resolution)]
	dx = (lrx - ulx) / prism_sz[1]
	dy = (uly - lry) / prism_sz[0]
	lon_prism = np.reshape(np.arange(ulx + dx/2, lrx, dx), (-1, prism_sz[1]))
	lon_prism = np.tile(lon_prism, [prism_sz[0], 1])
	lat_prism = np.reshape(np.arange(uly - dy/2, lry, -dy), (-1, prism_sz[0])).T
	lat_prism = np.tile(lat_prism, [1, prism_sz[1]])

	# Original Grid
	resolution_orig = 0.04166666667
	prism_sz_orig = [round((uly-lry) / resolution_orig), round((lrx-ulx) / resolution_orig)]
	dx_orig = (lrx - ulx) / prism_sz_orig[1]
	dy_orig = (uly - lry) / prism_sz_orig[0]
	lon_prism_orig = np.reshape(np.arange(ulx + dx_orig/2, lrx, dx_orig), (-1, prism_sz_orig[1]))
	lon_prism_orig = np.tile(lon_prism_orig, [prism_sz_orig[0], 1])
	lat_prism_orig = np.reshape(np.arange(uly - dy_orig/2, lry, -dy_orig), (-1, prism_sz_orig[0])).T
	lat_prism_orig = np.tile(lat_prism_orig, [1, prism_sz_orig[1]])

	for i in range(len(InFiles)):
		print('Reading data for ' + InFiles[i])

		# Read in the 4 km data
		ds = nc.Dataset(InFiles[i])
		for k in ds.variables:
			ds.variables[k].set_auto_mask(False)

		time = ds['time'][:]
		SWE = ds['SWE'][:]
		DEPTH = ds['DEPTH'][:]


		for d in range(len(SWE[:, 1, 1])):

			# Keep track of dates, etc
			ts = datetime(1900, 1, 1) + timedelta(float(time[d]))
			yyyy = str(ts.year)
			mm = str(ts.month)
			if len(mm) < 2:
				mm = '0' + mm
			dd = str(ts.day)
			if len(dd) < 2:
				dd = '0' + dd
			doy = ts.timetuple().tm_yday

			if (not os.path.exists(OutDir + '/' + yyyy + '/' + mm + '/' + dd + '/SWE.tif') and not os.path.exists(
					OutDir + '/' + yyyy + '/' + mm + '/' + dd + '/Depth.tif')) or overwrite:

				print('Downscaling Data for ' + ts.strftime("%b-%d-%Y"))

				# Figure out how to weight seasonal SWE multiplier maps
				frac_spring = max(0, min(1, (doy - 0) / (91 - 0)))
				frac_summer = max(0, min(1, (doy - 91) / (182 - 91)))
				frac_spring = frac_spring - frac_summer
				frac_fall = max(0, min(1, (doy - 182) / (274 - 182)))
				frac_summer = frac_summer - frac_fall
				frac_winter = max(0, min(1, (doy - 274) / (365 - 274)))
				frac_fall = frac_fall - frac_winter
				frac_winter = 1 - (frac_spring + frac_summer + frac_fall)
				SWEMult = frac_winter * md['SWEMult_winter'] + frac_spring * md['SWEMult_spring'] + frac_summer * md['SWEMult_summer'] + frac_fall * md['SWEMult_fall']

				# Extract the 4 km data for this date

				swe_4km = np.flipud(SWE[d, :, :].astype('float'))
				swe_4km[swe_4km == -999.] = np.nan
				depth_4km = np.flipud(DEPTH[d, :, :].astype('float'))
				depth_4km[depth_4km == -999.] = np.nan
				density_4km = swe_4km / depth_4km
				density_4km[depth_4km == 0] = np.nan

				# Downscale SWE data
				if np.nanmax(swe_4km) > 0:
					nanlocs = np.isnan(swe_4km)
					# To prevent missing data near edges 
					swe_4km_filled = interpolate_missing_pixels(image=swe_4km, mask=nanlocs)
				else:
					swe_4km_filled = np.zeros(swe_4km.shape)

				swe_800m = np.array(Image.fromarray(swe_4km_filled).resize(size=(lat_prism.shape[1], lat_prism.shape[0]))) * SWEMult
				swe_800m[swe_800m < 1] = 0
				swe_800m[mask > 0] = np.nan

				# Downscale depth data (through snow density)
				if np.nanmax(swe_4km) > 0:
					nanlocs = np.isnan(density_4km)
					density_4km_filled = interpolate_missing_pixels(image=density_4km, mask=nanlocs)
				else:
					density_4km_filled = np.ones(density_4km.shape) * 0.2

				density_800m = np.array(Image.fromarray(density_4km_filled).resize(size=(lat_prism.shape[1], lat_prism.shape[0])))
				density_800m[swe_800m == 0 + np.isnan(swe_800m)] = 0.2
				density_800m[density_800m > 0.5] = 0.5
				density_800m[density_800m < 0.1] = 0.1
				depth_800m = swe_800m / density_800m

				swe_800m[mask > 0] = -999
				depth_800m[mask > 0] = -999

				# Output Files
				if not os.path.exists(OutDir + '/' + yyyy + '/' + mm + '/' + dd):
					os.makedirs(OutDir + '/' + yyyy + '/' + mm + '/' + dd)

				[rows, cols] = swe_800m.shape
				driver = gdal.GetDriverByName("GTiff")
				co = ["COMPRESS=DEFLATE"]
				outdata = driver.Create(OutDir + '/' + yyyy + '/' + mm + '/' + dd + '/SWE.tif', cols, rows, 1, gdal.GDT_Float32, options=co)
				outdata.SetGeoTransform(geotransform)
				outdata.SetProjection(projection)
				outdata.GetRasterBand(1).WriteArray(swe_800m)
				outdata.GetRasterBand(1).SetNoDataValue(-999)
				outdata.FlushCache()
				outdata = None

				[rows, cols] = depth_800m.shape
				driver = gdal.GetDriverByName("GTiff")
				co = ["COMPRESS=DEFLATE"]
				outdata = driver.Create(OutDir + '/' + yyyy + '/' + mm + '/' + dd + '/Depth.tif', cols, rows, 1, gdal.GDT_Float32, options=co)
				outdata.SetGeoTransform(geotransform)
				outdata.SetProjection(projection)
				outdata.GetRasterBand(1).WriteArray(depth_800m)
				outdata.GetRasterBand(1).SetNoDataValue(-999)
				outdata.FlushCache()
				outdata = None
