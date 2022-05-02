
Downscales 4 km UA Snowpack data to create 800 m maps of SWE and Snow Depth 

USAGE: python Downscale.py InFiles OutDir overwrite 
	
INPUTS: InFiles - a list of input UA SWE files
        OutDir - directory where downscaled files will be written
        Overwrite - Flag whether or not to overwrite existing files
OUTPUT: Creates a series of georeferenced rasters (tif) in the OutDir

EXAMPLES:	python Downscale.py 4km/4km_SWE_Depth_WY2016_v01.nc 800m 1
			Converts the contents of 4km_SWE_Depth_WY2016_v01.nc to 800m SWE 
			maps, with the overwrite flag set to 1
			python Downscale.py 4km/4km_SWE_Depth_WY2015_v01.nc,4km/4km_SWE_Depth_WY2016_v01.nc 800m 0
			Converts the contents of 4km_SWE_Depth_WY2015_v01.nc and 
			4km_SWE_Depth_WY2016_v01.nc to 800 m SWE maps (no overwrite)
 
Note: Download 4 km UA Snowpack NetCDF files from https://nsidc.org/data/nsidc-0719
