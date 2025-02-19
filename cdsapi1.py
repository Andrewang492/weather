import cdsapi
from eccodes import *

import pygrib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import shiftgrid
import numpy as np

import csv

import argparse

import cfgrib
import xarray as xr

from netCDF4 import Dataset as NetCDFFile 
#import matplotlib.pyplot as plt
#import numpy as np
#from mpl_toolkits.basemap import Basemap

from request_params import request

target = "gribs/cdsdata2.grib" #TODO phase this out

def download(outputFilePath):
    dataset = "reanalysis-era5-land"
    request = {
        "variable": [        
            "10m_u_component_of_wind",
            "10m_v_component_of_wind",
            "surface_pressure",
            "total_precipitation"
        ],
        "year": "2025",
        "month": "02",
        "day": ["01"],
        "time": ["00:00"],
        "data_format": "grib",
        "download_format": "unarchived",
        "area": [-23, 141, -39, 154]
    }

    client = cdsapi.Client()
    client.retrieve(dataset, request, outputFilePath) #.download()

def read():
    # ----- Now read -----
    f = open(target, 'rb')
    gid = codes_grib_new_from_file(f)

    values = codes_get_values(gid)
    num_vals = len(values)
    for i in range(num_vals):
        print("%d %.10e" % (i + 1, values[i]))

    print('%d values found in %s' % (num_vals, target))

    for key in ('max', 'min', 'average'):
        print('%s=%.10e' % (key, codes_get(gid, key)))

    # Example of accessing specific elements from data values
    # Get first, middle and last elements
    indexes = [0, int(num_vals/2), num_vals-1]
    elems = codes_get_double_elements(gid, 'values', indexes)

    codes_release(gid)
    f.close()

def downloadCSV(gribFilePath, outputFilePath):
    # Open the GRIB file
    with open(gribFilePath, 'rb') as f:
        # Create a CSV file
        with open(outputFilePath, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)

            # Initialize header
            header = ['Latitude', 'Longitude']
            values_list = []

            while True:
                gid = codes_grib_new_from_file(f)
                if gid is None:
                    break

                print(f'gid: {gid}')
                # Get the values from the GRIB file
                values = codes_get_values(gid)
                num_vals = len(values)

                # Get the latitudes and longitudes
                lats = codes_get_double_array(gid, 'distinctLatitudes')
                lons = codes_get_double_array(gid, 'distinctLongitudes')

                # Add a new Value column to the header
                short_name = codes_get(gid, 'shortName')
                header.append(f'Value_{short_name}')
                values_list.append(values)

                # Release the GRIB file
                codes_release(gid)

            # Write the header
            csvwriter.writerow(header)

            # Write the data
            for i in range(num_vals):
                row = [lats[i % len(lats)], lons[i % len(lons)]]
                for values in values_list:
                    row.append(values[i])
                csvwriter.writerow(row)

            print(f'{num_vals} values written to {outputFilePath}')

def toNetCDF():
    ds = cfgrib.open_file(target)

    ds = xr.open_dataset(target, engine='cfgrib')
    ds.to_netcdf("ncs/era5.nc")

def plot3():
    # Plotting
    
    plt.figure(figsize=(12,8))
    
    grib = target # Set the file name of your input GRIB file
    grbs = pygrib.open(grib)
    
    grb = grbs.select()[0]
    data = grb.values
    
    # need to shift data grid longitudes from (0..360) to (-180..180)
    lons = np.linspace(float(grb['longitudeOfFirstGridPointInDegrees']), \
        float(grb['longitudeOfLastGridPointInDegrees']), int(grb['Ni']) )
    lats = np.linspace(float(grb['latitudeOfFirstGridPointInDegrees']), \
        float(grb['latitudeOfLastGridPointInDegrees']), int(grb['Nj']) )
    # data, lons = shiftgrid(180., data, lons, start=False)
    data, lons = shiftgrid(lons.max(), data, lons, start=False)
    grid_lon, grid_lat = np.meshgrid(lons, lats) #regularly spaced 2D grid
    
    print(grid_lon, grid_lat)

    m = Basemap(projection='cyl', llcrnrlon=lons.min(), \
        urcrnrlon=lons.max(),llcrnrlat=lats.min(),urcrnrlat=lats.max(), \
        resolution='c')
    
    x, y = m(grid_lon, grid_lat)
    
    cs = m.pcolormesh(x,y,data,shading='flat',cmap=plt.cm.gist_stern_r)
    
    m.drawcoastlines()
    m.drawmapboundary()
    m.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])
    
   # plt.colorbar(cs,orientation='vertical', shrink=0.5)
    plt.title('CAMS AOD forecast') # Set the name of the variable to plot
    #plt.savefig(grib+'.png') # Set the output file name
    plt.savefig('plots/basemap.png') # Set the output file name

def plot2():
    ds = cfgrib.open_file(target)

    ds = xr.open_dataset(target, engine='cfgrib')
    precip = ds.variables['tp']
    precipitation_da = xr.DataArray(precip)
    print(precipitation_da)
    
    fig = plt.figure(figsize=(10,6))
    precipitation_da.plot(
        cmap="Blues",
        cbar_kwargs={"label": "Total Precipitation [mm]"},
    )
    plt.savefig("plots/plot.png") # Set the output file name

def plot():
    # Derived from: https://www2.atmos.umd.edu/~cmartin/python/examples/netcdf_example1.html
    nc = NetCDFFile('ncs/era5.nc')

    lat = nc.variables['latitude'][:]
    lon = nc.variables['longitude'][:]
    time = nc.variables['time'][:]
    t2 = nc.variables['tp'][:] # precipitation
    #u = nc.variables['10u'][:] # 10m u-component of winds

    #Let's create a map centered over ___
    map = Basemap(projection='merc',llcrnrlon=141.,llcrnrlat=-39.,urcrnrlon=154.,urcrnrlat=-23.,resolution='i') # projection, lat/lon extents and resolution of polygons to draw
    # resolutions: c - crude, l - low, i - intermediate, h - high, f - full

    map.drawcoastlines()
    #map.drawrivers(color='lightblue')
    map.drawstates()
    map.drawcountries()
    map.drawlsmask(land_color='Linen', ocean_color='#CCFFFF') # can use HTML names or codes for colors

    parallels = np.arange(-20,-40,5.) # make latitude lines ever 5 degrees from 30N-50N
    meridians = np.arange(140, 155,5.) # make longitude lines every 5 degrees from 95W to 70W
    map.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

    lons,lats= np.meshgrid(lon,lat) #do not need to shift 180.
    x,y = map(lons,lats)

    plt.title('Precipitation (mm)')
    print(t2)
    print(f"Dimensions of t2: {t2.shape}")
    print(t2[:, :])
    temp = map.contourf(x, y, t2[:, :], cmap='coolwarm')
    cb = map.colorbar(temp,"right", size="5%", pad="2%")
    plt.title('Precip')
    cb.set_label('precipitation (mm)')
    
    plt.savefig('plots/plot3.png')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process GRIB data.')
    parser.add_argument('--download', metavar='OUTPUT_FILE', help='Download data to the specified file')
    parser.add_argument('--read', action='store_true', help='Read data')
    parser.add_argument('--csv', nargs=2, metavar=('GRIB_FILE', 'CSV_FILE'), help='Convert GRIB file to CSV')
    parser.add_argument('--netcdf', action='store_true', help='Convert GRIB file to NetCDF')
    parser.add_argument('--plot', action='store_true', help='Plot data')

    args = parser.parse_args()

    if args.download:
        download(args.download)
    if args.read:
        read()
    if args.csv:
        gribFilePath, outputFilePath = args.csv
        downloadCSV(gribFilePath, outputFilePath)
    if args.netcdf:
        toNetCDF()
    if args.plot:
        plot()