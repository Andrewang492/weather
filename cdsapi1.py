from cdsrequest import request # Your request input.
import cdsapi
from eccodes import *
import csv
import cfgrib
import xarray as xr
from netCDF4 import Dataset as NetCDFFile 
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import argparse


def download(outputFileName):
    dataset = "reanalysis-era5-land"
    client = cdsapi.Client()
    client.retrieve(dataset, request, f'gribs/{outputFileName}.grib') #.download()
    print(f'GRIB file saved to gribs/{outputFileName}.grib')

def downloadCSV(gribFileName, outputFileName):
    # Open the GRIB file
    with open(f'gribs/{gribFileName}.grib', 'rb') as f:
        # Create a CSV file
        with open(f'csvs/{outputFileName}.csv', 'w', newline='') as csvfile:
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

            print(f'{num_vals} values written to csvs/{outputFileName}.csv')

def toNetCDF(gribFileName, outputFileName):
    ds = cfgrib.open_file(f'gribs/{gribFileName}.grib')

    ds = xr.open_dataset(f'gribs/{gribFileName}.grib', engine='cfgrib')
    ds.to_netcdf(f'ncs/{outputFileName}.nc')
    print(f'NetCDF file saved to ncs/{outputFileName}.nc')

def plot(netcdfFileName, variableShortName, outputFileName ):
    # Derived from: https://www2.atmos.umd.edu/~cmartin/python/examples/netcdf_example1.html
    nc = NetCDFFile(f'ncs/{netcdfFileName}.nc')
    print(nc.variables.keys())
    lat = nc.variables['latitude'][:]
    lon = nc.variables['longitude'][:]
    time = nc.variables['time'][:]
    variable = nc.variables[variableShortName][:]

    #Let's create a map centered over ___
    latLower = min(lat)
    latUpper = max(lat)
    lonLower = min(lon)
    lonUpper = max(lon)

    map = Basemap(projection='merc',llcrnrlon=lonLower,llcrnrlat=latLower,urcrnrlon=lonUpper,urcrnrlat=latUpper,resolution='i') # projection, lat/lon extents and resolution of polygons to draw
    # Draw lines for background
    map.drawcoastlines()
    #map.drawrivers(color='lightblue')
    map.drawstates()
    map.drawcountries()
    map.drawlsmask(land_color='Linen', ocean_color='#CCFFFF')
    # Draw parallels and meridians
    parallels = np.arange(-90, 90,5.)
    meridians = np.arange(-180, 180,5.) 
    map.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

    lons,lats= np.meshgrid(lon,lat) #do not need to shift 180.
    x,y = map(lons,lats)

    plt.title(f'{variableShortName} (u)') #TODO make this dynamic
    temp = map.contourf(x, y, variable, cmap='coolwarm')
    cb = map.colorbar(temp,"right", size="5%", pad="2%")
    plt.title(f'{variableShortName} (u)')
    cb.set_label(f'{variableShortName} (u)')
    
    plt.savefig(f'plots/{outputFileName}.png', format='png', dpi=300)
    print(f'Plot saved to plots/{outputFileName}.png')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process GRIB data.')
    parser.add_argument('--download', metavar='OUTPUT_FILE', help='Download data to the specified file')
    parser.add_argument('--csv', nargs=2, metavar=('GRIB_FILE', 'CSV_FILE'), help='Convert GRIB file to CSV')
    parser.add_argument('--netcdf', nargs=2, metavar=('GRIB_FILE', 'NETCDF_FILE'), help='Convert GRIB file to NetCDF')
    parser.add_argument('--plot', nargs=3, metavar=('NETCDF_FILE', 'VAR_SHORT_NAME', 'PNG_FILE'), help='Plot data')

    args = parser.parse_args()

    if args.download:
        download(args.download)
    if args.csv:
        gribFileName, outputFileName = args.csv
        downloadCSV(gribFileName, outputFileName)
    if args.netcdf:
        gribFileName, outputFileName = args.netcdf
        toNetCDF(gribFileName, outputFileName)
    if args.plot:
        netcdfFileName, variableShortName, outputFileName = args.plot
        #print(variableShortName)

        plot(netcdfFileName, variableShortName, outputFileName)