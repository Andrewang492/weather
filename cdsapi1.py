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

target = "gribs/cdsdata2.grib"

def download():
    dataset = "reanalysis-era5-land"
    request = {
        "variable": [        
            "total_precipitation",
            "10m_u_component_of_wind",
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
    client.retrieve(dataset, request, target) #.download()

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



def downloadCSV():
    # Open the GRIB file
    with open(target, 'rb') as f:
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

            # Create a CSV file
            #csv_file = target.replace('.grib', '.csv')
            short_name = codes_get(gid, 'shortName')
            csv_file = f'csvs/cdsdata_{short_name}.csv'
            with open(csv_file, 'w', newline='') as csvfile:
                csvwriter = csv.writer(csvfile)

                # Write the header
                csvwriter.writerow(['Latitude', 'Longitude', 'Value'])

                # Write the data
                for i in range(num_vals):
                    lat = lats[i % len(lats)]
                    lon = lons[i % len(lons)]
                    value = values[i]
                    csvwriter.writerow([lat, lon, value])

            print(f'{num_vals} values written to {csv_file}')

            # Release the GRIB file
            codes_release(gid)


def plot():
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
    #data, lons = shiftgrid(180., data, lons, start=False)
    grid_lon, grid_lat = np.meshgrid(lons, lats) #regularly spaced 2D grid
    
    m = Basemap(projection='cyl', llcrnrlon=-180, \
        urcrnrlon=180.,llcrnrlat=lats.min(),urcrnrlat=lats.max(), \
        resolution='c')
    
    x, y = m(grid_lon, grid_lat)
    
    cs = m.pcolormesh(x,y,data,shading='flat',cmap=plt.cm.gist_stern_r)
    
    m.drawcoastlines()
    m.drawmapboundary()
    m.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])
    
    plt.colorbar(cs,orientation='vertical', shrink=0.5)
    plt.title('CAMS AOD forecast') # Set the name of the variable to plot
    plt.savefig(grib+'.png') # Set the output file name




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process GRIB data.')
    parser.add_argument('--download', action='store_true', help='Download data')
    parser.add_argument('--read', action='store_true', help='Read data')
    parser.add_argument('--csv', type=str, help='Convert data to CSV')
    parser.add_argument('--plot', action='store_true', help='Plot data')

    args = parser.parse_args()

    if args.download:
        download()
    if args.read:
        read()
    if args.csv:
        downloadCSV(args.csv) #TODO add argument to downloadCSV
    if args.plot:
        plot()