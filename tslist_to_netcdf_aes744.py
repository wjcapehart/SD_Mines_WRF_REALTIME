#!/usr/bin/env python
# coding: utf-8

# ![AES-744 Masthead](https://kyrill.ias.sdsmt.edu/wjc/eduresources/AES_744_Masthead.png)
# 
# # WRF Time Series File to netCDF
# 
# Code converts WRF Timeseries ASCII Files to netCDF files.
# 
# The documentation on WRF's TimeSeries configuration is available here:
# 
# *  https://github.com/wrf-model/WRF/blob/master/run/README.tslist
# 
# ## Code Overview
# 
# *  Reads *tslist* station file
# *  Cross-references available output timeseries ASCII files and for each station determines the highest/innermost] domain.  This domain will be used to create the combined netCDF file. 
#     *  Write optional excel file for future access
# *  Extracts sigma layers coordinates
# *  Loops through station files
#     * Read individual files where ???? is the *tslist* station files station abbreviations (colunn 2), and dxx is the WRF domain from which the time series is extracted.
#         * single layer fields (*????.dxx.TS*)
#         * geopotential height profile (*????.dxx.PH*)
#         * potential temperature profile (*????.dxx.TH*)
#         * specific humidity profile (*????.dxx.QV*)
#         * single layer fields (*????.dxx.TS*)
#         * eastward wind speed profile (*????.dxx.UU*)
#         * northward wind speed profile (*????.dxx.VV*)
#         * vertical wind speed profile (*????.dxx.WW*)
#         * pressure profile (*????.dxx.PR*)   
#     * Output files into one CF-compliant netCDF timeSeries feature file per station
#         * file output nameform is: *wrfout_dxx_YYYY-MM-DD_HH_????".nc*, where YYYY-MM-DD_HH is the model start time.
#     
# 
# ## External Software Requirements
# 
# This notebook will require access to the following external packages
# 
# *  [The netCDF Command Operators (NCO) toolkit](https://nco.sourceforge.net): a unix-based tooklit that manipulates and analyzes data stored in netCDF-accessible formats
#    *   [ncatted](https://nco.sourceforge.net/nco.html#ncatted-netCDF-Attribute-Editor): The NCO netCDF Attribute Editor
# 

# ## Python Libraries
# 
# This script requires the following libraries.
# 
# * On-board Python Interfaces
#     *  [os](https://docs.python.org/3/library/os.html): Miscellaneous operating system interfaces
#     *  [glob](https://docs.python.org/3/library/glob.html): Unix style pathname pattern expansion
#     *  [datetime](https://docs.python.org/3/library/datetime.html): Basic date and time types
# *  External Libraries
#     *  [numpy](https://numpy.org): The fundamental package for scientific computing with Python
#     *  [xarray](https://xarray.dev): N-D labeled arrays and datasets in Python
#     *  [pandas](https://pandas.pydata.org): Flexible and powerful data analysis / manipulation library for Python, providing labeled data structures similar to R data.frame objects, statistical functions, and much more
#     *  [pint_xarray](https://pint-xarray.readthedocs.io/en/stable/): A convenience wrapper for using [pint](https://pint.readthedocs.io/en/stable/) in [xarray](https://xarray.dev) objects.
#     *  [metpy](https://unidata.github.io/MetPy/latest/index.html): A collection of tools in Python for reading, visualizing, and performing calculations with weather data
#     *  [haversine](https://github.com/mapado/haversine): calculates spherical line distance between two points (i.e., the [haversine distance](https://en.wikipedia.org/wiki/Haversine_formula)).
# 
# x
# 
# *  [platform](https://docs.python.org/3/library/platform.html): Access to underlying platform’s identifying data
# *  [socket](https://docs.python.org/3/library/socket.html): Low-level networking interface
# 

# In[11]:


####################################################
####################################################
####################################################
#
# Libraries
#

import numpy           as np
import xarray          as xr
import pandas          as pd

import pint_xarray     as px

import metpy.calc      as mpcalc
from   metpy.units import units

import haversine       as haversine

import os              as os
import glob            as glob
import datetime        as dt

#import platform          as platform
#import socket            as socket

#
####################################################
####################################################
####################################################


# ---
# ## Begin User Modification Area
# 
# Most general use of this code should only require manipulation of the following fields below
# 
# ### File Organization, Paths, and Filenames
# 
# ### Model Start Date String 
# 
# ### Latitude and Longitude of WRF Model Area of Interest
# 

# In[23]:


####################################################
####################################################
####################################################
#
# BEGIN USER MODIFICATION AREA
#

#
# File Organization, Paths, and Filenames
#

PATH_TO_WRF_OUTPUT_FILES = "/Users/wjc/GitHub/SD_Mines_WRF_REALTIME/WRF4/WRF/test/em_real/"
TSLIST_FILE              = "/Users/wjc/GitHub/SD_Mines_WRF_REALTIME/WRF4/WRF/test/em_real/tslist"

USE_STATION_EXCEL_FILE   = False
STATION_EXCEL_FILE       = "~/tslist_station_files.xlsx"

NETCDF_OUTPUT_PATH       = "~/."

#
# Model Start Date String (format must be YYYY_MM_DD_HH w/ "-"s and "_")
#

MODEL_START_DATE_YYYY_MM_DD_HH = "2023-03-28_18" # YYYY_MM_DD_HH

#
# Latitude and Longitude of WRF Model Area of Interest
#

TARGET_LATITUDE  =   44.074915 # degrees north
TARGET_LONGITUDE = -103.206571 # degrees east

#
# END USER MODIFICATION AREA
#
####################################################
####################################################
####################################################


# ## End User Modification Area
# ---
# 
# ## Process Model Start Date
# 
# String manipulation is applied to the user-provided model start date to create a date format matching the WRF Model's *wrfout* files' nameform.
# 
# The model start time is also converted into a Python Datetime Object using the [datetime.datetime.strptime()](https://docs.python.org/3/library/datetime.html#datetime.datetime.strptime)
# 

# In[19]:


####################################################
####################################################
####################################################
#
# Model Start Date
#

model_start_date_YYYY_MM_DD_HH0000 = MODEL_START_DATE_YYYY_MM_DD_HH + ":00:00"

print("WRFOUT Model Time Field: ", model_start_date_YYYY_MM_DD_HH0000)
    
    
    
model_start_datetime = dt.datetime.strptime(model_start_date_YYYY_MM_DD_HH0000, 
                                            '%Y-%m-%d_%H:%M:%S')

print("  Model Simulation Date: ", model_start_datetime)


#
####################################################
####################################################
####################################################


# ## Read tslist file
# 
# The *tslist* file contains the station location information for fields needed to extract.

# ### Pull Available Files for Each Domain.
# 
# The user can use the script to generate a pandas-readable MS Excel file with of the stations IDs, Station Long-Name (editable), their latitude & longitudes, and a distance of each station from a single reference point (e.g., the center point of the event under scrutiny).
# https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html
# 
# https://pandas.pydata.org/docs/reference/api/pandas.read_excel.html
# 
# https://pandas.pydata.org/docs/reference/api/pandas.read_fwf.html
# 
# https://docs.python.org/3/library/glob.html#glob.glob
# 
# 
# 
# ### Template for *tslist* Station File 
# 
# ```
# 
# 0000000000111111111122222222223333333333444444444
# 0123456789012345678901234567890123456789012345678
# #-----------------------------------------------#
# # 24 characters for name | pfx |  LAT  |   LON  |
# #-----------------------------------------------#
# RAPID CITY NWS SD         KUNR  44.0727 -103.211
# RAPID CITY ARPT SD        KRAP  44.0430 -103.054
# 
# 0000000000111111111122222222223333333333444444444
# 0123456789012345678901234567890123456789012345678
# 
# Fortran I/O Format : (A25,1X,A5,1X,F7,1X,F8)
# ```

# In[43]:


####################################################
####################################################
####################################################
#
# Read TSLIST File or Excel File derived from TSLIST
#

if (USE_STATION_EXCEL_FILE):

    print("Read file from "+station_list_file)
    
    available_time_series_list = pd.read_excel(STATION_EXCEL_FILE,
                                               index_col=0)
    
else:
    
    import haversine as haversine

    # Pull Full Time Series List
    
    print("Generating time series list from original wrf time series files")
    
    print("Accessing TSLIST: ", TSLIST_FILE)

    full_time_series_list = pd.read_fwf(filepath_or_buffer =    TSLIST_FILE, 
                                        header             =              2,
                                        colspecs           = [     [ 0,25],
                                                                   [26,31],
                                                                   [31,39],
                                                                   [39,-1] ],
                                        names              = ['Station Name',
                                                                'Station ID',
                                                                  'Latitude',
                                                                 'Longitude'])    
    
          
    # Grep the Library for *.TS
    
    display(len(full_time_series_list['Longitude'].values))

    available_time_series_list = glob.glob(pathname =            "????.d??.TS",
                                           root_dir = PATH_TO_WRF_OUTPUT_FILES)
    
    # Trim *.TS Suffix

    available_time_series_list = {x.replace('.TS', '') for x in available_time_series_list}

    # Convert to Data Frame

    available_time_series_list = pd.DataFrame(available_time_series_list, columns=['a'])

    # Convert to Split Domains and Station IDs

    available_time_series_list = available_time_series_list["a"].str.split(pat    = ".d", 
                                                                           n      =    1, 
                                                                           expand = True).rename(columns={0:"Station ID",
                                                                                                          1:    "Domain"})
    # Convert Domain to an integer

    available_time_series_list['Domain'] = available_time_series_list['Domain'].astype(int)

    # Capture only innermost available domain for each station

    available_time_series_list = available_time_series_list.groupby('Station ID').agg({'Domain':'max'}).reset_index()

    # Left Join Available and Full tables

    available_time_series_list = pd.merge(available_time_series_list,
                                          full_time_series_list,
                                          on  = 'Station ID',
                                          how = 'left').sort_values(by=['Domain','Station ID'])
    
    lats = available_time_series_list['Latitude'].values
    lons = available_time_series_list['Longitude'].values
    dist = lats.copy()
    for i in range(len(lats)):
        dist[i] = haversine.haversine([lats[i],lons[i]],[ TARGET_LATITUDE, TARGET_LONGITUDE])
                                  


    available_time_series_list["Distane from SDMines"] = dist

    available_time_series_list.sort_values(by = "Distane from SDMines", inplace=True)

    
    available_time_series_list.to_excel(STATION_EXCEL_FILE)

    print("Available Time Series List")
    display(available_time_series_list)

#
####################################################
####################################################
####################################################


# In[31]:





# In[ ]:





# ## Read common values for sigma levels and top of model space

# In[ ]:


####################################################
####################################################
####################################################
#
# Read common values for sigma levels and top of model space
#

model_start_date_YYYY_MM_DD_HH0000 = model_start_datetime.strftime('%Y-%m-%d_%H:00:00')

wrf_file  = WRF_EXE  + "./wrfout_d01_" + model_start_date_YYYY_MM_DD_HH0000
    
#
# Extract Sigma Coordinates
#

sigma = xr.DataArray(xr.open_dataset(wrf_file, engine="netcdf4")["ZNU"][0][0:15].values, 
                          name  =  "sigma",
                          dims  = ["sigma"],
                          attrs = {"description"   : "vertical sigma coordinates on mass points",
                                   "long_name"     : "vertical sigma coordinates on mass points",
                                   "standard_name" : "atmosphere_sigma_coordinate",
                                   "formula_terms" : "sigma: sigma ps: p_sfc ptop: wrf_ptop",
                                   "comment1"      : "pressure(n,k) = ptop + sigma(k)*(ps(k)-ptop)",
                                   "positive"      : "down"})

sigma = sigma.assign_coords({"sigma":sigma.values})
                        
wrf_ptop = xr.DataArray(xr.open_dataset(wrf_file)["P_TOP"].values,
                        name  = "wrf_ptop",
                             dims=["wrf_ptop"],
                             attrs = {"description"   : "Top-most Model Pressure",
                                      "long_name"     : "Top-most Model Pressure",
                                      "standard_name" : "air_pressure",
                                      "units"         : "Pa",
                                      "positive"      : "down"})

wrf_ptop = wrf_ptop.assign_coords({"wrf_ptop":wrf_ptop.values})

#
####################################################
####################################################
####################################################


# ## Loop through all stations

# In[2]:


####################################################
####################################################
####################################################
#
# Read TSLIST File or Excel File derived from TSLIST
#

#
# Loop through stations
#

for station in available_time_series_list.iterrows():
    station_id     = station[1][0]
    grid_domain    = station[1][1]
    station_name   = station[1][2]
    station_lat    = station[1][3]
    station_lon    = station[1][4]
    
    print("")  
    
    #
    #  Pull Grid-Relative Information
    #
    
    ts_file = WRF_EXE + station_id + ".d" + str(grid_domain).zfill(2) + ".TS"
    print(ts_file)
    
    myfile = open(ts_file, "r")
    headerline = myfile.readline()
    
    grid_i         =   int(headerline[ 58:62])
    grid_j         =   int(headerline[ 63:67])
    grid_lat       = float(headerline[ 70:77])
    grid_lon       = float(headerline[ 78:86])
    grid_elevation = float(headerline[ 87:94])
    grid_time      = datetime.datetime.strptime(headerline[103:122], '%Y-%m-%d_%H:%M:%S')
    units_time     = grid_time.strftime('%Y-%m-%d_%H:%M:%S.00')

    print(units_time)

    #
    #  Pull Time Series
    #
    
    ts_input = pd.read_fwf(ts_file,
                           header   = None,
                           skiprows = 1, delim_whitespace=True,
                           infer_nrows = 38880,
                           names    = ["id", "ts_hour", "id_tsloc", "ix", "iy", "t", "q", "u", "v", "psfc", "glw", "gsw", "hfx", "lfx", "tsk", "tslb", "rainc", "rainnc", "clw"],
                           width    = [   2,        13,          5,    5,    5,  14,  14,  14,  14,     14,    14,    14,    14,   14,    14,     14,      14,       14,    14])


    
    time = xr.DataArray(ts_input["ts_hour"].to_numpy(), 
                        name="time",
                        dims=["time"],
                        attrs = {"description"   : "time",
                                 "long_name"     : "time",
                                 "standard_name" : "time",
                                 "units": "hours since "+units_time})
    time = time.assign_coords({"time":time.values})

    

    t_air = xr.DataArray(ts_input["t"].to_numpy(), 
                         coords=[time],
                         dims  = ["time"],
                         name  =  "air_temperature_2m",
                         attrs = {"description"   : "2-m Air Temperatrue",
                                  "long_name"     : "2-m Air Temperatrue",
                                  "standard_name" : "air_temperature",
                                  "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude wrf_grid_elevation",
                                  "units"         : "K"})

    q_air = xr.DataArray(ts_input["q"].to_numpy(), 
                         coords=[time],
                         dims  = ["time"],
                         name  =  "specific_humidity_2m",
                         attrs = {"description"   : "2-m Specific Humidity",
                                  "long_name"     : "2-m Specific Humidity",
                                  "standard_name" : "specific_humidity",
                                  "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude wrf_grid_elevation",
                                  "units"         : "kg3 kg-3"})

    p_sfc = xr.DataArray(ts_input["psfc"].to_numpy(), 
                         coords=[time],
                         dims  = ["time"],
                         name  =  "surface_air_pressure",
                         attrs = {"description"   : "Surface Air Pressure",
                                  "long_name"     : "Surface Air Pressure",
                                  "standard_name" : "surface_air_pressure",
                                  "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude wrf_grid_elevation",
                                  "units"         : "Pa"})  
    
    u_10m = xr.DataArray(ts_input["u"].to_numpy(), 
                         coords=[time],
                         dims  = ["time"],
                         name  =  "eastward_wind_10m",
                         attrs = {"description"   : "10-m Eastward Wind",
                                  "long_name"     : "10-m Eastward Wind",
                                  "standard_name" : "eastward_wind",
                                  "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude wrf_grid_elevation",
                                  "units"         : "m s-1"})

    v_10m = xr.DataArray(ts_input["v"].to_numpy(), 
                         coords=[time],
                         name  =  "northward_wind_10m",
                         dims  = ["time"],
                         attrs = {"description"   : "10-m Northward Wind",
                                  "long_name"     : "10-m Northward Wind",
                                  "standard_name" : "northward_wind",
                                  "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude wrf_grid_elevation",
                                  "units"         : "m s-1"})
  

    glw   = xr.DataArray(ts_input["glw"].to_numpy(), 
                         coords=[time],
                         dims  = ["time"],
                         name  =  "surface_net_downward_longwave_flux",
                         attrs = {"description"   : "Surface Net Downward Longwave Flux",
                                  "long_name"     : "Surface Net Downward Longwave Flux",
                                  "standard_name" : "surface_net_downward_longwave_flux",
                                  "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude wrf_grid_elevation",
                                  "units"         : "W m-2"})    

    gsw   = xr.DataArray(ts_input["gsw"].to_numpy(), 
                         coords=[time],
                         dims  = ["time"],
                         name  =  "surface_net_downward_shortwave_flux",
                         attrs = {"description"   : "Surface Net Downward Shortwave Flux",
                                  "long_name"     : "Surface Net Downward Shortwave Flux",
                                  "standard_name" : "surface_net_downward_shortwave_flux",
                                  "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude wrf_grid_elevation",
                                  "units"         : "W m-2"})    

    hfx   = xr.DataArray(ts_input["hfx"].to_numpy(), 
                         coords=[time],
                         dims  = ["time"],
                         name  =  "surface_upward_sensible_heat_flux",
                         attrs = {"description"   : "Surface Net Upward Sensible Heat Flux",
                                  "long_name"     : "Surface Net Upward Sensible Heat Flux",
                                  "standard_name" : "surface_upward_sensible_heat_flux",
                                  "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude wrf_grid_elevation",
                                  "units"         : "W m-2"})    
    
    lfx   = xr.DataArray(ts_input["lfx"].to_numpy(), 
                         coords = [time],
                         dims   =  ["time"],
                         name   =  "surface_upward_latent_heat_flux",
                         attrs  =  {"description"   : "Surface Net Upward Latent Heat Flux",
                                  "long_name"     : "Surface Net Upward Latent Heat Flux",
                                  "standard_name" : "surface_upward_latent_heat_flux",
                                  "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude wrf_grid_elevation",
                                  "units"         : "W m-2"})   
    
    tsk   = xr.DataArray(ts_input["tsk"].to_numpy(), 
                         coords = [time],
                         dims   = ["time"],
                         name   =  "surface_temperature",
                         attrs  = {"description"  : "Surface Skin Temperature",
                                  "long_name"     : "Surface Skin Temperature",
                                  "standard_name" : "surface_temperature",
                                  "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude wrf_grid_elevation",
                                  "units"         : "K"})       
        
    tslb  = xr.DataArray(ts_input["tslb"].to_numpy(), 
                         coords = [time],
                         dims   = ["time"],
                         name   =  "topmost_soil_temperature",
                         attrs  =  {"description" : "Temperature of Top-Most Soil Layer",
                                  "long_name"     : "Temperature of Top-Most Soil Layer",
                                  "standard_name" : "soil_temperature",
                                  "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude wrf_grid_elevation",
                                  "units"         : "K"})       

    rainnc  = xr.DataArray(ts_input["rainnc"].to_numpy(), 
                         coords = [time],
                         dims   = ["time"],
                         name   =  "stratiform_precipitation_amount",
                         attrs  =  {"description"   : "Grid-Scale Precipitation Amount",
                                    "long_name"     : "Grid-Scale Precipitation Amount",
                                    "standard_name" : "stratiform_precipitation_amount",
                                    "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude wrf_grid_elevation",
                                    "units"         : "kg m-2"})       
 
    rainc  = xr.DataArray(ts_input["rainc"].to_numpy(), 
                         coords = [time],
                         dims   = ["time"],
                         name   =  "convective_precipitation_amount",
                         attrs  =  {"description"   : "Convective Precipitation Amount",
                                    "long_name"     : "Convective Precipitation Amount",
                                    "standard_name" : "convective_precipitation_amount",
                                    "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude wrf_grid_elevation",
                                    "units"         : "kg m-2"})     
           
    clw  = xr.DataArray(ts_input["clw"].to_numpy(), 
                         coords = [time],
                         dims   = ["time"],
                         name   =  "atmosphere_mass_content_of_water",
                         attrs  =  {"description"   : "Total Integrated Water",
                                    "long_name"     : "Total Integrated Water",
                                    "standard_name" : "atmosphere_mass_content_of_water",
                                    "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude wrf_grid_elevation",
                                    "units"         : "kg m-2"})     
   

    td_2m =mpcalc.dewpoint_from_specific_humidity(p_sfc, t_air, q_air).pint.to("K")


    td_2m  = xr.DataArray(td_2m, 
                         coords = [time],
                         dims   = ["time"],
                         name   =  "dew_point_temperature_2m",
                         attrs  =  {"description"   : "2-m Dew Point Temperatrue",
                                    "long_name"     : "2-m Dew Point Temperatrue",
                                    "standard_name" : "dew_point_temperature",
                                    "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude wrf_grid_elevation",
                                    "units"         : "K"})     

    #
    #  Pull Geopotential Height
    #
    
    file_2d = WRF_EXE + station_id + ".d" + str(grid_domain).zfill(2) + ".PH"

    colnames = ["ts_hour","0","1","2","3","4","5","6","7","8","9","10","11","12","13","14"]
    
    ts_input = pd.read_fwf(file_2d,
                           header      = None,
                           skiprows    = 1,
                           infer_nrows = 38880,
                           names       = colnames)
    
    ts_input['time'] = grid_time + ts_input["ts_hour"].to_numpy() * datetime.timedelta(hours=1)

    ts_input = ts_input.set_index('time')
    ts_input = ts_input.drop(columns = ["ts_hour"])
    
    ph = xr.DataArray(ts_input.to_numpy(), 
                      coords=[time,sigma],
                      name  =  "geopotential_height",
                      dims  = ["time","sigma"],
                      attrs = {"description"   : "geopotential height on mass points",
                               "long_name"     : "geopotential height on mass points",
                               "standard_name" : "geopotential_height",
                               "units"         : "m",
                               "positive"      : "down",
                               "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude sigma"})
    
    #
    #  Pull Potential Temperature 
    #
    
    file_2d = WRF_EXE + station_id + ".d" + str(grid_domain).zfill(2) + ".TH"

    colnames = ["ts_hour","0","1","2","3","4","5","6","7","8","9","10","11","12","13","14"]
    
    ts_input = pd.read_fwf(file_2d,
                           header      = None,
                           skiprows    = 1,
                           infer_nrows = 38880,
                           names       = colnames)
    
    ts_input['time'] = grid_time + ts_input["ts_hour"].to_numpy() * datetime.timedelta(hours=1)

    ts_input = ts_input.set_index('time')
    ts_input = ts_input.drop(columns = ["ts_hour"])
    
    th = xr.DataArray(ts_input.to_numpy(), 
                      coords=[time,sigma],
                      name  =  "air_potential_temperature",
                      dims  = ["time","sigma"],
                      attrs = {"description"   : "Potential Air Temperature",
                               "long_name"     : "Potential Air Temperature",
                               "standard_name" : "air_potential_temperature",
                               "units"         : "K",
                               "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude sigma"})

    #
    #  Pull Specific Humidity
    #
    
    file_2d = WRF_EXE + station_id + ".d" + str(grid_domain).zfill(2) + ".QV"

    colnames = ["ts_hour","0","1","2","3","4","5","6","7","8","9","10","11","12","13","14"]
    
    ts_input = pd.read_fwf(file_2d,
                           header      = None,
                           skiprows    = 1,
                           infer_nrows = 38880,
                           names       = colnames)
    
    ts_input['time'] = grid_time + ts_input["ts_hour"].to_numpy() * datetime.timedelta(hours=1)

    ts_input = ts_input.set_index('time')
    ts_input = ts_input.drop(columns = ["ts_hour"])
    
    
    qv = xr.DataArray(ts_input.to_numpy(), 
                      coords=[time,sigma],
                      name  =  "specific_humidity",
                      dims  = ["time","sigma"],
                      attrs = {"description"   : "Specific Humidity",
                               "long_name"     : "Specific Humidity",
                               "standard_name" : "specific_humidity",
                               "units"         : "kg3 kg-3",
                               "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude sigma"})
    
    #
    #  Pull  "Eastward Wind"
    #
    
    file_2d = WRF_EXE + station_id + ".d" + str(grid_domain).zfill(2) + ".UU"

    colnames = ["ts_hour","0","1","2","3","4","5","6","7","8","9","10","11","12","13","14"]
    
    ts_input = pd.read_fwf(file_2d,
                           header      = None,
                           skiprows    = 1,
                           infer_nrows = 38880,
                           names       = colnames)
    
    ts_input['time'] = grid_time + ts_input["ts_hour"].to_numpy() * datetime.timedelta(hours=1)

    ts_input = ts_input.set_index('time')
    ts_input = ts_input.drop(columns = ["ts_hour"])
    
    
    uu = xr.DataArray(ts_input.to_numpy(), 
                      coords=[time,sigma],
                      name  =  "eastward_wind",
                      dims  = ["time","sigma"],
                      attrs = {"description"   : "Eastward Wind",
                               "long_name"     : "Eastward Wind",
                               "standard_name" : "eastward_wind",
                               "units"         : "m s-1",
                               "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude sigma"})

    #
    #  Pull  "Nortward Wind"
    #
    
    file_2d = WRF_EXE + station_id + ".d" + str(grid_domain).zfill(2) + ".VV"

    colnames = ["ts_hour","0","1","2","3","4","5","6","7","8","9","10","11","12","13","14"]
    
    ts_input = pd.read_fwf(file_2d,
                           header      = None,
                           skiprows    = 1,
                           infer_nrows = 38880,
                           names       = colnames)
    
    ts_input['time'] = grid_time + ts_input["ts_hour"].to_numpy() * datetime.timedelta(hours=1)

    ts_input = ts_input.set_index('time')
    ts_input = ts_input.drop(columns = ["ts_hour"])
    
    
    vv = xr.DataArray(ts_input.to_numpy(), 
                      coords=[time,sigma],
                      name  =  "northward_wind",
                      dims  = ["time","sigma"],
                      attrs = {"description"   : "Northward Wind",
                               "long_name"     : "Northward Wind",
                               "standard_name" : "northward_wind",
                               "units"         : "m s-1",
                               "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude sigma"})

    #
    #  Pull  "Vertical Velocities 
    #
    
    file_2d = WRF_EXE + station_id + ".d" + str(grid_domain).zfill(2) + ".WW"

    colnames = ["ts_hour","0","1","2","3","4","5","6","7","8","9","10","11","12","13","14"]
    
    ts_input = pd.read_fwf(file_2d,
                           header      = None,
                           skiprows    = 1,
                           infer_nrows = 38880,
                           names       = colnames)
    
    ts_input['time'] = grid_time + ts_input["ts_hour"].to_numpy() * datetime.timedelta(hours=1)

    ts_input = ts_input.set_index('time')
    ts_input = ts_input.drop(columns = ["ts_hour"])
    
    ww = xr.DataArray(ts_input.to_numpy(), 
                      coords=[time,sigma],
                      name  =  "upward_air_velocity",
                      dims  = ["time","sigma"],
                      attrs = {"description"   : "Vertical Velocity",
                               "long_name"     : "Vertical Velocity",
                               "standard_name" : "upward_air_velocity",
                               "units"         : "m s-1",
                               "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude sigma"})
    
    #
    #  Pull  "Pressure 
    #
    
    file_2d = WRF_EXE + station_id + ".d" + str(grid_domain).zfill(2) + ".PR"

    colnames = ["ts_hour","0","1","2","3","4","5","6","7","8","9","10","11","12","13","14"]
    
    ts_input = pd.read_fwf(file_2d,
                           header      = None,
                           skiprows    = 1,
                           infer_nrows = 38880,
                           names       = colnames)
    
    ts_input['time'] = grid_time + ts_input["ts_hour"].to_numpy() * datetime.timedelta(hours=1)

    ts_input = ts_input.set_index('time')
    ts_input = ts_input.drop(columns = ["ts_hour"])
    
    
    pr = xr.DataArray(ts_input.to_numpy(), 
                      coords=[time,sigma],
                      name  =  "air_pressure",
                      dims  = ["time","sigma"],
                      attrs = {"description"   : "Air Pressure",
                               "long_name"     : "Air Pressure",
                               "standard_name" : "air_pressure",
                               "units"         : "Pa",
                               "coordinates"   : "time wrf_grid_latitude wrf_grid_longitude sigma"})

    #
    # Single Value Arrays
    #
    
    wrf_grid_longitude = xr.DataArray(np.array([grid_lon]),
                                 name  = "wrf_grid_longitude",
                                 dims=["wrf_grid_longitude"],
                                 attrs = {"description"   : "WRF Grid Longitude",
                                          "long_name"     : "WRF Grid Longitude",
                                          "standard_name" : "longitude",
                                          "units"         : "degrees_east"})
    
    wrf_grid_longitude = wrf_grid_longitude.assign_coords({"wrf_grid_longitude":wrf_grid_longitude.values})


    wrf_grid_latitude = xr.DataArray(np.array([grid_lat]),
                                 name  = "wrf_grid_latitude",
                                 dims=["wrf_grid_latitude"],
                                 attrs = {"description"   : "WRF Grid Latitude",
                                          "long_name"     : "WRF Grid Latitude",
                                          "standard_name" : "latitude",
                                          "units"         : "degrees_north"})
    wrf_grid_latitude = wrf_grid_latitude.assign_coords({"wrf_grid_latitude":wrf_grid_latitude.values})


    wrf_grid_elevation = xr.DataArray(np.array([grid_elevation]),
                                 name  = "wrf_grid_elevation",
                                 dims=["wrf_grid_elevation"],
                                 attrs = {"description"   : "WRF Grid Elevation",
                                          "long_name"     : "WRF Grid Elevation",
                                          "standard_name" : "elevation",
                                          "units"         : "m"})
    wrf_grid_elevation = wrf_grid_elevation.assign_coords({"wrf_grid_elevation":wrf_grid_elevation.values})


    
    #
    # Group Dataset for Export
    #

    time_series = xr.Dataset(coords = {"time": time,
                                      "sigma":sigma,
                                      "wrf_grid_latitude"  : wrf_grid_latitude,
                                      "wrf_grid_longitude" : wrf_grid_longitude,
                                      "wrf_grid_elevation" : wrf_grid_elevation},
                            data_vars = {"air_temperature_2m"                  : t_air,
                                          "specific_humidity_2m"                : q_air,
                                          "dew_point_temperature_2m"            : td_2m.as_numpy(),
                                          "eastward_wind_10m"                   : u_10m,
                                          "northward_wind_10m"                  : v_10m,
                                          "surface_air_pressure"                : p_sfc,
                                          "surface_net_downward_shortwave_flux" : gsw,
                                          "surface_net_downward_longwave_flux"  : glw,
                                          "surface_upward_sensible_heat_flux"   : hfx,
                                          "surface_upward_latent_heat_flux"     : lfx,
                                          "surface_temperature"                 : tsk,
                                          "topmost_soil_temperature"            : tslb,
                                          "stratiform_precipitation_amount"     : rainnc,
                                          "convective_precipitation_amount"     : rainc,
                                          "atmosphere_mass_content_of_water"    : clw,
                                          "geopotential_height"                 : ph,
                                          "air_potential_temperature"           : th,
                                          "specific_humidity"                   : qv,
                                          "eastward_wind"                       : uu,
                                          "northward_wind"                      : vv,
                                          "upward_air_velocity"                 : ww,
                                          "air_pressure"                        : pr},
                            attrs = {"featureType"                : "timeSeries",
                                     "Conventions"                : "CF-1.6",
                                     "Station_Name"               : station_name, 
                                     "Station_ID"                 : station_id,
                                     "Station_Longitude"          : station_lon,
                                     "Station_Latitude"           : station_lat,
                                     "WRF_Domain"                 : grid_domain,
                                     "WRF_Grid_Longitude"         : grid_lon,
                                     "WRF_Grid_Latitude"          : grid_lat,
                                     "WRF_Grid_I"                 : grid_i,
                                     "WRF_Grid_J"                 : grid_j,
                                     "WRF_Grid_Surface_Elevation" : grid_elevation})
 
    netcdf_file_name = NETCDF_OUTPUT_PATH + "./wrfout_d"+str(grid_domain).zfill(2)+"_"+ MODEL_START_DATE_YYYY_MM_DD_HH+"_"+station_id+".nc"
    
    print(netcdf_file_name)
    
    time_series.to_netcdf(path   = netcdf_file_name, 
                          mode   =              'w', 
                          format =        "NETCDF4")
    
    os.system("ncatted -Oh --attribute _FillValue,sigma,d,, "+netcdf_file_name)

print("")    
print("--------------------")

#
####################################################
####################################################
####################################################
    


# 

# In[ ]:





# In[ ]:





# ## Departing TSLIST_to_NETCDF

# In[ ]:


####################################################
####################################################
####################################################
#
# End of Script
#

print("tslist_to_netcdf complete.")

#
####################################################
####################################################
####################################################


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




