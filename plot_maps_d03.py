#!/usr/bin/env python
# coding: utf-8

# # Plotting Maps

# ## Libraries

# In[ ]:


####################################################
####################################################
####################################################
#
# Libraries
#

import numpy             as np
import datetime          as datetime
import os                as os
import platform          as platform


import matplotlib        as mpl
import matplotlib.pyplot as plt
import pandas            as pd
import xarray            as xr
import pint_xarray       as px

import netCDF4           as nc4

import wrf               as wrf

import cartopy.crs       as ccrs
import cartopy.feature   as cfeature

import matplotlib.gridspec as gridspec

import metpy.calc        as mpcalc
from metpy.units import units, pandas_dataframe_to_unit_arrays
from metpy.plots import colortables,  USCOUNTIES


import metpy.calc  as mpcalc

from metpy.units import units

import seaborn           as sns

import timezonefinder    as tzf
import pytz as pytz

import socket as socket
import matplotlib.font_manager as fm
import matplotlib as mpl



#
####################################################
####################################################
####################################################


# ## File Organization

# In[ ]:


####################################################
####################################################
####################################################
#
# Directory Workspaces
#

beta_on       = 0
max_domains   = 3
chosen_domain = 3

if (socket.gethostname() == "kyrill"):
    WRF_OVERALL_DIR = "/projects/SD_Mines_WRF_REALTIME/"
else:
    if (platform.system() == "Darwin"):
         WRF_OVERALL_DIR = "/Users/wjc/GitHub/SD_Mines_WRF_REALTIME/"
    else:
         WRF_OVERALL_DIR = "/home/wjc/GitHub/SD_Mines_WRF_REALTIME/"

        



os.chdir(WRF_OVERALL_DIR)



print( "Current Working Directory is now " + os.getcwd() )

WPS_WORK    = WRF_OVERALL_DIR + "./WPS_PrepArea/"
WPS_EXE     = WRF_OVERALL_DIR + "./WRF4/WPS/"
WRF_EXE     = WRF_OVERALL_DIR + "./WRF4/WRF/test/em_real/"
WRF_ARCHIVE = WRF_OVERALL_DIR + "./ARCHIVE/"
WRF_IMAGES  = WRF_OVERALL_DIR + "./WEB_IMAGES/"

NCEP_FTP_URLROOT       = "ftp://ftpprd.ncep.noaa.gov/pub/data/nccf/com/nam/prod/nam."

NCEP_FTP_SERVER        = "ftpprd.ncep.noaa.gov"

NCEP_FTP_GRIB_DIR_ROOT = "/pub/data/nccf/com/nam/prod/nam."

#
####################################################
####################################################
####################################################
os.chdir(WRF_OVERALL_DIR)

print( "Current Working Directory is now " + os.getcwd() )
    
WPS_WORK    = WRF_OVERALL_DIR + "./WPS_PrepArea/"
WPS_EXE     = WRF_OVERALL_DIR + "./WRF4/WPS/"
WRF_EXE     = WRF_OVERALL_DIR + "./WRF4/WRF/test/em_real/"
WRF_ARCHIVE = WRF_OVERALL_DIR + "./ARCHIVE/"
WRF_IMAGES  = WRF_OVERALL_DIR + "./WEB_IMAGES/"

TS_DIR = WRF_EXE

station_list_file = WRF_OVERALL_DIR + "namelist_files_and_local_scripts/time_series_station_files_"+str(max_domains)+"_dom.xlsx"

os.chdir(WRF_EXE)

#
####################################################
####################################################
####################################################


# ## Time Control

# In[ ]:


####################################################
####################################################
####################################################
#
# Model Start Date
#

with open(WRF_OVERALL_DIR + "./current_run.txt") as f:
    model_start_date_YYYY_MM_DD_HH = f.readlines()

model_start_date_YYYY_MM_DD_HH     = model_start_date_YYYY_MM_DD_HH[0][0:13]

model_start_date_YYYY_MM_DD_HH0000 = model_start_date_YYYY_MM_DD_HH + ":00:00"
print(model_start_date_YYYY_MM_DD_HH0000)
    
model_start_datetime = datetime.datetime.strptime(model_start_date_YYYY_MM_DD_HH0000, '%Y-%m-%d_%H:%M:%S')
print("Model Simulation Date ", model_start_datetime)

model_end_datetime  = model_start_datetime + datetime.timedelta(hours=36)
current_datetime    = datetime.datetime.utcnow()
siphon_end_datetime = min(current_datetime,model_end_datetime)

print( "Current Working Directory is now " + os.getcwd() )
print( "         Model Start Datetime is " + model_start_datetime.strftime("%Y-%m-%d %H:00:00"))
print( "           Model End Datetime is " +   model_end_datetime.strftime("%Y-%m-%d %H:00:00"))
print( "             Current Datetime is " +     current_datetime.strftime("%Y-%m-%d %H:00:00"))
print( "          Siphon End Datetime is " +  siphon_end_datetime.strftime("%Y-%m-%d %H:00:00"))
print( "               Station List File " +    station_list_file)

wrf_skewt_time    = model_start_datetime.strftime("%Y-%m-%d %H UTC")




tf     = tzf.TimezoneFinder()
tz     = tf.certain_timezone_at(lng=-104, lat=44)

tzabbr = pytz.timezone(tz).localize(model_start_datetime)



#
####################################################
####################################################
####################################################


# In[ ]:


###################################################
#
# NWS Rainfall Color Table.
#

nws_precip_colors = [
    "#04e9e7",  # 0.01 - 0.10 inches
    "#019ff4",  # 0.10 - 0.25 inches
    "#0300f4",  # 0.25 - 0.50 inches
    "#02fd02",  # 0.50 - 0.75 inches
    "#01c501",  # 0.75 - 1.00 inches
    "#008e00",  # 1.00 - 1.50 inches
    "#fdf802",  # 1.50 - 2.00 inches
    "#e5bc00",  # 2.00 - 2.50 inches
    "#fd9500",  # 2.50 - 3.00 inches
    "#fd0000",  # 3.00 - 4.00 inches
    "#d40000",  # 4.00 - 5.00 inches
    "#bc0000",  # 5.00 - 6.00 inches
    "#f800fd",  # 6.00 - 8.00 inches
    "#9854c6",  # 8.00 - 10.00 inches
    "#fdfdfd"]  # 10.00+

precip_colormap = mpl.colors.ListedColormap(colors = nws_precip_colors)

precip_levels_in = [   0.01,   0.10,  0.25,   0.50, 
                       0.75,   1.00,  1.50,   2.00, 
                       2.50,   3.00,  4.00,   5.00,
                       6.00,   8.00, 10.00,  20.00] # in Inches!!!

precip_levels_mm = [  0.25,   2.50,   5.00,  10.00, 
                     20.00,  25.00,  40.00,  50.00, 
                     60.00,  75.00, 100.00, 125.00,
                    150.00, 200.00, 250.00, 500.00] # in mm

#
###################################################


# In[ ]:


###################################################
#
# BOM Temperature Color Table.
#

bom_temp_colors = ['#CFEAFD',
                   '#B8CCE3',
                   '#8D99BF',
                   '#8D97C2',
                   '#835BA1',
                   '#752078',
                   '#0F1D53',
                   '#172974',
                   '#28368E',
                   '#2B4999',
                   '#2F5CA3',
                   '#3C75B0',
                   '#478DBC',
                   '#51A6C2',
                   '#61B2C1',
                   '#7CC7C1',
                   '#91CBBD',
                   '#A6D9A8',
                   '#CDE9B7',
                   '#E0F0CA',
                   '#F0F8DE',
                   '#F6FBD8',
                   '#FCFED1',
                   '#FCFABD',
                   '#FAEDAA',
                   '#FAE495',
                   '#F6DA85',
                   '#F5C671',
                   '#F1B460',
                   '#EDA255',
                   '#EF914D',
                   '#EE7444',
                   '#E95B3C',
                   '#DF4431',
                   '#AF2B2F',
                   '#C1322D',
                   '#B02C2F',
                   '#93232D',
                   '#6A1721',
                   '#450B17',
                   '#17020A',
                   '#441293',
                   '#BE34CB',
                   '#FFFFFF']
 
temperature_colormap = mpl.colors.ListedColormap(colors = bom_temp_colors)

temperature_levels_degC = np.linspace(-32,54,44) # in degC
temperature_levels_degF = temperature_levels_degC * 9./5. + 32.


#
###################################################


# In[ ]:


####################################################
####################################################
####################################################
#
# MetPy Color Tables
#

nws_dbz_colors = colortables.get_colortable("NWSReflectivity")
nws_rainbow    = colortables.get_colortable("ir_rgbv")

clearsky_dbz_values = np.arange(-28, 28.1, 2)
stormy_dbz_values   = np.arange(  5, 75.1, 5)

#
####################################################
####################################################
####################################################


# ## Crack WRF Files

# In[ ]:


####################################################
####################################################
####################################################
#
# Rotate through Available Files
#

for domain in range(chosen_domain,chosen_domain+1):
    
    
    if (domain == 1): 
        figure_domain_size = (7,6)
    elif (domain ==2):
        figure_domain_size = (7,5.25)
    else:
        figure_domain_size = (7,6)
        
        

    
    graphics_directory = WRF_IMAGES + "/" + model_start_date_YYYY_MM_DD_HH + "/MAPS/d" +  str(domain).zfill(2) + "/"
    
    os.system("mkdir -pv " + graphics_directory + "SFCT")
    os.system("mkdir -pv " + graphics_directory + "WIND")
    
    os.system("mkdir -pv " + graphics_directory + "PBL")
    os.system("mkdir -pv " + graphics_directory + "DBZ")


    os.system("mkdir -pv " + graphics_directory + "SNOWH")

    os.system("mkdir -pv " + graphics_directory + "RAIN")
    os.system("mkdir -pv " + graphics_directory + "WEASD")

    os.system("mkdir -pv " + graphics_directory + "TRAIN")
    os.system("mkdir -pv " + graphics_directory + "TSNOW")

    wrf_file  = WRF_EXE  + "./wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH0000
    print(wrf_file)
    
    ncf = nc4.Dataset(filename = wrf_file)
    
    wrf_time_steps     = wrf.getvar(wrfin    =           ncf,
                                    varname  =       'times',
                                    timeidx  = wrf.ALL_TIMES)  
    
    nt                 = len(wrf_time_steps)
    
    #
    # Precipitation
    #
    
    rainc_maps         = wrf.getvar(wrfin    =           ncf,
                                    varname  =          'RAINC',
                                    timeidx  = wrf.ALL_TIMES) 
    
    rainnc_maps        = wrf.getvar(wrfin    =           ncf,
                                    varname  =          'RAINNC',
                                    timeidx  = wrf.ALL_TIMES)     
    
    rainsc_maps        = wrf.getvar(wrfin    =           ncf,
                                    varname  =          'RAINSH',
                                    timeidx  = wrf.ALL_TIMES)    
    
    ##########################################################
    # 
    # Pulling at the CRS Proj-4 fields.  
    #

    cart_proj = wrf.get_cartopy(wrfin = ncf)

    # Display the Proj-4 Information

    print("-- Proj4 Settings --")
    print("")
    #display("cart_proj.proj4_params",cart_proj.proj4_params)

    # Display the Flat Earth Orange Peel.  

    print(cart_proj)

    #
    ##########################################################

    rain_maps          = rainc_maps.copy()
    rain_maps.values   = rainc_maps.values + rainnc_maps.values + rainsc_maps.values
    
    rain_maps.values   = rain_maps.values / 25.4
    
    hrly_rain_maps     = rain_maps.copy()
    hrly_rain_maps.values[1:,:,:] = hrly_rain_maps.values[1:,:,:] - hrly_rain_maps.values[0:-1,:,:]
    
    hrly_rain_maps.values = hrly_rain_maps.values / 25.4
    
    #
    # DBZ
    #
    
    dbz_maps           = wrf.getvar(wrfin    =           ncf,
                                    varname  =          'mdbz',
                                    timeidx  = wrf.ALL_TIMES)
    #
    # Helicity
    #
    
    #srh_maps           = wrf.getvar(wrfin    =           ncf,
    #                                varname  =          'srh',
    #                                timeidx  = wrf.ALL_TIMES)
    
    #
    # T2M
    # 
    
    t2m_maps           = wrf.getvar(wrfin    =           ncf,
                                    varname  =          'T2',
                                    timeidx  = wrf.ALL_TIMES)
    
    t2m_maps.values         = (t2m_maps.values - 273.15) * 9./5. + 32.
    t2m_maps.attrs['units'] = 'degF'



    #
    # U10 / V10 / M10
    # 
        
    u10_maps           = wrf.getvar(wrfin    =           ncf,
                                    varname  =         'U10',
                                    timeidx  = wrf.ALL_TIMES) 
    
    v10_maps           = wrf.getvar(wrfin    =           ncf,
                                    varname  =         'V10',
                                    timeidx  = wrf.ALL_TIMES) 
    
    u10_maps.values         = u10_maps.values * 1.9438444924406     
    v10_maps.values         = v10_maps.values * 1.9438444924406     
    t2m_maps.attrs['units'] = 'kt'



    m10_maps, d10_maps = wrf.getvar(wrfin    =           ncf,
                                    varname  = 'wspd_wdir10',
                                    timeidx  = wrf.ALL_TIMES, 
                                    units    =          'kt') 

    #
    # MSLP
    # 
    
    mslp_maps          = wrf.getvar(wrfin    =           ncf,
                                    varname  =         'slp',
                                    timeidx  = wrf.ALL_TIMES, 
                                    units    =         'hPa')

    
    #
    # Boundary Layer Height
    # 

    pblh_maps          = wrf.getvar(wrfin    =           ncf,
                                    varname  =        'PBLH',
                                    timeidx  = wrf.ALL_TIMES)    
    #
    # Accumulated Snow Total
    # 

    snow_fall_maps     = wrf.getvar(wrfin    =           ncf,
                                    varname  =      'SNOWNC',
                                    timeidx  = wrf.ALL_TIMES) 

    hrly_snowfall_maps                = snow_fall_maps.copy()
    hrly_snowfall_maps.values[1:,:,:] = snow_fall_maps.values[1:,:,:] - snow_fall_maps.values[0:-1,:,:]


    snow_depth_map     = wrf.getvar(wrfin    =           ncf,
                                    varname  =        'SNOWH',
                                    timeidx  = wrf.ALL_TIMES) 
    
    snow_depth_map.values         = snow_depth_map.values * 0.1     
    snow_depth_map.attrs['units'] = 'cm'

   
    hrly_snowdepth_maps                = snow_depth_map.copy()
    hrly_snowdepth_maps.values[1:,:,:] = snow_depth_map.values[1:,:,:] - snow_depth_map.values[0:-1,:,:]
    
    #
    # Accumulated Snow Total
    # 

    lat2d, lon2d = wrf.latlon_coords(var = mslp_maps)
    
    # Extract Eastings and Northings at the grid's corners

    west_east_range   = wrf.cartopy_xlim(var       = mslp_maps, 
                                         timeidx   = 0)

    south_north_range = wrf.cartopy_ylim(var       = mslp_maps,
                                         timeidx   = 0)
    
    ####################################################
    ####################################################
    #
    # Rotate through stations Files
    #
    
    valid_time0         = pd.to_datetime(wrf_time_steps[0].values).tz_localize(tz="UTC").strftime("%Y-%m-%d %H")
    local_time0         = pd.to_datetime(wrf_time_steps[0].values).tz_localize(tz="UTC").tz_convert(tz=tz).strftime("%Y-%m-%d %H")

    valid_timef         = pd.to_datetime(wrf_time_steps[-1].values).tz_localize(tz="UTC").strftime("%Y-%m-%d %H %Z")
    local_timef         = pd.to_datetime(wrf_time_steps[-1].values).tz_localize(tz="UTC").tz_convert(tz=tz).strftime("%Y-%m-%d %H %Z")


        
    for t in range(nt) :

        ####################################################
        #
        # Establish Times and Labels
        #

        valid_time = pd.to_datetime(wrf_time_steps[t].values).tz_localize(tz="UTC").strftime("%Y-%m-%d %H %Z")
        local_time = pd.to_datetime(wrf_time_steps[t].values).tz_localize(tz="UTC").tz_convert(tz=tz).strftime("%Y-%m-%d %H %Z")



        model_run_label    = "Model Run " + wrf_skewt_time + "; WRF Domain " + str(domain).zfill(2)

        print(valid_time, "     -- " + model_run_label)

        #
        ####################################################
        
        
        ####################################################
        #
        # Temperature
        #

        v_name       = "SFCT"
        fig_dir_name = graphics_directory + "/" + v_name + "/"
        file_name    = "wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH + "_F" + str(t).zfill(2) + "_MAP_" + v_name + ".png"

        print(fig_dir_name + file_name)

        fig = plt.figure(figsize=figure_domain_size)

        fig.suptitle(model_run_label)

        ax1 = fig.add_subplot(1,  # nrows
                              1,  # ncols 
                              1,  # index of figure you're installing
                              projection = cart_proj) # cartopy CRS Projection

        ax1.set_xlim(  west_east_range)
        ax1.set_ylim(south_north_range)
        
        
        ax1.add_feature(cfeature.LAKES, 
                        linewidths =     0.5, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.RIVERS, 
                        linewidths =     0.5, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.BORDERS, 
                        linewidths =     1.0, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.COASTLINE, 
                        linewidths =     1.0, 
                        edgecolor  = 'black',
                        facecolor  =  'none')

        ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                     'admin_1_states_provinces_lines',
                                                     '50m',
                                                     linewidths = 0.75,
                                                     facecolor  = 'none',
                                                     edgecolor  = 'black'))
 
        if (domain > 1) :
            ax1.add_feature(USCOUNTIES, linewidths=0.5,edgecolor  = 'black',facecolor='none')
  

        ax1.set_title(valid_time + "  (" + local_time+")")


        line_contour_levels = np.arange(start = 900, 
                                        stop  = 1200,
                                        step  = 4)

        line_contour_map = ax1.contour(lon2d, 
                                      lat2d, 
                                      mslp_maps.isel(Time=t),
                                      transform  = ccrs.PlateCarree(),
                                      colors     = "white",
                                      linewidths = 1,
                                      levels     = line_contour_levels,
                                      extend     = "both")

        ax1.clabel(line_contour_map, inline=True, fontsize=8)

        filled_cm     = ax1.contourf(lon2d, 
                                     lat2d, 
                                     t2m_maps.isel(Time=t),
                                     transform = ccrs.PlateCarree(),
                                     cmap      = temperature_colormap,
                                     levels    = temperature_levels_degF)
        plt.colorbar(filled_cm, 
                     label  = r"2-m Temperature (°F)",
                     shrink = 0.8,
                     pad    = 0.012,
                     format = '%+d')
        gap = 10
        plt.barbs(wrf.to_np(lon2d[::gap,::gap]),
                  wrf.to_np(lat2d[::gap,::gap]),
                  wrf.to_np(u10_maps.isel(Time=t)[::gap,::gap]),
                  wrf.to_np(v10_maps.isel(Time=t)[::gap,::gap]),
                  transform = ccrs.PlateCarree(), 
                  length=5)


        # plt.show()
        plt.tight_layout()
        plt.subplots_adjust(top=0.90)


        fig.savefig(fig_dir_name + file_name)

        plt.close('all')



        #
        ####################################################

        ####################################################
        #
        # Winds
        #

        v_name       = "WIND"
        fig_dir_name = graphics_directory + "/" + v_name + "/"
        file_name    = "wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH + "_F" + str(t).zfill(2) + "_MAP_" + v_name + ".png"

        print(fig_dir_name + file_name)

        fig = plt.figure(figsize=(7,6))

        fig.suptitle(model_run_label)

        ax1 = fig.add_subplot(1,  # nrows
                              1,  # ncols 
                              1,  # index of figure you're installing
                              projection = cart_proj) # cartopy CRS Projection

        ax1.set_xlim(  west_east_range)
        ax1.set_ylim(south_north_range)

        
        ax1.add_feature(cfeature.LAKES, 
                        linewidths =     0.5, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.RIVERS, 
                        linewidths =     0.5, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.BORDERS, 
                        linewidths =     1.0, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.COASTLINE, 
                        linewidths =     1.0, 
                        edgecolor  = 'black',
                        facecolor  =  'none')

        ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                     'admin_1_states_provinces_lines',
                                                     '50m',
                                                     linewidths = 0.75,
                                                     facecolor  = 'none',
                                                     edgecolor  = 'black'))
 
        if (domain > 1) :
            ax1.add_feature(USCOUNTIES, 
                            linewidths=0.5,
                            edgecolor  = 'black',
                            facecolor='none')
            
        ax1.set_title(valid_time + "  (" + local_time+")")

        ax1.streamplot(wrf.to_np(lon2d),
                       wrf.to_np(lat2d), 
                       wrf.to_np(u10_maps.isel(Time=t)),
                       wrf.to_np(v10_maps.isel(Time=t)),
                       transform = ccrs.PlateCarree(),
                       color = "black")

        filled_cm     = ax1.contourf(lon2d, 
                                     lat2d, 
                                     m10_maps.isel(Time=t),
                                     transform = ccrs.PlateCarree(),
                                     cmap      = mpl.cm.Blues,
                                     levels    = np.arange(0,41,1))
        plt.colorbar(filled_cm, 
                     label  = r"10-m Wind Speed (kts)",
                     shrink = 0.8,
                     pad    = 0.012)



        # plt.show()
        plt.tight_layout()
        plt.subplots_adjust(top=0.90)


        fig.savefig(fig_dir_name + file_name)

        plt.close('all')


        #
        ####################################################

        ####################################################
        #
        # Hourly Rainfall
        #

        v_name       = "RAIN"
        fig_dir_name = graphics_directory + "/" + v_name + "/"
        file_name    = "wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH + "_F" + str(t).zfill(2) + "_MAP_" + v_name + ".png"

        print(fig_dir_name + file_name)

        fig = plt.figure(figsize=(7,6))

        fig.suptitle(model_run_label)

        ax1 = fig.add_subplot(1,  # nrows
                              1,  # ncols 
                              1,  # index of figure you're installing
                              projection = cart_proj) # cartopy CRS Projection

        ax1.set_xlim(  west_east_range)
        ax1.set_ylim(south_north_range)

        
        ax1.add_feature(cfeature.LAKES, 
                        linewidths =     0.5, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.RIVERS, 
                        linewidths =     0.5, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.BORDERS, 
                        linewidths =     1.0, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.COASTLINE, 
                        linewidths =     1.0, 
                        edgecolor  = 'black',
                        facecolor  =  'none')

        ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                     'admin_1_states_provinces_lines',
                                                     '50m',
                                                     linewidths = 0.75,
                                                     facecolor  = 'none',
                                                     edgecolor  = 'black'))
 
        if (domain > 1) :
            ax1.add_feature(USCOUNTIES, 
                            linewidths=0.5,
                            edgecolor  = 'black',
                            facecolor='none')
            
            
         
        
        ax1.set_title(valid_time + "  (" + local_time+")")

        rain_norm = mpl.colors.BoundaryNorm(boundaries = precip_levels_in, 
                                                ncolors    = 15)

        filled_cm     = ax1.contourf(lon2d, 
                                     lat2d, 
                                     hrly_rain_maps.isel(Time=t),
                                     transform = ccrs.PlateCarree(),
                                     norm      = rain_norm,
                                     cmap      = precip_colormap,
                                     extend    = 'max',
                                     levels    = precip_levels_in)
        plt.colorbar(filled_cm, 
                     label  = "Hourly Precip (in)",
                     shrink = 0.8,
                     ticks=precip_levels_in,
                     pad    = 0.012)



        # plt.show()
        plt.tight_layout()
        plt.subplots_adjust(top=0.90)


        fig.savefig(fig_dir_name + file_name)

        plt.close('all')


        #
        ####################################################
        

        ####################################################
        #
        # Hourly  Snowfall Equivalent
        #

        v_name       = "WEASD"
        fig_dir_name = graphics_directory + "/" + v_name + "/"
        file_name    = "wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH + "_F" + str(t).zfill(2) + "_MAP_" + v_name + ".png"

        print(fig_dir_name + file_name)


        fig = plt.figure(figsize=(7,6))

        fig.suptitle(model_run_label)

        ax1 = fig.add_subplot(1,  # nrows
                              1,  # ncols 
                              1,  # index of figure you're installing
                              projection = cart_proj) # cartopy CRS Projection

        ax1.set_xlim(  west_east_range)
        ax1.set_ylim(south_north_range)

        
        ax1.add_feature(cfeature.LAKES, 
                        linewidths =     0.5, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.RIVERS, 
                        linewidths =     0.5, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.BORDERS, 
                        linewidths =     1.0, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.COASTLINE, 
                        linewidths =     1.0, 
                        edgecolor  = 'black',
                        facecolor  =  'none')

        ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                     'admin_1_states_provinces_lines',
                                                     '50m',
                                                     linewidths = 0.75,
                                                     facecolor  = 'none',
                                                     edgecolor  = 'black'))
 
        if (domain > 1) :
            ax1.add_feature(USCOUNTIES, 
                            linewidths=0.5,
                            edgecolor  = 'black',
                            facecolor='none')
            
        ax1.set_title(valid_time + "  (" + local_time+")")

        rain_norm = mpl.colors.BoundaryNorm(boundaries = precip_levels_in, 
                                                ncolors    = 15)

        filled_cm     = ax1.contourf(lon2d, 
                                     lat2d, 
                                     hrly_snowfall_maps.isel(Time=t),
                                     transform = ccrs.PlateCarree(),
                                     norm      = rain_norm,
                                     cmap      = precip_colormap,
                                     extend    = 'max',
                                     levels    = precip_levels_in)
        plt.colorbar(filled_cm, 
                     label  = "Hourly Snow-Water Equivalent (in)",
                     shrink = 0.8,
                     ticks=precip_levels_in,
                     pad    = 0.012)



        # plt.show()
        plt.tight_layout()
        plt.subplots_adjust(top=0.90)


        fig.savefig(fig_dir_name + file_name)

        plt.close('all')


        #
        ####################################################

        ####################################################
        #
        # Hourly Accumulated Depth
        #

        v_name       = "SNOWH"
        fig_dir_name = graphics_directory + "/" + v_name + "/"
        file_name    = "wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH + "_F" + str(t).zfill(2) + "_MAP_" + v_name + ".png"

        print(fig_dir_name + file_name)


        fig = plt.figure(figsize=(7,6))

        fig.suptitle(model_run_label)

        ax1 = fig.add_subplot(1,  # nrows
                              1,  # ncols 
                              1,  # index of figure you're installing
                              projection = cart_proj) # cartopy CRS Projection

        ax1.set_xlim(  west_east_range)
        ax1.set_ylim(south_north_range)

        
        ax1.add_feature(cfeature.LAKES, 
                        linewidths =     0.5, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.RIVERS, 
                        linewidths =     0.5, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.BORDERS, 
                        linewidths =     1.0, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.COASTLINE, 
                        linewidths =     1.0, 
                        edgecolor  = 'black',
                        facecolor  =  'none')

        ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                     'admin_1_states_provinces_lines',
                                                     '50m',
                                                     linewidths = 0.75,
                                                     facecolor  = 'none',
                                                     edgecolor  = 'black'))
 
        if (domain > 1) :
            ax1.add_feature(USCOUNTIES, 
                            linewidths=0.5,
                            edgecolor  = 'black',
                            facecolor='none')
            
        ax1.set_title(valid_time + "  (" + local_time+")")



        rain_norm = mpl.colors.BoundaryNorm(boundaries = precip_levels_in, 
                                                ncolors    = 15)

        filled_cm     = ax1.contourf(lon2d, 
                                     lat2d, 
                                     snow_depth_map.isel(Time=t),
                                     transform = ccrs.PlateCarree(),
                                     norm      = rain_norm,
                                     cmap      = precip_colormap,
                                     extend    = 'max',
                                     levels    = precip_levels_in)
        plt.colorbar(filled_cm, 
                     label  = "Hourly Snow Depth (cm)",
                     shrink = 0.8,
                     ticks=precip_levels_in,
                     pad    = 0.012)



        # plt.show()
        plt.tight_layout()
        plt.subplots_adjust(top=0.90)


        fig.savefig(fig_dir_name + file_name)

        plt.close('all')


        #
        ####################################################

        ####################################################
        #
        # Hourly DBZ
        #

        v_name       = "DBZ"
        fig_dir_name = graphics_directory + "/" + v_name + "/"
        file_name    = "wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH + "_F" + str(t).zfill(2) + "_MAP_" + v_name + ".png"

        fig = plt.figure(figsize=(7,6))

        fig.suptitle(model_run_label)

        ax1 = fig.add_subplot(1,  # nrows
                              1,  # ncols 
                              1,  # index of figure you're installing
                              projection = cart_proj) # cartopy CRS Projection

        ax1.set_xlim(  west_east_range)
        ax1.set_ylim(south_north_range)

        
        ax1.add_feature(cfeature.LAKES, 
                        linewidths =     0.5, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.RIVERS, 
                        linewidths =     0.5, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.BORDERS, 
                        linewidths =     1.0, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.COASTLINE, 
                        linewidths =     1.0, 
                        edgecolor  = 'black',
                        facecolor  =  'none')

        ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                     'admin_1_states_provinces_lines',
                                                     '50m',
                                                     linewidths = 0.75,
                                                     facecolor  = 'none',
                                                     edgecolor  = 'black'))
 
        if (domain > 1) :
            ax1.add_feature(USCOUNTIES, 
                            linewidths=0.5,
                            edgecolor  = 'black',
                            facecolor='none')
            
        ax1.set_title(valid_time + "  (" + local_time+")")



        rain_norm = mpl.colors.BoundaryNorm(boundaries = precip_levels_in, 
                                                ncolors    = 15)

        filled_cm     = ax1.contourf(lon2d, 
                                     lat2d, 
                                     dbz_maps.isel(Time=t),
                                     transform = ccrs.PlateCarree(),
                                     levels    = stormy_dbz_values,
                                     cmap      =  nws_dbz_colors,
                                     extend    = 'max')
        plt.colorbar(filled_cm, 
                     label  = "Max Columnar Storm Reflectivity (dbZ)",
                     shrink = 0.8,
                     ticks = stormy_dbz_values,
                     pad    = 0.012)




        # plt.show()
        plt.tight_layout()
        plt.subplots_adjust(top=0.90)


        fig.savefig(fig_dir_name + file_name)

        plt.close('all')

        #
        ####################################################


        ####################################################
        #
        # PBLH
        #

        v_name       = "PBL"
        fig_dir_name = graphics_directory + "/" + v_name + "/"
        file_name    = "wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH + "_F" + str(t).zfill(2) + "_MAP_" + v_name + ".png"

        print(fig_dir_name + file_name)

        pbl_height_levels = np.arange(0,5000,100)

        fig = plt.figure(figsize=(7,6))

        fig.suptitle(model_run_label)

        ax1 = fig.add_subplot(1,  # nrows
                              1,  # ncols 
                              1,  # index of figure you're installing
                              projection = cart_proj) # cartopy CRS Projection

        ax1.set_xlim(  west_east_range)
        ax1.set_ylim(south_north_range)

        
        ax1.add_feature(cfeature.LAKES, 
                        linewidths =     0.5, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.RIVERS, 
                        linewidths =     0.5, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.BORDERS, 
                        linewidths =     1.0, 
                        edgecolor  = 'black',
                        facecolor  =  'none')
        ax1.add_feature(cfeature.COASTLINE, 
                        linewidths =     1.0, 
                        edgecolor  = 'black',
                        facecolor  =  'none')

        ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                     'admin_1_states_provinces_lines',
                                                     '50m',
                                                     linewidths = 0.75,
                                                     facecolor  = 'none',
                                                     edgecolor  = 'black'))
 
        if (domain > 1) :
            ax1.add_feature(USCOUNTIES, 
                            linewidths=0.5,
                            edgecolor  = 'black',
                            facecolor='none')
            
        ax1.set_title(valid_time + "  (" + local_time+")")




        filled_cm     = ax1.contourf(lon2d, 
                                     lat2d, 
                                     pblh_maps.isel(Time=t),
                                     transform = ccrs.PlateCarree(),
                                     levels    = pbl_height_levels,
                                     cmap      = mpl.cm.rainbow,
                                     extend    = 'max')
        plt.colorbar(filled_cm, 
                     label  = "PBL Height (m)",
                     shrink = 0.8,
                     #ticks = stormy_dbz_values,
                     pad    = 0.012)



 
        # plt.show()
        plt.tight_layout()
        plt.subplots_adjust(top=0.90)


        fig.savefig(fig_dir_name + file_name)

        plt.close('all')

        #
        ####################################################

    ####################################################
    #
    # making gifs
    #
    
    for v_name in ("DBZ", "PBL", "RAIN", "SFCT", "SNOWH", "WIND", "WEASD"):          
        fig_dir_name = graphics_directory + "/" + v_name + "/"
        png_file_name    = "wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH + "_F??_MAP_" + v_name + ".png"
        gif_file_name    = "wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH + "_Fxx_MAP_" + v_name + ".gif"
        os.system("convert -delay 25 " + fig_dir_name + png_file_name + " " + fig_dir_name + gif_file_name)
    
    #
    ####################################################    
    
        
    ####################################################
    #
    # 36-hr Total Rainfall
    #

    v_name       = "TRAIN"
    fig_dir_name = graphics_directory + "/" + v_name + "/"
    file_name    = "wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH + "_FXX_MAP_" + v_name + ".png"

    print(fig_dir_name + file_name)


    fig = plt.figure(figsize=(7,6))

    fig.suptitle(model_run_label)

    ax1 = fig.add_subplot(1,  # nrows
                          1,  # ncols 
                          1,  # index of figure you're installing
                          projection = cart_proj) # cartopy CRS Projection

    ax1.set_xlim(  west_east_range)
    ax1.set_ylim(south_north_range)

    
    ax1.add_feature(cfeature.LAKES, 
                    linewidths =     0.5, 
                    edgecolor  = 'black',
                    facecolor  =  'none')
    ax1.add_feature(cfeature.RIVERS, 
                    linewidths =     0.5, 
                    edgecolor  = 'black',
                    facecolor  =  'none')
    ax1.add_feature(cfeature.BORDERS, 
                    linewidths =     1.0, 
                    edgecolor  = 'black',
                    facecolor  =  'none')
    ax1.add_feature(cfeature.COASTLINE, 
                    linewidths =     1.0, 
                    edgecolor  = 'black',
                    facecolor  =  'none')

    ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                 'admin_1_states_provinces_lines',
                                                 '50m',
                                                 linewidths = 0.75,
                                                 facecolor  = 'none',
                                                 edgecolor  = 'black'))

    if (domain > 1) :
        ax1.add_feature(USCOUNTIES, 
                        linewidths=0.5,
                        edgecolor  = 'black',
                        facecolor='none')
            



    ax1.set_title(valid_time0 + " - " + valid_timef + "  (" + local_time0 + " - " + local_timef+")")

    rain_norm = mpl.colors.BoundaryNorm(boundaries = precip_levels_in, 
                                            ncolors    = 15)

    filled_cm     = ax1.contourf(lon2d, 
                                 lat2d, 
                                 rain_maps.isel(Time=36),
                                 transform = ccrs.PlateCarree(),
                                 norm      = rain_norm,
                                 cmap      = precip_colormap,
                                 extend    = 'max',
                                 levels    = precip_levels_in)
    plt.colorbar(filled_cm, 
                 label  = "36-hr Total Precip (in)",
                 shrink = 0.8,
                 ticks = precip_levels_in,
                 pad    = 0.012)




    # plt.show()
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)


    fig.savefig(fig_dir_name + file_name)

    plt.close('all')

    #
    ####################################################


    ####################################################
    #
    # Total Snowfall
    #

    v_name       = "TSNOW"
    fig_dir_name = graphics_directory + "/" + v_name + "/"
    file_name    = "wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH + "_FXX_MAP_" + v_name + ".png"

    print(fig_dir_name + file_name)

    fig = plt.figure(figsize=(7,6))

    fig.suptitle(model_run_label)

    ax1 = fig.add_subplot(1,  # nrows
                          1,  # ncols 
                          1,  # index of figure you're installing
                          projection = cart_proj) # cartopy CRS Projection

    ax1.set_xlim(  west_east_range)
    ax1.set_ylim(south_north_range)

    
    ax1.add_feature(cfeature.LAKES, 
                    linewidths =     0.5, 
                    edgecolor  = 'black',
                    facecolor  =  'none')
    ax1.add_feature(cfeature.RIVERS, 
                    linewidths =     0.5, 
                    edgecolor  = 'black',
                    facecolor  =  'none')
    ax1.add_feature(cfeature.BORDERS, 
                    linewidths =     1.0, 
                    edgecolor  = 'black',
                    facecolor  =  'none')
    ax1.add_feature(cfeature.COASTLINE, 
                    linewidths =     1.0, 
                    edgecolor  = 'black',
                    facecolor  =  'none')

    ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                 'admin_1_states_provinces_lines',
                                                 '50m',
                                                 linewidths = 0.75,
                                                 facecolor  = 'none',
                                                 edgecolor  = 'black'))

    if (domain > 1) :
        ax1.add_feature(USCOUNTIES, 
                        linewidths=0.5,
                        edgecolor  = 'black',
                        facecolor  =  'none')   
                        
    ax1.set_title(valid_time0 + " - " + valid_timef + "  (" + local_time0 + " - " + local_timef+")")

    rain_norm = mpl.colors.BoundaryNorm(boundaries = precip_levels_in, 
                                            ncolors    = 15)

    filled_cm     = ax1.contourf(lon2d, 
                                 lat2d, 
                                 snow_fall_maps.isel(Time=36),
                                 transform = ccrs.PlateCarree(),
                                 norm      = rain_norm,
                                 cmap      = precip_colormap,
                                 extend    = 'max',
                                 levels    = precip_levels_in)
    plt.colorbar(filled_cm, 
                 label  = "36-hr Total Snowfall (in)",
                 shrink = 0.8,
                 ticks = precip_levels_in,
                 pad    = 0.012)




    # plt.show()
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)


    fig.savefig(fig_dir_name + file_name)

    plt.close('all')

    #
    ####################################################

    #
    ####################################################
    ####################################################
            

print("done")
#
####################################################
####################################################
####################################################


# ## Ending Script

# In[ ]:


####################################################
####################################################
####################################################
#
# End of Script
#

print("End Sounding Plotting Script")

#
####################################################
####################################################
####################################################

