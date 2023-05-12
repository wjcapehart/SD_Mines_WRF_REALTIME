#!/usr/bin/env python
# coding: utf-8

# # Plotting Skew-Ts

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


import matplotlib              as mpl
import matplotlib.pyplot       as plt
import matplotlib.patches      as patches
import matplotlib.font_manager as fm
import matplotlib              as mpl
import matplotlib.gridspec     as gridspec
from   matplotlib.ticker       import (MultipleLocator, NullFormatter, ScalarFormatter)

import seaborn           as sns

import pandas            as pd
import xarray            as xr
import pint_xarray       as px

import scipy.interpolate as scintp 

import netCDF4           as nc4

import wrf               as wrf



import metpy.calc  as     mpcalc
from   metpy.plots import SkewT, Hodograph
from   metpy.units import units, pandas_dataframe_to_unit_arrays


import timezonefinder    as tzf
import pytz              as pytz
import socket            as socket

#
####################################################
####################################################
####################################################


# In[ ]:





# In[ ]:


####################################################
####################################################
####################################################
#
# Mines Colors and Fonts
#

Mines_Blue = "#002554"


plt.rcParams.update({'text.color'      : Mines_Blue,
                     'axes.labelcolor' : Mines_Blue,
					 'axes.edgecolor'  : Mines_Blue,
					 'xtick.color'     : Mines_Blue,
					 'ytick.color'     : Mines_Blue})


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
# Mines Colors and Fonts
#

Mines_Blue = "#002554"


plt.rcParams.update({'text.color'      : Mines_Blue,
                     'axes.labelcolor' : Mines_Blue,
					 'axes.edgecolor'  : Mines_Blue,
					 'xtick.color'     : Mines_Blue,
					 'ytick.color'     : Mines_Blue})


#
####################################################
####################################################
########################################################################################################
####################################################
####################################################
#
# File Organization
#
intel         = True
beta_on       = 0
max_domains   = 2
chosen_domain = 2

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

TS_DIR = WRF_EXE

station_list_file = WRF_OVERALL_DIR + "namelist_files_and_local_scripts/time_series_station_files_"+str(max_domains)+"_dom_all.xlsx"

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


print(model_start_date_YYYY_MM_DD_HH0000)

#
####################################################
####################################################
####################################################


# ## Read tslist excel file

# In[ ]:


####################################################
####################################################
####################################################
#
# Read TSLIST Excel File
#

print("read file from "+station_list_file)

available_time_series_list = pd.read_excel(station_list_file,
                                           index_col=0)

print(available_time_series_list)

#
####################################################
####################################################
####################################################


# ## Crack WRF Files

# In[ ]:





# In[ ]:





# In[ ]:


####################################################
####################################################
####################################################
#
# Rotate through Available Files
#

for domain in range(chosen_domain,chosen_domain+1):
    station_doms = available_time_series_list[available_time_series_list['Domain'] == domain]
    print(station_doms)
    
    wrf_file  = WRF_EXE  + "./wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH0000
    
    ncf = nc4.Dataset(filename = wrf_file)
    
    wrf_time_steps     = wrf.getvar(wrfin    =           ncf,
                                    varname  =       'times',
                                    timeidx  = wrf.ALL_TIMES)  
    
    nt                 = len(wrf_time_steps)
    
    temperature_4d     = wrf.getvar(wrfin    =           ncf,
                                    varname  =          'tc',
                                    timeidx  = wrf.ALL_TIMES)   * units.degC

    dew_point_4d       = wrf.getvar(wrfin    =           ncf,
                                    varname  =          'td',
                                    timeidx  = wrf.ALL_TIMES, 
                                    units    =        'degC')  * units.degC

    isobar_hgt_4d      = wrf.getvar(wrfin    =           ncf,
                                    varname  =        'pres',
                                    timeidx  = wrf.ALL_TIMES, 
                                    units    =         'hPa') * units.hPa

    uwind_4d, vwind_4d = wrf.getvar(wrfin    =           ncf,
                                    varname  =       'uvmet',
                                    timeidx  = wrf.ALL_TIMES, 
                                    units    =          'kt') * units.kt
    
    height_4d          = wrf.getvar(wrfin    =           ncf,
                                    varname  =      'height',
                                    timeidx  = wrf.ALL_TIMES, 
                                    units    =           'm') * units.m
    
    height_agl_4d      = wrf.getvar(wrfin        =           ncf,
                                    varname      =  'height_agl',
                                    timeidx      = wrf.ALL_TIMES, 
                                    units        =           'm') * units.m

    qv_4d              = wrf.getvar(wrfin        =           ncf,
                                    varname      =      'QVAPOR',
                                    timeidx      = wrf.ALL_TIMES)     
    qc_4d              = wrf.getvar(wrfin        =           ncf,
                                    varname      =      'QCLOUD',
                                    timeidx      = wrf.ALL_TIMES)    
    
    qi_4d              = wrf.getvar(wrfin        =           ncf,
                                    varname      =         'QICE',
                                    timeidx      = wrf.ALL_TIMES)    
    
    qr_4d              = wrf.getvar(wrfin        =           ncf,
                                    varname      =      'QRAIN',
                                    timeidx      = wrf.ALL_TIMES)    
    
    qs_4d              = wrf.getvar(wrfin        =           ncf,
                                    varname      =      'QSNOW',
                                    timeidx      = wrf.ALL_TIMES)    
    
    qg_4d              = wrf.getvar(wrfin        =           ncf,
                                    varname      =      'QGRAUP',
                                    timeidx      = wrf.ALL_TIMES)     
    
    z_700_500_4d = wrf.vinterp(wrfin         =           ncf,
                               field         = height_agl_4d,
                               vert_coord    =    "pressure",
                               interp_levels =     [700,500],
                               log_p         =          True,
                               timeidx       = wrf.ALL_TIMES) * units.m

    t_700_500_4d = wrf.vinterp(wrfin        =            ncf,
                              field         = temperature_4d,
                              vert_coord    =     "pressure",
                              interp_levels =      [700,500],
                              log_p         =           True,
                              timeidx       =  wrf.ALL_TIMES) * units.degC
        
    lat2d, lon2d = wrf.latlon_coords(var = height_4d)
    
    ####################################################
    ####################################################
    #
    # Rotate through stations Files
    #

    for station in station_doms.iterrows():
        station_id     = station[1][0]
        grid_domain    = station[1][1]
        station_name   = station[1][2]
        station_lat    = station[1][3]
        station_lon    = station[1][4]
        
        tf     = tzf.TimezoneFinder()
        tz     = tf.certain_timezone_at(lng=station_lon, lat=station_lat)
        tzabbr = pytz.timezone(tz).localize(model_start_datetime)
        
        #
        # Creating Graphics Directory
        #

        graphics_directory = WRF_IMAGES + "/" + model_start_date_YYYY_MM_DD_HH + "/SKEWTS/" + station_id + "/"

        print("Creating " + graphics_directory)

        os.system("mkdir -pv " + graphics_directory )
        
        wrf_x, wrf_y = wrf.ll_to_xy(wrfin     =         ncf,
                                    latitude  = station_lat, 
                                    longitude = station_lon, 
                                    timeidx   = 0, 
                                    squeeze   = True, 
                                    meta      = False, 
                                    stagger   = None, 
                                    as_int    = True)
        
        clouds_maxx = np.max([qv_4d[:, :, wrf_y, wrf_x],
                              qc_4d[:, :, wrf_y, wrf_x],
                              qi_4d[:, :, wrf_y, wrf_x],
                              qr_4d[:, :, wrf_y, wrf_x],
                              qs_4d[:, :, wrf_y, wrf_x],
                              qg_4d[:, :, wrf_y, wrf_x]]) * 1000.
        print('  Max q* = ',clouds_maxx, "g/kg")
               
        z_700_500 = z_700_500_4d[:,:,wrf_y,wrf_x]
        t_700_500 = t_700_500_4d[:,:,wrf_y,wrf_x]
        
        above_30        = np.arange(30,51,5)
        max_temperature = temperature_4d[:, :,wrf_y, wrf_x].values.max()
        max_axis_temp   = above_30[above_30 > max_temperature].min() 
        
        sounding_file_name_gif = "wrfout_dxx_" + model_start_date_YYYY_MM_DD_HH + "_SKEWT_" + station_id + ".gif"

        print(" - " + sounding_file_name_gif)
               
        ####################################################
        ### ### ### ### ### ### ### ### ### ### ### ### ####
        #
        # Rotate through Time Steps
        #
        
        for t in range(nt) :

            ####################################################
            #
            # Extract 3D fields for each timestep
            #

            
            
            #
            ####################################################

            ####################################################
            #
            # Extract Point Locations
            #
        

            valid_time = pd.to_datetime(wrf_time_steps[t].values).tz_localize(tz="UTC").strftime("%Y-%m-%d %H %Z")
            local_time = pd.to_datetime(wrf_time_steps[t].values).tz_localize(tz="UTC").tz_convert(tz=tz).strftime("%Y-%m-%d %H %Z")
            local_time_zone = pd.to_datetime(wrf_time_steps[t].values).tz_localize(tz="UTC").tz_convert(tz=tz).strftime("%Z")
            dow = pd.to_datetime(wrf_time_steps[t].values).tz_localize(tz="UTC").tz_convert(tz=tz).strftime("%a")


            time_for_clock = pd.to_datetime(wrf_time_steps[t].values).tz_localize(tz="UTC").tz_convert(tz=tz).time()

            hour   = time_for_clock.hour
            minute = time_for_clock.minute
            second = time_for_clock.second
            percent_done = (1.0 * t) / (nt-1.)


            if ((hour >= 6) and (hour < 18)):
                Clock_Color = Mines_Blue
                Clock_BgndC = "white"           
            else:
                Clock_Color = "white"
                Clock_BgndC = Mines_Blue        
                
            #
            # Walking Clock Parameters
            #

            circle_theta  = np.deg2rad(np.arange(0,360,0.01))
            circle_radius = circle_theta * 0 + 1

            if (hour > 12) :
                hour = hour - 12

            angles_h = 2*np.pi*hour/12+2*np.pi*minute/(12*60)+2*second/(12*60*60)
            angles_m = 2*np.pi*minute/60+2*np.pi*second/(60*60)
            
            
            #
            # Labels
            #

            model_run_label    = "Model Run " + wrf_skewt_time + "; WRF Domain " + str(domain).zfill(2)

            print(valid_time, "     -- " + model_run_label, " :: ", (percent_done*100) )


            sounding_file_name_png = "wrfout_dxx_" + model_start_date_YYYY_MM_DD_HH + "_F" + str(t).zfill(2) + "_SKEWT_" + station_id + ".png"
            
            print("     -- " + sounding_file_name_png)
    
            #
            # Make Sounding Dataframe
            #

            


            sounding_df = pd.DataFrame({"pressure"    :   isobar_hgt_4d[t, :, wrf_y, wrf_x],
                                        "temperature" :  temperature_4d[t, :, wrf_y, wrf_x],
                                        "dewpoint"    :    dew_point_4d[t, :, wrf_y, wrf_x],
                                        "u_wind"      :        uwind_4d[t, :, wrf_y, wrf_x],
                                        "v_wind"      :        vwind_4d[t, :, wrf_y, wrf_x],
                                        "height"      :       height_4d[t, :, wrf_y, wrf_x],
                                        "agl"         :   height_agl_4d[t, :, wrf_y, wrf_x],
                                        "qv"          :           qv_4d[t, :, wrf_y, wrf_x]*1000,
                                        "qc"          :           qc_4d[t, :, wrf_y, wrf_x]*1000,
                                        "qi"          :           qi_4d[t, :, wrf_y, wrf_x]*1000,
                                        "qr"          :           qr_4d[t, :, wrf_y, wrf_x]*1000,
                                        "qs"          :           qs_4d[t, :, wrf_y, wrf_x]*1000,
                                        "qg"          :           qg_4d[t, :, wrf_y, wrf_x]*1000})

            units_df = {"pressure"    :  "hPa",
                        "temperature" : "degC",
                        "dewpoint"    : "degC",
                        "u_wind"      :   "kt",
                        "v_wind"      :   "kt",
                        "height"      :    "m",
                        "agl"         :    "m",
                        "qv"          : "g/kg",
                        "qc"          : "g/kg",
                        "qi"          : "g/kg",
                        "qr"          : "g/kg",
                        "qs"          : "g/kg",
                        "qh"          : "g/kg",
                        "speed"       : "kt"}


            sounding_df = pandas_dataframe_to_unit_arrays(df           = sounding_df, 
                                                          column_units =    units_df)
            
            
            #################################################
            #
            # Make Hodo Data Frame
            #

            hodo_frame = pd.DataFrame.from_dict(sounding_df)[["agl","u_wind","v_wind"]].copy()
            hodo_frame["speed"] = np.sqrt((hodo_frame["u_wind"].values**2) + (hodo_frame["v_wind"].values**2))



            hodo_extra = pd.DataFrame({"agl":[500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000]})



            fx = scintp.interp1d(x    = hodo_frame[   "agl"].values, 
                                 y    = hodo_frame["u_wind"].values,
                                 kind = "cubic")


            fy = scintp.interp1d(x    = hodo_frame[   "agl"].values, 
                                 y    = hodo_frame["v_wind"].values,
                                 kind = "cubic")

            hodo_extra["u_wind"] = fx(hodo_extra["agl"].values)
            hodo_extra["v_wind"] = fy(hodo_extra["agl"].values)
            hodo_extra["speed"]  = np.sqrt((hodo_extra["u_wind"].values**2 + hodo_extra["v_wind"].values**2))

            hodo_extra[hodo_extra["speed"] <= 40]



            hodo_frame = pd.concat([hodo_frame, hodo_extra])
            hodo_frame = hodo_frame.sort_values(by = "agl",
                                               ignore_index = True)
            


            hodo_extra = pandas_dataframe_to_unit_arrays(df           = hodo_extra, 
                                                         column_units =   units_df)
            hodo_frame = pandas_dataframe_to_unit_arrays(df           = hodo_frame, 
                                                         column_units =   units_df)
            

            hodo_extra["u_wind"] = hodo_extra["u_wind"].to(units("m/s"))
            hodo_extra["v_wind"] = hodo_extra["v_wind"].to(units("m/s"))
            hodo_extra["v_wind"] = hodo_extra["v_wind"].to(units("m/s"))
            hodo_frame["u_wind"] = hodo_frame["u_wind"].to(units("m/s"))
            hodo_extra["speed"]  = hodo_extra["speed"].to(units("m/s"))
            hodo_frame["speed"]  = hodo_frame["speed"].to(units("m/s"))
            

                
            #
            #################################################

            # Calculate thermodynamics
            lcl_pressure, lcl_temperature = mpcalc.lcl(sounding_df["pressure"][0],
                                                       sounding_df["temperature"][0],
                                                       sounding_df["dewpoint"][0])

            lfc_pressure, lfc_temperature = mpcalc.lfc(sounding_df["pressure"],
                                                       sounding_df["temperature"],
                                                       sounding_df["dewpoint"])

            el_pressure, el_temperature = mpcalc.el(sounding_df["pressure"],
                                                    sounding_df["temperature"],
                                                    sounding_df["dewpoint"])

            parcel_profile = mpcalc.parcel_profile(sounding_df["pressure"],
                                                   sounding_df["temperature"][0],
                                                   sounding_df["dewpoint"][0])

            surface_cape, surface_cin = mpcalc.surface_based_cape_cin(sounding_df["pressure"],
                                                                      sounding_df["temperature"],
                                                                      sounding_df["dewpoint"])

            lcl_hgt = np.round(mpcalc.pressure_to_height_std(lcl_pressure), 
                               decimals=3).to(units.meter)
            
            lfc_hgt = np.round(mpcalc.pressure_to_height_std(lfc_pressure), 
                               decimals=3).to(units.meter)

            sb_cape, sb_cin = mpcalc.surface_based_cape_cin(sounding_df["pressure"], 
                                                            sounding_df["temperature"], 
                                                            sounding_df["dewpoint"])
            ml_cape, ml_cin = mpcalc.mixed_layer_cape_cin(sounding_df["pressure"], 
                                                          sounding_df["temperature"], 
                                                          sounding_df["dewpoint"])
            mu_cape, mu_cin = mpcalc.most_unstable_cape_cin(sounding_df["pressure"], 
                                                            sounding_df["temperature"], 
                                                            sounding_df["dewpoint"])

            t_700m500  =  (t_700_500[t,1].values- t_700_500[t,0].values)*units.delta_degC
            z_700m500  = ((z_700_500[t,1].values- z_700_500[t,0].values)/1000)*units.km
            lr_700_500 = np.round(-1 * np.divide( t_700m500, z_700m500),2)

            sbcape = np.round(sb_cape, 1)
            sbcin  = np.round(sb_cin,  1)
            mlcape = np.round(ml_cape, 1)
            mlcin  = np.round(ml_cin,  1)
            mucape = np.round(mu_cape, 1)

            u_shear01, v_shear01 = mpcalc.bulk_shear(sounding_df["pressure"], 
                                                     sounding_df["u_wind"].to(units.meter/units.second), 
                                                     sounding_df["v_wind"].to(units.meter/units.second), 
                                                     depth = 1000 * units.meter)
            
            shear01 = np.round((np.sqrt(u_shear01**2 + v_shear01**2)), 1)
            
            u_shear06, v_shear06 = mpcalc.bulk_shear(sounding_df["pressure"], 
                                                     sounding_df["u_wind"].to(units.meter/units.second), 
                                                     sounding_df["v_wind"].to(units.meter/units.second), 
                                                     depth = 6000 * units.meter)
            
            shear06 = np.round((np.sqrt(u_shear06**2 + v_shear06**2)), 1)
            
            rmover, lmover, mean = mpcalc.bunkers_storm_motion(sounding_df["pressure"], 
                                                               sounding_df["u_wind"].to(units.meter/units.second), 
                                                               sounding_df["v_wind"].to(units.meter/units.second), 
                                                               sounding_df["agl"])
            
            srh_01_pos, srh_01_neg, srh_01_tot = mpcalc.storm_relative_helicity(sounding_df["height"],
                                                                                sounding_df["u_wind"].to(units.meter/units.second), 
                                                                                sounding_df["v_wind"].to(units.meter/units.second),  
                                                                                depth   = 1000 * units.meter, 
                                                                                bottom  = sounding_df["agl"][0], 
                                                                                storm_u = lmover[0], 
                                                                                storm_v = lmover[1])
            srh_01 = np.round(srh_01_neg, 1)
            
            srh_03_pos, srh_03_neg, srh_03_tot = mpcalc.storm_relative_helicity(sounding_df["height"],
                                                                                sounding_df["u_wind"].to(units.meter/units.second), 
                                                                                sounding_df["v_wind"].to(units.meter/units.second), 
                                                                                depth   = 3000 * units.meter, 
                                                                                bottom  = sounding_df["agl"][0], 
                                                                                storm_u = lmover[0], storm_v = lmover[1])
            srh_03 = np.round(srh_03_neg, 1)

            ###################################################
            #
            # Generate SkewT
            #

            fig = plt.figure(figsize   = (9, 9))

            fig.suptitle(station_name + "; Model Run " + wrf_skewt_time + "; WRF Domain " + str(domain).zfill(2), 
                         fontsize=18,
                         verticalalignment = "center",
                         horizontalalignment = "center")

            # Grid for plots

            gs   = gridspec.GridSpec(nrows = 3, 
                                     ncols = 3)

            skew = SkewT(fig      =       fig, 
                         rotation =        45,  
                         rect     = [0.10, 0.00, 0.60, 0.95])

            # Plot the sounding using normal plotting functions, in this case using
            # log scaling in Y, as dictated by the typical meteorological plot

            skew.plot(sounding_df["pressure"], 
                      sounding_df["temperature"], 
                      color = 'red')
            
            skew.plot(sounding_df["pressure"], 
                      sounding_df["dewpoint"], 
                      color = 'green')
            
            skew.plot(sounding_df["pressure"], 
                      parcel_profile, 
                      color = Mines_Blue, 
                      linewidth = 0.5)

            # Mask barbs to be below 100 hPa only

            mask = (sounding_df["pressure"] >= 100 * units.hPa)
            for i in range(mask.size):
                if (sounding_df["pressure"][i] > 500 * units.hPa) :
                    if (i%2 == 1):
                        mask[i] = False

            skew.plot_barbs(sounding_df["pressure"][mask], 
                            sounding_df[  "u_wind"][mask], 
                            sounding_df[  "v_wind"][mask],
                            color=Mines_Blue)

            skew.ax.set_ylim=(1000, 100)

            # Add the relevant special lines

            skew.plot_dry_adiabats()
            skew.plot_moist_adiabats()
            skew.plot_mixing_lines()

            #
            ###################################################


            ###################################################
            #
            # Skew-T Thermodynamics
            #


            # Shade areas

            skew.shade_cin(sounding_df["pressure"], 
                           sounding_df["temperature"], 
                           parcel_profile, 
                           color = "lightcyan")
            
            skew.shade_cape(sounding_df["pressure"], 
                            sounding_df["temperature"], 
                            parcel_profile, 
                            color = "mistyrose")


            # Good bounds for aspect ratio

            skew.ax.set_xlim(-25, max_axis_temp)
            skew.ax.set_title(valid_time + "  (" + local_time+")", fontsize=15)


            if lcl_pressure:
                skew.ax.plot(lcl_temperature, 
                             lcl_pressure, 
                             marker          =     "_", 
                             color           = 'orange', 
                             markersize      =      30, 
                             markeredgewidth =       3)

            if lfc_pressure:
                skew.ax.plot(lfc_temperature, 
                             lfc_pressure, 
                             marker          =     "_", 
                             color           = 'brown', 
                             markersize      =      30, 
                             markeredgewidth =       3)

            if el_pressure:
                skew.ax.plot(el_temperature, 
                             el_pressure, 
                             marker          =    "_", 
                             color           = 'blue', 
                             markersize      =     30, 
                             markeredgewidth =      3)
                
            skew.ax.spines["right"].set_visible(False)
            skew.ax.spines[  "top"].set_color("white")
            skew.ax.set_xlabel('Temperature (\N{DEGREE CELSIUS})')
            skew.ax.set_ylabel("Isobaric Height (hPa)")

            skew_box          = skew.ax.get_position()
            skew_box_x_start  = skew_box.x0
            skew_box_y_start  = skew_box.y0
            skew_box_x_end    = skew_box.x1
            skew_box_y_end    = skew_box.y1
            skew_box_x_length = skew_box.x1 - skew_box.x0
            skew_box_y_length = skew_box.y1 - skew_box.y0
        

            #
            ###################################################



            
            ###################################################
            #
            # Create a Cloud Water Profile
            #
         
            if (True):
                
                axclouds = SkewT(fig      =       fig,
                                 rotation =         0, 
                                 aspect = 'auto',
                                 #aspect   =    (90./np.log10(1050.-100)) / (np.ceil(clouds_maxx)/np.log10(1050.-100)) ,
                                 rect     = [0.77, skew_box_y_start, 0.28, skew_box_y_length])
                
                axclouds.ax.set_title('Hodograph & Moisture Profile', fontsize=15)
                axclouds.ax.set_xlabel('Mixing Ratio (g/kg)')
                axclouds.ax.set_ylabel("")
                axclouds.ax.set_xlim(0, np.ceil(clouds_maxx))
                axclouds.ax.xaxis.set_major_locator(MultipleLocator(1))
                axclouds.ax.xaxis.set_units(units("g/kg"))
                axclouds.ax.spines["right"].set_visible(False)
                axclouds.ax.spines[  "top"].set_color("white")
                
                
                axclouds.plot(sounding_df["pressure"], sounding_df["qv"], color="greenyellow")
                axclouds.plot(sounding_df["pressure"], sounding_df["qc"], color="darkgrey")            
                axclouds.plot(sounding_df["pressure"], sounding_df["qi"], color="cyan")            
                axclouds.plot(sounding_df["pressure"], sounding_df["qr"], color="darkgreen")            
                axclouds.plot(sounding_df["pressure"], sounding_df["qs"], color="blue")            
                axclouds.plot(sounding_df["pressure"], sounding_df["qg"], color="red")  


                axclouds.ax.legend(["$q_v$","$q_c$","$q_i$","$q_r$","$q_s$","$q_g$"],
                                  loc = 'center right',
                                  frameon = False)
            
            axclouds_box          = axclouds.ax.get_position()
            axclouds_box_x_start  = axclouds_box.x0
            axclouds_box_y_start  = axclouds_box.y0
            axclouds_box_x_end    = axclouds_box.x1
            axclouds_box_y_end    = axclouds_box.y1
            axclouds_box_x_length = axclouds_box.x1 - axclouds_box.x0
            axclouds_box_y_length = axclouds_box.y1 - axclouds_box.y0


            #
            ###################################################             

        
            ###################################################
            #
            # Create a hodograph
            #
            
            do_hodo = True
            
            if (do_hodo):
                
                
                
                mask = hodo_frame["agl"] <= 10 * units.km

                axhodo = fig.add_axes(rect = [0.75, 
                                              skew_box_y_start+skew_box_y_length-0.32, 
                                              0.3, 
                                              0.3],
                                      zorder=1001)
                axhodo.zorder=1001

                h = Hodograph(ax              = axhodo, 
                              component_range = 40.)

                h.add_grid(increment=10)
                axhodo.spines["right"].set_visible(False)
                axhodo.spines[  "top"].set_visible(False)


                cmh = h.plot_colormapped(u    = hodo_frame["u_wind"][mask], 
                                         v    = hodo_frame["v_wind"][mask], 
                                         c    = hodo_frame["agl"][mask].to(units.km), 
                                         cmap = "jet")
                cmh.set_clim(0, 10)
                #cmh2 = h.plot(hodo_extra["u_wind"], hodo_extra["v_wind"], 
                #              color = Mines_Blue, 
                #              marker = "o",
                #              linestyle="None")
                

                fig.patches.extend([plt.Rectangle((0.73,skew_box_y_start+skew_box_y_length-0.35),0.3,0.355,
                                  fill=True, color='w', alpha=0.99, zorder=1000,
                                  transform=fig.transFigure, figure=fig)])

                cbhodo = plt.colorbar(mappable = cmh, 
                                      ax = axhodo,
                                      #labelpad = -2,
                                      orientation = "horizontal", 
                                      label       = 'Height (km AGL)')
                
                #cbhodo.ax.scatter(hodo_extra["agl"]/1000,
                #                 (np.zeros(len(hodo_extra["agl"]))+0.5),
                #                 color = "black")
            
                cbhodo.ax.zorder = 1001

                axhodo.set_xlabel("u (m s$^{-1}$)",labelpad=-2)
                axhodo.set_ylabel("v (m s$^{-1}$)",labelpad=-2)
            #
            ###################################################  
            
            
            #plt.tight_layout()

            #plt.subplots_adjust(top=0.93)    
            
            
            
            ###################################################
            #
            # Add Status Bars 
            #
            
            rect1 = patches.Rectangle(xy        = (0, 0-0.01/2),
                                      width     = percent_done,
                                      height    = 0.01, 
                                      edgecolor = Mines_Blue, 
                                      facecolor = Mines_Blue,
                                      transform = skew.ax.transAxes)
            skew.ax.add_patch(rect1)

            rect2 = patches.Rectangle(xy        = (0, 0-0.01/2),
                                      width     = percent_done,
                                      height    = 0.015, 
                                      edgecolor = Mines_Blue, 
                                      facecolor = Mines_Blue,
                                      transform = axhodo.transAxes, zorder = 1001)
            axhodo.add_patch(rect2)

            rect3 = patches.Rectangle(xy        = (0, 0-0.01/2),
                                      width     = percent_done,
                                      height    = 0.01, 
                                      edgecolor = Mines_Blue, 
                                      facecolor = Mines_Blue,
                                      transform = axclouds.ax.transAxes)
            axclouds.ax.add_patch(rect3)
            
            #######
            #
            # Walking Clock
            #
            #######
            
            plot_box          = skew.ax.get_position()
            plot_box_x_start  = plot_box.x0
            plot_box_y_start  = plot_box.y0
            plot_box_x_end    = plot_box.x1
            plot_box_y_end    = plot_box.y1
            plot_box_x_length = plot_box.x1 - plot_box.x0
            plot_box_y_length = plot_box.y1 - plot_box.y0

            size_of_clock = 0.05

            x_clock = percent_done*(plot_box_x_end-plot_box_x_start) + plot_box_x_start - size_of_clock/2

            y_clock = plot_box_y_start-size_of_clock/2
            
            x_dow = percent_done
            y_dow = size_of_clock    
            
            skew.ax.annotate(dow+"-"+local_time_zone, 
                             [x_dow,y_dow],
                             horizontalalignment = "center",
                             verticalalignment   = "center",
                             xycoords            = 'axes fraction')


            
            axins = fig.add_axes(rect     =    [x_clock,
                                                y_clock,
                                                size_of_clock,
                                                size_of_clock],
                                 projection =  "polar")
            
            plt.setp(axins.get_yticklabels(), visible=False)
            plt.setp(axins.get_xticklabels(), visible=False)
            axins.spines['polar'].set_visible(False)
            axins.set_ylim(0,1)
            axins.set_theta_zero_location('N')
            axins.set_theta_direction(-1)
            axins.set_facecolor(Clock_BgndC)
            axins.grid(False)

            axins.plot([angles_h,angles_h], [0,0.60], color=Clock_Color, linewidth=1)
            axins.plot([angles_m,angles_m], [0,0.95], color=Clock_Color, linewidth=1)
            axins.plot(circle_theta, circle_radius, color="darkgrey", linewidth=1)
            

            

            #
            ###################################################
    
            ###################################################
            #
            # Add Descriptive Statistics
            #   
            
            xtext_start =  0.22;  0.115
            ytext_start =  0.9
            text_deltay =  0.02
            xtext_numbr =  xtext_start + 0.005

            
            plt.figtext( xtext_start, ytext_start - text_deltay* 0, 'LCL Height:', horizontalalignment = "right")
            plt.figtext( xtext_numbr, ytext_start - text_deltay* 0, '{0:.1f} m'.format(lcl_hgt.magnitude))
            plt.figtext( xtext_start, ytext_start - text_deltay* 1, 'LFC Height:', horizontalalignment = "right")
            plt.figtext( xtext_numbr, ytext_start - text_deltay* 1, '{0:.1f} m'.format(lfc_hgt.magnitude))
            plt.figtext( xtext_start, ytext_start - text_deltay* 2, 'MLLR:', horizontalalignment = "right")
            plt.figtext( xtext_numbr, ytext_start - text_deltay* 2, '{0:.1f} K'.format(lr_700_500.magnitude))
            plt.figtext( xtext_start, ytext_start - text_deltay* 3, 'SBCAPE:', horizontalalignment = "right")
            plt.figtext( xtext_numbr, ytext_start - text_deltay* 3, '{0:.1f} J/kg'.format(sbcape.magnitude))
            plt.figtext( xtext_start, ytext_start - text_deltay* 4, 'SBCIN:', horizontalalignment = "right")
            plt.figtext( xtext_numbr, ytext_start - text_deltay* 4, '{0:.1f} J/kg'.format(sbcin.magnitude))
            plt.figtext( xtext_start, ytext_start - text_deltay* 5, 'MLCAPE:', horizontalalignment = "right")
            plt.figtext( xtext_numbr, ytext_start - text_deltay* 5, '{0:.1f} J/kg'.format(mlcape.magnitude))
            plt.figtext( xtext_start, ytext_start - text_deltay* 6, 'MLCIN:', horizontalalignment = "right")
            plt.figtext( xtext_numbr, ytext_start - text_deltay* 6, '{0:.1f} J/kg'.format(mlcin.magnitude))
            plt.figtext( xtext_start, ytext_start - text_deltay* 7, 'MUCAPE:', horizontalalignment = "right")
            plt.figtext( xtext_numbr, ytext_start - text_deltay* 7, '{0:.1f} J/kg'.format(mucape.magnitude))
            plt.figtext( xtext_start, ytext_start - text_deltay* 8, 'Shear 0-1 km:', horizontalalignment = "right")
            plt.figtext( xtext_numbr, ytext_start - text_deltay* 8, '{0:.1f} m/s'.format(shear01.magnitude))
            plt.figtext( xtext_start, ytext_start - text_deltay* 9, 'Shear 0-6 km:', horizontalalignment = "right")
            plt.figtext( xtext_numbr, ytext_start - text_deltay* 9, '{0:.1f} m/s'.format(shear06.magnitude))
            plt.figtext( xtext_start, ytext_start - text_deltay*10, 'SRH 0-1 km:', horizontalalignment = "right")
            plt.figtext( xtext_numbr, ytext_start - text_deltay*10, '{0:.1f} m\u00b2/s\u00b2'.format(srh_01.magnitude))
            plt.figtext( xtext_start, ytext_start - text_deltay*11, 'SRH 0-3 km:', horizontalalignment = "right")
            plt.figtext( xtext_numbr, ytext_start - text_deltay*11, '{0:.1f} m\u00b2/s\u00b2'.format(srh_03.magnitude))

            #
            ###################################################           
    
   
            ###################################################
            #
            # Close SkewT
            #
            
            plt.savefig(graphics_directory + sounding_file_name_png,
                        facecolor   = 'white', 
                        transparent =   False,
                        bbox_inches = 'tight', 
                        pad_inches  =       0)
            
            #print(1/0)
                        
            plt.close('all')
            
        
            #
            ###################################################

        #
        ### ### ### ### ### ### ### ### ### ### ### ### ####
        ####################################################
        
        
        ####################################################
        #
        # making gifs
        #

        png_file_name    = "wrfout_dxx_" + model_start_date_YYYY_MM_DD_HH + "_F??_SKEWT_" + station_id + ".png"
        gif_file_name    = "wrfout_dxx_" + model_start_date_YYYY_MM_DD_HH + "_Fxx_SKEWT_" + station_id + ".gif"

        
        print("creating " + WRF_OVERALL_DIR + "./processing_skewt"+ str(domain).zfill(1) + "_gif.sh")
        with open(WRF_OVERALL_DIR + "./processing_skewt"+ str(domain).zfill(1) + "_gif.sh", 'w') as f:
            print("#!/bin/bash", file =  f)
            print(". ~/.bashrc", file =  f)
            print("ulimit -s unlimited", file = f)
            if intel:
                print(". /opt/intel/oneapi/setvars.sh", file = f)
            print("cd " + WRF_OVERALL_DIR, file =  f) 
            print("convert -delay 50 " + graphics_directory + png_file_name + " " + graphics_directory + gif_file_name, file =  f) 
            print("echo MAIN:SKEWT"+ str(domain).zfill(1) + "::: We^re Outahere Like Vladimir", file =  f) 

        os.system("chmod a+x " + WRF_OVERALL_DIR + "./processing_skewt"+ str(domain).zfill(1) + "_gif.sh")
        os.system(WRF_OVERALL_DIR + "./processing_skewt"+ str(domain).zfill(1) + "_gif.sh > ./processing_skewt"+ str(domain).zfill(1) + "_gif." + model_start_date_YYYY_MM_DD_HH + ".LOG 2>&1 ")
        os.system("date")
        print()
    



        
        #
        ####################################################  
        
        

            
    #
    ####################################################
    ####################################################

#
####################################################
####################################################
####################################################


# ## Ending ScriptMines_Blue

# In[ ]:


####################################################
####################################################
####################################################
#
# End of Script
#

print("End Sounding Plotting Scriot")

#
####################################################
####################################################
####################################################


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




