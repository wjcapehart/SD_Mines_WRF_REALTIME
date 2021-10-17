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


import matplotlib.pyplot as plt
import pandas            as pd
import xarray            as xr
import pint_xarray       as px

import netCDF4           as nc4

import wrf               as wrf


import matplotlib.gridspec as gridspec

import metpy.calc        as mpcalc
from metpy.plots import SkewT, Hodograph
from metpy.units import units, pandas_dataframe_to_unit_arrays


import metpy.calc  as mpcalc

from metpy.units import units

import seaborn           as sns

import timezonefinder    as tzf
import pytz as pytz
import socket as socket

import matplotlib.font_manager as fm
import matplotlib as mpl

sns.set_theme(style="ticks")

if (platform.system() != "Darwin"):
    path = '/usr/share/fonts/truetype/open-sans/OpenSans-Regular.ttf'
    print("          Enabling OpenSans " + path)
    prop = fm.FontProperties(fname=path)
    print( prop.get_name() )
    mpl.rcParams['font.family'] = prop.get_name()[0]
    # print("mpl.rcParams['font.family'] " + mpl.rcParams['font.family'])

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
# File Organization
#

beta_on       = 0
max_domains   = 3
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


            sounding_file_name_png = "wrfout_dxx_" + model_start_date_YYYY_MM_DD_HH + "_F" + str(t).zfill(2) + "_SKEWT_" + station_id + ".png"
            
            print("     -- " + sounding_file_name_png)
    
            sounding_df = pd.DataFrame({"pressure"    :   isobar_hgt_4d[t, :, wrf_y, wrf_x],
                                        "temperature" :  temperature_4d[t, :, wrf_y, wrf_x],
                                        "dewpoint"    :    dew_point_4d[t, :, wrf_y, wrf_x],
                                        "u_wind"      :        uwind_4d[t, :, wrf_y, wrf_x],
                                        "v_wind"      :        vwind_4d[t, :, wrf_y, wrf_x],
                                        "height"      :       height_4d[t, :, wrf_y, wrf_x],
                                        "agl"         :   height_agl_4d[t, :, wrf_y, wrf_x]})

            units_df = {"pressure"    :  "hPa",
                        "temperature" : "degC",
                        "dewpoint"    : "degC",
                        "u_wind"      :   "kt",
                        "v_wind"      :   "kt",
                        "height"      :    "m",
                        "agl"         :    "m"}


            sounding_df = pandas_dataframe_to_unit_arrays(sounding_df, 
                                                          column_units = units_df)

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

            fig = plt.figure(figsize=(9, 9))

            fig.suptitle(station_name + "; Model Run " + wrf_skewt_time + "; WRF Domain " + str(domain).zfill(2), 
                         fontsize=20)

            # Grid for plots

            gs   = gridspec.GridSpec(3, 3)

            skew = SkewT(fig, 
                         rotation =        45,  
                         subplot  = gs[:, :2])

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
                      color = 'black', 
                      linewidth = 0.5)

            # Mask barbs to be below 100 hPa only

            mask = (sounding_df["pressure"] >= 100 * units.hPa)
            for i in range(mask.size):
                if (sounding_df["pressure"][i] > 500 * units.hPa) :
                    if (i%2 == 1):
                        mask[i] = False

            skew.plot_barbs(sounding_df["pressure"][mask], 
                            sounding_df[  "u_wind"][mask], 
                            sounding_df[  "v_wind"][mask])

            skew.ax.set_ylim(1000, 100)

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

            skew.ax.set_xlabel = "Temperature (C)"
            skew.ax.set_ylabel = "Isobaric Height (hPa)"


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

            #
            ###################################################

            ###################################################
            #
            # Add Heights
            #



            #
            ###################################################

            ###################################################
            #
            # Create a hodograph
            #

            mask = sounding_df["agl"] <= 10 * units.km

            ax = fig.add_subplot(gs[0, -1])

            h = Hodograph(ax, component_range=40.)

            h.add_grid(increment=10)

            cmh = h.plot_colormapped(sounding_df["u_wind"][mask], sounding_df["v_wind"][mask], sounding_df[   "agl"][mask].to(units.km), cmap = "jet")
            cmh.set_clim(0, 10)

            plt.colorbar(cmh, orientation = "horizontal", label       = 'Height (km AGL)')

            plt.figtext( 0.65, 0.60, 'LCL Height:')
            plt.figtext( 0.80, 0.60, '{}'.format(lcl_hgt))
            plt.figtext( 0.65, 0.58, 'LFC Height:')
            plt.figtext( 0.80, 0.58, '{}'.format(lfc_hgt))
            plt.figtext( 0.65, 0.56, 'MLLR:')
            plt.figtext( 0.80, 0.56, '{}'.format(lr_700_500))
            plt.figtext( 0.65, 0.54, 'SBCAPE:')
            plt.figtext( 0.80, 0.54, '{0} J/kg'.format(sbcape))
            plt.figtext( 0.65, 0.52, 'SBCIN:')
            plt.figtext( 0.80, 0.52, '{}'.format(sbcin))
            plt.figtext( 0.65, 0.50, 'MLCAPE:')
            plt.figtext( 0.80, 0.50, '{}'.format(mlcape))
            plt.figtext( 0.65, 0.48, 'MLCIN:')
            plt.figtext( 0.80, 0.48, '{}'.format(mlcin))
            plt.figtext( 0.65, 0.46, 'MUCAPE:')
            plt.figtext( 0.80, 0.46, '{}'.format(mucape))
            plt.figtext( 0.65, 0.44, 'Shear 0-1 km:')
            plt.figtext( 0.80, 0.44, '{}'.format(shear01))
            plt.figtext( 0.65, 0.42, 'Shear 0-6 km:')
            plt.figtext( 0.80, 0.42, '{}'.format(shear06))
            plt.figtext( 0.65, 0.40, 'SRH 0-1 km:')
            plt.figtext( 0.80, 0.40, '{}'.format(srh_01))
            plt.figtext( 0.65, 0.38, 'SRH 0-3 km:')
            plt.figtext( 0.80, 0.38, '{}'.format(srh_03))

            plt.tight_layout()

            plt.subplots_adjust(top=0.93)    
            
            plt.savefig(graphics_directory + sounding_file_name_png)
            
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
        os.system("convert -delay 25 " + graphics_directory + png_file_name + " " + graphics_directory + gif_file_name)

        #
        ####################################################  
        
        

            
    #
    ####################################################
    ####################################################

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

print("End Sounding Plotting Scriot")

#
####################################################
####################################################
####################################################


# In[ ]:





# In[ ]:




