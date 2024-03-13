#!/usr/bin/env python
# coding: utf-8

# # Plotting Maps

# ## Libraries

# In[1]:


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
import matplotlib.patches as patches

import pandas            as pd
import xarray            as xr
import pint_xarray       as px

import netCDF4           as nc4


import cartopy.crs       as ccrs
import cartopy.feature   as cfeature

import pyproj            as pyproj

import matplotlib.gridspec as gridspec

import metpy.calc        as mpcalc

from metpy.units import units, pandas_dataframe_to_unit_arrays
from metpy.plots import colortables,  USCOUNTIES


import metpy.calc  as mpcalc

from metpy.units import units

import timezonefinder    as tzf
import pytz as pytz

import socket as socket
import matplotlib.font_manager as fm
import matplotlib as mpl


from      joblib import Parallel, delayed

from metpy.plots import colortables
from datetime    import timezone



import wrf               as wrf


colorbar_pad    = 0.012
colorbar_shrink = 0.8
#
####################################################
####################################################
####################################################


# In[2]:


####################################################
####################################################
####################################################
#
# Mines Colors and Fonts
#

Mines_Blue = "#002554"

#plt.rcParams['font.family'] = 'OpenSans'

plt.rcParams.update({'text.color'      : Mines_Blue,
                     'axes.labelcolor' : Mines_Blue,
					 'axes.edgecolor'  : Mines_Blue,
					 'xtick.color'     : Mines_Blue,
					 'ytick.color'     : Mines_Blue})


#
####################################################
####################################################
####################################################


# In[3]:


####################################################
####################################################
####################################################
#
# MetPy Color Tables
#

nws_dbz_colors = colortables.get_colortable("NWSReflectivity")
nws_rainbow    = colortables.get_colortable("ir_rgbv")

## radar colormap

norm_radar, cmap_radar = colortables.get_with_steps("NWSStormClearReflectivity", 
                                        -20, 
                                        0.5)





clearsky_dbz_values = np.arange(-28, 28.1, 2)
stormy_dbz_values   = np.arange(  5, 75.1, 5)

#
####################################################
####################################################
####################################################


# In[4]:


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
    "#8B008B",  # 8.00 - 10.00 inches
    "#ff00ff",  # 6.00 - 8.00 inches
	"#FFCCFF",
    "#888888"]  # 10.00+

precip_colormap = mpl.colors.ListedColormap(colors = nws_precip_colors)

precip_levels_hrly = [   0.01,   0.05,   0.10,  0.25,    
                         0.50,   0.75,   1.00,  1.25,   
                         1.50,   1.75,   2.00,  2.25,  
                         2.50,   3.00,   4.00,  5.00] # in Inches!!!


precip_levels_full = [   0.01,   0.05,   0.10,  0.25,    
                         0.50,   0.75,   1.00,  1.50,   
                         2.00,   2.50,   3.00,  4.00,  
                         5.00,   6.00,   8.00, 10.00] # in Inches!!!

snow_levels_hrly   = [   0.10,   0.25,   0.50,  0.75,    
                         1.00,   1.25,   1.50,  1.75,   
                         2.00,   2.25,   2.50,  2.75,  
                         3.00,   4.00,   5.00,  6.00] # in Inches!!!

snow_levels_full   = [   0.10,   0.25,   0.50,  0.75,    
                         1.00,   2.00,   3.00,  4.00,   
                         5.00,   6.00,  12.00, 18.00,  
                        24.00,  30.00,  36.00, 48.00] # in Inches!!!

#
###################################################


# In[5]:


###################################################
#
# BOM Temperature Color Table.
#

###################################################
#
# BOM Temperature Color Table.
#


bom_temp_colors  = ['#ccebff', 
					'#b3cde3',
					'#a0b2d4', 
					'#8c96c6', 
					'#8856a7', 
					'#810f7c', 
					'#081d58', 
					'#132778', 
					'#253494', 
					'#23479e', 
					'#225ca7', 
					'#1f76b4', 
					'#1d91c0', 
					'#2ca7c5', 
					'#43b5c5', 
					'#63c8c5', 
					'#7fcdbb', 
					'#98dca6', 
					'#c7e9b4',
					'#dcf2c6', 
					'#edf8d9', 
					'#f5fcd3', 
					'#fcffcc', 
					'#fff9b6', 
					'#ffeda0', 
					'#ffe48b', 
					'#fed976', 
					'#fec761', 
					'#feb24c', 
					'#fea044', 
					'#fd8d3c', 
					'#fd6e33', 
					'#fc4e2a', 
					'#f23120', 
					'#e31a1c', 
					'#d20b20', 
					'#bd0026', 
					'#9f0027', 
					'#73001f', 
					'#4c0019', 
					'#19000d', 
					'#4c0099', 
					'#cc17cc',
				    '#ffffff']

bom_colormap = mpl.colors.ListedColormap(colors = bom_temp_colors)


temp_range_degF =  np.linspace(-20, 115, 44-2)
delta_t_scale   =  temp_range_degF[1]-temp_range_degF[0]
temp_range_degF =  np.linspace(-20-delta_t_scale, 115+delta_t_scale, 44)

bom_colorbar_ticks = np.arange(-20,116,5)

dewpoint_levels_degF = np.linspace(30,70,41) # in DegF


#
###################################################


# ## Parallel Code Area
# 
# ### Making Time Series Maps

# In[7]:


####################################################

def plot_time_series_maps_func(t):
	
	
    import matplotlib.font_manager as font_manager


    Mines_Blue = "#002554"


    plt.rcParams.update({'text.color'      : Mines_Blue,
                         'axes.labelcolor' : Mines_Blue,
                         'xtick.color'     : Mines_Blue,
                         'ytick.color'    : Mines_Blue})



    ####################################################
    #
    # Establish Times and Labels
    #

    valid_time      = pd.to_datetime(wrf_time_steps[t].values).tz_localize(tz="UTC").strftime("%Y-%m-%d %H %Z")
    local_time      = pd.to_datetime(wrf_time_steps[t].values).tz_localize(tz="UTC").tz_convert(tz=tz).strftime("%Y-%m-%d %H %Z")
    local_time_zone = pd.to_datetime(wrf_time_steps[t].values).tz_localize(tz="UTC").tz_convert(tz=tz).strftime("%Z")
    dow             = pd.to_datetime(wrf_time_steps[t].values).tz_localize(tz="UTC").tz_convert(tz=tz).strftime("%a")

    
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

    
    model_run_label    = "Model Run " + wrf_skewt_time + "; WRF Domain " + str(domain).zfill(2)

    print(valid_time, "     -- " + model_run_label, " :: ", (percent_done*100) )

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

    fig.suptitle(model_run_label,fontsize="x-large")

    ax1 = fig.add_subplot(1,  # nrows
                          1,  # ncols 
                          1,  # index of figure you're installing
                          projection = cart_proj) # cartopy CRS Projection

    ax1.set_xlim(  west_east_range)
    ax1.set_ylim(south_north_range)
    ax1.set_frame_on(False)
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)


    ax1.add_feature(cfeature.LAKES, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.RIVERS, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.BORDERS, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.COASTLINE, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')

    ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                 'admin_1_states_provinces_lines',
                                                 '50m',
                                                 linewidths = 0.75,
                                                 facecolor  = 'none',
                                                 edgecolor  = Mines_Blue))

    if (domain > 1) :
        ax1.add_feature(USCOUNTIES, linewidths=0.5,edgecolor  = Mines_Blue,facecolor='none')

    line_contour_levels = np.arange(start =  900, 
                                    stop  = 1200,
                                    step  =    4)

    line_contour_map = ax1.contour(lon2d, 
                                  lat2d, 
                                  mslp_maps.isel(Time=t),
                                  transform  = ccrs.PlateCarree(),
                                  colors     = "white",
                                  linewidths = 1,
                                  levels     = line_contour_levels,
                                  extend     = "both")

    ax1.clabel(line_contour_map, inline=True, fontsize=8)

    filled_cm = t2m_maps.isel(Time=t).plot.imshow(cmap          = bom_colormap,
												  alpha         = alpha2d,
												  interpolation = "bilinear",
												  vmin          = temp_range_degF.min(),
												  vmax          = temp_range_degF.max(),
												  add_colorbar  = False)
	
    cb = plt.colorbar(mappable = filled_cm, 
					  label    = r"2-m Temperature (°F)",
					  shrink   = colorbar_shrink, 
					  pad      = colorbar_pad,
					  ticks    = bom_colorbar_ticks,
					  format   = '%d')
	
    cb.outline.set_color(Mines_Blue)

    gap = 10
    plt.barbs(wrf.to_np(lon2d[::gap,::gap]),
              wrf.to_np(lat2d[::gap,::gap]),
              wrf.to_np(u10_maps.isel(Time=t)[::gap,::gap]),
              wrf.to_np(v10_maps.isel(Time=t)[::gap,::gap]),
              transform = ccrs.PlateCarree(), 
              length=5)
    
    
    ##################################
    # 
    # Progress Bar
    #
    
    plot_box         = ax1.get_position()
    plot_box_x_start = plot_box.x0
    plot_box_y_start = plot_box.y0
    plot_box_x_end   = plot_box.x1
    plot_box_y_end   = plot_box.y1



    rect1 = patches.Rectangle(xy        = (0, 0),
                              width     = percent_done,
                              height    = 0.01, 
                              edgecolor = Mines_Blue, 
                              facecolor = Mines_Blue,
                              transform = ax1.transAxes)
    ax1.add_patch(rect1)
    
    #
    # Walking Clock
    #
 
    circle_theta  = np.deg2rad(np.arange(0,360,0.01))
    circle_radius = circle_theta * 0 + 1

    if (hour > 12) :
        hour = hour - 12

    angles_h = 2*np.pi*hour/12+2*np.pi*minute/(12*60)+2*second/(12*60*60)
    angles_m = 2*np.pi*minute/60+2*np.pi*second/(60*60)



    size_of_clock = 0.07

    x_clock = percent_done*(plot_box_x_end-plot_box_x_start) + plot_box_x_start - size_of_clock/2
    y_clock = plot_box_y_start-size_of_clock/2+0.01

    x_dow = percent_done
    y_dow = size_of_clock

    ax1.annotate(dow+"-"+local_time_zone, 
                 [x_dow,y_dow],
                 horizontalalignment="center",
                verticalalignment="center",
                xycoords='axes fraction')


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
    axins.plot(circle_theta, circle_radius,   color=Mines_Blue,  linewidth=1)

    #
    ##################################
            
    # plt.show()
    ax1.set_frame_on(False)
    ax1.set_title(valid_time + "  (" + local_time+")")

    fig.savefig(fname       = fig_dir_name + file_name,
				facecolor   =                  'white', 
				transparent =                    False)

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

    fig = plt.figure(figsize=figure_domain_size)

    fig.suptitle(model_run_label,fontsize="x-large")

    ax1 = fig.add_subplot(1,  # nrows
                          1,  # ncols 
                          1,  # index of figure you're installing
                          projection = cart_proj) # cartopy CRS Projection

    ax1.set_xlim(  west_east_range)
    ax1.set_ylim(south_north_range)
    ax1.set_frame_on(False)
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)

    ax1.add_feature(cfeature.LAKES, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.RIVERS, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.BORDERS, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.COASTLINE, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')

    ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                 'admin_1_states_provinces_lines',
                                                 '50m',
                                                 linewidths = 0.75,
                                                 facecolor  = 'none',
                                                 edgecolor  = Mines_Blue))

    if (domain > 1) :
        ax1.add_feature(USCOUNTIES, 
                        linewidths=0.5,
                        edgecolor  = Mines_Blue,
                        facecolor='none')

    ax1.streamplot(wrf.to_np(lon2d),
                   wrf.to_np(lat2d), 
                   wrf.to_np(u10_maps.isel(Time=t)),
                   wrf.to_np(v10_maps.isel(Time=t)),
                   transform = ccrs.PlateCarree(),
                   color = Mines_Blue)

    filled_cm = m10_maps.isel(Time=t).plot.imshow(cmap          = mpl.cm.Blues,
												  alpha         = alpha2d,
												  interpolation = "bilinear",
												  vmin          = 0,
												  vmax          = 41,
												  add_colorbar  = False)

    cb = plt.colorbar(filled_cm, 
					  label  = r"10-m Wind Speed (kts)",
					  shrink = colorbar_shrink, 
					  pad    = colorbar_pad)
    cb.outline.set_color(Mines_Blue)
        
    
    
    ##################################
    # 
    # Progress Bar
    #
    
    plot_box         = ax1.get_position()
    plot_box_x_start = plot_box.x0
    plot_box_y_start = plot_box.y0
    plot_box_x_end   = plot_box.x1
    plot_box_y_end   = plot_box.y1



    rect1 = patches.Rectangle(xy        = (0, 0),
                              width     = percent_done,
                              height    = 0.01, 
                              edgecolor = Mines_Blue, 
                              facecolor = Mines_Blue,
                              transform = ax1.transAxes)
    ax1.add_patch(rect1)
    
    #
    # Walking Clock
    #

    circle_theta  = np.deg2rad(np.arange(0,360,0.01))
    circle_radius = circle_theta * 0 + 1

    if (hour > 12) :
        hour = hour - 12

    angles_h = 2*np.pi*hour/12+2*np.pi*minute/(12*60)+2*second/(12*60*60)
    angles_m = 2*np.pi*minute/60+2*np.pi*second/(60*60)



    size_of_clock = 0.07

    x_clock = percent_done*(plot_box_x_end-plot_box_x_start) + plot_box_x_start - size_of_clock/2
    y_clock = plot_box_y_start-size_of_clock/2+0.01

    x_dow = percent_done
    y_dow = size_of_clock

    ax1.annotate(dow+"-"+local_time_zone, 
                 [x_dow,y_dow],
                 horizontalalignment="center",
                verticalalignment="center",
                xycoords='axes fraction')


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
    axins.plot(circle_theta, circle_radius,   color=Mines_Blue,  linewidth=1)

    #
    ##################################
            
    # plt.show()
    ax1.set_frame_on(False)
    ax1.set_title(valid_time + "  (" + local_time+")")

    fig.savefig(fname       = fig_dir_name + file_name,
				facecolor   =                  'white', 
				transparent =                    False)

    plt.close('all')    

    #
    ####################################################

    ####################################################
    #
    # Dew Points
    #

    v_name       = "DEWP"
    fig_dir_name = graphics_directory + "/" + v_name + "/"
    file_name    = "wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH + "_F" + str(t).zfill(2) + "_MAP_" + v_name + ".png"

    print(fig_dir_name + file_name)

    fig = plt.figure(figsize=figure_domain_size)

    fig.suptitle(model_run_label,fontsize="x-large")

    ax1 = fig.add_subplot(1,  # nrows
                          1,  # ncols 
                          1,  # index of figure you're installing
                          projection = cart_proj) # cartopy CRS Projection

    ax1.set_xlim(  west_east_range)
    ax1.set_ylim(south_north_range)
    ax1.set_frame_on(False)
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)


    ax1.add_feature(cfeature.LAKES, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.RIVERS, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.BORDERS, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.COASTLINE, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')

    ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                 'admin_1_states_provinces_lines',
                                                 '50m',
                                                 linewidths = 0.75,
                                                 facecolor  = 'none',
                                                 edgecolor  = Mines_Blue))

    if (domain > 1) :
        ax1.add_feature(USCOUNTIES, 
                        linewidths=0.5,
                        edgecolor  = Mines_Blue,
                        facecolor='none')

    
    lw = 5 * m10_maps.isel(Time=t) / np.nanmax(m10_maps.isel(Time=t).values)
    ax1.streamplot(wrf.to_np(lon2d),
                   wrf.to_np(lat2d), 
                   wrf.to_np(u10_maps.isel(Time=t)),
                   wrf.to_np(v10_maps.isel(Time=t)),
                   linewidth = lw.values,
                   transform = ccrs.PlateCarree(),
                   color = Mines_Blue)

    filled_cm = td2m_maps.isel(Time=t).plot.imshow(cmap          = mpl.cm.YlGn,
												   alpha         = alpha2d,
												   interpolation = "bilinear",
												   vmin          = dewpoint_levels_degF.min(),
												   vmax          = dewpoint_levels_degF.max(),
												   add_colorbar  = False)

    cb = plt.colorbar(filled_cm, 
				 label  = r"2-m Dew Point Temperatuere (°F)",
				 shrink = colorbar_shrink, 
				 pad    = colorbar_pad)
    cb.outline.set_color(Mines_Blue)
 
    
    
    ##################################
    # 
    # Progress Bar
    #
    
    plot_box         = ax1.get_position()
    plot_box_x_start = plot_box.x0
    plot_box_y_start = plot_box.y0
    plot_box_x_end   = plot_box.x1
    plot_box_y_end   = plot_box.y1



    rect1 = patches.Rectangle(xy        = (0, 0),
                              width     = percent_done,
                              height    = 0.01, 
                              edgecolor = Mines_Blue, 
                              facecolor = Mines_Blue,
                              transform = ax1.transAxes)
    ax1.add_patch(rect1)
    
    #
    # Walking Clock
    #

    circle_theta  = np.deg2rad(np.arange(0,360,0.01))
    circle_radius = circle_theta * 0 + 1

    if (hour > 12) :
        hour = hour - 12

    angles_h = 2*np.pi*hour/12+2*np.pi*minute/(12*60)+2*second/(12*60*60)
    angles_m = 2*np.pi*minute/60+2*np.pi*second/(60*60)



    size_of_clock = 0.07

    x_clock = percent_done*(plot_box_x_end-plot_box_x_start) + plot_box_x_start - size_of_clock/2
    y_clock = plot_box_y_start-size_of_clock/2+0.01

    x_dow = percent_done
    y_dow = size_of_clock

    ax1.annotate(dow+"-"+local_time_zone, 
                 [x_dow,y_dow],
                 horizontalalignment="center",
                verticalalignment="center",
                xycoords='axes fraction')


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
    axins.plot(circle_theta, circle_radius,   color=Mines_Blue,  linewidth=1)

    #
    ##################################
            
    # plt.show()
    ax1.set_frame_on(False)
    ax1.set_title(valid_time + "  (" + local_time+")")

    fig.savefig(fname       = fig_dir_name + file_name,
				facecolor   =                  'white', 
				transparent =                    False)

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

    fig = plt.figure(figsize=figure_domain_size)

    fig.suptitle(model_run_label,fontsize="x-large")

    ax1 = fig.add_subplot(1,  # nrows
                          1,  # ncols 
                          1,  # index of figure you're installing
                          projection = cart_proj) # cartopy CRS Projection

    ax1.set_xlim(  west_east_range)
    ax1.set_ylim(south_north_range)
    ax1.set_frame_on(False)
    ax1.set_xlim(  west_east_range)
    ax1.set_ylim(south_north_range)
    ax1.set_frame_on(False)
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)

    ax1.add_feature(cfeature.LAKES, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.RIVERS, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.BORDERS, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.COASTLINE, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')

    ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                 'admin_1_states_provinces_lines',
                                                 '50m',
                                                 linewidths = 0.75,
                                                 facecolor  = 'none',
                                                 edgecolor  = Mines_Blue))

    if (domain > 1) :
        ax1.add_feature(USCOUNTIES, 
                        linewidths=0.5,
                        edgecolor  = Mines_Blue,
                        facecolor='none')

    rain_norm = mpl.colors.BoundaryNorm(boundaries = precip_levels_hrly, 
											ncolors    = 15)
	
    filled_cm = hrly_rain_maps.isel(Time=t).plot.imshow(cmap          = precip_colormap,
											  alpha         = alpha2d,
											  ax            = ax1,  
											  interpolation = "bilinear",
												   extend        = 'max',
												   norm          = rain_norm,
												   levels        = precip_levels_hrly,
												  add_colorbar  = False)

    cb = plt.colorbar(filled_cm, 
				 label  = "Hourly Water-Equivalent Snowfall (in)",
				 shrink = colorbar_shrink, 
				 pad    = colorbar_pad,
				 ticks  = precip_levels_hrly)
    cb.outline.set_color(Mines_Blue)

    contour_plot2 = ax1.contour(lon2d, 
								lat2d, 
								hrly_rain_maps.isel(Time=t),
								transform = ccrs.PlateCarree(),
								colors    =            "cyan",
								linewidths=1,
								levels    = np.array([0.002]))    

    
    
    ##################################
    # 
    # Progress Bar
    #
    
    plot_box         = ax1.get_position()
    plot_box_x_start = plot_box.x0
    plot_box_y_start = plot_box.y0
    plot_box_x_end   = plot_box.x1
    plot_box_y_end   = plot_box.y1



    rect1 = patches.Rectangle(xy        = (0, 0),
                              width     = percent_done,
                              height    = 0.01, 
                              edgecolor = Mines_Blue, 
                              facecolor = Mines_Blue,
                              transform = ax1.transAxes)
    ax1.add_patch(rect1)
    
    #
    # Walking Clock
    #

    circle_theta  = np.deg2rad(np.arange(0,360,0.01))
    circle_radius = circle_theta * 0 + 1

    if (hour > 12) :
        hour = hour - 12

    angles_h = 2*np.pi*hour/12+2*np.pi*minute/(12*60)+2*second/(12*60*60)
    angles_m = 2*np.pi*minute/60+2*np.pi*second/(60*60)



    size_of_clock = 0.07

    x_clock = percent_done*(plot_box_x_end-plot_box_x_start) + plot_box_x_start - size_of_clock/2
    y_clock = plot_box_y_start-size_of_clock/2+0.01

    x_dow = percent_done
    y_dow = size_of_clock

    ax1.annotate(dow+"-"+local_time_zone, 
                 [x_dow,y_dow],
                 horizontalalignment="center",
                verticalalignment="center",
                xycoords='axes fraction')


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
    axins.plot(circle_theta, circle_radius,   color=Mines_Blue,  linewidth=1)

    #
    ##################################
            
    # plt.show()
    ax1.set_frame_on(False)
    ax1.set_title(valid_time + "  (" + local_time+")")

    fig.savefig(fname       = fig_dir_name + file_name,
				facecolor   =                  'white', 
				transparent =                    False)

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


    fig = plt.figure(figsize=figure_domain_size)

    fig.suptitle(model_run_label,fontsize="x-large")

    ax1 = fig.add_subplot(1,  # nrows
                          1,  # ncols 
                          1,  # index of figure you're installing
                          projection = cart_proj) # cartopy CRS Projection

    ax1.set_xlim(  west_east_range)
    ax1.set_ylim(south_north_range)
    ax1.set_frame_on(False)
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)


    ax1.add_feature(cfeature.LAKES, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.RIVERS, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.BORDERS, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.COASTLINE, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')

    ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                 'admin_1_states_provinces_lines',
                                                 '50m',
                                                 linewidths = 0.75,
                                                 facecolor  = 'none',
                                                 edgecolor  = Mines_Blue))

    if (domain > 1) :
        ax1.add_feature(USCOUNTIES, 
                        linewidths=0.5,
                        edgecolor  = Mines_Blue,
                        facecolor='none')

    ax1.set_title(valid_time + "  (" + local_time+")")
	
    rain_norm = mpl.colors.BoundaryNorm(boundaries = precip_levels_hrly, 
											ncolors    = 15)
	
    filled_cm = hrly_snowfall_maps.isel(Time=t).plot.imshow(cmap          = precip_colormap,
											  alpha         = alpha2d,
											  ax            = ax1,  
												   extend        = 'max',
												   norm          = rain_norm,
												   levels        = precip_levels_hrly,
												  add_colorbar  = False)

    cb = plt.colorbar(filled_cm, 
				 label  = "Hourly Water-Equivalent Snowfall (in)",
				 shrink = colorbar_shrink, 
				 pad    = colorbar_pad,
				 ticks  = precip_levels_hrly)
    cb.outline.set_color(Mines_Blue)

    contour_plot2 = ax1.contour(lon2d, 
								lat2d, 
								hrly_snowfall_maps.isel(Time=t),
								transform = ccrs.PlateCarree(),
								colors    =            "cyan",
								linewidths=1,
								levels    = np.array([0.002]))  

    rect1 = patches.Rectangle(xy        = (0, 0),
                              width     = percent_done,
                              height    = 0.01, 
                              edgecolor = Mines_Blue, 
                              facecolor = Mines_Blue,
                              transform = ax1.transAxes)
	
    ax1.set_title(valid_time + "  (" + local_time+")")

        
    
    
    ##################################
    # 
    # Progress Bar
    #
    
    plot_box         = ax1.get_position()
    plot_box_x_start = plot_box.x0
    plot_box_y_start = plot_box.y0
    plot_box_x_end   = plot_box.x1
    plot_box_y_end   = plot_box.y1



    rect1 = patches.Rectangle(xy        = (0, 0),
                              width     = percent_done,
                              height    = 0.01, 
                              edgecolor = Mines_Blue, 
                              facecolor = Mines_Blue,
                              transform = ax1.transAxes)
    ax1.add_patch(rect1)
    
    #
    # Walking Clock
    #
 

    circle_theta  = np.deg2rad(np.arange(0,360,0.01))
    circle_radius = circle_theta * 0 + 1

    if (hour > 12) :
        hour = hour - 12

    angles_h = 2*np.pi*hour/12+2*np.pi*minute/(12*60)+2*second/(12*60*60)
    angles_m = 2*np.pi*minute/60+2*np.pi*second/(60*60)



    size_of_clock = 0.07

    x_clock = percent_done*(plot_box_x_end-plot_box_x_start) + plot_box_x_start - size_of_clock/2
    y_clock = plot_box_y_start-size_of_clock/2+0.01

    x_dow = percent_done
    y_dow = size_of_clock

    ax1.annotate(dow+"-"+local_time_zone, 
                 [x_dow,y_dow],
                 horizontalalignment="center",
                verticalalignment="center",
                xycoords='axes fraction')


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
    axins.plot(circle_theta, circle_radius,   color=Mines_Blue,  linewidth=1)

    #
    ##################################
            
    # plt.show()
    ax1.set_frame_on(False)
    ax1.set_title(valid_time + "  (" + local_time+")")

    fig.savefig(fname       = fig_dir_name + file_name,
				facecolor   =                  'white', 
				transparent =                    False)

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


    fig = plt.figure(figsize=figure_domain_size)

    fig.suptitle(model_run_label,fontsize="x-large")

    ax1 = fig.add_subplot(1,  # nrows
                          1,  # ncols 
                          1,  # index of figure you're installing
                          projection = cart_proj) # cartopy CRS Projection

    ax1.set_xlim(  west_east_range)
    ax1.set_ylim(south_north_range)
    ax1.set_frame_on(False)
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)


    ax1.add_feature(cfeature.LAKES, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.RIVERS, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.BORDERS, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.COASTLINE, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')

    ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                 'admin_1_states_provinces_lines',
                                                 '50m',
                                                 linewidths = 0.75,
                                                 facecolor  = 'none',
                                                 edgecolor  = Mines_Blue))

    if (domain > 1) :
        ax1.add_feature(USCOUNTIES, 
                        linewidths=0.5,
                        edgecolor  = Mines_Blue,
                        facecolor='none')
	
    rain_norm = mpl.colors.BoundaryNorm(boundaries = snow_levels_full, 
										ncolors    = 15)

    filled_cm = snow_depth_map.isel(Time=t).plot.imshow(cmap              = precip_colormap,
															alpha         = alpha2d,
															ax            = ax1,  
															extend        = 'max',
															norm          = rain_norm,
															levels        = snow_levels_full,
															add_colorbar  = False)

    cb = plt.colorbar(filled_cm, 
					  label  = "Hourly Snow Depth (in)",
					  shrink = colorbar_shrink, 
					  pad    = colorbar_pad,
					  ticks  = snow_levels_full)  


    ax1.set_title(valid_time + "  (" + local_time+")")

    
    
    ##################################
    # 
    # Progress Bar
    #
    
    plot_box         = ax1.get_position()
    plot_box_x_start = plot_box.x0
    plot_box_y_start = plot_box.y0
    plot_box_x_end   = plot_box.x1
    plot_box_y_end   = plot_box.y1



    rect1 = patches.Rectangle(xy        = (0, 0),
                              width     = percent_done,
                              height    = 0.01, 
                              edgecolor = Mines_Blue, 
                              facecolor = Mines_Blue,
                              transform = ax1.transAxes)
    ax1.add_patch(rect1)
    
    #
    # Walking Clock
    #

    circle_theta  = np.deg2rad(np.arange(0,360,0.01))
    circle_radius = circle_theta * 0 + 1

    if (hour > 12) :
        hour = hour - 12

    angles_h = 2*np.pi*hour/12+2*np.pi*minute/(12*60)+2*second/(12*60*60)
    angles_m = 2*np.pi*minute/60+2*np.pi*second/(60*60)



    size_of_clock = 0.07

    x_clock = percent_done*(plot_box_x_end-plot_box_x_start) + plot_box_x_start - size_of_clock/2
    y_clock = plot_box_y_start-size_of_clock/2+0.01

    x_dow = percent_done
    y_dow = size_of_clock

    ax1.annotate(dow+"-"+local_time_zone, 
                 [x_dow,y_dow],
                 horizontalalignment="center",
                verticalalignment="center",
                xycoords='axes fraction')


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
    axins.plot(circle_theta, circle_radius,   color=Mines_Blue,  linewidth=1)

    #
    ##################################
            
    # plt.show()
    ax1.set_frame_on(False)
    ax1.set_title(valid_time + "  (" + local_time+")")

    fig.savefig(fname       = fig_dir_name + file_name,
				facecolor   =                  'white', 
				transparent =                    False)

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

    fig = plt.figure(figsize=figure_domain_size)

    fig.suptitle(model_run_label,fontsize="x-large")

    ax1 = fig.add_subplot(1,  # nrows
                          1,  # ncols 
                          1,  # index of figure you're installing
                          projection = cart_proj) # cartopy CRS Projection

    ax1.set_xlim(  west_east_range)
    ax1.set_ylim(south_north_range)
    ax1.set_frame_on(False)
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)


    ax1.add_feature(cfeature.LAKES, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.RIVERS, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.BORDERS, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.COASTLINE, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')

    ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                 'admin_1_states_provinces_lines',
                                                 '50m',
                                                 linewidths = 0.75,
                                                 facecolor  = 'none',
                                                 edgecolor  = Mines_Blue))
    if (domain > 1) :
        ax1.add_feature(USCOUNTIES, 
                        linewidths=0.5,
                        edgecolor  = Mines_Blue,
                        facecolor='none')

    filled_cm =  dbz_maps.isel(Time=t).plot.imshow(cmap          = cmap_radar,
												   alpha         = alpha2d,
												   norm          = norm_radar,
												   levels        = stormy_dbz_values,
												   add_colorbar  = False)

    cb = plt.colorbar(filled_cm, 
					  label  = "Reflectivity (dbZ)",
					  shrink = colorbar_shrink, 
					  pad    = colorbar_pad)
    cb.outline.set_color(Mines_Blue)
   
    ax1.set_title(valid_time + "  (" + local_time+")")
         
    
    
    ##################################
    # 
    # Progress Bar
    #
    
    plot_box         = ax1.get_position()
    plot_box_x_start = plot_box.x0
    plot_box_y_start = plot_box.y0
    plot_box_x_end   = plot_box.x1
    plot_box_y_end   = plot_box.y1



    rect1 = patches.Rectangle(xy        = (0, 0),
                              width     = percent_done,
                              height    = 0.01, 
                              edgecolor = Mines_Blue, 
                              facecolor = Mines_Blue,
                              transform = ax1.transAxes)
    ax1.add_patch(rect1)
    
    #
    # Walking Clock
    #

    circle_theta  = np.deg2rad(np.arange(0,360,0.01))
    circle_radius = circle_theta * 0 + 1

    if (hour > 12) :
        hour = hour - 12

    angles_h = 2*np.pi*hour/12+2*np.pi*minute/(12*60)+2*second/(12*60*60)
    angles_m = 2*np.pi*minute/60+2*np.pi*second/(60*60)



    size_of_clock = 0.07

    x_clock = percent_done*(plot_box_x_end-plot_box_x_start) + plot_box_x_start - size_of_clock/2
    y_clock = plot_box_y_start-size_of_clock/2+0.01

    x_dow = percent_done
    y_dow = size_of_clock

    ax1.annotate(dow+"-"+local_time_zone, 
                 [x_dow,y_dow],
                 horizontalalignment="center",
                verticalalignment="center",
                xycoords='axes fraction')


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
    axins.plot(circle_theta, circle_radius,   color=Mines_Blue,  linewidth=1)

    #
    ##################################
            
    # plt.show()
    ax1.set_frame_on(False)
    ax1.set_title(valid_time + "  (" + local_time+")")

    fig.savefig(fname       = fig_dir_name + file_name,
				facecolor   =                  'white', 
				transparent =                    False)

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

    fig = plt.figure(figsize=figure_domain_size)

    fig.suptitle(model_run_label,fontsize="x-large")

    ax1 = fig.add_subplot(1,  # nrows
                          1,  # ncols 
                          1,  # index of figure you're installing
                          projection = cart_proj) # cartopy CRS Projection

    ax1.set_xlim(  west_east_range)
    ax1.set_ylim(south_north_range)
    ax1.set_frame_on(False)
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)

    ax1.add_feature(cfeature.LAKES, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.RIVERS, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.BORDERS, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.COASTLINE, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')

    ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                 'admin_1_states_provinces_lines',
                                                 '50m',
                                                 linewidths = 0.75,
                                                 facecolor  = 'none',
                                                 edgecolor  = Mines_Blue))

    if (domain > 1) :
        ax1.add_feature(USCOUNTIES, 
                        linewidths=0.5,
                        edgecolor  = Mines_Blue,
                        facecolor='none')

    filled_cm = pblh_maps.isel(Time=t).plot.imshow(cmap          = mpl.cm.rainbow,
												   alpha         = alpha2d,
												   vmin          = pbl_height_levels.min(),
												   vmax          = pbl_height_levels.max(),
												   add_colorbar  = False)

    cb = plt.colorbar(filled_cm, 
				 label  = "PBL Height (m)",
				 shrink = colorbar_shrink, 
				 pad    = colorbar_pad)
    cb.outline.set_color(Mines_Blue)
   
    ax1.set_title(valid_time + "  (" + local_time+")")
         
    
    
    ##################################
    # 
    # Progress Bar
    #
    
    plot_box         = ax1.get_position()
    plot_box_x_start = plot_box.x0
    plot_box_y_start = plot_box.y0
    plot_box_x_end   = plot_box.x1
    plot_box_y_end   = plot_box.y1



    rect1 = patches.Rectangle(xy        = (0, 0),
                              width     = percent_done,
                              height    = 0.01, 
                              edgecolor = Mines_Blue, 
                              facecolor = Mines_Blue,
                              transform = ax1.transAxes)
    ax1.add_patch(rect1)
    
    #
    # Walking Clock
    #
   

    circle_theta  = np.deg2rad(np.arange(0,360,0.01))
    circle_radius = circle_theta * 0 + 1

    if (hour > 12) :
        hour = hour - 12

    angles_h = 2*np.pi*hour/12+2*np.pi*minute/(12*60)+2*second/(12*60*60)
    angles_m = 2*np.pi*minute/60+2*np.pi*second/(60*60)



    size_of_clock = 0.07

    x_clock = percent_done*(plot_box_x_end-plot_box_x_start) + plot_box_x_start - size_of_clock/2
    y_clock = plot_box_y_start-size_of_clock/2+0.01

    x_dow = percent_done
    y_dow = size_of_clock

    ax1.annotate(dow+"-"+local_time_zone, 
                 [x_dow,y_dow],
                 horizontalalignment="center",
                verticalalignment="center",
                xycoords='axes fraction')


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
    axins.plot(circle_theta, circle_radius,   color=Mines_Blue,  linewidth=1)

    #
    ##################################
            
    # plt.show()
    ax1.set_frame_on(False)
    ax1.set_title(valid_time + "  (" + local_time+")")

    fig.savefig(fname       = fig_dir_name + file_name,
				facecolor   =                  'white', 
				transparent =                    False)

    plt.close('all')
    
    #
    ####################################################
    
    
    
    

    ####################################################
    #
    # SRH
    #

    v_name       = "SRH"
    fig_dir_name = graphics_directory + "/" + v_name + "/"
    file_name    = "wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH + "_F" + str(t).zfill(2) + "_MAP_" + v_name + ".png"

    print(fig_dir_name + file_name)


    fig = plt.figure(figsize=figure_domain_size)

    fig.suptitle(model_run_label,fontsize="x-large")

    ax1 = fig.add_subplot(1,  # nrows
                          1,  # ncols 
                          1,  # index of figure you're installing
                          projection = cart_proj) # cartopy CRS Projection

    ax1.set_xlim(  west_east_range)
    ax1.set_ylim(south_north_range)
    ax1.set_frame_on(False)
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)

    ax1.add_feature(cfeature.LAKES, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.RIVERS, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.BORDERS, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.COASTLINE, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')

    ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                 'admin_1_states_provinces_lines',
                                                 '50m',
                                                 linewidths = 0.75,
                                                 facecolor  = 'none',
                                                 edgecolor  = Mines_Blue))

    if (domain > 1) :
        ax1.add_feature(USCOUNTIES, 
                        linewidths=0.5,
                        edgecolor  = Mines_Blue,
                        facecolor='none')

    filled_cm = helic_maps.isel(Time=t).plot.imshow(cmap          = mpl.cm.bwr,
												   alpha         = alpha2d,
                                                   vmin          = np.fabs(helic_maps).max(),
												   vmax          = np.fabs(helic_maps).max(),
												   add_colorbar  = False)

    cb = plt.colorbar(filled_cm, 
				 label  = "Sfc-3-km Storm-Relative helicity (m² s⁻²)",
				 shrink = colorbar_shrink, 
				 pad    = colorbar_pad)
    cb.outline.set_color(Mines_Blue)
   
    ax1.set_title(valid_time + "  (" + local_time+")")
         
    
    
    ##################################
    # 
    # Progress Bar
    #
    
    plot_box         = ax1.get_position()
    plot_box_x_start = plot_box.x0
    plot_box_y_start = plot_box.y0
    plot_box_x_end   = plot_box.x1
    plot_box_y_end   = plot_box.y1



    rect1 = patches.Rectangle(xy        = (0, 0),
                              width     = percent_done,
                              height    = 0.01, 
                              edgecolor = Mines_Blue, 
                              facecolor = Mines_Blue,
                              transform = ax1.transAxes)
    ax1.add_patch(rect1)
    
    #
    # Walking Clock
    #
   

    circle_theta  = np.deg2rad(np.arange(0,360,0.01))
    circle_radius = circle_theta * 0 + 1

    if (hour > 12) :
        hour = hour - 12

    angles_h = 2*np.pi*hour/12+2*np.pi*minute/(12*60)+2*second/(12*60*60)
    angles_m = 2*np.pi*minute/60+2*np.pi*second/(60*60)



    size_of_clock = 0.07

    x_clock = percent_done*(plot_box_x_end-plot_box_x_start) + plot_box_x_start - size_of_clock/2
    y_clock = plot_box_y_start-size_of_clock/2+0.01

    x_dow = percent_done
    y_dow = size_of_clock

    ax1.annotate(dow+"-"+local_time_zone, 
                 [x_dow,y_dow],
                 horizontalalignment="center",
                verticalalignment="center",
                xycoords='axes fraction')


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
    axins.plot(circle_theta, circle_radius,   color=Mines_Blue,  linewidth=1)

    #
    ##################################
            
    # plt.show()
    ax1.set_frame_on(False)
    ax1.set_title(valid_time + "  (" + local_time+")")

    fig.savefig(fname       = fig_dir_name + file_name,
				facecolor   =                  'white', 
				transparent =                    False)

    plt.close('all')
    
    #
    ####################################################
    
  

    ####################################################
    #
    # SRH
    #

    v_name       = "SRH"
    fig_dir_name = graphics_directory + "/" + v_name + "/"
    file_name    = "wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH + "_F" + str(t).zfill(2) + "_MAP_" + v_name + ".png"

    print(fig_dir_name + file_name)


    fig = plt.figure(figsize=figure_domain_size)

    fig.suptitle(model_run_label,fontsize="x-large")

    ax1 = fig.add_subplot(1,  # nrows
                          1,  # ncols 
                          1,  # index of figure you're installing
                          projection = cart_proj) # cartopy CRS Projection

    ax1.set_xlim(  west_east_range)
    ax1.set_ylim(south_north_range)
    ax1.set_frame_on(False)
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)

    ax1.add_feature(cfeature.LAKES, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.RIVERS, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.BORDERS, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.COASTLINE, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')

    ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                 'admin_1_states_provinces_lines',
                                                 '50m',
                                                 linewidths = 0.75,
                                                 facecolor  = 'none',
                                                 edgecolor  = Mines_Blue))

    if (domain > 1) :
        ax1.add_feature(USCOUNTIES, 
                        linewidths=0.5,
                        edgecolor  = Mines_Blue,
                        facecolor='none')

    filled_cm = helic_maps.isel(Time=t).plot.imshow(cmap          = mpl.cm.bwr,
												   alpha         = alpha2d,
                                                   vmin          = -np.fabs(helic_maps).max(),
												   vmax          = np.fabs(helic_maps).max(),
												   add_colorbar  = False)

    cb = plt.colorbar(filled_cm, 
				 label  = "Sfc-3km Storm-Relative Helicity (m² s⁻²)",
				 shrink = colorbar_shrink, 
				 pad    = colorbar_pad)
    cb.outline.set_color(Mines_Blue)
   
    ax1.set_title(valid_time + "  (" + local_time+")")
         
    
    
    ##################################
    # 
    # Progress Bar
    #
    
    plot_box         = ax1.get_position()
    plot_box_x_start = plot_box.x0
    plot_box_y_start = plot_box.y0
    plot_box_x_end   = plot_box.x1
    plot_box_y_end   = plot_box.y1



    rect1 = patches.Rectangle(xy        = (0, 0),
                              width     = percent_done,
                              height    = 0.01, 
                              edgecolor = Mines_Blue, 
                              facecolor = Mines_Blue,
                              transform = ax1.transAxes)
    ax1.add_patch(rect1)
    
    #
    # Walking Clock
    #
   

    circle_theta  = np.deg2rad(np.arange(0,360,0.01))
    circle_radius = circle_theta * 0 + 1

    if (hour > 12) :
        hour = hour - 12

    angles_h = 2*np.pi*hour/12+2*np.pi*minute/(12*60)+2*second/(12*60*60)
    angles_m = 2*np.pi*minute/60+2*np.pi*second/(60*60)



    size_of_clock = 0.07

    x_clock = percent_done*(plot_box_x_end-plot_box_x_start) + plot_box_x_start - size_of_clock/2
    y_clock = plot_box_y_start-size_of_clock/2+0.01

    x_dow = percent_done
    y_dow = size_of_clock

    ax1.annotate(dow+"-"+local_time_zone, 
                 [x_dow,y_dow],
                 horizontalalignment="center",
                verticalalignment="center",
                xycoords='axes fraction')


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
    axins.plot(circle_theta, circle_radius,   color=Mines_Blue,  linewidth=1)

    #
    ##################################
            
    # plt.show()
    ax1.set_frame_on(False)
    ax1.set_title(valid_time + "  (" + local_time+")")

    fig.savefig(fname       = fig_dir_name + file_name,
				facecolor   =                  'white', 
				transparent =                    False)

    plt.close('all')
    
    #
    ####################################################
    
    
  

    ####################################################
    #
    # UHEL
    #

    v_name       = "UHEL"
    fig_dir_name = graphics_directory + "/" + v_name + "/"
    file_name    = "wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH + "_F" + str(t).zfill(2) + "_MAP_" + v_name + ".png"

    print(fig_dir_name + file_name)


    fig = plt.figure(figsize=figure_domain_size)

    fig.suptitle(model_run_label,fontsize="x-large")

    ax1 = fig.add_subplot(1,  # nrows
                          1,  # ncols 
                          1,  # index of figure you're installing
                          projection = cart_proj) # cartopy CRS Projection

    ax1.set_xlim(  west_east_range)
    ax1.set_ylim(south_north_range)
    ax1.set_frame_on(False)
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)

    ax1.add_feature(cfeature.LAKES, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.RIVERS, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.BORDERS, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.COASTLINE, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')

    ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                 'admin_1_states_provinces_lines',
                                                 '50m',
                                                 linewidths = 0.75,
                                                 facecolor  = 'none',
                                                 edgecolor  = Mines_Blue))

    if (domain > 1) :
        ax1.add_feature(USCOUNTIES, 
                        linewidths=0.5,
                        edgecolor  = Mines_Blue,
                        facecolor='none')

    filled_cm = uphelic_maps.isel(Time=t).plot.imshow(cmap          = mpl.cm.bwr,
												   alpha         = alpha2d,
                                                   vmin          = -np.fabs(uphelic_maps).max(),
												   vmax          = np.fabs(uphelic_maps).max(),
												   add_colorbar  = False)

    cb = plt.colorbar(filled_cm, 
				 label  = "2km-5km Upward Helicity (m² s⁻²)",
				 shrink = colorbar_shrink, 
				 pad    = colorbar_pad)
    cb.outline.set_color(Mines_Blue)
   
    ax1.set_title(valid_time + "  (" + local_time+")")
         
    
    
    ##################################
    # 
    # Progress Bar
    #
    
    plot_box         = ax1.get_position()
    plot_box_x_start = plot_box.x0
    plot_box_y_start = plot_box.y0
    plot_box_x_end   = plot_box.x1
    plot_box_y_end   = plot_box.y1



    rect1 = patches.Rectangle(xy        = (0, 0),
                              width     = percent_done,
                              height    = 0.01, 
                              edgecolor = Mines_Blue, 
                              facecolor = Mines_Blue,
                              transform = ax1.transAxes)
    ax1.add_patch(rect1)
    
    #
    # Walking Clock
    #
   

    circle_theta  = np.deg2rad(np.arange(0,360,0.01))
    circle_radius = circle_theta * 0 + 1

    if (hour > 12) :
        hour = hour - 12

    angles_h = 2*np.pi*hour/12+2*np.pi*minute/(12*60)+2*second/(12*60*60)
    angles_m = 2*np.pi*minute/60+2*np.pi*second/(60*60)



    size_of_clock = 0.07

    x_clock = percent_done*(plot_box_x_end-plot_box_x_start) + plot_box_x_start - size_of_clock/2
    y_clock = plot_box_y_start-size_of_clock/2+0.01

    x_dow = percent_done
    y_dow = size_of_clock

    ax1.annotate(dow+"-"+local_time_zone, 
                 [x_dow,y_dow],
                 horizontalalignment="center",
                verticalalignment="center",
                xycoords='axes fraction')


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
    axins.plot(circle_theta, circle_radius,   color=Mines_Blue,  linewidth=1)

    #
    ##################################
            
    # plt.show()
    ax1.set_frame_on(False)
    ax1.set_title(valid_time + "  (" + local_time+")")

    fig.savefig(fname       = fig_dir_name + file_name,
				facecolor   =                  'white', 
				transparent =                    False)

    plt.close('all')
    
    #
    ####################################################


# ### PNG to Animated GIF

# In[ ]:


# for v_name in ("DBZ", "PBL", "RAIN", "SFCT", "DEWP",  "SNOWH", "WIND", "WEASD"):   
        
        
def png_to_gif_func(v_name):
    fig_dir_name = graphics_directory + "/" + v_name + "/"
    png_file_name    = "wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH + "_F??_MAP_" + v_name + ".png"
    gif_file_name    = "wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH + "_Fxx_MAP_" + v_name + ".gif"


    print("creating " + WRF_OVERALL_DIR + "./processing_"+v_name+"2_gif.sh")
    with open(WRF_OVERALL_DIR + "./processing_"+v_name+"2_gif.sh", 'w') as f:
        print("#!/bin/bash", file =  f)
        print(". ~/.bashrc", file =  f)                                  
        print("ulimit -s unlimited", file = f)
        if intel:
            print(". /opt/intel/oneapi/setvars.sh", file = f)
        print("cd " + WRF_OVERALL_DIR, file =  f) 
        print("convert -delay 50 " + fig_dir_name + png_file_name + " " + fig_dir_name + gif_file_name, file =  f) 
        print("echo MAIN:MAPS_"+v_name+"2::: We^re Outahere Like Vladimir", file =  f) 

    os.system("chmod a+x " + WRF_OVERALL_DIR + "./processing_"+v_name+"2_gif.sh")
    os.system(WRF_OVERALL_DIR + "./processing_"+v_name+"2_gif.sh > ./processing_"+v_name+"2_gif." + model_start_date_YYYY_MM_DD_HH + ".LOG 2>&1 ")
    os.system("date")
    print()
        


# ## File Organization

# In[ ]:


####################################################
####################################################
####################################################
#
# Directory Workspaces
#
intel         = True
beta_on       = 0
max_domains   = 2
chosen_domain = -1

if (socket.gethostname() == "kyrill"):
    WRF_OVERALL_DIR = "/projects/SD_Mines_WRF_REALTIME/"
else:
    if (platform.system() == "Darwin"):
        intel = False
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
    
model_start_datetime = pd.to_datetime(datetime.datetime.strptime(model_start_date_YYYY_MM_DD_HH0000, '%Y-%m-%d_%H:%M:%S')).tz_localize(tz="UTC")
print("Model Simulation Date ", model_start_datetime)

model_end_datetime  = model_start_datetime + datetime.timedelta(hours=36)
current_datetime     = datetime.datetime.now(tz=timezone.utc)
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


# ## Crack WRF Files
# 
# 
# 

# In[ ]:


####################################################
####################################################
####################################################
#
# Rotate through Available Files
#

for domain in range(1,max_domains+1):
    
    
    if (domain == 1): 
        figure_domain_size = (7,6)
    elif (domain ==2):
        figure_domain_size = (7,5.25)
    else:
        figure_domain_size = (7,6)
        
        
    print("Domain")
    
    graphics_directory = WRF_IMAGES + "/" + model_start_date_YYYY_MM_DD_HH + "/MAPS/d" +  str(domain).zfill(2) + "/"
    
    os.system("mkdir -pv " + graphics_directory + "SFCT")
    os.system("mkdir -pv " + graphics_directory + "WIND")
    
    os.system("mkdir -pv " + graphics_directory + "PBL")
    os.system("mkdir -pv " + graphics_directory + "DBZ")

    os.system("mkdir -pv " + graphics_directory + "DEWP")
    os.system("mkdir -pv " + graphics_directory + "SRH")
    os.system("mkdir -pv " + graphics_directory + "UHEL")


    os.system("mkdir -pv " + graphics_directory + "SNOWH")

    os.system("mkdir -pv " + graphics_directory + "RAIN")
    os.system("mkdir -pv " + graphics_directory + "WEASD")

    os.system("mkdir -pv " + graphics_directory + "TRAIN")
    os.system("mkdir -pv " + graphics_directory + "TSNOW")

    wrf_file  = WRF_EXE  + "./wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH0000
    print(wrf_file)
    
    ncf = nc4.Dataset(filename = wrf_file)
    ds  = xr.open_dataset(filename_or_obj = wrf_file)
    
    ##########################################################
    # 
    # Pulling at the CRS Proj-4 fields.  
    #

    cart_proj = wrf.get_cartopy(wrfin = ncf)


    lat2d = wrf.getvar(wrfin = ncf, varname="lat", timeidx=0)
    lon2d = wrf.getvar(wrfin = ncf, varname="lon", timeidx=0)

    wrf_proj = pyproj.Proj(proj  =                             'lcc',
                           lat_1 =     np.round(ds.TRUELAT1*100)/100, 
                           lat_2 =     np.round(ds.TRUELAT2*100)/100, 
                           lat_0 = np.round(ds.MOAD_CEN_LAT*100)/100,
                           lon_0 =                      ds.STAND_LON,
                           a     =    wrf.Constants.WRF_EARTH_RADIUS,
                           b     =    wrf.Constants.WRF_EARTH_RADIUS)

    LambertConformal_Projection = xr.DataArray(data   = int(0),
                                               attrs  = dict(grid_mapping_name              = "lambert_conformal_conic",
                                                             latitude_of_projection_origin  = np.round(ds.MOAD_CEN_LAT*100)/100,
                                                             longitude_of_central_meridian  = ds.STAND_LON,
                                                             standard_parallel              = np.round([ds.TRUELAT1*100,ds.TRUELAT2*100])/100,
                                                             earth_radius                   = wrf.Constants.WRF_EARTH_RADIUS,
                                                             _CoordinateTransformType       = "Projection",
                                                             _CoordinateAxisTypes           = "GeoX GeoY"))

    # Display the Proj-4 Information

    print("-- Proj4 Settings --")

    print(wrf_proj.to_proj4())

    # Display the Proj-4 Information

    eastings,northings=wrf_proj(lon2d,lat2d)

    west_east = xr.DataArray(data   = eastings.mean(axis = 0),
                             dims   = ["west_east"],
                             coords = dict(west_east = eastings.mean(axis = 0)),
                             attrs  = dict(standard_name       = "projection_x_coordinate",
                                           description         = "Eastings",
                                           units               = "m",
                                           _CoordinateAxisType = "GeoX"))

    south_north = xr.DataArray(data   = northings.mean(axis = 1),
                               dims   = ["south_north"],
                               coords = dict(south_north         = northings.mean(axis = 1)),
                               attrs  = dict(standard_name       = "projection_y_coordinate",
                                             description         = "Eastings",
                                             units               = "m",
                                             _CoordinateAxisType = "GeoY"))
    #
    ###########################################

    
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
    

    rain_maps          = rainc_maps.copy()
    rain_maps.values   = rainc_maps.values + rainnc_maps.values + rainsc_maps.values
    
    rain_maps.values         = rain_maps.values / 25.4
    rain_maps.attrs['units'] = 'in'
    
    rain_maps = rain_maps.assign_coords(coords = dict(south_north = south_north,
                                                  west_east   =   west_east))

   
    hrly_rain_maps     = rain_maps.copy()



    
    hrly_rain_maps.values[1:,:,:] = hrly_rain_maps.values[1:,:,:] - hrly_rain_maps.values[0:-1,:,:]
    
    hrly_rain_maps = hrly_rain_maps.assign_coords(coords = dict(south_north = south_north,
                                                                west_east   =   west_east))

    hrly_rain_maps = hrly_rain_maps.where(hrly_rain_maps >= 0.01)
    rain_maps      =      rain_maps.where(     rain_maps >= 0.01)

    
    #
    # DBZ
    #
    
    dbz_maps           = wrf.getvar(wrfin    =           ncf,
                                    varname  =          'mdbz',
                                    timeidx  = wrf.ALL_TIMES)
    
    dbz_maps = dbz_maps.assign_coords(coords = dict(south_north = south_north,
                                                                west_east   =   west_east))
    #
    # helicity
    #
    
    helic_maps         = wrf.getvar(wrfin    =           ncf,
                                    varname  =         'helicity',
                                    timeidx  = wrf.ALL_TIMES) 
    helic_maps   = helic_maps.assign_coords(coords = dict(south_north = south_north,
                                                      west_east   =   west_east))

    #
    # helicity upeard
    #
    
    uphelic_maps         = wrf.getvar(wrfin    =           ncf,
                                    varname  =         'updraft_helicity',
                                    timeidx  = wrf.ALL_TIMES) 
    uphelic_maps   = uphelic_maps.assign_coords(coords = dict(south_north = south_north,
                                                      west_east   =   west_east))
    #
    # T2M
    # 
    
    t2m_maps           = wrf.getvar(wrfin    =           ncf,
                                    varname  =          'T2',
                                    timeidx  = wrf.ALL_TIMES)
    
    t2m_maps.values         = (t2m_maps.values - 273.15) * 9./5. + 32.
    t2m_maps.attrs['units'] = 'degF'
    t2m_maps = t2m_maps.assign_coords(coords = dict(south_north = south_north,
                                                    west_east   =   west_east))

    #
    # TD2M
    # 
    
    td2m_maps           = wrf.getvar(wrfin    =           ncf,
                                    varname  =          'td2',
                                    timeidx  = wrf.ALL_TIMES, 
                                    units    =          'degF')
    
    td2m_maps.attrs['units'] = 'degF'
    td2m_maps = td2m_maps.assign_coords(coords = dict(south_north = south_north,
                                                      west_east   =   west_east))

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
    u10_maps.attrs['units'] = 'kt'
    v10_maps.attrs['units'] = 'kt'
    
    u10_maps   = u10_maps.assign_coords(coords = dict(south_north = south_north,
                                                      west_east   =   west_east))
    v10_maps   = v10_maps.assign_coords(coords = dict(south_north = south_north,
                                                      west_east   =   west_east))

    m10_maps, d10_maps = wrf.getvar(wrfin    =           ncf,
                                    varname  = 'wspd_wdir10',
                                    timeidx  = wrf.ALL_TIMES, 
                                    units    =          'kt') 
    m10_maps   = m10_maps.assign_coords(coords = dict(south_north = south_north,
                                                      west_east   =   west_east))
    d10_maps   = d10_maps.assign_coords(coords = dict(south_north = south_north,
                                                      west_east   =   west_east))
    #
    # MSLP
    # 
    
    mslp_maps          = wrf.getvar(wrfin    =           ncf,
                                    varname  =         'slp',
                                    timeidx  = wrf.ALL_TIMES, 
                                    units    =         'hPa')
    mslp_maps   = mslp_maps.assign_coords(coords = dict(south_north = south_north,
                                                      west_east   =   west_east))

    
    #
    # Boundary Layer Height
    # 

    pblh_maps          = wrf.getvar(wrfin    =           ncf,
                                    varname  =        'PBLH',
                                    timeidx  = wrf.ALL_TIMES)
    pblh_maps   = pblh_maps.assign_coords(coords = dict(south_north = south_north,
                                                      west_east   =   west_east))  
    #
    # Accumulated Snow Total
    # 

    snow_fall_maps     = wrf.getvar(wrfin    =           ncf,
                                    varname  =      'SNOWNC',
                                    timeidx  = wrf.ALL_TIMES) 
    
    snow_fall_maps.values         = snow_fall_maps.values / 25.4    
    snow_fall_maps.attrs['units'] = 'in'
    
    snow_fall_maps   = snow_fall_maps.assign_coords(coords = dict(south_north = south_north,
                                                      west_east   =   west_east))  

    hrly_snowfall_maps                = snow_fall_maps.copy()
    hrly_snowfall_maps.values[1:,:,:] = snow_fall_maps.values[1:,:,:] - snow_fall_maps.values[0:-1,:,:]

    hrly_snowfall_maps   = hrly_snowfall_maps.assign_coords(coords = dict(south_north = south_north,
                                                      west_east   =   west_east))  
    
    snow_depth_map     = wrf.getvar(wrfin    =           ncf,
                                    varname  =        'SNOWH',
                                    timeidx  = wrf.ALL_TIMES) 
    
    snow_depth_map.values         = snow_depth_map.values * 39.3701 
    snow_depth_map.attrs['units'] = 'in'

    snow_depth_map   = snow_depth_map.assign_coords(coords = dict(south_north = south_north,
                                                      west_east   =   west_east))  
   
    hrly_snowdepth_maps                = snow_depth_map.copy()
    hrly_snowdepth_maps.values[1:,:,:] = snow_depth_map.values[1:,:,:] - snow_depth_map.values[0:-1,:,:]
    hrly_snowdepth_maps   = hrly_snowdepth_maps.assign_coords(coords = dict(south_north = south_north,
                                                      west_east   =   west_east))  



    hrly_snowfall_maps = hrly_snowfall_maps.where(hrly_snowfall_maps >= 0.01)
    snow_fall_maps     =     snow_fall_maps.where(    snow_fall_maps >= 0.01)
    snow_depth_map     =     snow_depth_map.where(    snow_depth_map >= 0.01)

    
   
    #
    # Accumulated Snow Total
    # 

    lat2d, lon2d = wrf.latlon_coords(var = mslp_maps)
    
    
    alpha_factor = 0.05

    ny = lon2d.shape[0]
    nx = lon2d.shape[1]      
    alpha2d = np.sqrt(np.outer(np.abs(np.hanning(ny)),np.abs(np.hanning(nx))))
    alpha2d = np.where(alpha2d>alpha_factor,alpha_factor,alpha2d)
    alpha2d = alpha2d / alpha_factor

    
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

    ###############################
    #
    # Parallel Block for Printing Maps
    #
    
    n_jobs = 16
    
    Parallel(n_jobs=n_jobs)(delayed(plot_time_series_maps_func)(MAPS_T) for MAPS_T in range(nt))

    #
    ###############################



    ####################################################
    #
    # making gifs
    #
    
    animation_list = ("DBZ", "PBL", "RAIN", "SFCT", "DEWP", "SNOWH", "SRH", "UHEL", "WIND", "WEASD")
    
    n_jobs = len(animation_list)
    
    Parallel(n_jobs=n_jobs)(delayed(png_to_gif_func)(V_NAME) for V_NAME in animation_list)
        
    #
    ####################################################    
    
        
    ####################################################
    #
    # 36-hr Total Rainfall
    #
    
    valid_time = pd.to_datetime(wrf_time_steps[-1].values).tz_localize(tz="UTC").strftime("%Y-%m-%d %H %Z")
    local_time = pd.to_datetime(wrf_time_steps[-1].values).tz_localize(tz="UTC").tz_convert(tz=tz).strftime("%Y-%m-%d %H %Z")

    model_run_label    = "Model Run " + wrf_skewt_time + "; WRF Domain " + str(domain).zfill(2)


    print(valid_time, "     -- " + model_run_label)



    v_name       = "TRAIN"
    fig_dir_name = graphics_directory + "/" + v_name + "/"
    file_name    = "wrfout_d" + str(domain).zfill(2) + "_" + model_start_date_YYYY_MM_DD_HH + "_FXX_MAP_" + v_name + ".png"

    print(fig_dir_name + file_name)


    fig = plt.figure(figsize=figure_domain_size)

    fig.suptitle(model_run_label,fontsize="x-large")
    
    ax1 = fig.add_subplot(1, 1, 1, projection = cart_proj) # cartopy CRS Projection

    ax1.set_xlim(  west_east_range)
    ax1.set_ylim(south_north_range)
    ax1.set_frame_on(False)

    ax1.add_feature(cfeature.LAKES, linewidths =     0.5, edgecolor  = Mines_Blue,facecolor  =  'none')
    ax1.add_feature(cfeature.RIVERS, linewidths =     0.5, edgecolor  = Mines_Blue, facecolor  =  'none')
    ax1.add_feature(cfeature.BORDERS, linewidths =     1.0, edgecolor  = Mines_Blue, facecolor  =  'none')
    ax1.add_feature(cfeature.COASTLINE, linewidths =     1.0, edgecolor  = Mines_Blue, facecolor  =  'none')

    ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                 'admin_1_states_provinces_lines',
                                                 '50m',
                                                 linewidths = 0.75,
                                                 facecolor  = 'none',
                                                 edgecolor  = Mines_Blue))


    if (domain > 1) :
        ax1.add_feature(USCOUNTIES, 
                        linewidths=0.5,
                        edgecolor  = Mines_Blue,
                        facecolor='none')
            
    rain_norm = mpl.colors.BoundaryNorm(boundaries = precip_levels_hrly, 
											ncolors    = 15)
	
    filled_cm = rain_maps[-1,:,:].plot.imshow(cmap          = precip_colormap,
											  alpha         = alpha2d,
											  ax            = ax1,  
											  interpolation = "bilinear",
												   extend        = 'max',
												   norm          = rain_norm,
												   levels        = precip_levels_hrly,
												  add_colorbar  = False)

    cb = plt.colorbar(filled_cm, 
				 label  = "36-hr Total Precip (in)",
				 shrink = colorbar_shrink, 
				 pad    = colorbar_pad,
				 ticks  = precip_levels_hrly)
    cb.outline.set_color(Mines_Blue)


    contour_plot2 = ax1.contour(lon2d, 
								 lat2d, 
								 rain_maps[-1,:,:],
								 transform = ccrs.PlateCarree(),
								 colors    =            "cyan",
								 linewidths=1,
								 levels    = np.array([0.002]))    

    ax1.set_title(valid_time0 + " - " + valid_timef + "  (" + local_time0 + " - " + local_timef+")")

    # plt.show()
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    
    ax1.set_frame_on(False)


    fig.savefig(fig_dir_name + file_name,
                        facecolor   = 'white', 
                        transparent =   False)

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

    fig = plt.figure(figsize=figure_domain_size, facecolor="white")

    fig.suptitle(model_run_label,fontsize="x-large")

    ax1 = fig.add_subplot(1,  # nrows
                          1,  # ncols 
                          1,  # index of figure you're installing
                          projection = cart_proj) # cartopy CRS Projection

    ax1.set_xlim(  west_east_range)
    ax1.set_ylim(south_north_range)
    ax1.set_frame_on(False)

    
    ax1.add_feature(cfeature.LAKES, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.RIVERS, 
                    linewidths =     0.5, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.BORDERS, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')
    ax1.add_feature(cfeature.COASTLINE, 
                    linewidths =     1.0, 
                    edgecolor  = Mines_Blue,
                    facecolor  =  'none')

    ax1.add_feature(cfeature.NaturalEarthFeature('cultural',
                                                 'admin_1_states_provinces_lines',
                                                 '50m',
                                                 linewidths = 0.75,
                                                 facecolor  = 'none',
                                                 edgecolor  = Mines_Blue))

    if (domain > 1) :
        ax1.add_feature(USCOUNTIES, 
                        linewidths=0.5,
                        edgecolor  = Mines_Blue,
                        facecolor  =  'none')   

    rain_norm = mpl.colors.BoundaryNorm(boundaries = snow_levels_hrly, 
										ncolors    = 15)
	
    filled_cm = snow_fall_maps[-1,:,:].plot.imshow(cmap          = precip_colormap,
												   alpha         = alpha2d,
												   ax            = ax1,  
												   interpolation = "bilinear",
												   extend        = 'max',
												   norm          = rain_norm,
												   levels        = snow_levels_hrly,
												   add_colorbar  = False)

    cb = plt.colorbar(filled_cm, 
				 label  = "36-hr Total Snowfall Liquid Water Equivalent (in)",
				 shrink = colorbar_shrink, 
				 pad    = colorbar_pad,
				 ticks  = snow_levels_hrly)
    cb.outline.set_color(Mines_Blue)

    contour_plot2 = ax1.contour(lon2d, 
								 lat2d, 
								 rain_maps[-1,:,:],
								 transform = ccrs.PlateCarree(),
								 colors    =            "cyan",
								 linewidths=1,
								 levels    = np.array([0.002]))    

    ax1.set_title(valid_time0 + " - " + valid_timef + "  (" + local_time0 + " - " + local_timef+")")

    


    # plt.show()
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)

    ax1.set_frame_on(False)

    fig.savefig(fig_dir_name + file_name,
                        facecolor   = 'white', 
                        transparent =   False)

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


# In[ ]:




