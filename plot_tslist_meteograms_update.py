#!/usr/bin/env python
# coding: utf-8

# # 
# Plot Time Series w/ Observations (To get siphcat to get the text files -- USE THIS ONE!)

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
import xarray            as xr
import pandas            as pd
import glob              as glob
import siphon.catalog    as siphcat  
import siphon.ncss       as siphncss
import seaborn           as sns
import matplotlib.pyplot as plt
import pint_xarray       as px
import pint              as pint
import matplotlib.dates  as mdates
import timezonefinder    as tzf
import pytz              as pytz
import haversine         as hs
import socket            as socket
import metpy.calc        as mpcalc
import metpy.units       as mpunits
import pathlib           as pathlib

import urllib.request
import shutil

import matplotlib.font_manager as fm
import matplotlib as mpl
import metpy.io          as mpio


from requests import HTTPError
from datetime import timezone


from metpy.units import units

import airportsdata as airpt

ureg = pint.UnitRegistry()
Q_ = ureg.Quantity


def haversine(row):
    lon1 = station_lon
    lat1 = station_lat
    lon2 = row['longitude']
    lat2 = row['latitude']
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) 
    km = 6367 * c
    return km
    
    
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


plt.rcParams.update({'text.color'      : Mines_Blue,
                     'axes.labelcolor' : Mines_Blue,
					 'axes.edgecolor'  :Mines_Blue,
					 'xtick.color'     : Mines_Blue,
					 'ytick.color'    : Mines_Blue})


#
####################################################
####################################################
####################################################


# ## File Organization

# In[3]:


####################################################
####################################################
####################################################
#
# File Organization
#



beta_on     = 0
max_domains = 2

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
METAR_DIR   = WRF_OVERALL_DIR + "./METARS/"



station_list_file = WRF_OVERALL_DIR + "namelist_files_and_local_scripts/time_series_station_files_"+str(max_domains)+"_dom_all.xlsx"

os.chdir(WRF_EXE)

#
####################################################
####################################################
####################################################


# In[4]:


####################################################
####################################################
####################################################
#
# Model Geographic Limits
#

f_geog = xr.open_dataset(filename_or_obj= WRF_ARCHIVE + "/GEOGRID_EM_FILES/geo_em.d01.nc")

lat2d_d01 = f_geog[ "XLAT_C"]
lon2d_d01 = f_geog["XLONG_C"]

geospatial_lat_min =  lat2d_d01.values.min()
geospatial_lat_max =  lat2d_d01.values.max()
geospatial_lon_min =  lon2d_d01.values.min()
geospatial_lon_max =  lon2d_d01.values.max()

#
####################################################
####################################################
####################################################


# In[ ]:





# ## Time Control

# In[5]:


with open(WRF_ARCHIVE  + "./current_complete_run/current_run.txt") as f:
    model_start_date_YYYY_MM_DD_HH = f.readlines()

model_start_date_YYYY_MM_DD_HH     = model_start_date_YYYY_MM_DD_HH[0][0:13]

model_start_date_YYYY_MM_DD_HH0000 = model_start_date_YYYY_MM_DD_HH + ":00:00"
print(model_start_date_YYYY_MM_DD_HH0000)
    
model_start_datetime = pd.to_datetime(datetime.datetime.strptime(model_start_date_YYYY_MM_DD_HH0000, '%Y-%m-%d_%H:%M:%S')).tz_localize(tz="UTC")
print("Model Simulation Date ", model_start_datetime)
    
model_end_datetime   = model_start_datetime + datetime.timedelta(hours=36)
current_datetime     = datetime.datetime.now(tz=timezone.utc)
siphon_end_datetime  = min(current_datetime, model_end_datetime)

print("         Model Start Datetime is ", model_start_datetime)
print("           Model End Datetime is ",   model_end_datetime)
print("             Current Datetime is ",     current_datetime)
print("          Siphon End Datetime is ",  siphon_end_datetime)


siphon_time_series       = pd.date_range(model_start_datetime, siphon_end_datetime,freq='h')
siphon_pulls_YYYYMMDD_HH = siphon_time_series.strftime("%Y%m%d_%H00")

print(siphon_pulls_YYYYMMDD_HH)

#
####################################################
####################################################
####################################################


# ## Read tslist excel file
# 
# 

# In[6]:


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

target_time_series_as_list = available_time_series_list["Station ID"].to_list()
print(target_time_series_as_list)

#
####################################################
####################################################
####################################################


# ## Get Station Information

# In[7]:


####################################################
####################################################
####################################################
#
# Airport Data
#

airport_database = airpt.load('ICAO')

#
####################################################
####################################################
####################################################


# ## Pull METARS from Siphon Services

# In[8]:


####################################################
####################################################
####################################################
#
# Pull METARS from UNIDATA NOAAPORT Experimental Site
#
# https://thredds-test.unidata.ucar.edu/thredds/fileServer/noaaport/text/metar/metar_20210924_0000.txt

#


print("Box = ",geospatial_lat_min, geospatial_lat_max)
print("      ",geospatial_lon_min, geospatial_lon_max)


first = True
for datehour in siphon_pulls_YYYYMMDD_HH:

    # URLs & Local Work File Names
    
    metar_url  = "https://thredds-dev.unidata.ucar.edu/thredds/fileServer/noaaport/text/metar/metar_"+datehour+".txt"
    metar_file = METAR_DIR + "./metar_"+datehour+".txt"
    
    path_to_file = pathlib.Path(metar_file)
    
    print(path_to_file, path_to_file.is_file())

    # Pull File 
    
    print("downloading "+ metar_url)
    with urllib.request.urlopen(metar_url) as response, open(metar_file, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)
            
    print("cracking "+metar_file)
    try:
        indata = mpio.metar.parse_metar_file(metar_file)
        indata = indata[(indata['latitude']  > geospatial_lat_min) & 
                        (indata['latitude']  < geospatial_lat_max) &
                        (indata['longitude'] > geospatial_lon_min) &
                        (indata['longitude'] < geospatial_lon_max) ]
        indata = indata.drop(["current_wx2_symbol",
                              "current_wx3_symbol",
                              "current_wx1",
                              "current_wx2",
                              "current_wx3",
                              "altimeter",
                              "wind_gust",
                              "visibility",
                              "cloud_coverage",
                              "current_wx1_symbol",
                              "air_pressure_at_sea_level",
                              "remarks",
                              "low_cloud_type", 
                              "low_cloud_level", 
                              "medium_cloud_type",
                              "medium_cloud_level", 
                              "high_cloud_type", 
                              "high_cloud_level",
                              "highest_cloud_type", 
                              "highest_cloud_level"], axis=1)
        if first:
            first = False
            metar_dataframe = indata
        else:
            metar_dataframe = pd.concat(objs = [metar_dataframe,indata],
                                        axis =                  "index")
    except ValueError:
        print("BALLS! Parse Error")
        error_404 = True
        pass

# Cookiecut the Radar Domain Metars and Clean Data Table





metar_dataframe = metar_dataframe.sort_values("date_time").reset_index(drop=True)




print(    "purging " + METAR_DIR + "./metar_*.txt")
os.system("rm -frv " + METAR_DIR + "./metar_*.txt")

print("Dropping Duplicate Records")

metar_dataframe = metar_dataframe.drop_duplicates()



try:
    metar_station_locs = metar_dataframe[["station_id","latitude","longitude"]].drop_duplicates()
except:
    print("Dangit: there is no metar stations to locate")

print("Storing Data in CSVs")
metar_station_locs.to_csv(WRF_OVERALL_DIR + "Metar_loc_dump.csv")
metar_dataframe.to_csv(   WRF_OVERALL_DIR + "Metar_Dump.csv")
print("Metar Extraction Complete")
#
####################################################
####################################################
####################################################


# ## Rotate through Available Files

# In[ ]:


####################################################
####################################################
####################################################
#
# Rotate through Available Files
#

file_time = model_start_datetime.strftime('%Y-%m-%d_%H')

TS_DIR    = WRF_ARCHIVE  + "./current_complete_run/STATION_TIME_SERIES/"

#
# Creating Graphics Directory
#

graphics_directory = WRF_IMAGES + "./current_complete_run/STATION_TIME_SERIES/"

print("Creating " + graphics_directory)

os.system("mkdir -pv " + graphics_directory )

#
# Start File Rotation
#


for station_row in range(len(available_time_series_list)):

    ###################################################################
    #
    # Pull Station Data 
    #



    station_id     = available_time_series_list.iloc[station_row]["Station ID"]
    grid_domain    = available_time_series_list.iloc[station_row]["Domain"]
    station_name   = available_time_series_list.iloc[station_row]["Station Name"]
    station_lat    = available_time_series_list.iloc[station_row]["Latitude"]
    station_lon    = available_time_series_list.iloc[station_row]["Longitude"]


    #    
    ###################################################################

    ###################################################################
    #
    # Pull WRF Time Series
    #
    
    netcdf_file_name = TS_DIR + "./wrfout_d"+str(grid_domain).zfill(2)+"_"+file_time+"_"+station_id+".nc"
    
    wrf_timeseries = xr.open_dataset(netcdf_file_name, 
                                     engine='netcdf4')
    
    print(netcdf_file_name)

    #
    ###################################################################

    ###################################################################
    #
    # Select Metars for Closest Station
    #

  

    ###################################################################
    #
    # Get METARS for Station 
    #
    


    #
    # Get Nearest Station to Requested Station Point
    #
    
    metar_station_locs['distance'] = metar_station_locs.apply(lambda row: haversine(row), axis=1)
    
    x=metar_station_locs[ metar_station_locs['distance'] == metar_station_locs['distance'].min() ]
    
    try:
        weather_station_name = airport_database[x.iloc[0]['station_id']]["name"] #+", "+airport_database[x['station_id'][0]]["subd"]
        print("Weatherstation from airport database ", weather_station_name)
    except KeyError:
        weather_station_name = x.iloc[0]['station_id']
        print("Weatherstation from short table ", weather_station_name)

               
    metar_data = metar_dataframe[metar_dataframe["station_id"]==x.iloc[0]['station_id']].set_index("date_time")
    

    
    
    print("lon metar/station/wrf:",metar_data.iloc[0]["longitude"],station_lon,wrf_timeseries["wrf_grid_longitude"].values)
    print("lat metar/station/wrf:",metar_data.iloc[0][ "latitude"],station_lat,wrf_timeseries["wrf_grid_latitude"].values)

   

    metar_to_sta_distance  = hs.haversine((metar_data.iloc[0][ "latitude"],
                                           metar_data.iloc[0]["longitude"]),
                                          (station_lat,
                                           station_lon))

    metar_to_wrf_distance = hs.haversine((metar_data.iloc[0][ "latitude"],
                                          metar_data.iloc[0]["longitude"]),
                                          (wrf_timeseries["wrf_grid_latitude" ].values[0],
                                           wrf_timeseries["wrf_grid_longitude"].values[0]))


    sta_to_wrf_distance   = hs.haversine((station_lat,
                                          station_lon),
                                         (wrf_timeseries["wrf_grid_latitude" ].values[0],
                                          wrf_timeseries["wrf_grid_longitude"].values[0]))    
    
    print("distance between  metar and tslist ",metar_to_sta_distance)
    print("distance between  metar and    wrf ",metar_to_wrf_distance)
    print("distance between tslist and    wrf ",sta_to_wrf_distance)

    
    
    #    
    ###################################################################    

        







    #
    ###################################################################    
    

    ###################################################################
    ###################################################################
    #
    # Create Meteogram
    #

    ###################################################################
    #
    # Time Axes
    #
                                         
    tf     = tzf.TimezoneFinder()
    tz     = tf.certain_timezone_at(lng=station_lon, lat=station_lat)
    
    tzabbr = pytz.timezone(tz)#.localize(model_start_datetime)

    wrf_times  = pd.to_datetime(wrf_timeseries["time"]).tz_localize(tz="UTC").tz_convert(tz=tz)
    
    ncss_times = pd.to_datetime(metar_data.index).tz_localize(tz="UTC").tz_convert(tz=tz)



    wrf_time_seconds  =  wrf_times.minute*60+wrf_times.second 
    on_the_hour       = np.where(wrf_time_seconds ==0)
    wrf_time_hrly     = wrf_times[on_the_hour]
    wrf_time_hrly_bar = wrf_times[on_the_hour]-datetime.timedelta(minutes=30)

    #
    ###################################################################


    ###################################################################
    #
    # Precip Prep
    #

    wrf_cum_prec      = wrf_timeseries["stratiform_precipitation_amount"].values + wrf_timeseries["convective_precipitation_amount"].values
    print("Total precip: ",np.sum(wrf_cum_prec)/25.4, " in")
    if (np.sum(wrf_cum_prec)/25.4 < 0.005):
        print("No Significant Rainfall")
        wrf_cum_prec[:] = 0.000
    
    wrf_cum_hrly_prec = wrf_cum_prec[on_the_hour]
    wrf_hrly_prec     = wrf_cum_hrly_prec.copy()

    wrf_hrly_prec[1:] = wrf_cum_hrly_prec[1:] - wrf_cum_hrly_prec[0:-1]
    
    wrf_hrly_prec     = wrf_hrly_prec     / 25.4
    wrf_cum_hrly_prec = wrf_cum_hrly_prec / 25.4
    wrf_cum_prec      = wrf_cum_prec      / 25.4

    #
    ###################################################################

    ###################################################################
    #
    # Wind Barb Prep
    #

    u_wrf = (wrf_timeseries["eastward_wind_10m"]*units("m")/units("s")).pint.to("knots")[on_the_hour]
    v_wrf = (wrf_timeseries["northward_wind_10m"]*units("m")/units("s")).pint.to("knots")[on_the_hour]

    spd_wrf = np.sqrt(wrf_timeseries[ "eastward_wind_10m"]**2 + 
                      wrf_timeseries["northward_wind_10m"]**2 )
    
    
    spd_wrf = (spd_wrf *units("m")/units("s")).pint.to("knots")                 

    obs_winddir   =  metar_data["wind_direction"].to_numpy() * units("deg")
    obs_windspeed = (metar_data["wind_speed"].to_numpy() * units("m")/units("s")).to("knots") 



    u_obs = (metar_data[ "eastward_wind"].to_numpy() * units("m")/units("s")).to("knots") 
    v_obs = (metar_data["northward_wind"].to_numpy() * units("m")/units("s")).to("knots") 

    #u_obs = xr.DataArray(data=u_obs)#, coords={"time":ncss_times})
    #v_obs = xr.DataArray(data=v_obs)#, coords={"time":ncss_times})
    #
    ###################################################################

    ###################################################################
    #
    # Plot Meteogram
    #

    fig, ax = plt.subplots(figsize = (12, 8),
                           nrows   =  2, 
                           ncols   =  2,
                           sharex  =  True)

    date_form = mdates.DateFormatter("%H %Z\n%d %b", tz=pytz.timezone(tz))
    xmajor = mdates.HourLocator(interval = 6)
    xminor = mdates.HourLocator(interval = 1)
    
    print(ncss_times)
    print(metar_data)
    
    

    #
    # Temperature and Humidity
    #
    
    ax[0,0].plot(wrf_times,
             (wrf_timeseries["air_temperature_2m"]*units("K")).pint.to("degF"),
              color = "red")
    ax[0,0].plot(wrf_times,
             (wrf_timeseries["dew_point_temperature_2m"]*units("K")).pint.to("degF"),
              color = "blue")
    ax[0,0].set_ylabel("Temperature/DewPoint (°F)")


    try:
        ax[0,0].plot(ncss_times,
                 (metar_data["air_temperature"].to_numpy()*units("degC")).to("degF"),
                 marker = "o",
                 color="magenta",
                linestyle = "None")
    except ValueError:
        print("balls: air_temperature plot error")

    try:
        ax[0,0].plot(ncss_times,
                 (metar_data["dew_point_temperature"].to_numpy()*units("degC")).to("degF"),
                 marker = "o",
                 color="cyan",
                linestyle = "None")
    except ValueError:
        print("balls: dew_point_temperature plot error")


    ax[0,1].set_title("Sta-to-WRF dist: " + str(round(metar_to_wrf_distance,1)) +" km")
    ax[0,0].set_title("Nearest Sta: "     + weather_station_name +" ("+x.iloc[0]['station_id'] +")")
        


    #
    # Total Atmos Column Water + Wind Speed
    #
    
    ax[0,1].plot(wrf_times,
            spd_wrf,
              color = "steelblue")
    ax[0,1].set_ylabel("WRF Wind Speed (kts)")
 

    ax01 = ax[0,1].twinx()
    
    ax01.set_ylim(0,1)
    ax01.set_yticks([1/3.,2/3.])
    ax01.set_yticklabels(["WRF","OBS"])

    
    ax01.barbs( wrf_time_hrly, 1/3.,  u_wrf, v_wrf, color=Mines_Blue )


    

    #try:
    ax01.barbs( ncss_times,    
                   2/3.,  
                   u_obs.astype('float64'),                  
                   v_obs.astype('float64'), 
                   color      = "blue")
    #except:
     #print("balls: wind plotting error")

   
    
    #
    # Surface Energy Budget
    #

    ax[1,0].plot(wrf_times,
                 wrf_timeseries["surface_net_downward_shortwave_flux"],
                 color = "goldenrod")
    ax[1,0].plot(wrf_times,
                 wrf_timeseries["surface_net_downward_longwave_flux"],
                 color = "magenta")
    ax[1,0].plot(wrf_times,
                 wrf_timeseries["surface_upward_sensible_heat_flux"],
                 color = "red")
    ax[1,0].plot(wrf_times,
                 wrf_timeseries["surface_upward_latent_heat_flux"],
                 color = "blue")
    ax[1,0].legend(["Solar↓","LongWave↓","Heat↑","Evap↑"],frameon=False)
    ax[1,0].set_ylabel("Surface Energy Flux (W/m²)")

    ax[1,0].axhline(y=0,color="grey", linewidth=0.5)
 
    #
    # Precipitation
    #

    ax[1,1].bar(wrf_time_hrly_bar,
                wrf_hrly_prec,
                linewidth=0,
                width=1/24, 
                color="lightgreen",
                edgecolor=None)
    ax11 = ax[1,1].twinx()
    ax11.plot(wrf_times,
              wrf_cum_prec, 
              color="darkgreen")
    ax11.set_ylabel("Cumulative Precipitation (in)")
    ax[1,1].set_ylabel("Hourly Precipitation (in)")
    ax[1,1].set_ylim(0., np.max([np.max(wrf_hrly_prec), 0.005]))
    ax11.set_ylim(   0., np.max([np.max( wrf_cum_prec), 0.005]))
    fig.suptitle(station_name+"; Model Run "+file_time+"; WRF Domain "+str(grid_domain),
                 fontsize=20)


    ax[1,0].set_xlim(model_start_datetime, model_end_datetime)
    ax[1,0].xaxis.set_major_formatter(date_form)
    ax[1,0].xaxis.set_major_locator(xmajor)
    ax[1,0].xaxis.set_minor_locator(xminor)
    ax[1,0].xaxis_date()

    ax[1,1].set_xlim(model_start_datetime, model_end_datetime)
    ax[1,1].xaxis.set_major_formatter(date_form)
    ax[1,1].xaxis.set_major_locator(xmajor)
    ax[1,1].xaxis.set_minor_locator(xminor)
    ax[1,1].xaxis_date()

    ax[0,0].set_xlim(model_start_datetime, model_end_datetime)
    ax[0,0].xaxis.set_major_formatter(date_form)
    ax[0,0].xaxis.set_major_locator(xmajor)
    ax[0,0].xaxis.set_minor_locator(xminor)
    ax[0,0].xaxis_date()

    ax[0,1].set_xlim(model_start_datetime, model_end_datetime)
    ax[0,1].xaxis.set_major_formatter(date_form)
    ax[0,1].xaxis.set_major_locator(xmajor)
    ax[0,1].xaxis.set_minor_locator(xminor)
    ax[0,1].xaxis_date()

    ax[0,0].spines["top"].set_visible(False)
    ax[1,0].spines["top"].set_visible(False)
    ax[0,1].spines["top"].set_visible(False)
    ax[1,1].spines["top"].set_visible(False)
    ax11.spines[   "top"].set_visible(False)
    ax01.spines[   "top"].set_visible(False)

    ax[0,0].spines["right"].set_visible(False)
    ax[1,0].spines["right"].set_visible(False)
    ax[0,1].spines["right"].set_visible(False)
    ax01.spines[   "right"].set_visible(False)

    product_label = "Generated "+current_datetime.strftime('%Y-%m-%d %H%M %Z')
    fig.text(0,0,product_label, fontsize="xx-small")
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)


    # plt.show()
    if (platform.system() != "Darwin"):
        fig.savefig(graphics_directory + "./wrfout_dxx_"+file_time+"_"+station_id+".png",
                        facecolor   = 'white', 
                        transparent =   False)
    else:
        plt.show()

    plt.close('all')


    #
    ###################################################################

    #
    ###################################################################
    ###################################################################
    
    print(" ")



# ## Depart 

# In[ ]:


####################################################
####################################################
####################################################
#
# End of Script
#

print("Ploting Meteogram Script complete.")

#
####################################################
####################################################
####################################################


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




