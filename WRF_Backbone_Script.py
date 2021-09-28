#!/usr/bin/env python
# coding: utf-8

# # WRF Realtime Script
# 
# ## Libraries

# In[ ]:


######################################
#
# Libraries
#

import numpy             as np
import matplotlib.pyplot as plt
import ftplib            as ftplib
import datetime          as datetime
import os                as os
import platform          as platform

#
######################################


# ## Directory Control

# In[ ]:


platform.system()


# In[ ]:


######################################
#
# Directory Workspaces
#

beta_test = 0

if (platform.system() == "Darwin"):
    WRF_OVERALL_DIR ="/Volumes/nfsdrives/ias_raid/projects/SD_Mines_WRF_REALTIME/"
else:
    WRF_OVERALL_DIR ="/projects/SD_Mines_WRF_REALTIME/"


os.chdir(WRF_OVERALL_DIR)

display( "Current Working Directory is now " + os.getcwd() )



WPS_WORK = WRF_OVERALL_DIR + "./WPS_PrepArea/"
WPS_EXE  = WRF_OVERALL_DIR + "./WRF4/WPS/"
WRF_EXE  = WRF_OVERALL_DIR + "./WRF4/WRF/test/em_real/"

NCEP_FTP_URLROOT = "ftp://ftpprd.ncep.noaa.gov/pub/data/nccf/com/nam/prod/nam."


NCEP_FTP_SERVER = "ftpprd.ncep.noaa.gov"

NCEP_FTP_GRIB_DIR_ROOT = "/pub/data/nccf/com/nam/prod/nam."


#
######################################


# ## Time Settings
# 
# ### Setting up Time Intervals

# In[ ]:


######################################
#
# Model Time Intervals
#

time_interval_between_runs            =  6 # hours

time_between_output_steps             =  1 # hours

total_sumulation_time                 = 36 # hours

time_between_boundary_condition_feeds =  3 # hours

#
######################################


# ### Timings for each run
# 
# The Realtime WRF is generated every 6 hr at best.  The model takes 3 hr to 
# 
# | Model Product Time (UTC) | Wallclock Start Time (UTC) |
# |:------------------------:|:--------------------------:|
# |        00 UTC            |        03 UTC              |
# |        06 UTC            |        09 UTC              |
# |        12 UTC            |        15 UTC              |
# |        18 UTC            |        21 UTC              |

# In[ ]:


######################################
#
# Identify Specific Run by Wall Clock Window
#

lag_hours = 3

current_datetime = datetime.datetime.utcnow()



if (not beta_test) :



    #current_datetime = datetime.datetime(year   = 2021,
    #                                     month  =    8, 
    #                                     day    =    4, 
    #                                     hour   =   22,
    #                                     minute =   15)

    current_datetime_lag3 = current_datetime - datetime.timedelta(hours=lag_hours)


    if (current_datetime.day == current_datetime_lag3.day):
        if (current_datetime_lag3.hour < 6):
            fx_hour =  0
        elif (current_datetime_lag3.hour < 12):
            fx_hour =  6
        elif (current_datetime_lag3.hour < 18):
            fx_hour = 12
        else:
            fx_hour = 18

        model_start_datetime = datetime.datetime(year  = current_datetime_lag3.year,
                                                 month = current_datetime_lag3.month, 
                                                 day   = current_datetime_lag3.day, 
                                                 hour  = fx_hour)     
    else:
        fx_hour = 18
        model_start_datetime = datetime.datetime(year  = current_datetime_lag3.year,
                                                 month = current_datetime_lag3.month, 
                                                 day   = current_datetime_lag3.day, 
                                                 hour  = fx_hour)
        
else:


    model_start_datetime = datetime.datetime(year  = 2021,
                                             month =    9,
                                             day   =    7,
                                             hour  =   12)
    
model_end_datetime = model_start_datetime + datetime.timedelta(hours=total_sumulation_time)

print("           Current Time ", current_datetime)
print("WRF Forecast Start Time ", model_start_datetime)
print("  WRF Forecast End Time ", model_end_datetime)


model_start_date_YYYYMMDDHH = model_start_datetime.strftime("%Y%m%d%H")

with open(WRF_OVERALL_DIR + "./current_run.txt", 'w') as f:
    print(model_start_date_YYYYMMDDHH, file =  f)

#
######################################


# ## WRF Preprocessing System Components (WPS)
# 
# 
# 
# ### Move to WPS Work Directory

# In[ ]:



print("Entering WPS Section")

os.chdir(WPS_WORK)

display( "Current Working Directory is now " + os.getcwd() )


# ### NCEP Model Boundary Condition Times

# In[ ]:


## Enter Working WPS Direcotry


ncep_boundary_condition_hour = np.arange(0, 
                                         total_sumulation_time + 1, 
                                         time_between_boundary_condition_feeds)

ncep_boundary_condition_hour = [str(x).zfill(2) for x in ncep_boundary_condition_hour]

print(ncep_boundary_condition_hour)

len(ncep_boundary_condition_hour)


# ### Times for First Guess NCEP Runs
# 
# 

# In[ ]:


model_start_YYYYMMDD = model_start_datetime.strftime("%Y%m%d")
model_start_HH       = model_start_datetime.strftime("%H")

ftp_directory = NCEP_FTP_URLROOT + model_start_YYYYMMDD

print(ftp_directory)

# nam.tHHz.conusnest.hiresfFF.tm00.grib2

ncep_ftp_address = [ftp_directory + "/nam.t" + model_start_HH + "z.conusnest.hiresf" + x + ".tm00.grib2" for x in ncep_boundary_condition_hour]

# NCEP File Name Template
# /pub/data/nccf/com/nam/prod/nam.20210904/nam.t00z.conusnest.hiresf03.tm00.grib2
# ncep_ftp_file = NCEP_FTP_GRIB_DIR_ROOT + model_start_YYYYMMDD 

ncep_ftp_dir = NCEP_FTP_GRIB_DIR_ROOT + model_start_YYYYMMDD 

print(ncep_ftp_dir)

ncep_ftp_file = ["./nam.t" + model_start_HH + "z.conusnest.hiresf" + x + ".tm00.grib2" for x in ncep_boundary_condition_hour]

local_ftp_file = ["./ncep_first_guess_grib_" + x + ".grib2" for x in ncep_boundary_condition_hour]


display(ncep_ftp_file)

display(local_ftp_file)


# ### FTP GRIB Files from NCEP
# 

# In[ ]:


os.system("rm -frv " + WPS_WORK + "./ungriblog.txt")
os.system("rm -frv " + WPS_WORK + "./ungrib.log")


os.system("rm -frv " + WPS_WORK + "./GRIBFILE*")
os.system("rm -frv " + WPS_WORK + "./met_em*.nc")
os.system("rm -frv " + WPS_WORK + "./FILE:????-??-??_??")

os.system("rm -frv " + WPS_WORK + "./namelist.wps")
os.system("rm -frv " + WPS_WORK + "./NAMELIST_WPS_SHARE.TXT")


if (not beta_test):
    
    os.system("rm -frv " + WPS_WORK + "./ncep_first_guess_grib_*.grib2")


    ftp = ftplib.FTP(NCEP_FTP_SERVER)
    ftp.login('Anonymous','William.Capehart@sdsmt.edu')

    for file in range(0,len(ncep_boundary_condition_hour)):
        print("Downloading..." + ncep_ftp_file[file] + " to " + local_ftp_file[file])
        ftp.retrbinary("RETR " + ncep_ftp_dir + "/" + ncep_ftp_file[file],
                       open(WPS_WORK + local_ftp_file[file], 'wb').write)

    ftp.close()
    

os.system(WPS_EXE +"./link_grib.csh " +  WPS_WORK + "./ncep_first_guess_grib_*.grib2")


# ### Create Namelist.WPS File
# 
# Creates the following template and merges with a root value.
# 
# ```
# &share
#  wrf_core = 'ARW',
#  max_dom = 3,
#  start_date = '2019-09-04_12:00:00','2019-09-04_12:00:00',
#  end_date   = '2019-09-06_00:00:00','2019-09-04_12:00:00',
#  interval_seconds = 10800
# /
# ```

# In[ ]:


print("Create Namelist.WPS File")

model_start_wpsdate = model_start_datetime.strftime("%Y-%m-%d_%H:%M:00")
model_end_wpsdate   =   model_end_datetime.strftime("%Y-%m-%d_%H:%M:00")



with open(WPS_WORK + "./NAMELIST_WPS_SHARE.TXT", 'w') as f:
    print("&share", file =  f)
    print(" wrf_core         = 'ARW',", file =  f)
    print(" max_dom          =     3,", file =  f)
    print(" start_date       = '" + model_start_wpsdate + "', '"
                                  + model_start_wpsdate + "', '"
                                  + model_start_wpsdate + "',", file =  f)
    print(" end_date         = '" + model_end_wpsdate   + "', '"
                                  + model_end_wpsdate   + "', '"
                                  + model_end_wpsdate   + "',", file =  f)
    print(" interval_seconds = "  + str(time_between_boundary_condition_feeds*3600) + ",", file =  f)
    print("/", file =  f)

os.system("cat  NAMELIST_WPS_SHARE.TXT "+ WRF_OVERALL_DIR +"./namelist_files_and_local_scripts/NAMELIST_WPS_ROOT.TXT > namelist.wps")


# ### Execute UNGRIB.EXE, 

# In[ ]:


print("Executing UnGrib.exe")
os.system("time " + WPS_EXE +"./ungrib.exe >& ./ungriblog.txt ")


# ### Execute METGRID.EXE

# In[ ]:


print("Executing MetGrid.exe")

os.system("time " + WPS_EXE +"./metgrid.exe ")


# ## WRF Section

# In[ ]:



print("Entering WRF Section")

os.chdir(WRF_EXE)

display( "Current Working Directory is now " + os.getcwd() )


# ### Cleaning WRF

# In[ ]:


print("Cleaning WRF Directory")

os.system("rm -frv " + WRF_EXE + "./met_em*.nc")
os.system("rm -frv " + WRF_EXE + "./NAMELIST_TIME_CONTROL.TXT")
os.system("rm -frv " + WRF_EXE + "./namelist.input")
os.system("rm -frv " + WRF_EXE + "./wrfbdy_*")
os.system("rm -frv " + WRF_EXE + "./wrfinput_*")
os.system("rm -frv " + WRF_EXE + "./wrfout*")


# ### Pulling Metgrid files over from WPS Section

# In[ ]:


print("Copying MetGRID Files from WPS to WRF Area")

os.system("mv -v " + WPS_WORK + "./met_em.d??.*.nc  ./")


# ### Write WRF Input File Area
# 
# ```
#  &time_control
#  run_days                            = 0,
#  run_hours                           = 36,
#  run_minutes                         = 0,
#  run_seconds                         = 0,
#  start_year                          = 2019, 2019,
#  start_month                         = 09,   09,
#  start_day                           = 04,   04,
#  start_hour                          = 12,   12,
#  end_year                            = 2019, 2019,
#  end_month                           = 09,   09,
#  end_day                             = 06,   06,
#  end_hour                            = 00,   00,
#  interval_seconds                    = 10800
#  input_from_file                     = .true.,.true.,
#  history_interval                    = 60,  60,
#  frames_per_outfile                  = 1, 1,
#  restart                             = .false.,
#  restart_interval                    = 7200,
#  io_form_history                     = 2
#  io_form_restart                     = 2
#  io_form_input                       = 2
#  io_form_boundary                    = 2
#  /
# ```

# In[ ]:


print("Create Namelist.Input File")


with open(WRF_EXE + "./NAMELIST_TIME_CONTROL.TXT", 'w') as f:
    print("&time_control", file =  f)
    print(" run_days           =  0,", file =  f)
    print(" run_hours          = "+str(total_sumulation_time)+",", file =  f)
    print(" run_minutes        =  0,", file =  f)
    print(" run_seconds        =  0,", file =  f)
    print(" start_year         = ",model_start_datetime.strftime("%Y, %Y, %Y,"), file =  f)
    print(" start_month        = ",model_start_datetime.strftime("  %m,   %m,   %m,"), file =  f)
    print(" start_day          = ",model_start_datetime.strftime("  %d,   %d,   %d,"), file =  f)
    print(" start_hour         = ",model_start_datetime.strftime("  %H,   %H,   %H,"), file =  f)
    print(" end_year           = ",  model_end_datetime.strftime("%Y, %Y, %Y,"), file =  f)
    print(" end_month          = ",  model_end_datetime.strftime("  %m,   %m,   %m,"), file =  f)
    print(" end_day            = ",  model_end_datetime.strftime("  %d,   %d,   %d,"), file =  f)
    print(" end_hour           = ",  model_end_datetime.strftime("  %H,   %H,   %H,"), file =  f)
    print(" interval_seconds   = " + str(time_between_boundary_condition_feeds*3600) + ",", file =  f)
    print(" history_interval   = " + str(time_between_output_steps*60) + ",", str(time_between_output_steps*60) + ",", str(time_between_output_steps*60) + ",", file =  f)
    print(" frames_per_outfile =  100,100,100,", file =  f)
    print(" input_from_file    = .true.,.true.,.true.,", file =  f)
    print(" restart            = .false.,", file =  f)
    print(" restart_interval   = 999999,", file =  f)
    print(" io_form_history    = 2,", file =  f)
    print(" io_form_restart    = 2,", file =  f)
    print(" io_form_input      = 2,", file =  f)
    print(" io_form_boundary   = 2,", file =  f)
    print("/", file =  f)

os.system("cat  ./NAMELIST_TIME_CONTROL.TXT "+ WRF_OVERALL_DIR +"./namelist_files_and_local_scripts/NAMELIST_WRF_ROOT.TXT > namelist.input")

os.system("cp -frv "+ WRF_OVERALL_DIR +"./namelist_files_and_local_scripts/tslist ./")


# ### Run Real

# In[ ]:


print("Executing Real")

os.system("time ./real.exe")


# ### Run WRF

# In[ ]:


print("Executing WRF")

os.system("time ./wrf.exe")


# In[ ]:




