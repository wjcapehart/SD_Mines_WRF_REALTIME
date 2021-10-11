#!/usr/bin/env python
# coding: utf-8

# # WRF Backbone Script

# ## Libraries

# In[ ]:


####################################################
####################################################
####################################################
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
####################################################
####################################################
####################################################


# ## Directory Control

# In[ ]:


####################################################
####################################################
####################################################
#
# Directory Workspaces
#

beta_on = 1

if (platform.system() == "Darwin"):
    WRF_OVERALL_DIR = "/Users/wjc/GitHub/SD_Mines_WRF_REALTIME/"
else:
    WRF_OVERALL_DIR = "/projects/SD_Mines_WRF_REALTIME/"


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


# ## Time Settings
# 
# ### Setting up Time Intervals

# In[ ]:


####################################################
####################################################
####################################################
#
# Model Time Intervals
#

time_interval_between_runs            =  6 # hours

time_between_output_steps             =  1 # hours

total_sumulation_time                 = 36 # hours

time_between_boundary_condition_feeds =  3 # hours

#
####################################################
####################################################
####################################################


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


####################################################
####################################################
####################################################
#
# Identify Specific Run by Wall Clock Window
#

lag_hours = 3

current_datetime = datetime.datetime.utcnow()

if (not beta_on) :

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
                                             day   =   28,
                                             hour  =   18)
    
model_end_datetime = model_start_datetime + datetime.timedelta(hours=total_sumulation_time)

print("           Current Time ", current_datetime)
print("WRF Forecast Start Time ", model_start_datetime)
print("  WRF Forecast End Time ", model_end_datetime)

#
# Burn Current Time to File in WRF Root Directory
#

model_start_date_YYYYMMDDHH = model_start_datetime.strftime("%Y%m%d%H")

file_time = model_start_datetime.strftime("%Y-%m-%d_%H")

with open(WRF_OVERALL_DIR + "./current_run.txt", 'w') as f:
    print(model_start_date_YYYYMMDDHH, file =  f)
    
    

#
####################################################
####################################################
####################################################


# ## WRF Preprocessing System (WPS) Workflow Area

# In[ ]:


####################################################
####################################################
####################################################
#
# WRF Preprocessing System (WPS) Workflow Area
#


# ### Move to WPS Work Directory

# In[ ]:


####################################################
####################################################
#
# Move to WPS Directory
#

print("Entering WPS Section")

os.chdir(WPS_WORK)

print( "Current Working Directory is now " + os.getcwd() )

#
####################################################
####################################################


# ### NCEP Model Boundary Condition Times

# In[ ]:


####################################################
####################################################
#
# Enter Working WPS Direcotry
#

ncep_boundary_condition_hour = np.arange(0, 
                                         total_sumulation_time + 1, 
                                         time_between_boundary_condition_feeds)

ncep_boundary_condition_hour = [str(x).zfill(2) for x in ncep_boundary_condition_hour]

print(ncep_boundary_condition_hour)

len(ncep_boundary_condition_hour)

#
####################################################
####################################################


# ### Pull Times and Create URLs for First Guess NCEP Files

# In[ ]:


####################################################
####################################################
#
# Pull Times and Create URLs for First Guess NCEP Files
#

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


print(ncep_ftp_file)

print(local_ftp_file)

#
####################################################
####################################################


# ### FTP GRIB Files from NCEP
# 

# In[ ]:


####################################################
####################################################
#
# FTP GRIB Files from NCEP
#

os.system("rm -frv " + WPS_WORK + "./ungriblog.txt")
os.system("rm -frv " + WPS_WORK + "./ungrib.log")


os.system("rm -frv " + WPS_WORK + "./GRIBFILE*")
os.system("rm -frv " + WPS_WORK + "./met_em*.nc")
os.system("rm -frv " + WPS_WORK + "./FILE:????-??-??_??")

os.system("rm -frv " + WPS_WORK + "./namelist.wps")
os.system("rm -frv " + WPS_WORK + "./NAMELIST_WPS_SHARE.TXT")


if (not beta_on):
    
    os.system("rm -frv " + WPS_WORK + "./ncep_first_guess_grib_*.grib2")


    ftp = ftplib.FTP(NCEP_FTP_SERVER)
    ftp.login('Anonymous','William.Capehart@sdsmt.edu')

    for file in range(0,len(ncep_boundary_condition_hour)):
        print("Downloading..." + ncep_ftp_file[file] + " to " + local_ftp_file[file])
        ftp.retrbinary("RETR " + ncep_ftp_dir + "/" + ncep_ftp_file[file],
                       open(WPS_WORK + local_ftp_file[file], 'wb').write)

    ftp.close()

# Link Files 

os.system(WPS_EXE +"./link_grib.csh " +  WPS_WORK + "./ncep_first_guess_grib_*.grib2")

#
####################################################
####################################################


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


####################################################
####################################################
#
# Create Namelist.WPS File
#

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

#
####################################################
####################################################


# ### Execute UNGRIB.EXE 

# In[ ]:


####################################################
####################################################
#
# Execute UNGRIB.EXE
#

print("Executing UnGrib.exe")

os.system("nohup time  " + WPS_EXE +"./ungrib.exe >& ./ungriblog.txt ")

#
####################################################
####################################################


# ### Execute METGRID.EXE

# In[ ]:


####################################################
####################################################
#
# Execute METGRID.EXE
#

print("Executing MetGrid.exe")

os.system("nohup time " + WPS_EXE +"./metgrid.exe  >& ./metgridlog.txt")

#
####################################################
####################################################


# ### End WPS Section

# In[ ]:


#
# End WPS Section
#
####################################################
####################################################
####################################################


# ## WRF Forecast Model Workflow Area

# In[ ]:


####################################################
####################################################
####################################################
#
# Entering WRF Forecast Model Workflow Area
#


# ### Enter WRF Section of Workflow and Burn Model Runtime to WRF Dir

# In[ ]:


####################################################
####################################################
#
# Enter WRF Section of Workflow and 
#   Burn Model Runtime to WRF Dir
#

print("Entering WRF Section")

os.chdir(WRF_EXE)

print( "Current Working Directory is now " + os.getcwd() )

with open(WRF_EXE + "./current_run.txt", 'w') as f:
    print(model_start_date_YYYYMMDDHH, file =  f)
    
#
####################################################
####################################################


# ### Cleaning WRF

# In[ ]:


####################################################
####################################################
#
# Cleaning WRF
#

print("Cleaning WRF Directory")

os.system("rm -frv " + WRF_EXE + "./met_em*.nc")
os.system("rm -frv " + WRF_EXE + "./NAMELIST_TIME_CONTROL.TXT")
os.system("rm -frv " + WRF_EXE + "./namelist.input")
os.system("rm -frv " + WRF_EXE + "./wrfbdy_*")
os.system("rm -frv " + WRF_EXE + "./wrfinput_*")
os.system("rm -frv " + WRF_EXE + "./wrfout*")
    
#
####################################################
####################################################


# ### Moving Metgrid files over from WPS Section

# In[ ]:


####################################################
####################################################
#
# Moving Metgrid files over from WPS Section
#

print("Copying MetGRID Files from WPS to WRF Area")

os.system("mv -v " + WPS_WORK + "./met_em.d??.*.nc  ./")


#
####################################################
####################################################


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


####################################################
####################################################
#
# Write WRF Namelist File
#

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
os.system("cp -frv "+ WRF_OVERALL_DIR +"./namelist_files_and_local_scripts/bk.sh ./")


#
####################################################
####################################################


# ### Run REAL.EXE

# In[ ]:


####################################################
####################################################
#
# Run REAL.EXE
#

print("Executing Real")

os.system("nohup time mpiexec -machinefile ~wjc/nodes.wrf -np 24 ./real.exe >& reallog.txt")

#
####################################################
####################################################


# ### Run WRF.EXE

# In[ ]:


####################################################
####################################################
#
# Run REAL.EXE
#

print("Executing WRF")

os.system("nohup time mpiexec -machinefile ~wjc/nodes.wrf -np 24 ./wrf.exe  >& wrflog.txt")

#
####################################################
####################################################


# ### End WRF Section

# In[ ]:


#
# End WRF Section
#
####################################################
####################################################
####################################################


# ## Post Processing Workflow Area

# In[ ]:


####################################################
####################################################
####################################################
#
# Post Processing Workflow Area
#


# ### Compose Final Post Processing

# In[ ]:


####################################################
####################################################
#
# TSLIST Translation to NetCDF
#


netcdf_archive_directory = WRF_ARCHIVE + "/" + file_time + "/NETCDF/"

print("Creating " + netcdf_archive_directory)

os.system("mkdir -pv " + netcdf_archive_directory )


print("Creating " + WRF_IMAGES + file_time  )

os.system("mkdir -pv " + WRF_IMAGES + file_time )




os.chdir(WRF_OVERALL_DIR)

print("creating " + WRF_OVERALL_DIR + "./wrf_post_processing.sh")
with open(WRF_OVERALL_DIR + "./wrf_post_processing.sh", 'w') as f:
    print("#!/bin/bash", file =  f)
    print("source ~/.bashrc", file =  f)
    print("cd " + WRF_OVERALL_DIR, file =  f)
    print("( python " + WRF_OVERALL_DIR+ "./tslist_to_netcdf.py >& TS2NC."+file_time+".LOG  && python " + WRF_OVERALL_DIR + "./plot_tslist_meteograms.py >& METOGRAMS."+file_time+".LOG ) & ", file =  f) 
    print("( python " + WRF_OVERALL_DIR+ "./plot_maps_d01.py >& MAPS_D01."+file_time+".LOG ) & ", file =  f) 
    print("( python " + WRF_OVERALL_DIR+ "./plot_maps_d02.py >& MAPS_D02."+file_time+".LOG ) & ", file =  f) 
    print("( python " + WRF_OVERALL_DIR+ "./plot_maps_d03.py >& MAPS_D03."+file_time+".LOG ) & ", file =  f) 
    print("( python " + WRF_OVERALL_DIR+ "./plot_skewt_d03.py >& SKEWT_D01."+file_time+".LOG ) & ", file =  f) 
    print("( python " + WRF_OVERALL_DIR+ "./plot_skewt_d03.py >& SKEWT_D02."+file_time+".LOG ) & ", file =  f) 
    print("( python " + WRF_OVERALL_DIR+ "./plot_skewt_d03.py >& SKEWT_D03."+file_time+".LOG ) & ", file =  f) 
    print("( " + WRF_OVERALL_DIR + "namelist_files_and_local_scripts/run_unipost_frames_SDMines >& UUP."+file_time+".LOG ) & ", file =  f) 
    print("wait", file =  f) 
    print("echo We're Outahere Like Vladimir", file =  f) 
    
# os.system("rm -frv " + WRF_OVERALL_DIR + "./MAPS_D??.LOG ")
# os.system("rm -frv " + WRF_OVERALL_DIR + "./SKEWT_D??01.LOG ")
# os.system("rm -frv " + WRF_OVERALL_DIR + "./UUP.LOG")

# os.system("rm -frv " + WRF_OVERALL_DIR + "./METOGRAMS.LOG")
# os.system("rm -frv " + WRF_OVERALL_DIR + "./TS2NC.LOG")


os.system("chmod a+x " + WRF_OVERALL_DIR + "./wrf_post_processing.sh")
os.system(WRF_OVERALL_DIR + "./wrf_post_processing.sh")


#
####################################################
####################################################


# ### Push NetCDF Files to Final Destination on Archive Drive

# In[ ]:


####################################################
####################################################
#
# Graphics Production for Webpage
#


for domain in range(1, 3+1):
    print("cp -v " + WRF_EXE + "wrfout_d" + str(domain).zfill(2) + "_" + file_time + ":00:00  " + netcdf_archive_directory + "wrfout_d" + str(domain).zfill(2) + "_" + file_time + ".nc")
    os.system("cp -v " + WRF_EXE + "wrfout_d" + str(domain).zfill(2) + "_" + file_time + ":00:00  " + netcdf_archive_directory + "wrfout_d" + str(domain).zfill(2) + "_" + file_time + ".nc")

#
# Lock in Final 
#

with open(WRF_OVERALL_DIR + "/current_complete_run.txt", 'w') as f:
    print(file_time, file =  f)


with open(WRF_IMAGES + file_time + "/current_run.txt", 'w') as f:
    print(model_start_date_YYYYMMDDHH, file =  f)

with open(WRF_ARCHIVE + file_time + "/current_run.txt", 'w') as f:
    print(model_start_date_YYYYMMDDHH, file =  f)

os.chdir(WRF_IMAGES)

os.system("rm -fv " + WRF_IMAGES + "current_complete_run")
os.system("ln -sv " + WRF_IMAGES + file_time  + " " + WRF_IMAGES + "current_complete_run")

os.chdir(WRF_ARCHIVE)

          
          
os.system("rm -fv " + WRF_ARCHIVE + "current_complete_run")
os.system("ln -sv " + WRF_ARCHIVE + file_time  + " " + WRF_ARCHIVE + "current_complete_run")



#
####################################################
####################################################


# ### End Post Processing Workflow Area

# In[ ]:


#
# End Post Processing Workflow Area
#
####################################################
####################################################
####################################################


# ## Clean Up.

# In[ ]:


####################################################
####################################################
####################################################
#
# Clean-up
#


# ### End Cleanup

# In[ ]:


#
# End Post Processing Workflow Area
#
####################################################
####################################################
####################################################


# ## End of Code
