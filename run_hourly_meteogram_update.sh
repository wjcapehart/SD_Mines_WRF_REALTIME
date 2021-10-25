#!/bin/bash
echo on
echo
echo "Starting Plotting Script"
echo
echo "Date: `date`"
echo
echo "Getting Date Stamps"
echo
echo $DATESTRING
echo
echo "Moving to Working Directory"
echo
cd /home/wjc/GitHub/SD_Mines_WRF_REALTIME
echo
echo "Starting Python Script"
echo
/home/wjc/miniconda3/bin/python /home/wjc/GitHub/SD_Mines_WRF_REALTIME/plot_tslist_meteograms_update.py 2>&1 /home/wjc/GitHub/SD_Mines_WRF_REALTIME/wrf_metogramupdate_${DATESTRING}.log
echo
echo "Ending Script"
echo
echo "Completed: `date`"
echo
echo off
