#!/bin/bash
echo on
echo
echo "Starting Plotting Script"
echo
echo "Date: `date`"
echo
echo "Getting Date Stamps"
echo
DATESTRING=`date +"%Y-%m-%d_%H%M"`
echo
echo "Moving to Working Directory"
echo
cd /home/wjc/GitHub/SD_Mines_WRF_REALTIME
echo
echo "Starting Python Script"
echo
. source ~/.bashrc ; /home/wjc/miniconda3/bin/python /home/wjc/GitHub/SD_Mines_WRF_REALTIME/plot_tslist_meteograms_update_alt.py > /home/wjc/GitHub/SD_Mines_WRF_REALTIME/wrf_metogramupdate_${DATESTRING}.log 2>&1
echo
echo "Ending Script"
echo
echo "Completed: `date`"
echo
echo off
