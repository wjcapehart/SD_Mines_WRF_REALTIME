#!/bin/bash
echo on
echo
echo "Starting WRF Realtime Execution"
echo "Date: `date`"
echo
echo
DATESTRING=`date +"%Y-%m-%d_%H%M"`
echo
echo
echo "Entering Working Directory"
echo
cd /home/wjc/GitHub/SD_Mines_WRF_REALTIME
echo
echo  "Firing Things Up!"
echo
.  /home/wjc/.bashrc ; /home/wjc/miniconda3/bin/python /home/wjc/GitHub/SD_Mines_WRF_REALTIME/WRF_Backbone_Script.py > /home/wjc/GitHub/SD_Mines_WRF_REALTIME/wrf_log_${DATESTRING}.log 2>&1
.  /home/wjc/.bashrc ; /home/wjc/miniconda3/bin/python /home/wjc/GitHub/SD_Mines_WRF_REALTIME/plot_tslist_meteograms_update.py > /home/wjc/GitHub/SD_Mines_WRF_REALTIME/wrf_metogramupdate_${DATESTRING}.log 2>&1
echo "Ending Script"
echo
echo "Completed: `date`"
echo
echo off
