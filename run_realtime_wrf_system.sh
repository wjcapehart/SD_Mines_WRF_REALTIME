#!/bin/bash
echo on
echo
echo "Starting WRF Realtime Execution"
echo "Date: `date`"
echo
echo
DATESTRING=`date +"%Y-%m-%d%H"`
echo
echo Running BashRC
echo
. ~/.bashrc
echo
echo "Checking for MPIRUN and MPIEXEX"
echo
which mpirun
which mpiexec
echo
echo "Entering Working Directory"
echo
cd /home/wjc/GitHub/SD_Mines_WRF_REALTIME
echo
echo  "Firing Things Up!"
echo
/home/wjc/miniconda3/bin/python /home/wjc/GitHub/SD_Mines_WRF_REALTIME/WRF_Backbone_Script.py 2>&1 /home/wjc/GitHub/SD_Mines_WRF_REALTIME/wrf_log_${DATESTRING}.log
echo
echo "Ending Script"
echo
echo "Completed: `date`"
echo
echo off
