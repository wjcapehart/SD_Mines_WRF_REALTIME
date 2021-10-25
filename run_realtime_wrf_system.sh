#!/bin/bash
echo on
echo
echo "Starting WRF Realtime Execution"
echo "Date: `date`"
echo
echo
DATESTRING=`date +"%Y-%m-%d_%H%M"`
echo
echo Running BashRC
echo
echo
echo "Checking for MPIRUN and MPIEXEX"
echo
which mpirun > /home/wjc/wheres_my_mpirun.txt 2>&1
which mpiexec > /home/wjc/wheres_my_mpiexec.txt 2>&1 
echo
echo "Entering Working Directory"
echo
cd /home/wjc/GitHub/SD_Mines_WRF_REALTIME
echo
echo  "Firing Things Up!"
echo
/home/wjc/miniconda3/bin/python /home/wjc/GitHub/SD_Mines_WRF_REALTIME/WRF_Backbone_Script.py > /home/wjc/GitHub/SD_Mines_WRF_REALTIME/wrf_log_${DATESTRING}.log 2>&1
echo
echo "Ending Script"
echo
echo "Completed: `date`"
echo
echo off
