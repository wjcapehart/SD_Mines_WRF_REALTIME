#!/bin/bash
echo on
echo "Date: `date`"


DATESTRING=`date +"%Y-%m-%d%H"`
cd /home/wjc/GitHub/SD_Mines_WRF_REALTIME
/home/wjc/miniconda3/bin/python /home/wjc/GitHub/SD_Mines_WRF_REALTIME/WRF_Backbone_Script.py 2>&1 /home/wjc/GitHub/SD_Mines_WRF_REALTIME/wrf_log_${DATESTRING}.log
echo "Completed: `date`"
echo
echo off

