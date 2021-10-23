#!/bin/bash
(
echo "Date: `date`"
source /home/wjc/.bashrc
DATESTRING=`date +"%Y-%m-%d"`
cd /home/wjc/GitHub/SD_Mines_WRF_REALTIME
python /home/wjc/GitHub/SD_Mines_WRF_REALTIME/WRF_Backbone_Script.py 2>&1 /home/wjc/GitHub/SD_Mines_WRF_REALTIME/wrf_log_${DATESTRING}.log
echo "Completed: `date`"
echo
) >>${LOGFILE:=log_wrf_realtime} 2>&1 &


