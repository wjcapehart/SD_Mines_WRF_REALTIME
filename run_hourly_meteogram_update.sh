#!/bin/bash
(
echo "Date: `date`"
DATESTRING=`date +"%Y-%m-%d"`
source /home/wjc/.bashrc
cd /home/wjc/GitHub/SD_Mines_WRF_REALTIME
python /home/wjc/GitHub/SD_Mines_WRF_REALTIME/plot_tslist_meteograms_update.py 2>&1 /home/wjc/GitHub/SD_Mines_WRF_REALTIME/wrf_metogramupdate_${DATESTRING}.log
echo "Completed: `date`"
echo
) >>${LOGFILE:=wrf_metogramupdate_log} 2>&1 &

