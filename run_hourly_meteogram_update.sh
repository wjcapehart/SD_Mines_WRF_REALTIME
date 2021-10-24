#!/bin/bash
echo on
echo "Date: `date`"


echo $DATESTRING
cd /home/wjc/GitHub/SD_Mines_WRF_REALTIME
/home/wjc/miniconda3/bin/python /home/wjc/GitHub/SD_Mines_WRF_REALTIME/plot_tslist_meteograms_update.py 2>&1 /home/wjc/GitHub/SD_Mines_WRF_REALTIME/wrf_metogramupdate_${DATESTRING}.log
echo "Completed: `date`"
echo
echo off


