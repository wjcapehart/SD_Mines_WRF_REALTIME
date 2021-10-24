#!/bin/bash
echo on
echo "Date: `date`"

__conda_setup="$('/home/wjc/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/wjc/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/wjc/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/wjc/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup

conda activate base

DATESTRING=`date +"%Y-%m-%d%H"`
cd /home/wjc/GitHub/SD_Mines_WRF_REALTIME
python /home/wjc/GitHub/SD_Mines_WRF_REALTIME/WRF_Backbone_Script.py 2>&1 /home/wjc/GitHub/SD_Mines_WRF_REALTIME/wrf_log_${DATESTRING}.log
echo "Completed: `date`"
echo
echo off

