#!/bin/bash
source /home/wjc/.bash_profile
cd /Users/wjc/GitHub/SD_Mines_WRF_REALTIME/
( python /Users/wjc/GitHub/SD_Mines_WRF_REALTIME/./tslist_to_netcdf.py  && python /Users/wjc/GitHub/SD_Mines_WRF_REALTIME/./plot_tslist_meteograms.py ) &
( python /Users/wjc/GitHub/SD_Mines_WRF_REALTIME/./plot_maps.py  ) &
( python /Users/wjc/GitHub/SD_Mines_WRF_REALTIME/./plot_skewt.py  ) &
( /Users/wjc/GitHub/SD_Mines_WRF_REALTIME/namelist_files_and_local_scripts/run_unipost_frames_SDMines >& UUP.LOG ) &
wait
echo We're Outahere Like Vladimir
