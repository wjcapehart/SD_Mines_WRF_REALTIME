#!/bin/bash
. ~/.bashrc
ulimit -s unlimited
cd /Users/wjc/GitHub/SD_Mines_WRF_REALTIME/
convert -delay 50 /Users/wjc/GitHub/SD_Mines_WRF_REALTIME/./WEB_IMAGES//2023-01-01_12/MAPS/d02//DEWP/wrfout_d02_2023-01-01_12_F??_MAP_DEWP.png /Users/wjc/GitHub/SD_Mines_WRF_REALTIME/./WEB_IMAGES//2023-01-01_12/MAPS/d02//DEWP/wrfout_d02_2023-01-01_12_Fxx_MAP_DEWP.gif
echo MAIN:MAPS_DEWP2::: We^re Outahere Like Vladimir
