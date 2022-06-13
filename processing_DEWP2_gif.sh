#!/bin/bash
. ~/.bashrc
ulimit -s unlimited
cd /Users/wjc/GitHub/SD_Mines_WRF_REALTIME/
convert -delay 50 /Users/wjc/GitHub/SD_Mines_WRF_REALTIME/./WEB_IMAGES//2022-05-29_00/MAPS/d02//DEWP/wrfout_d02_2022-05-29_00_F??_MAP_DEWP.png /Users/wjc/GitHub/SD_Mines_WRF_REALTIME/./WEB_IMAGES//2022-05-29_00/MAPS/d02//DEWP/wrfout_d02_2022-05-29_00_Fxx_MAP_DEWP.gif
echo MAIN:MAPS_DEWP2::: We^re Outahere Like Vladimir
