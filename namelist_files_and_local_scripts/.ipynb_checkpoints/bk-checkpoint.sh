#!/bin/sh
#
# @(#)$Id: bk.sh,v 1.9 2008/06/25 16:43:25 jleffler Exp $"
#
# Run process in background
# Immune from logoffs -- output to file log
# Ussage Example:
#
# bk mpirun -np 8 solverName -parallel
#
source ~/.bashrc
(
echo "Date: `date`"
echo "Command: $*"
nohup "$@"
echo "Completed: `date`"
echo
) >>${LOGFILE:=log} 2>&1 &


