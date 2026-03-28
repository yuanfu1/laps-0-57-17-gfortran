#! /bin/bash

export LAPS_SRC_ROOT=`pwd`
export LAPSINSTALLROOT=/home/ywLaps/LAPS/LAPS

export CPP_INCLUDE_PATH=/home/ywLaps/software/grib2/include
./configure --prefix=$LAPSINSTALLROOT --netcdf=/home/ywLaps/software/netcdf4 --cc=icc --fc=ifort

make >& log.make
make install
#make install_lapsplot

