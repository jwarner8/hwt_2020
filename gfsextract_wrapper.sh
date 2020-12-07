#!/bin/bash

# Script using the wgrib2 library to extract required grib messages from GFS data file by performing
# a filtered search then extract, and then running the python script to processes these grib
# messages.
# GFS atmosphere files can be downloaded from https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/

# Make sure you export the following shell environment variables (export var=/some/path);
# $wgrib2path (location of wgrib2 executable) and
# $pythonpath (path to python env)
# $outputpath (path where you want files to be built, NOTE make sure its an empty dir)

# Run as ./gfsextract_wrapper.sh /path/to/gfs/file

$wgrib2path $1 -s | grep ':TMP:' | grep 'mb:' | $wgrib2path -i $1 -grib $outputpath/extract_3D_T.grib2
$wgrib2path $1 -s | grep ':RH:' | grep 'mb:' | $wgrib2path -i $1 -grib $outputpath/extract_3D_RH.grib2
$wgrib2path $1 -s | grep ':UGRD:' | grep 'mb:' | $wgrib2path -i $1 -grib $outputpath/extract_3D_U.grib2
$wgrib2path $1 -s | grep ':VGRD:' | grep 'mb:' | $wgrib2path -i $1 -grib $outputpath/extract_3D_V.grib2
$wgrib2path $1 -s | grep ':HGT:surface' | $wgrib2path -i $1 -grib $outputpath/extract_2D_O.grib2
$wgrib2path $1 -s | grep ':PRES:surface' | $wgrib2path -i $1 -grib $outputpath/extract_2D_P.grib2
$wgrib2path $1 -s | grep ':TMP:surface' | $wgrib2path -i $1 -grib $outputpath/extract_2D_T.grib2
$wgrib2path $1 -s | grep ':LANDN:' | $wgrib2path -i $1 -grib $outputpath/extract_2D_L.grib2
$wgrib2path $1 -s | grep ':TSOIL:' | $wgrib2path -i $1 -grib $outputpath/extract_SOIL_T.grib2
$wgrib2path $1 -s | grep ':SOILW:' | $wgrib2path -i $1 -grib $outputpath/extract_SOIL_V.grib2
$wgrib2path $1 -s | grep ':SNOD:' | $wgrib2path -i $1 -grib $outputpath/extract_2D_SNOW.grib2
$wgrib2path $1 -s | grep ':ICEC:' | $wgrib2path -i $1 -grib $outputpath/extract_2D_ICE.grib2

$pythonpath conv_gfs_grib2_to_grib.py

mv $outputpath/construct.grib $outputpath/ec_grib_X.t+X
rm $outputpath/extract*
