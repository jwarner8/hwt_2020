"""
   Python script that works with the wgrib2 shell wrapper and reads in the extracted GFS grib
   messages, tweaks them using the gribapi, and saves them. It also generates required fields for
   the UM start dump such as specific humidity, derived from relative humidity.
   The file produced is compatible with the ec_recon tasks used to generate the model start dump.
   TODO: Currently mapping to ec soil levels. How suitable is this?
   TODO: Can we optimise/test the conversion of RH to Specific Humidity?
"""

__author__ = 'James Warner'
__email__  = 'james.warner@metoffice.gov.uk'

import iris
import iris_grib
import gribapi
import numpy as np
import sys
import shutil
import os


# grib message mappings for ec data, using table 128
# (https://apps.ecmwf.int/codes/grib/param-db?&filter=grib1&table=128)
ec_grib_centre = {"TA":[("centre",98), ("grib2code",0,0,0), ("edition",1), ("shortnameECMF","t"),
                       ("indicatorOfParameter",130)],
                  "UW":[("centre",98), ("grib2code",0,2,2), ("edition",1), ("shortnameECMF","u"),
                        ("indicatorOfParameter",131)],
                  "VW":[("centre",98), ("grib2code",0,2,3), ("edition",1), ("shortnameECMF","v"),
                        ("indicatorOfParameter",132)],
                  "SH":[("centre",98), ("grib2code",0,1,0), ("edition",1), ("shortnameECMF","q"),
                        ("indicatorOfParameter",133)],
                  "ST":[("centre",98), ("grib2code",0,0,17), ("edition",1), ("shortnameECMF","skt"),
                        ("indicatorOfParameter",235)],
                  "LM":[("centre",98), ("grib2code",2,0,0), ("edition",1), ("shortnameECMF","lsm"),
                        ("indicatorOfParameter",172)],
                  "SG":[("centre",98), ("grib2code",0,3,4), ("edition",1), ("shortnameECMF","z"),
                        ("indicatorOfParameter",129)],
                  "SP":[("centre",98), ("grib2code",0,3,0), ("edition",1), ("shortnameECMF","sp"),
                        ("indicatorOfParameter",134)],
                  "T1":[("centre",98), ("grib2code",192,128,139), ("edition",1),
                        ("shortnameECMF","stl1"), ("indicatorOfParameter",139),
                        ("table2Version",128), ("topLevel",0), ("bottomLevel",10)],
                  "T2":[("centre",98), ("grib2code",192,128,170), ("edition",1),
                        ("shortnameECMF","stl2"), ("indicatorOfParameter",170),
                        ("table2Version",128), ("topLevel",10), ("bottomLevel",40)],
                  "T3":[("centre",98), ("grib2code",192,128,183), ("edition",1),
                        ("shortnameECMF","stl3"), ("indicatorOfParameter",183),
                        ("table2Version",128), ("topLevel",40), ("bottomLevel",100)],
                  "T4":[("centre",98), ("grib2code",192,128,236), ("edition",1),
                        ("shortnameECMF","stl4"), ("indicatorOfParameter",236),
                        ("table2Version",128), ("topLevel",100), ("bottomLevel",200)],
                  "V1":[("centre",98), ("grib2code",192,128,39), ("edition",1),
                        ("shortnameECMF","swvl1"), ("indicatorOfParameter",39),
                        ("table2Version",128), ("topLevel",0), ("bottomLevel",10)],
                  "V2":[("centre",98), ("grib2code",192,128,40), ("edition",1),
                        ("shortnameECMF","swvl2"), ("indicatorOfParameter",40),
                        ("table2Version",128), ("topLevel",10), ("bottomLevel",40)],
                  "V3":[("centre",98), ("grib2code",192,128,41), ("edition",1),
                        ("shortnameECMF","swvl3"), ("indicatorOfParameter",41),
                        ("table2Version",128), ("topLevel",40), ("bottomLevel",100)],
                  "V4":[("centre",98), ("grib2code",192,128,42), ("edition",1),
                        ("shortnameECMF","swvl4"),("indicatorOfParameter",42),
                        ("table2Version",128), ("topLevel",100), ("bottomLevel",200)],
                  "SD":[("centre",98), ("grib2code",0,1,11), ("edition",1), ("shortnameECMF","sd"),
                        ("indicatorOfParameter",141)],
                  "CI":[("centre",98), ("grib2code",10,2,0), ("edition",1), ("shortnameECMF","ci"),
                        ("indicatorOfParameter",31)]}


def rh2sph(rh,t,p):

    """
       Function that converts relative humidity to specific humidity. There are alternative versions of this approximation.
       Sourced from https://github.com/hendrikwout/meteo/blob/master/meteo/humidity.py
       Some verification has been performed of these equations vs known values, but needs to be explored further.
    """

    Mw=18.0160 # molecular weight of water
    Md=28.9660 # molecular weight of dry air
    R =  8.31432E3 # gas constant
    Rd = R/Md # specific gas constant for dry air
    Rv = R/Mw # specific gas constant for vapour
    Lv = 2.5e6 # heat release for condensation of water vapour [J kg-1]
    eps = Mw/Md

    def esat(T):
        """ get saturation pressure (units [Pa]) for a given air temperature (units [K]). Other ways to compute this include Bolton 1980, but they compare well."""
        TK = 273.15
        e1 = 101325.0
        logTTK = np.log10(T/TK)
        esat =  e1*10**(10.79586*(1-TK/T)-5.02808*logTTK+ 1.50474*1e-4*(1.-10**(-8.29692*(T/TK-1)))+ 0.42873*1e-3*(10**(4.76955*(1-TK/T))-1)-2.2195983)
        return esat

    def esat2(T):
        """ a simpler form for the saturation pressure (units [Pa]) for a given air temperature (units [K]), based on clausius-claperyon"""
        return 611.*np.exp(-Lv/Rv*(1./T - 1./273.16))

    def rh2mixr(RH,T,p):
        """purpose: conversion relative humidity (unitless) to mixing ratio [kg/kg]"""
        es = esat(T)

        return Mw/Md*RH*es/(p-RH*es)

    def mixr2sh(W):
        """conversion from mixing ratio (units [kg/kg]) to specific humidity (units also [kg/kg])"""

        return W/(1.+W)

    w = rh2mixr(rh,t,p)
    q = mixr2sh(w)

    return q


def tweak_grib_msg(variable,cube):

    """ Iterate through cubes and pair with grib messages (using reference table) before saving"""

    for cube, grib_message in iris_grib.save_pairs_from_cube(cube):

        gribapi.grib_set_long(grib_message, "centre", ec_grib_centre[variable][0][1])
        gribapi.grib_set_long(grib_message, "discipline", ec_grib_centre[variable][1][1])
        gribapi.grib_set_long(grib_message, "parameterCategory", ec_grib_centre[variable][1][2])
        gribapi.grib_set_long(grib_message, "parameterNumber", ec_grib_centre[variable][1][3])

        gribapi.grib_set_long(grib_message, "edition", ec_grib_centre[variable][2][1])
        gribapi.grib_set(grib_message, "shortNameECMF", ec_grib_centre[variable][3][1])
        gribapi.grib_set_long(grib_message, "indicatorOfParameter", ec_grib_centre[variable][4][1])

        if (
            variable == 'V1'
            or variable == 'V2'
            or variable == 'V3'
            or variable == 'V4'
            or variable == 'T1'
            or variable == 'T2'
            or variable == 'T3'
            or variable == 'T4'
            ):

            gribapi.grib_set_long(grib_message, "table2Version", ec_grib_centre[variable][5][1])
            gribapi.grib_set_long(grib_message, "topLevel", ec_grib_centre[variable][6][1])
            gribapi.grib_set_long(grib_message, "bottomLevel", ec_grib_centre[variable][7][1])

    yield grib_message


def proc_temperature():

    """ subsample pressure levels so all variables are on common pressure grid (remove T > 10hPa) """

    cube   = iris.load(os.environ['outputpath']+'/extract_3D_T.grib2')[0]
    cube   = cube[[1,2,3,4,5,6,8,9,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33]]

    for i in range(0,31):
        iris_grib.save_messages(tweak_grib_msg("TA",cube[i]), os.environ['outputpath']+'/construct.grib',append=True)


def proc_uwind():

    cube   = iris.load(os.environ['outputpath']+'/extract_3D_U.grib2')[0]
    for i in range(0,31):
        iris_grib.save_messages(tweak_grib_msg("UW",cube[i]), os.environ['outputpath']+'/construct.grib',append=True)


def proc_vwind():

    cube   = iris.load(os.environ['outputpath']+'/extract_3D_V.grib2')[0]
    for i in range(0,31):
        iris_grib.save_messages(tweak_grib_msg("VW",cube[i]), os.environ['outputpath']+'/construct.grib',append=True)


def proc_spech():

    """ derive q using relative humidity, temperature and pressure """

    cube   = iris.load(os.environ['outputpath']+'/extract_3D_RH.grib2')[0]
    T_3D   = iris.load(os.environ['outputpath']+'/extract_3D_T.grib2')[0]
    T_3D = T_3D[[1,2,3,4,5,6,8,9,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33]]

    qcubedata = np.zeros((cube.data.shape))

    for y in range(cube.data.shape[1]):
        for x in range(cube.data.shape[2]):

            qcubedata[:,y,x] = rh2sph(cube.data[:,y,x]/100.,T_3D.data[:,y,x],cube.coord(axis='z').points)

    qcube = cube
    qcube.data = qcubedata[:] # fill with new data
    qcube.units='kg kg-1'
    qcube.rename('specific_humidity')

    for i in range(0,31):
        iris_grib.save_messages(tweak_grib_msg("SH",qcube[i]), os.environ['outputpath']+'/construct.grib',append=True)


def proc_skint():

    cube   = iris.load(os.environ['outputpath']+'/extract_2D_T.grib2')[0]
    iris_grib.save_messages(tweak_grib_msg("ST",cube), os.environ['outputpath']+'/construct.grib',append=True)


def proc_land():

    cube   = iris.load(os.environ['outputpath']+'/extract_2D_L.grib2')[0]
    cube.units='1'
    cube.rename('land_binary_mask')

    iris_grib.save_messages(tweak_grib_msg("LM",cube), os.environ['outputpath']+'/construct.grib',append=True)


def proc_orog():

    """ convert orography to surface geopotential using constant g """

    cube   = iris.load(os.environ['outputpath']+'/extract_2D_O.grib2')[0]

    cube = cube * 9.80665

    cube.units='m2 s-2'
    cube.rename('geopotential')

    iris_grib.save_messages(tweak_grib_msg("SG",cube), os.environ['outputpath']+'/construct.grib',append=True)


def proc_surfp():

    cube   = iris.load(os.environ['outputpath']+'/extract_2D_P.grib2')[0]

    iris_grib.save_messages(tweak_grib_msg("SP",cube), os.environ['outputpath']+'/construct.grib',append=True)


def proc_soilt():

    cube = iris.load(os.environ['outputpath']+'/extract_SOIL_T.grib2')

    depth = iris.coords.AuxCoord(0.05*100., bounds=[[0*100.,0.1*100.]],standard_name='depth',long_name='depth', units='cm')
    cube[0].add_aux_coord(depth)
    depth = iris.coords.AuxCoord(0.225*100., bounds=[[0*100.,0.35*100.]],standard_name='depth',long_name='depth', units='cm')
    cube[1].add_aux_coord(depth)
    depth = iris.coords.AuxCoord(0.675*100., bounds=[[0.35*100.,1.0*100.]],standard_name='depth',long_name='depth', units='cm')
    cube[2].add_aux_coord(depth)
    depth = iris.coords.AuxCoord(2*100., bounds=[[1.0*100.,3.0*100.]],standard_name='depth',long_name='depth', units='cm')
    cube[3].add_aux_coord(depth)
    
    merged_cube = cube.merge()[0]

    iris_grib.save_messages(tweak_grib_msg("T1",merged_cube[0]), os.environ['outputpath']+'/construct.grib',append=True)
    iris_grib.save_messages(tweak_grib_msg("T2",merged_cube[1]), os.environ['outputpath']+'/construct.grib',append=True)
    iris_grib.save_messages(tweak_grib_msg("T3",merged_cube[2]), os.environ['outputpath']+'/construct.grib',append=True)
    iris_grib.save_messages(tweak_grib_msg("T4",merged_cube[3]), os.environ['outputpath']+'/construct.grib',append=True)


def proc_soilv():

    cube = iris.load(os.environ['outputpath']+'/extract_SOIL_V.grib2')

    depth = iris.coords.AuxCoord(0.05*100., bounds=[[0*100.,0.1*100.]],standard_name='depth',long_name='depth', units='cm')
    cube[0].add_aux_coord(depth)
    depth = iris.coords.AuxCoord(0.225*100., bounds=[[0*100.,0.35*100.]],standard_name='depth',long_name='depth', units='cm')
    cube[1].add_aux_coord(depth)
    depth = iris.coords.AuxCoord(0.675*100., bounds=[[0.35*100.,1.0*100.]],standard_name='depth',long_name='depth', units='cm')
    cube[2].add_aux_coord(depth)
    depth = iris.coords.AuxCoord(2*100., bounds=[[1.0*100.,3.0*100.]],standard_name='depth',long_name='depth', units='cm')
    cube[3].add_aux_coord(depth)
    
    merged_cube = cube.merge()[0]

    merged_cube.units='m3 m-3'

    iris_grib.save_messages(tweak_grib_msg("V1",merged_cube[0]*10.), os.environ['outputpath']+'/construct.grib',append=True)
    iris_grib.save_messages(tweak_grib_msg("V2",merged_cube[1]*10.), os.environ['outputpath']+'/construct.grib',append=True)
    iris_grib.save_messages(tweak_grib_msg("V3",merged_cube[2]*10.), os.environ['outputpath']+'/construct.grib',append=True)
    iris_grib.save_messages(tweak_grib_msg("V4",merged_cube[3]*10.), os.environ['outputpath']+'/construct.grib',append=True)


def proc_snow():

    cube   = iris.load(os.environ['outputpath']+'/extract_2D_SNOW.grib2')[0]

    iris_grib.save_messages(tweak_grib_msg("SD",cube), os.environ['outputpath']+'/construct.grib',append=True)


def proc_ice():

    cube   = iris.load(os.environ['outputpath']+'/extract_2D_ICE.grib2')[0]

    iris_grib.save_messages(tweak_grib_msg("CI",cube), os.environ['outputpath']+'/construct.grib',append=True)


def main():

    # Individual tasks to process each extracted variable.
    proc_temperature()
    proc_uwind()
    proc_vwind()
    proc_spech()
    proc_skint()
    proc_land()
    proc_orog()
    proc_surfp()
    proc_soilt()
    proc_soilv()
    proc_snow()
    proc_ice()


if __name__ == '__main__':

   main()
