

"""
   Script that produces some initial analysis of the corrected HWT runs - pulling the data, and
   plotting.
"""

__author__ = "James Warner"
__email__  = "james.warner@metoffice.gov.uk"

import numpy as np
from matplotlib import pyplot as plt
import iris
import os
import matplotlib.gridspec as gridspec
import datetime as dt
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import scipy as sp
from scipy.stats import pearsonr
import subprocess
import glob


def get_model_data(outpath, suite):

    """
    Create queryfile for MASS to pull soil moisture, temperature profiles, along with 1.5 surface
    temperature and dewpoint.
    """

    print('Creating filtered retrieval file')
    queryfile = open(outpath+'/qf', 'w')
    queryfile.write('begin\n')
    queryfile.write('stash=(8225,8223,3236,3250)'+'\n') # soil temperature and moisture (a), 1.5m t,d (b)
    queryfile.write('pp_file=(')

    # Iterate over months, days and initialisation times. It doesn't matter if a file doesn't
    # exist - in this case experiments start on the 28th May. The pver stash streams are stored
    # in blocks of 12h for both experiments.
    for month in ['04','05']:
        for day in range(1,32):
            for hrblock in ['000','012','024','036','048']:

                queryfile.write('"2020'+month+str(day).zfill(2)+'T0000Z_HWT_km2p2_RA2M_pvera'+\
                                  hrblock+'.pp", ')

                queryfile.write('"2020'+month+str(day).zfill(2)+'T0000Z_HWT_km2p2_RA2M_pverb'+\
                                   hrblock+'.pp", ')

    # Write latest contents to file
    queryfile.close()

    # Remove last ', ' and replace with a bracket/new line to close the queryfile
    with open(outpath+'/qf', 'rb+') as filehandle:
        filehandle.seek(-2, os.SEEK_END)
        filehandle.truncate()

    queryfile = open(outpath+'/qf', 'a')
    queryfile.write(')\n')
    queryfile.write('end')
    queryfile.close()

    # Get data from MASS
    print('Pulling data from MASS')
    subprocess.check_output(['/opt/moose-client-wrapper/bin/moo', 'select',outpath+'/qf',
                           'moose:/devfc/'+suite+'/field.pp/',
                           outpath+'/'])


def proc_data(outpath, suite):

    """
    Function to read in fields pulled from mass in outpath, subselect correct variable from cube,
    and concatenate all timeslices into one cube before writing to disk and deleting original files.
    """

    year = '2020'

    for month in ['04','05']:
        for day in range(1,32):

            try:
                print('Try proc '+suite+', M '+month+' D '+str(day).zfill(2)+' (1.5m)')

                cubes = iris.load(outpath+'/2020'+month+str(day).zfill(2)+'T0000Z_HWT_km2p2_RA2M_pvera*.pp')

                # Create filelist of original files to delete.
                filelist = glob.glob(outpath+'/2020'+month+str(day).zfill(2)+'T0000Z_HWT_km2p2_RA2M_pvera*.pp')

                # Save new concatenated cubes
                iris.save(cubes,outpath+'/'+suite+'_2020'+month+str(day).zfill(2)+'T0000Z_T_DP_1p5m.pp')

                # Remove files
                for file in filelist:
                    os.remove(file)
                print('-> Success proc Y '+year+' M '+month+' D '+str(day).zfill(2))

                print('Try proc '+suite+', M '+month+' D '+str(day).zfill(2)+' (Soil)')

                cubes = iris.load(outpath+'/2020'+month+str(day).zfill(2)+'T0000Z_HWT_km2p2_RA2M_pverb*.pp')

                # Create filelist of original files to delete.
                filelist = glob.glob(outpath+'/2020'+month+str(day).zfill(2)+'T0000Z_HWT_km2p2_RA2M_pverb*.pp')

                # Save new concatenated cubes
                iris.save(cubes,outpath+'/'+suite+'_2020'+month+str(day).zfill(2)+'T0000Z_ST_SMC.pp')

                # Remove files
                for file in filelist:
                    os.remove(file)
                print('-> Success proc Y '+year+' M '+month+' D '+str(day).zfill(2))

            except OSError:
                print('-> Cant proc '+suite+', M '+month+' D '+str(day).zfill(2))
                pass


def analysis(outpath):

    """
    This function takes the pre-processed HWT data and produces some analysis.
    Plots that I want from this:

    - Anomaly from control with lead time (all days), for all variables (domain average), or just
      plotting the control alongside the faulty and corrected experiment.

    Can run locally but stretches resources (4G RAM), best to run on SPICE using:
    sbatch --qos=normal --mem=8G --ntasks=2 --time=120
    """

    # EXP (u-bt240, u-bv107, u-bz435), days, lead time, layers (if appropriate)
    SMC = np.zeros((3, 43, 37, 4))
    ST  = np.zeros((3, 43, 37, 4))
    TS  = np.zeros((3, 43, 37))
    DP  = np.zeros((3, 43, 37))

    # Control Simulation - u-bt240 (all UM)
    index = 0
    for month in ['04','05']:
        for day in range(1,32):

            try:
                print('Loading u-bt240, M '+month+' D '+str(day).zfill(2)+' (SMC)')
                cubes = iris.load(outpath+'/u-bt240_2020'+month+str(day).zfill(2)+'T0000Z_ST_SMC.pp',
                                  iris.NameConstraint(standard_name='moisture_content_of_soil_layer'))[0]

                SMC[0, index, :, :] = np.mean(np.mean(cubes.data,axis=3),axis=2)

                print('Loading u-bt240, M '+month+' D '+str(day).zfill(2)+' (ST)')
                cubes = iris.load(outpath+'/u-bt240_2020'+month+str(day).zfill(2)+'T0000Z_ST_SMC.pp',
                                  iris.NameConstraint(standard_name='soil_temperature'))[0]

                ST[0, index, :, :] = np.mean(np.mean(cubes.data,axis=3),axis=2)

                print('Loading u-bt240, M '+month+' D '+str(day).zfill(2)+' (TS)')
                cubes = iris.load(outpath+'/u-bt240_2020'+month+str(day).zfill(2)+'T0000Z_T_DP_1p5m.pp',
                                  iris.NameConstraint(standard_name='air_temperature'))[0]

                TS[0, index, :] = np.mean(np.mean(cubes.data,axis=2),axis=1)

                print('Loading u-bt240, M '+month+' D '+str(day).zfill(2)+' (TS)')
                cubes = iris.load(outpath+'/u-bt240_2020'+month+str(day).zfill(2)+'T0000Z_T_DP_1p5m.pp',
                                  iris.NameConstraint(standard_name='dew_point_temperature'))[0]

                DP[0, index, :] = np.mean(np.mean(cubes.data,axis=2),axis=1)

                index = index+1

            except OSError:
                pass

    # Poorly Simulation - u-bv107
    index = 0
    for month in ['04','05']:
        for day in range(1,32):

            try:
                print('Loading u-bv107, M '+month+' D '+str(day).zfill(2)+' (SMC)')
                cubes = iris.load(outpath+'/u-bv107_2020'+month+str(day).zfill(2)+'T0000Z_ST_SMC.pp',
                                  iris.NameConstraint(standard_name='moisture_content_of_soil_layer'))[0]

                SMC[1, index, :, :] = np.mean(np.mean(cubes.data,axis=3),axis=2)

                print('Loading u-bv107, M '+month+' D '+str(day).zfill(2)+' (ST)')
                cubes = iris.load(outpath+'/u-bv107_2020'+month+str(day).zfill(2)+'T0000Z_ST_SMC.pp',
                                  iris.NameConstraint(standard_name='soil_temperature'))[0]

                ST[1, index, :, :] = np.mean(np.mean(cubes.data,axis=3),axis=2)

                print('Loading u-bv107, M '+month+' D '+str(day).zfill(2)+' (TS)')
                cubes = iris.load(outpath+'/u-bv107_2020'+month+str(day).zfill(2)+'T0000Z_T_DP_1p5m.pp',
                                  iris.NameConstraint(standard_name='air_temperature'))[0]

                TS[1, index, :] = np.mean(np.mean(cubes.data,axis=2),axis=1)

                print('Loading u-bv107, M '+month+' D '+str(day).zfill(2)+' (TS)')
                cubes = iris.load(outpath+'/u-bv107_2020'+month+str(day).zfill(2)+'T0000Z_T_DP_1p5m.pp',
                                  iris.NameConstraint(standard_name='dew_point_temperature'))[0]

                DP[1, index, :] = np.mean(np.mean(cubes.data,axis=2),axis=1)

                index = index+1

            except OSError:
                pass

    # Corrected Simulation - u-bz435
    index = 0
    for month in ['04','05']:
        for day in range(1,32):

            try:
                print('Loading u-bz435, M '+month+' D '+str(day).zfill(2)+' (SMC)')
                cubes = iris.load(outpath+'/u-bz435_2020'+month+str(day).zfill(2)+'T0000Z_ST_SMC.pp',
                                  iris.NameConstraint(standard_name='moisture_content_of_soil_layer'))[0]

                SMC[2, index, :, :] = np.mean(np.mean(cubes.data,axis=3),axis=2)

                print('Loading u-bz435, M '+month+' D '+str(day).zfill(2)+' (ST)')
                cubes = iris.load(outpath+'/u-bz435_2020'+month+str(day).zfill(2)+'T0000Z_ST_SMC.pp',
                                  iris.NameConstraint(standard_name='soil_temperature'))[0]

                ST[2, index, :, :] = np.mean(np.mean(cubes.data,axis=3),axis=2)

                print('Loading u-bz435, M '+month+' D '+str(day).zfill(2)+' (TS)')
                cubes = iris.load(outpath+'/u-bz435_2020'+month+str(day).zfill(2)+'T0000Z_T_DP_1p5m.pp',
                                  iris.NameConstraint(standard_name='air_temperature'))[0]

                TS[2, index, :] = np.mean(np.mean(cubes.data,axis=2),axis=1)

                print('Loading u-bz435, M '+month+' D '+str(day).zfill(2)+' (TS)')
                cubes = iris.load(outpath+'/u-bz435_2020'+month+str(day).zfill(2)+'T0000Z_T_DP_1p5m.pp',
                                  iris.NameConstraint(standard_name='dew_point_temperature'))[0]

                DP[2, index, :] = np.mean(np.mean(cubes.data,axis=2),axis=1)

                index = index+1

            except OSError:
                pass

    # Save postprocessed output
    np.save(outpath+'/SMC.npy', SMC)
    np.save(outpath+'/ST.npy', ST)
    np.save(outpath+'/TS.npy', TS)
    np.save(outpath+'/DP.npy', DP)


def plots(outpath):

    """
    Function that reads in the .npy processed files and plots them.
    """

    # Load processed files
    SMC = np.load(outpath+'/SMC.npy')
    ST  = np.load(outpath+'/ST.npy')
    TS  = np.load(outpath+'/TS.npy')
    DP  = np.load(outpath+'/DP.npy')

    # Create and plot figure.
    gs = gridspec.GridSpec(400,200)

    ax = plt.subplot(gs[10:90,10:90])
    ax.plot(np.arange(0,37),np.mean(TS[0,:,:],axis=0), color='black',marker='o', label='u-bt240 (UM->UM)')
    ax.plot(np.arange(0,37),np.mean(TS[1,:,:],axis=0), color='red',marker='o', label='u-bv107 (GFS->UM) [bug]')
    ax.plot(np.arange(0,37),np.mean(TS[2,:,:],axis=0), color='blue',marker='o', label='u-bz435 (GFS->UM) [corrected]')
    ax.set_title('1.5m temperature [K]')
    ax.set_xlim(0,36)
    ax.legend(loc='lower right')

    ax = plt.subplot(gs[10:90,110:190])
    ax.plot(np.arange(0,37),np.mean(DP[0,:,:],axis=0), color='black',marker='o', label='u-bt240 (UM->UM)')
    ax.plot(np.arange(0,37),np.mean(DP[1,:,:],axis=0), color='red',marker='o', label='u-bv107 (GFS->UM) [bug]')
    ax.plot(np.arange(0,37),np.mean(DP[2,:,:],axis=0), color='blue',marker='o', label='u-bz435 (GFS->UM) [corrected]')
    ax.set_title('1.5m dewpoint [K]')
    ax.set_xlim(0,36)

    ax = plt.subplot(gs[110:190,10:90])
    ax.plot(np.arange(0,37),np.mean(ST[0,:,:,0],axis=0), color='black',marker='o', label='u-bt240 (UM->UM)')
    ax.plot(np.arange(0,37),np.mean(ST[1,:,:,0],axis=0), color='red',marker='o', label='u-bv107 (GFS->UM) [bug]')
    ax.plot(np.arange(0,37),np.mean(ST[2,:,:,0],axis=0), color='blue',marker='o', label='u-bz435 (GFS->UM) [corrected]')
    ax.set_title('Top soil T ~5cm [K]')
    ax.set_xlim(0,36)

    ax = plt.subplot(gs[110:190,110:190])
    ax.plot(np.arange(0,37),np.mean(SMC[0,:,:,0],axis=0), color='black',marker='o', label='u-bt240 (UM->UM)')
    ax.plot(np.arange(0,37),np.mean(SMC[1,:,:,0],axis=0), color='red',marker='o', label='u-bv107 (GFS->UM) [bug]')
    ax.plot(np.arange(0,37),np.mean(SMC[2,:,:,0],axis=0), color='blue',marker='o', label='u-bz435 (GFS->UM) [corrected]')
    ax.set_title('Top soil SMC ~5cm [kgm-2]')
    ax.set_xlim(0,36)

    ax = plt.subplot(gs[210:290,10:90])
    ax.plot(np.arange(0,37),np.mean(ST[0,:,:,1],axis=0), color='black',marker='o', label='u-bt240 (UM->UM)')
    ax.plot(np.arange(0,37),np.mean(ST[1,:,:,1],axis=0), color='red',marker='o', label='u-bv107 (GFS->UM) [bug]')
    ax.plot(np.arange(0,37),np.mean(ST[2,:,:,1],axis=0), color='blue',marker='o', label='u-bz435 (GFS->UM) [corrected]')
    ax.set_title('2nd soil T ~22.5cm [K]')
    ax.set_xlim(0,36)

    ax = plt.subplot(gs[210:290,110:190])
    ax.plot(np.arange(0,37),np.mean(SMC[0,:,:,1],axis=0), color='black',marker='o', label='u-bt240 (UM->UM)')
    ax.plot(np.arange(0,37),np.mean(SMC[1,:,:,1],axis=0), color='red',marker='o', label='u-bv107 (GFS->UM) [bug]')
    ax.plot(np.arange(0,37),np.mean(SMC[2,:,:,1],axis=0), color='blue',marker='o', label='u-bz435 (GFS->UM) [corrected]')
    ax.set_title('2nd soil SMC ~22.5cm [kgm-2]')
    ax.set_xlim(0,36)

    ax = plt.subplot(gs[310:390,10:90])
    ax.plot(np.arange(0,37),np.mean(ST[0,:,:,3],axis=0), color='black',marker='o', label='u-bt240 (UM->UM)')
    ax.plot(np.arange(0,37),np.mean(ST[1,:,:,3],axis=0), color='red',marker='o', label='u-bv107 (GFS->UM) [bug]')
    ax.plot(np.arange(0,37),np.mean(ST[2,:,:,3],axis=0), color='blue',marker='o', label='u-bz435 (GFS->UM) [corrected]')
    ax.set_title('Lowest soil T ~2m [K]')
    ax.set_xlim(0,36)
    ax.set_xlabel('Lead time (hours)')

    ax = plt.subplot(gs[310:390,110:190])
    ax.plot(np.arange(0,37),np.mean(SMC[0,:1,:,3],axis=0), color='black',marker='o', label='u-bt240 (UM->UM)')
    ax.plot(np.arange(0,37),np.mean(SMC[1,:1,:,3],axis=0), color='red',marker='o', label='u-bv107 (GFS->UM) [bug]')
    ax.plot(np.arange(0,37),np.mean(SMC[2,:1,:,3],axis=0), color='blue',marker='o', label='u-bz435 (GFS->UM) [corrected]')
    ax.set_title('Lowest soil SMC ~2m [kgm-2]')
    ax.set_xlim(0,36)
    ax.set_xlabel('Lead time (hours)')

    plt.suptitle('Whole domain average, averaged over 20th April - 31st May, HWT', weight='bold',y=0.92)

    fig = plt.gcf()
    fig.set_size_inches(14,14)
    plt.savefig('hwt_soil_spinup.png',dpi=300,bbox_inches='tight')


def main():

    # Get control data, u-bt240
    get_model_data(os.environ['SCRATCH']+'/hwt_test', 'u-bt240')
    proc_data(os.environ['SCRATCH']+'/hwt_test', 'u-bt240')

    # Get poorly data, u-bv107
    get_model_data(os.environ['SCRATCH']+'/hwt_test2', 'u-bv107')
    proc_data(os.environ['SCRATCH']+'/hwt_test2', 'u-bv107')

    # Get corrected data, u-bz435
    get_model_data(os.environ['SCRATCH']+'/hwt_test3', 'u-bz435')
    proc_data(os.environ['SCRATCH']+'/hwt_test3', 'u-bz435')

    # Run with sbatch, i.e. sbatch --qos=normal --mem=8G --ntasks=2 --time=120
    analysis(os.environ['SCRATCH']+'/hwt_test')

    plots(os.environ['SCRATCH']+'/hwt_test')


if __name__ == '__main__':

    main()
