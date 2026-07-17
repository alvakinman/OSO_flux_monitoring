"""
VGOS-OP flux calibration script.
Computes flux densities from a VO fitsidi file and a calibration scan fitsidi file.

Created by Alva Kinman 2025-11-21.
Cleaned-up version created 2026-07-16.

Usage: casa -c cal_vo_kinman2026.py vo5043 [-d]
The flag -d means use default options.
The default options assume we start from a fitsidi and antab file and want to do all processing steps.

Input files:
FITSIDI: vo5043.fits. Contains raw visibilities from Oe, Ow and potentially other stations (including autocorrelations)
CALIBRATION FITSIDI: ca5043.fits. Contains raw visibilities from flux calibration scans observed by Oe and Ow after the VO session.
ANTAB: vo5043oe+ow.antab. Contains system temeratures from Oe and Ow. Typically clipped at 500 K.
CALIBRATION ANTAB: ca5043oe+ow.antab. Like above, but for the flux calibration scans.
FLAG FILE (optional): vo5043_flags.txt. Contains flag commands in CASA format, to flag e.g. off-source time or frequency bands to remove.
CALIBRATION FLAG FILE (optional): ca5043_flags.txt. Like above but for the flux calibration scans.

NOTE: MUST BE RUN WITH CASA 6.7!

NOTE: To import gain data, a script from github.com/jive-vlbi/casa-vlbi.git is needed. This directory
 will be automatically downloaded if not present.
"""

import sys, os
import numpy as np
from time import time
from datetime import datetime

#Astropy is used to compute elevation of potential fringefit caibration scans, checking that they are above 30 deg.
# The script can run without this functionality.
# If astropy is not available, switch the function altaz_coord() with the dummy version to allow the script to run.

from astropy.coordinates import EarthLocation, SkyCoord, AltAz
from astropy.time import Time
from astropy import units as u

########################################################
# Defining names of files 
########################################################

e = sys.argv[1]

default = ''
try:
    default = sys.argv[2]
    if default == '-d':
        print('Default options will be used in this run!')
except:
    pass

if not default=='-d':
    print('The default flag "-d" was not detected! Will use the options currently in the script file.')

#Input file and main measurement sets
infitsidi=e + '.fits'
infits_calsource= 'ca' + e[2:] + '.fits'
antabfile=e + 'oe+ow.antab'
antabfile_cal = 'ca' + e[2:] + 'oe+ow.antab'
flagfile = e + '_flags.txt'
flagfile_cal = 'ca' + e[2:] + '_flags.txt'
msname_full= e + '_full.ms' #including all antennas in the fits file
msname=e+'_OTT.ms' #including only Oe and Ow
msname_calsource=e + '_cal.ms'

#Calibration tables, VO measurement set
tsystab=e+'.tsys'
tsystab_smooth = e+'.tsys_smooth'
gaindata=e+'oe+ow.gain'
gaintab=e+'.gain'
accortab=e+'.accor'
fringetab=e+'.sb_fringefit'
fringetab2 = e+'.mb_fringefit'

#Calibration tables, calibrator measurement set
tsystab_cal=e+'_cal.tsys'
tsystab_smooth_cal = e+'.tsys_smooth_cal'
gaintab_cal=e+'_cal.gain'
accortab_cal=e+'_cal.accor'
fringetab_cal=e+'_cal.fringefit'
bptab_cal=e+'_cal.bandpass'

#Averaged measurement sets
ms_avg=e + '_avg.ms'
ms_avg_cal=e + '_avg_cal.ms'

#directory to save flux data and textfile to save scale factors
fluxdir = 'fluxdata_avgamp_nooutl'
scalefile = 'scaling_factors_avgamp_nooutl'

#Outlier handling (relevant for averaging scans together) -----
#define median absolute deviation
def MAD(x):
    sample_median = np.nanmedian(x)
    return np.nanmedian(np.abs(sample_median-x))

mad_num = 5 #how many MADs away from the median does a point need to be to be an outlier?
#--------------------------------------------------------------


########################################################
# OPTIONS
########################################################
#Controls which steps to do and which to skip.
#If we have run some of the steps before, we don't need to repeat them.

options={'read_tsysdata':True, #Reads tsys data from the antab file and puts it into measurement set.
        'read_gaindata':True, #Reads gain data from the antab file and creates the gaindata catalog to be used in gencal.
        'create_ms':'fitsidi', #can be 'none' (assume {e}_OTT.ms already exists), 'split' (split the OTT baselines from {e}_full.ms), or 'fitsidi' (read infitsidi file and create {e}_full.ms, then split OTT baselines).
        'prepare_calsource':True, #Run all analysis steps on the MS containing only the flux calibration sources.
        'restore_flags':True, #Only relevant if create_ms=none. Restores MS to the chosen flag version below.
        'flagversion':'original', #Sensible options are 'original' (no flags), 'rflag1' (after earlyflags), 'before_applycal'. If we plan to run all steps, we should pick 'original'!
        'flag_spw_01':True, #Flag OTT spectral windows 0 and 1, which always have very strong unwanted radiation.
        'custom_flags':True, #if the file {e}_flags.txt exists, aply the flags in it to the main MS.
        'earlyflags':True, #do edge flagging, pcal flagging and the first rflag. Should be True unless restore_flags is False.
        'remove_gaincurve':True, #SHOULD BE TRUE! Only relevant if read_data!=none. Removes the gaincurve from the MS.
        'calc_tsys':True, #Generate calibration table for tsys. If False, assume the table already exists.
        'smooth_tsys':True, #Smooth the raw Tsys values in time to make them less noisy.
        'apply_smooth_tsys': True, #Use the smoothed Tsys values instead of the raw values.
        'calc_gain':True, #Generate calibration table for gain. If False, assume the table already exists.
        'calc_accor':True, #Generate autocorrelation correction table. If False, assume the table already exists.
        'calc_fringefit':True, #Generate instrumental fringefit table ("manual phase cal"). If False, assume the table already exists.
        'calc_fringefit2':True, #Generate multi-band fringefit table. If False, assume the table already exists. 
        'fringefit2_type':'calibrators', #relevant for calculation and application of multiband fringefit. If fringefit2_type = 'all', do multi-band fringefit on all scans. If 'calibrators', do it on all fringefit calibrator scans (listed above). If 'highsnr', do multiband fringefit on the scans with a multiband SNR above the given value below. WARNING: 'highsnr' will flag all scans that do not reach the threshold!
        'fringefit_snr':50, #only relevant if 'fringefit2_type'=='highsnr'
        'do_applycal':True, #Apply the above calibration tables.
        'apply_fringefit2':True, #apply the multi-band fringefit in Applycal.
        'rflag2':True, #do the last round of rflag after applycal.
        'set_weights':True, #set statistical weights based on stddev of visibilities.
        'make_avg_ms':True, #make a final MS where freq and time is averaged over each spw and scan. THIS NEEDS TO EXIST!

        'save_stokesi':False, #Save stokes I values in addition to the XX, YY etc. fluxes.
        'save_unscaled':True #Save unscaled fluxes in addition to the scaled ones.
        }


#THESE ARE THE DEFAULT OPTIONS, DO NOT CHANGE! ###################################
default_options={'read_tsysdata':True, #Reads tsys data from the antab file and puts it into measurement set.
        'read_gaindata':True, #Reads gain data from the antab file and creates the gaindata catalog to be used in gencal.
        'create_ms':'fitsidi', #can be 'none' (assume {e}_OTT.ms already exists), 'split' (split the OTT baselines from {e}_full.ms), or 'fitsidi' (read infitsidi file and create {e}_full.ms, then split OTT baselines).
        'prepare_calsource':True, #Run all analysis steps on the MS containing only the flux calibration sources.
        'restore_flags':True, #Only relevant if create_ms=none. Restores MS to the chosen flag version below.
        'flagversion':'original', #Sensible options are 'original' (no flags), 'rflag1' (after earlyflags), 'before_applycal'. If we plan to run all steps, we should pick 'original'!
        'flag_spw_01':True, #Flag OTT spectral windows 0 and 1, which always have very strong unwanted radiation.
        'custom_flags':True, #if the file {e}_flags.txt exists, aply the flags in it to the main MS.
        'earlyflags':True, #do edge flagging, pcal flagging and the first rflag. Should be True unless restore_flags is False.
        'remove_gaincurve':True, #SHOULD BE TRUE! Only relevant if read_data!=none. Removes the gaincurve from the MS.
        'calc_tsys':True, #Generate calibration table for tsys. If False, assume the table already exists.
        'smooth_tsys':True, #Smooth the raw Tsys values in time to make them less noisy.
        'apply_smooth_tsys': True, #Use the smoothed Tsys values instead of the raw values.
        'calc_gain':True, #Generate calibration table for gain. If False, assume the table already exists.
        'calc_accor':True, #Generate autocorrelation correction table. If False, assume the table already exists.
        'calc_fringefit':True, #Generate instrumental fringefit table ("manual phase cal"). If False, assume the table already exists.
        'calc_fringefit2':True, #Generate multi-band fringefit table. If False, assume the table already exists. 
        'fringefit2_type':'calibrators', #relevant for calculation and application of multiband fringefit. If fringefit2_type = 'all', do multi-band fringefit on all scans. If 'calibrators', do it on all fringefit calibrator scans (listed above). If 'highsnr', do multiband fringefit on the scans with a multiband SNR above the given value below. WARNING: 'highsnr' will flag all scans that do not reach the threshold!
        'fringefit_snr':50, #only relevant if 'fringefit2_type'=='highsnr'
        'do_applycal':True, #Apply the above calibration tables.
        'apply_fringefit2':True, #apply the multi-band fringefit in Applycal.
        'rflag2':True, #do the last round of rflag after applycal.
        'set_weights':True, #set statistical weights based on stddev of visibilities.
        'make_avg_ms':True, #make a final MS where freq and time is averaged over each spw and scan. THIS NEEDS TO EXIST!

        'save_stokesi':False, #Save stokes I values in addition to the XX, YY etc. fluxes.
        'save_unscaled':True #Save unscaled fluxes in addition to the scaled ones.
        }

# THE ABOVE ARE DEFAULT OPTIONS, DO NOT CHANGE! ################################

if default=='-d':
    options = default_options

casaversion = casalog.version()
scriptpath = os.path.abspath(__file__)
current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
print(f'\n Analysis done with script {scriptpath}, {casaversion}, at {current_time}\n')


#Print the current options
print('OPTIONS FOR THIS RUN ARE:')

print('read_tsysdata is ', options['read_tsysdata'])
print('read_gaindata is ', options['read_gaindata'])
print('create_ms is ', options['create_ms'])
print('prepare_calsource is ',options['prepare_calsource'])
if options['create_ms']=='none':
    print('restore_flags is ', options['restore_flags'])
    if options['restore_flags']:
        print('flagversion is ', options['flagversion'])
print('remove_gaincurve is ', options['remove_gaincurve'])
print('calc_tsys is ',options['calc_tsys'])
print('smooth_tsys is ',options['smooth_tsys'])
print('apply_smooth_tsys is ',options['apply_smooth_tsys'])
print('calc_gain is ',options['calc_gain'])
print('calc_accor is ',options['calc_accor'])
print('calc_fringefit is ',options['calc_fringefit'])
print('calc_fringefit2 is ',options['calc_fringefit2'])
if options['calc_fringefit2']:
    print('fringefit2_type is ',options['fringefit2_type'])
    if options['fringefit2_type']=='highsnr':
        print('\t fringefit_snr is',options['fringefit_snr'])
print('do_applycal is ',options['do_applycal'])
if options['do_applycal']:
    print('apply_fringefit2 is ', options['apply_fringefit2'])
print('rflag2 is ',options['rflag2'])
print('set_weights is ',options['set_weights'])
print('make_avg_ms is ',options['make_avg_ms'])
print('save_stokesi is ',options['save_stokesi'])
print('save_unscaled is ',options['save_unscaled'])


########################################################
# DEFINE PCAL FLAG FUNCTION
########################################################
#Since the phasecal for Oe and Ow is correlated, we need to flag out the peaks every 5 MHz.
#There are also spurious peaks at integer MHz frequencies.

def pcal_flag(msname, N_CHANNELS, N_SPW):
    CHANNEL_WIDTH = 32/N_CHANNELS
    #Flag pcal signal
    flagcmd = ""
    tb.open(msname+'/SPECTRAL_WINDOW')
    freqlist = tb.getcol("CHAN_FREQ") #Frequencies for each channel
    tb.close()
    for spw in range(N_SPW):
        chans = freqlist[:,spw]
        for ch,f in enumerate(chans):
            ispcal = False
            fm = f/1e6 #frequency in MHz
            # Check if channel is less than 0.1 MHz from pcal frequency (0,5,10,... MHz)
            #if (abs(fm-5*round(fm/5)) < 0.1):
            # Check if if channel is less than 0.1 MHz from integer MHz frequency (0,1,...). 
            # This should not be needed, every 5 MHz should be enough, but it seems we almost always
            # have spurious peaks at ever 5 MHz +- 1 MHz as well. Maybe Pcal malfunction.

            if N_CHANNELS >= 160: #If narrow channels, we need to flag the spurious peaks.
                if (abs(fm-1*round(fm/1)) < CHANNEL_WIDTH):
                    ispcal = True
            elif N_CHANNELS < 160: #If broad channels, the spurious peaks are less important (and we don't want to flag everything!)
                if (abs(fm-5*round(fm/5)) < CHANNEL_WIDTH):
                    ispcal = True
                elif (abs(fm-1*round(fm/1)) < CHANNEL_WIDTH) and int(round(fm))%5 == 1: #But keep flagging the channels 1 MHz above even 5:s, they need it!
                    ispcal = True
            if ispcal:
                flagcmd += ","+str(int(spw))+':'+str(int(ch))
    flagcmd = flagcmd[1:] # Remove initial comma
    flagdata(vis=msname, spw=flagcmd, mode="manual", antenna="OE&&OW")


########################################################
# SANITY CHECKS
########################################################
#Checking that all neccessary files exist, to avoid the script crashing later.

if '6.7' not in casalog.version():
    print(f"This script needs CASA 6.7, but the detected CASA version was {casalog.version()}. Cannot contiune.")
    #Note: The only task that strictly requires CASA 6.7 or later is appendantab. If you skip that task, you can disable this check (on your own risk!)
    sys.exit(1)

if (options['read_tsysdata'] or options['read_gaindata']) and not os.path.exists(antabfile):
    print(f"antab file {antabfile} is missing, cannot continue.")
    sys.exit(1)
if (not options['read_gaindata']) and (not os.path.exists(gaindata)) and options['calc_gain']:
    print(f"gaindata folder {gaindata} is missing and read_gaindata is False, cannot continue.")
    sys.exit(1)
if (options['create_ms']=='split') and (not os.path.exists(msname_full)):
    print("full measurement set is missing, cannot split data from it.")
    sys.exit(1) 
if (options['create_ms']=='none') and (not os.path.exists(msname)):
    print(f"Measurement set {msname} is missing and will not be created. Cannot continue.")
    sys.exit(1)
if options['prepare_calsource'] and not os.path.exists(infits_calsource):
    print(f"Calibrator fitsidi file {infits_calsource} is missing! Cannot continue.")
    sys.exit(1) 
if not options['make_avg_ms'] and not os.path.exists(ms_avg):
    print('Average MS does not exist and will not be created. Cannot continue.')
    sys.exit(1)
if not options['make_avg_ms'] and (options['do_applycal'] or options['rflag2']):
    print('Applycal or rflag2 will be done, but will not be saved in the averaged MS. This looks like a mistake, exiting.')
    sys.exit(1)

#Check that calibration tables exist, if the option to create them is False.
calc_options = [options['calc_tsys'], options['calc_gain'], options['calc_accor'], 
                options['calc_fringefit'], options['prepare_calsource']]
caltables = [tsystab, gaintab, accortab, fringetab, bptab_cal]
for option, table in zip(calc_options, caltables):
    if (not option) and (not os.path.exists(table)):
        print(f'Table {table} does not exist and is not set to be created. Cannot continue.')
        sys.exit(1)
if not os.path.exists(fringetab2) and not options['calc_fringefit2'] and options['do_applycal'] and options['apply_fringefit2']:
    print(f'Table {fringetab2} is set to be used, but it does not exist and is not set to be created. Cannot continue.')


#TIME THE SCRIPT
starttime = time()

########################################################
# IMPORT AND/OR PREPARE MEASUREMENT SET
########################################################

if options['create_ms']=='none':
    print(f'Using existing measurement set {msname}.')


elif options['create_ms']=='fitsidi':
    print(f'Loading complete vo data from {infitsidi}.')
    start_ms=time()
    if os.path.exists(msname_full):
        print(f'full measurement set {msname_full} already exists. Will not overwrite it.')
        sys.exit(1)
    importfitsidi(fitsidifile=infitsidi, vis=msname_full, scanreindexgap_s=10)
    end_ms=time()
    print(f'Measurement set created in {(end_ms-start_ms)/60:.2f} minutes = {(end_ms-start_ms)/3600:.2f} hours.')

    listobs(vis=msname_full,listfile=msname_full + '.listobs',overwrite=True)


#Split out only the data for Oe and Ow
if (options['create_ms']=='split' or options['create_ms']=='fitsidi'):

    #---- FIND THE SCANS THAT HAVE BOTH OE AND OW -------------------------
    #We only want to process the scans that have the OE-OW baseline, and ignore those that have only OE or only OW
    relevant_scans=''
    msmd.open(msname_full) 
    all_scans=msmd.scannumbers()
    for i in all_scans: #check if each scan contains both OE and OW
        antennas=msmd.antennasforscan(i) #ID numbers of antennas participating in the scan
        antnames=msmd.antennanames(antennas) #Convert to antenna names
        if 'OE' in antnames and 'OW' in antnames:
            relevant_scans+= str(i) + ','        
    msmd.done()
    relevant_scans=relevant_scans[:-1] #remove the last comma

    print('The scans containing the OE-OW baseline are:')
    if len(relevant_scans) < 100:
        print(relevant_scans)
    else:
        print(relevant_scans[:50] + '...' +  relevant_scans[-50:])
    scanfile=open('relevant_scans.txt','w')
    print(relevant_scans,file=scanfile)
    scanfile.close()
    #----------------------------------------------------------------------

    print(f'Splitting OTT baselines from measurement set {msname_full}, creating {msname}.')
    rmtables(msname)
    os.system('rm -r ' + msname + '.flagversions')
    split(vis=msname_full, outputvis=msname, datacolumn='data',scan=relevant_scans, antenna='OE&&OE;OE&&OW;OW&&OW')

    flagmanager(vis=msname, versionname='original',mode='save') #save the original flags


#The MS may contain flags form previous runs.  If 'flagversion'=='original', remove all flags.
if options['restore_flags'] and options['flagversion']=='original':
    try:
        flagmanager(vis=msname,versionname='original', mode='restore')
        os.system("rm -r " + msname+".flagversions") #remove old flags
        print(f'Flags were restored to original!')
        flagmanager(vis=msname, versionname='original',mode='save')
    except:
        flagdata(vis=msname,mode='unflag')
        flagmanager(vis=msname, versionname='original',mode='save')

#Restore flags to another state than original (e.g. if we did early flags before and want to restore to the state after rflag1).
elif options['restore_flags']:
    flagmanager(vis=msname,versionname=options['flagversion'], mode='restore')
    print(f'Flags were restored to version {options["flagversion"]}!')


#List observation info to file
listobs(vis=msname, listfile=e+'_OTT.listobs',overwrite=True)

########################################################
# IMPORT TSYS AND GAIN DATA
########################################################

#Tsys: Read from the antab file and append to the measurement set.
if options['read_tsysdata']:

    rmtables(msname+'.withtsys')
    appendantab(vis=msname,outvis=msname+'.withtsys',antab = antabfile, append_tsys = True, append_gc = False, overwrite = True)
    #overwrite the old msname, we don't want to keep two measurement sets!
    rmtables(msname)
    os.system('mv ' + msname + '.withtsys ' + msname)

    print(f'Tsys was read from file {antabfile}!')

#Read gain data from an external file. Appendantab does not work with my gain data!
if options['read_gaindata']:
    if not os.path.exists("../../casa-vlbi.git"):
        os.system("git clone https://github.com/jive-vlbi/casa-vlbi.git ../../casa-vlbi.git")
    os.system("casa --nologger -c ../../casa-vlbi.git/gc.py "+antabfile+" "+gaindata)

    if os.path.exists(gaindata):
        print("\nSuccessfully read gain data from antab file.\n")
    else:
        print("Failed to read gain data from antab file. Cannot continue.")
        sys.exit(1)


#Remove potential incorrect gaincurve from MS, allowing the data in gaindata to be used instead.
if options['remove_gaincurve']:
    tb.open(msname, nomodify=False)
    try:
        tb.removekeyword('GAIN_CURVE')
        tb.flush()
        print('Default gain curve was removed from MS.')
    except:
        print('Gain curve already removed!')
    tb.done()

########################################################
# FULL PROCESSING OF CALIBRATION SCANS
########################################################
if options['prepare_calsource']:
    start_cal=time()

    rmtables(msname_calsource)
    os.system("rm -r " + msname_calsource+".flagversions") #remove old flags
    print('Importing calibration scan measurement set...')
    importfitsidi(fitsidifile=infits_calsource, vis=msname_calsource,scanreindexgap_s=10)
    print('MS with calibration scans successfully imported.')
    listobs(vis=msname_calsource, listfile=msname_calsource + '.listobs',overwrite=True)

    ## Import Tsys

    rmtables(msname_calsource+'.withtsys')
    appendantab(vis=msname_calsource,outvis=msname_calsource+'.withtsys',antab = antabfile_cal, append_tsys = True, append_gc = False, overwrite = True)
    rmtables(msname_calsource)
    os.system('mv ' + msname_calsource + '.withtsys ' + msname_calsource)

    print(f'Tsys was read from file {antabfile_cal}!')

    #Remove potential incorrect gaincurve from MS, allowing the data in gaindata to be used instead.
    tb.open(msname_calsource, nomodify=False)
    try:
        tb.removekeyword('GAIN_CURVE')
        tb.flush()
        print('Default gain curve was removed from calibrator MS.')
    except:
        print('Gain curve already removed from calibrator MS!')
    tb.done()

    #Check the number of channels in calibration scan data
    tb.open(msname_calsource + '/SPECTRAL_WINDOW')
    channels=tb.getcol('CHAN_FREQ') # an array of shape (N_channels, N_spw)
    spw_freq_cal=(channels[0,:] + channels[-1,:])/2  #spw freq = (first_ch + last_ch)/2
    tb.close()

    NUM_CHANNELS_CAL=channels.shape[0]
    NUM_SPW_CAL = channels.shape[1]
    print(f'Calibration scan has {NUM_CHANNELS_CAL} channels per spw.')
    if spw_freq_cal[0]>spw_freq_cal[1]:
        print('Calibration scan seems to follow the conventional VO spw order.')
    else:
        print('WARNING! Calibration scan does not follow conventional VO spw order!')


    #Use gencal to make Tsys and gain tables for calibrator. Note: we can use the same gaindata table that we read for the main measurement set.
    gencal(vis=msname_calsource,caltype='tsys',caltable=tsystab_cal)
    # Smooth tsys table!
    smoothcal(vis=msname_calsource,tablein=tsystab_cal,caltable=tsystab_smooth_cal,smoothtype='median',smoothtime=60.0)
    if options['apply_smooth_tsys']:
        tsystab_cal = tsystab_smooth_cal

    gencal(vis=msname_calsource, caltype='gc', caltable=gaintab_cal, infile=gaindata)

    if options['flag_spw_01']:
        #Flag spw 0 and 1 (highest freq in band A) as they are always affected by UER (Unwanted Electromagnetic Radiation).
        flagdata(vis=msname_calsource, spw='0,1', mode="manual")

    if options['custom_flags']:
        #Use flags from flagfile. Can be e.g. times when a telescope was off-source.
        if os.path.exists(flagfile_cal):
            try:
                flagdata(vis=msname_calsource,mode='list',inpfile=flagfile_cal)
            except:
                print(f'WARNING! custom flags failed. Probably because none of the flags in {flagfile_cal} are within the timerange of the experiment.')
        else:
            print(f'WARNING! {flagfile_cal} does not exist. Will not apply custom flags to flux calibration sources!')
         

    #Do early flagging of calibration data
    pcal_flag(msname_calsource,NUM_CHANNELS_CAL,NUM_SPW_CAL)

    #save these flags
    flagmanager(vis=msname_calsource, mode='save', versionname='pcal')

    #Flag edges of all spws
    if NUM_CHANNELS_CAL==320:
        flagdata(vis=msname_calsource, mode="manual", spw="*:0~31;304~319", antenna="OE&&OW")
    elif NUM_CHANNELS_CAL==160:
        flagdata(vis=msname_calsource, mode="manual", spw="*:0~16;152~159", antenna="OE&&OW")
    elif NUM_CHANNELS==128:
        flagdata(vis=msname_calsource, mode="manual", spw="*:0~12;121~127", antenna="OE&&OW")
    elif NUM_CHANNELS==64:
        flagdata(vis=msname_calsource, mode="manual", spw="*:0~5;57~59", antenna="OE&&OW")
    else:
        print('Unexpected number of channels, cannot do edgeflag of calibration source!')

    flagmanager(vis=msname_calsource, mode='save', versionname='edgeflag')

    #Make autocorrelation correction table.
    #Note: We get one solution per scan, spw and polarization.
    accor(vis=msname_calsource, caltable = accortab_cal, corrdepflags=True, solint="inf")

    #Rflag
    flagdata(vis=msname_calsource, mode="rflag", freqdevscale=5.0, timedevscale=5.0, antenna="OE&OW", datacolumn="data", extendflags=True, ntime="scan")
        
    flagmanager(vis=msname_calsource, mode='save', versionname='rflag1')

    #Fringefit: Do instrumental fringefit on all the calibrator sources.
    #First, pick out all the flux calibrators observed in this experiment.
    #Also pick out the observation date.
    msmd.open(msname_calsource)
    fieldnames_cal = msmd.fieldnames()
    time_mjd = msmd.timerangeforobs(0)['begin']['m0']['value']
    msmd.done()
    fringefit_sources = []
    for source in fieldnames_cal:
        if source in ['3C286','3C147','3C295']:
            fringefit_sources.append(source)

    #If this experiment was observed after 2025-06-30, it was correlated with an empirical clock offset of 0.11 us. This causes smaller delays than previously!
    if time_mjd < 60856:
        delaywindow = [100,160]
    else:
        delaywindow = [-10,50]
    if e == 'vo5134':
        delaywindow = [-10,50]

    solint='180s'
    fringefit(vis=msname_calsource, field=','.join(fringefit_sources), refant="OE", caltable=fringetab_cal, gaintable=[gaintab_cal,tsystab_cal,accortab_cal], interp = ['',',nearest',''], solint=solint, globalsolve=True, corrdepflags=True, delaywindow = delaywindow, zerorates=True, antenna="OE&&OW")

    #bandpass
    bandpass(vis=msname_calsource, gaintable=[gaintab_cal,tsystab_cal,accortab_cal,fringetab_cal], interp=['',',nearest','',''], scan='1', refant="OE", solnorm=True, caltable=bptab_cal, solint = "180s", fillgaps=10, minblperant=1, corrdepflags=True)

    flagmanager(vis=msname_calsource, mode='save', versionname='before_first_applycal') 
    #save the flags here so we can restore to this point after.


    #Applycal
    applycal(vis = msname_calsource, gaintable = [gaintab_cal, tsystab_cal, accortab_cal, fringetab_cal, bptab_cal],
        interp=["", ",nearest", "","linearperscan",""]) #note: we may have a weak source in the calibration file, which should not be used! "linearperscan" makes sure it is not used.

    flagmanager(vis=msname_calsource, mode='save', versionname='after_applycal')

    #rflag2
    flagdata(vis=msname_calsource, mode="rflag", freqdevscale=7.0, timedevscale=7.0, antenna="OE&OW", datacolumn="corrected", extendflags=True, ntime="scan")

    flagmanager(vis=msname_calsource, mode='save', versionname='rflag2')

    #Set weights based on variance
    statwt(vis=msname_calsource,timebin='10s')

    #split: average over time and channel. ALSO remove autocorrelations
    rmtables(ms_avg_cal)
    split(vis=msname_calsource, outputvis = ms_avg_cal, width=NUM_CHANNELS_CAL, timebin = "3600s", datacolumn="corrected", antenna="OE&OW", keepflags=False)

    end_cal=time()
    print(f'Processing of calibration scan completed in {(end_cal-start_cal)/60} minutes.')
####################################################################
##End of calibration scan processing

    

########################################################
# INITIAL DATA CHECK ON MAIN MEASUREMENT SET
########################################################


tb.open(msname + '/SPECTRAL_WINDOW')
channels=tb.getcol('CHAN_FREQ') # an array of shape (N_channels, N_spw)
tb.close()

#Check number of spw and number of channels!
NUM_SPW=channels.shape[1]
NUM_CHANNELS=channels.shape[0]
CHANNEL_WIDTH = 32/NUM_CHANNELS #MHz

print(f'Dataset has {NUM_SPW} spectral windows, each with {NUM_CHANNELS} channels. The channel width is {CHANNEL_WIDTH:.2f} MHz.')

spw_freq=(channels[0,:] + channels[-1,:])/2  #spw freq = (first_ch + last_ch)/2
print('\nCenter frequencies of spectral windows:')
print(spw_freq, '\n')
#Save spw frequencies to file
np.save('frequencies.npy',spw_freq)

#Check the ordering of the spws.
#In VO data, it can be either (from low to high freq) 0,1,2,3,...,31  or  7,6,5,.., 15,14,13,..., 23,22,21..., 31,30,29...24.
#Tsys values in antab file assumes the latter option, so if the MS has the former, we need to apply a spw-map in Tsys calibration.
if spw_freq[1]<spw_freq[0]: #conventional case, 7,6,5,.., 15,14,13,..., 23,22,21..., 31,30,29...24.
    spwmap_tsys=[]
    OTT_order=True
else:
    OTT_order=False
    print('Unconventional spw order detected! creating spw map for Tsys calibration.')
    #ex: spwmap=[0,0,1,1] means apply the caltable solutions from spw = 0 to the spw 0,1 and spw 1 to spw 2,3.
    spwmap_tsys=[7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8,23,22,21,20,19,18,17,16,31,30,29,28,27,26,25,24]
    print(spwmap_tsys)
    print('Sanity check:')
    print(f'spw 0 has freq {spw_freq[0]/1e9:.5f} GHz and gets the Tsys values for spw {spwmap_tsys[0]}.')
    print(f'spw 1 has freq {spw_freq[1]/1e9:.5f} GHz and gets the Tsys values for spw {spwmap_tsys[1]}.')
    print(f'spw 7 has freq {spw_freq[7]/1e9:.5f} GHz and gets the Tsys values for spw {spwmap_tsys[7]}.')
    print(f'spw 8 has freq {spw_freq[8]/1e9:.5f} GHz and gets the Tsys values for spw {spwmap_tsys[8]}.')
    print(f'spw 15 has freq {spw_freq[15]/1e9:.5f} GHz and gets the Tsys values for spw {spwmap_tsys[15]}.')
    print(f'spw 16 has freq {spw_freq[16]/1e9:.5f} GHz and gets the Tsys values for spw {spwmap_tsys[16]}.')
    print(f'spw 23 has freq {spw_freq[23]/1e9:.5f} GHz and gets the Tsys values for spw {spwmap_tsys[23]}.')
    print(f'spw 24 has freq {spw_freq[24]/1e9:.5f} GHz and gets the Tsys values for spw {spwmap_tsys[24]}.')
    print(f'spw 30 has freq {spw_freq[30]/1e9:.5f} GHz and gets the Tsys values for spw {spwmap_tsys[30]}.')
    print(f'spw 31 has freq {spw_freq[31]/1e9:.5f} GHz and gets the Tsys values for spw {spwmap_tsys[31]}.')


########################################################
# PICK REFERENCE SCANS FOR FRINGEFIT
########################################################
#Function to compute az-el coordinates of the source. If the source has too low elevation, we do not want to use it for instrumental fringe fit.
def altaz_coord(time_mjd, ant_coord, source_coord):
    """
    Calculates horizontal coordinates of a source as seen from the given antenna location.
    
    Parameters
    ----------
    time_mjd : float
        Time in modified Julian days
    ant_coord: tuple
        (lat,lon) of antenna in unit of radians
    source_coord: tuple
        (RA,dec) of source in unit of radians

    Returns
    -------
    (az,el) as floats (unit degrees).

    """
    #Convert indata to useful formats
    time_dt = Time(time_mjd,format='mjd').datetime
    telescope_location = EarthLocation(lat=ant_coord[0]*u.rad,lon=ant_coord[1]*u.rad)
    source_location = SkyCoord(ra=source_coord[0]*u.rad,dec=source_coord[1]*u.rad)
    
    altaz_frame = AltAz(location=telescope_location,obstime=time_dt)
    
    altaz_coord = source_location.transform_to(altaz_frame)
    
    return altaz_coord.az.value,altaz_coord.alt.value


#def altaz_coord(time_mjd,ant_coord,source_coord):
#"""
#Dummy function to use if astropy is not available. Will return a high elevation value, so that no scan is filtered out based on low elevation.
#"""
#    return (0,999)
#--------------------------------------------------------


msmd.open(msname) 
fieldnames = msmd.fieldnames() #get all filenames

#get antenna position, to be used for elevation calculation
pos=msmd.antennaposition('OE')
lon_rad = pos['m0']['value'] #longitude in radians
lat_rad=pos['m1']['value'] #latitude in radians
antenna_coord = (lat_rad,lon_rad)
print(f'OE has lat={lat_rad*180/np.pi} deg and lon={lon_rad*180/np.pi} deg.')

#get MS summary, which contains all sky coordinates of sources
ms.open(msname)
mssum = ms.summary()
ms.close()

#Pick a bright, long scan for instrumental fringefit
#Check each source on the list. If we find a source in the measurement set, check how long each scan is. If the scan is <29 s, ignore it.
ok_refsources= ['OJ287','0552+598','0133+476','1156+295','0059+581'] #priority list of reference sources.

ref_source=''
frscan_instr = -1
for s in ok_refsources:
    if s in fieldnames:
        ref_source=s #We found a possible reference source!

        #Now find a long enough scan.
        ref_scans=msmd.scansforfield(ref_source)

        for scan in ref_scans:
            times = msmd.timesforscan(scan) #get time in MJDseconds for every integration time in this scan
            duration = times[-1]-times[0]

            if duration < 28:
                #ignore any short scans
                continue 
                
            #Compute elevation of source
            begin_time = mssum['scan_'+str(scan)]['0']['BeginTime'] #time in MJD
            field_id = mssum['scan_'+str(scan)]['0']['FieldId']
            coord = mssum['field_'+str(field_id)]['direction']
            ra_rad = coord['m0']['value']
            dec_rad = coord['m1']['value']
            source_coord = (ra_rad,dec_rad)
            
            (az,el) = altaz_coord(begin_time,antenna_coord,source_coord) #azimuth, elevation in degrees

            if el < 30:
                #ignore any low-elevation scans
                print(f'Will not use scan {scan} of {s} as instrumental fringefit scan, as elevation is only {el:.1f} degrees.')
                continue
                
            frscan_instr = scan
            print(f"Reference scan for instrumental fringefit is scan {scan}, source {ref_source}. Duration: {duration+1:.1f} seconds. Elevation: {el:.1f} degrees.")
            break
            #stop looping through scans when we found a good one
        if frscan_instr==-1:
            print(f"No acceptable fringe finder scans found for source {s}!")
        else: 
            break
            #stop looping through sources when we found one good fringefinder scan

if ref_source=='':
    raise Exception('Source list contains no valid reference source:', ok_refsources)
elif frscan_instr == -1:
    raise Exception("No appropriate calibration scan found! All scans of the possible fringe finders are either too short of have too low elevation.")


#Pick out all long scans to use for multiband fringefit.
if options['fringefit2_type']=='calibrators':

    #Check each scan in the MS. If it is long (2 minutes), it is a calibration scan.
    all_scans=msmd.scannumbers()

    ref_source=''
    all_ref_sources = []
    frscans=[]

    for scan in all_scans:
        times = msmd.timesforscan(scan) #get time in MJDseconds for every integration time in this scan

        duration = times[-1]-times[0]
        if duration > 89:
            frscans.append(scan)
            ref_source = fieldnames[msmd.fieldsforscan(scan)[0]]
            previous_scan_time = times[-1]
            if not ref_source in all_ref_sources:
                all_ref_sources.append(ref_source)
            print(f"Reference scan for fringefit found: scan {scan}, source {ref_source}. Duration: {duration+1:.1f} seconds.")

    if frscans == []:
        raise Exception("No appropriate calibration scan found! All scans are less than 90 seconds.")

    #Also pick the first and last scan that have data
    for i in range(100):
        firstscan = all_scans[i]
        nrows = mssum['scan_'+str(firstscan)]['0']['nRow']
        if nrows > 500:
            break
    for i in range(100):
        lastscan = all_scans[-i]
        nrows = mssum['scan_'+str(lastscan)]['0']['nRow']
        if nrows > 500:
            break

    print("Using all long scans of ", all_ref_sources, "to do multi-band fringefit, i.e. frscans=", frscans)
    print(f"Also using first and last scan with data: {firstscan} and {lastscan}.")

    scan_strings=[str(f) for f in frscans]
    frscans_string = ','.join(scan_strings)
    frscans_string = str(firstscan) + ',' + frscans_string + ',' + str(lastscan) #add the first and last scan as well, to get interpolated solutions also in the first and last hour

elif options['fringefit2_type']=='highsnr':
    print(f'Multi-band fringefit will be done on all scans with an SNR above {options["fringefit_snr"]}. Other scans will be flagged!')

########################################################
# FIRST ROUND OF FLAGGING
########################################################

if options['flag_spw_01']:
    #Flag spw 0 and 1 (highest freq in band A) as they are always affected by unwanted electromagnetic radiation.
    if len(spwmap_tsys)==0:
        flagdata(vis=msname, spw='0,1', mode="manual")
    else: #If spws are in unconventional order, the bad spws are 6 and 7.
        flagdata(vis=msname, spw='6,7', mode="manual")

if options['custom_flags']:
    #Use flags from flagfile. Can be e.g. times when a telescope was off-source.
    if os.path.exists(flagfile):
        try:
            flagdata(vis=msname,mode='list',inpfile=flagfile)
        except:
            print(f'WARNING! custom flags failed. Probably because none of the flags in {flagfile} are within the timerange of the experiment.')
    else:
        print(f'WARNING! {flagfile} does not exist. Will not apply custom flags!')


if options['earlyflags']:
    #Flag phasecal peaks
    pcal_flag(msname, NUM_CHANNELS, NUM_SPW)

    flagmanager(vis=msname, mode='save', versionname='pcal')

    #Flag edges of all spws
    if NUM_CHANNELS==320:
        flagdata(vis=msname, mode="manual", spw="*:0~31;304~319", antenna="OE&&OW")
    elif NUM_CHANNELS==160:
        flagdata(vis=msname, mode="manual", spw="*:0~16;152~159", antenna="OE&&OW")
    elif NUM_CHANNELS==128:
        flagdata(vis=msname, mode="manual", spw="*:0~12;121~127", antenna="OE&&OW")
    elif NUM_CHANNELS==64:
        flagdata(vis=msname, mode="manual", spw="*:0~5;57~59", antenna="OE&&OW")
    else:
        print(f'Unexpected number of channels: {NUM_CHANNELS}. Cannot do edgeflag on main measurement set!!')
        sys.exit(1)

    flagmanager(vis=msname, mode='save', versionname='edgeflag')


if options['earlyflags']:
    #Rflag: flag outliers
    flagdata(vis=msname, mode="rflag", freqdevscale=5.0, timedevscale=5.0, antenna="OE&OW", datacolumn="data", extendflags=True, ntime="scan")
        
    flagmanager(vis=msname, mode='save', versionname='rflag1')


########################################################
# MAKE CALIBRATION TABLES
########################################################

#Use gencal to make Tsys and gain tables
if options['calc_tsys']:
    gencal(vis=msname,caltype='tsys',caltable=tsystab)

if options['smooth_tsys']:
    #Smooth Tsys in time to remove noise
    smoothcal(vis=msname,tablein=tsystab,caltable=tsystab_smooth,smoothtype='median',smoothtime=60.0)

if options['apply_smooth_tsys']:
    tsystab = tsystab_smooth
    #The smooth Tsys table will be used.

if options['calc_gain']:
    gencal(vis=msname, caltype='gc', caltable=gaintab, infile=gaindata)

#Make autocorrelation correction table. One solution per scan, polarization and spw.
if options['calc_accor']:
    accor(vis=msname, caltable = accortab, corrdepflags=True, solint="inf")


#Instrumental fringefit: 
#note: delaywindow is set to [100,160] as correlators typically keep the 0.11 us offset between Oe and Ow.
if options['calc_fringefit']:
    fringefit(vis=msname, scan=str(frscan_instr), refant="OE", gaintable = [gaintab, tsystab, accortab], interp = ["",",nearest", ""],spwmap=[[],spwmap_tsys,[]], caltable=fringetab, solint="inf", globalsolve=True, corrdepflags=True, delaywindow = [100,160], zerorates=True, antenna="OE&&OW")


#Multi-band fringefit:
if options['calc_fringefit2']:
    if options['fringefit2_type']=='calibrators':
        fringefit(vis=msname, scan=frscans_string, refant="OE", caltable=fringetab2, gaintable=[gaintab, tsystab, accortab,fringetab],interp = ["",",nearest", "",""],spwmap=[[],spwmap_tsys,[],[]], solint="inf", globalsolve=True, corrdepflags=True, delaywindow = [-5,5], zerorates=True, antenna="OE&&OW",combine='spw',corrcomb='stokes') #combine X and Y polarization!

    elif options['fringefit2_type']=='all':
        fringefit(vis=msname, scan='', refant="OE", caltable=fringetab2, gaintable=[gaintab, tsystab, accortab,fringetab],interp = ["",",nearest", "",""],spwmap=[[],spwmap_tsys,[],[]], solint="inf", globalsolve=True, corrdepflags=True, delaywindow = [-10,10], zerorates=True, antenna="OE&&OW",combine='spw',corrcomb='stokes') #combine X and Y polarization!

    elif options['fringefit2_type']=='highsnr':
        fringefit(vis=msname, scan='', refant="OE", caltable=fringetab2, gaintable=[gaintab, tsystab, accortab,fringetab],interp = ["",",nearest", "",""],spwmap=[[],spwmap_tsys,[],[]], solint="inf", globalsolve=True, corrdepflags=True, delaywindow = [-5,5], zerorates=True, antenna="OE&&OW",combine='spw',corrcomb='stokes',minsnr = options['fringefit_snr']) #combine X and Y polarization!



#Note: We don't sove for bandpass here. See below.


########################################################
# APPLY CALIBRATION
########################################################

#Final applycal! NOTE: We use the bandpass solution from the calibrator, to make sure the scalefactors are consistent with the bandpass solution.

if options['do_applycal']:

    flagmanager(vis=msname, mode='save', versionname='before_applycal') 
    #save the flags here so we can restore to this point after.

    if options['apply_fringefit2']:
        applycal(vis = msname, gaintable = [gaintab, tsystab, accortab, fringetab,fringetab2, bptab_cal],
            interp=["", ",nearest", "","","",""],spwmap=[[],spwmap_tsys,[],[],[0]*NUM_SPW,spwmap_tsys])
    else:
        #Only instrumental fringefit!
        applycal(vis = msname, gaintable = [gaintab, tsystab, accortab, fringetab, bptab_cal],
            interp=["", ",nearest", "","",""],spwmap=[[],spwmap_tsys,[],[],spwmap_tsys])


    flagmanager(vis=msname, mode='save', versionname='after_applycal')  

########################################################
# SECOND ROUND OF FLAGGING
########################################################

#rflag2
if options['rflag2']:
    flagdata(vis=msname, mode="rflag", freqdevscale=7.0, timedevscale=7.0, antenna="OE&OW", datacolumn="corrected", extendflags=True, ntime="scan")

    flagmanager(vis=msname, mode='save', versionname='rflag2')

########################################################
# PREPARE AND AVERAGE DATA
########################################################

if options['set_weights']:
    #set statistical weights based on standard devation of visibilities
    statwt(vis=msname,timebin='10s')

#split: average over time and channel. Also remove autocorrelations
if options['make_avg_ms']:
    rmtables(ms_avg)
    split(vis=msname, outputvis = ms_avg, width=NUM_CHANNELS, timebin = "3600s", datacolumn="corrected", antenna="OE&OW", keepflags=False)


########################################################
# CALCULATE FLUX DENSITIES AND SAVE DATA
########################################################

if not os.path.exists(fluxdir):
    os.makedirs(fluxdir)

vis = ms_avg
msmd.open(vis) 
fieldnames = msmd.fieldnames()
print('fieldnames:',fieldnames)
msmd.done()
if 'fm' not in e:
    msmd.open(ms_avg_cal)
    fieldnames_cal=msmd.fieldnames()
    print('Calibrators: ',fieldnames_cal)
else:
    fieldnames_cal = [] #for a FM experiment, we don't have any sourcesin the additional measurement set
msmd.done()

#Make list of the MS where to find each source:
ms_list=[ms_avg]*len(fieldnames) + [ms_avg_cal]*len(fieldnames_cal)
fieldname_list=list(fieldnames)  + list(fieldnames_cal)   

data = {} #saves stokes I data
corr_data={} #saves XX, XY, YX, YY data
unflagged={} #stores which spws are unflagged for each field.

flagged_file=open('flagged_spws.txt','w') #a file that logs which spws are totally flagged.

spw_OTT=list(range(0,32))
if spwmap_tsys==[]: #If no reordering of spws needed:
    spw_VO=spw_OTT
else:
    spw_VO=spwmap_tsys #[7,6,5,4,3,2,1,0,15,14,13...]

outlier_percent = []

for field, vis in zip(fieldname_list, ms_list):
    phase_removed=False
    idata = []
    full_data=[] #including all corrs!
    unflagged_spw=[] #keep track of which spws actually have data.
    flagged_spw=[]
    print("Processing field " + field)
    print(field, file=flagged_file)

    #Some explanation of the spw numbering:
    #spw1 decides what the spw should be named in the outfile.
    #spw2 is the corresponding name for that spw in the MS. For example, we want the highest frequency spw in band A to be called 0,
    #but it may be stored as '7' in the MS.
    #spw_OTT=[0,1,2,3...], spw_VO=[7,6,5,4,3,2,1,0,15,14,13,...] if the order is unconventional.
    for spw1, spw2 in zip(spw_OTT,spw_VO):

        if field in fieldnames_cal: #if this is a calibrator source:
            spw2=spw1 #spw order should not be permutated. It is already in the OTT frequency order!

        ms.open(vis)
        ms.selectinit(datadescid=spw2) #select data from one spw only
        staql={'field':field} #selct data from one field only
        ms.msselect(staql)
        msum =ms.summary()
        if msum:
        # If we got some data (all could be flagged)
            unflagged_spw.append(spw1)
            start = msum['BeginTime']
            stop = msum['EndTime']
            avgtime = start+0.5*(start-stop)
            d = ms.getdata(["amplitude","weight","axis_info"], ifraxis=True) # ifraxis to reformat into baseline structure
            f = d["axis_info"]["freq_axis"]["chan_freq"][0][0] #spw frequency
            corrs = list(d["axis_info"]["corr_axis"])
            weights = d["weight"] #shape e.g. (4,1,3): (pol, baseline, scan)
            xx = d["amplitude"][corrs.index("XX")][0][0] # Pol XX, frequency 0, baseline 0 (only one freq and one baseline ramains in these data)
            xy = d["amplitude"][corrs.index("XY")][0][0] #array structure: (pol, channel, baseline, time)
            yx = d["amplitude"][corrs.index("YX")][0][0]
            yy = d["amplitude"][corrs.index("YY")][0][0]

            try:
                len(xx)
            except: #if there is only one scan of a source, xx will be a float. Convert it to an array to avoid issues later.
                xx = np.array([xx])
                xy = np.array([xy])
                yx = np.array([yx])
                yy = np.array([yy])

            #remove points that deviate too much from the median amplitude
            is_outlier = np.zeros(len(xx),dtype='bool')
            for v in [xx,yy]:
                if len(xx)<2:
                    continue #cannot do outlier filtering on one datapoint!
                median_val = np.nanmedian(v)
                mad_val = MAD(v)

                for n in range(len(is_outlier)):
                    if np.abs(v[n]-median_val)>mad_num*mad_val:
                        is_outlier[n]=True

            outlier_percent.append(np.round(sum(is_outlier)/len(xx)*100,2))

            xx = xx[~is_outlier]
            xy = xy[~is_outlier]
            yx = yx[~is_outlier]
            yy = yy[~is_outlier]

            weights_reshaped = weights.reshape((4,int(weights.size/4))) #shape (4,X) where X is the number of scans


            #weighted average of amplitudes
            xx_avg = np.sum(xx*weights_reshaped[corrs.index("XX"),~is_outlier])/np.sum(weights_reshaped[corrs.index("XX"),~is_outlier])
            xy_avg = np.sum(xy*weights_reshaped[corrs.index("XY"),~is_outlier])/np.sum(weights_reshaped[corrs.index("XY"),~is_outlier])
            yx_avg = np.sum(yx*weights_reshaped[corrs.index("YX"),~is_outlier])/np.sum(weights_reshaped[corrs.index("YX"),~is_outlier])
            yy_avg = np.sum(yy*weights_reshaped[corrs.index("YY"),~is_outlier])/np.sum(weights_reshaped[corrs.index("YY"),~is_outlier])

            #calculate new weights
            weights_nooutliers = weights_reshaped[:,~is_outlier]
            weights_total = np.sum(weights_nooutliers,axis=1)
            final_errors=1/np.sqrt(weights_total) #error on the averaged visibilities

            #Calculate and store Stokes I
            stokesi =  0.5*(np.abs(xx_avg)+np.abs(yy_avg))
            idata.append([f, stokesi, avgtime, spw1]) #Note: Here we have changed SPW order to the OTT order!

            #Store XX,XY,YX,YY values and their errors
            full_data.append([avgtime,f,np.abs(xx_avg),np.abs(xy_avg),np.abs(yx_avg),np.abs(yy_avg), final_errors[0],final_errors[1],final_errors[2],final_errors[3], spw1])

            if spw1%8==0 and field=='OJ287': #sanity check, make sure the spw mapping is correct
                print(f'We read spw {spw2} from the MS with freq={f/1e9:.5f} GHz, and store it as spw {spw1}.')
        else:
            idata.append([]) #append an empty list to show that a spectral window is missing.
            full_data.append([])
            flagged_spw.append(spw1) #Note: This is in OTT order

        ms.done()

    #end of "for spw1, spw2 in zip(spw_OTT,spw_VO)"
        
    data[field] = idata
    corr_data[field]=full_data
    unflagged[field] = unflagged_spw
    print('Unflagged spw for field ', field, ': ', unflagged_spw)
    print(flagged_spw, file=flagged_file)
    

    #SAVE DATA (WITHOUT SCALEFACTORS) TO FILE ###########
    if options['save_unscaled'] and options['save_stokesi']:
    
        outfile=open(f'{fluxdir}/flux_{field}.txt','w')
        print('freq stokesI time spw', file=outfile)
        for row in fdata:
            if len(row)!=0: #If there is data for this spw
                print(f'{row[0]:16.1f} {row[1]:8.3f} {row[2]:14.5f} {row[3]}',file=outfile)
        outfile.close()

    if options['save_unscaled']:
        #Save the individual polarizations to file
        outfile=open(f'{fluxdir}/corrs_{field}.txt','w')
        print('time           \t freq             \t XX       XY       YX       YY       \t XXerr    XYerr    YXerr    YYerr    \t spw', file=outfile)
        for row in full_data:
            if len(row)!=0: #If there is data for this spw
                print(f'{row[0]:14.5f} {row[1]:16.1f} {row[2]:8.3f} {row[3]:8.3f} {row[4]:8.3f} {row[5]:8.3f} \t {row[6]:8.3f} {row[7]:8.3f} {row[8]:8.3f} {row[9]:8.3f} \t {row[10]}',file=outfile)
        outfile.close()
    ####################################################
    
#end of "for field, vis in zip(fieldname_list, ms_list)"


########################################################
# CALCULATE SCALE FACTORS
########################################################

#First, check if there is a calibration source. If not, no scaled flux files are created.
ampcals=["3C286", "3C147", "3C295"]
used_calsources=[] #Calibration sources observed in this exp
for s in ampcals:
    if s in fieldname_list:
        used_calsources.append(s)

if len(used_calsources)==0:
    print('No calibration source was observed. Cannot calculate scale factors.')

    endtime=time()
    print(f'Script finished in {(endtime-starttime)/60:.2f} minutes = {(endtime-starttime)/3600:.2f} hours.')
    sys.exit(0)
    
        
# define Perley-Butler models, x in GHz #####################
def pb_3c286(x):
    """Spectrum for 3c286, Perley and Butler (2017) """
    return 10**(1.2481-0.4507*np.log10(x)-0.1798*(np.log10(x))**2+0.0357*(np.log10(x))**3)
def pb_3c147(x):
    """Spectrum for 3c147, Perley and Butler (2017) """
    return 10**(1.4516-0.6961*np.log10(x)-0.2007*(np.log10(x))**2+0.0640*(np.log10(x))**3-0.0464*(np.log10(x))**4+0.0289*(np.log10(x))**5)
def pb_3c295(x):
    """Spectrum for 3c295, Perley and Butler (2017) """
    return 10**(1.4701-0.7658*np.log10(x)-0.2780*(np.log10(x))**2+0.0347*(np.log10(x))**3-0.0399*(np.log10(x))**4)
 ############################################################
        
ampscales = []
for cal in ampcals:
    if cal in data.keys():
        calscale = []

        for spw in range(0,32):

            row=data[cal][spw]
            if len(row)==0: #if the data row is empty, i.e. there is no data for this spw
                calscale.append(float('nan'))
                continue
            
            frq = row[0] / 1e9 # converted to GHz
            amp = row[1]
            
            if cal == "3C286":
                mod = pb_3c286(frq)
            elif cal == "3C147":
                mod = pb_3c147(frq)
            elif cal == "3C295":
                mod = pb_3c295(frq)
            calscale.append(mod/amp)
            
        ampscales.append(calscale)

ampscales = np.array(ampscales)
scalemean = np.nanmean(ampscales, axis=0)
scalestd = np.nanstd(ampscales, axis=0)

#Save detailed scalefactors to file
scaleinfo = open(f"{scalefile}_all.txt", 'w')

header='spw '
formatstring='{0} '
values=[list(range(0,32))] #data to be written to the file. In the beginning, it contains only spectral window indices.
#We may have either 1, 2 or 3 flux calibrators. Add a column to the file format string for each calibrator we observed.
for k,s in enumerate(used_calsources):
    header= header + s + ' '
    values.append(ampscales[k]) #all 32 scale factors calculated from source s
    formatstring = formatstring + '{'+ str(k+1)+ ':.4f} ' #Add a place for individual scale factor

header = header + 'scalemean scalestd'
values.append(scalemean)
values.append(scalestd)
formatstring = formatstring + '{' + str(k+2) +':.4f} {' + str(k+3) + ':.4f}' #Add the places for std and mean scale factor

print(header, file=scaleinfo)
for i in range(32):
    values_thisrow=[v[i] for v in values]
    print(formatstring.format(*values_thisrow),file=scaleinfo)
scaleinfo.close()
np.save(f'{scalefile}_all.npy',ampscales) #save the same data as a .npy file also.

########################################################
# APPLY SCALE FACTORS
########################################################
    

#Correct the stokes I amplitudes using the scalefactors and save them to file.
if options['save_stokesi']:
    for field in data.keys():
        print("Writing data for field " + field)
        
        outfile=open(f'{fluxdir}/scaled_flux_{field}.txt','w')
        print('freq stokesI time', file=outfile)
        for spw in unflagged[field]: #only loop through the spws that have data, making sure corrections are applied in the right place!
            row=data[field][spw]
            if len(row)==0: #this check is probably unnecessary given that we don't loop though all
                continue
            
            flux=row[1]
            corr_flux=flux*scalemean[spw]
            print(f'{row[0]:16.1f} {corr_flux:8.3f} {row[2]:14.5f}',file=outfile)
        outfile.close()
    
    
#Correct the individual corrs and save to file.
#Note: errors are first scaled up (or down) using the scale factor, then the uncertainty of the scalefactor itself is added.
for field in data.keys():
    
    outfile=open(f'{fluxdir}/scaled_corrs_{field}.txt','w')
    print('time           \t freq             \t XX       \t XY       \t YX       \t YY       \t XXerr    \t XYerr    \t YXerr    \t YYerr    \t spw', file=outfile)
    for spw in unflagged[field]: #only loop through the spws that have data, making sure corrections are applied in the right place!
        updated_row = []
        row=corr_data[field][spw]
        updated_row = updated_row + row[:2] #copy time and freq
        for num in row[2:10]:
            updated_row.append(scalemean[spw]*num) #scale fluxes and errors by scalefactor

        rel_scale_err = (scalestd[spw]/np.sqrt(len(used_calsources)))/scalemean[spw] #relative scale error
        for i in range(4):
            flux = updated_row[2+i]
            err = updated_row[6+i]
            new_err = np.sqrt(err**2 + (flux*rel_scale_err)**2)
            updated_row[6+i] = new_err
        
        print(f'{updated_row[0]:14.5f} {updated_row[1]:16.1f} {updated_row[2]:8.3f} {updated_row[3]:8.3f} {updated_row[4]:8.3f} {updated_row[5]:8.3f} {updated_row[6]:8.3f} {updated_row[7]:8.3f} {updated_row[8]:8.3f} {updated_row[9]:8.3f} {row[10]}',file=outfile)
    outfile.close()  


#If the scaling factors are too large, warn!
bad_scale = scalemean > 2
bad_spw = [i for i in range(32) if bad_scale[i]]
if len(bad_spw)>0:
    print('WARNING! SCALE FACTORS > 2 ENCOUNTERED IN SPW ', bad_spw )
    print(scalemean[bad_scale])

endtime=time()

print('Outlier percentages:')
print(f'Median: {np.nanmedian(outlier_percent):.1f} %, max: {np.nanmax(outlier_percent):.1f} %')

print(f'Script finished in {(endtime-starttime)/60:.2f} minutes = {(endtime-starttime)/3600:.2f} hours.')


