# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 16:05:17 2015

@author: dongyu
"""
import os
import utm
import hydro_wrapper
import oilspill_wrapper
from DownloadTool import downloadROMS
import pdb



def HyosPy(starttime, endtime, period, mode=1):
    """
    input:
        starttime='2014-02-21'
        endtime='2014-02-22'
        period=24         oil spill simulation duration;
        mode = 1
        
        Mode 1:
            download ROMS output from the server and run GNOME, have a quicker look (takes a few seconds)
        Mode 2:
            download river, wind, tidal, ROMS initial data, run SUNTANS and GNOME, have an accurate look
            (takes a few minutes)
        Mode 3: 
            a combination of mode 1 and 2, use the output of both ROMS and SUNTANS, duplicate the particles
            that enters the overlap region
    
        Some suggestions: 
            1) if the initial location is out of SUNTANS domain, run both mode 1 and mode 3
            2) if the initial location is in SUNTANS domain, simply run mode 2
    """
    #starttime='2014-02-21'
    #endtime='2014-02-22'
    dire=os.getcwd()
    
    if mode==1:
        
        downloadROMS(starttime,endtime)
        ROMS_file=dire+'/DATA/'+'txla_subset_HIS.nc'
        ROMS_out=dire+'/GNOME/'+'hiroms_ss_rho.nc'        
        oilspill_wrapper.init_model(opt='ROMS')
        oilspill_wrapper.HIROMS(ROMS_file,ROMS_out)
        oilspill_wrapper.run_mul_GNOME(311584.1,3113650.2,starttime,endtime,period,900, opt='ROMS')
        oilspill_wrapper.GNOME_GM_visualization(opt='ROMS')
        oilspill_wrapper.GNOME_GE_animation(13,starttime,opt='ROMS')
    
    elif mode==2:
    
        infile=dire+'/CoarseTri/rundata/'+'GalvCoarse_0000.nc'
        outfile=dire+'/GNOME/'+'txsuntans.nc'
        hydro_wrapper.runSUNTANS(starttime, endtime)
        oilspill_wrapper.init_model()
        oilspill_wrapper.Tx_SUNTANS(infile,outfile) ##input: SUNTANS input; output: GNOME input
        oilspill_wrapper.run_mul_GNOME(321947.94,3256260.05,starttime,endtime,period,900)
        oilspill_wrapper.GNOME_GM_visualization()
        oilspill_wrapper.GNOME_GE_animation(13,starttime)
        print 'end'
        
    elif mode==3:
        
        print "under developing!!\n"
        ####use the ROMS current data to run GNOME####
        downloadROMS(starttime,endtime)
        ROMS_file=dire+'/DATA/'+'txla_subset_HIS.nc'
        ROMS_out=dire+'/GNOME/'+'hiroms_ss_rho.nc'
        oilspill_wrapper.init_model(opt='ROMS') #The option is 'ROMS', since will run GNOME using ROMS output first
        oilspill_wrapper.HIROMS(ROMS_file,ROMS_out)                
        #oilspill_wrapper.run_mul_GNOME(311584.1,3113650.2,starttime,endtime,period,900, opt='ROMS')
        #(utm_x,utm_y)=utm.from_latlon(28.600572, -94.728385)[0:2]
        (utm_x,utm_y)=utm.from_latlon(28.353786, -95.315109)[0:2]
        oilspill_wrapper.run_mul_GNOME(utm_x,utm_y,starttime,endtime,period,900, opt='ROMS')
        oilspill_wrapper.GNOME_GM_visualization(opt='ROMS')
        ####duplicate the particles####
        oilspill_wrapper.duplicate(starttime,endtime,900)
        oilspill_wrapper.visualization2()
        
        
                                        
    else:
        print "There is no such mode, check HyosPy input mode!!!\n"
        
        
        

    
        
        
starttime='2014-08-01'
endtime='2014-08-02'
HyosPy(starttime, endtime, 24, mode=2)
