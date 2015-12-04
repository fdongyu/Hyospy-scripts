# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 11:16:08 2015

@author: dongyu
"""


import os
import subprocess
#subprocess.call(["chmod", "+x", "bashrc"])
#subprocess.Popen('./bashrc', shell=True )

from DownloadTool import downloadUSGSriver, downloadTCOONTide, updateWind, downloadROMS
from FileTool import write2db, increaseT1, increaseT2, increaseT3, convertFormat, reviseDatFile
from suntans_driver import generateBI
import pdb


def runSUNTANS(starttime,endtime):
    """
    Specify the time for a model run and import envir var
    for example: starttime='2014-08-21', endtime='2014-08-23'
    """
    
    """
    specify the time period to download the data. 
    Note that this period should be longer and include the run period
    For example the start='2014-12-01' and end='2014-12-10' 
    This can be done by employing the function increaseT()
    """
      
    (start,end)=increaseT1(starttime,endtime)
    (romst1,romst2)=increaseT3(starttime,endtime)
    downloadROMS(romst1,romst2)    #download initial condition data  
    
    downloadUSGSriver(start,end)              	#download the river inflow data
    staid='022' 					#specify the station id to download tide data 
    downloadTCOONTide(start,end,staid) 			#download the tide
    #Update and write the resulting data into the database file
    dbfile = 'GalvestonObs.db'
    write2db(dbfile)
    
    """
    download the wind data from TCOON and insert the data into the existing file
    Galveston_NARR_20122016 in folder DATA/ which can be read by SUNTANS
    """
    windID='207' #specify the wind data location
    updateWind(windID)
    """
    run the shell script: buildFolder to build the rundata folder 
    and copy the necessary files to start a new run
    """

    subprocess.Popen('./CoarseTri/buildFolder',shell=False)

    '''
    generate the boundary condition and initial condition, 
    the starttime and endtime should be in the format that
    start = '20141202.000000'
    end = '20141206.000000'
    '''

    (t1,t2)=increaseT2(starttime,endtime)
    t1=convertFormat(t1)+'.000000'
    t2=convertFormat(t2)+'.000000'
    generateBI(t1,t2)
    
    '''
    revise the file suntans.dat to specify the timeperiod and starttime
    Note that the time format for starttime and endtime is like: 
    starttime='2014-12-03', endtime='2014-12-08'
    which is also the run time
    '''
    reviseDatFile(starttime,endtime)

    ##########################run the model############################
    dire=os.getcwd()
    os.chdir(dire+'/CoarseTri')
    #subprocess.call('make test &> suntans.out & ',shell=True)
    subprocess.call('make test ',shell=True)
    
    os.chdir(dire)
    
###########################################test only#################################################
#starttime='2014-08-21'
#endtime='2014-08-23'
#runSUNTANS(starttime, endtime)
