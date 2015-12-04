# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 16:05:51 2015

@author: dongyu
"""
import urllib2
import string
import datetime
import numpy as np
import os
import netcdfio
import urllib
import shutil
import tarfile
from contextlib import closing
import math
from netCDF4 import Dataset
import time as tim
from othertime import MinutesSince
from romsio import roms_subset
from FileTool import copy_file
import pdb




def downloadUSGSriver(starttime,endtime):
    from getUSGSnwis import getUSGSnwis
    stationids = ['08066500',\
                '08078000',\
                '08067500',\
                '08073600',\
                '08042558',\
                '08031000',\
                '08067525']
          
    #starttime = '2014-11-01'
    #endtime = '2014-11-10'
    ncfile = 'USGS_Rivers.nc'
    getUSGSnwis(stationids,starttime,endtime,ncfile)    
    files=['USGS_Rivers.dbf','USGS_Rivers.nc','USGS_Rivers.shp','USGS_Rivers.shx']
    dst_dir=os.getcwd()+'/DATA/'
    for src_file in files:
        copy_file(src_file, dst_dir)


def downloadTCOONTide(starttime,endtime,staid):
    data=downloadTide(starttime,endtime,staid)
    dire=os.path.dirname(__file__)
    dire=os.getcwd()
    ncfile = dire+'/DATA/TCOONTide.nc'
    globalatts = {'title':'TCOON oceanographic observation data'}
    netcdfio.writePointData2Netcdf(ncfile,data,globalatts)
    print 'end writing tide data...'
    

def downloadTide(starttime,endtime,staid):
    """ Main function to construct data in a uniform format"""
    lon=-94.7933
    lat=29.31
    ID=8771450
    nn='Galveston Pier 21'
    vv='waterlevel'
    # Build up the output data as a list of dictionaries
    meta=[]
    coords = [{'Name':'longitude','Value':lon,'units':'degrees East'},\
        {'Name':'latitude','Value':lat,'units':'degrees North'},\
        {'Name':'time','Value':[],'units':'minutes since 1970-01-01 00:00:00'}]
    attribs = {'StationID':str(ID),'StationName':nn,'Data':[],'coordinates':'time, longitude, latitude','coords':coords} 
    meta.append(attribs)

    data=[]
    tmp = {vv:meta[0]}
    ctr=0

    #staid='014' #TCOON station ID
    output,t,atts = getTCOON(starttime,endtime,staid)
    if ctr==0:                
        tmp[vv].update(atts)
    # Append the data to the list array
    tmp[vv]['Data'] += output
    # Append the time data
    ctr=-1
    for cc in tmp[vv]['coords']:
        ctr+=1
        if cc['Name']=='time':
            tmp[vv]['coords'][ctr]['Value'] += t 
    if np.size(tmp[vv]['Data']) > 0:
            data.append(tmp)
            
    return data

def getTCOON(starttime,endtime,staid,choice=True):
    startt=starttime[5:7]+'.'+starttime[8:10]+'.'+starttime[0:4]
    endt=endtime[5:7]+'.'+endtime[8:10]+'.'+endtime[0:4]
    url = 'http://lighthouse.tamucc.edu/pd?stnlist='+staid+'&serlist=pwl%2Charmwl&when='+startt+'-'+endt+'&whentz=UTC0&-action=c&unit=metric&elev=msl'
    try:
        print 'Opening: %s'%url
        f=urllib2.urlopen(url)
    except:
        raise Exception, 'cannot open url:\n%s'%url
    data1=[]
    data0=[] #Note data0 is different from data that data0 includes the possible string 'NA'
    data=[] #To store the tidal elevation and time
    attribs = {'long_name':'Water surface elevation','units':'m'}
    for s in f:
        if '#' not in s:
            line0=s.split()
            data0.append(line0)
            if 'NA' not in s:
                line=s.split()
                data1.append(line)
    if choice:
        # calculate the mean difference between pwl and harmwl with N = 30 (3 hours)        
        sum=0
        #pdb.set_trace()
        for n in range(len(data1),len(data1)-30,-1):
            if data1[n-1][1]== 'RM': data1[n-1][1] = data1[n-1][2]
            aa=string.atof(data1[n-1][1]);bb=string.atof(data1[n-1][2])
            diff=aa-bb
            sum=sum+diff
        x=sum/len(data1)
        # find the index of the first 'NA'
        column2=[]
        for line in data0:
            column2.append(line[1])
        #pdb.set_trace()
        if 'NA' in column2:
            y=column2.index('NA')
        else:
            y=len(data0)
        # add x with harmwl to creat the forecast part of pwl, set the return interval as 72 hours (72*60/6=720)
        for ii in range(len(data0)):
            if data0[ii][1]=='RM':
                data0[ii][1]=data0[ii][2]
        for i in range(y,len(data0)):
            data0[i][1]=string.atof(data0[i][2])+x*(1-i/720)
        for k in range(y):
            data0[k][1]=string.atof(data0[k][1])
        for m in range(len(data0)):
            del data0[m][2]
        for line in data0:
            data.append(line)
    else:
        for line in data0:
            if line[1]=='NA':
                line[1]=0
            del line[2]
            data.append(line)
    tt=[] #time
    for kk in range(len(data)):
        outTime=parseTime(data[kk][0])
        tt.append(outTime)
    elev=[] #water elevation
    for ii in range(len(data)):
        tem=data[ii][1]
        elev.append(tem)
    
    print 'finish downloading data from TCOON...'
#    pdb.set_trace()
    return elev, tt, attribs
    
 
def parseTime(inTime):
    basetime=datetime.datetime(2014,1,1,00,00)
    yy=(int(inTime[0:4])-2014)*365
    dd=inTime[4:7]
    hh=inTime[8:10]
    mm=inTime[10:]
    delta=datetime.timedelta(days = yy+ int(dd)-1,hours=int(hh),minutes=int(mm))
    outTime=basetime+delta
    return MinutesSince(outTime)[0]  


####Below is the example    
#starttime='2014-03-01'
#endtime='2014-03-02'
#staid='014'
#downTCOONTide(start,end,staid)
##output,t,atts=getTCOON(starttime,endtime,staid)
##data=downloadTide(starttime,endtime)
##dire=os.path.dirname(__file__)
##ncfile = dire[0:56]+'DATA/TCOONTide.nc'
##globalatts = {'title':'TCOON oceanographic observation data'}
##netcdfio.writePointData2Netcdf(ncfile,data,globalatts)
##print 'end writing tide data...'

def updateWind(windID):
    '''
    This function updates the wind file, insert the downloaded TCOON wind data
    into the file and replace the original wind velocity
    '''
    tt,data = getWindData(windID) #tt is time series, data is wind data
    filename=os.getcwd()+'/DATA/Galveston_NARR_20122016'
    wind=Dataset(filename+'.nc','r')
    time=wind.variables['Time'][:]
    #the unit of time is 'minutes since 1970', next step is to change the time
    # format to something like "19901231.120000"
    otherTime=[]
    timeformat='%Y%m%d.%H%M%S'
    for ii in range(len(time)):
        timeArray=tim.localtime(time[ii])
        otherTime.append(tim.strftime(timeformat, timeArray))
    #Next is to adjust the time, which is increasing the time to real time
    a=datetime.datetime.strptime(otherTime[0],timeformat)
    b=datetime.datetime.strptime('20111231.120000',timeformat)
    diff=b-a
    inTime=[]
    for jj in range(len(otherTime)):
        counter=datetime.datetime.strptime(otherTime[jj],timeformat)+datetime.timedelta(diff.days)
        inTime.append(counter)
    print "The time period is from %s to %s"%(inTime[0],inTime[-1])
    # transform the resulting format (2011-12-31 12:00:00) to (20141101)
    oldTime=[]
    for kk in range(len(inTime)):
        counter=str(inTime[kk])[0:4]+str(inTime[kk])[5:7]+str(inTime[kk])[8:10]
        oldTime.append(counter)
    print 'The time period of the new format is %s to %s'%(oldTime[0],oldTime[-1])
    #determine the start index and the end index of gl.days in oldTime
    oldStart=oldTime.index(tt[0])
    oldEnd=oldTime.index(tt[-1])+8
    wind.close()
    nc=Dataset(filename+'.nc','a')
    UUwind=[]
    VVwind=[]
    for ii in range(len(data)):
        UUwind.append(data[ii][0])
        VVwind.append(data[ii][1])
    for pp in range(16):
        nc.variables['Uwind'][oldStart:oldEnd,pp]=UUwind
        nc.variables['Vwind'][oldStart:oldEnd,pp]=VVwind
    nc.close
    print "end updating wind data..."



def getWindData(windID,setGNOME=True):
    '''
    This function download the wind data from TCOON and generate the wind file for GNOME
    '''
    # Download wind forecast data from "http://seawater.tamu.edu/tglopu/twdb_lc.tar" and save it at "wind/twdb_1c.tar".
    dire=os.getcwd()+'/'
    
    if os.path.isdir(dire+'wind'):
        print 'wind directory exists, delete it and generate a new one'
        shutil.rmtree(dire+'wind')
            
        
    os.mkdir(dire+'wind')
    url = 'http://seawater.tamu.edu/tglopu/twdb_lc.tar'
    path = dire+'wind/twdb_1c.tar'
    try:
        print 'Opening: %s'%url
        data = urllib.urlopen(url).read()
    except:
        raise Exception, 'cannot open url:\n%s'%url
    f = file(path,'wb')
    f.write(data)
    f.close()
    # decide whether the download process is over and unzip the tar file to 'D:/wind' when the download is over.
    path = dire+'wind/twdb_1c.tar'
    if os.path.isfile(dire+'wind/twdb_1c.tar'):
       os.mkdir(dire+'wind/wind_data')
       with closing(tarfile.open(dire+'wind/twdb_1c.tar','r')) as t:
            t.extractall(dire+'wind/wind_data')
    
    data=[] # data to store wind data
    #choose the file of the target wind station    
    a=open(dire+'wind/wind_data/twdb'+windID+'.wndq', 'r').readlines()
    for s in a:
        if '*' not in s:
            if 'days' not in s:
                line=s.split()
                data.append(line)

    ######  create wind.WND for PyGNOME ######
    if setGNOME is True:
        row,col=np.shape(data)
        data0=np.zeros((row,col+1)).tolist() # to store the data for constructing GNOME wind file      
        #for i in range(len(data0)):
        #    data0[i].append(i)    
        for kk in range(len(data0)):
        ## wind direction
            data0[kk][6]=data[kk][5]
        ## wind magnitude
            data0[kk][5]=str(float(data[kk][4])*0.447)+','
        ## time
            data0[kk][4]='00,'
            data0[kk][3]=data[kk][3]+','
            data0[kk][1]=str('%01d' % int(data[kk][1]))+','
            stock=data[kk][0]
            data0[kk][0]=str('%01d' % int(data[kk][2]))+','
            data0[kk][2]=stock[-2:]+','  
        
        if os.path.isfile(dire+'GNOME/wind.WND'):
            print 'The GNOME wind data file already exists, generate new file!!!'
            os.remove(dire+'GNOME/wind.WND')
        else:
            print 'Generating GNOME wind file!!'
        
        ff=open(dire+'GNOME/wind.WND','w')

        ff.write('twdb'+windID+'\n')
        ff.write('27.798 -97.311'+'\n')
        ff.write('mps'+'\n')
        ff.write('LTime'+'\n')
        ff.write('0,0,0,0,0,0,0,0'+'\n')
        for i in (data0):
            k='   '.join([str(j) for j in i])
            ff.write(k+'\n')
        ff.close
    # make record of what time period of wind data are required
    time=[] #To store the time of wind data '20141201'
    for ii in range(len(data)):
        time.append(data[ii][0]+data[ii][1]+data[ii][2])
    # delete the time info in data
    for jj in range(len(data)):
         data[jj][0:4]=[]
    for nn in range(len(data)):
        aa=string.atof(data[nn][0]);bb=string.atof(data[nn][1])
        data[nn][0]=0.447*aa*math.sin(bb*math.pi/180);data[nn][1]=0.447*aa*math.cos(bb*math.pi/180)

    if os.path.isdir(dire+'wind'):
        shutil.rmtree(dire+'wind')

    return time, data

def downloadROMS(starttime,endtime):
    """
    This function is used to download the ROMS data from http://barataria.tamu.edu:8080/thredds/catalog.html
    The ROMS data is also subsetted to a certain area.
    The resulting data is in the format of txla_subset_HIS.nc, which is used to generate initial condition
    input:
        starttime = '2014-08-21'
        endtime   = '2014-08-22'
        timelims: ('20140801000000','20140831000000')        
    """
    timelims = (starttime.replace("-", "")+'000000', endtime.replace("-", "")+'000000')
       
    grdfile = 'http://barataria.tamu.edu:8080/thredds/dodsC/txla_nesting6_grid/txla_grd_v4_new.nc'
    #grdfile = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    
    #######specify the ROMS file on the server: http://barataria.tamu.edu:8080/thredds/ #######
    #ncfiles=['http://barataria.tamu.edu:8080/thredds/dodsC/txla/roms/2014/ocean_his_'+starttime[5:7]+'.nc']
    #ncfiles=['http://barataria.tamu.edu:8080/thredds/dodsC/oofv2/out/files/2014/ocean_his_0016.nc',\
    #		'http://barataria.tamu.edu:8080/thredds/dodsC/oofv2/out/files/2014/ocean_his_0017.nc']
    
    def _int2str(num):
            """
            convert the month format
            input a integer;
            output a string;
            """
            if num<10:
                return '0%s'%str(num)
            else:
                return '%s'%str(num)
        
    basedir='http://barataria.tamu.edu:8080/thredds/dodsC/oofv2/out/files'
    ncfiles=[]
    start=2*string.atoi(starttime[5:7])-1
    end=2*string.atoi(endtime[5:7])+1
    if starttime[0:4]==endtime[0:4]:
        print "beginning downloading ROMS data in %s"%starttime[0:4]
    else:
        raise IOError('Error downloading ROMS from the server, starttime endtime not in the same year')
        
    for i in range(start,end+1):
        filestr='%s/%s/ocean_his_00%s.nc'%(basedir,starttime[0:4],_int2str(i))
        ncfiles.append(filestr)
        
    #pdb.set_trace()
    #timelims = ('20140801000000','20140831000000')
    bbox = [-95.53,-94.25,28.17,30.0]
    roms = roms_subset(ncfiles,bbox,timelims,gridfile=grdfile)
    outfile = os.getcwd()+'/DATA/txla_subset_HIS.nc'
    roms.Writefile(outfile)
    roms.Go()
