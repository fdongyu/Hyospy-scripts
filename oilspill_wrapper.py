# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 11:34:40 2015

@author: dongyu
"""

import os
import numpy as np
import subprocess
import shutil
import math
from netCDF4 import Dataset
from datetime import datetime, timedelta
from urllib2 import HTTPError
import utm
import FileTool
import nctools
import hydro_wrapper
from gnome import scripting
from gnome.basic_types import datetime_value_2d
from gnome.utilities.remote_data import get_datafile
from gnome.model import Model
from gnome.map import MapFromBNA
from gnome.environment import Wind
from gnome.spill import point_line_release_spill
from gnome.movers import RandomMover, WindMover, GridCurrentMover
from gnome.outputters import Renderer
from gnome.outputters import NetCDFOutput
import pdb



def Tx_SUNTANS(infile,outfile):
    """
    step 1: 
    interpolate the data from SUNTANS output to GNOME input
    Input:
        infile='GalvCoarse_0000.nc'
        outfile='txsuntans_example.nc'
    """
    '''
    Sample script to retrieve data from unstructured grid netcdf "file" (can be
    OPeNDAP url), generate necessary grid topology (boundary info), and write 
    GNOME compatible output.
    The boundary file is saved to the data files directory so it only needs 
    to be generated once (unless you are subsetting the grid).
    To process multiple files (urls) either
    a) pass the filenames/urls in as a list -- this creates a netcdf4 MFDataset and is
    a good option for not too many files (all output is written to one nc file for GNOME 
    in this case)
    b) add a file list loop -- in this case put it after the grid topo vars are loaded (as
    this only has to be done once). See NGOFS_multifile_example.py
    '''
    
    data_files_dir=os.path.dirname(__file__)
       
    # specify local file or opendap url 
    data_file = os.path.join(data_files_dir,infile)
    
    # the utools class requires a mapping of specific model variable names (values)
    # to common names (keys) so that the class methods can work with FVCOM, SELFE,
    # and ADCIRC which have different variable names
    # (This seemed easier than finding them by CF long_names etc)
    #!!!!!!!!txsuntans output on server does not include eles_surrounding_ele info
    #I have it saved as a netcdf file included in libgoods data_files directory
    var_map = { 
                'time':'time',\
                'u_velocity':'uc', \
                'v_velocity':'vc', \
                'nodes_surrounding_ele':'cells',\
                'eles_surrounding_ele':'nbe',\
                'edge_node_connectivity':'edges',\
                }  

    # class instantiation creates a netCDF Dataset object as an attribute
    txsuntans = nctools.ugrid(data_file)

    # get longitude, latitude, and time variables
    print 'Downloading data dimensions'
    txsuntans.get_dimensions(var_map)
    
    # UTM coordinates -- calculate lat/lon
    x = txsuntans.Dataset.variables['xp'][:]
    y = txsuntans.Dataset.variables['yp'][:]
    lon = np.ones_like(x); lat = np.ones_like(x)
    for ii in range(len(x)):
        lat[ii], lon[ii] = nctools.utmToLatLng(15,x[ii],y[ii])
    txsuntans.data['lon'] = lon
    txsuntans.data['lat'] = lat
    txsuntans.atts['lon'] = {'long_name': 'longitude'}
    txsuntans.atts['lat'] = {'long_name': 'latitude'}
    
    # get grid topo variables (nbe, nv)
    print 'Downloading grid topo variables'
    txsuntans.get_grid_topo(var_map)

    edge_types = txsuntans.Dataset.variables['mark'][:].tolist()
    #"0 - computational; 1 - closed; 2 flux BC; 3 - stage BC; 4 - other BC; 5 - interproc; 6 - ghost
    #Based on personal communication we use 1,2, and 3
    bound_id, bound_type = zip(*[(i,x) for i,x in enumerate(edge_types) if x>0 and x<4])
    bound_segs = txsuntans.data['edges'][bound_id,:] + 1 #Node numbering starts at 0, GNOME expects it to start at 1 -- adjust nv and bnd
    bound_type_gnome = (np.array(bound_type)<=2).choose(1,0)
    
    txsuntans.order_boundary(bound_segs.tolist(),list(bound_type_gnome))
    txsuntans.atts['bnd'] = {'long_name':'Boundary segment information required for GNOME model'} 
    
    #Node numbering starts at 0, GNOME expects it to start at 1 -- adjust nv and bnd
    txsuntans.data['nv'] = txsuntans.data['nv'] + 1
    
    ## GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
    txsuntans.atts['nbe']['order'] = 'ccw'

    '''
    !!!!!!!!!!!!!All the stuff above here only has to be done once -- if you want to 
    process multiple files, I'd put a loop here and just keep overwriting 
    txsuntans.data['u'] and ['v'] and incrementing the output file name
    Also need to change txsuntans.data['time'] appropriately
    '''
    # get the data
    print 'Loading u/v'
    txsuntans.get_data(var_map,zindex=0) 
    
    
    print 'Writing to GNOME file'
    txsuntans.write_unstruc_grid(os.path.join(data_files_dir, outfile))
    
    txsuntans.Dataset.close()
    

def HIROMS(infile,outfile,subset=1):
    
    '''
    Sample script to retrieve data from Arakawa c-grid type model
    '''     
    
    data_files_dir=os.path.dirname(__file__)
        
    #url = 'http://barataria.tamu.edu:8080/thredds/dodsC/txla/roms/2014/ocean_his_03.nc'
    #url='txla_subset_HIS_Aug2014.nc'
    #url=infile
    url = os.path.join(data_files_dir,infile)
     
    
    #HI ROMS output does not include psi grid -- create from rho grid
    var_map = { 'time':'ocean_time',
               }  
    hiroms = nctools.roms(url)
    hiroms.get_dimensions(var_map)
    
    # try subset
    #subset = 1
    
    #Only download last five timesteps
    #ti=[len(hiroms.data['time'])-5,len(hiroms.data['time']),1]
    ti=[0,len(hiroms.data['time']),1]
    #pdb.set_trace()
    
    if subset:
        # this case interpolates to rho grid and reduces lon/lat size 
        # (compatible with GNOME at present)
        hiroms.subset([28.17,-95.53,30.0,-94.25],lat='lat_psi',lon='lon_psi')
        hiroms.get_grid_info(xindex=hiroms.x,yindex=hiroms.y)
        hiroms.get_data(var_map,tindex=ti,xindex=hiroms.x,yindex=hiroms.y)
        hiroms.reduce_latlon_mesh_for_GNOME()
        #ofn = os.path.join(data_files_dir,'hiroms_ss_rho_reduced.nc')
        ofn = os.path.join(data_files_dir,outfile)
        
        hiroms.data['lon_ss'] = hiroms.data['lon_psi_ss']
        hiroms.data['lat_ss'] = hiroms.data['lat_psi_ss']
        hiroms.write_nc(ofn,is3d=False)
    else:
        #u/v interpolated to rho grid, u/v and lat/lon the same size (works in current GNOME)
        print 'interp and reduce'
        hiroms.get_grid_info()
        hiroms.get_data(var_map,tindex=ti)
        hiroms.reduce_latlon_mesh_for_GNOME()
        ofn = os.path.join(data_files_dir,'hiroms_rho_reduced.nc')
        hiroms.data['lon'] = hiroms.data['lon_psi']
        hiroms.data['lat'] = hiroms.data['lat_psi']
        hiroms.write_nc(ofn,is3d=False)
        
        #u/v interpolated to rho grid, lat/lon on psi grid larger than u/v
        print 'interp only'
        hiroms.get_dimensions(var_map)
        hiroms.get_grid_info()
        hiroms.get_data(var_map,tindex=ti)
        ofn = os.path.join(data_files_dir,'hiroms_rho.nc')
        hiroms.data['lon'] = hiroms.data['lon_psi']
        hiroms.data['lat'] = hiroms.data['lat_psi']
        hiroms.write_nc(ofn,is3d=False)
    
        #u/v on native grids
        print 'native'
        hiroms.get_dimensions(var_map)
        hiroms.get_grid_info()
        hiroms.get_data(var_map,tindex=ti,interp=False)
        ofn = os.path.join(data_files_dir,'hiroms_native.nc')
        hiroms.write_nc_native(ofn,is3d=False)
    

    
def init_model(opt='SUNTANS'):
    """
    This function is used to move the necessary files for a GNOME run into the GNOME folder
    gird file: coast.bna
    KML file:  GE_animation.txt
    """
    
    data_dir=os.getcwd()+'/DATA/'
    GNOME_dir=os.getcwd()+'/GNOME/' 
    
    if os.path.exists(GNOME_dir):
        shutil.rmtree(GNOME_dir)
    os.makedirs(GNOME_dir)
    
    filelist=['coast.bna', 'GE_animation.txt', 'javascript.txt']
    try:
        for ff in filelist:
            FileTool.move_file(data_dir+ff, GNOME_dir)
    except:
        raise Exception, 'No target directory called:\n%s or necessary file missing in folder /DATA/'%GNOME_dir
    
    FileTool.move_file(os.getcwd()+'/DATA/wind.WND',GNOME_dir)

    if opt=='SUNTANS':
        row1='NetCDF Files'
        row2='[File]    ./txsuntans.nc'
        f=open(GNOME_dir+'current_SUNTANS.txt','w')
        f.write(row1+'\n')
        f.write(row2+'\n')
        f.close()
    elif opt=='ROMS':
        row1='NetCDF Files'
        row2='[File]    ./hiroms_ss_rho.nc'
        f=open(GNOME_dir+'current_ROMS.txt','w')
        f.write(row1+'\n')
        f.write(row2+'\n')
        f.close()
    else:
        raise IOError('No such option, check input for init_model')
        
    
#def init_ROMS():
#    data_dir=os.getcwd()+'/DATA/'
#    ROMS_dir=os.getcwd()+'/ROMS_GNOME/'
#    filelist=['coast.bna', 'GE_animation.txt', 'javascript.txt']
#    FileTool.move_file(os.getcwd()+'/DATA/wind.WND',ROMS_dir)
#    try:
#        for ff in filelist:
#            FileTool.move_file(data_dir+ff, ROMS_dir)
#    except:
#        raise Exception, 'No target directory called:\n%s or necessary file missing in folder /DATA/'%ROMS_dir
#    row1='NetCDF Files'    
#    row2='[File]    ./hiroms_ss_rho.nc'
#    f=open(ROMS_dir+'current.txt','w')
#    f.write(row1+'\n')
#    f.write(row2+'\n')
#    f.close()
#    
#def init_ROMS_SUN():
#    
#    data_dir=os.getcwd()+'/DATA/'
#    RS_dir=os.getcwd()+'/ROMS_SUN_GNOME/'
#    if not os.path.exists(RS_dir):
#        os.makedirs(RS_dir)
#    filelist=['coast.bna', 'GE_animation.txt', 'javascript.txt']
#    FileTool.move_file(os.getcwd()+'/DATA/wind.WND',RS_dir)
#    try:
#        for ff in filelist:
#            FileTool.move_file(data_dir+ff, RS_dir)
#    except:
#        raise Exception, 'No target directory called:\n%s or necessary file missing in folder /DATA/'%RS_dir
#    row1='NetCDF Files'    
#    row2='[File]    ./hiroms_ss_rho.nc'
#    f=open(RS_dir+'current.txt','w')
#    f.write(row1+'\n')
#    f.write(row2+'\n')
#    f.close()
    
        
def make_model(coor, yr, month, day, period=46, dt=900 ,opt='SUNTANS', images_dir=os.getcwd()+"/images"):
    '''
    Run multiple GNOME;
    yr,month,day---oil spill start time;
    period---oil spill simulation duration;
    dt--oil spill time step in second
    '''    
    
    print "initializing the model"

    base_dir=os.getcwd()
    #start_time = datetime(2014,8,21,0)
    start_time = datetime(yr,month,day,0)
    model = Model(start_time = start_time,
                              duration = timedelta(hours=period),
                              time_step =dt,
                              uncertain = False,
                              )
    
    mapfile = os.path.join(base_dir, './coast.bna')
    print "adding the map"
    gnome_map = MapFromBNA(mapfile, refloat_halflife=6)  # hours
    
    print "adding renderer" 
    model.outputters += Renderer(mapfile, images_dir, size=(1800, 1600))

    
    print "adding a wind mover from a time-series"
    ## this is wind
    wind_file=get_datafile(os.path.join(base_dir, 'wind.WND'))
    wind = Wind(filename=wind_file)
    w_mover = WindMover(wind)
    model.movers += w_mover
    
    print "adding a current mover:"
    ## this is currents
    curr_file = get_datafile(os.path.join(base_dir, 'current_'+opt+'.txt'))
    model.movers += GridCurrentMover(curr_file)

    #model.movers += RandomMover(1000)
    ##
    ## Add some spills (sources of elements)
    ##
    print "adding 13 points in a cluster that has some small initial separation as the source of spill"
    
    for i in range(len(coor)):
        
        xcoor,ycoor=nctools.utmToLatLng(15,coor[i][0],coor[i][1],northernHemisphere=True)
        model.spills += point_line_release_spill(num_elements=1,
                                                start_position = (ycoor,xcoor, 0.0),
                                                release_time = start_time,
                                                )

    print "adding netcdf output"
    netcdf_output_file = os.path.join(base_dir,'GNOME_'+opt+'.nc')
    scripting.remove_netcdf(netcdf_output_file)
    model.outputters += NetCDFOutput(netcdf_output_file, which_data='all')
    
    return model
        
def run_mul_GNOME(x,y,starttime,endtime,period,dt,opt='SUNTANS'):
    
    '''
    number--set up how many GNOME to run
    Run multiple GNOME;
    x,y--UTM coordinates for oil spill location;
    starttime='2014-08-21'
    endtime='2014-08-22'
    yr,month,day---oil spill start time;
    period---oil spill simulation duration;
    dt--oil spill time step in second
    opt--'SUNTANS' or 'ROMS' or 'both'
    '''
    yr=int(starttime[0:4]);
    month=int(starttime[5:7].lstrip("0").replace("0", " "))
    day=int(starttime[8:].lstrip("0").replace("0", " "))
    timeformat='%Y-%m-%d'
    diff=datetime.strptime(endtime,timeformat)-datetime.strptime(starttime,timeformat)

    if diff.days*24 < period:
        raise IOError('Input period is longer than the data length')
            
    utm_x=x ; utm_y=y
    coor=[[utm_x,utm_y],[utm_x-10,utm_y],[utm_x-20,utm_y],[utm_x+10,utm_y],[utm_x+20,utm_y],[utm_x,utm_y-10],\
     [utm_x,utm_y-20],[utm_x,utm_y+10],[utm_x,utm_y+20],[utm_x-10,utm_y+10],[utm_x-10,utm_y-10],[utm_x+10,utm_y+10],[utm_x+10,utm_y-10]]
    
    work_dir=os.getcwd()
    
    os.chdir(work_dir+'/GNOME/')
            
    scripting.make_images_dir()
    model = make_model(coor,yr,month,day,period,dt,opt,images_dir=os.getcwd()+"/images")
    wdd=[0.01,0.04]
    model.full_run(wdd,logger=True)
    print 'end'
    os.chdir(work_dir) 
     
             
    
def GNOME_GM_visualization(opt='SUNTANS'):
    
    """
    Visualize different GNOME oil spill trajectories on 2D Google Map GIS;
    number --- set up how many set of GNOME outputs to visualize
    """
    
    
    locations=[]

    print "Generating different oil spill tracks (GNOME) on Google Map"
    
    
    base_dir=os.getcwd()+'/GNOME/'
    
        
    a=Dataset(base_dir+'GNOME_'+opt+'.nc','r')
    
    for j in range(len(a.variables[u'longitude'][:])):
        
        locations.append([a.variables[u'latitude'][j],a.variables[u'longitude'][j]])
    
    if opt=='ROMS' or opt=='both':
        f=open(base_dir+'/javascript.txt','r')
        content=f.readlines()
        os.remove(base_dir+'/javascript.txt')
        open(base_dir+'/javascript.txt','w').write('')
        h=open(base_dir+'/javascript.txt','w')
        for line in content:
            h.write(line.replace('zoom: 12','zoom: 9').replace('(29.51, -94.84)','(29.031161, -94.574199)'))
        h.close()
    
    file1 = open(base_dir+'/javascript.txt','r')     # open the javascript
    row=[]
    for s in file1.readlines():        
        row.append(s)

    row[11]='    '+'var'+' '+'locations'+'='+str(locations) + '\n'  
    
    h=open(base_dir+'/GNOME_Google_map.html','w')
    for l in row:
        g=''.join([str(j) for j in l])     
        h.write(g+'\n')
    h.close
           


def GNOME_GE_animation(np,starttime,opt='SUNTANS'):
    
    '''
    Animate different GNOME oil spill trajectories on 3D Google Earth GIS;
    number --- set up how many set of GNOME outputs to visualize;
    np is the total particle number;
    starttime='2014-08-21'
    yr,month,day---oil spill start time
    '''
    yr=int(starttime[0:4]);
    month=int(starttime[5:7].lstrip("0").replace("0", " "))
    day=int(starttime[8:].lstrip("0").replace("0", " "))
    

    base_dir=os.getcwd()+'/GNOME/'
    
    a=Dataset(base_dir+'GNOME_'+opt+'.nc','r')
        
    file1 = open(base_dir+'/GE_animation.txt','r')     
    row=[]
    for s in file1.readlines():     
        row.append(s)
        
    nt=len(a.variables[u'time'][:])   #  nt is the total time steps

    print "Generating multiple oil spill tracks (GNOME) on 3D Google Earth"
    
    for n in range(nt):   
    
        for j in range(np*n,np*(n+1)):                       
            
            #row.append('  <Placemark>'+'\n'+'    <TimeStamp>'+'\n'+'      <when>'+str(yr)+'-'+str('%02d' % month)+'-'+str('%02d' % (day+n/96))+'T'+ str("%02d" % ((n/4)%24))+':' +str("%02d" %((n*15)%60))+':00Z</when>'+'\n'+'    </TimeStamp>'+'\n'+'    <styleUrl>#</styleUrl>'+'\n'+'    <Point>'+'\n'+'      <coordinates>'+str(a.variables[u'longitude'][j])+','+str(a.variables[u'latitude'][j])+'</coordinates>'+'\n'+'    </Point>'+'\n'+'  </Placemark>'+'\n')
            row.append('  <Placemark>'+'\n'+'    <TimeStamp>'+'\n'+'      <when>'+str(yr)+'-'+str('%02d' % month)+'-'+str('%02d' % (day+n/96))+'T'+ str("%02d" % ((n/4)%24))+':' +str("%02d" %((n*15)%60))+':00Z</when>'+'\n'+'    </TimeStamp>'+'\n'+'    <styleUrl>#'+str(1)+'</styleUrl>'+'\n'+'    <Point>'+'\n'+'      <coordinates>'+str(a.variables[u'longitude'][j])+','+str(a.variables[u'latitude'][j])+'</coordinates>'+'\n'+'    </Point>'+'\n'+'  </Placemark>'+'\n')    
    row.append('</Document>'+'\n'+'</kml>')    
      
    h=open(base_dir+'/GNOME_GE.kml','w')
    for l in row:
        g=''.join([str(j) for j in l])
        h.write(g+'\n')
    h.close         
    
def spill_locations(infile):
    """
    input:
        infile--GNOME output file 
        example: 'GNOME_output.nc'
    output:
        loc
        oil spill particle locations
        loc is divided into n sub-list, n is the number of oil spill particles
        in each sub-list is the trajectory of the specified particle
    """
    base_dir = os.path.dirname(__file__)
    filename=base_dir+'/GNOME/'+infile
    a=Dataset(filename,'r')
    print "Writing the particle points into the locations"
    
    locations=[]
    for j in range(len(a.variables[u'longitude'][:])):        
        locations.append([a.variables[u'latitude'][j],a.variables[u'longitude'][j]])
    
    tlen=len(a.variables['time'])
    td=len(a.variables['spill_num'])
    nn=td/tlen
    
    group=[]    
    for ii in range(tlen):
        group.append(locations[ii*nn:(ii+1)*nn])      
    
    one_track=[]
    for jj in range(nn):
        pp=[]
        for ii in range(tlen):
            pp.append(group[ii][jj])
        one_track.append(pp)

    loc=[]
    for i in range(len(one_track)):
        loc.append(one_track[i])
    
    return loc
    
        
def duplicate(starttime,endtime,dt):
    """
    
    Traverse the SUNTANS boundary grid to locate the points that ROMS result enters
    and run GNOME
    
    """
    print "Trying to duplicate the oil spill particles at the SUNTANS boundary!!\n"
    base_dir = os.path.dirname(__file__)
    
    try:
        ncfile=Dataset(base_dir+'/CoarseTri/rundata/'+'GalvCoarse_0000.nc','r')
    except IOError:
        print 'No SUNTANS output available'
        
    mark=ncfile.variables['mark'][:]
    xe=ncfile.variables['xe'][:]
    ye=ncfile.variables['ye'][:]

    tide=[] #tide boundary mark
    river=[] #river boundary mark
    for ii in mark:
        if ii==3:
            tide.append(ii)
        if ii==2:
            river.append(ii)

    tide_index=[]
    for ind,mk in enumerate(mark):
        if mk==3:
            tide_index.append(ind)
    ###########Output the cooridinate of SUNTANS open boundary points#########
    open_boundary=[]
    for kk in tide_index:
        open_boundary.append([xe[kk],ye[kk]])
    boundary=[]
    for i in range(len(open_boundary)):
        boundary.append(nctools.utmToLatLng(15,open_boundary[i][0],open_boundary[i][1]))
        
    ###########ROMS output############
    roms_spill=spill_locations('GNOME_ROMS.nc')
    ##################################
    cbp=[]            ###cbp is cross SUNTANS boundary point 
    SUN13=[]
    for nn in range(len(roms_spill)):
        particle=roms_spill[nn]
        for i in range(len(particle)):
            find=False
            for j in range(len(boundary)):                                
                if abs(boundary[j][0]-particle[i][0])<0.0035 \
                and abs(boundary[j][1]-particle[i][1])<0.0035:
                    SUN13.append([boundary[j][0],boundary[j][1]])
                    find=True
                    break
            if find==True:
                cbp.append(i) ###cbp is cross SUNTANS boundary point 
                break     
         
    if SUN13==[]:
        print '#########\n ROMS output particles will not cross SUNTANS boundary!! \n#########'
        print '#########\n check the Google map for details!! \n#########'
        raise RuntimeError('adjust initial locations and time period for GNOME run \
                            or you can change to mode 1 and rerun')
    else:
        print '#########\n %s particles generating at SUNTANS boundary!! \n#########'%str(len(SUN13))
        
    timeformat='%Y-%m-%d'
    SUN_starttime=datetime.strptime(starttime,timeformat)+timedelta(math.floor(min(cbp)/(24*3600./dt))+1)
    SUN_starttime=SUN_starttime.strftime(timeformat) ##SUNTANS starttime before particle cross ROMS boundary
    diff=datetime.strptime(endtime,timeformat)-datetime.strptime(SUN_starttime,timeformat) 
    period=diff.days*24

    ####running SUNTANS from the SUN_starttime to endtime 
    hydro_wrapper.runSUNTANS(SUN_starttime, endtime) 
    ####generating input of GNOME from SUNTANS output
    infile=base_dir+'/CoarseTri/rundata/'+'GalvCoarse_0000.nc'
    outfile=base_dir+'/GNOME/'+'txsuntans.nc'
    Tx_SUNTANS(infile,outfile)
    
    coor=[]
    for i in range(len(SUN13)):
        coor.append(utm.from_latlon(SUN13[i][0],SUN13[i][1])[0:2])
    
    
    yr=int(SUN_starttime[0:4]);
    month=int(SUN_starttime[5:7].lstrip("0").replace("0", " "))
    day=int(SUN_starttime[8:].lstrip("0").replace("0", " "))
        
    row1='NetCDF Files'
    row2='[File]    ./txsuntans.nc'
    f=open(os.getcwd()+'/GNOME/'+'current_SUNTANS.txt','w')
    f.write(row1+'\n')
    f.write(row2+'\n')
    f.close()
    
    work_dir=os.getcwd()
    os.chdir(work_dir+'/GNOME/')
    scripting.make_images_dir()         
    model = make_model(coor,yr,month,day,period,dt,images_dir=os.getcwd()+"/images")
    wdd=[0.01,0.04]
    model.full_run(wdd,logger=True)
    os.chdir(work_dir) 
    
def visualization2():
    
    '''
    Visualize different GNOME oil spill trajectories on 2D Google Map;
    number --- set up how many set of GNOME outputs to visualize
    '''
    base_dir=os.getcwd()+'/GNOME/'
        
    locations=[]

    print "Generating different oil spill tracks (GNOME) on Google Map"
    filelist=['GNOME_ROMS.nc','GNOME_SUNTANS.nc']
    #filelist=['GNOME_SUNTANS.nc']
    
    for ii in range(len(filelist)):
        filename=filelist[ii]
    
        a=Dataset(base_dir+filename,'r')
    
        for j in range(len(a.variables[u'longitude'][:])):
        
            locations.append([a.variables[u'latitude'][j],a.variables[u'longitude'][j]])
            
        
    file1 = open(base_dir+'/javascript.txt','r')     # open the javascript
    row=[]
    for s in file1.readlines():        
        row.append(s)

    row[11]='    '+'var'+' '+'locations'+'='+str(locations)   
    
    h=open(base_dir+'/GNOME_Google_map2.html','w')
    for l in row:
        g=''.join([str(j) for j in l])
        h.write(g+'\n')
    h.close        
    
    
    
    
    

    


    
    
    


    
