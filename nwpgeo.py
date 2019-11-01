# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 09:25:19 2019
Numerical Weather Prediction Reader and Processing for ICON Model and COSMO Model
ICON BASIC INFORMATION
Die Vorhersagen des Regionalmodell ICON-EU werden aus den vier Modellläufen um 00, 06, 12, und 18 UTC 
bis +120 Stunden bereitgestellt und aus den Modelläufen um 03, 09, 15 und 21 UTC bis +30 Stunden. 
Für den Vorhersagezeitraum bis + 78 Stunden sind einstündige Zeitschritte verfügbar, 
von +81 Stunden bis +120 Stunden dreistündige Zeitschritte

COSMO BASIC INFORMATION:
    
Aufgrund der Zielstellung des COSMO-D2 (und des kleinen Modellgebiets)
sind auch nur relativ kurze Vorhersagezeiten sinnvoll. 
Vom COSMO-D2 werden alle 3 Stunden, ausgehend von 00, 06, 09, 12, 15, 18 und 21 UTC, 
neue +27-Stunden-Vorhersagen bereitgestellt. Der Modelllauf von 03 UTC liefert 
sogar einen Vorhersagezeitraum von +45 Stunden

@author: nixdorf
"""

#%% load libraries
import datetime
import numpy as np
from ftplib import FTP
import time
from io import BytesIO
import bz2
from rasterio.io import MemoryFile
# import required geometrytools from geotools
import geotools.geotools as gs
import re
import sys
import fiona
#%% First we write the ICON_IO downloader
def nwp_io(nwp_model='cosmo-d2',bufferhours=4,feature='tot_prec',shapefile='.\Examples\\Mueglitz_Basin.shp',number_frmt=np.float16):
    """
    icon_io checks the current time, downloads the datasets and converts to standard cells
    
    
    """
    print('Start retrieving NWP data for modell',nwp_model)
    initDf = True
    # Get the current time in order to choose appropriate forecast directory
    local_utctime=datetime.datetime.utcnow()
    #available runs 
    if nwp_model=='icon-eu':
        runtms=np.array((0,6,12,18))
    else:
        if nwp_model=='cosmo-d2':
           runtms=np.array((0,6,9,12,15,18,21))
        else:
           sys.exit('Unknown climate model have been chosen, only icon-eu and cosmo-d2 are integrated')
    #calculate nearest run with a buffer for simulation time of model
    delta_hours=list(runtms-(local_utctime.hour-bufferhours))
    smallest_hour=max(i for i in delta_hours if i < 0)
    runtm=runtms[delta_hours.index(smallest_hour)]
    
    # Connect to ftp
    
    server='opendata.dwd.de'
    connected = False
    
    while not connected:
        try:
            ftp = FTP(server)
            ftp.login()
            connected = True
    
        except:
            time.sleep(5)
            print('Reconnect to Server')
            pass
    if len(str(runtm))<2:
        ftp.cwd('weather/nwp/' + nwp_model + '/grib/'+'0'+str(runtm)+'/'+feature+'/')
    else:
        ftp.cwd('weather/nwp/' + nwp_model + '/grib/'+str(runtm)+'/'+feature+'/')   
    files_raw = ftp.nlst()
    #get the grid data only
    files=[i for i in files_raw if 'regular' in i]
    
    
    
    for file in files:
        print('Retrieving {}...'.format(file))
        retrieved=False
        archive = BytesIO()
        # try to retrieve file
        while not retrieved:
            try:
                ftp.retrbinary("RETR " + file, archive.write)
                retrieved=True
            except:
                print('reconnect to ftp')
                ftp = FTP(server)
                ftp.login()
                if len(str(runtm))<2:
                    ftp.cwd('weather/nwp/' + nwp_model + '/grib/'+'0'+str(runtm)+'/tot_prec/')
                else:
                    ftp.cwd('weather/nwp/' + nwp_model + '/grib/'+str(runtm)+'/tot_prec/')  
                
        archive.seek(0)
        grib_data = bz2.decompress(archive.read())
        #Get raster Data
        with MemoryFile(grib_data) as memfile:
            dataset=memfile.open()
            data_array = dataset.read()[0].astype(number_frmt)
            data_proj=dataset.crs
            if initDf:
                NaN_Value = dataset.nodata
                afn_transform = dataset.transform
                dataset_transform = (afn_transform[2],
                                  afn_transform[0],
                                  0,
                                  afn_transform[5],
                                  0,
                                  afn_transform[4])
                # do the complicated buffer clipping
                # if a shapefile exist
                if shapefile is not None:
                    data_clip, data_clip_transform, cols, rows = gs.buffered_raster_clipping(
                        data_array,
                        shape_inpt=shapefile,
                        raster_transfrm=
                        dataset_transform,
                        raster_proj=data_proj)
                else:
                    data_clip = data_array
                    data_clip_transform = dataset_transform
                    rows = [data_array.shape[0], 0]
                    cols = [0, data_array.shape[1]]
        
                #generate the footprint cells
                datacells = gs.create_footprint_cells(
                    transform=data_clip_transform,
                    data_size=data_clip.shape,
                    proj_crs=data_proj)
                #initialize the merged dataset
                data_stacked = data_clip
                #get the correct date attached to the simulation
                # get the numbers from filename
                filenm_nmbrs=re.findall(r'\d+', file)
                # if cosmo-2, we remove first number
                if nwp_model=='cosmo-d2':
                    del(filenm_nmbrs[0])
                sim_starttime=datetime.datetime.strptime(filenm_nmbrs[0], '%Y%m%d%H')
                sim_dates=[sim_starttime+ datetime.timedelta(hours=int(filenm_nmbrs[1]))]
                initDf = False
            # if we initialised already, computation is easy
            else:
                rado_clip_data = data_array[rows[
                    1]:rows[0], cols[0]:cols[1]]
                try:
                    data_stacked = np.dstack(
                        (data_stacked,
                         rado_clip_data))
                except Exception as e:
                    print(e)
                    sys.exit(
                        'Memory Error :-(, buy more RAM'
                    )
               # get the numbers from filename
                filenm_nmbrs=re.findall(r'\d+', file)
                # if cosmo-2, we remove first number
                if nwp_model=='cosmo-d2':
                    del(filenm_nmbrs[0])
                sim_dates.append(sim_starttime + datetime.timedelta(hours=int(filenm_nmbrs[1])))
        print('Processing {}...finished'.format(file))
        sim_dates = sorted(sim_dates)
    try:
        ftp.quit()
    except Exception as e:
        print(e)
    # repair the radocell crs
    datacells.crs = fiona.crs.to_string(datacells.crs)
    print('IO Operation of',nwp_model,' conducted successfully')
    return data_stacked, sim_dates, datacells

def nwpgeo(nwp_model='cosmo-d2',shape_inpt='.\Examples\Mueglitz_Basin.shp',feature='tot_prec',
              shape_integration=True,
              outpt_proj='epsg:25833',
              Output=True,
              outpt_nm='nwp_icon',number_frmt=np.float16):
    """
    Couples the four main function for the entire workflow for icon data processing
    """

    # first we start with icon_io to retrieve the data
    data_stacked, sim_dates, datacells = nwp_io(nwp_model=nwp_model,bufferhours=4,
                                                 feature=feature,
                                                 shapefile=shape_inpt,
                                                 number_frmt=number_frmt)
    #map data on vectorgrid
    cell_precip, precip_col_nms = gs.map_arraystack_on_cellgrd(
        data_stacked, sim_dates, datacells, numerator=1,
        Output=True, outpt_proj=outpt_proj,number_frmt=number_frmt)
    # delete stacked numpy array
    del data_stacked
    # Integrate values to boundary shape, if needed
    if shape_integration:
        precipcell_clip, gdf_shape = gs.map_cellgrd_on_polyg(
            cell_precip, shape_inpt=shape_inpt, outpt_proj=outpt_proj)
        
        gs.compute_polyg_values(precipcell_clip, gdf_shape, 
            header=nwp_model+'_rainfall_mm_h',
            datacol_type='AllDigits',
            Output=Output,
            outpt_proj=outpt_proj,
            outpt_nm=outpt_nm)
    else:
        print('No integration of data on geometry requested')
    return None

#%% Next part is the COSMO Dataset
    

nwpgeo()