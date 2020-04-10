#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 17:16:03 2018

A collection of Python functions to aid with 
receiver function analysis in rf

@author: Brook Keats
"""

from __future__ import print_function
from obspy import read
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
from obspy.taup import TauPyModel
import rf
import glob
import os
from tqdm import tqdm
from geopy.distance import vincenty

# OBTAIN EARTHQUAKE CATALOG FOR RECEIVER FUNCTION ANALYSIS
def gen_eq_catalog(starttime, endtime, centre_lat, centre_lon, min_mag=5.5, min_r=30, max_r=90, plot=True, savefile=None):
    """
    Obtain an earthquake catalog for receiver function analysis
    
    :param starttime: start time for earthquake search (UTCDateTime object)
    :param endtime: end time for earthquake search (UTCDateTime object)
    :param centre_lat: latitude of central point for earthquake search (e.g. station latitude)
    :param centre_lon: longitude of central point for earthquake search (e.g. station longitude)
    :param min_mag: minimum earthquake magnitude for earthquake search
    :param min_r: minimum earthquake search radius in degrees
    :param max_r: maximum earthquake search radius in degrees
    :param plot: generate a plot of the resulting earthquake catalog (True or False)
    :param savefile: save earthquake catalog to XML file (savepath/filename)
    
    """
    # SELECT DATA PROVIDER
    client = Client("USGS")
    
    # CREATE CATALOG
    cat = client.get_events(starttime=starttime, endtime=endtime, latitude=centre_lat, longitude=centre_lon, minradius=min_r, maxradius=max_r, minmagnitude=min_mag)

    # PLOT CATALOG
    if plot == True:
        cat.plot()
    
    # WRITE CATALOG TO FILE
    if savefile != None:
        cat.write(savefile, format="QUAKEML")
    else:
        return cat


# GENERATE STATION DICT TO CALCULATE RFSTATS
def get_station_info(stat, inv, xml=False):
    """
    
    :param stat: station code (string)
    :param inv: obspy network inventory object containing stations
    :param xml: return station as obspy inventory object instead of dict
    """
    found=0
    for netw in inv:
        #print(netw)
        for inv_sta in netw:
            #print(inv_sta)
            if inv_sta.code == stat:
                op_stat = inv_sta
                stat_dict = dict(latitude=inv_sta.latitude, longitude=inv_sta.longitude, elevation=inv_sta.elevation)
                found=1
                if xml:
                    return op_stat
                else:
                    return stat_dict
                break
    if not found:
        print("Station {0} not found in NCMS network".format(stat))


# CALCULATE EXPECTED ARRIVAL TIME AT STATION FROM GLOBAL VELOCITY MODELS FOR EARTHQUAKES IN CATALOG
def get_tt_delay(sta_lat, sta_lon, ev_lat, ev_lon, ev_depth, phase_list=['P']):
    """
    :param sta_lat: latitude of seismic station
    :param sta_lon: longitude of seismic station
    :param ev_lat: latitude of earthquake event
    :param ev_lon: longitude of earthquake event
    :param ev_depth: depth of earthquake event in kilometers
    :param phase_list: array containing list of phases to calculate arrivals for    
    """
    km_per_deg = 111.195
    dist_deg=(vincenty((sta_lat, sta_lon), (ev_lat, ev_lon)).km)/km_per_deg
    #print("Event is {0} degrees from station at a depth of {1} km".format(round(dist_deg,4), ev_depth))
    model = TauPyModel(model='iasp91')
    arrivals = model.get_travel_times(source_depth_in_km=ev_depth, distance_in_degree=dist_deg, phase_list=phase_list)
    #arr = arrivals[0]
    #delay = arr.time
    #return delay
    return arrivals
    #print(arrivals)
    
    
# READ WAVEFORM FILES FROM MARAID2 DATABASE
def read_passive(op_stat, event, st=False, filt_kws=None, **kwargs):

    # GET EVENT INFO
    ev_datetime = event.origins[0].time
    ev_lat = event.origins[0].latitude
    ev_lon = event.origins[0].longitude
    ev_dep = event.origins[0].depth/1000.
    print("Event origin time", ev_datetime)

    # GET STATION INFO
    sta_lat = op_stat.latitude
    sta_lon = op_stat.longitude
    sta_code = op_stat.code

    # CALCULATE EXPECTED ARRIVAL TIME
    arrivals = get_tt_delay(sta_lat, sta_lon, ev_lat, ev_lon, ev_dep, **kwargs)
    delay = arrivals[0].time    
    exp_arr_time = ev_datetime + delay
    print("Expected arrival time", exp_arr_time)

    # FIND FILE (def_path is a function to find the correct miniseed file for this event)
    seis_file = def_path(sta_code, exp_arr_time)
    # IF EQ FALLS WITHIN 1 MIN OF THE HOUR FETCH PREV/NEXT HOURS WAVEFORM    
    if seis_file == None:
        return
    if exp_arr_time.minute == 0:
        print("Fetching previous hours waveform as well")
        prev_file = def_path(sta_code, exp_arr_time-60)
        if prev_file == None:
            return
        else:
            seis_file += prev_file
    elif exp_arr_time.minute == 59 or exp_arr_time.minute == 58:
        print("Fetching next hours waveform as well")
        next_file = def_path(sta_code, exp_arr_time+60)
        if next_file == None:
            return
        else:
            seis_file += next_file

    # READ DATA WITH OBSPY
    print("Using file(s) {0}".format(seis_file))    
    for i, seis in enumerate(seis_file):
        if i == 0:
            st = read(seis)
        else:
            st += read(seis)
            
    if len(st) > 0:
        st = st.select(channel="HH?")
        if len(st) == 0:
            print("No HH? traces in stream...")
            return
    else:
        print("No traces in stream...")
        return
    
    filt_st = st.copy()
    filt_st.merge(method=1, fill_value='interpolate', interpolation_samples=-1)
    
    # DOWNSAMPLE TO 25 HZ
    if filt_kws != None:
        filt_st.filter(**filt_kws)
    samp_rate = st[0].stats.sampling_rate
    # CHECK SAMPLE RATE IS THE SAME ON ALL TRACES
    for tr in filt_st:
        tr_sr = tr.stats.sampling_rate        
        assert samp_rate == tr_sr, "Sampling rate differs across components"
    
    dec_factor = int(samp_rate/25.)
    filt_st = filt_st.detrend(type='simple')
    filt_st = filt_st.detrend(type='demean')    
    filt_st.decimate(dec_factor)  
    
    # CUT AROUND EXPECTED ARRIVAL TIME
    filt_st = filt_st.slice(starttime=exp_arr_time-30, endtime=exp_arr_time+60)
    
    #filt_st = filt_st.filter('bandpass', freqmin=1/30., freqmax=1/2.5)
    
    # return stream object if requested
    if st:
        return filt_st


# DEFINE PATHS TO DATA FILES (date in DD/MM/YYYY and time in HH:MM:SS.SS)
def def_path(stat, datetime):
    """
    
    :param stat: station code (string)
    :param datetime: UTCDateTime object with the projected 
                     arrival time of the primary P phase
    
    """
    
    # SET UP BASE DIRECTORIES
    seis_base = "/home/maraid2/brookk/"
    ncms_data = seis_base + "NCMS_data/Raw_data/"
    ncms_dataless = ncms_data + "meta data/Dataless/"

    # SET UP LISTS WITH STATION CODES OF THE NCMS NETWORKS
    ncms_ae = ["AJN", "ALN", "JRN", "MSF", "SHM", "SRB", "SAK", "UMQ", "MZR"]
    ncms_ue = ["UMZA", "MZWR", "GHWR", "SLWR"]
    ncms_dn = ["HAT", "FAQ", "ASU", "NAZ"]
    ncms_om = ["BAN", "HOQ", "MDH", "ASH", "SOH", "BID", "SHA"]
    
    year = str(datetime.year)
    month = str(datetime.month).zfill(2)
    day = str(datetime.day).zfill(2)
    
    hour = str(datetime.hour).zfill(2)
    minute = str(datetime.minute).zfill(2)
    second = str(datetime.second).zfill(2)

    if stat in ncms_ae:
        data_path = ncms_data + "AE-Network/" + year + "/" + month + "/" + day + "/"
        data_file = glob.glob(data_path + '*' + stat + '*_{0}0000.*seed'.format(hour))
        
    elif stat in ncms_ue:
        data_path = ncms_data + "UE-Network/" + year + "/" + month + "/" + day + "/"
        data_file = glob.glob(data_path + '*' + stat + '*_{0}0000.*seed'.format(hour))
        
    elif stat in ncms_dn:
        data_path = ncms_data + "OM&DN-Network/" + year + "/" + month + "/" + day + "/"
        data_file = glob.glob(data_path + '*' + stat + '*_{0}0000.*seed'.format(hour))
        
    elif stat in ncms_om:
        data_path = ncms_data + "OM&DN-Network/" + year + "/" + month + "/" + day + "/"
        data_file = glob.glob(data_path + '*' + stat + '*_{0}0000.*seed'.format(hour))
        
    
    else:
        data_file = []
        print("Station not found...")
        return
        #sys.exit("Station not found... aborting script.")
    
    if not os.path.isdir(data_path):
        print("Directory {0} does not exist...".format(data_path))
        return        
        
    #print("Searching directory", data_path, "for file", data_file)
    if len(data_file) < 1:
        print("Data file does not exist...")
        return
        #sys.exit("No data file found... aborting script.")
    else:
        datafile = data_file
        #print("Using data file(s)", datafile)
        return datafile



# GET WAVEFORMS FOR RF ANALYSIS FROM SEISMIC DATABASE FOR STATIONS IN STATION_LIST  
def get_eq_waveforms(station_list, sta_inv, eq_cat, filt_kws=None, **kwargs):
    """
    
    :param station_list: list of stations to retrieve data from
    :param sta_inv: obspy network inventory object containing station information
    :param eq_cat: obspy inventory object containing earthquake catalog
    
    """
    ev_cnt = 0
    for stat in station_list:
        op_stat = get_station_info(stat, sta_inv, xml=True)
        stat_dict = get_station_info(stat, sta_inv)    
        for i, event in enumerate(eq_cat):        
            cat_id = event.resource_id
            test_id = str(cat_id).split("=")[1].split("&")[0]
            print("\nSearching for event {0} at station {1}".format(test_id, stat))
            waveforms = read_passive(op_stat, event, filt_kws=filt_kws, **kwargs)
            if waveforms != None:
                if len(waveforms) == 3:
                    #print(len(waveforms))
                    if ev_cnt == 0:
                        #op_stream = waveforms
                        eq_stream = rf.RFStream(waveforms)
                        #print kwargs
                        if 'phase_list' in kwargs:
                            phase = kwargs['phase_list'][0]
                        else:
                            phase = 'P'
                            
                        stats = rf.rfstats(station=stat_dict, event=event, phase=phase, dist_range=(30,90))
                        for tr in eq_stream:
                            tr.stats.update(stats)
                    else:
                        #op_stream += waveforms
                        temp_stream = rf.RFStream(waveforms)
                        stats = rf.rfstats(station=stat_dict, event=event, phase=phase, dist_range=(30,90))                        
                        #print(stats)
                        #print(stat_dict, event)
                        if stats == None:
                            print("No rfstats calculated... skipping event")
                            continue
                        else:
                            for tr in temp_stream:
                                tr.stats.update(stats)
                            eq_stream.extend(temp_stream)
                    ev_cnt += 1
                else:
                    print("Imported stream does not have 3 traces ({0})... skipping event".format(len(waveforms)))
    
    if ev_cnt != 0:
        return eq_stream
    else:
        print("No earthquake waveforms found...")
        return None

# CALCULATE RECEIVER FUNCTIONS
def calculate_rfs(eq_stream, filt_kw=None, deconvolve='time', moveout=True, savefile=None, **kwargs):
    """
    
    :param eq_stream: obspy stream object containing earthquake waveforms
    :param filt_kw: dict object containing parameters for an obspy filter
    :param deconvolve: string to select deconvolution method ('time' or 'freq')
    :param moveout: correct time delays of receiver function results for moveout (True or False)
    :param savefile: save result to hdf5 file (string with the full file path to file)
    """
    rf_stream = rf.RFStream()
    working_stream = eq_stream.copy()
    for stream3c in tqdm(rf.IterMultipleComponents(working_stream, 'onset', number_components=3)):
        bad_npts = 0
        bad_start_t = 0
        for i, tr in enumerate(stream3c):
            samp_rate = tr.stats.sampling_rate
            #if samp_rate != 25:
            #    warnings.warn("Sampling rate for tr is {0}".format(samp_rate))
            npts = tr.stats.npts
            thr_npts = samp_rate*90    # stream should have sr*streamlength datapoints
            if (npts - thr_npts) == 1:
                thr_npts = thr_npts + 1
            if thr_npts != npts:
                print("Expected npts {0}, actual number {1}... skipping event".format(thr_npts, npts))
                bad_npts = 1
            if i == 0:
                ref_t = tr.stats.starttime
            else:
                start_t = tr.stats.starttime
                if start_t != ref_t:
                    print("Inconsistent start times ({0} and {1}) in traces... skipping event".format(ref_t, start_t))
                    bad_start_t = 1
        
        if not bad_npts and not bad_start_t:
            if len(stream3c) != 3:
                continue

            stream3c.trim2(-25, 75, 'onset')
            stream3c.rotate('ZNE->LQT')
            #stream3c.deconvolve(method=deconvolve)
            stream3c.rf(deconvolve=deconvolve, filter=filt_kw, **kwargs)
            if moveout:                
                stream3c.moveout()
            #print(stream3c)
            rf_stream.extend(stream3c)

    if savefile != None:
        rf_stream.write(savefile, 'H5')
    
    return rf_stream
