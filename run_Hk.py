'''
Script to calculate receiver functions and run Hk stacking
'''

import numpy as np
import matplotlib.pyplot as plt
from obspy import taup
import pandas as pd
import rf
import geopy
from geopy.distance import VincentyDistance, distance
from shapely.geometry import LineString, Point
from rf_tools import *
from latexify import latexify

# IMPORT PASSIVE SEISMIC DATA
#eqdir='/Users/brookkeats/Dropbox/Documents/Oxford_DPhil/Misc_data_files/RF_data/'
eqdir='/Users/brookkeats/Documents/DPhil/Data/RF_data/'

#eqfile='PI_netw_M5.5+v2_raw_data.h5'
#eqfile='PI_netw_M5.5+v2_instr_corr_data.h5'
eqfile='NCMS_ASU_M5.5+instr_corr_data.h5'

eq_st = rf.read_rf(eqdir+eqfile)

# PRE PROCESS DATA
eq_st.filter('bandpass', freqmin=1./20, freqmax=4., corners=2, zerophase=True)
#eq_st.filter('lowpass', freq=4., corners=2, zerophase=True)

# FILTER EVENTS BY SNR
SNR_llim = 2
snr_eq_st = SNR_st_filter(eq_st, SNR_llim=SNR_llim)

# CALCULATE RECEIVER FUNCTIONS
fmax=1
bp_filt = dict(type='bandpass', freqmin=1/20., freqmax=fmax)
freq_kws = dict(waterlevel=0.05, gauss=2., normalize=0)
mtc_kws = dict(module='pmtm', fc=fmax)

#rf_st = calculate_rfs(snr_eq_st, filt_kw = bp_filt, deconvolve="time")
#rf_st = calculate_rfs(snr_eq_st, deconvolve="freq", **freq_kws)
rf_st = calculate_rfs(snr_eq_st, deconvolve='MTC', **mtc_kws)
#rf_st = calculate_rfs(snr_eq_st, deconvolve='ET-MTC', **mtc_kws)

# AUTO QC RESULTS
qc_rf_st = qc_rfs(rf_st, qt_ll=2, pp_ll=2, return_st=True)

qc_rf_qst = qc_rf_st.select(component="Q")

# PLOT RF STACK
kw = {'trim': (-5,25), 'fillcolors': ('black', 'gray'), 'trace_height': 0.1}

qc_rf_qst.sort(['back_azimuth']).plot_rf(**kw)
#rffig = "{0}_qc_rf_results.pdf".format(stat)
#rf_qst.sort(['back_azimuth']).plot_rf(fname=savedir+rffig, **kw)
plt.show()

# FILTER RESULTS BY BACK AZIMUTH/RAY PARAMETER
try:
    baz_min, baz_max = input("Select back azimuth range (e.g [0, 360]) : ")
except SyntaxError:
    baz_min, baz_max = [0,360]
filt_back_az(qc_rf_st, baz_min=baz_min, baz_max=baz_max)

print("\nDisplaying selected events : ")
#stat_qst = stat_st.select(component="Q")
filt_qst = qc_rf_st.select(component="Q")
filt_qst.sort(['distance']).plot_rf(**kw)
plt.show()

try:
    deg_min, deg_max = input("Select epicentre distance range (e.g [30, 98]) : ") or [30,90]
except SyntaxError:
    deg_min, deg_max = [30,98]
filt_epi_dist(qc_rf_st, deg_min=deg_min, deg_max=deg_max)

# CHECK STREAM CONTAINS ONLY A SINGLE STATION
stat = qc_rf_st[0].stats.station
for tr in eq_st:
    assert tr.stats.station == stat

# SAVE EQ CATALOGUES TO FILE
cat_dir = '/Users/brookkeats/Documents/DPhil/Data/RF_data/EQ_catalogues/'
full_cat = '{}_full_rf_eq_cat.txt'.format(stat)
used_cat = '{}_used_rf_eq_cat.txt'.format(stat)
get_st_cat(rf_st, savefile=cat_dir+full_cat)
get_st_cat(qc_rf_st, savefile=cat_dir+used_cat)


# PLOT HK STACK
H_min, H_max = [5,50]
k_min, k_max = [1.6,1.9]
vp_av = 6.5

w1, w2, w3 = 0.5, 0.3, 0.2
savedir = "/Users/brookkeats/Dropbox/Documents/Oxford_DPhil/Figures/Thesis/Receiver_functions/"

hkfig = "{0}_{1}-{2}baz_{3}-{4}deg_{5}-{6}km_qc_Hk_stack.pdf".format(stat,baz_min, baz_max, deg_min, deg_max, H_min, H_max)
H, k = Hk_stack(qc_rf_st, Vp_av=vp_av, w1=w1, w2=w2, w3=w3, H_min=H_min, H_max=H_max, k_min=k_min, k_max=k_max, plot=True)
#H, k = Hk_stack(rf_st, Vp_av=vp_av, w1=0.6, w2=0.3, w3=0.1, H_min=H_min, H_max=H_max, plot=True, savefigure=savedir+hkfig)

# BOOTSTRAP HK
B=50
hkbsfig = "{}_{}_{}-{}baz_{}-{}deg_{}-{}km_qc_Hk_bs_stack.pdf".format(stat, fmax, baz_min, baz_max, deg_min, deg_max, H_min, H_max)
#H, H_std, k, k_std, angle = Hk_stack(qc_rf_st, Vp_av=vp_av, w1=w1, w2=w2, w3=w3, bootstrap=B, H_min=H_min, H_max=H_max, k_min=k_min, k_max=k_max, plot=True, savefigure=savedir+hkbsfig)

print H, k

multfig = "{0}_{1}-{2}baz_{3}-{4}deg_{5}-{6}km_qc_rf_multiples.pdf".format(stat,baz_min, baz_max, deg_min, deg_max, H_min, H_max)
plot_rf_multiples(qc_rf_st, H=H, k=k, vp_av=vp_av)
#plot_rf_multiples(qc_rf_st, H=H, k=k, vp_av=vp_av, sort='back_azimuth', savefigure=savedir+multfig)
#H, k = Hk_stack(rf_st, Vp_av=vp_av, w1=0.5, w2=0.3, w3=0.2, H_min=H_min, H_max=H_max, plot=True, bootstrap=B)

