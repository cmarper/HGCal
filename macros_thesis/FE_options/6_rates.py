# C. Martin Perez cmartinp@cern.ch, Nov. 2019

###########

import os
import pandas as pd
import numpy as np
import root_pandas
import pickle
import xgboost as xgb
from matplotlib import pyplot as plt
import matplotlib
import sys

###########

#algo = 'A'
algo = 'B'

PUBDTWP = '99'
#PUBDTWP = '95'
#PUBDTWP = '90'

###########

if algo == 'A':
  print 'Algorithm A'

elif algo == 'B':
  print 'Algorithm B'

else:
  print 'Error: No algorithm chosen!'

###########

if PUBDTWP == '99':
  text = 'WP99'
  print 'PU BDT WP99'

elif PUBDTWP == '95':
  text = 'WP95'
  print 'PU BDT WP95'

elif PUBDTWP == '90':
  text = 'WP90'
  print 'PU BDT WP90'

else:
  print 'Error: No PU BDT WP chosen!'

###########

fe_names = {}

fe_names[0] = 'Threshold'
fe_names[1] = 'STC'
fe_names[2] = 'BestChoice'
fe_names[3] = 'BestChoiceCoarse'
fe_names[4] = 'MixedBC+STC'

###########

file_in_nu = {}

if algo == 'A':

  dir_in = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoA/DMid/'

  file_in_nu[0] = dir_in+'ntuple_NuGun_THRESHOLD_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  file_in_nu[1] = dir_in+'ntuple_NuGun_SUPERTRIGGERCELL_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  file_in_nu[2] = dir_in+'ntuple_NuGun_BESTCHOICE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  file_in_nu[3] = dir_in+'ntuple_NuGun_BESTCHOICECOARSE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  file_in_nu[4] = dir_in+'ntuple_NuGun_MIXEDBCSTC_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'

if algo == 'B':

  dir_in = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoB/DMid/'

  file_in_nu[0] = dir_in+'ntuple_NuGun_THRESHOLD_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  file_in_nu[1] = dir_in+'ntuple_NuGun_SUPERTRIGGERCELL_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  file_in_nu[2] = dir_in+'ntuple_NuGun_BESTCHOICE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  file_in_nu[3] = dir_in+'ntuple_NuGun_BESTCHOICECOARSE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  file_in_nu[4] = dir_in+'ntuple_NuGun_MIXEDBCSTC_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'

###########

file_in_mapping = {}

if algo == 'A':

  dir_in = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoA/mapping/'

  file_in_mapping[0] = dir_in+'mapping_95eff_THRESHOLD_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.pkl'
  file_in_mapping[1] = dir_in+'mapping_95eff_SUPERTRIGGERCELL_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.pkl'
  file_in_mapping[2] = dir_in+'mapping_95eff_BESTCHOICE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.pkl'
  file_in_mapping[3] = dir_in+'mapping_95eff_BESTCHOICECOARSE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.pkl'
  file_in_mapping[4] = dir_in+'mapping_95eff_MIXEDBCSTC_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.pkl'

if algo == 'B':

  dir_in = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoB/mapping/'

  file_in_mapping[0] = dir_in+'mapping_95eff_THRESHOLD_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.pkl'
  file_in_mapping[1] = dir_in+'mapping_95eff_SUPERTRIGGERCELL_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.pkl'
  file_in_mapping[2] = dir_in+'mapping_95eff_BESTCHOICE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.pkl'
  file_in_mapping[3] = dir_in+'mapping_95eff_BESTCHOICECOARSE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.pkl'
  file_in_mapping[4] = dir_in+'mapping_95eff_MIXEDBCSTC_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.pkl'

pt_95s = {}

for name in file_in_mapping:

  with open(file_in_mapping[name], 'rb') as f:
    pt_95s[name] = pickle.load(f)

###########

df_nu = {}

for name in file_in_nu:

  store_nu = pd.HDFStore(file_in_nu[name], mode='r')
  df_nu[name] = store_nu['df_nu_PU200']
  store_nu.close()

###########

events_total = {}
events_pubdtwp = {}

for name in file_in_nu: 

  events_total[name] = np.unique(df_nu[name].reset_index()['event']).shape[0]

  if PUBDTWP == '99':
    sel = df_nu[name]['cl3d_pubdt_passWP99'] == True

  elif PUBDTWP == '95':
    sel = df_nu[name]['cl3d_pubdt_passWP95'] == True

  elif PUBDTWP == '90':
    sel = df_nu[name]['cl3d_pubdt_passWP90'] == True

  df_nu[name] = df_nu[name][sel]

  df_nu[name].set_index('event', inplace=True)

  group = df_nu[name].groupby('event')
  df_nu[name]['leading_cl3d_pt_c3'] = group['cl3d_pt_c3'].max()

  df_nu[name]['cl3d_isleading'] = df_nu[name]['leading_cl3d_pt_c3'] == df_nu[name]['cl3d_pt_c3']

  sel = df_nu[name]['cl3d_isleading'] == True
  df_nu[name] = df_nu[name][sel]

  events_pubdtwp[name] = np.unique(df_nu[name].reset_index()['event']).shape[0]

  print fe_names[name], float(events_pubdtwp[name])/float(events_total[name])


rates = {}
rates_rescale = {}

#scale=2760*11246./1000  # bunches * frequency , /1000 to pass to kHz
scale = 1

from decimal import *
getcontext().prec = 30

for name in file_in_nu:

  tmp = df_nu[name][df_nu[name].cl3d_pt_c3 > 10]

  tmp['cl3d_pt_95'] = tmp.cl3d_pt_c3.apply(lambda x : np.interp(x, pt_95s[name].threshold, pt_95s[name].pt95))

  rate = np.arange( float(tmp.shape[0]*scale)/events_total[name], 0., (-1.*scale)/events_total[name])
  rate = rate[rate>0]

  rates[name] = ( np.sort(tmp.cl3d_pt_c3), rate )
  #print rates[name][0],rates[name][1]
  rates_rescale[name] = ( np.sort(tmp.cl3d_pt_95), rate )
  #print rates_rescale[name][0],rates_rescale[name][1]

#####

colors = {}
colors[0] = 'blue'
colors[1] = 'red'
colors[2] = 'green'
colors[3] = 'orange'
colors[4] = 'fuchsia'

legends = {}
legends[0] = 'Threshold'
legends[1] = 'STC'
legends[2] = 'BC'
legends[3] = 'BC Coarse'
legends[4] = 'Mixed BC+STC'

#####

if algo == 'A':
  plotdir = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoA/plots/rates/'

elif algo == 'B':
  plotdir = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoB/plots/rates/'

#####

import geeksw.plotting.cmsplot as plt
from geeksw.plotting.root_colors import *

plt.matplotlib.font_manager._rebuild()

import matplotlib.lines as mlines

'''
matplotlib.rcParams.update({'font.size': 22})

interp_points = np.arange(10, 80, 0.5)
rate_ref = rates[0]
rate_ref_interp = np.interp(interp_points, rate_ref[0], rate_ref[1]/scale)
fig, axs = plt.subplots(2, 1, figsize=(15,20))

for name in file_in_nu:

  rate = rates[name]
  #plt.plot(rate[0],rate[1]*scale,label=legends[name], linewidth=2, color=colors[name])
  #plt.plot(rate[0],rate[1],label=legends[name], linewidth=2, color=colors[name])
  axs[0].plot(rate[0], rate[1]/scale, label=legends[name], linewidth=2, color=colors[name])
  rate_interp = np.interp(interp_points, rate[0], rate[1]/scale)
  ratio =  rate_interp / rate_ref_interp
  axs[1].plot(interp_points, ratio, label=legends[name], linewidth=2, color=colors[name])

axs[0].legend(loc = 'upper right', fontsize=22)
axs[0].set_yscale("log")
axs[0].set_xlim(10, 80)
axs[0].set_ylim(6e-4, 1)
axs[0].grid()
#axs[0].set_xlabel('L1 threshold [GeV]')
axs[0].set_ylabel('Rate [a.u.]')
axs[1].set_ylim(0.6,2.0)
axs[1].set_xlim(10, 80)
axs[1].set_xlabel('L1 threshold [GeV]')
#axs[1].set_ylabel('Ratio to STC')
axs[1].set_ylabel('Ratio to Threshold')
axs[1].grid()
plt.savefig(plotdir+'rate_L1_algo'+algo+'_WP'+PUBDTWP+'_v2.png')
plt.savefig(plotdir+'rate_L1_algo'+algo+'_WP'+PUBDTWP+'_v2.pdf')
'''

#plt.legend(loc = 'upper right', fontsize=16)
#plt.yscale('log')
#plt.xlabel('L1 threshold [GeV]')
#plt.ylabel('Rate [a.u.]')
#plt.xlim(20, 80)
#plt.ylim(0.0007, 0.3)
#plt.grid()
#plt.savefig(plotdir+'rate_L1_algo'+algo+'_WP'+PUBDTWP+'.png')
#plt.savefig(plotdir+'rate_L1_algo'+algo+'_WP'+PUBDTWP+'.pdf')
#plt.show()

'''
plt.figure(figsize=(10,10))

for name in file_in_nu:

  rate = rates_rescale[name]
  #plt.plot(rate[0],rate[1]*scale,label=legends[name], linewidth=2, color=colors[name])
  plt.plot(rate[0],rate[1],label=legends[name], linewidth=2, color=colors[name])

plt.legend(loc = 'upper right', fontsize=16)
plt.yscale('log')
plt.xlabel('Offline threshold [GeV]')
plt.ylabel('Rate [a.u.]')
plt.xlim(30, 90)
plt.ylim(0.1, 0.6)
plt.grid()
plt.savefig(plotdir+'rate_offline_algo'+algo+'_WP'+PUBDTWP+'.png')
plt.savefig(plotdir+'rate_offline_algo'+algo+'_WP'+PUBDTWP+'.pdf')
plt.show()
'''

interp_points = np.arange(10, 90, 0.5)
rate_ref = rates_rescale[0]
rate_ref_interp = np.interp(interp_points, rate_ref[0], rate_ref[1]/scale)
fig, axs = plt.subplots(2, 1, figsize=(8,8))

for name in file_in_nu:

  rate = rates_rescale[name]
  #plt.plot(rate[0],rate[1]*scale,label=legends[name], linewidth=2, color=colors[name])
  #plt.plot(rate[0],rate[1],label=legends[name], linewidth=2, color=colors[name])
  axs[0].plot(rate[0], rate[1]/scale, label=legends[name], linewidth=2, color=colors[name])
  rate_interp = np.interp(interp_points, rate[0], rate[1]/scale)
  ratio =  rate_interp / rate_ref_interp
  axs[1].plot(interp_points, ratio, label=legends[name], linewidth=2, color=colors[name])

axs[0].legend(loc = 'upper right', fontsize=14)
axs[0].set_yscale("log")
axs[0].set_xlim(20, 90)
axs[0].set_ylim(6e-4, 1)
axs[0].grid()
#axs[0].set_xlabel('Offline threshold [GeV]')
axs[0].set_ylabel('Rate [a.u.]')
axs[1].set_ylim(0.6,2.0)
axs[1].set_xlim(20, 90)
axs[1].set_xlabel('Offline threshold [GeV]')
#axs[1].set_ylabel('Ratio to STC')
axs[1].set_ylabel('Ratio to Threshold')
axs[1].grid()
plt.savefig(plotdir+'rate_offline_algo'+algo+'_WP'+PUBDTWP+'_v2.png')
plt.savefig(plotdir+'rate_offline_algo'+algo+'_WP'+PUBDTWP+'_v2.pdf')
