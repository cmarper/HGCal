# C. Martin Perez cmartinp@cern.ch, Nov. 2019

###########

import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import pickle

###########

#algo = 'A'
algo = 'B'

#PUBDTWP = '99'
PUBDTWP = '95'
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
  bdtcut = 'cl3d_pubdt_passWP99'
  text = 'WP99'
  print 'PU BDT WP99'

elif PUBDTWP == '95':
  bdtcut = 'cl3d_pubdt_passWP95'
  text = 'WP95'
  print 'PU BDT WP95'

elif PUBDTWP == '90':
  bdtcut = 'cl3d_pubdt_passWP90'
  text = 'WP90'
  print 'PU BDT WP90'

else:
  print 'Error: No PU BDT WP chosen!'

###########

fe_names = {}

fe_names[0] = 'Thr'
fe_names[1] = 'STC'
fe_names[2] = 'BC'
fe_names[3] = 'BCCoarse'
fe_names[4] = 'BC+STC'

###########

file_in_tau = {}

if algo == 'A':

  dir_in = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoA/DMid/'

  file_in_tau[0] = dir_in+'ntuple_TauGun_THRESHOLD_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  file_in_tau[1] = dir_in+'ntuple_TauGun_SUPERTRIGGERCELL_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  file_in_tau[2] = dir_in+'ntuple_TauGun_BESTCHOICE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  file_in_tau[3] = dir_in+'ntuple_TauGun_BESTCHOICECOARSE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  file_in_tau[4] = dir_in+'ntuple_TauGun_MIXEDBCSTC_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'

if algo == 'B':

  dir_in = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoB/DMid/'

  file_in_tau[0] = dir_in+'ntuple_TauGun_THRESHOLD_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  file_in_tau[1] = dir_in+'ntuple_TauGun_SUPERTRIGGERCELL_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  file_in_tau[2] = dir_in+'ntuple_TauGun_BESTCHOICE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  file_in_tau[3] = dir_in+'ntuple_TauGun_BESTCHOICECOARSE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  file_in_tau[4] = dir_in+'ntuple_TauGun_MIXEDBCSTC_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'

###########

file_out_mapping = {}

if algo == 'A':

  dir_out = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoA/mapping/'

  file_out_mapping[0] = dir_out+'mapping_95eff_THRESHOLD_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.pkl'
  file_out_mapping[1] = dir_out+'mapping_95eff_SUPERTRIGGERCELL_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.pkl'
  file_out_mapping[2] = dir_out+'mapping_95eff_BESTCHOICE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.pkl'
  file_out_mapping[3] = dir_out+'mapping_95eff_BESTCHOICECOARSE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.pkl'
  file_out_mapping[4] = dir_out+'mapping_95eff_MIXEDBCSTC_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.pkl'

if algo == 'B':

  dir_out = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoB/mapping/'

  file_out_mapping[0] = dir_out+'mapping_95eff_THRESHOLD_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.pkl'
  file_out_mapping[1] = dir_out+'mapping_95eff_SUPERTRIGGERCELL_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.pkl'
  file_out_mapping[2] = dir_out+'mapping_95eff_BESTCHOICE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.pkl'
  file_out_mapping[3] = dir_out+'mapping_95eff_BESTCHOICECOARSE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.pkl'
  file_out_mapping[4] = dir_out+'mapping_95eff_MIXEDBCSTC_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.pkl'

###########

if PUBDTWP == '99':
  bdtcut = 'cl3d_pubdt_passWP99'

elif PUBDTWP == '95':
  bdtcut = 'cl3d_pubdt_passWP95'

elif PUBDTWP == '90':
  bdtcut = 'cl3d_pubdt_passWP90'

###########

df_tau = {}

for name in file_in_tau:

  store_tau = pd.HDFStore(file_in_tau[name], mode='r')
  df_tau[name]  = store_tau['df_tau_PU200']
  store_tau.close()

###########

for name in df_tau:

  sel = df_tau[name]['gentau_vis_pt'] > 1
  df_tau[name] = df_tau[name][sel]
  
  sel = np.abs(df_tau[name]['gentau_vis_eta']) > 1.6
  df_tau[name] = df_tau[name][sel]
  
  sel = np.abs(df_tau[name]['gentau_vis_eta']) < 2.9
  df_tau[name] = df_tau[name][sel]
  
  sel = df_tau[name]['cl3d_isbestmatch'] == True
  df_tau[name] = df_tau[name][sel]

  sel = df_tau[name]['cl3d_pt'] > 4
  df_tau[name] = df_tau[name][sel]
  
  sel = df_tau[name][bdtcut] == True
  df_tau[name] = df_tau[name][sel]
  
  sel = df_tau[name]['gentau_decayMode'] >= 0
  df_tau[name] = df_tau[name][sel]

###########

if algo == 'A':
  plotdir = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoA/plots/turnons/'

elif algo == 'B':
  plotdir = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoB/plots/turnons/'

ptcut = 1
etamin = 1.6

for name in df_tau:

  df_tau[name]['gentau_vis_abseta'] = np.abs(df_tau[name]['gentau_vis_eta'])
  df_tau[name]['gentau_bin_eta'] = ((df_tau[name]['gentau_vis_abseta'] - etamin)/0.1).astype('int32')
  df_tau[name]['gentau_bin_pt']  = ((df_tau[name]['gentau_vis_pt'] - ptcut)/1).astype('int32')

def efficiency(group, threshold):

  tot = group.shape[0]
  sel = group[(group.cl3d_pt_c3 > threshold)].shape[0]
  return float(sel)/float(tot)

efficiencies_vs_pt = {}

for name in df_tau:

  efficiencies_vs_pt[name] = {}
  efficiencies_vs_pt[name] = df_tau[name].groupby('gentau_bin_pt').mean()

turnon_thresholds = range(1, 100, 1)

for name in df_tau:
  for threshold in turnon_thresholds:
    eff = df_tau[name].groupby('gentau_bin_pt').apply(lambda x : efficiency(x, threshold))
    eff_smooth = eff.rolling(window=15, win_type='triang', center=True).mean()
    #eff_smooth.fillna(0., inplace=True)
    efficiencies_vs_pt[name]['efficiency_{}'.format(threshold)] = eff
    efficiencies_vs_pt[name]['efficiency_smooth_{}'.format(threshold)] = eff_smooth

pt_95s = {}

for name, df in efficiencies_vs_pt.items():
  mapping = {'threshold':[], 'pt95':[]}
  for threshold in turnon_thresholds:
    eff_smooth = df['efficiency_smooth_{}'.format(threshold)]
    pt_95 = np.interp(0.95, eff_smooth, df.gentau_vis_pt)
    mapping['threshold'].append(threshold)
    mapping['pt95'].append(pt_95)
    #print threshold, pt_95
  pt_95s[name] = pd.DataFrame(mapping)

for name in pt_95s:
  with open(file_out_mapping[name], 'wb') as f:
    pickle.dump(pt_95s[name], f)

##############

colors = {}
colors[0] = 'blue'
colors[1] = 'red'
colors[2] = 'olive'
colors[3] = 'orange'
colors[4] = 'fuchsia'

legends = {}
legends[0] = 'Threshold 1.35 mipT'
legends[1] = 'STC4+16'
legends[2] = 'BC Decentral'
legends[3] = 'BC Coarse 2x2 TC'
legends[4] = 'Mixed BC + STC'

matplotlib.rcParams.update({'font.size': 22})

plt.figure(figsize=(10,10))
for name in df_tau: 
  df = efficiencies_vs_pt[name]
  eff = df['efficiency_40']
  eff_smooth = df['efficiency_smooth_40']
  plt.plot(df.gentau_vis_pt, eff_smooth, label=legends[name], linewidth=2, color=colors[name])
plt.ylim(0., 1.01)
plt.xlim(28, 90)
plt.legend(loc = 'lower right', fontsize=16)
plt.xlabel(r'Gen. tau $p_{T}\,[GeV]$')
plt.ylabel('Efficiency')
plt.grid()
plt.savefig(plotdir+'eff_vs_pt_L1_40_WP_algo'+algo+'_'+PUBDTWP+'.png')
plt.savefig(plotdir+'eff_vs_pt_L1_40_WP_algo'+algo+'_'+PUBDTWP+'.pdf')

plt.figure(figsize=(10,10))
for name in df_tau: 
  df = efficiencies_vs_pt[name]
  eff = df['efficiency_60']
  eff_smooth = df['efficiency_smooth_60']
  plt.plot(df.gentau_vis_pt, eff_smooth, label=legends[name], linewidth=2, color=colors[name])
plt.ylim(0., 1.01)
plt.xlim(28, 90)
plt.legend(loc = 'lower right', fontsize=16)
plt.xlabel(r'Gen. tau $p_{T}\,[GeV]$')
plt.ylabel('Efficiency')
plt.grid()
plt.savefig(plotdir+'eff_vs_pt_L1_60_WP_algo'+algo+'_'+PUBDTWP+'.png')
plt.savefig(plotdir+'eff_vs_pt_L1_60_WP_algo'+algo+'_'+PUBDTWP+'.pdf')

plt.figure(figsize=(10,10))
for name in df_tau:
  df = pt_95s[name]
  plt.plot(df.threshold, df.pt95, label=legends[name], linewidth=2, color=colors[name])
plt.legend(loc = 'lower right', fontsize=16)
plt.xlabel('L1 Threshold [GeV]')
plt.ylabel('Offline threshold [GeV]')
plt.xlim(10, 80)
plt.ylim(10, 100)
plt.grid()
plt.savefig(plotdir+'L1_to_offline_WP_algo'+algo+'_'+PUBDTWP+'.png')
plt.savefig(plotdir+'L1_to_offline_WP_algo'+algo+'_'+PUBDTWP+'.pdf')
