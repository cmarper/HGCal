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
#fe_names[1] = 'STC'
#fe_names[2] = 'BC'
#fe_names[3] = 'BCCoarse'
#fe_names[4] = 'BC+STC'

###########

file_in_tau = {}

if algo == 'A':

  dir_in = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoA/DMid/'

  file_in_tau[0] = dir_in+'ntuple_TauGun_THRESHOLD_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  #file_in_tau[1] = dir_in+'ntuple_TauGun_SUPERTRIGGERCELL_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  #file_in_tau[2] = dir_in+'ntuple_TauGun_BESTCHOICE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  #file_in_tau[3] = dir_in+'ntuple_TauGun_BESTCHOICECOARSE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  #file_in_tau[4] = dir_in+'ntuple_TauGun_MIXEDBCSTC_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'

if algo == 'B':

  dir_in = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoB/DMid/'

  file_in_tau[0] = dir_in+'ntuple_TauGun_THRESHOLD_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  #file_in_tau[1] = dir_in+'ntuple_TauGun_SUPERTRIGGERCELL_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  #file_in_tau[2] = dir_in+'ntuple_TauGun_BESTCHOICE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  #file_in_tau[3] = dir_in+'ntuple_TauGun_BESTCHOICECOARSE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  #file_in_tau[4] = dir_in+'ntuple_TauGun_MIXEDBCSTC_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'

###########

file_out_mapping = {}

if algo == 'A':

  dir_out = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoA/mapping/'

  file_out_mapping[0] = dir_out+'mapping_95eff_THRESHOLD_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.pkl'
  #file_out_mapping[1] = dir_out+'mapping_95eff_SUPERTRIGGERCELL_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.pkl'
  #file_out_mapping[2] = dir_out+'mapping_95eff_BESTCHOICE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.pkl'
  #file_out_mapping[3] = dir_out+'mapping_95eff_BESTCHOICECOARSE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.pkl'
  #file_out_mapping[4] = dir_out+'mapping_95eff_MIXEDBCSTC_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.pkl'

if algo == 'B':

  dir_out = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoB/mapping/'

  file_out_mapping[0] = dir_out+'mapping_95eff_THRESHOLD_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.pkl'
  #file_out_mapping[1] = dir_out+'mapping_95eff_SUPERTRIGGERCELL_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.pkl'
  #file_out_mapping[2] = dir_out+'mapping_95eff_BESTCHOICE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.pkl'
  #file_out_mapping[3] = dir_out+'mapping_95eff_BESTCHOICECOARSE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.pkl'
  #file_out_mapping[4] = dir_out+'mapping_95eff_MIXEDBCSTC_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.pkl'

###########

if PUBDTWP == '99':
  bdtcut = 'cl3d_pubdt_passWP99'

elif PUBDTWP == '95':
  bdtcut = 'cl3d_pubdt_passWP95'

elif PUBDTWP == '90':
  bdtcut = 'cl3d_pubdt_passWP90'

###########

df_tau = {}
df_tau_DM01 = {}
df_tau_DM45 = {}

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

  seldm01 = (df_tau[name]['gentau_decayMode'] == 0) | (df_tau[name]['gentau_decayMode'] == 1)
  df_tau_DM01[name] = df_tau[name][seldm01]

  seldm45 = (df_tau[name]['gentau_decayMode'] == 2)
  df_tau_DM45[name] = df_tau[name][seldm45]

###########

if algo == 'A':
  plotdir = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoA/plots/turnons/'

elif algo == 'B':
  plotdir = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoB/plots/turnons/'

ptcut = 1
etamin = 1.6

for name in df_tau:

  df_tau[name]['gentau_vis_abseta'] = np.abs(df_tau[name]['gentau_vis_eta'])
  df_tau[name]['gentau_bin_eta'] = ((df_tau[name]['gentau_vis_abseta'] - etamin)/0.1).astype('int32')
  df_tau[name]['gentau_bin_pt']  = ((df_tau[name]['gentau_vis_pt'] - ptcut)/2).astype('int32')

  df_tau_DM01[name]['gentau_vis_abseta'] = np.abs(df_tau_DM01[name]['gentau_vis_eta'])
  df_tau_DM01[name]['gentau_bin_eta'] = ((df_tau_DM01[name]['gentau_vis_abseta'] - etamin)/0.1).astype('int32')
  df_tau_DM01[name]['gentau_bin_pt']  = ((df_tau_DM01[name]['gentau_vis_pt'] - ptcut)/2).astype('int32')

  df_tau_DM45[name]['gentau_vis_abseta'] = np.abs(df_tau_DM45[name]['gentau_vis_eta'])
  df_tau_DM45[name]['gentau_bin_eta'] = ((df_tau_DM45[name]['gentau_vis_abseta'] - etamin)/0.1).astype('int32')
  df_tau_DM45[name]['gentau_bin_pt']  = ((df_tau_DM45[name]['gentau_vis_pt'] - ptcut)/2).astype('int32')

def efficiency(group, threshold):

  tot = group.shape[0]
  sel = group[(group.cl3d_pt_c3 > threshold)].shape[0]
  return float(sel)/float(tot)

efficiencies_vs_pt = {}
efficiencies_vs_pt_DM01 = {}
efficiencies_vs_pt_DM45 = {}

for name in df_tau:

  efficiencies_vs_pt[name] = {}
  efficiencies_vs_pt[name] = df_tau[name].groupby('gentau_bin_pt').mean()

  efficiencies_vs_pt_DM01[name] = {}
  efficiencies_vs_pt_DM01[name] = df_tau_DM01[name].groupby('gentau_bin_pt').mean()

  efficiencies_vs_pt_DM45[name] = {}
  efficiencies_vs_pt_DM45[name] = df_tau_DM45[name].groupby('gentau_bin_pt').mean()

turnon_thresholds = range(1, 100, 1)
#turnon_thresholds = range(50, 60, 10)

for name in df_tau:

  for threshold in turnon_thresholds:

    eff = df_tau[name].groupby('gentau_bin_pt').apply(lambda x : efficiency(x, threshold))
    eff_DM01 = df_tau_DM01[name].groupby('gentau_bin_pt').apply(lambda x : efficiency(x, threshold))
    eff_DM45 = df_tau_DM45[name].groupby('gentau_bin_pt').apply(lambda x : efficiency(x, threshold))

    eff_smooth = eff.rolling(window=7, win_type='triang', center=True).mean()
    eff_smooth_DM01 = eff_DM01.rolling(window=7, win_type='triang', center=True).mean()
    eff_smooth_DM45 = eff_DM45.rolling(window=7, win_type='triang', center=True).mean()

    eff_smooth.fillna(0., inplace=True)
    eff_smooth_DM01.fillna(0., inplace=True)
    eff_smooth_DM45.fillna(0., inplace=True)

    efficiencies_vs_pt[name]['efficiency_{}'.format(threshold)] = eff
    efficiencies_vs_pt_DM01[name]['efficiency_{}'.format(threshold)] = eff_DM01
    efficiencies_vs_pt_DM45[name]['efficiency_{}'.format(threshold)] = eff_DM45

    efficiencies_vs_pt[name]['efficiency_smooth_{}'.format(threshold)] = eff_smooth
    efficiencies_vs_pt_DM01[name]['efficiency_smooth_{}'.format(threshold)] = eff_smooth_DM01
    efficiencies_vs_pt_DM45[name]['efficiency_smooth_{}'.format(threshold)] = eff_smooth_DM45

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
#colors[1] = 'red'
#colors[2] = 'olive'
#colors[3] = 'orange'
#colors[4] = 'fuchsia'

legends = {}
legends[0] = 'Threshold 1.35 mipT'
#legends[1] = 'STC4+16'
#legends[2] = 'BC Decentral'
#legends[3] = 'BC Coarse 2x2 TC'
#legends[4] = 'Mixed BC + STC'

import matplotlib.lines as mlines

import geeksw.plotting.cmsplot as plt
from geeksw.plotting.root_colors import *

plt.matplotlib.font_manager._rebuild()

plt.rcParams['legend.numpoints'] = 1


plt.figure(figsize=(8,8))

lab_40 = r"$E_{T}^{L1,\tau}$ > 40 GeV"
lab_50 = r"$E_{T}^{L1,\tau}$ > 50 GeV"
lab_60 = r"$E_{T}^{L1,\tau}$ > 60 GeV"

for name in df_tau:
  df = efficiencies_vs_pt[name]
  y_eff_40 = df['efficiency_40']
  y_eff_50 = df['efficiency_50']
  y_eff_60 = df['efficiency_60']
  y_eff_smooth_40 = df['efficiency_smooth_40']
  y_eff_smooth_50 = df['efficiency_smooth_50']
  y_eff_smooth_60 = df['efficiency_smooth_60']
  x = df.gentau_vis_pt

plt.errorbar(x,y_eff_40,xerr=1,ls='None',label=lab_40,color='blue',lw=2,marker='o',mec='blue')
plt.errorbar(x,y_eff_50,xerr=1,ls='None',label=lab_50,color='green',lw=2,marker='o',mec='green')
plt.errorbar(x,y_eff_60,xerr=1,ls='None',label=lab_60,color='red',lw=2,marker='o',mec='red')
plt.plot(x,y_eff_smooth_40,label=lab_40,color='blue',lw=1.5)
plt.plot(x,y_eff_smooth_50,label=lab_50,color='green',lw=1.5)
plt.plot(x,y_eff_smooth_60,label=lab_60,color='red',lw=1.5)
plt.xlim(20, 91.5)
plt.ylim(0., 1.10)
blue_line = mlines.Line2D([], [], color='blue',markersize=15, label=lab_40,lw=2)
green_line = mlines.Line2D([], [], color='green',markersize=15, label=lab_50,lw=2)
red_line = mlines.Line2D([], [], color='red',markersize=15, label=lab_60,lw=2)
plt.legend(loc = 'lower right', fontsize=16, handles=[blue_line,green_line,red_line])
#plt.legend(loc = 'lower right', fontsize=16)
plt.xlabel(r'$p_{T}^{gen,\tau}\ [GeV]$')
plt.ylabel('Efficiency')
plt.grid()
plt.subplots_adjust(bottom=0.12)
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("PU=200","HGCAL")
plt.savefig(plotdir+'eff_vs_pt_TDRv1.png')
plt.savefig(plotdir+'eff_vs_pt_TDRv1.pdf')


import matplotlib.lines as mlines
plt.figure(figsize=(8,8))

lab_50 = r"$E_{T}^{L1,\tau}$ > 50 GeV"

for name in df_tau:

  df_DM01 = efficiencies_vs_pt_DM01[name]
  df_DM45 = efficiencies_vs_pt_DM45[name]

  eff_DM01 = df_DM01['efficiency_50']
  eff_DM45 = df_DM45['efficiency_50']

  eff_smooth_DM01 = df_DM01['efficiency_smooth_50']
  eff_smooth_DM45 = df_DM45['efficiency_smooth_50']

  x_DM01 = efficiencies_vs_pt_DM01[name].gentau_vis_pt
  x_DM45 = efficiencies_vs_pt_DM45[name].gentau_vis_pt

  plt.errorbar(x_DM01,eff_DM01,xerr=1,ls='None',label='1-prong (+ $\pi^{0}$\'s)',color='red',lw=2,marker='o',mec='red')
  plt.errorbar(x_DM45,eff_DM45,xerr=1,ls='None',label='3-prongs (+ $\pi^{0}$\'s)',color='blue',lw=2,marker='o',mec='blue')

  plt.plot(x_DM01,eff_smooth_DM01,color='red',lw=1.5)
  plt.plot(x_DM45,eff_smooth_DM45,color='blue',lw=1.5)

lab_50 = r"$E_{T}^{L1,\tau}$ > 50 GeV"

plt.xlim(20, 91.5)
plt.ylim(0., 1.10)

red_line = mlines.Line2D([], [], color='red',markersize=15, label='1-prong (+ $\pi^{0}$\'s)',lw=2)
blue_line = mlines.Line2D([], [], color='blue',markersize=15, label='3-prong (+ $\pi^{0}$\'s)',lw=2)

plt.legend(loc = 'lower right', fontsize=18, handles=[red_line,blue_line])
txt = (r'Gen. $\tau$ decay mode:')
t = plt.text(63,0.20, txt, ha='left', wrap=True, fontsize=18)
t.set_bbox(dict(facecolor='white', edgecolor='white'))
txt2 = (r'$E_{T}^{L1,\tau}$ > 50 GeV')
t2 = plt.text(26,0.83, txt2, ha='left', wrap=True, fontsize=18)
t2.set_bbox(dict(facecolor='white', edgecolor='white'))
plt.xlabel(r'$p_{T}^{gen,\tau}\ [GeV]$')
plt.ylabel('Efficiency')
plt.grid()
plt.subplots_adjust(bottom=0.12)
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("PU=200","HGCAL")
plt.savefig(plotdir+'eff_vs_pt_DM_TDRv1.png')
plt.savefig(plotdir+'eff_vs_pt_DM_TDRv1.pdf')


'''
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
'''

'''
plt.figure(figsize=(8,8))
for name in df_tau:
  df = pt_95s[name]
  plt.plot(df.threshold, df.pt95, label=legends[name], linewidth=2, color=colors[name])
#plt.legend(loc = 'lower right', fontsize=16)
plt.xlabel('L1 Threshold [GeV]')
plt.ylabel('Offline threshold [GeV]')
plt.xlim(10, 80)
plt.ylim(10, 100)
plt.grid()
plt.subplots_adjust(bottom=0.12)
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("PU=200","HGCAL")
plt.savefig(plotdir+'L1_to_offline_WP_algo'+algo+'_'+PUBDTWP+'.png')
plt.savefig(plotdir+'L1_to_offline_WP_algo'+algo+'_'+PUBDTWP+'.pdf')
plt.show()
'''