# C. Martin Perez cmartinp@cern.ch, Nov. 2019

###########

import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import xgboost as xgb
import matplotlib
import pickle
from sklearn import metrics
from sklearn.model_selection import train_test_split
import root_pandas

###########

def train_xgb(df, inputs, output, hyperparams, test_fraction=0.4):

    X_train, X_test, y_train, y_test = train_test_split(df[inputs], df[output], test_size=test_fraction)
    train = xgb.DMatrix(data=X_train,label=y_train, feature_names=inputs)
    test = xgb.DMatrix(data=X_test,label=y_test,feature_names=inputs)
    full = xgb.DMatrix(data=df[inputs],label=df[output],feature_names=inputs)
    booster = xgb.train(hyperparams, full, num_boost_round=hyperparams['num_trees'])
    df['bdt_output'] = booster.predict(full)
    fpr, tpr, threshold = metrics.roc_curve(df[output], df['bdt_output'])

    return booster, fpr, tpr, threshold

def efficiency(group, cut):
    tot = group.shape[0]
    sel = group[group.cl3d_pubdt_score > cut].shape[0]
    return float(sel)/float(tot)

###########

fe_names = {}

fe_names[0] = 'Thr'
fe_names[1] = 'STC'
fe_names[2] = 'BC'
fe_names[3] = 'BCCoarse'
fe_names[4] = 'BC+STC'

###########

dir_in = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoB/calibrated/'

file_in_tau = {}

file_in_tau[0] = dir_in+'ntuple_TauGun_THRESHOLD_calibrated_algoB.hdf5'
file_in_tau[1] = dir_in+'ntuple_TauGun_SUPERTRIGGERCELL_calibrated_algoB.hdf5'
file_in_tau[2] = dir_in+'ntuple_TauGun_BESTCHOICE_calibrated_algoB.hdf5'
file_in_tau[3] = dir_in+'ntuple_TauGun_BESTCHOICECOARSE_calibrated_algoB.hdf5'
file_in_tau[4] = dir_in+'ntuple_TauGun_MIXEDBCSTC_calibrated_algoB.hdf5'

file_in_nu = {}

file_in_nu[0] = dir_in+'ntuple_NuGun_THRESHOLD_calibrated_algoB.hdf5'
file_in_nu[1] = dir_in+'ntuple_NuGun_SUPERTRIGGERCELL_calibrated_algoB.hdf5'
file_in_nu[2] = dir_in+'ntuple_NuGun_BESTCHOICE_calibrated_algoB.hdf5'
file_in_nu[3] = dir_in+'ntuple_NuGun_BESTCHOICECOARSE_calibrated_algoB.hdf5'
file_in_nu[4] = dir_in+'ntuple_NuGun_MIXEDBCSTC_calibrated_algoB.hdf5'

###########

dir_out = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoB/pubdt/'

file_out_tau = {}

file_out_tau[0] = dir_out+'ntuple_TauGun_THRESHOLD_calibrated_pubdt_algoB.hdf5'
file_out_tau[1] = dir_out+'ntuple_TauGun_SUPERTRIGGERCELL_calibrated_pubdt_algoB.hdf5'
file_out_tau[2] = dir_out+'ntuple_TauGun_BESTCHOICE_calibrated_pubdt_algoB.hdf5'
file_out_tau[3] = dir_out+'ntuple_TauGun_BESTCHOICECOARSE_calibrated_pubdt_algoB.hdf5'
file_out_tau[4] = dir_out+'ntuple_TauGun_MIXEDBCSTC_calibrated_pubdt_algoB.hdf5'

file_out_nu = {}

file_out_nu[0] = dir_out+'ntuple_NuGun_THRESHOLD_calibrated_pubdt_algoB.hdf5'
file_out_nu[1] = dir_out+'ntuple_NuGun_SUPERTRIGGERCELL_calibrated_pubdt_algoB.hdf5'
file_out_nu[2] = dir_out+'ntuple_NuGun_BESTCHOICE_calibrated_pubdt_algoB.hdf5'
file_out_nu[3] = dir_out+'ntuple_NuGun_BESTCHOICECOARSE_calibrated_pubdt_algoB.hdf5'
file_out_nu[4] = dir_out+'ntuple_NuGun_MIXEDBCSTC_calibrated_pubdt_algoB.hdf5'

###########

dir_out_model = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoB/models/pubdt/'

file_out_model = {}

file_out_model[0] = dir_out_model+'model_pubdt_algoB_threshold.pkl'
file_out_model[1] = dir_out_model+'model_pubdt_algoB_supertriggercell.pkl'
file_out_model[2] = dir_out_model+'model_pubdt_algoB_bestchoice.pkl'
file_out_model[3] = dir_out_model+'model_pubdt_algoB_bestchoicecoarse.pkl'
file_out_model[4] = dir_out_model+'model_pubdt_algoB_mixedbcstc.pkl'

###########

plotdir = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoB/plots/pubdt/'

###########

df_tau = {}

for name in file_in_tau:

  store_tau = pd.HDFStore(file_in_tau[name], mode='r')
  df_tau[name] = store_tau['df_tau_PU200'] 
  store_tau.close()

df_nu = {}

for name in file_in_nu:
  
  store_nu = pd.HDFStore(file_in_nu[name], mode='r')
  df_nu[name] = store_nu['df_nu_PU200']
  store_nu.close()

###########

df_merged_train = {}

df_tau_train = {}
df_nu_train = {}

events_total_tau = {}
events_stored_tau = {}

events_total_nu = {}
events_stored_nu = {}

for name in df_tau:

  dfs = []

  # SIGNAL

  events_total_tau[name] = np.unique(df_tau[name].reset_index()['event']).shape[0]

  df_tau_train[name] = df_tau[name]

  df_tau_train[name]['gentau_pid'] = 1
  
  sel = df_tau_train[name]['gentau_vis_pt'] > 20
  df_tau_train[name] = df_tau_train[name][sel]
  
  sel = np.abs(df_tau_train[name]['gentau_vis_eta']) > 1.6
  df_tau_train[name] = df_tau_train[name][sel]
  
  sel = np.abs(df_tau_train[name]['gentau_vis_eta']) < 2.9
  df_tau_train[name] = df_tau_train[name][sel]
  
  sel = df_tau_train[name]['cl3d_isbestmatch'] == True
  df_tau_train[name] = df_tau_train[name][sel]

  sel = df_tau_train[name]['cl3d_pt_c3'] > 4
  df_tau_train[name] = df_tau_train[name][sel]

  events_stored_tau[name] = np.unique(df_tau_train[name].reset_index()['event']).shape[0]
  
  dfs.append(df_tau_train[name])
  
  # BACKGROUND

  events_total_nu[name] = np.unique(df_nu[name].reset_index()['event']).shape[0]

  df_nu_train[name] = df_nu[name]
  
  df_nu_train[name]['gentau_pid'] = 0
  
  sel = df_nu_train[name]['cl3d_pt_c3'] > 20
  df_nu_train[name] = df_nu_train[name][sel]

  events_stored_nu[name] = np.unique(df_nu_train[name].reset_index()['event']).shape[0]
  
  dfs.append(df_nu_train[name])
  
  # MERGE
  
  df_merged_train[name] = pd.concat(dfs)

  print ' ',fe_names[name],' #cl3d bkg: ', df_nu_train[name].shape[0], ', #cl3d signal: ', df_tau_train[name].shape[0]

print ' '


###########

# PU CLUSTERS PER EVENT

nclusters_per_event = {}

for name in df_merged_train:

  sel = df_merged_train[name]['gentau_pid']==0
  df_pu = df_merged_train[name][sel]
  clusters_per_event = df_pu.groupby('event').count()
  nclusters_per_event[name] = np.mean(clusters_per_event.gentau_pid)

for name, nclusters in nclusters_per_event.items():

  nclusters_per_event[name] = nclusters*events_stored_nu[name]/events_total_nu[name]
  print ' ', fe_names[name],' #events stored ',events_stored_nu[name], ', #events total ',events_total_nu[name], ', clusters/event: ',nclusters_per_event[name]

print ' '

###########

# Training

for name in df_merged_train:

  df_merged_train[name]['cl3d_abseta'] = np.abs(df_merged_train[name].cl3d_eta)

inputs = ['cl3d_abseta','cl3d_showerlength','cl3d_coreshowerlength',
   'cl3d_firstlayer', 'cl3d_maxlayer',
   'cl3d_szz', 'cl3d_seetot', 'cl3d_spptot', 'cl3d_srrtot', 'cl3d_srrmean',
   'cl3d_hoe', 'cl3d_meanz',
   'cl3d_layer10', 'cl3d_layer50', 'cl3d_layer90',
   'cl3d_ntc67', 'cl3d_ntc90']

output = 'gentau_pid'

param = {}

param['nthread']          	= 10  # limit number of threads
param['eta']              	= 0.2 # learning rate
param['max_depth']        	= 4  # maximum depth of a tree
param['subsample']        	= 0.8 # fraction of events to train tree on
param['colsample_bytree'] 	= 0.8 # fraction of features to train tree on
param['silent'] 			      = True
param['objective']   		    = 'binary:logistic' # objective function
param['num_trees'] 			    = 81  # number of trees to make

model = {}
fpr = {}
fpr_ncluster = {}
tpr = {}
threshold = {}

for name in df_merged_train:

  model[name], fpr[name], tpr[name], threshold[name] = train_xgb(df_merged_train[name], inputs, output, param, test_fraction=0.4)
  fpr_ncluster[name] = fpr[name]*nclusters_per_event[name]

for name in df_merged_train:

  with open(file_out_model[name], 'wb') as f:
    pickle.dump(model[name], f)

bdt_thr_WP99 = {}
bdt_thr_WP95 = {}
bdt_thr_WP90 = {}

for name in df_merged_train:

  full = xgb.DMatrix(data=df_merged_train[name][inputs], label=df_merged_train[name][output], feature_names=inputs)
  df_merged_train[name]['cl3d_pubdt_score'] = model[name].predict(full)

  bdt_thr_WP99[name] = np.interp(0.99, tpr[name], threshold[name])
  bdt_thr_WP95[name] = np.interp(0.95, tpr[name], threshold[name])
  bdt_thr_WP90[name] = np.interp(0.90, tpr[name], threshold[name])

  print ' ',fe_names[name], ' threshold WP99 ',bdt_thr_WP99[name],', threshold WP95 ',bdt_thr_WP95[name],', threshold WP90 ',bdt_thr_WP90[name]

print ' '

###########

# Application

for name in df_tau:

  full = xgb.DMatrix(data=df_tau[name][inputs], label=df_tau[name][output], feature_names=inputs)
  df_tau[name]['cl3d_pubdt_score'] = model[name].predict(full)

  df_tau[name]['cl3d_pubdt_passWP99'] = df_tau[name]['cl3d_pubdt_score'] > bdt_thr_WP99[name]
  df_tau[name]['cl3d_pubdt_passWP95'] = df_tau[name]['cl3d_pubdt_score'] > bdt_thr_WP95[name]
  df_tau[name]['cl3d_pubdt_passWP90'] = df_tau[name]['cl3d_pubdt_score'] > bdt_thr_WP90[name]

for name in df_nu:

  full = xgb.DMatrix(data=df_nu[name][inputs], label=df_nu[name][output], feature_names=inputs)
  df_nu[name]['cl3d_pubdt_score'] = model[name].predict(full)

  df_nu[name]['cl3d_pubdt_passWP99'] = df_nu[name]['cl3d_pubdt_score'] > bdt_thr_WP99[name]
  df_nu[name]['cl3d_pubdt_passWP95'] = df_nu[name]['cl3d_pubdt_score'] > bdt_thr_WP95[name]
  df_nu[name]['cl3d_pubdt_passWP90'] = df_nu[name]['cl3d_pubdt_score'] > bdt_thr_WP90[name]

###########

# Save files

for name in df_tau:

  store_tau = pd.HDFStore(file_out_tau[name], mode='w')
  store_tau['df_tau_PU200'] = df_tau[name]
  store_tau.close()

for name in df_nu:

  store_nu = pd.HDFStore(file_out_nu[name], mode='w')
  store_nu['df_nu_PU200'] = df_nu[name]
  store_nu.close()

###########

for name in df_nu:

  #print np.unique(df_nu[name].reset_index()['event']).shape[0]

  sel = df_nu[name]['cl3d_pt_c3'] > 20
  df_nu[name] = df_nu[name][sel]

  #print np.unique(df_nu[name].reset_index()['event']).shape[0]

  sel = df_nu[name]['cl3d_pubdt_passWP99'] == True
  df_nu[name] = df_nu[name][sel]

  #print np.unique(df_nu[name].reset_index()['event']).shape[0]
  #print '---'


###########

import geeksw.plotting.cmsplot as plt
from geeksw.plotting.root_colors import *

plt.matplotlib.font_manager._rebuild()

import matplotlib.lines as mlines


# PLOTTING

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

# ROC
'''
matplotlib.rcParams.update({'font.size': 22})

plt.figure(figsize=(15,10))

for name in df_tau:

  roc_auc = metrics.auc(fpr[name], tpr[name])
  print ' ',fe_names[name], 'AUC ',roc_auc
  plt.plot(tpr[name],fpr[name],label=legends[name], color=colors[name],lw=2)

plt.grid(which='both')
plt.legend(loc = 'upper left', fontsize=22)
plt.xlim(0.9,1.001)
plt.yscale('log')
plt.ylim(0.0,1)
plt.xlabel('Signal efficiency')
plt.ylabel('Background efficiency') 
plt.savefig(plotdir+'pubdt_roc.png')
plt.savefig(plotdir+'pubdt_roc.pdf')
'''

plt.figure(figsize=(8,8))

for name in df_tau:

  roc_auc = metrics.auc(fpr_ncluster[name], tpr[name])
  plt.plot(tpr[name],fpr_ncluster[name],label=legends[name], color=colors[name],lw=2)

plt.grid(which='both')
plt.legend(loc = 'upper left', fontsize=18)
plt.xlim(0.9,1.001)
plt.yscale('log')
plt.ylim(0.04,0.5)
plt.xlabel('Signal efficiency')
plt.ylabel('PU clusters / event') 
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("PU=200"," ")
plt.subplots_adjust(bottom=0.12)
plt.savefig(plotdir+'pubdt_roc_nclusters_FE.png')
plt.savefig(plotdir+'pubdt_roc_nclusters_FE.pdf')

# FEATURE IMPORTANCES
'''
matplotlib.rcParams.update({'font.size': 16})

for name in df_tau:

  plt.figure(figsize=(15,10))
  xgb.plot_importance(model[name], grid=False, importance_type='gain',lw=2)
  plt.subplots_adjust(left=0.28, right=0.85, top=0.9, bottom=0.1)
  plt.savefig(plotdir+'pubdt_importances_'+fe_names[name]+'.png')
  plt.savefig(plotdir+'pubdt_importances_'+fe_names[name]+'.pdf')

for name in df_tau:

  plt.figure(figsize=(12,10))
  plt.hist(df_tau[name]['cl3d_pubdt_score'], bins=np.arange(-0.2, 1.2, 0.02), normed=True, color='red', histtype='step', lw=2, label='Tau PU=200')
  plt.hist(df_nu[name]['cl3d_pubdt_score'],  bins=np.arange(-0.2, 1.2, 0.02), normed=True, color='blue', histtype='step', lw=2, label='Nu PU=200')
  plt.legend(loc = 'upper right', fontsize=22)
  plt.xlabel(r'PU BDT score')
  plt.ylabel(r'Entries')
  plt.savefig(plotdir+'pubdt_scores_'+fe_names[name]+'.png')
  plt.savefig(plotdir+'pubdt_scores_'+fe_names[name]+'.pdf')
'''

###########

# EFFICIENCIES
'''
matplotlib.rcParams.update({'font.size': 22})

# Cuts for plotting

for name in df_tau:

  df_tau[name]['gentau_vis_abseta'] = np.abs(df_tau[name]['gentau_vis_eta'])

  sel = df_tau[name]['gentau_vis_pt'] > 20
  df_tau[name] = df_tau[name][sel]
  
  sel = np.abs(df_tau[name]['gentau_vis_eta']) > 1.6
  df_tau[name] = df_tau[name][sel]
  
  sel = np.abs(df_tau[name]['gentau_vis_eta']) < 2.9
  df_tau[name] = df_tau[name][sel]
  
  sel = df_tau[name]['cl3d_isbestmatch'] == True
  df_tau[name] = df_tau[name][sel]

  sel = df_tau[name]['cl3d_pt_c3'] > 4
  df_tau[name] = df_tau[name][sel]

ptmin = 20
etamin = 1.6

for name in df_tau:

  df_tau[name]['gentau_bin_eta'] = ((df_tau[name]['gentau_vis_abseta'] - etamin)/0.1).astype('int32')
  df_tau[name]['gentau_bin_pt']  = ((df_tau[name]['gentau_vis_pt'] - ptmin)/5).astype('int32')

# EFFICIENCY VS ETA

efficiencies_vs_eta = {}

for name in df_tau:

  efficiencies_vs_eta[name] = df_tau[name].groupby('gentau_bin_eta').mean()

for name in df_tau:

  efficiencies_vs_eta[name]['efficiency99'] = df_tau[name].groupby('gentau_bin_eta').apply(lambda x : efficiency(x, bdt_thr_WP99[name]))
  efficiencies_vs_eta[name]['efficiency95'] = df_tau[name].groupby('gentau_bin_eta').apply(lambda x : efficiency(x, bdt_thr_WP95[name]))
  efficiencies_vs_eta[name]['efficiency90'] = df_tau[name].groupby('gentau_bin_eta').apply(lambda x : efficiency(x, bdt_thr_WP90[name]))

plt.figure(figsize=(15,10))

for name in df_tau:

  df = efficiencies_vs_eta[name]
  plt.plot(df.gentau_vis_abseta, df.efficiency99, label=legends[name], color=colors[name],lw=2)

plt.ylim(0.8, 1.01)
plt.legend(loc = 'lower left', fontsize=16)
plt.xlabel(r'$|\eta|$')
plt.ylabel('Efficiency')
plt.grid()
plt.savefig(plotdir+'pubdt_eff99_eta.png')
plt.savefig(plotdir+'pubdt_eff99_eta.pdf')

plt.figure(figsize=(15,10))

for name in df_tau:

  df = efficiencies_vs_eta[name]
  plt.plot(df.gentau_vis_abseta, df.efficiency95, label=legends[name], color=colors[name],lw=2)

plt.ylim(0.8, 1.01)
plt.legend(loc = 'lower left', fontsize=16)
plt.xlabel(r'$|\eta|$')
plt.ylabel('Efficiency')
plt.grid()
plt.savefig(plotdir+'pubdt_eff95_eta.png')
plt.savefig(plotdir+'pubdt_eff95_eta.pdf')

plt.figure(figsize=(15,10))

for name in df_tau:

  df = efficiencies_vs_eta[name]
  plt.plot(df.gentau_vis_abseta, df.efficiency90, label=legends[name], color=colors[name],lw=2)

plt.ylim(0.7, 1.01)
plt.legend(loc = 'lower left', fontsize=16)
plt.xlabel(r'$|\eta|$')
plt.ylabel('Efficiency')
plt.grid()
plt.savefig(plotdir+'pubdt_eff90_eta.png')
plt.savefig(plotdir+'pubdt_eff90_eta.pdf')

# EFFICIENCY VS PT

efficiencies_vs_pt = {}

for name in df_tau:

  efficiencies_vs_pt[name] = df_tau[name].groupby('gentau_bin_pt').mean()

for name in df_tau:

  efficiencies_vs_pt[name]['efficiency99'] = df_tau[name].groupby('gentau_bin_pt').apply(lambda x : efficiency(x, bdt_thr_WP99[name]))
  efficiencies_vs_pt[name]['efficiency95'] = df_tau[name].groupby('gentau_bin_pt').apply(lambda x : efficiency(x, bdt_thr_WP95[name]))
  efficiencies_vs_pt[name]['efficiency90'] = df_tau[name].groupby('gentau_bin_pt').apply(lambda x : efficiency(x, bdt_thr_WP90[name]))

plt.figure(figsize=(15,10))

for name in df_tau:

  df = efficiencies_vs_pt[name]
  plt.plot(df.gentau_vis_pt, df.efficiency99, label=legends[name], color=colors[name],lw=2)

plt.ylim(0.75, 1.01)
plt.legend(loc = 'lower right', fontsize=16)
plt.xlabel(r'$p_{T}\,[GeV]$')
plt.ylabel('Efficiency')
plt.grid()
plt.savefig(plotdir+'pubdt_eff99_pt.png')
plt.savefig(plotdir+'pubdt_eff99_pt.pdf')

plt.figure(figsize=(15,10))

for name in df_tau:

  df = efficiencies_vs_pt[name]
  plt.plot(df.gentau_vis_pt, df.efficiency95, label=legends[name], color=colors[name],lw=2)

plt.ylim(0.75, 1.01)
plt.legend(loc = 'lower right', fontsize=16)
plt.xlabel(r'$p_{T}\,[GeV]$')
plt.ylabel('Efficiency')
plt.grid()
plt.savefig(plotdir+'pubdt_eff95_pt.png')
plt.savefig(plotdir+'pubdt_eff95_pt.pdf')

plt.figure(figsize=(15,10))

for name in df_tau:

  df = efficiencies_vs_pt[name]
  plt.plot(df.gentau_vis_pt, df.efficiency90, label=legends[name], color=colors[name],lw=2)

plt.ylim(0.7, 1.01)
plt.legend(loc = 'lower right', fontsize=16)
plt.xlabel(r'$p_{T}\,[GeV]$')
plt.ylabel('Efficiency')
plt.grid()
plt.savefig(plotdir+'pubdt_eff90_pt.png')
plt.savefig(plotdir+'pubdt_eff90_pt.pdf')
'''
