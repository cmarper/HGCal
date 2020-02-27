# C. Martin Perez cmartinp@cern.ch, Sep. 2019

###########

import os
from glob import glob
import itertools
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import xgboost as xgb
import matplotlib
import pickle
from sklearn import metrics
from sklearn.model_selection import train_test_split

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

dir_in = '/data_CMS/cms/mperez/HGCal_data/Sep19/clustered/'

filein_tau_PU200 = dir_in+'RelValDiTau_Pt20To100_v10_PU200_matched.hdf5'
filein_nu_PU200	 = dir_in+'Nu_E10_v10_PU200_matched.hdf5'

store_tau_PU200 = pd.HDFStore(filein_tau_PU200, mode='r')
df_tau_PU200	= store_tau_PU200['df_tau_PU200']
store_tau_PU200.close()

store_nu_PU200 	= pd.HDFStore(filein_nu_PU200, mode='r')
df_nu_PU200 	= store_nu_PU200['df_nu_PU200']
store_tau_PU200.close()

###########

dir_out = '/data_CMS/cms/mperez/HGCal_data/Sep19/pubdt/'

fileout_tau_PU200 = dir_out+'RelValDiTau_Pt20To100_v10_PU200_pubdt.hdf5'
fileout_nu_PU200  = dir_out+'Nu_E10_v10_PU200_pubdt.hdf5'

###########

dfs = []

ptcut 	= 20
etamin 	= 1.6
etamax 	= 2.9

# SIGNAL

df_tau_PU200['gentau_pid'] = 1

sel = df_tau_PU200['gentau_vis_pt'] > ptcut
df_tau_PU200 = df_tau_PU200[sel]

sel = np.abs(df_tau_PU200['gentau_vis_eta']) > etamin
df_tau_PU200 = df_tau_PU200[sel]

sel = np.abs(df_tau_PU200['gentau_vis_eta']) < etamax
df_tau_PU200 = df_tau_PU200[sel]

sel = df_tau_PU200['cl3d_isbestmatch'] == True
df_tau_PU200 = df_tau_PU200[sel]

dfs.append(df_tau_PU200)

# BACKGROUND

df_nu_PU200['gentau_pid'] = 0

sel = df_nu_PU200['cl3d_pt'] > ptcut
df_nu_PU200 = df_nu_PU200[sel]

#sel = np.abs(df_nu_PU200['cl3d_eta']) > etamin
#df_nu_PU200 = df_nu_PU200[sel]

#sel = np.abs(df_nu_PU200['cl3d_eta']) < etamax
#df_nu_PU200 = df_nu_PU200[sel]

dfs.append(df_nu_PU200)

# MERGE

df_merged = pd.concat(dfs)

print 'Bkg events: ', df_merged[df_merged.gentau_pid==0].shape[0], ', signal events =', df_merged[df_merged.gentau_pid==1].shape[0]

###########

# TRAINING

df_merged['cl3d_abseta'] = np.abs(df_merged.cl3d_eta)

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

model, fpr, tpr, threshold = train_xgb(df_merged, inputs, output, param, test_fraction=0.4)

# Save model
with open('/data_CMS/cms/mperez/HGCal_data/Sep19/models/model_pubdt.pkl', 'wb') as f:
  pickle.dump(model, f)

###########

# APPLICATION

full = xgb.DMatrix(data=df_merged[inputs], label=df_merged[output], feature_names=inputs)
df_merged['cl3d_pubdt_score'] = model.predict(full)

bdt_thr_WP99 = np.interp(0.99, tpr, threshold)
bdt_thr_WP95 = np.interp(0.95, tpr, threshold)
bdt_thr_WP90 = np.interp(0.90, tpr, threshold)

#print bdt_thr_WP99 #0.108599005491
#print bdt_thr_WP95 #0.44313165769
#print bdt_thr_WP90 #0.642932797968

df_merged['cl3d_pubdt_passWP99'] = df_merged['cl3d_pubdt_score'] > bdt_thr_WP99
df_merged['cl3d_pubdt_passWP95'] = df_merged['cl3d_pubdt_score'] > bdt_thr_WP95
df_merged['cl3d_pubdt_passWP90'] = df_merged['cl3d_pubdt_score'] > bdt_thr_WP90


###########

# SAVE DFS

# Separate tau and PU dfs

df_tau_PU200  = df_merged.query('gentau_pid==1')
df_nu_PU200  = df_merged.query('gentau_pid==0')

# Save files
store_tau_PU200 = pd.HDFStore(fileout_tau_PU200, mode='w')
store_tau_PU200['df_tau_PU200'] = df_tau_PU200
store_tau_PU200.close()

store_nu_PU200 = pd.HDFStore(fileout_nu_PU200, mode='w')
store_nu_PU200['df_nu_PU200'] = df_nu_PU200
store_nu_PU200.close()

###########

# EFFICIENCIES

df_tau_PU200['gentau_vis_abseta'] = np.abs(df_tau_PU200['gentau_vis_eta'])

ptcut=20
etamin=1.6
etamax=2.9

# eta

df_tau_PU200['gentau_bin_eta'] = ((df_tau_PU200['gentau_vis_abseta'] - etamin)/0.1).astype('int32')

efficiencies_vs_eta = df_tau_PU200.groupby('gentau_bin_eta').mean()

efficiencies_vs_eta['efficiency99'] = df_tau_PU200.groupby('gentau_bin_eta').apply(lambda x : efficiency(x, bdt_thr_WP99))
efficiencies_vs_eta['efficiency95'] = df_tau_PU200.groupby('gentau_bin_eta').apply(lambda x : efficiency(x, bdt_thr_WP95))
efficiencies_vs_eta['efficiency90'] = df_tau_PU200.groupby('gentau_bin_eta').apply(lambda x : efficiency(x, bdt_thr_WP90))

# pt

df_tau_PU200['gentau_bin_pt']  = ((df_tau_PU200['gentau_vis_pt'] - ptcut)/5).astype('int32')

efficiencies_vs_pt = df_tau_PU200.groupby('gentau_bin_pt').mean()

efficiencies_vs_pt['efficiency99'] = df_tau_PU200.groupby('gentau_bin_pt').apply(lambda x : efficiency(x, bdt_thr_WP99))
efficiencies_vs_pt['efficiency95'] = df_tau_PU200.groupby('gentau_bin_pt').apply(lambda x : efficiency(x, bdt_thr_WP95))
efficiencies_vs_pt['efficiency90'] = df_tau_PU200.groupby('gentau_bin_pt').apply(lambda x : efficiency(x, bdt_thr_WP90))

###########

# PLOTTING

plotdir = '/home/llr/cms/mperez/HGCal/v10_geometry/tau_algorithm/plots/'

import geeksw.plotting.cmsplot as plt
from geeksw.plotting.root_colors import *

plt.matplotlib.font_manager._rebuild()

plt.figure(figsize=(8,7))
plt.plot(tpr,fpr,lw=2)
plt.xlabel('Signal efficiency')
plt.ylabel('Background efficiency') 
plt.xlim(0.89,1.001)
plt.ylim(0.007,0.51)
plt.grid()
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("200 PU","HGCal calorimeter-only")
plt.savefig(plotdir+'pubdt_roc_TDR2.png')
plt.savefig(plotdir+'pubdt_roc_TDR2.pdf')
#plt.show()

# ROC

'''matplotlib.rcParams.update({'font.size': 22})
roc_auc = metrics.auc(fpr, tpr)
plt.figure(figsize=(20,10))
plt.plot(tpr,fpr)
plt.grid(which='both')
#plt.xlim(0.9,1.001)
#plt.yscale('log')
#plt.ylim(0.007,1)
plt.xlabel('Signal efficiency')
plt.ylabel('Background efficiency') 
plt.savefig(plotdir+'pubdt_roc.png')
plt.savefig(plotdir+'pubdt_roc.pdf')
'''

# FEATURE IMPORTANCES
'''
matplotlib.rcParams.update({'font.size': 13})
plt.figure(figsize=(20,10))
xgb.plot_importance(model, grid=False, importance_type='gain')
plt.subplots_adjust(left=0.28, right=0.85, top=0.9, bottom=0.1)
plt.savefig(plotdir+'pubdt_importances.png')
plt.savefig(plotdir+'pubdt_importances.pdf')
'''

# EFFICIENCY VS ETA
'''
matplotlib.rcParams.update({'font.size': 22})
plt.figure(figsize=(15,10))
df = efficiencies_vs_eta
plt.plot(df.gentau_vis_abseta, df.efficiency99, label='WP 99%', color='blue')
plt.plot(df.gentau_vis_abseta, df.efficiency95, label='WP 95%', color='red')
plt.plot(df.gentau_vis_abseta, df.efficiency90, label='WP 90%', color='olive')
plt.ylim(0.5, 1.01)
plt.legend(loc = 'lower left', fontsize=16)
plt.xlabel(r'$|\eta|$')
plt.ylabel('Efficiency')
plt.grid()
plt.savefig(plotdir+'pubdt_eff_eta.png')
plt.savefig(plotdir+'pubdt_eff_eta.pdf')

matplotlib.rcParams.update({'font.size': 22})
plt.figure(figsize=(15,10))
df = efficiencies_vs_pt
plt.plot(df.gentau_vis_pt, df.efficiency99, label='WP 99%', color='blue')
plt.plot(df.gentau_vis_pt, df.efficiency95, label='WP 95%', color='red')
plt.plot(df.gentau_vis_pt, df.efficiency90, label='WP 90%', color='olive')
plt.ylim(0.5, 1.01)
plt.legend(loc = 'lower left', fontsize=16)
plt.xlabel(r'$p_{T}\,[GeV]$')
plt.ylabel('Efficiency')
plt.grid()
plt.savefig(plotdir+'pubdt_eff_pt.png')
plt.savefig(plotdir+'pubdt_eff_pt.pdf')

matplotlib.rcParams.update({'font.size': 22})
plt.figure(figsize=(12,10))
plt.hist(df_tau_PU200['cl3d_pubdt_score'], bins=np.arange(-0.2, 1.2, 0.02), normed=True, color='red', histtype='step', lw=2, label='Tau PU=200')
plt.hist(df_nu_PU200['cl3d_pubdt_score'],  bins=np.arange(-0.2, 1.2, 0.02), normed=True, color='blue', histtype='step', lw=2, label='Nu PU=200')
plt.legend(loc = 'upper right', fontsize=22)
plt.xlabel(r'PU BDT score')
plt.ylabel(r'Entries')
plt.savefig(plotdir+'pubdt_scores.png')
plt.savefig(plotdir+'pubdt_scores.pdf')
'''
