# C. Martin Perez cmartinp@cern.ch, Nov. 2019

###########

import os
import pandas as pd
import pickle
import numpy as np
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import LinearRegression
from matplotlib import pyplot as plt
import matplotlib

###########

fe_names = {}

fe_names[0] = 'Threshold'
#fe_names[1] = 'STC'
#fe_names[2] = 'BestChoice'
#fe_names[3] = 'BestChoiceCoarse'
#fe_names[4] = 'MixedBC+STC'

###########

dir_in = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/matched/'

file_in_tau = {}

file_in_tau[0] = dir_in+'ntuple_TauGun_THRESHOLD_matched.hdf5'
#file_in_tau[1] = dir_in+'ntuple_TauGun_SUPERTRIGGERCELL_matched.hdf5'
#file_in_tau[2] = dir_in+'ntuple_TauGun_BESTCHOICE_matched.hdf5'
#file_in_tau[3] = dir_in+'ntuple_TauGun_BESTCHOICECOARSE_matched.hdf5'
#file_in_tau[4] = dir_in+'ntuple_TauGun_MIXEDBCSTC_matched.hdf5'

file_in_nu = {}

file_in_nu[0] = dir_in+'ntuple_NuGun_THRESHOLD_matched.hdf5'
#file_in_nu[1] = dir_in+'ntuple_NuGun_SUPERTRIGGERCELL_matched.hdf5'
#file_in_nu[2] = dir_in+'ntuple_NuGun_BESTCHOICE_matched.hdf5'
#file_in_nu[3] = dir_in+'ntuple_NuGun_BESTCHOICECOARSE_matched.hdf5'
#file_in_nu[4] = dir_in+'ntuple_NuGun_MIXEDBCSTC_matched.hdf5'

###########

dir_out = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoB/calibrated/'

file_out_tau = {}

file_out_tau[0] = dir_out+'ntuple_TauGun_THRESHOLD_calibrated_algoB.hdf5'
#file_out_tau[1] = dir_out+'ntuple_TauGun_SUPERTRIGGERCELL_calibrated_algoB.hdf5'
#file_out_tau[2] = dir_out+'ntuple_TauGun_BESTCHOICE_calibrated_algoB.hdf5'
#file_out_tau[3] = dir_out+'ntuple_TauGun_BESTCHOICECOARSE_calibrated_algoB.hdf5'
#file_out_tau[4] = dir_out+'ntuple_TauGun_MIXEDBCSTC_calibrated_algoB.hdf5'

file_out_nu = {}

file_out_nu[0] = dir_out+'ntuple_NuGun_THRESHOLD_calibrated_algoB.hdf5'
#file_out_nu[1] = dir_out+'ntuple_NuGun_SUPERTRIGGERCELL_calibrated_algoB.hdf5'
#file_out_nu[2] = dir_out+'ntuple_NuGun_BESTCHOICE_calibrated_algoB.hdf5'
#file_out_nu[3] = dir_out+'ntuple_NuGun_BESTCHOICECOARSE_calibrated_algoB.hdf5'
#file_out_nu[4] = dir_out+'ntuple_NuGun_MIXEDBCSTC_calibrated_algoB.hdf5'

###########

dir_out_model = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoB/models/calibration/'

file_out_model_c1 = {}

file_out_model_c1[0] = dir_out_model+'model_c1_algoB_threshold.pkl'
#file_out_model_c1[1] = dir_out_model+'model_c1_algoB_supertriggercell.pkl'
#file_out_model_c1[2] = dir_out_model+'model_c1_algoB_bestchoice.pkl'
#file_out_model_c1[3] = dir_out_model+'model_c1_algoB_bestchoicecoarse.pkl'
#file_out_model_c1[4] = dir_out_model+'model_c1_algoB_mixedbcstc.pkl'

file_out_model_c2 = {}

file_out_model_c2[0] = dir_out_model+'model_c2_algoB_threshold.pkl'
#file_out_model_c2[1] = dir_out_model+'model_c2_algoB_supertriggercell.pkl'
#file_out_model_c2[2] = dir_out_model+'model_c2_algoB_bestchoice.pkl'
#file_out_model_c2[3] = dir_out_model+'model_c2_algoB_bestchoicecoarse.pkl'
#file_out_model_c2[4] = dir_out_model+'model_c2_algoB_mixedbcstc.pkl'

file_out_model_c3 = {}

file_out_model_c3[0] = dir_out_model+'model_c3_algoB_threshold.pkl'
#file_out_model_c3[1] = dir_out_model+'model_c3_algoB_supertriggercell.pkl'
#file_out_model_c3[2] = dir_out_model+'model_c3_algoB_bestchoice.pkl'
#file_out_model_c3[3] = dir_out_model+'model_c3_algoB_bestchoicecoarse.pkl'
#file_out_model_c3[4] = dir_out_model+'model_c3_algoB_mixedbcstc.pkl'

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

# Selection for training

df_tau_train = {}

for name in df_tau:

  df_tau[name]['cl3d_abseta'] = np.abs(df_tau[name]['cl3d_eta'])
  df_tau[name]['cl3d_response'] = df_tau[name]['cl3d_pt']/df_tau[name]['gentau_vis_pt']

  df_tau_train[name] = df_tau[name]

  sel = df_tau_train[name]['gentau_vis_pt'] > 20
  df_tau_train[name] = df_tau_train[name][sel]
  
  sel = np.abs(df_tau_train[name]['gentau_vis_eta']) > 1.6
  df_tau_train[name] = df_tau_train[name][sel]
  
  sel = np.abs(df_tau_train[name]['gentau_vis_eta']) < 2.9
  df_tau_train[name] = df_tau_train[name][sel]
  
  sel = df_tau_train[name]['cl3d_isbestmatch'] == True
  df_tau_train[name] = df_tau_train[name][sel]

  sel = df_tau_train[name]['cl3d_pt'] > 4
  df_tau_train[name] = df_tau_train[name][sel]

for name in df_nu:

  df_nu[name]['cl3d_abseta'] = np.abs(df_nu[name]['cl3d_eta'])
  df_nu[name]['n_matched_cl3d'] = 1

###########

# Training calibration 1
print 'Training calibration 1...'

model_c1 = {}

for name in df_tau_train:

  input_c1 = df_tau_train[name][['cl3d_abseta']]
  target_c1 = df_tau_train[name].gentau_vis_pt - df_tau_train[name].cl3d_pt
  model_c1[name] = LinearRegression().fit(input_c1, target_c1)

for name in df_tau_train:

  with open(file_out_model_c1[name], 'wb') as f:
    pickle.dump(model_c1[name], f)

for name in df_tau_train:

  df_tau_train[name]['cl3d_c1'] = model_c1[name].predict(df_tau_train[name][['cl3d_abseta']])
  df_tau_train[name]['cl3d_pt_c1'] = df_tau_train[name].cl3d_c1 + df_tau_train[name].cl3d_pt

###########

# Training calibration 2
print 'Training calibration 2...'

features = ['n_matched_cl3d', 'cl3d_abseta', 
  'cl3d_showerlength', 'cl3d_coreshowerlength', 
  'cl3d_firstlayer', 'cl3d_maxlayer', 
  'cl3d_szz', 'cl3d_seetot', 'cl3d_spptot', 'cl3d_srrtot', 'cl3d_srrmean',
  'cl3d_hoe', 'cl3d_meanz', 
  'cl3d_layer10', 'cl3d_layer50', 'cl3d_layer90', 
  'cl3d_ntc67', 'cl3d_ntc90']

model_c2 = {}

for name in df_tau_train:

  input_c2 = df_tau_train[name][features]
  target_c2 = df_tau_train[name].gentau_vis_pt / df_tau_train[name].cl3d_pt_c1
  model_c2[name] = GradientBoostingRegressor(n_estimators=1000, learning_rate=0.1, max_depth=2, random_state=0, loss='huber').fit(input_c2, target_c2)

for name in df_tau_train:

  with open(file_out_model_c2[name], 'wb') as f:
    pickle.dump(model_c2[name], f)

for name in df_tau_train:

  df_tau_train[name]['cl3d_c2'] = model_c2[name].predict(df_tau_train[name][features])
  df_tau_train[name]['cl3d_pt_c2'] = df_tau_train[name].cl3d_c2 * df_tau_train[name].cl3d_pt_c1
  df_tau_train[name]['cl3d_response_c2'] = df_tau_train[name].cl3d_pt_c2 / df_tau_train[name].gentau_vis_pt

'''
for name in df_tau:

  feature_importance = model_c2[name].feature_importances_
  sorted_idx = np.argsort(feature_importance)
  pos = np.arange(sorted_idx.shape[0]) + .5
  fig = plt.figure(figsize=(12, 6))
  plt.subplot(1, 2, 1)
  plt.barh(pos, feature_importance[sorted_idx], align='center')
  plt.yticks(pos, np.array(features)[sorted_idx])
  plt.title('Feature Importance (MDI)')
  plt.show()
'''


###########

# Training calibration 3
print 'Training calibration 3...'

ptmin = 20
etamin = 1.6

for name in df_tau_train:

  df_tau_train[name]['gentau_vis_abseta'] = np.abs(df_tau_train[name]['gentau_vis_eta'])
  df_tau_train[name]['gentau_vis_bin_eta'] = ((df_tau_train[name]['gentau_vis_abseta'] - etamin)/0.1).astype('int32')
  df_tau_train[name]['gentau_vis_bin_pt']  = ((df_tau_train[name]['gentau_vis_pt'] - ptmin)/5).astype('int32')

vars = ['gentau_vis_pt', 'gentau_vis_bin_pt', 'cl3d_pt_c2', 'cl3d_response_c2']

pt_means_train = {}
pt_rmss_train = {}

for name in df_tau_train:

  pt_means_train[name] = df_tau_train[name][vars].groupby('gentau_vis_bin_pt').mean() 
  pt_rmss_train[name] = df_tau_train[name][vars].groupby('gentau_vis_bin_pt').std() 

for name in df_tau_train:

  pt_means_train[name]['logpt1'] = np.log(pt_means_train[name]['cl3d_pt_c2'])
  pt_means_train[name]['logpt2'] = pt_means_train[name].logpt1**2
  pt_means_train[name]['logpt3'] = pt_means_train[name].logpt1**3
  pt_means_train[name]['logpt4'] = pt_means_train[name].logpt1**4

model_c3 = {}

for name in df_tau_train:
  input_c3 = pt_means_train[name][['logpt1', 'logpt2', 'logpt3', 'logpt4']]
  target_c3 = pt_means_train[name]['cl3d_response_c2']
  model_c3[name] = LinearRegression().fit(input_c3, target_c3)

for name in df_tau_train:

  with open(file_out_model_c3[name], 'wb') as f:
    pickle.dump(model_c3[name], f)

for name in df_tau_train:

  logpt1 = np.log(abs(df_tau_train[name]['cl3d_pt_c2']))
  logpt2 = logpt1**2
  logpt3 = logpt1**3
  logpt4 = logpt1**4
  
  df_tau_train[name]['cl3d_c3'] = model_c3[name].predict(np.vstack([logpt1, logpt2, logpt3, logpt4]).T)
  df_tau_train[name]['cl3d_pt_c3'] = df_tau_train[name].cl3d_pt_c2 / df_tau_train[name].cl3d_c3


###########

# Application calibration 1

for name in df_tau:

  df_tau[name]['cl3d_c1'] = model_c1[name].predict(df_tau[name][['cl3d_abseta']])
  df_tau[name]['cl3d_pt_c1'] = df_tau[name].cl3d_c1 + df_tau[name].cl3d_pt
  df_tau[name]['cl3d_response_c1'] = df_tau[name].cl3d_pt_c1 / df_tau[name].gentau_vis_pt

for name in df_nu:

  df_nu[name]['cl3d_c1'] = model_c1[name].predict(df_nu[name][['cl3d_abseta']])
  df_nu[name]['cl3d_pt_c1'] = df_nu[name].cl3d_c1 + df_nu[name].cl3d_pt

###########

# Application calibration 2

for name in df_tau:

  df_tau[name]['cl3d_c2'] = model_c2[name].predict(df_tau[name][features])
  df_tau[name]['cl3d_pt_c2'] = df_tau[name].cl3d_c2 * df_tau[name].cl3d_pt_c1
  df_tau[name]['cl3d_response_c2'] = df_tau[name].cl3d_pt_c2 / df_tau[name].gentau_vis_pt

for name in df_nu:

  df_nu[name]['cl3d_c2'] = model_c2[name].predict(df_nu[name][features])
  df_nu[name]['cl3d_pt_c2'] = df_nu[name].cl3d_c2 * df_nu[name].cl3d_pt_c1

###########

# Application calibration 3

for name in df_tau:

  logpt1 = np.log(abs(df_tau[name]['cl3d_pt_c2']))
  logpt2 = logpt1**2
  logpt3 = logpt1**3
  logpt4 = logpt1**4
  
  df_tau[name]['cl3d_c3'] = model_c3[name].predict(np.vstack([logpt1, logpt2, logpt3, logpt4]).T)
  df_tau[name]['cl3d_pt_c3'] = df_tau[name].cl3d_pt_c2 / df_tau[name].cl3d_c3
  df_tau[name]['cl3d_response_c3'] = df_tau[name].cl3d_pt_c3 / df_tau[name].gentau_vis_pt

for name in df_nu:

  logpt1 = np.log(abs(df_nu[name]['cl3d_pt_c2']))
  logpt2 = logpt1**2
  logpt3 = logpt1**3
  logpt4 = logpt1**4

  df_nu[name]['cl3d_c3'] = model_c3[name].predict(np.vstack([logpt1, logpt2, logpt3, logpt4]).T)
  df_nu[name]['cl3d_pt_c3'] = df_nu[name].cl3d_pt_c2 / df_nu[name].cl3d_c3

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

# BINNED RESPONSES AND RESOLUTIONS

def effrms(df, c=0.68):
  """Compute half-width of the shortest interval
  containing a fraction 'c' of items in a 1D array.
  """
  out = {}
  for col in df:
      x = df[col]
      x = np.sort(x, kind="mergesort")
      m = int(c * len(x)) + 1
      out[col] = [np.min(x[m:] - x[:-m]) / 2.0]
  return pd.DataFrame(out).iloc[0]

###########

# Selection for plotting

plot_var = ['gentau_vis_pt', 'gentau_vis_abseta',
          'gentau_vis_bin_eta', 'gentau_vis_bin_pt',
          'cl3d_pt', 'cl3d_response', 'cl3d_abseta',
          'cl3d_pt_c1', 'cl3d_response_c1',
          'cl3d_pt_c2', 'cl3d_response_c2',
          'cl3d_pt_c3', 'cl3d_response_c3']

etameans = {}
etarmss = {}
etaeffrmss = {}
ptmeans = {}
ptrmss = {}
pteffrmss = {}

for name in df_tau:

  sel = df_tau[name]['gentau_vis_pt'] > 20
  df_tau[name] = df_tau[name][sel]
  
  sel = np.abs(df_tau[name]['gentau_vis_eta']) > 1.6
  df_tau[name] = df_tau[name][sel]
  
  sel = np.abs(df_tau[name]['gentau_vis_eta']) < 2.9
  df_tau[name] = df_tau[name][sel]
  
  sel = df_tau[name]['cl3d_isbestmatch'] == True
  df_tau[name] = df_tau[name][sel]

  sel = df_tau[name]['cl3d_pt'] > 4
  df_tau[name] = df_tau[name][sel]

  df_tau[name]['gentau_vis_abseta'] = np.abs(df_tau[name]['gentau_vis_eta'])

  df_tau[name]['gentau_vis_bin_eta'] = ((df_tau[name]['gentau_vis_abseta'] - etamin)/0.1).astype('int32')
  df_tau[name]['gentau_vis_bin_pt']  = ((df_tau[name]['gentau_vis_pt'] - ptmin)/5).astype('int32')

  etameans[name]    = df_tau[name][plot_var].groupby('gentau_vis_bin_eta').mean()
  etarmss[name]     = df_tau[name][plot_var].groupby('gentau_vis_bin_eta').std()
  etaeffrmss[name]  = df_tau[name][plot_var].groupby('gentau_vis_bin_eta').apply(effrms)
  ptmeans[name]     = df_tau[name][plot_var].groupby('gentau_vis_bin_pt').mean()
  ptrmss[name]      = df_tau[name][plot_var].groupby('gentau_vis_bin_pt').std()
  pteffrmss[name]   = df_tau[name][plot_var].groupby('gentau_vis_bin_pt').apply(effrms)

###########

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

plotdir = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoB/plots/calibration/'

matplotlib.rcParams.update({'font.size': 22})

###########

# Print values

'''
print 'RAW: '

for name in df_tau:

  print ' ',fe_names[name],(': mean={0}, rms={1}, rms/mean={2}'.format(
      df_tau[name]['cl3d_response'].mean(),
      df_tau[name]['cl3d_response'].std(),
      df_tau[name]['cl3d_response'].std()/df_tau[name]['cl3d_response'].mean()
  ))

print ' '

print 'CALIBRATED 1: '

for name in df_tau:

  print ' ',fe_names[name],(': mean={0}, rms={1}, rms/mean={2}'.format(
      df_tau[name]['cl3d_response_c1'].mean(),
      df_tau[name]['cl3d_response_c1'].std(),
      df_tau[name]['cl3d_response_c1'].std()/df_tau[name]['cl3d_response_c1'].mean()
  ))

print ' '

print 'CALIBRATED 2: '

for name in df_tau:

  print ' ',fe_names[name],(': mean={0}, rms={1}, rms/mean={2}'.format(
      df_tau[name]['cl3d_response_c2'].mean(),
      df_tau[name]['cl3d_response_c2'].std(),
      df_tau[name]['cl3d_response_c2'].std()/df_tau[name]['cl3d_response_c2'].mean()
  ))

print ' '

print 'CALIBRATED 3: '

for name in df_tau:

  print ' ',fe_names[name], (': mean={0}, rms={1}, rms/mean={2}'.format(
      df_tau[name]['cl3d_response_c3'].mean(),
      df_tau[name]['cl3d_response_c3'].std(),
      df_tau[name]['cl3d_response_c3'].std()/df_tau[name]['cl3d_response_c3'].mean()
  ))
print ' '
'''
###########

# Inclusive responses

import geeksw.plotting.cmsplot as plt
from geeksw.plotting.root_colors import *

plt.matplotlib.font_manager._rebuild()

import matplotlib.lines as mlines


for name in df_tau:

  plt.figure(figsize=(8,8))
  plt.hist(df_tau[name]['cl3d_response'],    bins=np.arange(0., 2., 0.03),  label='Uncalibrated, mean: 0.69, RMS: 0.22',      color='red',    histtype='step', lw=2)
  plt.hist(df_tau[name]['cl3d_response_c1'], bins=np.arange(0., 2., 0.03),  label='Calib. 1, mean: 1.05, RMS: 0.23', color='blue',    histtype='step', lw=2)
  plt.hist(df_tau[name]['cl3d_response_c2'], bins=np.arange(0., 2., 0.03),  label='Calib. 2, mean: 1.02, RMS: 0.17', color='limegreen',  histtype='step', lw=2)
  plt.hist(df_tau[name]['cl3d_response_c3'], bins=np.arange(0., 2., 0.03),  label='Calib. 3, mean: 1.00, RMS: 0.18', color='black',   histtype='step', lw=2)
  red_line = mlines.Line2D([], [], color='red',markersize=15, label='Uncalib., mean: 0.69, RMS: 0.22',lw=2)
  blue_line = mlines.Line2D([], [], color='blue',markersize=15, label='Calib. 1, mean: 1.05, RMS: 0.23',lw=2)
  green_line = mlines.Line2D([], [], color='limegreen',markersize=15, label='Calib. 2, mean: 1.02, RMS: 0.17',lw=2)
  black_line = mlines.Line2D([], [], color='black',markersize=15, label='Calib. 3, mean: 1.00, RMS: 0.18',lw=2)
  plt.legend(loc = 'upper right',handles=[red_line,blue_line,green_line,black_line],fontsize=15)
  #plt.legend(loc = 'upper right', fontsize=22)
  plt.grid()
  plt.xlabel(r'$E_{T}^{L1,\tau}\ /\ p_{T}^{gen,\tau}$')
  plt.ylabel(r'a. u.')
  plt.ylim(0, 1800)
  plt.cmstext("CMS"," Phase-2 Simulation")
  plt.lumitext("PU=200","HGCAL")
  plt.subplots_adjust(bottom=0.12)
  txt = (r'1.6 < | $\eta_{gen,\tau}$ | < 2.9')
  t = plt.text(70,0.32, txt, ha='left', wrap=True, fontsize=16)
  t.set_bbox(dict(facecolor='white', edgecolor='white'))
  plt.savefig(plotdir+'calib_responses_TDRv2.png')
  plt.savefig(plotdir+'calib_responses_TDRv2.pdf')
  plt.show()


'''
myvariable = 'cl3d_spptot'
myxmin = 0
myxmax = 0.08
mybinsize = 0.005

df_tau_DM0 = {}
df_tau_DM1 = {}
df_tau_DM4 = {}
df_tau_DM5 = {}

for name in df_tau:

  seldm0 = (df_tau[name]['gentau_decayMode'] == 0)
  df_tau_DM0[name] = df_tau[name][seldm0]

  seldm1 = (df_tau[name]['gentau_decayMode'] == 1)
  df_tau_DM1[name] = df_tau[name][seldm1]

  seldm4 = (df_tau[name]['gentau_decayMode'] == 4)
  df_tau_DM4[name] = df_tau[name][seldm4]

  seldm5 = (df_tau[name]['gentau_decayMode'] == 5)
  df_tau_DM5[name] = df_tau[name][seldm5]

  plt.figure(figsize=(8,8))
  #plt.hist(df_tau[name]['cl3d_response_c2'], bins=np.arange(0., 2., 0.03),  label='Calib. 2', color='limegreen',  histtype='step', lw=2)
  plt.hist(df_tau_DM0[name][myvariable], normed=True, bins=np.arange(myxmin,myxmax,mybinsize),  label='1-prong', color='limegreen',  histtype='step', lw=2)
  plt.hist(df_tau_DM1[name][myvariable], normed=True, bins=np.arange(myxmin,myxmax,mybinsize),  label='1-prong + $\pi^{0}$\'s', color='blue',  histtype='step', lw=2)
  plt.hist(df_tau_DM4[name][myvariable], normed=True, bins=np.arange(myxmin,myxmax,mybinsize),  label='3-prongs', color='red',  histtype='step', lw=2)
  plt.hist(df_tau_DM5[name][myvariable], normed=True, bins=np.arange(myxmin,myxmax,mybinsize),  label='3-prongs + $\pi^{0}$\'s', color='fuchsia',  histtype='step', lw=2)
  green_line = mlines.Line2D([], [], color='limegreen',markersize=15, label='1-prong',lw=2)
  blue_line = mlines.Line2D([], [], color='blue',markersize=15, label='1-prong + $\pi^{0}$\'s',lw=2)
  red_line = mlines.Line2D([], [], color='red',markersize=15, label='3-prongs',lw=2)
  fuchsia_line = mlines.Line2D([], [], color='fuchsia',markersize=15, label='3-prongs + $\pi^{0}$\'s', lw=2)
  plt.legend(loc = 'upper right',handles=[green_line,blue_line,red_line,fuchsia_line])
  #plt.legend(loc = 'upper right', fontsize=18)
  plt.grid()
  plt.xlabel(r'Total $\sigma_{\phi\phi}$')
  plt.ylabel(r'Normalized events')
  plt.ylim(0, 50)
  plt.cmstext("CMS"," Phase-2 Simulation")
  plt.lumitext("PU=200","HGCAL")
  plt.subplots_adjust(bottom=0.12)
  #txt = (r'1.6 < | $\eta_{gen,\tau}$ | < 2.9')
  #t = plt.text(70,0.32, txt, ha='left', wrap=True, fontsize=16)
  #t.set_bbox(dict(facecolor='white', edgecolor='white'))
  plt.savefig(plotdir+'calib2_'+myvariable+'_DM_TDRv1.png')
  plt.savefig(plotdir+'calib2_'+myvariable+'_DM_TDRv1.pdf')
  plt.show()
'''

###########

# Model c3
'''
matplotlib.rcParams.update({'font.size': 22})
plt.figure(figsize=(8,8))
pt = np.arange(25, 100, 1)

logpt1 = np.log(pt)
logpt2 = logpt1**2
logpt3 = logpt1**3
logpt4 = logpt1**4

for name in df_tau:

  ptmeans[name]['logpt1'] = np.log(ptmeans[name]['cl3d_pt_c2'])
  ptmeans[name]['logpt2'] = ptmeans[name].logpt1**2
  ptmeans[name]['logpt3'] = ptmeans[name].logpt1**3
  ptmeans[name]['logpt4'] = ptmeans[name].logpt1**4

for name in df_tau:

  plt.plot(ptmeans[name]['cl3d_pt_c2'], ptmeans[name]['cl3d_response_c2'], marker='s', markersize=9, ls='None', color=colors[name], label='Observed')
  plt.plot(pt, model_c3[name].predict(np.vstack([logpt1, logpt2, logpt3, logpt4]).T), ls='--', color=colors[name], label='Predicted')

plt.grid()
plt.legend(loc = 'upper right', fontsize=16)
plt.xlabel(r'$\langle E_{T,calib2}^{L1}\rangle\, [GeV]$')
plt.ylabel(r'$\langle E_{T,calib2}^{L1}/p_{T}^{gen}\rangle$')
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("PU=200","HGCAL")
plt.subplots_adjust(bottom=0.12)
txt = (r'1.6 < | $\eta_{gen,\tau}$ | < 2.9')
t = plt.text(70,0.32, txt, ha='left', wrap=True, fontsize=16)
t.set_bbox(dict(facecolor='white', edgecolor='white'))
plt.savefig(plotdir+'calib_modelc3_TDRv2.png')
plt.savefig(plotdir+'calib_modelc3_TDRv2.pdf')
#plt.show()
'''

###########

# Responses and resolutions
'''
plt.figure(figsize=(15,10))

for name in df_tau:

  df = etameans[name]
  plt.plot(df['gentau_vis_abseta'], etaeffrmss[name]['cl3d_response']/df['cl3d_response'], label="Raw",color='red',lw=2)
  plt.plot(df['gentau_vis_abseta'], etaeffrmss[name]['cl3d_response_c3']/df['cl3d_response_c3'], label="Calibrated",color='blue',lw=2)
  y_array_raw = (etaeffrmss[name]['cl3d_response']/df['cl3d_response']).values
  y_array_calib = (etaeffrmss[name]['cl3d_response_c3']/df['cl3d_response_c3']).values

plt.errorbar(x_array,y_array_raw,xerr=2.5,yerr=None,marker="o",mec='red',ls='None',label='Raw',color='red',lw=2)
plt.errorbar(x_array,y_array_calib,xerr=2.5,yerr=None,marker="o",mec='blue',ls='None',label='Calibrated',color='blue',lw=2)
plt.xlabel(r'$|\eta|$')
plt.ylabel(r'$RMS_{eff}(E_{T}^{L1}/p_{T}^{gen})$')
plt.grid()
plt.legend(loc = 'upper right', fontsize=16)
plt.savefig(plotdir+'calib_mean_vs_eta.png')
plt.savefig(plotdir+'calib_mean_vs_eta.pdf')
'''

'''
matplotlib.rcParams.update({'font.size': 22})

fig, axs = plt.subplots(2, 1, figsize=(15,15))

for name in df_tau:

  df = etameans[name]
  axs[0].plot(df['gentau_vis_abseta'], etarmss[name]['cl3d_response_c3'],color=colors[name], label=legends[name],lw=2)
  ratio = etarmss[name]['cl3d_response_c3']/etarmss[0]['cl3d_response_c3']
  axs[1].plot(df['gentau_vis_abseta'], ratio, color=colors[name], label=legends[name],lw=2)

axs[0].legend(loc = 'upper right', fontsize=18)
axs[0].set_xlim(1.6, 3.0)
axs[0].set_ylim(0.10,0.35)
axs[0].grid()
#axs[0].set_xlabel(r'$|\eta|$')
axs[0].set_ylabel(r'$RMS(p_{T}^{L1}/p_{T}^{gen})$')
axs[1].set_xlim(1.6, 3.0)
axs[1].set_ylim(0.95,1.20)
axs[1].set_xlabel(r'$|\eta|$')
axs[1].set_ylabel('Ratio to Threshold')
axs[1].grid()
plt.savefig(plotdir+'calib_rms_vs_eta_ratio.png')
plt.savefig(plotdir+'calib_rms_vs_eta_ratio.pdf')
'''

'''
matplotlib.rcParams.update({'font.size': 22})

fig, axs = plt.subplots(2, 1, figsize=(15,15))

for name in df_tau:

  df = etameans[name]
  axs[0].plot(df['gentau_vis_abseta'], etaeffrmss[name]['cl3d_response_c3'],color=colors[name], label=legends[name],lw=2)
  ratio = etaeffrmss[name]['cl3d_response_c3']/etaeffrmss[0]['cl3d_response_c3']
  axs[1].plot(df['gentau_vis_abseta'], ratio, color=colors[name], label=legends[name],lw=2)

axs[0].legend(loc = 'upper right', fontsize=18)
axs[0].set_xlim(1.6, 3.0)
axs[0].set_ylim(0.10,0.30)
axs[0].grid()
#axs[0].set_xlabel(r'$|\eta|$')
axs[0].set_ylabel(r'$RMS_{eff}(p_{T}^{L1}/p_{T}^{gen})$')
axs[1].set_xlim(1.6, 3.0)
axs[1].set_ylim(0.95,1.25)
axs[1].set_xlabel(r'$|\eta|$')
axs[1].set_ylabel('Ratio to Threshold')
axs[1].grid()
plt.savefig(plotdir+'calib_effrms_vs_eta_ratio.png')
plt.savefig(plotdir+'calib_effrms_vs_eta_ratio.pdf')
'''

'''
plt.figure(figsize=(15,10))

for name in df_tau:

  df = ptmeans[name]
  plt.plot(df['gentau_vis_pt'], df['cl3d_response_c3'],color=colors[name], label=legends[name],lw=2)

plt.ylim(0.9, 1.1)
plt.xlabel(r'$p_{T}\, [GeV]$', fontsize=18)
plt.ylabel(r'$\langle p_{T}^{L1}/p_{T}^{gen}\rangle$')
plt.grid()
plt.legend(loc = 'upper right', fontsize=16)
plt.ylim(0.95, 1.05)
plt.savefig(plotdir+'calib_mean_vs_pt.png')
plt.savefig(plotdir+'calib_mean_vs_pt.pdf')
'''

'''
plt.figure(figsize=(15,10))

for name in df_tau:

  df = ptmeans[name]
  plt.plot(df['gentau_vis_pt'], ptrmss[name]['cl3d_response_c3'],color=colors[name], label=legends[name],lw=2)

plt.ylim(0., 0.3)
plt.xlabel(r'$p_{T}\, [GeV]$')
plt.ylabel(r'$RMS(p_{T}^{L1}/p_{T}^{gen})$')
plt.grid()
plt.legend(loc = 'upper right', fontsize=16)
plt.ylim(0.0, 0.4)
plt.savefig(plotdir+'calib_rms_vs_pt.png')
plt.savefig(plotdir+'calib_rms_vs_pt.pdf')
'''
'''
matplotlib.rcParams.update({'font.size': 22})

fig, axs = plt.subplots(2, 1, figsize=(15,15))

for name in df_tau:

  df = ptmeans[name]
  axs[0].plot(df['gentau_vis_pt'], ptrmss[name]['cl3d_response_c3'],color=colors[name], label=legends[name],lw=2)
  ratio = ptrmss[name]['cl3d_response_c3']/ptrmss[0]['cl3d_response_c3']
  axs[1].plot(df['gentau_vis_pt'], ratio, color=colors[name], label=legends[name],lw=2)

axs[0].legend(loc = 'upper right', fontsize=18)
axs[0].set_xlim(20., 100.)
axs[0].set_ylim(0.0, 0.4)
axs[0].grid()
#axs[0].set_xlabel(r'$RMS(p_{T}^{L1}/p_{T}^{gen})$')
axs[0].set_ylabel(r'$RMS(p_{T}^{L1}/p_{T}^{gen})$')
axs[1].set_xlim(20., 100.)
axs[1].set_ylim(0.95,1.40)
axs[1].set_xlabel(r'$p_{T}\, [GeV]$')
axs[1].set_ylabel('Ratio to Threshold')
axs[1].grid()
plt.savefig(plotdir+'calib_rms_vs_pt_ratio.png')
plt.savefig(plotdir+'calib_rms_vs_pt_ratio.pdf')
'''
'''
plt.figure(figsize=(8,9))

for name in df_tau:

  df = ptmeans[name]
  plt.plot(df['gentau_vis_pt'], pteffrmss[name]['cl3d_response']/df['cl3d_response'], label="Raw",color='red',lw=2)
  plt.plot(df['gentau_vis_pt'], pteffrmss[name]['cl3d_response_c3']/df['cl3d_response_c3'], label="Calibrated",color='blue',lw=2)

plt.ylim(0.01, 0.45)
plt.legend(loc = 'lower left')
plt.xlabel(r'$p_{T}^{gen,\tau}\ [GeV]$')
plt.ylabel(r'$RMS_{eff}\ /\ <p_{T}^{L1,\tau}\ /\ p_{T}^{gen,\tau}>$')
plt.grid()
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("PU=200","HGCAL")
plt.subplots_adjust(bottom=0.12)
plt.subplots_adjust(left=0.14)
#plt.ylim(0.0, 0.5)
plt.savefig(plotdir+'calib_effrms_vs_pt_TDR.png')
plt.savefig(plotdir+'calib_effrms_vs_pt_TDR.pdf')
plt.show()
'''

'''
x_array = [22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 82.5, 87.5, 92.5, 97.5]

for name in df_tau:
  df = ptmeans[name]
  y_array_raw = (pteffrmss[name]['cl3d_response']/df['cl3d_response']).values
  y_array_calib = (pteffrmss[name]['cl3d_response_c3']/df['cl3d_response_c3']).values

plt.figure(figsize=(8,8))
plt.errorbar(x_array,y_array_raw,xerr=2.5,yerr=None,marker="o",mec='red',ls='None',label='Raw',color='red',lw=2)
plt.errorbar(x_array,y_array_calib,xerr=2.5,yerr=None,marker="o",mec='blue',ls='None',label='Calibrated',color='blue',lw=2)
plt.ylim(0.01, 0.45)
red_line = mlines.Line2D([], [], color='red',markersize=15, label='Raw',lw=2)
blue_line = mlines.Line2D([], [], color='blue',markersize=15, label='Calibrated',lw=2)
plt.legend(loc = 'lower left',handles=[red_line,blue_line])
plt.xlabel(r'$p_{T}^{gen,\tau}\ [GeV]$')
plt.ylabel(r'$RMS_{eff}\ /\ <E_{T}^{L1,\tau}\ /\ p_{T}^{gen,\tau}>$')
plt.grid()
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("PU=200","HGCAL")
plt.subplots_adjust(bottom=0.12)
plt.subplots_adjust(left=0.14)
#txt = (r'1.6 < | $\eta_{gen,\tau}$ | < 2.9')
#t = plt.text(70,0.32, txt, ha='left', wrap=True, fontsize=16)
#t.set_bbox(dict(facecolor='white', edgecolor='white'))
#plt.ylim(0.0, 0.5)
plt.savefig(plotdir+'calib_effrms_vs_pt_TDRv1.png')
plt.savefig(plotdir+'calib_effrms_vs_pt_TDRv1.pdf')
plt.show()
'''
'''
x_array = [1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85]

for name in df_tau:
  df = etameans[name]
  y_array_raw = (etaeffrmss[name]['cl3d_response']/df['cl3d_response']).values
  y_array_calib = (etaeffrmss[name]['cl3d_response_c3']/df['cl3d_response_c3']).values

plt.figure(figsize=(8,8))
plt.errorbar(x_array,y_array_raw,xerr=0.05,yerr=None,marker="o",mec='red',ls='None',label='Raw',color='red',lw=2)
plt.errorbar(x_array,y_array_calib,xerr=0.05,yerr=None,marker="o",mec='blue',ls='None',label='Calibrated',color='blue',lw=2)
red_line = mlines.Line2D([], [], color='red',markersize=15, label='Raw',lw=2)
blue_line = mlines.Line2D([], [], color='blue',markersize=15, label='Calibrated',lw=2)
plt.legend(loc = 'lower left',handles=[red_line,blue_line])
plt.xlabel(r'$|\eta^{gen,\tau}|$')
plt.ylabel(r'$RMS_{eff}\ /\ <E_{T}^{L1,\tau}\ /\ p_{T}^{gen,\tau}>$')
plt.grid()
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("PU=200","HGCAL")
plt.subplots_adjust(bottom=0.12)
plt.subplots_adjust(left=0.14)
#txt = (r'1.6 < | $\eta_{gen,\tau}$ | < 2.9')
#t = plt.text(70,0.32, txt, ha='left', wrap=True, fontsize=16)
#t.set_bbox(dict(facecolor='white', edgecolor='white'))
plt.savefig(plotdir+'calib_effrms_vs_eta_TDRv1.png')
plt.savefig(plotdir+'calib_effrms_vs_eta_TDRv1.pdf')
plt.show()
'''
'''
matplotlib.rcParams.update({'font.size': 22})

fig, axs = plt.subplots(2, 1, figsize=(15,15))

for name in df_tau:

  df = ptmeans[name]
  axs[0].plot(df['gentau_vis_pt'], pteffrmss[name]['cl3d_response_c3'],color=colors[name], label=legends[name],lw=2)
  ratio = pteffrmss[name]['cl3d_response_c3']/pteffrmss[0]['cl3d_response_c3']
  axs[1].plot(df['gentau_vis_pt'], ratio, color=colors[name], label=legends[name],lw=2)

axs[0].legend(loc = 'upper right', fontsize=18)
axs[0].set_xlim(20., 100.)
axs[0].set_ylim(0., 0.35)
axs[0].grid()
#axs[0].set_xlabel(r'$RMS(p_{T}^{L1}/p_{T}^{gen})$')
axs[0].set_ylabel(r'$RMS_{eff}(p_{T}^{L1}/p_{T}^{gen})$')
axs[1].set_xlim(20., 100.)
axs[1].set_ylim(0.95,1.40)
axs[1].set_xlabel(r'$p_{T}\, [GeV]$')
axs[1].set_ylabel('Ratio to Threshold')
axs[1].grid()
plt.savefig(plotdir+'calib_effrms_vs_pt_ratio.png')
plt.savefig(plotdir+'calib_effrms_vs_pt_ratio.pdf')
'''
'''
plt.figure(figsize=(8,7))

for name in df_tau:
  plt.hist(df_tau[name]['gentau_vis_eta']-df_tau[name]['cl3d_eta'], bins=np.arange(-0.15, 0.15, 0.005), color='blue',lw=2,histtype='step')

#plt.ylim(0.01, 0.4)
plt.xlabel(r'$\eta_{gen,\tau}$ - $\eta_{L1,\tau}$')
plt.ylabel(r'a. u. ')
#plt.grid()
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("PU=200","HGCAL")
plt.subplots_adjust(bottom=0.12)
plt.subplots_adjust(left=0.14)
plt.grid()
#txt = (r'1.6 < | $\eta_{gen,\tau}$ | < 2.9')
#t = plt.text(-0.13,1400, txt, ha='left', wrap=True, fontsize=16)
#t.set_bbox(dict(facecolor='white', edgecolor='white'))
#txt2 = ('Mean: 0.00, RMS: 0.02')
#t2 = plt.text(0.025,1700, txt2, ha='left', wrap=True, fontsize=16)
#t2.set_bbox(dict(facecolor='white', edgecolor='white'))
#plt.ylim(0.0, 0.5)
plt.savefig(plotdir+'res_eta_TDRv1.png')
plt.savefig(plotdir+'res_eta_TDRv1.pdf')

plt.figure(figsize=(8,7))

for name in df_tau:
  plt.hist(df_tau[name]['gentau_vis_phi']-df_tau[name]['cl3d_phi'], bins=np.arange(-0.15, 0.15, 0.005), color='blue',lw=2,histtype='step')

#plt.ylim(0.01, 0.4)
plt.xlabel(r'$\phi_{gen,\tau}$ - $\phi_{L1,\tau}$')
plt.ylabel(r'a. u. ')
#plt.grid()
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("PU=200","HGCAL")
plt.subplots_adjust(bottom=0.12)
plt.subplots_adjust(left=0.14)
plt.grid()
#txt = (r'1.6 < | $\eta_{gen,\tau}$ | < 2.9')
#t = plt.text(-0.13,1400, txt, ha='left', wrap=True, fontsize=16)
#t.set_bbox(dict(facecolor='white', edgecolor='white'))
#txt2 = ('Mean: 0.00, RMS: 0.02')
#t2 = plt.text(0.025,1700, txt2, ha='left', wrap=True, fontsize=16)
#t2.set_bbox(dict(facecolor='white', edgecolor='white'))
#plt.ylim(0.0, 0.5)
plt.savefig(plotdir+'res_phi_TDRv1.png')
plt.savefig(plotdir+'res_phi_TDRv1.pdf')
'''
