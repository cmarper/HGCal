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

# INPUTS

print 'Algorithm A'

PUBDTWP = 99
#PUBDTWP = 95
#PUBDTWP = 90

if PUBDTWP == 99:
  bdtcut = 'cl3d_pubdt_passWP99'
  text = 'WP99'
  print 'PU BDT WP99'

elif PUBDTWP == 95:
  bdtcut = 'cl3d_pubdt_passWP95'
  text = 'WP95'
  print 'PU BDT WP95'

elif PUBDTWP == 90:
  bdtcut = 'cl3d_pubdt_passWP90'
  text = 'WP90'
  print 'PU BDT WP90'

else:
  print 'Error: No PU BDT WP chosen!'

print ' '

###########

fe_names = {}

fe_names[0] = 'Threshold'
fe_names[1] = 'STC'
fe_names[2] = 'BestChoice'
fe_names[3] = 'BestChoiceCoarse'
fe_names[4] = 'MixedBC+STC'

###########

dir_in = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoA/pubdt/'

file_in_tau = {}

file_in_tau[0] = dir_in+'ntuple_TauGun_THRESHOLD_pubdt_algoA.hdf5'
file_in_tau[1] = dir_in+'ntuple_TauGun_SUPERTRIGGERCELL_pubdt_algoA.hdf5'
file_in_tau[2] = dir_in+'ntuple_TauGun_BESTCHOICE_pubdt_algoA.hdf5'
file_in_tau[3] = dir_in+'ntuple_TauGun_BESTCHOICECOARSE_pubdt_algoA.hdf5'
file_in_tau[4] = dir_in+'ntuple_TauGun_MIXEDBCSTC_pubdt_algoA.hdf5'

file_in_nu = {}

file_in_nu[0] = dir_in+'ntuple_NuGun_THRESHOLD_pubdt_algoA.hdf5'
file_in_nu[1] = dir_in+'ntuple_NuGun_SUPERTRIGGERCELL_pubdt_algoA.hdf5'
file_in_nu[2] = dir_in+'ntuple_NuGun_BESTCHOICE_pubdt_algoA.hdf5'
file_in_nu[3] = dir_in+'ntuple_NuGun_BESTCHOICECOARSE_pubdt_algoA.hdf5'
file_in_nu[4] = dir_in+'ntuple_NuGun_MIXEDBCSTC_pubdt_algoA.hdf5'

###########

dir_out = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoA/calibrated/'

file_out_tau = {}

file_out_tau[0] = dir_out+'ntuple_TauGun_THRESHOLD_pubdt_calibrated'+text+'_algoA.hdf5'
file_out_tau[1] = dir_out+'ntuple_TauGun_SUPERTRIGGERCELL_pubdt_calibrated'+text+'_algoA.hdf5'
file_out_tau[2] = dir_out+'ntuple_TauGun_BESTCHOICE_pubdt_calibrated'+text+'_algoA.hdf5'
file_out_tau[3] = dir_out+'ntuple_TauGun_BESTCHOICECOARSE_pubdt_calibrated'+text+'_algoA.hdf5'
file_out_tau[4] = dir_out+'ntuple_TauGun_MIXEDBCSTC_pubdt_calibrated'+text+'_algoA.hdf5'

file_out_nu = {}

file_out_nu[0] = dir_out+'ntuple_NuGun_THRESHOLD_pubdt_calibrated'+text+'_algoA.hdf5'
file_out_nu[1] = dir_out+'ntuple_NuGun_SUPERTRIGGERCELL_pubdt_calibrated'+text+'_algoA.hdf5'
file_out_nu[2] = dir_out+'ntuple_NuGun_BESTCHOICE_pubdt_calibrated'+text+'_algoA.hdf5'
file_out_nu[3] = dir_out+'ntuple_NuGun_BESTCHOICECOARSE_pubdt_calibrated'+text+'_algoA.hdf5'
file_out_nu[4] = dir_out+'ntuple_NuGun_MIXEDBCSTC_pubdt_calibrated'+text+'_algoA.hdf5'

###########

dir_out_model = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoA/models/calibration/'

file_out_model_c1 = {}

file_out_model_c1[0] = dir_out_model+'model_c1_algoA_'+text+'_threshold.pkl'
file_out_model_c1[1] = dir_out_model+'model_c1_algoA_'+text+'_supertriggercell.pkl'
file_out_model_c1[2] = dir_out_model+'model_c1_algoA_'+text+'_bestchoice.pkl'
file_out_model_c1[3] = dir_out_model+'model_c1_algoA_'+text+'_bestchoicecoarse.pkl'
file_out_model_c1[4] = dir_out_model+'model_c1_algoA_'+text+'_mixedbcstc.pkl'

file_out_model_c2 = {}

file_out_model_c2[0] = dir_out_model+'model_c2_algoA_'+text+'_threshold.pkl'
file_out_model_c2[1] = dir_out_model+'model_c2_algoA_'+text+'_supertriggercell.pkl'
file_out_model_c2[2] = dir_out_model+'model_c2_algoA_'+text+'_bestchoice.pkl'
file_out_model_c2[3] = dir_out_model+'model_c2_algoA_'+text+'_bestchoicecoarse.pkl'
file_out_model_c2[4] = dir_out_model+'model_c2_algoA_'+text+'_mixedbcstc.pkl'

file_out_model_c3 = {}

file_out_model_c3[0] = dir_out_model+'model_c3_algoA_'+text+'_threshold.pkl'
file_out_model_c3[1] = dir_out_model+'model_c3_algoA_'+text+'_supertriggercell.pkl'
file_out_model_c3[2] = dir_out_model+'model_c3_algoA_'+text+'_bestchoice.pkl'
file_out_model_c3[3] = dir_out_model+'model_c3_algoA_'+text+'_bestchoicecoarse.pkl'
file_out_model_c3[4] = dir_out_model+'model_c3_algoA_'+text+'_mixedbcstc.pkl'

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

  sel = df_tau_train[name][bdtcut] == True
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

  sel = df_tau[name][bdtcut] == True
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

plotdir = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoA/plots/calibration/'

matplotlib.rcParams.update({'font.size': 22})

###########

# Print values

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

###########

# Inclusive responses

for name in df_tau:

  plt.figure(figsize=(15,10))
  plt.hist(df_tau[name]['cl3d_response'], bins=np.arange(0., 2., 0.02), alpha=0.3, label='Before calib', color='blue')
  plt.hist(df_tau[name]['cl3d_response_c1'], bins=np.arange(0., 2., 0.02), alpha=0.3, label='After $\eta$ corr', color='red')
  plt.hist(df_tau[name]['cl3d_response_c2'], bins=np.arange(0., 2., 0.02), alpha=0.3, label='After BDT', color='green')
  plt.hist(df_tau[name]['cl3d_response_c3'], bins=np.arange(0., 2., 0.02), label='After $p_T$ corr', color='black', histtype='step', lw=2)
  plt.legend(loc = 'upper left', fontsize=22)
  plt.xlabel(r'$p_{T}^{L1}/p_{T}^{gen}$')
  plt.ylabel(r'Entries')
  plt.ylim(0, 1000)
  name = 'calib_responses_'+fe_names[name]
  plt.savefig(plotdir+name+'_'+text+'.png')
  plt.savefig(plotdir+name+'_'+text+'.pdf')

###########

# Model c3

plt.figure(figsize=(15,10))
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

  plt.plot(pt, model_c3[name].predict(np.vstack([logpt1, logpt2, logpt3, logpt4]).T), ls='--', color=colors[name], label=legends[name])
  plt.plot(ptmeans[name]['cl3d_pt_c2'], ptmeans[name]['cl3d_response_c2'], marker='s', markersize=8, ls='None', color=colors[name], label=legends[name])

plt.legend(loc = 'upper right', fontsize=16)
plt.grid()
plt.legend(loc = 'upper right', fontsize=16)
plt.xlabel(r'$\langle p_{T}^{L1}\rangle\, [GeV]$')
plt.ylabel(r'$\langle p_{T}^{L1}/p_{T}^{gen}\rangle$')
plt.savefig(plotdir+'calib_modelc3_'+text+'.png')
plt.savefig(plotdir+'calib_modelc3_'+text+'.pdf')

###########

# Responses and resolutions

plt.figure(figsize=(15,10))

for name in df_tau:

  df = etameans[name]
  plt.plot(df['gentau_vis_abseta'], df['cl3d_response_c3'],color=colors[name], label=legends[name],lw=2)

plt.ylim(0.9, 1.2)
plt.xlabel(r'$|\eta|$')
plt.ylabel(r'$\langle p_{T}^{L1}/p_{T}^{gen}\rangle$')
plt.grid()
plt.legend(loc = 'upper right', fontsize=16)
plt.ylim(0.95, 1.05)
plt.savefig(plotdir+'calib_mean_vs_eta_'+text+'.png')
plt.savefig(plotdir+'calib_mean_vs_eta_'+text+'.pdf')

plt.figure(figsize=(15,10))

for name in df_tau:

  df = etameans[name]
  plt.plot(df['gentau_vis_abseta'], etarmss[name]['cl3d_response_c3'],color=colors[name], label=legends[name],lw=2)

plt.ylim(0., 0.4)
plt.xlabel(r'$|\eta|$')
plt.ylabel(r'$RMS(p_{T}^{L1}/p_{T}^{gen})$')
plt.grid()
plt.legend(loc = 'upper right', fontsize=16)
plt.ylim(0.1, 0.3)
plt.savefig(plotdir+'calib_rms_vs_eta_'+text+'.png')
plt.savefig(plotdir+'calib_rms_vs_eta_'+text+'.pdf')
  
plt.figure(figsize=(15,10))

for name in df_tau:

  df = etameans[name]
  plt.plot(df['gentau_vis_abseta'], etaeffrmss[name]['cl3d_response_c3'],color=colors[name], label=legends[name],lw=2)

plt.ylim(0., 0.4)
plt.xlabel(r'$|\eta|$')
plt.ylabel(r'$RMS_{eff}(p_{T}^{L1}/p_{T}^{gen})$')
plt.grid()
plt.legend(loc = 'upper right', fontsize=16)
plt.ylim(0.0, 0.35)
plt.savefig(plotdir+'calib_effrms_vs_eta_'+text+'.png')
plt.savefig(plotdir+'calib_effrms_vs_eta_'+text+'.pdf')
  
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
plt.savefig(plotdir+'calib_mean_vs_pt_'+text+'.png')
plt.savefig(plotdir+'calib_mean_vs_pt_'+text+'.pdf')

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
plt.savefig(plotdir+'calib_rms_vs_pt_'+text+'.png')
plt.savefig(plotdir+'calib_rms_vs_pt_'+text+'.pdf')
  
plt.figure(figsize=(15,10))

for name in df_tau:

  df = ptmeans[name]
  plt.plot(df['gentau_vis_pt'], pteffrmss[name]['cl3d_response_c3'],color=colors[name], label=legends[name],lw=2)

plt.ylim(0., 0.3)
plt.xlabel(r'$p_{T}\, [GeV]$')
plt.ylabel(r'$RMS_{eff}(p_{T}^{L1}/p_{T}^{gen})$')
plt.grid()
plt.legend(loc = 'upper right', fontsize=16)
plt.ylim(0.0, 0.35)
plt.savefig(plotdir+'calib_effrms_vs_pt_'+text+'.png')
plt.savefig(plotdir+'calib_effrms_vs_pt_'+text+'.pdf')

