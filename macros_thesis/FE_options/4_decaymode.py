# C. Martin Perez cmartinp@cern.ch, Nov. 2019

###########

import os
import pandas as pd
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from sklearn import metrics
from sklearn.metrics import confusion_matrix
from sklearn.ensemble import RandomForestClassifier
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

fe_names[0] = 'Threshold'
fe_names[1] = 'STC'
fe_names[2] = 'BestChoice'
fe_names[3] = 'BestChoiceCoarse'
fe_names[4] = 'MixedBC+STC'

###########

file_in_tau = {}

if algo == 'A':

  dir_in = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoA/calibrated/'

  file_in_tau[0] = dir_in+'ntuple_TauGun_THRESHOLD_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'
  file_in_tau[1] = dir_in+'ntuple_TauGun_SUPERTRIGGERCELL_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'
  file_in_tau[2] = dir_in+'ntuple_TauGun_BESTCHOICE_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'
  file_in_tau[3] = dir_in+'ntuple_TauGun_BESTCHOICECOARSE_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'
  file_in_tau[4] = dir_in+'ntuple_TauGun_MIXEDBCSTC_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'

if algo == 'B':

  dir_in = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoB/pubdt/'

  file_in_tau[0] = dir_in+'ntuple_TauGun_THRESHOLD_calibrated_pubdt_algoB.hdf5'
  file_in_tau[1] = dir_in+'ntuple_TauGun_SUPERTRIGGERCELL_calibrated_pubdt_algoB.hdf5'
  file_in_tau[2] = dir_in+'ntuple_TauGun_BESTCHOICE_calibrated_pubdt_algoB.hdf5'
  file_in_tau[3] = dir_in+'ntuple_TauGun_BESTCHOICECOARSE_calibrated_pubdt_algoB.hdf5'
  file_in_tau[4] = dir_in+'ntuple_TauGun_MIXEDBCSTC_calibrated_pubdt_algoB.hdf5'

###########

file_out_tau = {}

if algo == 'A':

  dir_out = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoA/DMid/'

  file_out_tau[0] = dir_out+'ntuple_TauGun_THRESHOLD_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  file_out_tau[1] = dir_out+'ntuple_TauGun_SUPERTRIGGERCELL_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  file_out_tau[2] = dir_out+'ntuple_TauGun_BESTCHOICE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  file_out_tau[3] = dir_out+'ntuple_TauGun_BESTCHOICECOARSE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  file_out_tau[4] = dir_out+'ntuple_TauGun_MIXEDBCSTC_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'

if algo == 'B':

  dir_out = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoB/DMid/'

  file_out_tau[0] = dir_out+'ntuple_TauGun_THRESHOLD_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  file_out_tau[1] = dir_out+'ntuple_TauGun_SUPERTRIGGERCELL_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  file_out_tau[2] = dir_out+'ntuple_TauGun_BESTCHOICE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  file_out_tau[3] = dir_out+'ntuple_TauGun_BESTCHOICECOARSE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  file_out_tau[4] = dir_out+'ntuple_TauGun_MIXEDBCSTC_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'

###########

file_in_nu = {}

if algo == 'A':

  dir_in = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoA/calibrated/'

  file_in_nu[0] = dir_in+'ntuple_NuGun_THRESHOLD_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'
  file_in_nu[1] = dir_in+'ntuple_NuGun_SUPERTRIGGERCELL_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'
  file_in_nu[2] = dir_in+'ntuple_NuGun_BESTCHOICE_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'
  file_in_nu[3] = dir_in+'ntuple_NuGun_BESTCHOICECOARSE_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'
  file_in_nu[4] = dir_in+'ntuple_NuGun_MIXEDBCSTC_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'

if algo == 'B':

  dir_in = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoB/pubdt/'

  file_in_nu[0] = dir_in+'ntuple_NuGun_THRESHOLD_calibrated_pubdt_algoB.hdf5'
  file_in_nu[1] = dir_in+'ntuple_NuGun_SUPERTRIGGERCELL_calibrated_pubdt_algoB.hdf5'
  file_in_nu[2] = dir_in+'ntuple_NuGun_BESTCHOICE_calibrated_pubdt_algoB.hdf5'
  file_in_nu[3] = dir_in+'ntuple_NuGun_BESTCHOICECOARSE_calibrated_pubdt_algoB.hdf5'
  file_in_nu[4] = dir_in+'ntuple_NuGun_MIXEDBCSTC_calibrated_pubdt_algoB.hdf5'

###########

file_out_nu = {}

if algo == 'A':

  dir_out = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoA/DMid/'

  file_out_nu[0] = dir_out+'ntuple_NuGun_THRESHOLD_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  file_out_nu[1] = dir_out+'ntuple_NuGun_SUPERTRIGGERCELL_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  file_out_nu[2] = dir_out+'ntuple_NuGun_BESTCHOICE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  file_out_nu[3] = dir_out+'ntuple_NuGun_BESTCHOICECOARSE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  file_out_nu[4] = dir_out+'ntuple_NuGun_MIXEDBCSTC_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'

if algo == 'B':

  dir_out = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoB/DMid/'

  file_out_nu[0] = dir_out+'ntuple_NuGun_THRESHOLD_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  file_out_nu[1] = dir_out+'ntuple_NuGun_SUPERTRIGGERCELL_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  file_out_nu[2] = dir_out+'ntuple_NuGun_BESTCHOICE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  file_out_nu[3] = dir_out+'ntuple_NuGun_BESTCHOICECOARSE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  file_out_nu[4] = dir_out+'ntuple_NuGun_MIXEDBCSTC_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'

###########

file_out_model = {}

if algo == 'A':

  dir_out_model = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoA/models/DMid/'

  file_out_model[0] = dir_out_model+'model_DMid_algoA_WP'+PUBDTWP+'_threshold.pkl'
  file_out_model[1] = dir_out_model+'model_DMid_algoA_WP'+PUBDTWP+'_supertriggercell.pkl'
  file_out_model[2] = dir_out_model+'model_DMid_algoA_WP'+PUBDTWP+'_bestchoice.pkl'
  file_out_model[3] = dir_out_model+'model_DMid_algoA_WP'+PUBDTWP+'_bestchoicecoarse.pkl'
  file_out_model[4] = dir_out_model+'model_DMid_algoA_WP'+PUBDTWP+'_mixedbcstc.pkl'

if algo == 'B':

  dir_out_model = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoB/models/DMid/'

  file_out_model[0] = dir_out_model+'model_DMid_algoB_WP'+PUBDTWP+'_threshold.pkl'
  file_out_model[1] = dir_out_model+'model_DMid_algoB_WP'+PUBDTWP+'_supertriggercell.pkl'
  file_out_model[2] = dir_out_model+'model_DMid_algoB_WP'+PUBDTWP+'_bestchoice.pkl'
  file_out_model[3] = dir_out_model+'model_DMid_algoB_WP'+PUBDTWP+'_bestchoicecoarse.pkl'
  file_out_model[4] = dir_out_model+'model_DMid_algoB_WP'+PUBDTWP+'_mixedbcstc.pkl'


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

df_nu = {}

for name in file_in_nu:

  store_nu = pd.HDFStore(file_in_nu[name], mode='r')
  df_nu[name]  = store_nu['df_nu_PU200']
  store_nu.close()

###########

df_tau_train = {}

for name in df_tau:

  df_tau[name]['gentau_decayMode'].replace([4,5],2, inplace=True)

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

  sel = df_tau_train[name]['gentau_decayMode'] >= 0
  df_tau_train[name] = df_tau_train[name][sel]


###########

print '-- Decay mode ID --'

# Training

features = [ 'cl3d_pt_c3',
    'cl3d_abseta', 'cl3d_showerlength',
    'cl3d_coreshowerlength', 'cl3d_firstlayer', 'cl3d_maxlayer', 'cl3d_szz',
    'cl3d_seetot', 'cl3d_spptot', 'cl3d_srrtot', 'cl3d_srrmean',
    'cl3d_hoe', 'cl3d_meanz', 'cl3d_layer10',
    'cl3d_layer50', 'cl3d_layer90', 'cl3d_ntc67', 'cl3d_ntc90'
  ]

model = {}

for name in df_tau_train:

  inputs = df_tau_train[name][features]
  target = df_tau_train[name]['gentau_decayMode']
  model[name] = RandomForestClassifier(n_jobs=10,
                             random_state=0,
                             class_weight='balanced',
                             max_depth=2,
                             n_estimators=1000).fit(inputs,target)

for name in df_tau:

  with open(file_out_model[name], 'wb') as f:
    pickle.dump(model[name], f)

# Application

# Tau

for name in df_tau:

  df_tau[name]['cl3d_predDM'] = model[name].predict(df_tau[name][features])

for name in df_tau:

  probas = model[name].predict_proba(df_tau[name][features])

  df_tau[name]['cl3d_probDM0'] = probas[:,0]
  df_tau[name]['cl3d_probDM1'] = probas[:,1]
  df_tau[name]['cl3d_probDM2'] = probas[:,2]

decaymodes = {}

for name in df_tau:

  decaymodes[name] = df_tau[name].groupby('gentau_decayMode')['gentau_vis_pt'].count()

# Nu

for name in df_nu:

  df_nu[name]['cl3d_predDM'] = model[name].predict(df_nu[name][features])

for name in df_nu:

  probas = model[name].predict_proba(df_nu[name][features])

  df_nu[name]['cl3d_probDM0'] = probas[:,0]
  df_nu[name]['cl3d_probDM1'] = probas[:,1]
  df_nu[name]['cl3d_probDM2'] = probas[:,2]

###########

# SAVE DFS

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

# Confusion matrix

cm = {}
ns = {}

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

  sel = df_tau[name]['gentau_decayMode'] >= 0
  df_tau[name] = df_tau[name][sel]

for name in df_tau:
  cm[name] = confusion_matrix(df_tau[name].gentau_decayMode, df_tau[name].cl3d_predDM)
  cm[name] = cm[name].astype('float') / cm[name].sum(axis=1)[:, np.newaxis]
  print ' ',fe_names[name]
  print cm[name]
  ns[name] = decaymodes[name]
  print ((cm[name][0,0]*ns[name][0]+cm[name][1,1]*ns[name][1]+cm[name][2,2]*ns[name][2])/np.sum(ns[name]))
  print ' '

# PLOTTING
'''
if algo == 'A':

  plotdir = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoA/plots/DMid/'

if algo == 'B':

  plotdir = '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Nov19/algoB/plots/DMid/' 

fe_names = {}

fe_names[0] = 'Threshold'
fe_names[1] = 'STC'
fe_names[2] = 'BestChoice'
fe_names[3] = 'BestChoiceCoarse'
fe_names[4] = 'MixedBC+STC'

matplotlib.rcParams.update({'font.size': 16})

for name in df_tau:

  plt.figure(figsize=(15,10))
  fig, ax = plt.subplots()
  im = ax.imshow(cm[name], interpolation='nearest', cmap=plt.cm.Blues)
  ax.figure.colorbar(im, ax=ax)
  plt.xticks([0,1,2], ('1-prong', r'1-prong$+\pi_0$', '3-prongs($+\pi_0$)'))
  plt.yticks([0,1,2], ('1-prong', r'1-prong$+\pi_0$', '3-prongs($+\pi_0$)'))
  plt.xlim(-0.5, 2.5)
  plt.ylim(-0.5, 2.5),
  plt.ylabel('True label')
  plt.xlabel('Predicted label')
  plt.setp(ax.get_xticklabels(), rotation=30, ha="right",
         rotation_mode="anchor")
  fmt = '.2f'
  thresh = cm[name].max() / 2.
  for i in range(cm[name].shape[0]):
      for j in range(cm[name].shape[1]):
          ax.text(j, i, format(cm[name][i, j], fmt),
                  ha="center", va="center",
                  color="white" if cm[name][i, j] > thresh else "black")
  plt.tight_layout()
  plt.savefig(plotdir+'dm_cm_'+fe_names[name]+'_algo'+algo+'_WP'+PUBDTWP+'.png')
  plt.savefig(plotdir+'dm_cm_'+fe_names[name]+'_algo'+algo+'_WP'+PUBDTWP+'.pdf')

# Predicted probabilities

for name in df_tau:

  plt.figure(figsize=(15,10))
  plt.hist(df_tau[name][df_tau[name].gentau_decayMode==0].cl3d_probDM0,
           bins=np.arange(0., 1., 0.05), histtype='step', lw=2, label='1-prong', normed=True)
  plt.hist(df_tau[name][df_tau[name].gentau_decayMode==1].cl3d_probDM0,
           bins=np.arange(0., 1., 0.05), histtype='step', lw=2, label=r'1-prong$+\pi_0$', normed=True)
  plt.hist(df_tau[name][df_tau[name].gentau_decayMode==2].cl3d_probDM0,
           bins=np.arange(0., 1., 0.05), histtype='step', lw=2, label='3-prongs($+\pi_0$)', normed=True)
  plt.legend(fontsize=22)
  plt.xlabel(r'Probability 1-prong')
  plt.ylabel(r'Entries')
  plt.savefig(plotdir+'dm_probDM0'+fe_names[name]+'_algo'+algo+'_WP'+PUBDTWP+'.png')
  plt.savefig(plotdir+'dm_probDM0'+fe_names[name]+'_algo'+algo+'_WP'+PUBDTWP+'.pdf')
  
  plt.figure(figsize=(15,10))
  plt.hist(df_tau[name][df_tau[name].gentau_decayMode==0].cl3d_probDM1,
           bins=np.arange(0., 1., 0.05), histtype='step', lw=2, label='1-prong', normed=True)
  plt.hist(df_tau[name][df_tau[name].gentau_decayMode==1].cl3d_probDM1,
           bins=np.arange(0., 1., 0.05), histtype='step', lw=2, label=r'1-prong$+\pi_0$', normed=True)
  plt.hist(df_tau[name][df_tau[name].gentau_decayMode==2].cl3d_probDM1,
           bins=np.arange(0., 1., 0.05), histtype='step', lw=2, label='3-prongs($+\pi_0$)', normed=True)
  plt.legend(fontsize=22)
  plt.xlabel(r'Probability 1-prong$+\pi_0$')
  plt.ylabel(r'Entries')
  plt.savefig(plotdir+'dm_probDM1'+fe_names[name]+'_algo'+algo+'_WP'+PUBDTWP+'.png')
  plt.savefig(plotdir+'dm_probDM1'+fe_names[name]+'_algo'+algo+'_WP'+PUBDTWP+'.pdf')
  
  plt.figure(figsize=(15,10))
  plt.hist(df_tau[name][df_tau[name].gentau_decayMode==0].cl3d_probDM2,
           bins=np.arange(0., 1., 0.05), histtype='step', lw=2, label='1-prong', normed=True)
  plt.hist(df_tau[name][df_tau[name].gentau_decayMode==1].cl3d_probDM2,
           bins=np.arange(0., 1., 0.05), histtype='step', lw=2, label=r'1-prong$+\pi_0$', normed=True)
  plt.hist(df_tau[name][df_tau[name].gentau_decayMode==2].cl3d_probDM2,
           bins=np.arange(0., 1., 0.05), histtype='step', lw=2, label='3-prongs($+\pi_0$)', normed=True)
  plt.legend(fontsize=22)
  plt.xlabel(r'Probability 3-prongs($+\pi_0$)')
  plt.ylabel(r'Entries')
  plt.savefig(plotdir+'dm_probDM2'+fe_names[name]+'_algo'+algo+'_WP'+PUBDTWP+'.png')
  plt.savefig(plotdir+'dm_probDM2'+fe_names[name]+'_algo'+algo+'_WP'+PUBDTWP+'.pdf')
'''
