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
#fe_names[1] = 'STC'
#fe_names[2] = 'BestChoice'
#fe_names[3] = 'BestChoiceCoarse'
#fe_names[4] = 'MixedBC+STC'

###########

file_in_tau = {}

if algo == 'A':

  dir_in = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoA/calibrated/'

  file_in_tau[0] = dir_in+'ntuple_TauGun_THRESHOLD_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'
  #file_in_tau[1] = dir_in+'ntuple_TauGun_SUPERTRIGGERCELL_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'
  #file_in_tau[2] = dir_in+'ntuple_TauGun_BESTCHOICE_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'
  #file_in_tau[3] = dir_in+'ntuple_TauGun_BESTCHOICECOARSE_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'
  #file_in_tau[4] = dir_in+'ntuple_TauGun_MIXEDBCSTC_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'

if algo == 'B':

  dir_in = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoB/pubdt/'

  file_in_tau[0] = dir_in+'ntuple_TauGun_THRESHOLD_calibrated_pubdt_algoB.hdf5'
  #file_in_tau[1] = dir_in+'ntuple_TauGun_SUPERTRIGGERCELL_calibrated_pubdt_algoB.hdf5'
  #file_in_tau[2] = dir_in+'ntuple_TauGun_BESTCHOICE_calibrated_pubdt_algoB.hdf5'
  #file_in_tau[3] = dir_in+'ntuple_TauGun_BESTCHOICECOARSE_calibrated_pubdt_algoB.hdf5'
  #file_in_tau[4] = dir_in+'ntuple_TauGun_MIXEDBCSTC_calibrated_pubdt_algoB.hdf5'

###########

file_out_tau = {}

if algo == 'A':

  dir_out = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoA/DMid/'

  file_out_tau[0] = dir_out+'ntuple_TauGun_THRESHOLD_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  #file_out_tau[1] = dir_out+'ntuple_TauGun_SUPERTRIGGERCELL_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  #file_out_tau[2] = dir_out+'ntuple_TauGun_BESTCHOICE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  #file_out_tau[3] = dir_out+'ntuple_TauGun_BESTCHOICECOARSE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  #file_out_tau[4] = dir_out+'ntuple_TauGun_MIXEDBCSTC_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'

if algo == 'B':

  dir_out = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoB/DMid/'

  file_out_tau[0] = dir_out+'ntuple_TauGun_THRESHOLD_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  #file_out_tau[1] = dir_out+'ntuple_TauGun_SUPERTRIGGERCELL_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  #file_out_tau[2] = dir_out+'ntuple_TauGun_BESTCHOICE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  #file_out_tau[3] = dir_out+'ntuple_TauGun_BESTCHOICECOARSE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  #file_out_tau[4] = dir_out+'ntuple_TauGun_MIXEDBCSTC_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'

###########

file_in_nu = {}

if algo == 'A':

  dir_in = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoA/calibrated/'

  file_in_nu[0] = dir_in+'ntuple_NuGun_THRESHOLD_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'
  #file_in_nu[1] = dir_in+'ntuple_NuGun_SUPERTRIGGERCELL_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'
  #file_in_nu[2] = dir_in+'ntuple_NuGun_BESTCHOICE_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'
  #file_in_nu[3] = dir_in+'ntuple_NuGun_BESTCHOICECOARSE_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'
  #file_in_nu[4] = dir_in+'ntuple_NuGun_MIXEDBCSTC_pubdt_calibratedWP'+PUBDTWP+'_algoA.hdf5'

if algo == 'B':

  dir_in = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoB/pubdt/'

  file_in_nu[0] = dir_in+'ntuple_NuGun_THRESHOLD_calibrated_pubdt_algoB.hdf5'
  #file_in_nu[1] = dir_in+'ntuple_NuGun_SUPERTRIGGERCELL_calibrated_pubdt_algoB.hdf5'
  #file_in_nu[2] = dir_in+'ntuple_NuGun_BESTCHOICE_calibrated_pubdt_algoB.hdf5'
  #file_in_nu[3] = dir_in+'ntuple_NuGun_BESTCHOICECOARSE_calibrated_pubdt_algoB.hdf5'
  #file_in_nu[4] = dir_in+'ntuple_NuGun_MIXEDBCSTC_calibrated_pubdt_algoB.hdf5'

###########

file_out_nu = {}

if algo == 'A':

  dir_out = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoA/DMid/'

  file_out_nu[0] = dir_out+'ntuple_NuGun_THRESHOLD_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  #file_out_nu[1] = dir_out+'ntuple_NuGun_SUPERTRIGGERCELL_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  #file_out_nu[2] = dir_out+'ntuple_NuGun_BESTCHOICE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  #file_out_nu[3] = dir_out+'ntuple_NuGun_BESTCHOICECOARSE_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'
  #file_out_nu[4] = dir_out+'ntuple_NuGun_MIXEDBCSTC_pubdt_calibratedWP'+PUBDTWP+'_DMid_algoA.hdf5'

if algo == 'B':

  dir_out = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoB/DMid/'

  file_out_nu[0] = dir_out+'ntuple_NuGun_THRESHOLD_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  #file_out_nu[1] = dir_out+'ntuple_NuGun_SUPERTRIGGERCELL_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  #file_out_nu[2] = dir_out+'ntuple_NuGun_BESTCHOICE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  #file_out_nu[3] = dir_out+'ntuple_NuGun_BESTCHOICECOARSE_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'
  #file_out_nu[4] = dir_out+'ntuple_NuGun_MIXEDBCSTC_calibrated_pubdt'+PUBDTWP+'_DMid_algoB.hdf5'

###########

file_out_model = {}

if algo == 'A':

  dir_out_model = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoA/models/DMid/'

  file_out_model[0] = dir_out_model+'model_DMid_algoA_WP'+PUBDTWP+'_threshold.pkl'
  #file_out_model[1] = dir_out_model+'model_DMid_algoA_WP'+PUBDTWP+'_supertriggercell.pkl'
  #file_out_model[2] = dir_out_model+'model_DMid_algoA_WP'+PUBDTWP+'_bestchoice.pkl'
  #file_out_model[3] = dir_out_model+'model_DMid_algoA_WP'+PUBDTWP+'_bestchoicecoarse.pkl'
  #file_out_model[4] = dir_out_model+'model_DMid_algoA_WP'+PUBDTWP+'_mixedbcstc.pkl'

if algo == 'B':

  dir_out_model = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoB/models/DMid/'

  file_out_model[0] = dir_out_model+'model_DMid_algoB_WP'+PUBDTWP+'_threshold.pkl'
  #file_out_model[1] = dir_out_model+'model_DMid_algoB_WP'+PUBDTWP+'_supertriggercell.pkl'
  #file_out_model[2] = dir_out_model+'model_DMid_algoB_WP'+PUBDTWP+'_bestchoice.pkl'
  #file_out_model[3] = dir_out_model+'model_DMid_algoB_WP'+PUBDTWP+'_bestchoicecoarse.pkl'
  #file_out_model[4] = dir_out_model+'model_DMid_algoB_WP'+PUBDTWP+'_mixedbcstc.pkl'


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

# PLOTTING

import geeksw.plotting.cmsplot as plt
from geeksw.plotting.root_colors import *

plt.matplotlib.font_manager._rebuild()

import matplotlib.lines as mlines

df_tau_DM0 = {}
df_tau_DM1 = {}
df_tau_DM4 = {}
df_tau_DM5 = {}
df_tau_DM45 = {}

plotdir = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoB/plots/DMid/' 

for name in df_tau:

  seldm0 = (df_tau_train[name]['gentau_decayMode'] == 0)
  df_tau_DM0[name] = df_tau_train[name][seldm0]

  seldm1 = (df_tau_train[name]['gentau_decayMode'] == 1)
  df_tau_DM1[name] = df_tau_train[name][seldm1]

  seldm4 = (df_tau_train[name]['gentau_decayMode'] == 4)
  df_tau_DM4[name] = df_tau_train[name][seldm4]

  seldm5 = (df_tau_train[name]['gentau_decayMode'] == 5)
  df_tau_DM5[name] = df_tau_train[name][seldm5]

  seldm45 = (df_tau_train[name]['gentau_decayMode'] > 3)
  df_tau_DM45[name] = df_tau_train[name][seldm45]
  print df_tau_DM45[name]['gentau_decayMode']

###########

for name in df_tau:

  df_tau[name]['gentau_decayMode'].replace([4,5],2, inplace=True)

myvariable = 'cl3d_hoe'
myxmin = 0.
myxmax = 2.5
mybinsize = 0.1

for name in df_tau:

  plt.figure(figsize=(8,8))

  plt.hist(df_tau_DM0[name][myvariable], normed=True, bins=np.arange(myxmin,myxmax,mybinsize),  label='1-prong', color='limegreen',  histtype='step', lw=2)
  plt.hist(df_tau_DM1[name][myvariable], normed=True, bins=np.arange(myxmin,myxmax,mybinsize),  label='1-prong + $\pi^{0}$\'s', color='blue',  histtype='step', lw=2)
  plt.hist(df_tau_DM45[name][myvariable], normed=True, bins=np.arange(myxmin,myxmax,mybinsize),  label='3-prongs (+ $\pi^{0}$\'s)', color='red',  histtype='step', lw=2)
  green_line = mlines.Line2D([], [], color='limegreen',markersize=15, label='1-prong',lw=2)
  blue_line = mlines.Line2D([], [], color='blue',markersize=15, label='1-prong + $\pi^{0}$\'s',lw=2)
  red_line = mlines.Line2D([], [], color='red',markersize=15, label='3-prongs (+ $\pi^{0}$\'s)',lw=2)
  plt.legend(loc = 'upper right',handles=[green_line,blue_line,red_line])
  plt.grid()
  plt.xlabel('Energy in CE-H / Energy in CE-E')
  plt.ylabel(r'Normalized events')
  #plt.ylim(0.0, 0.07)
  plt.cmstext("CMS"," Phase-2 Simulation")
  plt.lumitext("PU=200"," ")
  plt.subplots_adjust(bottom=0.12)
  plt.savefig(plotdir+'DMid_'+myvariable+'_DM_merged.png')
  plt.savefig(plotdir+'DMid_'+myvariable+'_DM_merged.pdf')
  plt.show()

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


'''for name in df_tau:

  feature_importance = model[name].feature_importances_
  sorted_idx = np.argsort(feature_importance)
  pos = np.arange(sorted_idx.shape[0]) + .5
  fig = plt.figure(figsize=(12, 6))
  plt.subplot(1, 2, 1)
  plt.barh(pos, feature_importance[sorted_idx], align='center')
  plt.yticks(pos, np.array(features)[sorted_idx])
  plt.title('Feature Importance (MDI)')
  plt.show()
'''

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

if algo == 'A':

  plotdir = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoA/plots/DMid/'

if algo == 'B':

  plotdir = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoB/plots/DMid/' 

fe_names = {}

fe_names[0] = 'Threshold'
#fe_names[1] = 'STC'
#fe_names[2] = 'BestChoice'
#fe_names[3] = 'BestChoiceCoarse'
#fe_names[4] = 'MixedBC+STC'

'''
for name in df_tau:

  plt.figure(figsize=(8,8))
  fig, ax = plt.subplots()
  im = ax.imshow(cm[name], interpolation='nearest', cmap=plt.cm.Blues)
  ax.figure.colorbar(im, ax=ax)
  plt.xticks([0,1,2], ('1-prong', r'1-prong$+\pi_0$', '3-prongs($+\pi_0$)'))
  plt.yticks([0,1,2], ('1-prong', r'1-prong$+\pi_0$', '3-prongs($+\pi_0$)'))
  plt.xlim(-0.5, 2.5)
  plt.ylim(-0.5, 2.5),
  plt.ylabel('Generated decay mode')
  plt.xlabel('Predicted decay mode')
  plt.setp(ax.get_xticklabels(), rotation=30, ha="right",
         rotation_mode="anchor")
  plt.setp(ax.get_yticklabels(), rotation=30, ha="right",
         rotation_mode="anchor")
  fmt = '.2f'
  thresh = cm[name].max() / 2.
  for i in range(cm[name].shape[0]):
      for j in range(cm[name].shape[1]):
          ax.text(j, i, format(cm[name][i, j], fmt),
                  ha="center", va="center", fontsize=18,
                  color="white" if cm[name][i, j] > thresh else "black")
  plt.tight_layout()
  plt.cmstext("CMS","    Phase-2 Simulation")
  plt.lumitext("PU=200"," ")
  plt.subplots_adjust(bottom=0.21)
  plt.subplots_adjust(top=0.9)
  #plt.subplots_adjust(left=0.01)
  plt.savefig(plotdir+'dm_cm_'+fe_names[name]+'_algo'+algo+'_WP'+PUBDTWP+'.png')
  plt.savefig(plotdir+'dm_cm_'+fe_names[name]+'_algo'+algo+'_WP'+PUBDTWP+'.pdf')
  plt.show()
'''



# Predicted probabilities
'''
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
