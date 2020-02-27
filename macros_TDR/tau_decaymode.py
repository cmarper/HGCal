# C. Martin Perez cmartinp@cern.ch, Sep. 2019

###########

import os
from glob import glob
import itertools
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn import metrics
import xgboost as xgb
import matplotlib
from scipy.optimize import lsq_linear
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LinearRegression
import pickle
from sklearn.metrics import confusion_matrix

###########

dir_in = '/data_CMS/cms/mperez/HGCal_data/Sep19/calibrated/'
filein_tau_PU200 = dir_in+'RelValDiTau_Pt20To100_v10_PU200_calibrated.hdf5'

store_tau_PU200 = pd.HDFStore(filein_tau_PU200, mode='r')
df_tau_PU200	= store_tau_PU200['df_tau_PU200']
store_tau_PU200.close()

dir_out = '/data_CMS/cms/mperez/HGCal_data/Sep19/DMid/'
fileout_tau_PU200 = dir_out+'RelValDiTau_Pt20To100_v10_PU200_DMid.hdf5'

###########

ptcut 	= 20
etamin 	= 1.6
etamax 	= 2.9

sel = df_tau_PU200['gentau_vis_pt'] > ptcut
df_tau_PU200 = df_tau_PU200[sel]

sel = np.abs(df_tau_PU200['gentau_vis_eta']) > etamin
df_tau_PU200 = df_tau_PU200[sel]

sel = np.abs(df_tau_PU200['gentau_vis_eta']) < etamax
df_tau_PU200 = df_tau_PU200[sel]

sel = df_tau_PU200['cl3d_isbestmatch'] == True
df_tau_PU200 = df_tau_PU200[sel]

sel = df_tau_PU200['cl3d_pubdt_passWP99'] == True
df_tau_PU200 = df_tau_PU200[sel]

sel = df_tau_PU200['gentau_decayMode'] >= 0
df_tau_PU200 = df_tau_PU200[sel]

df_tau_PU200['gentau_decayMode'].replace([4,5],2, inplace=True)


###########

print '-- Decay Mode id --'

# Training

features = [ 'cl3d_pt_c3',
    'cl3d_abseta', 'cl3d_showerlength',
    'cl3d_coreshowerlength', 'cl3d_firstlayer', 'cl3d_maxlayer', 'cl3d_szz',
    'cl3d_seetot', 'cl3d_spptot', 'cl3d_srrtot', 'cl3d_srrmean',
    'cl3d_hoe', 'cl3d_meanz', 'cl3d_layer10',
    'cl3d_layer50', 'cl3d_layer90', 'cl3d_ntc67', 'cl3d_ntc90'
    ]

inputs = df_tau_PU200[features]
target = df_tau_PU200['gentau_decayMode']

model = RandomForestClassifier(n_jobs=10,
                             random_state=0,
                             class_weight='balanced',
                             max_depth=2,
                             n_estimators=1000).fit(inputs,target)

with open('/data_CMS/cms/mperez/HGCal_data/Sep19/models/model_DMid.pkl', 'wb') as f:
  pickle.dump(model, f)

# Predicted DM
df_tau_PU200['cl3d_predDM'] = model.predict(df_tau_PU200[features])

# Predicted probabilities
probas = model.predict_proba(df_tau_PU200[features])
df_tau_PU200['cl3d_probDM0'] = probas[:,0]
df_tau_PU200['cl3d_probDM1'] = probas[:,1]
df_tau_PU200['cl3d_probDM2'] = probas[:,2]

decaymodes = df_tau_PU200.groupby('gentau_decayMode')['gentau_vis_pt'].count()

###########

# SAVE DFS

# Save files
store_tau_PU200 = pd.HDFStore(fileout_tau_PU200, mode='w')
store_tau_PU200['df_tau_PU200'] = df_tau_PU200
store_tau_PU200.close()

###########

# PLOTTING

plotdir = '/home/llr/cms/mperez/HGCal/v10_geometry/tau_algorithm/plots/'
import geeksw.plotting.cmsplot as plt
from geeksw.plotting.root_colors import *

plt.matplotlib.font_manager._rebuild()


# Confusion matrix

cm = confusion_matrix(df_tau_PU200.gentau_decayMode, df_tau_PU200.cl3d_predDM)
cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
print(cm)
ns = decaymodes
print(ns)
print((cm[0,0]*ns[0]+cm[1,1]*ns[1]+cm[2,2]*ns[2])/np.sum(ns))

'''matplotlib.rcParams.update({'font.size': 17})
plt.figure(figsize=(15,10))
fig, ax = plt.subplots()
im = ax.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
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
thresh = cm.max() / 2.
for i in range(cm.shape[0]):
    for j in range(cm.shape[1]):
        ax.text(j, i, format(cm[i, j], fmt),
                ha="center", va="center",
                color="white" if cm[i, j] > thresh else "black")
plt.tight_layout()
plt.savefig(plotdir+'DMid_confusionmatrix.png')
plt.savefig(plotdir+'DMid_confusionmatrix.pdf')
'''

plt.figure(figsize=(8,8))
fig, ax = plt.subplots()
im = ax.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
ax.figure.colorbar(im, ax=ax)
plt.xticks([0,1,2], ('1-prong', r'1-prong$+\pi_0$', '3-prongs($+\pi_0$)'))
plt.yticks([0,1,2], ('1-prong', r'1-prong$+\pi_0$', '3-prongs($+\pi_0$)'))
plt.xlim(-0.5, 2.5)
plt.ylim(-0.5, 2.5),
plt.ylabel('True decay mode')
plt.xlabel('Predicted decay mode')
plt.setp(ax.get_xticklabels(), rotation=30, ha="right",
         rotation_mode="anchor")
plt.setp(ax.get_yticklabels(), rotation=30, ha="right",
         rotation_mode="anchor")
fmt = '.2f'
thresh = cm.max() / 2.
for i in range(cm.shape[0]):
    for j in range(cm.shape[1]):
        ax.text(j, i, format(cm[i, j], fmt),
                ha="center", va="center", fontsize=18,
                color="white" if cm[i, j] > thresh else "black")
plt.tight_layout()
plt.cmstext("CMS","    Phase-2 Simulation")
plt.lumitext("200 PU","HGCal calorimeter-only")
plt.subplots_adjust(bottom=0.21)
plt.subplots_adjust(top=0.9)
#plt.subplots_adjust(left=0.01)
plt.savefig(plotdir+'DMid_confusionmatrix_TDR2.png')
plt.savefig(plotdir+'DMid_confusionmatrix_TDR2.pdf')
plt.show()

# Predicted probabilities
'''
plt.figure(figsize=(15,10))
plt.hist(df_tau_PU200[df_tau_PU200.gentau_decayMode==0].cl3d_probDM0,
         bins=np.arange(0., 1., 0.05), histtype='step', lw=2, label='1-prong', normed=True)
plt.hist(df_tau_PU200[df_tau_PU200.gentau_decayMode==1].cl3d_probDM0,
         bins=np.arange(0., 1., 0.05), histtype='step', lw=2, label=r'1-prong$+\pi_0$', normed=True)
plt.hist(df_tau_PU200[df_tau_PU200.gentau_decayMode==2].cl3d_probDM0,
         bins=np.arange(0., 1., 0.05), histtype='step', lw=2, label='3-prongs($+\pi_0$)', normed=True)
plt.legend(fontsize=22)
plt.xlabel(r'Probability 1-prong')
plt.ylabel(r'Entries')
plt.savefig(plotdir+'DMid_probDM0.png')
plt.savefig(plotdir+'DMid_probDM0.pdf')

plt.figure(figsize=(15,10))
plt.hist(df_tau_PU200[df_tau_PU200.gentau_decayMode==0].cl3d_probDM1,
         bins=np.arange(0., 1., 0.05), histtype='step', lw=2, label='1-prong', normed=True)
plt.hist(df_tau_PU200[df_tau_PU200.gentau_decayMode==1].cl3d_probDM1,
         bins=np.arange(0., 1., 0.05), histtype='step', lw=2, label=r'1-prong$+\pi_0$', normed=True)
plt.hist(df_tau_PU200[df_tau_PU200.gentau_decayMode==2].cl3d_probDM1,
         bins=np.arange(0., 1., 0.05), histtype='step', lw=2, label='3-prongs($+\pi_0$)', normed=True)
plt.legend(fontsize=22)
plt.xlabel(r'Probability 1-prong$+\pi_0$')
plt.ylabel(r'Entries')
plt.savefig(plotdir+'DMid_probDM1.png')
plt.savefig(plotdir+'DMid_probDM1.pdf')

plt.figure(figsize=(15,10))
plt.hist(df_tau_PU200[df_tau_PU200.gentau_decayMode==0].cl3d_probDM2,
         bins=np.arange(0., 1., 0.05), histtype='step', lw=2, label='1-prong', normed=True)
plt.hist(df_tau_PU200[df_tau_PU200.gentau_decayMode==1].cl3d_probDM2,
         bins=np.arange(0., 1., 0.05), histtype='step', lw=2, label=r'1-prong$+\pi_0$', normed=True)
plt.hist(df_tau_PU200[df_tau_PU200.gentau_decayMode==2].cl3d_probDM2,
         bins=np.arange(0., 1., 0.05), histtype='step', lw=2, label='3-prongs($+\pi_0$)', normed=True)
plt.legend(fontsize=22)
plt.xlabel(r'Probability 3-prongs($+\pi_0$)')
plt.ylabel(r'Entries')
plt.savefig(plotdir+'DMid_probDM2.png')
plt.savefig(plotdir+'DMid_probDM2.pdf')
'''
