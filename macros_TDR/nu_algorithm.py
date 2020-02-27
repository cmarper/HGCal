# C. Martin Perez cmartinp@cern.ch, Sep. 2019

###########

import os
import pandas as pd
import numpy as np
import root_pandas
import pickle
import xgboost as xgb

###########
###########

dir_in 	= '/data_CMS_upgrade/mperez/HGCal_data/Sep19/skimmed/'
dir_out = '/data_CMS_upgrade/mperez/HGCal_data/Sep19/rates/'

filein_nu_PU200 = dir_in+'ntuple_Nu_E10_v10_PU200_skimmed.root'

fileout_rates_inclusive = dir_out+'ntuple_Nu_E10_v10_PU200_leading.root'
fileout_rates_DM0 = dir_out+'ntuple_Nu_E10_v10_PU200_rates_leading_DM0.root'
fileout_rates_DM1 = dir_out+'ntuple_Nu_E10_v10_PU200_rates_leading_DM1.root'
fileout_rates_DM2 = dir_out+'ntuple_Nu_E10_v10_PU200_rates_leading_DM2.root'

fileout_rates_inclusive_sublead = dir_out+'ntuple_Nu_E10_v10_PU200_subleading.root'

treename = 'SkimmedTree'

branches_event_cl3d 	= ['event','cl3d_pt','cl3d_eta','cl3d_phi','cl3d_showerlength','cl3d_coreshowerlength',
	'cl3d_firstlayer','cl3d_maxlayer','cl3d_seetot','cl3d_spptot','cl3d_szz', 'cl3d_srrtot', 'cl3d_srrmean',
	'cl3d_hoe', 'cl3d_meanz', 'cl3d_layer10', 'cl3d_layer50', 'cl3d_layer90', 'cl3d_ntc67', 'cl3d_ntc90']
branches_cl3d 		= ['cl3d_pt','cl3d_eta','cl3d_phi','cl3d_showerlength','cl3d_coreshowerlength',
	'cl3d_firstlayer','cl3d_maxlayer','cl3d_seetot','cl3d_spptot','cl3d_szz', 'cl3d_srrtot', 'cl3d_srrmean',
	'cl3d_hoe', 'cl3d_meanz', 'cl3d_layer10', 'cl3d_layer50', 'cl3d_layer90', 'cl3d_ntc67', 'cl3d_ntc90']

df_nu_PU200_cl3d 	= root_pandas.read_root( filein_nu_PU200, 	key=treename, columns=branches_event_cl3d, flatten=branches_cl3d )

###########
###########

# PU BDT

print 'Applying PU BDT...'

df_nu_PU200_cl3d['cl3d_abseta'] = np.abs(df_nu_PU200_cl3d.cl3d_eta)

inputs = ['cl3d_abseta','cl3d_showerlength','cl3d_coreshowerlength',
   'cl3d_firstlayer', 'cl3d_maxlayer',
   'cl3d_szz', 'cl3d_seetot', 'cl3d_spptot', 'cl3d_srrtot', 'cl3d_srrmean',
   'cl3d_hoe', 'cl3d_meanz',
   'cl3d_layer10', 'cl3d_layer50', 'cl3d_layer90',
   'cl3d_ntc67', 'cl3d_ntc90']

filename = '/data_CMS_upgrade/mperez/HGCal_data/Sep19/models/model_pubdt.pkl'
model = pickle.load(open(filename, 'rb'))

full = xgb.DMatrix(data=df_nu_PU200_cl3d[inputs], feature_names=inputs) 
df_nu_PU200_cl3d['cl3d_pubdt_score'] = model.predict(full)

bdt_thr_WP99 = 0.108599005491
bdt_thr_WP95 = 0.44313165769
bdt_thr_WP90 = 0.642932797968

df_nu_PU200_cl3d['cl3d_pubdt_passWP99'] = df_nu_PU200_cl3d['cl3d_pubdt_score'] > bdt_thr_WP99
df_nu_PU200_cl3d['cl3d_pubdt_passWP99'] = df_nu_PU200_cl3d['cl3d_pubdt_score'] > bdt_thr_WP95
df_nu_PU200_cl3d['cl3d_pubdt_passWP90'] = df_nu_PU200_cl3d['cl3d_pubdt_score'] > bdt_thr_WP90

sel = df_nu_PU200_cl3d['cl3d_pubdt_passWP99'] == True
df_nu_PU200_cl3d = df_nu_PU200_cl3d[sel]

print '...Done'
print ' '

###########
###########

# CALIBRATION

print 'Applying calibration...'

df_nu_PU200_cl3d['n_matched_cl3d'] = 1

print('Raw: mean={0}, rms={1}'.format(
    df_nu_PU200_cl3d['cl3d_pt'].mean(),
    df_nu_PU200_cl3d['cl3d_pt'].std(),
))

###########

filename_c1 = '/data_CMS_upgrade/mperez/HGCal_data/Sep19/models/model_calib1.pkl'
model_c1 = pickle.load(open(filename_c1, 'rb'))

# Application

df_nu_PU200_cl3d['cl3d_c1'] = model_c1.predict(df_nu_PU200_cl3d[['cl3d_abseta']])
df_nu_PU200_cl3d['cl3d_pt_c1'] = df_nu_PU200_cl3d.cl3d_c1 + df_nu_PU200_cl3d.cl3d_pt

print('Calibration 1: mean={0}, rms={1}'.format(
    df_nu_PU200_cl3d['cl3d_pt_c1'].mean(),
    df_nu_PU200_cl3d['cl3d_pt_c1'].std(),
))

###########

filename_c2 = '/data_CMS_upgrade/mperez/HGCal_data/Sep19/models/model_calib2.pkl'
model_c2 = pickle.load(open(filename_c2, 'rb'))

# Application

features = ['n_matched_cl3d', 'cl3d_abseta', 
  'cl3d_showerlength', 'cl3d_coreshowerlength', 
  'cl3d_firstlayer', 'cl3d_maxlayer', 
  'cl3d_szz', 'cl3d_seetot', 'cl3d_spptot', 'cl3d_srrtot', 'cl3d_srrmean',
  'cl3d_hoe', 'cl3d_meanz', 
  'cl3d_layer10', 'cl3d_layer50', 'cl3d_layer90', 
  'cl3d_ntc67', 'cl3d_ntc90']

df_nu_PU200_cl3d['cl3d_c2'] = model_c2.predict(df_nu_PU200_cl3d[features])
df_nu_PU200_cl3d['cl3d_pt_c2'] = df_nu_PU200_cl3d.cl3d_c2 * df_nu_PU200_cl3d.cl3d_pt_c1

print('Calibration 2: mean={0}, rms={1}'.format(
    df_nu_PU200_cl3d['cl3d_pt_c2'].mean(),
    df_nu_PU200_cl3d['cl3d_pt_c2'].std(),
))

###########

filename_c3 = '/data_CMS_upgrade/mperez/HGCal_data/Sep19/models/model_calib3.pkl'
model_c3 = pickle.load(open(filename_c3, 'rb'))

# Application

logpt1 = np.log(abs(df_nu_PU200_cl3d['cl3d_pt_c2']))
logpt2 = logpt1**2
logpt3 = logpt1**3
logpt4 = logpt1**4

df_nu_PU200_cl3d['cl3d_c3'] = model_c3.predict(np.vstack([logpt1, logpt2, logpt3, logpt4]).T)
df_nu_PU200_cl3d['cl3d_pt_c3'] = df_nu_PU200_cl3d.cl3d_pt_c2 / df_nu_PU200_cl3d.cl3d_c3

print('Calibration 3: mean={0}, rms={1}'.format(
    df_nu_PU200_cl3d['cl3d_pt_c3'].mean(),
    df_nu_PU200_cl3d['cl3d_pt_c3'].std(),
))

print '...Done'
print ' '

###########

# Decay mode

print 'Applying DM identification'

filename_dm = '/data_CMS_upgrade/mperez/HGCal_data/Sep19/models/model_DMid.pkl'
model_dm = pickle.load(open(filename_dm, 'rb'))

features_dm = [ 'cl3d_pt_c3',
    'cl3d_abseta', 'cl3d_showerlength',
    'cl3d_coreshowerlength', 'cl3d_firstlayer', 'cl3d_maxlayer', 'cl3d_szz',
    'cl3d_seetot', 'cl3d_spptot', 'cl3d_srrtot', 'cl3d_srrmean',
    'cl3d_hoe', 'cl3d_meanz', 'cl3d_layer10',
    'cl3d_layer50', 'cl3d_layer90', 'cl3d_ntc67', 'cl3d_ntc90' ]

inputs = df_nu_PU200_cl3d[features_dm]

df_nu_PU200_cl3d['cl3d_predDM'] = model_dm.predict(inputs)

# Predicted probabilities
probas = model_dm.predict_proba(inputs)
df_nu_PU200_cl3d['cl3d_probDM0'] = probas[:,0]
df_nu_PU200_cl3d['cl3d_probDM1'] = probas[:,1]
df_nu_PU200_cl3d['cl3d_probDM2'] = probas[:,2]

print '...Done'
print ' '

###########
###########

# RATES

df_nu_PU200_cl3d.set_index('event', inplace=True)
group = df_nu_PU200_cl3d.groupby('event')
#print 'all',np.shape(df_nu_PU200_cl3d['cl3d_pt_c3'])
#print df_nu_PU200_cl3d['cl3d_pt_c3']

df_nu_PU200_cl3d['leading_cl3d_pt_c3'] = group['cl3d_pt_c3'].max()
#print df_nu_PU200_cl3d['leading_cl3d_pt_c3']
df_nu_PU200_cl3d['cl3d_isleading'] = df_nu_PU200_cl3d['leading_cl3d_pt_c3'] == df_nu_PU200_cl3d['cl3d_pt_c3']
df_nu_PU200_cl3d['cl3d_isnotleading'] = df_nu_PU200_cl3d['leading_cl3d_pt_c3'] != df_nu_PU200_cl3d['cl3d_pt_c3']


print df_nu_PU200_cl3d['cl3d_pt_c3'].head(1000)
print df_nu_PU200_cl3d['cl3d_isleading'].head(1000)
print df_nu_PU200_cl3d['cl3d_isnotleading'].head(1000)

#sel = df_nu_PU200_cl3d['cl3d_isleading']==True
#df_nu_PU200_lead_cl3d = df_nu_PU200_cl3d[sel]
#print 'leading',np.shape(df_nu_PU200_lead_cl3d['cl3d_pt_c3'])

#sel2 = df_nu_PU200_cl3d['cl3d_isleading']==False
#df_nu_PU200_nonlead_cl3d = df_nu_PU200_cl3d[sel2]
#print 'nonleading ',df_nu_PU200_nonlead_cl3d['cl3d_pt_c3']
#df_nu_PU200_nonlead_cl3d['subleading_cl3d_pt_c3'] = group['cl3d_pt_c3'].max()
#print 'subleading ',np.shape(df_nu_PU200_nonlead_cl3d['cl3d_pt_c3'])


'''
sel_dm0 = df_nu_PU200_lead_cl3d['cl3d_predDM']==0.0
sel_dm1 = df_nu_PU200_lead_cl3d['cl3d_predDM']==1.0
sel_dm2 = df_nu_PU200_lead_cl3d['cl3d_predDM']==2.0

df_nu_PU200_lead_cl3d_dm0 = df_nu_PU200_lead_cl3d[sel_dm0]
df_nu_PU200_lead_cl3d_dm1 = df_nu_PU200_lead_cl3d[sel_dm1]
df_nu_PU200_lead_cl3d_dm2 = df_nu_PU200_lead_cl3d[sel_dm2]

df_nu_PU200_lead_cl3d.reset_index()['event']
df_nu_PU200_lead_cl3d_dm0.reset_index()['event']
df_nu_PU200_lead_cl3d_dm1.reset_index()['event']
df_nu_PU200_lead_cl3d_dm2.reset_index()['event']

df_nu_PU200_lead_cl3d.to_root(fileout_rates_inclusive, key=treename, store_index=True, mode='a') #events stored in branch __index__event
df_nu_PU200_lead_cl3d_dm0.to_root(fileout_rates_DM0, key=treename, store_index=True, mode='a')
df_nu_PU200_lead_cl3d_dm1.to_root(fileout_rates_DM1, key=treename, store_index=True, mode='a')
df_nu_PU200_lead_cl3d_dm2.to_root(fileout_rates_DM2, key=treename, store_index=True, mode='a')
'''
