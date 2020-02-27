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
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import LinearRegression
import pickle

###########

dir_in = '/data_CMS_upgrade/mperez/HGCal_data/Sep19/pubdt/'
filein_tau_PU200 = dir_in+'RelValDiTau_Pt20To100_v10_PU200_pubdt.hdf5'

store_tau_PU200 = pd.HDFStore(filein_tau_PU200, mode='r')
df_tau_PU200	= store_tau_PU200['df_tau_PU200']
store_tau_PU200.close()

dir_out = '/data_CMS_upgrade/mperez/HGCal_data/Sep19/calibrated/'
fileout_tau_PU200 = dir_out+'RelValDiTau_Pt20To100_v10_PU200_calibrated.hdf5'

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

###########

df_tau_PU200['cl3d_response'] = df_tau_PU200['cl3d_pt']/df_tau_PU200['gentau_vis_pt']

print '-- Raw --'

print('Raw: mean={0}, rms={1}, rms/mean={2}'.format(
    df_tau_PU200['cl3d_response'].mean(),
    df_tau_PU200['cl3d_response'].std(),
    df_tau_PU200['cl3d_response'].std()/df_tau_PU200['cl3d_response'].mean()
))

print ' '

###########

print '-- Calibration 1 --'

# Training

input_c1 = df_tau_PU200[['cl3d_abseta']]
target_c1 = df_tau_PU200.gentau_vis_pt - df_tau_PU200.cl3d_pt

model_c1 = LinearRegression().fit(input_c1, target_c1)

with open('/data_CMS_upgrade/mperez/HGCal_data/Sep19/models/model_calib1.pkl', 'wb') as f:
  pickle.dump(model_c1, f)

# Application

df_tau_PU200['cl3d_c1'] = model_c1.predict(df_tau_PU200[['cl3d_abseta']])
df_tau_PU200['cl3d_pt_c1'] = df_tau_PU200.cl3d_c1 + df_tau_PU200.cl3d_pt
df_tau_PU200['cl3d_response_c1'] = df_tau_PU200.cl3d_pt_c1 / df_tau_PU200.gentau_vis_pt

print('Calib1: mean={0}, rms={1}, rms/mean={2}'.format(
    df_tau_PU200['cl3d_response_c1'].mean(),
    df_tau_PU200['cl3d_response_c1'].std(),
    df_tau_PU200['cl3d_response_c1'].std()/df_tau_PU200['cl3d_response_c1'].mean()
))

print ' '

###########

print '-- Calibration 2 --'

# Training

features = ['n_matched_cl3d', 'cl3d_abseta', 
	'cl3d_showerlength', 'cl3d_coreshowerlength', 
	'cl3d_firstlayer', 'cl3d_maxlayer', 
	'cl3d_szz', 'cl3d_seetot', 'cl3d_spptot', 'cl3d_srrtot', 'cl3d_srrmean',
    'cl3d_hoe', 'cl3d_meanz', 
    'cl3d_layer10', 'cl3d_layer50', 'cl3d_layer90', 
    'cl3d_ntc67', 'cl3d_ntc90']

input_c2 = df_tau_PU200[features]
target_c2 = df_tau_PU200.gentau_vis_pt / df_tau_PU200.cl3d_pt_c1

model_c2 = GradientBoostingRegressor(n_estimators=1000, learning_rate=0.1,
                                max_depth=2, random_state=0, loss='huber').fit(input_c2, target_c2)

with open('/data_CMS_upgrade/mperez/HGCal_data/Sep19/models/model_calib2.pkl', 'wb') as f:
  pickle.dump(model_c2, f)

# Application
df_tau_PU200['cl3d_c2'] = model_c2.predict(df_tau_PU200[features])
df_tau_PU200['cl3d_pt_c2'] = df_tau_PU200.cl3d_c2 * df_tau_PU200.cl3d_pt_c1
df_tau_PU200['cl3d_response_c2'] = df_tau_PU200.cl3d_pt_c2 / df_tau_PU200.gentau_vis_pt

print('Calib2: mean={0}, rms={1}, rms/mean={2}'.format(
    df_tau_PU200['cl3d_response_c2'].mean(),
    df_tau_PU200['cl3d_response_c2'].std(),
    df_tau_PU200['cl3d_response_c2'].std()/df_tau_PU200['cl3d_response_c2'].mean()
))

print ' '

###########

print '-- Calibration 3 --'

# Training

df_tau_PU200['gentau_vis_abseta'] = np.abs(df_tau_PU200['gentau_vis_eta'])

ptcut=20
etamin=1.6
etamax=2.9

df_tau_PU200['gentau_vis_bin_eta'] = ((df_tau_PU200['gentau_vis_abseta'] - etamin)/0.1).astype('int32')
df_tau_PU200['gentau_vis_bin_pt']  = ((df_tau_PU200['gentau_vis_pt'] - ptcut)/5).astype('int32')

# Training

vars = ['gentau_vis_pt', 'gentau_vis_bin_pt',
        'cl3d_pt_c2', 'cl3d_response_c2']

pt_means = df_tau_PU200[vars].groupby('gentau_vis_bin_pt').mean() 
pt_rmss = df_tau_PU200[vars].groupby('gentau_vis_bin_pt').std() 

pt_means['logpt1'] = np.log(pt_means['cl3d_pt_c2'])
pt_means['logpt2'] = pt_means.logpt1**2
pt_means['logpt3'] = pt_means.logpt1**3
pt_means['logpt4'] = pt_means.logpt1**4

input_c3 = pt_means[['logpt1', 'logpt2', 'logpt3', 'logpt4']]
target_c3 = pt_means['cl3d_response_c2']

model_c3 = LinearRegression().fit(input_c3, target_c3)

with open('/data_CMS_upgrade/mperez/HGCal_data/Sep19/models/model_calib3.pkl', 'wb') as f:
  pickle.dump(model_c3, f)

# Application

logpt1 = np.log(abs(df_tau_PU200['cl3d_pt_c2']))
logpt2 = logpt1**2
logpt3 = logpt1**3
logpt4 = logpt1**4

df_tau_PU200['cl3d_c3'] = model_c3.predict(np.vstack([logpt1, logpt2, logpt3, logpt4]).T)
df_tau_PU200['cl3d_pt_c3'] = df_tau_PU200.cl3d_pt_c2 / df_tau_PU200.cl3d_c3
df_tau_PU200['cl3d_response_c3'] = df_tau_PU200.cl3d_pt_c3 / df_tau_PU200.gentau_vis_pt

print('Calib3: mean={0}, rms={1}, rms/mean={2}'.format(
    df_tau_PU200['cl3d_response_c3'].mean(),
    df_tau_PU200['cl3d_response_c3'].std(),
    df_tau_PU200['cl3d_response_c3'].std()/df_tau_PU200['cl3d_response_c3'].mean()
))

print ' '

###########

# SAVE DFS

# Save files
#store_tau_PU200 = pd.HDFStore(fileout_tau_PU200, mode='w')
#store_tau_PU200['df_tau_PU200'] = df_tau_PU200
#store_tau_PU200.close()

###########

# PLOTTING

plotdir = '/home/llr/cms/mperez/HGCal/v10_geometry/tau_algorithm/plots/'

import geeksw.plotting.cmsplot as plt
from geeksw.plotting.root_colors import *

plt.matplotlib.font_manager._rebuild()

# INCLUSIVE RESPONSES
'''
matplotlib.rcParams.update({'font.size': 22})
plt.figure(figsize=(15,10))
plt.hist(df_tau_PU200['cl3d_response'], bins=np.arange(0., 2., 0.02), alpha=0.3, label='Before calib', color='blue')
plt.hist(df_tau_PU200['cl3d_response_c1'], bins=np.arange(0., 2., 0.02), alpha=0.3, label='After $\eta$ corr', color='red')
plt.hist(df_tau_PU200['cl3d_response_c2'], bins=np.arange(0., 2., 0.02), alpha=0.3, label='After BDT', color='green')
plt.hist(df_tau_PU200['cl3d_response_c3'], bins=np.arange(0., 2., 0.02), label='After $p_T$ corr', color='black', histtype='step', lw=2)
plt.legend(loc = 'upper left', fontsize=22)
plt.xlabel(r'$p_{T}^{L1}/p_{T}^{gen}$')
plt.ylabel(r'Entries')
plt.ylim(0, 1000)
plt.savefig(plotdir+'calib_responses.png')
plt.savefig(plotdir+'calib_responses.pdf')
'''
'''
#TDR
plt.figure(figsize=(8,7))
plt.hist(df_tau_PU200['cl3d_response'], bins=np.arange(0., 1.8, 0.03), label='Raw', color='red', histtype='step', lw=2)
plt.hist(df_tau_PU200['cl3d_response_c3'], bins=np.arange(0., 1.8, 0.03), label='Calibrated', color='blue', histtype='step', lw=2)
plt.legend(loc = 'upper right')
plt.grid()
plt.xlabel(r'$p_{T}^{L1,\tau}\ /\ p_{T}^{gen,\tau}$')
plt.ylabel(r'a. u.')
#plt.ylim(0, 1000)
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("200 PU","HGCal calorimeter-only")
plt.subplots_adjust(bottom=0.12)
plt.savefig(plotdir+'calib_responses_TDR2.png')
plt.savefig(plotdir+'calib_responses_TDR2.pdf')
plt.show()
'''

# MODEL CALIBRATION 3
'''
matplotlib.rcParams.update({'font.size': 22})
plt.figure(figsize=(15,10))
pt = np.arange(25, 100, 1)
logpt1 = np.log(pt)
logpt2 = logpt1**2
logpt3 = logpt1**3
logpt4 = logpt1**4	
plt.plot(pt_means['cl3d_pt_c2'], pt_means['cl3d_response_c2'], marker='s', markersize=8, ls='None', color='blue')
plt.plot(pt, model_c3.predict(np.vstack([logpt1, logpt2, logpt3, logpt4]).T), ls='--', color='blue')	
plt.grid()
plt.xlabel(r'$\langle p_{T}^{L1}\rangle\, [GeV]$')
plt.ylabel(r'$\langle p_{T}^{L1}/p_{T}^{gen}\rangle$')
plt.savefig(plotdir+'calib_modelc3.png')
plt.savefig(plotdir+'calib_modelc3.pdf')
'''
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


plot_var = ['gentau_vis_pt', 'gentau_vis_abseta',
	        'gentau_vis_bin_eta', 'gentau_vis_bin_pt',
	        'cl3d_pt', 'cl3d_response', 'cl3d_abseta',
	        'cl3d_pt_c1', 'cl3d_response_c1',
	        'cl3d_pt_c2', 'cl3d_response_c2',
	        'cl3d_pt_c3', 'cl3d_response_c3']

ptcut   = 20
etamin  = 1.6
etamax  = 2.9

df_tau_PU200['gentau_vis_bin_eta'] = ((df_tau_PU200['gentau_vis_abseta'] - etamin)/0.1).astype('int32')
df_tau_PU200['gentau_vis_bin_pt']  = ((df_tau_PU200['gentau_vis_pt'] - ptcut)/5).astype('int32')

etameans    = df_tau_PU200[plot_var].groupby('gentau_vis_bin_eta').mean()
etarmss     = df_tau_PU200[plot_var].groupby('gentau_vis_bin_eta').std()
etaeffrmss  = df_tau_PU200[plot_var].groupby('gentau_vis_bin_eta').apply(effrms)
ptmeans     = df_tau_PU200[plot_var].groupby('gentau_vis_bin_pt').mean()
ptrmss      = df_tau_PU200[plot_var].groupby('gentau_vis_bin_pt').std()
pteffrmss   = df_tau_PU200[plot_var].groupby('gentau_vis_bin_pt').apply(effrms)
'''
plt.figure(figsize=(20,10))
df = etameans
plt.plot(df['gentau_vis_abseta'], df['cl3d_response_c3'])
plt.ylim(0.9, 1.2)
plt.legend(loc = 'upper right', fontsize=16)
plt.xlabel(r'$|\eta|$')
plt.ylabel(r'$\langle p_{T}^{L1}/p_{T}^{gen}\rangle$')
plt.grid()
plt.ylim(0.9, 1.1)
plt.savefig(plotdir+'calib_mean_vs_eta.png')
plt.savefig(plotdir+'calib_mean_vs_eta.pdf')

plt.figure(figsize=(15,10))
df = etameans
plt.plot(df['gentau_vis_abseta'], etarmss['cl3d_response_c3'])
plt.ylim(0., 0.4)
plt.legend(loc = 'upper right', fontsize=16)
plt.xlabel(r'$|\eta|$')
plt.ylabel(r'$RMS(p_{T}^{L1}/p_{T}^{gen})$')
plt.grid()
plt.ylim(0.0, 0.4)
plt.savefig(plotdir+'calib_rms_vs_eta.png')
plt.savefig(plotdir+'calib_rms_vs_eta.pdf')
	
plt.figure(figsize=(15,10))
df = etameans
plt.plot(df['gentau_vis_abseta'], etaeffrmss['cl3d_response_c3'])
plt.ylim(0., 0.4)
plt.legend(loc = 'upper right', fontsize=16)
plt.xlabel(r'$|\eta|$')
plt.ylabel(r'$RMS_{eff}(p_{T}^{L1}/p_{T}^{gen})$')
plt.grid()
plt.ylim(0.0, 0.4)
plt.savefig(plotdir+'calib_effrms_vs_eta.png')
plt.savefig(plotdir+'calib_effrms_vs_eta.pdf')
	
plt.figure(figsize=(15,10))
df = ptmeans
plt.plot(df['gentau_vis_pt'], df['cl3d_response_c3'])
plt.ylim(0.9, 1.1)
plt.legend(loc = 'upper right', fontsize=16)
plt.xlabel(r'$p_{T}\, [GeV]$', fontsize=18)
plt.ylabel(r'$\langle p_{T}^{L1}/p_{T}^{gen}\rangle$')
plt.grid()
plt.ylim(0.9, 1.1)
plt.savefig(plotdir+'calib_mean_vs_pt.png')
plt.savefig(plotdir+'calib_mean_vs_pt.pdf')

plt.figure(figsize=(15,10))
df = ptmeans
plt.plot(df['gentau_vis_pt'], ptrmss['cl3d_response_c3'])
plt.ylim(0., 0.3)
#plt.xlim(5., 100)
plt.legend(loc = 'upper right', fontsize=16)
plt.xlabel(r'$p_{T}\, [GeV]$')
plt.ylabel(r'$RMS(p_{T}^{L1}/p_{T}^{gen})$')
plt.grid()
plt.ylim(0.0, 0.4)
plt.savefig(plotdir+'calib_rms_vs_pt.png')
plt.savefig(plotdir+'calib_rms_vs_pt.pdf')
	
plt.figure(figsize=(15,10))
df = ptmeans
plt.plot(df['gentau_vis_pt'], pteffrmss['cl3d_response_c3'])
plt.ylim(0., 0.3)
plt.legend(loc = 'upper right', fontsize=16)
plt.xlabel(r'$p_{T}\, [GeV]$')
plt.ylabel(r'$RMS_{eff}(p_{T}^{L1}/p_{T}^{gen})$')
plt.grid()
plt.ylim(0.0, 0.4)
plt.savefig(plotdir+'calib_effrms_vs_pt.png')
plt.savefig(plotdir+'calib_effrms_vs_pt.pdf')
'''
#TDR
df = ptmeans

'''plt.figure(figsize=(8,7))
plt.plot(df['gentau_vis_pt'], pteffrmss['cl3d_response']/df['cl3d_response'], label="Raw",color='red',lw=2)
plt.plot(df['gentau_vis_pt'], pteffrmss['cl3d_response_c3']/df['cl3d_response_c3'], label="Calibrated",color='blue',lw=2)
plt.ylim(0.01, 0.4)
plt.legend(loc = 'lower left')
plt.xlabel(r'$p_{T}^{gen,\tau}\ [GeV]$')
plt.ylabel(r'$RMS_{eff}\ /\ <p_{T}^{L1,\tau}\ /\ p_{T}^{gen,\tau}>$')
plt.grid()
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("200 PU","HGCal calorimeter-only")
plt.subplots_adjust(bottom=0.12)
plt.subplots_adjust(left=0.14)
#plt.ylim(0.0, 0.5)
plt.savefig(plotdir+'calib_effrms_vs_pt_TDR2.png')
plt.savefig(plotdir+'calib_effrms_vs_pt_TDR2.pdf')
'''


import matplotlib.lines as mlines

plt.figure(figsize=(8,7))

x_array = [22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 82.5, 87.5, 92.5, 97.5]
y_array_raw = (pteffrmss['cl3d_response']/df['cl3d_response']).values
y_array_calib = (pteffrmss['cl3d_response_c3']/df['cl3d_response_c3']).values

plt.figure(figsize=(8,7))
plt.errorbar(x_array,y_array_raw,xerr=2.5,yerr=None,marker="o",mec='red',ls='None',label='Raw',color='red',lw=2)
plt.errorbar(x_array,y_array_calib,xerr=2.5,yerr=None,marker="o",mec='blue',ls='None',label='Calibrated',color='blue',lw=2)
plt.ylim(0.01, 0.4)
red_line = mlines.Line2D([], [], color='red',markersize=15, label='Raw',lw=2)
blue_line = mlines.Line2D([], [], color='blue',markersize=15, label='Calibrated',lw=2)
plt.legend(loc = 'lower left',handles=[red_line,blue_line])
plt.xlabel(r'$p_{T}^{gen,\tau}\ [GeV]$')
plt.ylabel(r'$RMS_{eff}\ /\ <p_{T}^{L1,\tau}\ /\ p_{T}^{gen,\tau}>$')
plt.grid()
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("200 PU","HGCal calorimeter-only")
plt.subplots_adjust(bottom=0.12)
plt.subplots_adjust(left=0.14)
txt = (r'1.6 < | $\eta_{gen,\tau}$ | < 2.9')
t = plt.text(70,0.32, txt, ha='left', wrap=True, fontsize=16)
t.set_bbox(dict(facecolor='white', edgecolor='white'))
#plt.ylim(0.0, 0.5)
plt.savefig(plotdir+'calib_effrms_vs_pt_TDR4.png')
plt.savefig(plotdir+'calib_effrms_vs_pt_TDR4.pdf')

'''
plt.figure(figsize=(8,7))
plt.hist(df_tau_PU200['gentau_vis_eta']-df_tau_PU200['cl3d_eta'], bins=np.arange(-0.15, 0.15, 0.005), color='blue',lw=2,histtype='step')
#plt.ylim(0.01, 0.4)
plt.xlabel(r'$\eta_{gen,\tau}$ - $\eta_{L1,\tau}$')
plt.ylabel(r'a. u. ')
#plt.grid()
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("200 PU","HGCal calorimeter-only")
plt.subplots_adjust(bottom=0.12)
plt.subplots_adjust(left=0.14)
txt = (r'1.6 < | $\eta_{gen,\tau}$ | < 2.9')
t = plt.text(-0.13,1400, txt, ha='left', wrap=True, fontsize=16)
t.set_bbox(dict(facecolor='white', edgecolor='white'))
#plt.ylim(0.0, 0.5)
plt.savefig(plotdir+'res_eta_TDR4.png')
plt.savefig(plotdir+'res_eta_TDR4.pdf')

plt.figure(figsize=(8,7))
plt.hist(df_tau_PU200['gentau_vis_phi']-df_tau_PU200['cl3d_phi'], bins=np.arange(-0.15, 0.15, 0.005), color='blue',lw=2,histtype='step')
plt.ylim(0, 2200)
plt.xlabel(r'$\phi_{gen,\tau}$ - $\phi_{L1,\tau}$')
plt.ylabel(r'a. u. ')
#plt.grid()
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("200 PU","HGCal calorimeter-only")
plt.subplots_adjust(bottom=0.12)
plt.subplots_adjust(left=0.14)
t = plt.text(-0.13,1400, txt, ha='left', wrap=True, fontsize=16)
t.set_bbox(dict(facecolor='white', edgecolor='white'))
#plt.ylim(0.0, 0.5)
plt.savefig(plotdir+'res_phi_TDR4.png')
plt.savefig(plotdir+'res_phi_TDR4.pdf')
'''
