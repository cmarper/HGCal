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

dir_in = '/data_CMS_upgrade/mperez/HGCal_data/Sep19/calibrated/'
filein_tau_PU200 = dir_in+'RelValDiTau_Pt20To100_v10_PU200_calibrated.hdf5'

store_tau_PU200 = pd.HDFStore(filein_tau_PU200, mode='r')
df_tau_PU200	= store_tau_PU200['df_tau_PU200']
store_tau_PU200.close()

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

seldm0 = df_tau_PU200['gentau_decayMode'] == 0
df_tau_PU200_DM0 = df_tau_PU200[seldm0]

seldm1 = df_tau_PU200['gentau_decayMode'] == 1
df_tau_PU200_DM1 = df_tau_PU200[seldm1]

seldm4 = df_tau_PU200['gentau_decayMode'] == 4
df_tau_PU200_DM4 = df_tau_PU200[seldm4]

seldm5 = df_tau_PU200['gentau_decayMode'] == 5
df_tau_PU200_DM5 = df_tau_PU200[seldm5]

seldm01 = (df_tau_PU200['gentau_decayMode'] == 0) | (df_tau_PU200['gentau_decayMode'] == 1)
df_tau_PU200_DM01 = df_tau_PU200[seldm01]

seldm45 = (df_tau_PU200['gentau_decayMode'] == 4) | (df_tau_PU200['gentau_decayMode'] == 5)
df_tau_PU200_DM45 = df_tau_PU200[seldm45]

###########

plotdir = '/home/llr/cms/mperez/HGCal/v10_geometry/tau_algorithm/plots/'
import geeksw.plotting.cmsplot as plt
from geeksw.plotting.root_colors import *

plt.matplotlib.font_manager._rebuild()

df_tau_PU200['gentau_vis_abseta'] = np.abs(df_tau_PU200['gentau_vis_eta'])
df_tau_PU200['gentau_bin_eta'] = ((df_tau_PU200['gentau_vis_abseta'] - etamin)/0.1).astype('int32')
df_tau_PU200['gentau_bin_pt']  = ((df_tau_PU200['gentau_vis_pt'] - ptcut)/2).astype('int32')

df_tau_PU200_DM0['gentau_vis_abseta'] = np.abs(df_tau_PU200_DM0['gentau_vis_eta'])
df_tau_PU200_DM0['gentau_bin_eta'] = ((df_tau_PU200_DM0['gentau_vis_abseta'] - etamin)/0.1).astype('int32')
df_tau_PU200_DM0['gentau_bin_pt']  = ((df_tau_PU200_DM0['gentau_vis_pt'] - ptcut)/2).astype('int32')

df_tau_PU200_DM1['gentau_vis_abseta'] = np.abs(df_tau_PU200_DM1['gentau_vis_eta'])
df_tau_PU200_DM1['gentau_bin_eta'] = ((df_tau_PU200_DM1['gentau_vis_abseta'] - etamin)/0.1).astype('int32')
df_tau_PU200_DM1['gentau_bin_pt']  = ((df_tau_PU200_DM1['gentau_vis_pt'] - ptcut)/2).astype('int32')

df_tau_PU200_DM4['gentau_vis_abseta'] = np.abs(df_tau_PU200_DM4['gentau_vis_eta'])
df_tau_PU200_DM4['gentau_bin_eta'] = ((df_tau_PU200_DM4['gentau_vis_abseta'] - etamin)/0.1).astype('int32')
df_tau_PU200_DM4['gentau_bin_pt']  = ((df_tau_PU200_DM4['gentau_vis_pt'] - ptcut)/2).astype('int32')

df_tau_PU200_DM5['gentau_vis_abseta'] = np.abs(df_tau_PU200_DM5['gentau_vis_eta'])
df_tau_PU200_DM5['gentau_bin_eta'] = ((df_tau_PU200_DM5['gentau_vis_abseta'] - etamin)/0.1).astype('int32')
df_tau_PU200_DM5['gentau_bin_pt']  = ((df_tau_PU200_DM5['gentau_vis_pt'] - ptcut)/2).astype('int32')

df_tau_PU200_DM01['gentau_vis_abseta'] = np.abs(df_tau_PU200_DM01['gentau_vis_eta'])
df_tau_PU200_DM01['gentau_bin_eta'] = ((df_tau_PU200_DM01['gentau_vis_abseta'] - etamin)/0.1).astype('int32')
df_tau_PU200_DM01['gentau_bin_pt']  = ((df_tau_PU200_DM01['gentau_vis_pt'] - ptcut)/2).astype('int32')

df_tau_PU200_DM45['gentau_vis_abseta'] = np.abs(df_tau_PU200_DM45['gentau_vis_eta'])
df_tau_PU200_DM45['gentau_bin_eta'] = ((df_tau_PU200_DM45['gentau_vis_abseta'] - etamin)/0.1).astype('int32')
df_tau_PU200_DM45['gentau_bin_pt']  = ((df_tau_PU200_DM45['gentau_vis_pt'] - ptcut)/2).astype('int32')


def efficiency(group, threshold):
	tot = group.shape[0]
	sel = group[(group.cl3d_pt_c3 > threshold)].shape[0]
	#print 'thr ',threshold, 'eff ', float(sel)/float(tot)
	return float(sel)/float(tot)


efficiencies_vs_pt = {}
efficiencies_vs_pt = df_tau_PU200.groupby('gentau_bin_pt').mean()

efficiencies_vs_pt_DM0 = {}
efficiencies_vs_pt_DM0 = df_tau_PU200_DM0.groupby('gentau_bin_pt').mean()

efficiencies_vs_pt_DM1 = {}
efficiencies_vs_pt_DM1 = df_tau_PU200_DM1.groupby('gentau_bin_pt').mean()

efficiencies_vs_pt_DM4 = {}
efficiencies_vs_pt_DM4 = df_tau_PU200_DM4.groupby('gentau_bin_pt').mean()

efficiencies_vs_pt_DM5 = {}
efficiencies_vs_pt_DM5 = df_tau_PU200_DM5.groupby('gentau_bin_pt').mean()

efficiencies_vs_pt_DM01 = {}
efficiencies_vs_pt_DM01 = df_tau_PU200_DM01.groupby('gentau_bin_pt').mean()

efficiencies_vs_pt_DM45 = {}
efficiencies_vs_pt_DM45 = df_tau_PU200_DM45.groupby('gentau_bin_pt').mean()

turnon_thresholds = range(10, 100, 10)

for threshold in turnon_thresholds:

	eff = df_tau_PU200.groupby('gentau_bin_pt').apply(lambda x : efficiency(x, threshold))
        eff_DM0 = df_tau_PU200_DM0.groupby('gentau_bin_pt').apply(lambda x : efficiency(x, threshold))
        eff_DM1 = df_tau_PU200_DM1.groupby('gentau_bin_pt').apply(lambda x : efficiency(x, threshold))
        eff_DM4 = df_tau_PU200_DM4.groupby('gentau_bin_pt').apply(lambda x : efficiency(x, threshold))
        eff_DM5 = df_tau_PU200_DM5.groupby('gentau_bin_pt').apply(lambda x : efficiency(x, threshold))
        eff_DM01 = df_tau_PU200_DM01.groupby('gentau_bin_pt').apply(lambda x : efficiency(x, threshold))
        eff_DM45 = df_tau_PU200_DM45.groupby('gentau_bin_pt').apply(lambda x : efficiency(x, threshold))        

	eff_smooth = eff.rolling(window=5, win_type='triang', center=True).mean()
        eff_smooth_DM0 = eff_DM0.rolling(window=5, win_type='triang', center=True).mean()
        eff_smooth_DM1 = eff_DM1.rolling(window=5, win_type='triang', center=True).mean()
        eff_smooth_DM4 = eff_DM4.rolling(window=5, win_type='triang', center=True).mean()
        eff_smooth_DM5 = eff_DM5.rolling(window=5, win_type='triang', center=True).mean()
        eff_smooth_DM01 = eff_DM01.rolling(window=7, win_type='triang', center=True).mean()
        eff_smooth_DM45 = eff_DM45.rolling(window=7, win_type='triang', center=True).mean()

	eff_smooth.fillna(0., inplace=True)
        eff_smooth_DM0.fillna(0., inplace=True)
        eff_smooth_DM1.fillna(0., inplace=True)
        eff_smooth_DM4.fillna(0., inplace=True)
        eff_smooth_DM5.fillna(0., inplace=True)
        eff_smooth_DM01.fillna(0., inplace=True)
        eff_smooth_DM45.fillna(0., inplace=True)

	efficiencies_vs_pt['efficiency_{}'.format(threshold)] = eff
        efficiencies_vs_pt_DM0['efficiency_{}'.format(threshold)] = eff_DM0
        efficiencies_vs_pt_DM1['efficiency_{}'.format(threshold)] = eff_DM1
        efficiencies_vs_pt_DM4['efficiency_{}'.format(threshold)] = eff_DM4
        efficiencies_vs_pt_DM5['efficiency_{}'.format(threshold)] = eff_DM5
        efficiencies_vs_pt_DM01['efficiency_{}'.format(threshold)] = eff_DM01
        efficiencies_vs_pt_DM45['efficiency_{}'.format(threshold)] = eff_DM45

	efficiencies_vs_pt['efficiency_smooth_{}'.format(threshold)] = eff_smooth
        efficiencies_vs_pt_DM0['efficiency_smooth_{}'.format(threshold)] = eff_smooth_DM0
        efficiencies_vs_pt_DM1['efficiency_smooth_{}'.format(threshold)] = eff_smooth_DM1
        efficiencies_vs_pt_DM4['efficiency_smooth_{}'.format(threshold)] = eff_smooth_DM4
        efficiencies_vs_pt_DM5['efficiency_smooth_{}'.format(threshold)] = eff_smooth_DM5
	efficiencies_vs_pt_DM01['efficiency_smooth_{}'.format(threshold)] = eff_smooth_DM01
	efficiencies_vs_pt_DM45['efficiency_smooth_{}'.format(threshold)] = eff_smooth_DM45

'''
matplotlib.rcParams.update({'font.size': 22})
plt.figure(figsize=(20,15))
df = efficiencies_vs_pt
for threshold in turnon_thresholds:
	eff = df['efficiency_{}'.format(threshold)]
	eff_smooth = df['efficiency_smooth_{}'.format(threshold)]
	line, = plt.plot(df.gentau_vis_pt, eff_smooth, label=threshold, ls='-', linewidth=2)
	#plt.plot(df.gentau_vis_pt, eff, color=line.get_color(), ls='--')
plt.ylim(0., 1.01)
plt.xlim(20, 95)
plt.legend(loc = 'lower right', title='L1 thresholds', fontsize=16)
plt.xlabel(r'Gen. tau $p_{T}\,[GeV]$')
plt.ylabel('Efficiency')
plt.grid()
plt.savefig(plotdir+'eff_vs_pt.png')
plt.savefig(plotdir+'eff_vs_pt.pdf')
'''
#TDR
'''
plt.figure(figsize=(8,7))
df = efficiencies_vs_pt
for threshold in turnon_thresholds:
        lab = "$E_{T}^{L1}$ > %s GeV" % threshold
        eff = df['efficiency_{}'.format(threshold)]
        eff_smooth = df['efficiency_smooth_{}'.format(threshold)]
        line, = plt.plot(df.gentau_vis_pt, eff_smooth, label=lab, ls='-', linewidth=2)
        #plt.plot(df.gentau_vis_pt, eff, color=line.get_color(), ls='--')
plt.ylim(0., 1.01)
plt.xlim(20, 95)
plt.legend(loc = 'lower right', fontsize=16)
plt.xlabel(r'$p_{T}^{gen,\tau}\ [GeV]$')
plt.ylabel('Efficiency')
plt.grid()
plt.ylim(0.,1.14)
plt.subplots_adjust(bottom=0.12)
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("200 PU","HGCal calorimeter-only")
plt.savefig(plotdir+'eff_vs_pt_TDR2.png')
plt.savefig(plotdir+'eff_vs_pt_TDR2.pdf')
'''
'''
plt.rcParams['legend.numpoints'] = 1

import matplotlib.lines as mlines
plt.figure(figsize=(8,7))
df = efficiencies_vs_pt
lab_40 = r"$E_{T}^{L1,\tau}$ > 40 GeV"
lab_50 = r"$E_{T}^{L1,\tau}$ > 50 GeV"
lab_60 = r"$E_{T}^{L1,\tau}$ > 60 GeV"
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
plt.xlim(20, 94.5)
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
txt = (r'1.6 < | $\eta_{gen,\tau}$ | < 2.9')
t = plt.text(67,0.45, txt, ha='left', wrap=True, fontsize=16)
t.set_bbox(dict(facecolor='white', edgecolor='white'))
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("200 PU","HGCal calorimeter-only")
plt.savefig(plotdir+'eff_vs_pt_TDR4.png')
plt.savefig(plotdir+'eff_vs_pt_TDR4.pdf')
'''
'''
import matplotlib.lines as mlines
plt.figure(figsize=(8,7))

df_DM0 = efficiencies_vs_pt_DM0
df_DM1 = efficiencies_vs_pt_DM1
df_DM4 = efficiencies_vs_pt_DM4
df_DM5 = efficiencies_vs_pt_DM5

lab_50 = r"$E_{T}^{L1,\tau}$ > 50 GeV"

eff_DM0 = df_DM0['efficiency_50']
eff_DM1 = df_DM1['efficiency_50']
eff_DM4 = df_DM4['efficiency_50']
eff_DM5 = df_DM5['efficiency_50']

eff_smooth_DM0 = df_DM0['efficiency_smooth_50']
eff_smooth_DM1 = df_DM1['efficiency_smooth_50']
eff_smooth_DM4 = df_DM4['efficiency_smooth_50']
eff_smooth_DM5 = df_DM5['efficiency_smooth_50']

x_DM0 = efficiencies_vs_pt_DM0.gentau_vis_pt
x_DM1 = efficiencies_vs_pt_DM1.gentau_vis_pt
x_DM4 = efficiencies_vs_pt_DM4.gentau_vis_pt
x_DM5 = efficiencies_vs_pt_DM5.gentau_vis_pt

plt.errorbar(x_DM0,eff_DM0,xerr=1,ls='None',label='1-prong',color='red',lw=2,marker='o',mec='red')
plt.errorbar(x_DM1,eff_DM1,xerr=1,ls='None',label='1-prong + $\pi^{0}$\'s',color='blue',lw=2,marker='o',mec='blue')
plt.errorbar(x_DM4,eff_DM4,xerr=1,ls='None',label='3-prongs',color='green',lw=2,marker='o',mec='green')
plt.errorbar(x_DM5,eff_DM5,xerr=1,ls='None',label='3-prongs + $\pi^{0}$\'s',color='fuchsia',lw=2,marker='o',mec='fuchsia')

plt.plot(x_DM0,eff_smooth_DM0,color='red',lw=1.5)
plt.plot(x_DM1,eff_smooth_DM1,color='blue',lw=1.5)
plt.plot(x_DM4,eff_smooth_DM4,color='green',lw=1.5)
plt.plot(x_DM5,eff_smooth_DM5,color='fuchsia',lw=1.5)

plt.xlim(20, 94.5)
plt.ylim(0., 1.10)

red_line = mlines.Line2D([], [], color='red',markersize=15, label='1-prong',lw=2)
blue_line = mlines.Line2D([], [], color='blue',markersize=15, label='1-prong + $\pi^{0}$\'s',lw=2)
green_line = mlines.Line2D([], [], color='green',markersize=15, label='3-prongs',lw=2)
fuchsia_line = mlines.Line2D([], [], color='fuchsia',markersize=15, label='3-prongs + $\pi^{0}$\'s',lw=2)

plt.legend(title=r'$\tau$ decay modes',loc = 'lower right', fontsize=16, handles=[red_line,blue_line,green_line,fuchsia_line])
#plt.legend(loc = 'lower right', fontsize=16)
txt = (r'Gen. $\tau$ decay mode:')
t = plt.text(67,0.35, txt, ha='left', wrap=True, fontsize=16)
t.set_bbox(dict(facecolor='white', edgecolor='white'))
txt2 = (r'$E_{T}^{L1,\tau}$ > 50 GeV')
t2 = plt.text(28,0.83, txt2, ha='left', wrap=True, fontsize=16)
t2.set_bbox(dict(facecolor='white', edgecolor='white'))
txt3 = (r'1.6 < | $\eta_{gen,\tau}$ | < 2.9')
t3 = plt.text(-0.13,1400, txt3, ha='left', wrap=True, fontsize=16)
t3.set_bbox(dict(facecolor='white', edgecolor='white'))
plt.xlabel(r'$p_{T}^{gen,\tau}\ [GeV]$')
plt.ylabel('Efficiency')
plt.grid()
plt.subplots_adjust(bottom=0.12)
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("200 PU","HGCal calorimeter-only")
plt.savefig(plotdir+'eff_vs_pt_DM_TDR4.png')
plt.savefig(plotdir+'eff_vs_pt_DM_TDR4.pdf')
'''
'''
plt.figure(figsize=(8,7))
df_DM0 = efficiencies_vs_pt_DM0
df_DM1 = efficiencies_vs_pt_DM1
df_DM4 = efficiencies_vs_pt_DM4
df_DM5 = efficiencies_vs_pt_DM5
lab = "$E_{T}^{L1}$ > %s GeV" % threshold
eff_DM0 = df_DM0['efficiency_{}'.format(threshold)]
eff_smooth_DM0 = df_DM0['efficiency_smooth_{}'.format(threshold)]
plt.plot(df_DM0.gentau_vis_pt, eff_smooth_DM0, label='1 prong', ls='-', linewidth=2, color='red')
eff_DM1 = df_DM1['efficiency_{}'.format(threshold)]
eff_smooth_DM1 = df_DM1['efficiency_smooth_{}'.format(threshold)]
plt.plot(df_DM1.gentau_vis_pt, eff_smooth_DM1, label='1 prong + $\pi^{0}$\'s', ls='-', linewidth=2, color='blue')
eff_DM4 = df_DM4['efficiency_{}'.format(threshold)]
eff_smooth_DM4 = df_DM4['efficiency_smooth_{}'.format(threshold)]
plt.plot(df_DM4.gentau_vis_pt, eff_smooth_DM4, label='3 prongs', ls='-', linewidth=2, color='green')
eff_DM5 = df_DM5['efficiency_{}'.format(threshold)]
eff_smooth_DM5 = df_DM5['efficiency_smooth_{}'.format(threshold)]
plt.plot(df_DM5.gentau_vis_pt, eff_smooth_DM5, label='3 prongs + $\pi^{0}$\'s', ls='-', linewidth=2, color='fuchsia')
plt.ylim(0., 1.01)
plt.xlim(20, 90)
plt.legend(loc = 'lower right', fontsize=16)
plt.xlabel(r'$p_{T}^{gen,\tau}\ [GeV]$')
plt.ylabel('Efficiency')
plt.grid()
plt.ylim(0.,1.14)
plt.subplots_adjust(bottom=0.12)
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("200 PU","HGCal calorimeter-only", "$E_{T}^{L1}$ > 50 GeV")
plt.savefig(plotdir+'eff_vs_pt_DM_TDR2.png')
plt.savefig(plotdir+'eff_vs_pt_DM_TDR2.pdf')
'''

import matplotlib.lines as mlines
plt.figure(figsize=(8,7))

df_DM01 = efficiencies_vs_pt_DM01
df_DM45 = efficiencies_vs_pt_DM45

lab_50 = r"$E_{T}^{L1,\tau}$ > 50 GeV"

eff_DM01 = df_DM01['efficiency_50']
eff_DM45 = df_DM45['efficiency_50']

eff_smooth_DM01 = df_DM01['efficiency_smooth_50']
eff_smooth_DM45 = df_DM45['efficiency_smooth_50']

x_DM01 = efficiencies_vs_pt_DM01.gentau_vis_pt
x_DM45 = efficiencies_vs_pt_DM45.gentau_vis_pt

plt.errorbar(x_DM01,eff_DM01,xerr=1,ls='None',label='1-prong',color='red',lw=2,marker='o',mec='red')
plt.errorbar(x_DM45,eff_DM45,xerr=1,ls='None',label='1-prong + $\pi^{0}$\'s',color='blue',lw=2,marker='o',mec='blue')

plt.plot(x_DM01,eff_smooth_DM01,color='red',lw=1.5)
plt.plot(x_DM45,eff_smooth_DM45,color='blue',lw=1.5)

plt.xlim(20, 91.5)
plt.ylim(0., 1.10)

red_line = mlines.Line2D([], [], color='red',markersize=15, label='1-prong (+ $\pi^{0}$\'s)',lw=2)
blue_line = mlines.Line2D([], [], color='blue',markersize=15, label='3-prong (+ $\pi^{0}$\'s)',lw=2)

plt.legend(loc = 'lower right', fontsize=16, handles=[red_line,blue_line])
#plt.legend(loc = 'lower right', fontsize=16)
txt = (r'Gen. $\tau$ decay mode:')
t = plt.text(66,0.21, txt, ha='left', wrap=True, fontsize=16)
t.set_bbox(dict(facecolor='white', edgecolor='white'))
txt2 = (r'$E_{T}^{L1,\tau}$ > 50 GeV')
t2 = plt.text(28,0.83, txt2, ha='left', wrap=True, fontsize=16)
t2.set_bbox(dict(facecolor='white', edgecolor='white'))
txt3 = (r'1.6 < | $\eta_{gen,\tau}$ | < 2.9')
t3 = plt.text(23,0.72, txt3, ha='left', wrap=True, fontsize=16)
t3.set_bbox(dict(facecolor='white', edgecolor='white'))
plt.xlabel(r'$p_{T}^{gen,\tau}\ [GeV]$')
plt.ylabel('Efficiency')
plt.grid()
plt.subplots_adjust(bottom=0.12)
plt.cmstext("CMS"," Phase-2 Simulation")
plt.lumitext("200 PU","HGCal calorimeter-only")
plt.savefig(plotdir+'eff_vs_pt_DM_TDR4.png')
plt.savefig(plotdir+'eff_vs_pt_DM_TDR4.pdf')

