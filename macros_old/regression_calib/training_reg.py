# C. Martin Perez cmartinp@cern.ch
# May. 2019
# Training of regression BDT for L1 taus HGCal calibration

import sys

# Training version, to run:
# python training.py v1
#version = sys.argv[1]
#print version

import xgboost as xgb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from numpy import loadtxt, histogram
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
#from sklearn.cross_validation import ShuffleSplit, train_test_split
#from sklearn.grid_search import GridSearchCV


# Load data
#v1
#features = ['eta_cl1', 'pTfraction_cl1', 'pTfraction_cl2', 'pTfraction_cl3', 'BDTeg_cl1', 'BDTeg_cl2', 'BDTeg_cl3', 'showerlength_cl1', 'showerlength_cl2', 'showerlength_cl3', 'maxlayer_cl1', 'maxlayer_cl2', 'maxlayer_cl3']
#v2
features = ['pT_tot_cl','eta_cl1', 'pTfraction_cl1', 'pTfraction_cl2', 'pTfraction_cl3', 'BDTeg_cl1', 'BDTeg_cl2', 'BDTeg_cl3', 'showerlength_cl1', 'showerlength_cl2', 'showerlength_cl3', 'maxlayer_cl1', 'maxlayer_cl2', 'maxlayer_cl3']
#v3
#features = ['pT_tot_cl','eta_cl1', 'BDTeg_cl1', 'BDTeg_cl2', 'BDTeg_cl3', 'showerlength_cl1', 'showerlength_cl2', 'showerlength_cl3', 'maxlayer_cl1', 'maxlayer_cl2', 'maxlayer_cl3']
#v4
#features = ['pT_tot_cl','eta_cl1', 'pTfraction_cl1', 'BDTeg_cl1', 'showerlength_cl1', 'maxlayer_cl1']


n_features = len(features)

#datasetname = './data/BDTv1_calib_vbles.csv'
datasetname = './data/BDTv2_calib_vbles.csv'
#datasetname = './data/BDTv3_calib_vbles.csv'
#datasetname = './data/BDTv4_calib_vbles.csv'
dataset = loadtxt(datasetname, delimiter=",")

#print ' '
print 'Loading data, with',n_features,'features:\n',features
print ' '

# Split data into X and Y
X = dataset[:,0:n_features] # features
Y = dataset[:,n_features] # target: pT_gentau_vis/pT_tot_cl3d

pT_genvistau = dataset[:,n_features+1]
pT_totcl3d = dataset[:,n_features+2]
eta_gentau = dataset[:,n_features+3]
DM_gentau = dataset[:,n_features+4]

# Split data into train and test sets
seed = 99
test_size = 0.3
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=test_size, random_state=seed)
pT_genvistau_train, pT_genvistau_test, pT_totcl3d_train, pT_totcl3d_test = train_test_split(pT_genvistau, pT_totcl3d, test_size=test_size, random_state=seed)
eta_gentau_train, eta_gentau_test, DM_gentau_train, DM_gentau_test = train_test_split(eta_gentau, DM_gentau, test_size=test_size, random_state=seed)

'''def GradientBooster(param_grid, n_jobs):
  estimator = GradientBoostingRegressor()
  cv = ShuffleSplit(X_train.shape[0], n_iter=10, test_size=0.2)
  classifier = GridSearchCV(estimator=estimator, cv=cv, param_grid=param_grid, n_jobs=n_jobs)
  classifier.fit(X_train, y_train)

  print "Best Estimator learned through GridSearch" 
  print classifier.best_estimator_ 

  return cv, classifier.best_estimator_


param_grid={'n_estimators':[100], 
            'learning_rate': [0.1],# 0.05, 0.02, 0.01],
            'max_depth':[6],#4,6], 
            'min_samples_leaf':[3],#,5,9,17], 
            'max_features':[1.0],#,0.3]#,0.1]
            }

n_jobs=4
cv,best_est=GradientBooster(param_grid, n_jobs)
'''

# Train the random forest classifier
#gbrt=GradientBoostingRegressor(loss='huber',n_estimators=5000)
gbrt=GradientBoostingRegressor(alpha=0.9, criterion='friedman_mse', init=None,
             learning_rate=0.1, loss='huber', max_depth=6, max_features=1.0,
             max_leaf_nodes=None, min_impurity_decrease=0.0,
             min_impurity_split=None, min_samples_leaf=3,
             min_samples_split=2, min_weight_fraction_leaf=0.0,
             n_estimators=1000, presort='auto', random_state=None,
             subsample=1.0, verbose=0, warm_start=False)

gbrt.fit(X_train, y_train)

# Get the predictions
preds = gbrt.predict(X_test)

# Correct for divergences
preds_corr = []
for i in xrange(len(preds)):
  if preds[i] < -20 or preds[i] > 20:
    preds_corr.append(1)
  else:
    preds_corr.append(preds[i])

#print preds[0:100]
#print preds_corr[0:100]

# Get the feature importances and sort them
feat_importances = gbrt.feature_importances_;
indices = np.argsort(feat_importances)[::-1]
ordered_feats = np.empty(n_features, dtype=object)

#print "Feature importance"
for f in xrange(n_features):
    ordered_feats[f] = features[indices[f]]

index = np.empty(n_features, dtype=object)
for i in xrange(n_features):
    index[i] = i

# Plot the feature importances
plt.figure()
plt.bar(index, feat_importances[indices])
plt.xticks(index, ordered_feats)
plt.xticks(rotation=45)
plt.ylabel('Importances')
plt.xlabel('Features',fontsize=12)
plt.title(' ')
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.3)
nameFI = "plots/feature_importance.png"
plt.savefig(nameFI)

#Let's print the R-squared value for train/test. This explains how much of the variance in the data our model is #able to decipher.
#print "R-squared for Train: %.2f" %gbrt.score(X_train, y_train)
#print "R-squared for Test: %.2f" %gbrt.score(X_test, y_test)

from matplotlib.colors import LogNorm
# Plot calibration factors
bins = np.linspace(0, 10, 50)
plt.figure(figsize=(10, 10))
plt.hist2d(y_test,preds,bins,norm=LogNorm())
plt.title(" ")
plt.xlabel("True calibration factor")
plt.ylabel("Predicted calibration factor")
plt.plot((0, 0), (10, 10))
namecalib = "./plots/calib_factors.png"
plt.savefig(namecalib)
#plt.show()

# Plot responses
uncalib_pt_test = pT_totcl3d_test
calib_pt_test = preds_corr*pT_totcl3d_test

#print 'uncalib ',uncalib_pt_test[0:10]
#print 'preds ', preds[0:10]
#print 'calib ',calib_pt_test[0:10]

res_uncalib_test = uncalib_pt_test/pT_genvistau_test
res_calib_test = calib_pt_test/pT_genvistau_test

mean_res_uncalib_test = np.mean(res_uncalib_test)
rms_res_uncalib_test = np.std(res_uncalib_test)
rms_o_mean_res_uncalib_test = rms_res_uncalib_test/mean_res_uncalib_test

mean_res_calib_test = np.mean(res_calib_test)
rms_res_calib_test = np.std(res_calib_test)
rms_o_mean_res_calib_test = rms_res_calib_test/mean_res_calib_test

res_uncalib_center_test = res_uncalib_test/mean_res_uncalib_test
res_calib_center_test = res_calib_test/mean_res_calib_test

print 'Inclusive'
print 'Uncalib: mean ',mean_res_uncalib_test,'rms ',rms_res_uncalib_test,'rms/mean ',rms_o_mean_res_uncalib_test
print 'Calib: mean ',mean_res_calib_test,'rms ',rms_res_calib_test,'rms/mean ',rms_o_mean_res_calib_test

bins = np.linspace(0, 3.5, 50)
plt.figure()
plt.hist([res_uncalib_test,res_calib_test],bins,histtype='step',linewidth=1.5,fill=False, stacked=False,color=['red','limegreen'], label=['Uncalibrated','Calibrated'])
plt.axvline(res_uncalib_test.mean(), color='red', linestyle='dashed', linewidth=1.2)
plt.axvline(res_calib_test.mean(), color='limegreen', linestyle='dashed', linewidth=1.2)
plt.legend(loc='upper right')
plt.title(" ")
plt.xlabel("p$_T$ (matched cl3d) / p$_T$ (gen. tau vis.)")
plt.ylabel("Entries")
namecalib2 = "./plots/response.png"
plt.savefig(namecalib2)
#plt.show()

bins = np.linspace(0, 4, 50)
plt.figure()
plt.hist([res_uncalib_center_test,res_calib_center_test],bins,histtype='step',linewidth=1.5,fill=False, stacked=False,color=['red','limegreen'], label=['Uncalibrated','Calibrated'])
plt.legend(loc='upper right')
plt.title(" ")
plt.xlabel("[p$_T$ (matched cl3d) / p$_T$ (gen. tau vis.)] / mean")
plt.ylabel("Entries")
namecalib3 = "./plots/response_centered.png"
plt.savefig(namecalib3)
#plt.show()


# Response in pT bins

res_uncalib_test_pT0to25 = []
res_uncalib_test_pT25to50 = []
res_uncalib_test_pT50to75 = []
res_uncalib_test_pT75 = []

res_calib_test_pT0to25 = []
res_calib_test_pT25to50 = []
res_calib_test_pT50to75 = []
res_calib_test_pT75 = []

for i in xrange(len(res_uncalib_test)):
  if (pT_genvistau_test[i] <= 25):
     res_uncalib_test_pT0to25.append(res_uncalib_test[i])
     res_calib_test_pT0to25.append(res_calib_test[i])
  elif (pT_genvistau_test[i] > 25 and pT_genvistau_test[i] <= 50):
     res_uncalib_test_pT25to50.append(res_uncalib_test[i])
     res_calib_test_pT25to50.append(res_calib_test[i])
  elif (pT_genvistau_test[i] > 50 and pT_genvistau_test[i] <= 75):
     res_uncalib_test_pT50to75.append(res_uncalib_test[i])
     res_calib_test_pT50to75.append(res_calib_test[i])
  elif (pT_genvistau_test[i] > 75):
     res_uncalib_test_pT75.append(res_uncalib_test[i])
     res_calib_test_pT75.append(res_calib_test[i])


mean_res_uncalib_test_pT0to25 = np.mean(res_uncalib_test_pT0to25)
rms_res_uncalib_test_pT0to25 = np.std(res_uncalib_test_pT0to25)
rms_o_mean_res_uncalib_test_pT0to25 = rms_res_uncalib_test_pT0to25/mean_res_uncalib_test_pT0to25
mean_res_calib_test_pT0to25 = np.mean(res_calib_test_pT0to25)
rms_res_calib_test_pT0to25 = np.std(res_calib_test_pT0to25)
rms_o_mean_res_calib_test_pT0to25 = rms_res_calib_test_pT0to25/mean_res_calib_test_pT0to25
print ' '
print 'pT0to25'
print 'Uncalib: mean ',mean_res_uncalib_test_pT0to25,'rms ',rms_res_uncalib_test_pT0to25,'rms/mean ',rms_o_mean_res_uncalib_test_pT0to25
print 'Calib: mean ',mean_res_calib_test_pT0to25,'rms ',rms_res_calib_test_pT0to25,'rms/mean ',rms_o_mean_res_calib_test_pT0to25

mean_res_uncalib_test_pT25to50 = np.mean(res_uncalib_test_pT25to50)
rms_res_uncalib_test_pT25to50 = np.std(res_uncalib_test_pT25to50)
rms_o_mean_res_uncalib_test_pT25to50 = rms_res_uncalib_test_pT25to50/mean_res_uncalib_test_pT25to50
mean_res_calib_test_pT25to50 = np.mean(res_calib_test_pT25to50)
rms_res_calib_test_pT25to50 = np.std(res_calib_test_pT25to50)
rms_o_mean_res_calib_test_pT25to50 = rms_res_calib_test_pT25to50/mean_res_calib_test_pT25to50
print ' '
print 'pT25to50'
print 'Uncalib: mean ',mean_res_uncalib_test_pT25to50,'rms ',rms_res_uncalib_test_pT25to50,'rms/mean ',rms_o_mean_res_uncalib_test_pT25to50
print 'Calib: mean ',mean_res_calib_test_pT25to50,'rms ',rms_res_calib_test_pT25to50,'rms/mean ',rms_o_mean_res_calib_test_pT25to50

mean_res_uncalib_test_pT50to75 = np.mean(res_uncalib_test_pT50to75)
rms_res_uncalib_test_pT50to75 = np.std(res_uncalib_test_pT50to75)
rms_o_mean_res_uncalib_test_pT50to75 = rms_res_uncalib_test_pT50to75/mean_res_uncalib_test_pT50to75
mean_res_calib_test_pT50to75 = np.mean(res_calib_test_pT50to75)
rms_res_calib_test_pT50to75 = np.std(res_calib_test_pT50to75)
rms_o_mean_res_calib_test_pT50to75 = rms_res_calib_test_pT50to75/mean_res_calib_test_pT50to75
print ' '
print 'pT50to75'
print 'Uncalib: mean ',mean_res_uncalib_test_pT50to75,'rms ',rms_res_uncalib_test_pT50to75,'rms/mean ',rms_o_mean_res_uncalib_test_pT50to75
print 'Calib: mean ',mean_res_calib_test_pT50to75,'rms ',rms_res_calib_test_pT50to75,'rms/mean ',rms_o_mean_res_calib_test_pT50to75

mean_res_uncalib_test_pT75 = np.mean(res_uncalib_test_pT75)
rms_res_uncalib_test_pT75 = np.std(res_uncalib_test_pT75)
rms_o_mean_res_uncalib_test_pT75 = rms_res_uncalib_test_pT75/mean_res_uncalib_test_pT75
mean_res_calib_test_pT75 = np.mean(res_calib_test_pT75)
rms_res_calib_test_pT75 = np.std(res_calib_test_pT75)
rms_o_mean_res_calib_test_pT75 = rms_res_calib_test_pT75/mean_res_calib_test_pT75
print ' '
print 'pT75'
print 'Uncalib: mean ',mean_res_uncalib_test_pT75,'rms ',rms_res_uncalib_test_pT75,'rms/mean ',rms_o_mean_res_uncalib_test_pT75
print 'Calib: mean ',mean_res_calib_test_pT75,'rms ',rms_res_calib_test_pT75,'rms/mean ',rms_o_mean_res_calib_test_pT75


bins = np.linspace(0, 3.5, 35)
fig = plt.figure(figsize=(10, 10))

plt.subplot(2, 2, 1)
plt.hist([res_uncalib_test_pT0to25,res_calib_test_pT0to25],bins,histtype='step',linewidth=1.5,fill=False, stacked=False,color=['red','limegreen'], label=['Uncalibrated','Calibrated'])
plt.axvline(mean_res_uncalib_test_pT0to25, color='red', linestyle='dashed', linewidth=1.2)
plt.axvline(mean_res_calib_test_pT0to25, color='limegreen', linestyle='dashed', linewidth=1.2)
plt.legend(loc='upper right')
plt.title("p$_T$ < 25 GeV")
#plt.xlabel("p$_T$ (matched cl3d) / p$_T$ (gen. tau vis.)")
plt.ylabel("Entries")

plt.subplot(2, 2, 2)
plt.hist([res_uncalib_test_pT25to50,res_calib_test_pT25to50],bins,histtype='step',linewidth=1.5,fill=False, stacked=False,color=['red','limegreen'], label=['Uncalibrated','Calibrated'])
plt.axvline(mean_res_uncalib_test_pT25to50, color='red', linestyle='dashed', linewidth=1.2)
plt.axvline(mean_res_calib_test_pT25to50, color='limegreen', linestyle='dashed', linewidth=1.2)
plt.legend(loc='upper right')
plt.title("25 GeV < p$_T$ < 50 GeV")
#plt.xlabel("p$_T$ (matched cl3d) / p$_T$ (gen. tau vis.)")
#plt.ylabel("Entries")

plt.subplot(2, 2, 3)
plt.hist([res_uncalib_test_pT50to75,res_calib_test_pT50to75],bins,histtype='step',linewidth=1.5,fill=False, stacked=False,color=['red','limegreen'], label=['Uncalibrated','Calibrated'])
plt.axvline(mean_res_uncalib_test_pT50to75, color='red', linestyle='dashed', linewidth=1.2)
plt.axvline(mean_res_calib_test_pT50to75, color='limegreen', linestyle='dashed', linewidth=1.2)
plt.legend(loc='upper right')
plt.title("50 GeV < p$_T$ < 75 GeV")
plt.xlabel("p$_T$ (matched cl3d) / p$_T$ (gen. tau vis.)")
plt.ylabel("Entries")

plt.subplot(2, 2, 4)
plt.hist([res_uncalib_test_pT75,res_calib_test_pT75],bins,histtype='step',linewidth=1.5,fill=False, stacked=False,color=['red','limegreen'], label=['Uncalibrated','Calibrated'])
plt.axvline(mean_res_uncalib_test_pT75, color='red', linestyle='dashed', linewidth=1.2)
plt.axvline(mean_res_calib_test_pT75, color='limegreen', linestyle='dashed', linewidth=1.2)
plt.legend(loc='upper right')
plt.title("75 GeV < p$_T$")
plt.xlabel("p$_T$ (matched cl3d) / p$_T$ (gen. tau vis.)")
#plt.ylabel("Entries")

namecalib4 = "./plots/response_ptbins.png"
plt.savefig(namecalib4)
#plt.show()

# Response in eta bins

res_uncalib_test_eta1p7to2p1 = []
res_uncalib_test_eta2p1to2p4 = []
res_uncalib_test_eta2p4to2p8 = []

res_calib_test_eta1p7to2p1 = []
res_calib_test_eta2p1to2p4 = []
res_calib_test_eta2p4to2p8 = []

for i in xrange(len(res_uncalib_test)):
  if (abs(eta_gentau_test[i]) >= 1.7 and abs(eta_gentau_test[i]) < 2.1):
     res_uncalib_test_eta1p7to2p1.append(res_uncalib_test[i])
     res_calib_test_eta1p7to2p1.append(res_calib_test[i])
  elif (abs(eta_gentau_test[i]) >= 2.1 and abs(eta_gentau_test[i]) < 2.4):
     res_uncalib_test_eta2p1to2p4.append(res_uncalib_test[i])
     res_calib_test_eta2p1to2p4.append(res_calib_test[i])
  elif (abs(eta_gentau_test[i]) >= 2.4 and abs(eta_gentau_test[i]) < 2.8):
     res_uncalib_test_eta2p4to2p8.append(res_uncalib_test[i])
     res_calib_test_eta2p4to2p8.append(res_calib_test[i])


mean_res_uncalib_test_eta1p7to2p1 = np.mean(res_uncalib_test_eta1p7to2p1)
rms_res_uncalib_test_eta1p7to2p1 = np.std(res_uncalib_test_eta1p7to2p1)
rms_o_mean_res_uncalib_test_eta1p7to2p1 = rms_res_uncalib_test_eta1p7to2p1/mean_res_uncalib_test_eta1p7to2p1
mean_res_calib_test_eta1p7to2p1 = np.mean(res_calib_test_eta1p7to2p1)
rms_res_calib_test_eta1p7to2p1 = np.std(res_calib_test_eta1p7to2p1)
rms_o_mean_res_calib_test_eta1p7to2p1 = rms_res_calib_test_eta1p7to2p1/mean_res_calib_test_eta1p7to2p1
print ' '
print 'eta1p7to2p1'
print 'Uncalib: mean ',mean_res_uncalib_test_eta1p7to2p1,'rms ',rms_res_uncalib_test_eta1p7to2p1,'rms/mean ',rms_o_mean_res_uncalib_test_eta1p7to2p1
print 'Calib: mean ',mean_res_calib_test_eta1p7to2p1,'rms ',rms_res_calib_test_eta1p7to2p1,'rms/mean ',rms_o_mean_res_calib_test_eta1p7to2p1

mean_res_uncalib_test_eta2p1to2p4 = np.mean(res_uncalib_test_eta2p1to2p4)
rms_res_uncalib_test_eta2p1to2p4 = np.std(res_uncalib_test_eta2p1to2p4)
rms_o_mean_res_uncalib_test_eta2p1to2p4 = rms_res_uncalib_test_eta2p1to2p4/mean_res_uncalib_test_eta2p1to2p4
mean_res_calib_test_eta2p1to2p4 = np.mean(res_calib_test_eta2p1to2p4)
rms_res_calib_test_eta2p1to2p4 = np.std(res_calib_test_eta2p1to2p4)
rms_o_mean_res_calib_test_eta2p1to2p4 = rms_res_calib_test_eta2p1to2p4/mean_res_calib_test_eta2p1to2p4
print ' '
print 'eta2p1to2p4'
print 'Uncalib: mean ',mean_res_uncalib_test_eta2p1to2p4,'rms ',rms_res_uncalib_test_eta2p1to2p4,'rms/mean ',rms_o_mean_res_uncalib_test_eta2p1to2p4
print 'Calib: mean ',mean_res_calib_test_eta2p1to2p4,'rms ',rms_res_calib_test_eta2p1to2p4,'rms/mean ',rms_o_mean_res_calib_test_eta2p1to2p4

mean_res_uncalib_test_eta2p4to2p8 = np.mean(res_uncalib_test_eta2p4to2p8)
rms_res_uncalib_test_eta2p4to2p8 = np.std(res_uncalib_test_eta2p4to2p8)
rms_o_mean_res_uncalib_test_eta2p4to2p8 = rms_res_uncalib_test_eta2p4to2p8/mean_res_uncalib_test_eta2p4to2p8
mean_res_calib_test_eta2p4to2p8 = np.mean(res_calib_test_eta2p4to2p8)
rms_res_calib_test_eta2p4to2p8 = np.std(res_calib_test_eta2p4to2p8)
rms_o_mean_res_calib_test_eta2p4to2p8 = rms_res_calib_test_eta2p4to2p8/mean_res_calib_test_eta2p4to2p8
print ' '
print 'eta2p4to2p8'
print 'Uncalib: mean ',mean_res_uncalib_test_eta2p4to2p8,'rms ',rms_res_uncalib_test_eta2p4to2p8,'rms/mean ',rms_o_mean_res_uncalib_test_eta2p4to2p8
print 'Calib: mean ',mean_res_calib_test_eta2p4to2p8,'rms ',rms_res_calib_test_eta2p4to2p8,'rms/mean ',rms_o_mean_res_calib_test_eta2p4to2p8

bins = np.linspace(0, 3.5, 35)
fig = plt.figure(figsize=(15, 5))

plt.subplot(1, 3, 1)
plt.hist([res_uncalib_test_eta1p7to2p1,res_calib_test_eta1p7to2p1],bins,histtype='step',linewidth=1.5,fill=False, stacked=False,color=['red','limegreen'], label=['Uncalibrated','Calibrated'])
plt.axvline(mean_res_uncalib_test_eta1p7to2p1, color='red', linestyle='dashed', linewidth=1.2)
plt.axvline(mean_res_calib_test_eta1p7to2p1, color='limegreen', linestyle='dashed', linewidth=1.2)
plt.legend(loc='upper right')
plt.title("1.7 < abs(eta) < 2.1")
plt.xlabel("p$_T$ (matched cl3d) / p$_T$ (gen. tau vis.)")
plt.ylabel("Entries")

plt.subplot(1, 3, 2)
plt.hist([res_uncalib_test_eta2p1to2p4,res_calib_test_eta2p1to2p4],bins,histtype='step',linewidth=1.5,fill=False, stacked=False,color=['red','limegreen'], label=['Uncalibrated','Calibrated'])
plt.axvline(mean_res_uncalib_test_eta2p1to2p4, color='red', linestyle='dashed', linewidth=1.2)
plt.axvline(mean_res_calib_test_eta2p1to2p4, color='limegreen', linestyle='dashed', linewidth=1.2)
plt.legend(loc='upper right')
plt.title("2.1 < abs(eta) < 2.4")
plt.xlabel("p$_T$ (matched cl3d) / p$_T$ (gen. tau vis.)")

plt.subplot(1, 3, 3)
plt.hist([res_uncalib_test_eta2p4to2p8,res_calib_test_eta2p4to2p8],bins,histtype='step',linewidth=1.5,fill=False, stacked=False,color=['red','limegreen'], label=['Uncalibrated','Calibrated'])
plt.axvline(mean_res_uncalib_test_eta2p4to2p8, color='red', linestyle='dashed', linewidth=1.2)
plt.axvline(mean_res_calib_test_eta2p4to2p8, color='limegreen', linestyle='dashed', linewidth=1.2)
plt.legend(loc='upper right')
plt.title("2.4 < abs(eta) < 2.8")
plt.xlabel("p$_T$ (matched cl3d) / p$_T$ (gen. tau vis.)")

namecalib5 = "./plots/response_etabins.png"
plt.savefig(namecalib5)
#plt.show()

# Response per DM

res_uncalib_test_DM0 = []
res_uncalib_test_DM1 = []
res_uncalib_test_DM4 = []
res_uncalib_test_DM5 = []

res_calib_test_DM0 = []
res_calib_test_DM1 = []
res_calib_test_DM4 = []
res_calib_test_DM5 = []

for i in xrange(len(res_uncalib_test)):
  if (DM_gentau_test[i] == 0):
     res_uncalib_test_DM0.append(res_uncalib_test[i])
     res_calib_test_DM0.append(res_calib_test[i])
  elif (DM_gentau_test[i] == 1):
     res_uncalib_test_DM1.append(res_uncalib_test[i])
     res_calib_test_DM1.append(res_calib_test[i])
  elif (DM_gentau_test[i] == 4):
     res_uncalib_test_DM4.append(res_uncalib_test[i])
     res_calib_test_DM4.append(res_calib_test[i])
  elif (DM_gentau_test[i] == 5):
     res_uncalib_test_DM5.append(res_uncalib_test[i])
     res_calib_test_DM5.append(res_calib_test[i])

mean_res_uncalib_test_DM0 = np.mean(res_uncalib_test_DM0)
rms_res_uncalib_test_DM0 = np.std(res_uncalib_test_DM0)
rms_o_mean_res_uncalib_test_DM0 = rms_res_uncalib_test_DM0/mean_res_uncalib_test_DM0
mean_res_calib_test_DM0 = np.mean(res_calib_test_DM0)
rms_res_calib_test_DM0 = np.std(res_calib_test_DM0)
rms_o_mean_res_calib_test_DM0 = rms_res_calib_test_DM0/mean_res_calib_test_DM0
print ' '
print 'DM0'
print 'Uncalib: mean ',mean_res_uncalib_test_DM0,'rms ',rms_res_uncalib_test_DM0,'rms/mean ',rms_o_mean_res_uncalib_test_DM0
print 'Calib: mean ',mean_res_calib_test_DM0,'rms ',rms_res_calib_test_DM0,'rms/mean ',rms_o_mean_res_calib_test_DM0

mean_res_uncalib_test_DM1 = np.mean(res_uncalib_test_DM1)
rms_res_uncalib_test_DM1 = np.std(res_uncalib_test_DM1)
rms_o_mean_res_uncalib_test_DM1 = rms_res_uncalib_test_DM1/mean_res_uncalib_test_DM1
mean_res_calib_test_DM1 = np.mean(res_calib_test_DM1)
rms_res_calib_test_DM1 = np.std(res_calib_test_DM1)
rms_o_mean_res_calib_test_DM1 = rms_res_calib_test_DM1/mean_res_calib_test_DM1
print ' '
print 'DM1'
print 'Uncalib: mean ',mean_res_uncalib_test_DM1,'rms ',rms_res_uncalib_test_DM1,'rms/mean ',rms_o_mean_res_uncalib_test_DM1
print 'Calib: mean ',mean_res_calib_test_DM1,'rms ',rms_res_calib_test_DM1,'rms/mean ',rms_o_mean_res_calib_test_DM1

mean_res_uncalib_test_DM4 = np.mean(res_uncalib_test_DM4)
rms_res_uncalib_test_DM4 = np.std(res_uncalib_test_DM4)
rms_o_mean_res_uncalib_test_DM4 = rms_res_uncalib_test_DM4/mean_res_uncalib_test_DM4
mean_res_calib_test_DM4 = np.mean(res_calib_test_DM4)
rms_res_calib_test_DM4 = np.std(res_calib_test_DM4)
rms_o_mean_res_calib_test_DM4 = rms_res_calib_test_DM4/mean_res_calib_test_DM4
print ' '
print 'DM4'
print 'Uncalib: mean ',mean_res_uncalib_test_DM4,'rms ',rms_res_uncalib_test_DM4,'rms/mean ',rms_o_mean_res_uncalib_test_DM4
print 'Calib: mean ',mean_res_calib_test_DM4,'rms ',rms_res_calib_test_DM4,'rms/mean ',rms_o_mean_res_calib_test_DM4

mean_res_uncalib_test_DM5 = np.mean(res_uncalib_test_DM5)
rms_res_uncalib_test_DM5 = np.std(res_uncalib_test_DM5)
rms_o_mean_res_uncalib_test_DM5 = rms_res_uncalib_test_DM5/mean_res_uncalib_test_DM5
mean_res_calib_test_DM5 = np.mean(res_calib_test_DM5)
rms_res_calib_test_DM5 = np.std(res_calib_test_DM5)
rms_o_mean_res_calib_test_DM5 = rms_res_calib_test_DM5/mean_res_calib_test_DM5
print ' '
print 'DM5'
print 'Uncalib: mean ',mean_res_uncalib_test_DM5,'rms ',rms_res_uncalib_test_DM5,'rms/mean ',rms_o_mean_res_uncalib_test_DM5
print 'Calib: mean ',mean_res_calib_test_DM5,'rms ',rms_res_calib_test_DM5,'rms/mean ',rms_o_mean_res_calib_test_DM5


bins = np.linspace(0, 3.5, 35)
fig = plt.figure(figsize=(10, 10))

plt.subplot(2, 2, 1)
plt.hist([res_uncalib_test_DM0,res_calib_test_DM0],bins,histtype='step',linewidth=1.5,fill=False, stacked=False,color=['red','limegreen'], label=['Uncalibrated','Calibrated'])
plt.axvline(mean_res_uncalib_test_DM0, color='red', linestyle='dashed', linewidth=1.2)
plt.axvline(mean_res_calib_test_DM0, color='limegreen', linestyle='dashed', linewidth=1.2)
plt.legend(loc='upper right')
plt.title("1prong")
#plt.xlabel("p$_T$ (matched cl3d) / p$_T$ (gen. tau vis.)")
plt.ylabel("Entries")

plt.subplot(2, 2, 2)
plt.hist([res_uncalib_test_DM1,res_calib_test_DM1],bins,histtype='step',linewidth=1.5,fill=False, stacked=False,color=['red','limegreen'], label=['Uncalibrated','Calibrated'])
plt.axvline(mean_res_uncalib_test_DM1, color='red', linestyle='dashed', linewidth=1.2)
plt.axvline(mean_res_calib_test_DM1, color='limegreen', linestyle='dashed', linewidth=1.2)
plt.legend(loc='upper right')
plt.title("1prong+pi0")
#plt.xlabel("p$_T$ (matched cl3d) / p$_T$ (gen. tau vis.)")
#plt.ylabel("Entries")

plt.subplot(2, 2, 3)
plt.hist([res_uncalib_test_DM4,res_calib_test_DM4],bins,histtype='step',linewidth=1.5,fill=False, stacked=False,color=['red','limegreen'], label=['Uncalibrated','Calibrated'])
plt.axvline(mean_res_uncalib_test_DM4, color='red', linestyle='dashed', linewidth=1.2)
plt.axvline(mean_res_calib_test_DM4, color='limegreen', linestyle='dashed', linewidth=1.2)
plt.legend(loc='upper right')
plt.title("3prong")
plt.xlabel("p$_T$ (matched cl3d) / p$_T$ (gen. tau vis.)")
plt.ylabel("Entries")

plt.subplot(2, 2, 4)
plt.hist([res_uncalib_test_DM5,res_calib_test_DM5],bins,histtype='step',linewidth=1.5,fill=False, stacked=False,color=['red','limegreen'], label=['Uncalibrated','Calibrated'])
plt.axvline(mean_res_uncalib_test_DM5, color='red', linestyle='dashed', linewidth=1.2)
plt.axvline(mean_res_calib_test_DM5, color='limegreen', linestyle='dashed', linewidth=1.2)
plt.legend(loc='upper right')
plt.title("3prong+pi0")
plt.xlabel("p$_T$ (matched cl3d) / p$_T$ (gen. tau vis.)")
#plt.ylabel("Entries")

namecalib6 = "./plots/response_DMbins.png"
plt.savefig(namecalib6)
#plt.show()

