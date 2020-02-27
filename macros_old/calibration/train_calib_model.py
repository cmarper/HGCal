# C. Martin Perez cmartinp@cern.ch, Aug. 2019
# Training of regression BDT for L1 taus HGCal calibration

import sys

import numpy
from numpy import loadtxt, shape

from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.externals import joblib

# Load data
features = ['pt_tot', 'eta_Eweighted','firstlayer', 'maxlayer', 'layer10', 'layer50', 'layer90', 'showerlength', 'hoe', 'meanz']
n_features = len(features)

datasetname = '../../data/calibration/input_train_calib_10vars_etaPlus_PUcut_PU200.csv'
dataset = loadtxt(datasetname, delimiter=",")

print 'Loading data, with',n_features,'features:\n',features
print ' '

# Split data into X (features) and Y (target)
X = dataset[:,0:n_features] # features
Y = dataset[:,n_features] # target: pT_gentau_vis/pT_tot_cl3d

# Split data into train and test sets
seed = 99
test_size = 0.0
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=test_size, random_state=seed)

# Train the random forest classifier
gbrt=GradientBoostingRegressor(alpha=0.9, criterion='friedman_mse', init=None,
             learning_rate=0.1, loss='huber', max_depth=6, max_features=1.0,
             max_leaf_nodes=None, min_impurity_decrease=0.0,
             min_impurity_split=None, min_samples_leaf=3,
             min_samples_split=2, min_weight_fraction_leaf=0.0,
             n_estimators=1000, presort='auto', random_state=None,
             subsample=1.0, verbose=0, warm_start=False)

print X_train
print numpy.shape(X_train)

# Fit the model
gbrt.fit(X_train, y_train)

# Save the model
filename = '../../data/calibration/model_train_calib_10vars_etaPlus_PUcut_PU200.sav'
joblib.dump(gbrt, filename)
