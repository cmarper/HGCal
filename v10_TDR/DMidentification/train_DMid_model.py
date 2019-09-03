# C. Martin Perez cmartinp@cern.ch, Aug. 2019
# Training of multiclassifier BDT for L1 taus HGCal DM identification

import sys

import numpy
from numpy import loadtxt, shape

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib

# Load data
features = ['firstlayer', 'maxlayer', 'layer10', 'layer50', 'layer90', 'showerlength', 'hoe', 'meanz']
n_features = len(features)

datasetname = '../../data/DMid/input_train_DMid_8vars_etaPlus_PUcut.csv'
dataset = loadtxt(datasetname, delimiter=",")

print 'Loading data, with',n_features,'features:\n',features
print ' '

# Split data into X (features) and Y (target)
X = dataset[:,0:n_features] # features
Y = dataset[:,n_features] # target: DM (0,1,4,5)

# Split data into train and test sets
seed = 99
test_size = 0.0
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=test_size, random_state=seed)

# Train the random forest classifier
clf = RandomForestClassifier(n_jobs=2, random_state=0, class_weight='balanced',n_estimators=3000)

# Fit the model
clf.fit(X_train, y_train)

# Save the model
filename = '../../data/DMid/model_train_DMid_8vars_etaPlus_PUcut.sav'
joblib.dump(clf, filename)
