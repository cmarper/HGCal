# C. Martin Perez cmartinp@cern.ch, Aug. 2019
# Training of binary BDT for L1 taus vs PU discrimination

execfile('xgboost_to_tmva.py')

import sys

import numpy
from numpy import loadtxt, shape

import xgboost as xgb

from sklearn.model_selection import train_test_split
from sklearn.externals import joblib
from sklearn import metrics
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score

import matplotlib.pyplot as plt

# Load data
features = ['cl3d_abseta','cl3d_showerlength',
       'cl3d_coreshowerlength', 'cl3d_firstlayer', 'cl3d_maxlayer', 'cl3d_szz',
       'cl3d_seetot', 'cl3d_spptot', 'cl3d_srrtot', 'cl3d_srrmean']

n_features = len(features)

# Merge the two csv files with 
# cat *puBDT*.csv > input_train_puBDT_10vars_etaPlus.csv
datasetname = '../../data/input_train_puBDT_10vars_etaPlus.csv'
dataset = loadtxt(datasetname, delimiter=",")

print 'Loading data, with',n_features,'features:\n',features
print ' '

# Split data into X (features) and Y (target)
X = dataset[:,0:n_features] # features
Y = dataset[:,n_features] # target: isTau

# Split data into train and test sets
seed = 99
test_size = 0.3
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=test_size, random_state=seed)

dtrain = xgb.DMatrix(data=X_train,label=y_train, feature_names=features)
dtest  = xgb.DMatrix(data=X_test,label=y_test,feature_names=features)

hyperparam = {
    'nthread' : 10,  # limit number of threads
    'eta' : 0.2, # learning rate
    'max_depth' : 4,  # maximum depth of a tree
    'subsample' : 0.8, # fraction of events to train tree on
    'colsample_bytree' : 0.8, # fraction of features to train tree on
    'silent' : True,
    'objective' : 'binary:logistic' }  

#hyperparam = {
#    'silent' : True,
#    'objective' : 'binary:logistic', # binary, output probability
#    'max_depth': 2,  # the maximum depth of each tree
#    'eta': 0.1 }  # the training step for each iteration

num_round = 100  # the number of training iterations

model = xgb.train(hyperparam, dtrain, num_round)

# Generate predictions
preds = model.predict(dtest)
print 'preds ', preds

# ROC curves
fpr, tpr, _ = roc_curve(y_test, preds, pos_label=1)
#fpr/tpr: false/true positive rate

plt.figure()
plt.plot(tpr, fpr, label='BDT')
plt.title("Signal: Tau / Background: PU")
plt.xlabel("Signal efficiency")
plt.ylabel("Background efficiency")
plt.legend(loc='upper left')
plt.savefig("../../plots/png/PU_BDT/roc_curve.png")
plt.savefig("../../plots/png/PU_BDT/roc_curve.pdf")

# Area under ROC curve

auc = roc_auc_score(y_test, preds)
print 'AUC ',auc

# Importance
xgb.plot_importance(model)
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.1)
plt.savefig("../../plots/png/PU_BDT/ranking.png")
plt.savefig("../../plots/png/PU_BDT/ranking.pdf")

## Output file with xml weights

model_dump = model.get_dump()
convert_model(model_dump,input_variables=[('cl3d_abseta','F'),('cl3d_showerlength','I'),('cl3d_coreshowerlength','I'),
	('cl3d_firstlayer','I'),('cl3d_maxlayer','I'),('cl3d_szz','F'),('cl3d_seetot','F'),('cl3d_spptot','F'),
	('cl3d_srrtot','F'),('cl3d_srrmean','F')],output_xml='../../data/xgboost_weights_puBDT.xml')
