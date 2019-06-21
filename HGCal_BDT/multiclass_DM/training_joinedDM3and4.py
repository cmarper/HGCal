# C. Martin Perez cmartinp@cern.ch
# May. 2019
# Training of multiclassifier BDT for L1 taus HGCal DM identification

import sys

# Training version, to run:
# python training.py v1
version = sys.argv[1]
#print version

import xgboost as xgb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from numpy import loadtxt, histogram
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from sklearn.utils.multiclass import unique_labels

execfile('xgboost_to_tmva.py')

# Load data #v6 and v7 have same features
features_v6 = ['pTfraction_cl1', 'pTfraction_cl2', 'pTfraction_cl3', 'BDTeg_cl1', 'BDTeg_cl2', 'BDTeg_cl3', 'showerlength_cl1', 'showerlength_cl2', 'showerlength_cl3', 'maxlayer_cl1', 'maxlayer_cl2', 'maxlayer_cl3']

if version=='v7':
  features = features_v6
else:
  print "Error: Wrong version given!"
n_features = len(features)

datasetname = './data/BDT'+version+'_vbles.csv'
dataset = loadtxt(datasetname, delimiter=",")

#print ' '
print 'Loading',version,'with',n_features,'features:\n',features
print ' '

# Split data into X and Y
X = dataset[:,0:n_features] # features
Y = dataset[:,n_features] # target: DM (0,1,3,4)

# Split data into train and test sets
seed = 99
test_size = 0.3
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=test_size, random_state=seed)

# Train the random forest classifier
clf = RandomForestClassifier(n_jobs=2, random_state=0, class_weight='balanced')
clf.fit(X_train, y_train)

# Get the predictions
preds = clf.predict(X_test)
#print(preds[0:20])

# Get the predicted probabilities
preds_proba = clf.predict_proba(X_test)
classes = clf.classes_
#print classes

#realDM = []

#probDM0 = []
#probDM1 = []
#probDM4 = []
#probDM5 = []

probDM0_realDM0 = []
probDM0_realDM4 = []
probDM0_realDM6 = []

probDM4_realDM0 = []
probDM4_realDM4 = []
probDM4_realDM6 = []

probDM6_realDM0 = []
probDM6_realDM4 = []
probDM6_realDM6 = []

for i in xrange(len(preds_proba)):

   if (y_test[i] == 0):
      probDM0_realDM0.append(preds_proba[i][0])
      probDM4_realDM0.append(preds_proba[i][1])
      probDM6_realDM0.append(preds_proba[i][2])

   elif (y_test[i] == 4):
      probDM0_realDM4.append(preds_proba[i][0])
      probDM4_realDM4.append(preds_proba[i][1])
      probDM6_realDM4.append(preds_proba[i][2])

   elif (y_test[i] == 6):
      probDM0_realDM6.append(preds_proba[i][0])
      probDM4_realDM6.append(preds_proba[i][1])
      probDM6_realDM6.append(preds_proba[i][2])

# Plot the probs of each DM

weights_probDM0_realDM0 = np.ones_like(probDM0_realDM0)/float(len(probDM0_realDM0))
weights_probDM0_realDM4 = np.ones_like(probDM0_realDM4)/float(len(probDM0_realDM4))
weights_probDM0_realDM6 = np.ones_like(probDM0_realDM6)/float(len(probDM0_realDM6))

weights_probDM0 = [weights_probDM0_realDM0,weights_probDM0_realDM4,weights_probDM0_realDM6]

bins = np.linspace(0, 1, 11)
plt.figure()
plt.hist([probDM0_realDM0,probDM0_realDM4,probDM0_realDM6],bins,stacked=True,color=['red','blue','forestgreen'],weights=weights_probDM0, label=['True 1prong','True 3prong','True 1/3prong+pi0'])
plt.legend(loc='upper right')
plt.title("v7")
plt.xlabel("Predicted probability of decay mode 1prong")
plt.ylabel("Normalized entries")
nameDM0 = "./plots/probDM0_stacked_"+version+".png"
plt.savefig(nameDM0)

plt.figure()
plt.hist([probDM0_realDM0,probDM0_realDM4,probDM0_realDM6],bins, histtype='step',linewidth=1.5,fill=False, stacked=False,color=['red','blue','forestgreen'],weights=weights_probDM0, label=['True 1prong','True 3prong','True 1/3prong+pi0'])
plt.legend(loc='upper right')
plt.title("v7")
plt.xlabel("Predicted probability of decay mode 1prong")
plt.ylabel("Normalized entries")
nameDM0 = "./plots/probDM0_"+version+".png"
plt.savefig(nameDM0)

weights_probDM4_realDM0 = np.ones_like(probDM4_realDM0)/float(len(probDM4_realDM0))
weights_probDM4_realDM4 = np.ones_like(probDM4_realDM4)/float(len(probDM4_realDM4))
weights_probDM4_realDM6 = np.ones_like(probDM4_realDM6)/float(len(probDM4_realDM6))

weights_probDM4 = [weights_probDM4_realDM0,weights_probDM4_realDM4,weights_probDM4_realDM6]

bins = np.linspace(0, 1, 11)
plt.figure()
plt.hist([probDM4_realDM0,probDM4_realDM4,probDM4_realDM6],bins,stacked=True,color=['red','blue','forestgreen'],weights=weights_probDM4, label=['True 1prong','True 3prong','True 1/3prong+pi0'])
plt.legend(loc='upper right')
plt.title("v7")
plt.xlabel("Predicted probability of decay mode 3prong")
plt.ylabel("Normalized entries")
nameDM0 = "./plots/probDM4_stacked_"+version+".png"
plt.savefig(nameDM0)

plt.figure()
plt.hist([probDM4_realDM0,probDM4_realDM4,probDM4_realDM6],bins, histtype='step',linewidth=1.5,fill=False, stacked=False,color=['red','blue','forestgreen'],weights=weights_probDM4, label=['True 1prong','True 3prong','True 1/3prong+pi0'])
plt.legend(loc='upper right')
plt.title("v7")
plt.xlabel("Predicted probability of decay mode 3prong")
plt.ylabel("Normalized entries")
nameDM0 = "./plots/probDM4_"+version+".png"
plt.savefig(nameDM0)

weights_probDM6_realDM0 = np.ones_like(probDM6_realDM0)/float(len(probDM6_realDM0))
weights_probDM6_realDM4 = np.ones_like(probDM6_realDM4)/float(len(probDM6_realDM4))
weights_probDM6_realDM6 = np.ones_like(probDM6_realDM6)/float(len(probDM6_realDM6))

weights_probDM6 = [weights_probDM6_realDM0,weights_probDM6_realDM4,weights_probDM6_realDM6]

bins = np.linspace(0, 1, 11)
plt.figure()
plt.hist([probDM6_realDM0,probDM6_realDM4,probDM6_realDM6],bins,stacked=True,color=['red','blue','forestgreen'],weights=weights_probDM6, label=['True 1prong','True 3prong','True 1/3prong+pi0'])
plt.legend(loc='upper right')
plt.title("v7")
plt.xlabel("Predicted probability of decay mode 1/3prong+pi0")
plt.ylabel("Normalized entries")
nameDM0 = "./plots/probDM6_stacked_"+version+".png"
plt.savefig(nameDM0)

plt.figure()
plt.hist([probDM6_realDM0,probDM6_realDM4,probDM6_realDM6],bins, histtype='step',linewidth=1.5,fill=False, stacked=False,color=['red','blue','forestgreen'],weights=weights_probDM6, label=['True 1prong','True 3prong','True 1/3prong+pi0'])
plt.legend(loc='upper right')
plt.title("v7")
plt.xlabel("Predicted probability of decay mode 1/3prong+pi0")
plt.ylabel("Normalized entries")
nameDM0 = "./plots/probDM6_"+version+".png"
plt.savefig(nameDM0)

#plt.show()

# Create the confusion matrix
conf_matrix = pd.crosstab(y_test, preds, rownames=['Real DM'], colnames=['Predicted DM'])

# Plot the confusion matrix
def plot_confusion_matrix(y_true, y_pred, classes,
                          normalize=False,
                          title=None,
                          cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    # Only use the labels that appear in the data
    #classes = classes[unique_labels(y_true, y_pred)]
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    #print(cm)
    #print " "

    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           xticklabels=classes, yticklabels=classes,
           title=title,
           ylabel='True decay mode',
           xlabel='Predicted decay mode')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    fig.tight_layout()
    return ax

np.set_printoptions(precision=2)

plot_confusion_matrix(y_test, preds, classes = ["1prong","3prong","1/3prong+pi0"],
                      title='Confusion matrix, without normalization')
plt.title(version)
nameCM = "plots/conf_matrix_"+version+".png"
plt.savefig(nameCM)

plot_confusion_matrix(y_test, preds, classes= ["1prong","3prong","1/3prong+pi0"], normalize=True,
                      title='Normalized confusion matrix')
plt.title(version)
nameCMN = "plots/conf_matrix_norm_"+version+".png"
plt.savefig(nameCMN)

# Get the feature importances and sort them
feat_importances = clf.feature_importances_;
indices = np.argsort(feat_importances)[::-1]
ordered_feats = np.empty(n_features, dtype=object)

#print "Feature importance"
for f in xrange(n_features):
    ordered_feats[f] = features[indices[f]]
    #print f+1,'.',ordered_feats[f],'->',feat_importances[indices[f]]

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
plt.title(version)
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.3)
nameFI = "plots/feature_importance_"+version+".png"
plt.savefig(nameFI)

#print ' '
