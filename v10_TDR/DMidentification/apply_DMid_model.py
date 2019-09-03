# C. Martin Perez cmartinp@cern.ch, Aug. 2019
# Training of multiclassifier BDT for L1 taus HGCal DM identification

from sklearn.externals import joblib
import pandas as pd
import numpy as np
import uproot

# Load the input variables from pandas
# from root file ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_clustered_flat.root
# inputs ['firstlayer', 'maxlayer', 'layer10', 'layer50', 'layer90', 'showerlength', 'hoe', 'meanz']

file = uproot.open('/data_CMS/cms/mperez/HGCal_data/Aug19/flat/ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU200_flat_PUcut_calibs.root')
file.keys()

tree = file["FlatTree"]
tree.keys

inputs = tree.pandas.df(["gentau_matchedSCL3D_firstlayer_seed","gentau_matchedSCL3D_maxlayer_seed","gentau_matchedSCL3D_layer10_seed", 
	"gentau_matchedSCL3D_layer50_seed","gentau_matchedSCL3D_layer90_seed","gentau_matchedSCL3D_showerlength_seed",
	"gentau_matchedSCL3D_hoe_seed","gentau_matchedSCL3D_meanz_seed"])

#print np.shape(inputs)

# Load the model trained previously
model = joblib.load('/home/llr/cms/mperez/HGCal/v10_geometry/data/DMid/model_train_DMid_8vars_etaPlus_PUcut.sav')

# Get predictions -> DM with the highest score (0,1,4,5)
DM = model.predict(inputs)
#print 'DM ', DM
#print 'shape DM ',np.shape(DM)

# Get probabilities -> probability for a specific DM
probs_DM = model.predict_proba(inputs)
#print 'preds_proba ', probs_DM
#print 'shape preds_proba ',np.shape(probs_DM)

#Correct for un-matched cases
ismatched = tree["gentau_isMatchedtoSCL3D"].array()
DM_corr = []
probDM0_corr = []
probDM1_corr = []
probDM4_corr = []
probDM5_corr = []

for i in range(0,len(ismatched)):
	if ismatched[i] == False:
		DM_corr.append(-999)
		probDM0_corr.append(-999)
		probDM1_corr.append(-999)
		probDM4_corr.append(-999)
		probDM5_corr.append(-999)
	elif ismatched[i] == True:
		DM_corr.append(int(DM[i]))
		probDM0_corr.append(probs_DM[i][0])
		probDM1_corr.append(probs_DM[i][1])
		probDM4_corr.append(probs_DM[i][2])
		probDM5_corr.append(probs_DM[i][3])

#for i in range(0,10):
#	print 'Pred DM: ',DM_corr[i],'; with probs',probDM0_corr[i],',',probDM1_corr[i],',',probDM4_corr[i],',',probDM5_corr[i]

#Dump the predicted DM (DM_corr) and probs of each DM (probDM*_corr) to a root file
from ROOT import TFile, TTree, gRandom
from array import array

f_in = TFile("/data_CMS/cms/mperez/HGCal_data/Aug19/flat/ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU200_flat_PUcut_calibs.root", 'read')
t_in = f_in.Get("FlatTree")

f_out = TFile("/data_CMS/cms/mperez/HGCal_data/Aug19/flat/ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU200_flat_PUcut_calibs_DM.root", 'recreate')
t_out = t_in.CloneTree()

N = t_out.GetEntries()

gentau_matchedSCL3D_predDM = array('i',[0])
gentau_matchedSCL3D_probDM0 = array('f',[0])
gentau_matchedSCL3D_probDM1 = array('f',[0])
gentau_matchedSCL3D_probDM4 = array('f',[0])
gentau_matchedSCL3D_probDM5 = array('f',[0])

b_predDM = t_out.Branch( "gentau_matchedSCL3D_predDM", gentau_matchedSCL3D_predDM, 'gentau_matchedSCL3D_predDM/I' )
b_probDM0 = t_out.Branch( "gentau_matchedSCL3D_probDM0", gentau_matchedSCL3D_probDM0, 'gentau_matchedSCL3D_probDM0/F' )
b_probDM1 = t_out.Branch( "gentau_matchedSCL3D_probDM1", gentau_matchedSCL3D_probDM1, 'gentau_matchedSCL3D_probDM1/F' )
b_probDM4 = t_out.Branch( "gentau_matchedSCL3D_probDM4", gentau_matchedSCL3D_probDM4, 'gentau_matchedSCL3D_probDM4/F' )
b_probDM5 = t_out.Branch( "gentau_matchedSCL3D_probDM5", gentau_matchedSCL3D_probDM5, 'gentau_matchedSCL3D_probDM5/F' )

for i in xrange(N):
	t_out.GetEntry(i)
	gentau_matchedSCL3D_predDM[0] = DM_corr[i]
	gentau_matchedSCL3D_probDM0[0] = probDM0_corr[i]
	gentau_matchedSCL3D_probDM1[0] = probDM1_corr[i]
	gentau_matchedSCL3D_probDM4[0] = probDM4_corr[i]
	gentau_matchedSCL3D_probDM5[0] = probDM5_corr[i]
	b_predDM.Fill()
	b_probDM0.Fill()
	b_probDM1.Fill()
	b_probDM4.Fill()
	b_probDM5.Fill()

t_out.Write()
f_out.Write("",TFile.kOverwrite)
f_in.Close()
f_out.Close()
