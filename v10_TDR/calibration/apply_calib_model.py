# C. Martin Perez cmartinp@cern.ch, Aug. 2019
# Training of regression BDT for L1 taus HGCal calibration

from sklearn.externals import joblib
import pandas as pd
import numpy as np
import uproot


# Load the input variables from pandas
# from root file ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_clustered.root
# inputs ['pt_tot', 'eta_Eweighted','firstlayer', 'maxlayer', 'layer10', 'layer50', 'layer90', 'showerlength', 'hoe', 'meanz']

file = uproot.open('/data_CMS/cms/mperez/HGCal_data/Aug19/flat/ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_flat.root')
file.keys()
tree = file["FlatTree"]
tree.keys

inputs = tree.pandas.df(["gentau_matchedSC_pt_tot","gentau_matchedSC_eta_Eweighted",
	"gentau_matchedSC_firstlayer_seed","gentau_matchedSC_maxlayer_seed","gentau_matchedSC_layer10_seed", 
	"gentau_matchedSC_layer50_seed","gentau_matchedSC_layer90_seed","gentau_matchedSC_showerlength_seed",
	"gentau_matchedSC_hoe_seed","gentau_matchedSC_meanz_seed"])

#print inputs
#print np.shape(inputs)

# Load the model trained previously
model = joblib.load('/home/llr/cms/mperez/HGCal/v10_geometry/data/model_train_calib_10vars_etaPlus.sav')

# Predict the calibration factor
calib_factors = model.predict(inputs)

# Correct for un-matched cases
ismatched = tree["gentau_isMatched"].array()
pt_raw = tree["gentau_matchedSC_pt_tot"].array()
calib_factors_corr = []
pt_calib = []

for i in range(0,len(ismatched)):
	if ismatched[i] == False:
		calib_factors_corr.append(-999)
		pt_calib.append(-999)
	elif ismatched[i] == True:
		calib_factors_corr.append(calib_factors[i])
		pt_calib.append(pt_raw[i]*calib_factors[i])

#print np.shape(pt_raw), np.shape(calib_factors), np.shape(pt_calib)
#for i in range(0,10):
#	print pt_raw[i],'*',calib_factors_corr[i],'=',pt_calib[i]

# Dump the calibration factors (calib_factors_corr) and calibrated pt (pt_calib) to a root file
from ROOT import TFile, TTree, gRandom
from array import array

f_in = TFile("/data_CMS/cms/mperez/HGCal_data/Aug19/flat/ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_flat.root", 'read')
t_in = f_in.Get("FlatTree")

f_out = TFile("/data_CMS/cms/mperez/HGCal_data/Aug19/flat/ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_flat_calibs.root", 'recreate')
t_out = t_in.CloneTree()

N = t_out.GetEntries()
#print N

gentau_matchedSC_calib_factor = array('f',[0])
gentau_matchedSC_pt_tot_calib = array('f',[0])

b_calib_factor = t_out.Branch( "gentau_matchedSC_calib_factor", gentau_matchedSC_calib_factor, 'gentau_matchedSC_calib_factor/F' )
b_pt_tot_calib = t_out.Branch( "gentau_matchedSC_pt_tot_calib", gentau_matchedSC_pt_tot_calib, 'gentau_matchedSC_pt_tot_calib/F' )

for i in xrange(N):
	t_out.GetEntry(i)
	gentau_matchedSC_calib_factor[0] = calib_factors_corr[i]
	gentau_matchedSC_pt_tot_calib[0] = pt_calib[i]
	b_calib_factor.Fill()
	b_pt_tot_calib.Fill()

t_out.Write()
f_out.Write("",TFile.kOverwrite)
f_in.Close()
f_out.Close()
