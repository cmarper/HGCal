# C. Martin Perez cmartinp@cern.ch, Sep. 2019

###########

import os
import pandas as pd
import numpy as np
import root_pandas

###########

def deltar( df ):

    df['deta_cl3d_gentau'] = df['cl3d_eta'] - df['gentau_vis_eta']
    df['dphi_cl3d_gentau'] = np.abs(df['cl3d_phi'] - df['gentau_vis_phi'])
    sel = df['dphi_cl3d_gentau'] > np.pi
    df['dphi_cl3d_gentau'] -= sel*(2*np.pi)

    return ( np.sqrt(df['dphi_cl3d_gentau']*df['dphi_cl3d_gentau']+df['deta_cl3d_gentau']*df['deta_cl3d_gentau']) )

def taumatching( df_cl3d, df_gen, deta, dphi ):

	df_cl3d_plus 	= df_cl3d.query('cl3d_eta>0')
	df_cl3d_minus 	= df_cl3d.query('cl3d_eta<=0')

	df_gen_plus 	= df_gen.query('gentau_vis_eta>0')	
	df_gen_minus 	= df_gen.query('gentau_vis_eta<=0')

	df_cl3d_plus.set_index('event', inplace=True)
	df_cl3d_minus.set_index('event', inplace=True)

	df_gen_plus.set_index('event', inplace=True)
	df_gen_minus.set_index('event', inplace=True)

	df_merged_plus 	= df_gen_plus.join(df_cl3d_plus, how='left', rsuffix='_cl3d')
	df_merged_minus = df_gen_minus.join(df_cl3d_minus, how='left', rsuffix='_cl3d')

	df_merged_plus['deltar_cl3d_gentau'] 	= deltar(df_merged_plus)
	df_merged_minus['deltar_cl3d_gentau']	= deltar(df_merged_minus)

	sel_plus = np.abs(df_merged_plus['deta_cl3d_gentau']) < (deta/2)
	sel_minus = np.abs(df_merged_minus['deta_cl3d_gentau']) < (deta/2)

	df_merged_plus 	= df_merged_plus[sel_plus]
	df_merged_minus = df_merged_minus[sel_minus]

	sel_plus 	= np.abs(df_merged_plus['dphi_cl3d_gentau']) < (dphi/2)
	sel_minus 	= np.abs(df_merged_minus['dphi_cl3d_gentau']) < (dphi/2)

	df_merged_plus 	= df_merged_plus[sel_plus]
	df_merged_minus = df_merged_minus[sel_minus]

	group_plus 	= df_merged_plus.groupby('event')
	group_minus = df_merged_minus.groupby('event')

	n_cl3d_plus 	= group_plus['cl3d_pt'].size()
	n_cl3d_minus 	= group_minus['cl3d_pt'].size()

	df_merged_plus['n_matched_cl3d'] 	= n_cl3d_plus
	df_merged_minus['n_matched_cl3d'] 	= n_cl3d_minus

	df_merged_plus['bestmatch_pt'] 	= group_plus['cl3d_pt'].max()
	df_merged_minus['bestmatch_pt'] = group_minus['cl3d_pt'].max()

	df_merged_plus['cl3d_isbestmatch'] 	= df_merged_plus['bestmatch_pt'] == df_merged_plus['cl3d_pt']
	df_merged_minus['cl3d_isbestmatch'] = df_merged_minus['bestmatch_pt'] == df_merged_minus['cl3d_pt']

	df_merged = pd.concat([df_merged_plus, df_merged_minus], sort=False).sort_values('event')

	return df_merged


###########

dir_in 	= '/data_CMS/cms/mperez/HGCal_data/Sep19/skimmed/'
dir_out = '/data_CMS/cms/mperez/HGCal_data/Sep19/clustered/'

filein_tau_PU0 		= dir_in+'ntuple_RelValDiTau_Pt20To100_v10_PU0_skimmed.root'
filein_tau_PU200 	= dir_in+'ntuple_RelValDiTau_Pt20To100_v10_PU200_skimmed.root'
filein_nu_PU200 	= dir_in+'ntuple_Nu_E10_v10_PU200_skimmed.root'

fileout_tau_PU0		= dir_out+'RelValDiTau_Pt20To100_v10_PU0_matched.hdf5'
fileout_tau_PU200	= dir_out+'RelValDiTau_Pt20To100_v10_PU200_matched.hdf5'
fileout_nu_PU200 	= dir_out+'Nu_E10_v10_PU200_matched.hdf5'

deta_matching = 0.1
dphi_matching = 0.2

###########

treename = 'SkimmedTree'

branches_event_cl3d 	= ['event','cl3d_pt','cl3d_eta','cl3d_phi','cl3d_showerlength','cl3d_coreshowerlength',
	'cl3d_firstlayer','cl3d_maxlayer','cl3d_seetot','cl3d_spptot','cl3d_szz', 'cl3d_srrtot', 'cl3d_srrmean',
	'cl3d_hoe', 'cl3d_meanz', 'cl3d_layer10', 'cl3d_layer50', 'cl3d_layer90', 'cl3d_ntc67', 'cl3d_ntc90']
branches_cl3d 		= ['cl3d_pt','cl3d_eta','cl3d_phi','cl3d_showerlength','cl3d_coreshowerlength',
	'cl3d_firstlayer','cl3d_maxlayer','cl3d_seetot','cl3d_spptot','cl3d_szz', 'cl3d_srrtot', 'cl3d_srrmean',
	'cl3d_hoe', 'cl3d_meanz', 'cl3d_layer10', 'cl3d_layer50', 'cl3d_layer90', 'cl3d_ntc67', 'cl3d_ntc90']

branches_event_gentau = ['event', 'gentau_pt', 'gentau_eta', 'gentau_phi', 'gentau_energy', 'gentau_mass',
	'gentau_vis_pt', 'gentau_vis_eta', 'gentau_vis_phi', 'gentau_vis_energy', 'gentau_vis_mass',
	'gentau_decayMode']
branches_gentau 	= ['gentau_pt', 'gentau_eta', 'gentau_phi', 'gentau_energy', 'gentau_mass',
	'gentau_vis_pt', 'gentau_vis_eta', 'gentau_vis_phi', 'gentau_vis_energy', 'gentau_vis_mass',
	'gentau_decayMode']

df_tau_PU0_cl3d 	= root_pandas.read_root( filein_tau_PU0, 	key=treename, columns=branches_event_cl3d, flatten=branches_cl3d )
df_tau_PU0_gentau 	= root_pandas.read_root( filein_tau_PU0, 	key=treename, columns=branches_event_gentau, flatten=branches_gentau )

df_tau_PU200_cl3d 	= root_pandas.read_root( filein_tau_PU200, 	key=treename, columns=branches_event_cl3d, flatten=branches_cl3d )
df_tau_PU200_gentau = root_pandas.read_root( filein_tau_PU200, 	key=treename, columns=branches_event_gentau, flatten=branches_gentau )

df_nu_PU200_cl3d 	= root_pandas.read_root( filein_nu_PU200, 	key=treename, columns=branches_event_cl3d, flatten=branches_cl3d )

###########

df_tau_PU0 	 = taumatching( df_tau_PU0_cl3d, 	 df_tau_PU0_gentau,   deta_matching, dphi_matching )
df_tau_PU200 = taumatching( df_tau_PU200_cl3d, df_tau_PU200_gentau, deta_matching, dphi_matching )

#print 'PU=0: '
#print df_tau_PU0.head(5)
#print ' '
#print 'PU=200: '
#print df_tau_PU200.head(5)

###########

#save files to savedir
store_tau_PU0 = pd.HDFStore(fileout_tau_PU0, mode='w')
store_tau_PU0['df_tau_PU0'] = df_tau_PU0
store_tau_PU0.close()

store_tau_PU200 = pd.HDFStore(fileout_tau_PU200, mode='w')
store_tau_PU200['df_tau_PU200'] = df_tau_PU200
store_tau_PU200.close()

store_nu_PU200 = pd.HDFStore(fileout_nu_PU200, mode='w')
store_nu_PU200['df_nu_PU200'] = df_nu_PU200_cl3d
store_nu_PU200.close()

