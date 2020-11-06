# C. Martin Perez cmartinp@cern.ch, Nov. 2019

###########

import os
import numpy as np
import pandas as pd
import root_pandas

###########

# INPUTS

deta_matching = 0.1
dphi_matching = 0.2

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

dir_in 	= '/data_CMS_upgrade/mperez/HGCal_data/FE_options_Oct19/skimmed/'

file_in_tau = {}

file_in_tau[0] = dir_in+'ntuple_TauGun_Pt20_100_PU200_RelVal1060p2_decentralized_xyseed_pt2_THRESHOLD_skimmed.root'
file_in_tau[1] = dir_in+'ntuple_TauGun_Pt20_100_PU200_RelVal1060p2_decentralized_xyseed_pt2_SUPERTRIGGERCELL_skimmed.root'
file_in_tau[2] = dir_in+'ntuple_TauGun_Pt20_100_PU200_RelVal1060p2_decentralized_xyseed_pt2_BESTCHOICE_skimmed.root'
file_in_tau[3] = dir_in+'ntuple_TauGun_Pt20_100_PU200_RelVal1060p2_decentralized_xyseed_pt2_BESTCHOICECOARSE_skimmed.root'
file_in_tau[4] = dir_in+'ntuple_TauGun_Pt20_100_PU200_RelVal1060p2_decentralized_xyseed_pt2_MIXEDBCSTC_skimmed.root'

file_in_nu = {}

file_in_nu[0] = dir_in+'ntuple_NuGun_PU200_L1TSpring19_decentralized_xyseed_pt10_THRESHOLD_skimmed.root'
file_in_nu[1] = dir_in+'ntuple_NuGun_PU200_L1TSpring19_decentralized_xyseed_pt10_SUPERTRIGGERCELL_skimmed.root'
file_in_nu[2] = dir_in+'ntuple_NuGun_PU200_L1TSpring19_decentralized_xyseed_pt10_BESTCHOICE_skimmed.root'
file_in_nu[3] = dir_in+'ntuple_NuGun_PU200_L1TSpring19_decentralized_xyseed_pt10_BESTCHOICECOARSE_skimmed.root'
file_in_nu[4] = dir_in+'ntuple_NuGun_PU200_L1TSpring19_decentralized_xyseed_pt10_MIXEDBCSTC_skimmed.root'

###########

dir_out = '/data_CMS_upgrade/mperez/HGCal_data/thesis/matched/'

file_out_tau = {}

file_out_tau[0] = dir_out+'ntuple_TauGun_THRESHOLD_matched.hdf5'
file_out_tau[1] = dir_out+'ntuple_TauGun_SUPERTRIGGERCELL_matched.hdf5'
file_out_tau[2] = dir_out+'ntuple_TauGun_BESTCHOICE_matched.hdf5'
file_out_tau[3] = dir_out+'ntuple_TauGun_BESTCHOICECOARSE_matched.hdf5'
file_out_tau[4] = dir_out+'ntuple_TauGun_MIXEDBCSTC_matched.hdf5'

file_out_nu = {}

file_out_nu[0] = dir_out+'ntuple_NuGun_THRESHOLD_matched.hdf5'
file_out_nu[1] = dir_out+'ntuple_NuGun_SUPERTRIGGERCELL_matched.hdf5'
file_out_nu[2] = dir_out+'ntuple_NuGun_BESTCHOICE_matched.hdf5'
file_out_nu[3] = dir_out+'ntuple_NuGun_BESTCHOICECOARSE_matched.hdf5'
file_out_nu[4] = dir_out+'ntuple_NuGun_MIXEDBCSTC_matched.hdf5'

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

###########

df_tau = {}

for name in file_in_tau:

	df_tau_cl3d = root_pandas.read_root(file_in_tau[name], key=treename, columns=branches_event_cl3d, flatten=branches_cl3d)
	df_tau_gentau = root_pandas.read_root(file_in_tau[name], key=treename, columns=branches_event_gentau, flatten=branches_gentau)
	df_tau[name] = taumatching(df_tau_cl3d, df_tau_gentau, deta_matching, dphi_matching)


df_nu = {}

for name in file_in_nu:

	df_nu[name] = root_pandas.read_root(file_in_nu[name], key=treename, columns=branches_event_cl3d, flatten=branches_cl3d)

########

# Save files

for name in file_in_tau:

  store_tau = pd.HDFStore(file_out_tau[name], mode='w')
  store_tau['df_tau_PU200'] = df_tau[name]
  store_tau.close()

for name in file_in_nu:
	
  store_nu = pd.HDFStore(file_out_nu[name], mode='w')
  store_nu['df_nu_PU200'] = df_nu[name]
  store_nu.close()

########
