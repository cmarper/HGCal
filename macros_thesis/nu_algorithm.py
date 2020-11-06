# C. Martin Perez cmartinp@cern.ch, Sep. 2019

###########

import os
import pandas as pd
import numpy as np
import root_pandas

###########
###########

dir_in 	= '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoB/DMid/'
file_in = dir_in+'ntuple_NuGun_THRESHOLD_calibrated_pubdt99_DMid_algoB.hdf5'

store_nu = pd.HDFStore(file_in, mode='r')
df_nu  = store_nu['df_nu_PU200']
store_nu.close()

dir_out = '/data_CMS_upgrade/mperez/HGCal_data/thesis/algoB/rates/'
fileout_root = dir_out+'ntuple_NuGun_passPUbdt_leading.root'

fileout_pandas = dir_out+'ntuple_NuGun_passPUbdt_notleading.root'

###########
###########

# RATES

#Pass PU WP
sel1 = (df_nu['cl3d_pubdt_passWP99']==True)
df_nu_passpu = df_nu[sel1]

# Group by events
df_nu_passpu.set_index('event', inplace=True)
group = df_nu_passpu.groupby('event')
#print 'all',np.shape(df_nu['cl3d_pt_c3'])
#print df_nu
#print '...'

#Is leading
df_nu_passpu['leading_cl3d_pt_c3'] = group['cl3d_pt_c3'].max()
df_nu_passpu['cl3d_isleading'] = (df_nu_passpu['leading_cl3d_pt_c3'] == df_nu_passpu['cl3d_pt_c3'])

sel2 = (df_nu_passpu['cl3d_isleading']==True)
df_nu_passpu_leading = df_nu_passpu[sel2]

print df_nu_passpu_leading

#To root
df_nu_passpu_leading.reset_index()['event']
df_nu_passpu_leading.to_root(fileout_root, key='Tree', store_index=True, mode='a') #events stored in branch __index__event

