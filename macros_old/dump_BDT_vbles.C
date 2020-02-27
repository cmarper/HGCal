/////////////////////////////////////////////////////////
///// HGCal L1 taus, C. Martin Perez, LLR, May 2019 /////
/////////////////////////////////////////////////////////

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <iostream>
#include <fstream>

using namespace std;


void dump_variables_classifierDM(int nevents = -1){
	
	ofstream out_file_v1;

	out_file_v1.open("/home/llr/cms/mperez/HGCal/v10_geometry/data/v10geometry_DMclassifier_Aug19_v1.csv");

	TChain * in_tree = new TChain("ClusteredTree");	
	in_tree->Add("/data_CMS/cms/mperez/HGCal_data/Aug19/ntuple_RelValDiTau_Pt20To100_clustered.root");

	Long64_t nentries = in_tree->GetEntries();
	cout<<"nentries="<<in_tree->GetEntries()<<endl;
	if (nevents != -1) 
		nentries = nevents;

	// old branches used

	int _gentau_n;

	vector<bool> 	*_gentau_isMatched;

    vector<int>    *_gentau_matchedSC_n_cl3d;
    vector<float>  *_gentau_matchedSC_pt_tot;
    vector<float>  *_gentau_matchedSC_pt_seed;
    vector<float>  *_gentau_matchedSC_eta_Eweighted;
    vector<float>  *_gentau_matchedSC_eta_seed;
    vector<float>  *_gentau_matchedSC_phi_Eweighted;
    vector<float>  *_gentau_matchedSC_phi_seed;

    vector<int>    *_gentau_matchedSC_showerlength_seed;
    vector<int>    *_gentau_matchedSC_coreshowerlength_seed;
    vector<int>    *_gentau_matchedSC_firstlayer_seed;
    vector<int>    *_gentau_matchedSC_maxlayer_seed;
    vector<float>  *_gentau_matchedSC_seetot_seed;
    vector<float>  *_gentau_matchedSC_seemax_seed;
    vector<float>  *_gentau_matchedSC_spptot_seed;
    vector<float>  *_gentau_matchedSC_sppmax_seed;
    vector<float>  *_gentau_matchedSC_szz_seed;
    vector<float>  *_gentau_matchedSC_srrtot_seed;
    vector<float>  *_gentau_matchedSC_srrmax_seed;
    vector<float>  *_gentau_matchedSC_srrmean_seed;
    vector<float>  *_gentau_matchedSC_emaxe_seed;
    vector<float>  *_gentau_matchedSC_hoe_seed;
    vector<float>  *_gentau_matchedSC_meanz_seed;
    vector<float>  *_gentau_matchedSC_layer10_seed;
    vector<float>  *_gentau_matchedSC_layer50_seed;
    vector<float>  *_gentau_matchedSC_layer90_seed;
    vector<float>  *_gentau_matchedSC_ntc67_seed;
    vector<float>  *_gentau_matchedSC_ntc90_seed;
    vector<float>  *_gentau_matchedSC_bdteg_seed;
    vector<int>    *_gentau_matchedSC_quality_seed;

	vector<int>   *_gentau_decayMode;

	in_tree->SetBranchAddress("gentau_n",	&_gentau_n);
	in_tree->SetBranchAddress("gentau_isMatched",	&_gentau_isMatched);

    in_tree->SetBranchAddress("gentau_matchedSC_n_cl3d", 		&_gentau_matchedSC_n_cl3d);
    in_tree->SetBranchAddress("gentau_matchedSC_pt_tot", 		&_gentau_matchedSC_pt_tot);
    in_tree->SetBranchAddress("gentau_matchedSC_pt_seed", 		&_gentau_matchedSC_pt_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_eta_Eweighted", &_gentau_matchedSC_eta_Eweighted);
    in_tree->SetBranchAddress("gentau_matchedSC_eta_seed", 		&_gentau_matchedSC_eta_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_phi_Eweighted", &_gentau_matchedSC_phi_Eweighted);
    in_tree->SetBranchAddress("gentau_matchedSC_phi_seed", 		&_gentau_matchedSC_phi_seed);

    in_tree->SetBranchAddress("gentau_matchedSC_showerlength_seed", 	&_gentau_matchedSC_showerlength_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_coreshowerlength_seed", &_gentau_matchedSC_coreshowerlength_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_firstlayer_seed", 		&_gentau_matchedSC_firstlayer_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_maxlayer_seed", 		&_gentau_matchedSC_maxlayer_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_seetot_seed", 			&_gentau_matchedSC_seetot_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_seemax_seed", 			&_gentau_matchedSC_seemax_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_spptot_seed", 			&_gentau_matchedSC_spptot_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_sppmax_seed", 			&_gentau_matchedSC_sppmax_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_szz_seed", 				&_gentau_matchedSC_szz_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_srrtot_seed", 			&_gentau_matchedSC_srrtot_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_srrmax_seed", 			&_gentau_matchedSC_srrmax_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_srrmean_seed", 			&_gentau_matchedSC_srrmean_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_emaxe_seed", 			&_gentau_matchedSC_emaxe_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_hoe_seed", 				&_gentau_matchedSC_hoe_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_meanz_seed", 			&_gentau_matchedSC_meanz_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_layer10_seed", 			&_gentau_matchedSC_layer10_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_layer50_seed", 			&_gentau_matchedSC_layer50_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_layer90_seed", 			&_gentau_matchedSC_layer90_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_ntc67_seed", 			&_gentau_matchedSC_ntc67_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_ntc90_seed", 			&_gentau_matchedSC_ntc90_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_bdteg_seed", 			&_gentau_matchedSC_bdteg_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_quality_seed", 			&_gentau_matchedSC_quality_seed);

    in_tree->SetBranchAddress("gentau_decayMode",	&_gentau_decayMode);


	for (int i=0;i<nentries;i++) {

		if(i%1000==0) cout<<"i="<<i<<endl;

		//old branches

		_gentau_n = 0;

		_gentau_isMatched = 0;

        _gentau_matchedSC_n_cl3d = 0;
        _gentau_matchedSC_pt_tot = 0;
        _gentau_matchedSC_pt_seed = 0;
        _gentau_matchedSC_eta_Eweighted = 0;
        _gentau_matchedSC_eta_seed = 0;
        _gentau_matchedSC_phi_Eweighted = 0;
        _gentau_matchedSC_phi_seed = 0;

        _gentau_matchedSC_showerlength_seed = 0;
        _gentau_matchedSC_coreshowerlength_seed = 0;
        _gentau_matchedSC_firstlayer_seed = 0;
        _gentau_matchedSC_maxlayer_seed = 0;
        _gentau_matchedSC_seetot_seed = 0;
        _gentau_matchedSC_seemax_seed = 0;
        _gentau_matchedSC_spptot_seed = 0;
        _gentau_matchedSC_sppmax_seed = 0;
        _gentau_matchedSC_szz_seed = 0;
        _gentau_matchedSC_srrtot_seed = 0;
        _gentau_matchedSC_srrmax_seed = 0;
        _gentau_matchedSC_srrmean_seed = 0;
        _gentau_matchedSC_emaxe_seed = 0;
        _gentau_matchedSC_hoe_seed = 0;
        _gentau_matchedSC_meanz_seed = 0;
        _gentau_matchedSC_layer10_seed = 0;
        _gentau_matchedSC_layer50_seed = 0;
        _gentau_matchedSC_layer90_seed = 0;
        _gentau_matchedSC_ntc67_seed = 0;
        _gentau_matchedSC_ntc90_seed = 0;
        _gentau_matchedSC_bdteg_seed = 0;
        _gentau_matchedSC_quality_seed = 0;

		_gentau_decayMode = 0;

		int entry_ok = in_tree->GetEntry(i);	
		if(entry_ok<0) 
			continue;

		for(int i_gentau=0; i_gentau<_gentau_n; i_gentau++){

			if (!(*_gentau_isMatched)[i_gentau]) continue;

			if ((*_gentau_decayMode)[i_gentau] == 0 || (*_gentau_decayMode)[i_gentau] == 1 || (*_gentau_decayMode)[i_gentau] == 4 || (*_gentau_decayMode)[i_gentau] == 5){


				int    f_gentau_matchedSC_n_cl3d = (*_gentau_matchedSC_n_cl3d)[i_gentau];
    			float  f_gentau_matchedSC_pt_tot = (*_gentau_matchedSC_pt_tot)[i_gentau];
    			float  f_gentau_matchedSC_pt_seed = (*_gentau_matchedSC_pt_seed)[i_gentau];
    			float  f_gentau_matchedSC_eta_Eweighted = (*_gentau_matchedSC_eta_Eweighted)[i_gentau];
    			float  f_gentau_matchedSC_eta_seed = (*_gentau_matchedSC_eta_seed)[i_gentau];
    			float  f_gentau_matchedSC_phi_Eweighted = (*_gentau_matchedSC_phi_Eweighted)[i_gentau];
    			float  f_gentau_matchedSC_phi_seed = (*_gentau_matchedSC_phi_seed)[i_gentau];

    			int    f_gentau_matchedSC_showerlength_seed = (*_gentau_matchedSC_showerlength_seed)[i_gentau];
    			int    f_gentau_matchedSC_coreshowerlength_seed = (*_gentau_matchedSC_coreshowerlength_seed)[i_gentau];
    			int    f_gentau_matchedSC_firstlayer_seed = (*_gentau_matchedSC_firstlayer_seed)[i_gentau];
    			int    f_gentau_matchedSC_maxlayer_seed = (*_gentau_matchedSC_maxlayer_seed)[i_gentau];
    			float  f_gentau_matchedSC_seetot_seed = (*_gentau_matchedSC_seetot_seed)[i_gentau];
    			float  f_gentau_matchedSC_seemax_seed = (*_gentau_matchedSC_seemax_seed)[i_gentau];
    			float  f_gentau_matchedSC_spptot_seed = (*_gentau_matchedSC_spptot_seed)[i_gentau];
    			float  f_gentau_matchedSC_sppmax_seed = (*_gentau_matchedSC_sppmax_seed)[i_gentau];
    			float  f_gentau_matchedSC_szz_seed = (*_gentau_matchedSC_szz_seed)[i_gentau];
    			float  f_gentau_matchedSC_srrtot_seed = (*_gentau_matchedSC_srrtot_seed)[i_gentau];
    			float  f_gentau_matchedSC_srrmax_seed = (*_gentau_matchedSC_srrmax_seed)[i_gentau];
    			float  f_gentau_matchedSC_srrmean_seed = (*_gentau_matchedSC_srrmean_seed)[i_gentau];
    			float  f_gentau_matchedSC_emaxe_seed = (*_gentau_matchedSC_emaxe_seed)[i_gentau];
    			float  f_gentau_matchedSC_hoe_seed = (*_gentau_matchedSC_hoe_seed)[i_gentau];
    			float  f_gentau_matchedSC_meanz_seed = (*_gentau_matchedSC_meanz_seed)[i_gentau];
    			float  f_gentau_matchedSC_layer10_seed = (*_gentau_matchedSC_layer10_seed)[i_gentau];
    			float  f_gentau_matchedSC_layer50_seed = (*_gentau_matchedSC_layer50_seed)[i_gentau];
    			float  f_gentau_matchedSC_layer90_seed = (*_gentau_matchedSC_layer90_seed)[i_gentau];
    			float  f_gentau_matchedSC_ntc67_seed = (*_gentau_matchedSC_ntc67_seed)[i_gentau];
    			float  f_gentau_matchedSC_ntc90_seed = (*_gentau_matchedSC_ntc90_seed)[i_gentau];
    			float  f_gentau_matchedSC_bdteg_seed = (*_gentau_matchedSC_bdteg_seed)[i_gentau];
    			int    f_gentau_matchedSC_quality_seed = (*_gentau_matchedSC_quality_seed)[i_gentau];

    			int f_gentau_decayMode = (*_gentau_decayMode)[i_gentau];

				int f_gentau_decayMode_joined;

				if (f_gentau_decayMode == 1 || f_gentau_decayMode == 5) f_gentau_decayMode_joined = 6;											
				else f_gentau_decayMode_joined = f_gentau_decayMode;

				TString s_gentau_matchedSC_n_cl3d = Form("%.0i", f_gentau_matchedSC_n_cl3d);
				TString s_gentau_matchedSC_pt_tot = Form("%.5f", f_gentau_matchedSC_pt_tot);
				TString s_gentau_matchedSC_pt_seed = Form("%.5f", f_gentau_matchedSC_pt_seed);
				TString s_gentau_matchedSC_eta_Eweighted = Form("%.5f", f_gentau_matchedSC_eta_Eweighted);
				TString s_gentau_matchedSC_eta_seed = Form("%.5f", f_gentau_matchedSC_eta_seed);
				TString s_gentau_matchedSC_phi_Eweighted = Form("%.5f", f_gentau_matchedSC_phi_Eweighted);
				TString s_gentau_matchedSC_phi_seed = Form("%.5f", f_gentau_matchedSC_phi_seed);

				TString s_gentau_matchedSC_showerlength_seed = Form("%.1i", f_gentau_matchedSC_showerlength_seed);
				TString s_gentau_matchedSC_coreshowerlength_seed = Form("%.1i", f_gentau_matchedSC_coreshowerlength_seed);
				TString s_gentau_matchedSC_firstlayer_seed = Form("%.1i", f_gentau_matchedSC_firstlayer_seed);
				TString s_gentau_matchedSC_maxlayer_seed  = Form("%.1i", f_gentau_matchedSC_maxlayer_seed);
				TString s_gentau_matchedSC_seetot_seed = Form("%.5f", f_gentau_matchedSC_seetot_seed);
    			TString s_gentau_matchedSC_seemax_seed = Form("%.5f", f_gentau_matchedSC_seemax_seed);
    			TString s_gentau_matchedSC_spptot_seed = Form("%.5f", f_gentau_matchedSC_spptot_seed);
    			TString s_gentau_matchedSC_sppmax_seed = Form("%.5f", f_gentau_matchedSC_sppmax_seed);
    			TString s_gentau_matchedSC_szz_seed = Form("%.5f", f_gentau_matchedSC_szz_seed);
    			TString s_gentau_matchedSC_srrtot_seed = Form("%.5f", f_gentau_matchedSC_srrtot_seed);
    			TString s_gentau_matchedSC_srrmax_seed = Form("%.5f", f_gentau_matchedSC_srrmax_seed);
    			TString s_gentau_matchedSC_srrmean_seed = Form("%.5f", f_gentau_matchedSC_srrmean_seed);
    			TString s_gentau_matchedSC_emaxe_seed = Form("%.5f", f_gentau_matchedSC_emaxe_seed);
    			TString s_gentau_matchedSC_hoe_seed = Form("%.5f", f_gentau_matchedSC_hoe_seed);
    			TString s_gentau_matchedSC_meanz_seed = Form("%.5f", f_gentau_matchedSC_meanz_seed);
    			TString s_gentau_matchedSC_layer10_seed = Form("%.5f", f_gentau_matchedSC_layer10_seed);
    			TString s_gentau_matchedSC_layer50_seed = Form("%.5f", f_gentau_matchedSC_layer50_seed);
    			TString s_gentau_matchedSC_layer90_seed = Form("%.5f", f_gentau_matchedSC_layer90_seed);
    			TString s_gentau_matchedSC_ntc67_seed = Form("%.5f", f_gentau_matchedSC_ntc67_seed);
    			TString s_gentau_matchedSC_ntc90_seed = Form("%.5f", f_gentau_matchedSC_ntc90_seed);
    			TString s_gentau_matchedSC_bdteg_seed = Form("%.5f", f_gentau_matchedSC_bdteg_seed);
    			TString s_gentau_matchedSC_quality_seed = Form("%.1i", f_gentau_matchedSC_quality_seed);

				TString s_gentau_decayMode = Form("%.1i", f_gentau_decayMode);
				TString s_gentau_decayMode_joined = Form("%.1i", f_gentau_decayMode_joined);
	
				//8+1
				TString dumper_v1 = s_gentau_matchedSC_firstlayer_seed+","+s_gentau_matchedSC_maxlayer_seed+","+s_gentau_matchedSC_layer10_seed+","+s_gentau_matchedSC_layer50_seed+","+s_gentau_matchedSC_layer90_seed+","+s_gentau_matchedSC_showerlength_seed+","+s_gentau_matchedSC_hoe_seed+","+s_gentau_matchedSC_meanz_seed+","+s_gentau_decayMode;
				out_file_v1<<dumper_v1<<'\n';
				
			}
		}
	}

	out_file_v1.close();

	return;
}


void dump_variables_regressorCalib(int nevents = -1){
	
	ofstream out_file_v2;

	out_file_v2.open("/home/llr/cms/mperez/HGCal/v10_geometry/data/input_train_calib_10vars_etaPlus.csv");

	TChain * in_tree = new TChain("ClusteredTree");	
	in_tree->Add("/data_CMS/cms/mperez/HGCal_data/Aug19/ntuple_RelValDiTau_Pt20To100_Eta1p6To2p9_clustered.root");

	Long64_t nentries = in_tree->GetEntries();
	cout<<"nentries="<<in_tree->GetEntries()<<endl;
	if (nevents != -1) 
		nentries = nevents;

	// old branches used

	int _gentau_n;

	vector<bool> 	*_gentau_isMatched;

    vector<int>    *_gentau_matchedSC_n_cl3d;
    vector<float>  *_gentau_matchedSC_pt_tot;
    vector<float>  *_gentau_matchedSC_pt_seed;
    vector<float>  *_gentau_matchedSC_eta_Eweighted;
    vector<float>  *_gentau_matchedSC_eta_seed;
    vector<float>  *_gentau_matchedSC_phi_Eweighted;
    vector<float>  *_gentau_matchedSC_phi_seed;

    vector<int>    *_gentau_matchedSC_showerlength_seed;
    vector<int>    *_gentau_matchedSC_coreshowerlength_seed;
    vector<int>    *_gentau_matchedSC_firstlayer_seed;
    vector<int>    *_gentau_matchedSC_maxlayer_seed;
    vector<float>  *_gentau_matchedSC_seetot_seed;
    vector<float>  *_gentau_matchedSC_seemax_seed;
    vector<float>  *_gentau_matchedSC_spptot_seed;
    vector<float>  *_gentau_matchedSC_sppmax_seed;
    vector<float>  *_gentau_matchedSC_szz_seed;
    vector<float>  *_gentau_matchedSC_srrtot_seed;
    vector<float>  *_gentau_matchedSC_srrmax_seed;
    vector<float>  *_gentau_matchedSC_srrmean_seed;
    vector<float>  *_gentau_matchedSC_emaxe_seed;
    vector<float>  *_gentau_matchedSC_hoe_seed;
    vector<float>  *_gentau_matchedSC_meanz_seed;
    vector<float>  *_gentau_matchedSC_layer10_seed;
    vector<float>  *_gentau_matchedSC_layer50_seed;
    vector<float>  *_gentau_matchedSC_layer90_seed;
    vector<float>  *_gentau_matchedSC_ntc67_seed;
    vector<float>  *_gentau_matchedSC_ntc90_seed;
    vector<float>  *_gentau_matchedSC_bdteg_seed;
    vector<int>    *_gentau_matchedSC_quality_seed;

	vector<int>   *_gentau_decayMode;
    vector<float> *_gentau_vis_pt;

	in_tree->SetBranchAddress("gentau_n",	&_gentau_n);
	in_tree->SetBranchAddress("gentau_isMatched",	&_gentau_isMatched);

    in_tree->SetBranchAddress("gentau_matchedSC_n_cl3d", 		&_gentau_matchedSC_n_cl3d);
    in_tree->SetBranchAddress("gentau_matchedSC_pt_tot", 		&_gentau_matchedSC_pt_tot);
    in_tree->SetBranchAddress("gentau_matchedSC_pt_seed", 		&_gentau_matchedSC_pt_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_eta_Eweighted", &_gentau_matchedSC_eta_Eweighted);
    in_tree->SetBranchAddress("gentau_matchedSC_eta_seed", 		&_gentau_matchedSC_eta_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_phi_Eweighted", &_gentau_matchedSC_phi_Eweighted);
    in_tree->SetBranchAddress("gentau_matchedSC_phi_seed", 		&_gentau_matchedSC_phi_seed);

    in_tree->SetBranchAddress("gentau_matchedSC_showerlength_seed", 	&_gentau_matchedSC_showerlength_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_coreshowerlength_seed", &_gentau_matchedSC_coreshowerlength_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_firstlayer_seed", 		&_gentau_matchedSC_firstlayer_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_maxlayer_seed", 		&_gentau_matchedSC_maxlayer_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_seetot_seed", 			&_gentau_matchedSC_seetot_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_seemax_seed", 			&_gentau_matchedSC_seemax_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_spptot_seed", 			&_gentau_matchedSC_spptot_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_sppmax_seed", 			&_gentau_matchedSC_sppmax_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_szz_seed", 				&_gentau_matchedSC_szz_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_srrtot_seed", 			&_gentau_matchedSC_srrtot_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_srrmax_seed", 			&_gentau_matchedSC_srrmax_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_srrmean_seed", 			&_gentau_matchedSC_srrmean_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_emaxe_seed", 			&_gentau_matchedSC_emaxe_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_hoe_seed", 				&_gentau_matchedSC_hoe_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_meanz_seed", 			&_gentau_matchedSC_meanz_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_layer10_seed", 			&_gentau_matchedSC_layer10_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_layer50_seed", 			&_gentau_matchedSC_layer50_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_layer90_seed", 			&_gentau_matchedSC_layer90_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_ntc67_seed", 			&_gentau_matchedSC_ntc67_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_ntc90_seed", 			&_gentau_matchedSC_ntc90_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_bdteg_seed", 			&_gentau_matchedSC_bdteg_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_quality_seed", 			&_gentau_matchedSC_quality_seed);

    in_tree->SetBranchAddress("gentau_decayMode",	&_gentau_decayMode);
    in_tree->SetBranchAddress("gentau_vis_pt", &_gentau_vis_pt);


	for (int i=0;i<nentries;i++) {

		if(i%1000==0) cout<<"i="<<i<<endl;

		//old branches

		_gentau_n = 0;

		_gentau_isMatched = 0;

        _gentau_matchedSC_n_cl3d = 0;
        _gentau_matchedSC_pt_tot = 0;
        _gentau_matchedSC_pt_seed = 0;
        _gentau_matchedSC_eta_Eweighted = 0;
        _gentau_matchedSC_eta_seed = 0;
        _gentau_matchedSC_phi_Eweighted = 0;
        _gentau_matchedSC_phi_seed = 0;

        _gentau_matchedSC_showerlength_seed = 0;
        _gentau_matchedSC_coreshowerlength_seed = 0;
        _gentau_matchedSC_firstlayer_seed = 0;
        _gentau_matchedSC_maxlayer_seed = 0;
        _gentau_matchedSC_seetot_seed = 0;
        _gentau_matchedSC_seemax_seed = 0;
        _gentau_matchedSC_spptot_seed = 0;
        _gentau_matchedSC_sppmax_seed = 0;
        _gentau_matchedSC_szz_seed = 0;
        _gentau_matchedSC_srrtot_seed = 0;
        _gentau_matchedSC_srrmax_seed = 0;
        _gentau_matchedSC_srrmean_seed = 0;
        _gentau_matchedSC_emaxe_seed = 0;
        _gentau_matchedSC_hoe_seed = 0;
        _gentau_matchedSC_meanz_seed = 0;
        _gentau_matchedSC_layer10_seed = 0;
        _gentau_matchedSC_layer50_seed = 0;
        _gentau_matchedSC_layer90_seed = 0;
        _gentau_matchedSC_ntc67_seed = 0;
        _gentau_matchedSC_ntc90_seed = 0;
        _gentau_matchedSC_bdteg_seed = 0;
        _gentau_matchedSC_quality_seed = 0;

		_gentau_decayMode = 0;
        _gentau_vis_pt = 0;

		int entry_ok = in_tree->GetEntry(i);	
		if(entry_ok<0) 
			continue;

		for(int i_gentau=0; i_gentau<_gentau_n; i_gentau++){

			if (!(*_gentau_isMatched)[i_gentau]) continue;

			if ((*_gentau_decayMode)[i_gentau] == 0 || (*_gentau_decayMode)[i_gentau] == 1 || (*_gentau_decayMode)[i_gentau] == 4 || (*_gentau_decayMode)[i_gentau] == 5){


				int    f_gentau_matchedSC_n_cl3d = (*_gentau_matchedSC_n_cl3d)[i_gentau];
    			float  f_gentau_matchedSC_pt_tot = (*_gentau_matchedSC_pt_tot)[i_gentau];
    			float  f_gentau_matchedSC_pt_seed = (*_gentau_matchedSC_pt_seed)[i_gentau];
    			float  f_gentau_matchedSC_eta_Eweighted = (*_gentau_matchedSC_eta_Eweighted)[i_gentau];
    			float  f_gentau_matchedSC_eta_seed = (*_gentau_matchedSC_eta_seed)[i_gentau];
    			float  f_gentau_matchedSC_phi_Eweighted = (*_gentau_matchedSC_phi_Eweighted)[i_gentau];
    			float  f_gentau_matchedSC_phi_seed = (*_gentau_matchedSC_phi_seed)[i_gentau];

    			int    f_gentau_matchedSC_showerlength_seed = (*_gentau_matchedSC_showerlength_seed)[i_gentau];
    			int    f_gentau_matchedSC_coreshowerlength_seed = (*_gentau_matchedSC_coreshowerlength_seed)[i_gentau];
    			int    f_gentau_matchedSC_firstlayer_seed = (*_gentau_matchedSC_firstlayer_seed)[i_gentau];
    			int    f_gentau_matchedSC_maxlayer_seed = (*_gentau_matchedSC_maxlayer_seed)[i_gentau];
    			float  f_gentau_matchedSC_seetot_seed = (*_gentau_matchedSC_seetot_seed)[i_gentau];
    			float  f_gentau_matchedSC_seemax_seed = (*_gentau_matchedSC_seemax_seed)[i_gentau];
    			float  f_gentau_matchedSC_spptot_seed = (*_gentau_matchedSC_spptot_seed)[i_gentau];
    			float  f_gentau_matchedSC_sppmax_seed = (*_gentau_matchedSC_sppmax_seed)[i_gentau];
    			float  f_gentau_matchedSC_szz_seed = (*_gentau_matchedSC_szz_seed)[i_gentau];
    			float  f_gentau_matchedSC_srrtot_seed = (*_gentau_matchedSC_srrtot_seed)[i_gentau];
    			float  f_gentau_matchedSC_srrmax_seed = (*_gentau_matchedSC_srrmax_seed)[i_gentau];
    			float  f_gentau_matchedSC_srrmean_seed = (*_gentau_matchedSC_srrmean_seed)[i_gentau];
    			float  f_gentau_matchedSC_emaxe_seed = (*_gentau_matchedSC_emaxe_seed)[i_gentau];
    			float  f_gentau_matchedSC_hoe_seed = (*_gentau_matchedSC_hoe_seed)[i_gentau];
    			float  f_gentau_matchedSC_meanz_seed = (*_gentau_matchedSC_meanz_seed)[i_gentau];
    			float  f_gentau_matchedSC_layer10_seed = (*_gentau_matchedSC_layer10_seed)[i_gentau];
    			float  f_gentau_matchedSC_layer50_seed = (*_gentau_matchedSC_layer50_seed)[i_gentau];
    			float  f_gentau_matchedSC_layer90_seed = (*_gentau_matchedSC_layer90_seed)[i_gentau];
    			float  f_gentau_matchedSC_ntc67_seed = (*_gentau_matchedSC_ntc67_seed)[i_gentau];
    			float  f_gentau_matchedSC_ntc90_seed = (*_gentau_matchedSC_ntc90_seed)[i_gentau];
    			float  f_gentau_matchedSC_bdteg_seed = (*_gentau_matchedSC_bdteg_seed)[i_gentau];
    			int    f_gentau_matchedSC_quality_seed = (*_gentau_matchedSC_quality_seed)[i_gentau];

                int f_gentau_decayMode = (*_gentau_decayMode)[i_gentau];

                int f_gentau_decayMode_joined;

                float f_gentau_vis_pt = (*_gentau_vis_pt)[i_gentau];
                float f_target = f_gentau_vis_pt/f_gentau_matchedSC_pt_tot;

				if (f_gentau_decayMode == 1 || f_gentau_decayMode == 5) f_gentau_decayMode_joined = 6;											
				else f_gentau_decayMode_joined = f_gentau_decayMode;

				TString s_gentau_matchedSC_n_cl3d = Form("%.0i", f_gentau_matchedSC_n_cl3d);
				TString s_gentau_matchedSC_pt_tot = Form("%.5f", f_gentau_matchedSC_pt_tot);
				TString s_gentau_matchedSC_pt_seed = Form("%.5f", f_gentau_matchedSC_pt_seed);
				TString s_gentau_matchedSC_eta_Eweighted = Form("%.5f", f_gentau_matchedSC_eta_Eweighted);
				TString s_gentau_matchedSC_eta_seed = Form("%.5f", f_gentau_matchedSC_eta_seed);
				TString s_gentau_matchedSC_phi_Eweighted = Form("%.5f", f_gentau_matchedSC_phi_Eweighted);
				TString s_gentau_matchedSC_phi_seed = Form("%.5f", f_gentau_matchedSC_phi_seed);

				TString s_gentau_matchedSC_showerlength_seed = Form("%.1i", f_gentau_matchedSC_showerlength_seed);
				TString s_gentau_matchedSC_coreshowerlength_seed = Form("%.1i", f_gentau_matchedSC_coreshowerlength_seed);
				TString s_gentau_matchedSC_firstlayer_seed = Form("%.1i", f_gentau_matchedSC_firstlayer_seed);
				TString s_gentau_matchedSC_maxlayer_seed  = Form("%.1i", f_gentau_matchedSC_maxlayer_seed);
				TString s_gentau_matchedSC_seetot_seed = Form("%.5f", f_gentau_matchedSC_seetot_seed);
    			TString s_gentau_matchedSC_seemax_seed = Form("%.5f", f_gentau_matchedSC_seemax_seed);
    			TString s_gentau_matchedSC_spptot_seed = Form("%.5f", f_gentau_matchedSC_spptot_seed);
    			TString s_gentau_matchedSC_sppmax_seed = Form("%.5f", f_gentau_matchedSC_sppmax_seed);
    			TString s_gentau_matchedSC_szz_seed = Form("%.5f", f_gentau_matchedSC_szz_seed);
    			TString s_gentau_matchedSC_srrtot_seed = Form("%.5f", f_gentau_matchedSC_srrtot_seed);
    			TString s_gentau_matchedSC_srrmax_seed = Form("%.5f", f_gentau_matchedSC_srrmax_seed);
    			TString s_gentau_matchedSC_srrmean_seed = Form("%.5f", f_gentau_matchedSC_srrmean_seed);
    			TString s_gentau_matchedSC_emaxe_seed = Form("%.5f", f_gentau_matchedSC_emaxe_seed);
    			TString s_gentau_matchedSC_hoe_seed = Form("%.5f", f_gentau_matchedSC_hoe_seed);
    			TString s_gentau_matchedSC_meanz_seed = Form("%.5f", f_gentau_matchedSC_meanz_seed);
    			TString s_gentau_matchedSC_layer10_seed = Form("%.5f", f_gentau_matchedSC_layer10_seed);
    			TString s_gentau_matchedSC_layer50_seed = Form("%.5f", f_gentau_matchedSC_layer50_seed);
    			TString s_gentau_matchedSC_layer90_seed = Form("%.5f", f_gentau_matchedSC_layer90_seed);
    			TString s_gentau_matchedSC_ntc67_seed = Form("%.5f", f_gentau_matchedSC_ntc67_seed);
    			TString s_gentau_matchedSC_ntc90_seed = Form("%.5f", f_gentau_matchedSC_ntc90_seed);
    			TString s_gentau_matchedSC_bdteg_seed = Form("%.5f", f_gentau_matchedSC_bdteg_seed);
    			TString s_gentau_matchedSC_quality_seed = Form("%.1i", f_gentau_matchedSC_quality_seed);

				TString s_gentau_decayMode = Form("%.1i", f_gentau_decayMode);
				TString s_gentau_decayMode_joined = Form("%.1i", f_gentau_decayMode_joined);

				TString s_gentau_vis_pt = Form("%.5f", f_gentau_vis_pt);

                TString s_target = Form("%.5f", f_target);
	
				//8+1
				/*TString dumper_v1 = s_gentau_matchedSC_firstlayer_seed+","+s_gentau_matchedSC_maxlayer_seed+","+s_gentau_matchedSC_layer10_seed+","+s_gentau_matchedSC_layer50_seed+","+s_gentau_matchedSC_layer90_seed+","+s_gentau_matchedSC_showerlength_seed+","+s_gentau_matchedSC_hoe_seed+","+s_gentau_matchedSC_meanz_seed+","+s_target+","+s_gentau_vis_pt+","+s_gentau_matchedSC_pt_tot+","+s_gentau_matchedSC_eta_Eweighted+","+s_gentau_decayMode;
				out_file_v1<<dumper_v1<<'\n';*/

				//10+1
				TString dumper_v2 = s_gentau_matchedSC_pt_tot+","+s_gentau_matchedSC_eta_Eweighted+","+s_gentau_matchedSC_firstlayer_seed+","+s_gentau_matchedSC_maxlayer_seed+","+s_gentau_matchedSC_layer10_seed+","+s_gentau_matchedSC_layer50_seed+","+s_gentau_matchedSC_layer90_seed+","+s_gentau_matchedSC_showerlength_seed+","+s_gentau_matchedSC_hoe_seed+","+s_gentau_matchedSC_meanz_seed+","+s_target;
				out_file_v2<<dumper_v2<<'\n';
				
			}
		}
	}

	out_file_v2.close();

	return;
}
