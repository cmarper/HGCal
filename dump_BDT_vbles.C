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
	ofstream out_file_v2;
	ofstream out_file_v3;
	ofstream out_file_v4;
	ofstream out_file_v5;
	ofstream out_file_v6;
	ofstream out_file_v7;

	out_file_v1.open("/home/llr/cms/mperez/BDT_training/CMSSW_10_2_1/xgboost/HGCal_BDT/data/BDTv1_DM_vbles.csv");
	out_file_v2.open("/home/llr/cms/mperez/BDT_training/CMSSW_10_2_1/xgboost/HGCal_BDT/data/BDTv2_DM_vbles.csv");
	out_file_v3.open("/home/llr/cms/mperez/BDT_training/CMSSW_10_2_1/xgboost/HGCal_BDT/data/BDTv3_DM_vbles.csv");
	out_file_v4.open("/home/llr/cms/mperez/BDT_training/CMSSW_10_2_1/xgboost/HGCal_BDT/data/BDTv4_DM_vbles.csv");
	out_file_v5.open("/home/llr/cms/mperez/BDT_training/CMSSW_10_2_1/xgboost/HGCal_BDT/data/BDTv5_DM_vbles.csv");
	out_file_v6.open("/home/llr/cms/mperez/BDT_training/CMSSW_10_2_1/xgboost/HGCal_BDT/data/BDTv6_DM_vbles.csv");
	out_file_v7.open("/home/llr/cms/mperez/BDT_training/CMSSW_10_2_1/xgboost/HGCal_BDT/data/BDTv7_DM_vbles.csv");

	TChain * in_tree = new TChain("SortedTree");	
	in_tree->Add("/data_CMS/cms/mperez/HGCal_data/May19/sorted/NTuple_ZTT_PU0_sorted0p3.root");

	Long64_t nentries = in_tree->GetEntries();
	cout<<"nentries="<<in_tree->GetEntries()<<endl;
	if (nevents != -1) 
		nentries = nevents;

	// old branches used

	int _gentau_n;

	vector<bool> 	*_gentau_isMatched;

	vector<int> 	*_gentau_numberMatchedcl3d;

	vector<float>   *_gentau_pT1cl3d_pTfraction;
	vector<float>   *_gentau_pT2cl3d_pTfraction;
	vector<float>   *_gentau_pT3cl3d_pTfraction;

	vector<float>   *_gentau_pT1cl3d_BDTeg;
	vector<float>   *_gentau_pT2cl3d_BDTeg;
	vector<float>   *_gentau_pT3cl3d_BDTeg;

	vector<int>	*_gentau_pT1cl3d_coreshowerlength;
	vector<int>	 *_gentau_pT2cl3d_coreshowerlength;
	vector<int>	 *_gentau_pT3cl3d_coreshowerlength;

	vector<int>	 *_gentau_pT1cl3d_maxlayer;
	vector<int>	 *_gentau_pT2cl3d_maxlayer;
	vector<int>	 *_gentau_pT3cl3d_maxlayer;

	vector<int>   *_gentau_decayMode;

	in_tree->SetBranchAddress("gentau_n",	&_gentau_n);
	in_tree->SetBranchAddress("gentau_isMatched",	&_gentau_isMatched);

	in_tree->SetBranchAddress("gentau_numberMatchedcl3d",	&_gentau_numberMatchedcl3d);

	in_tree->SetBranchAddress("gentau_pT1cl3d_pTfraction",	&_gentau_pT1cl3d_pTfraction);
	in_tree->SetBranchAddress("gentau_pT2cl3d_pTfraction",	&_gentau_pT2cl3d_pTfraction);
	in_tree->SetBranchAddress("gentau_pT3cl3d_pTfraction",	&_gentau_pT3cl3d_pTfraction);

	in_tree->SetBranchAddress("gentau_pT1cl3d_BDTeg", 		&_gentau_pT1cl3d_BDTeg);
	in_tree->SetBranchAddress("gentau_pT2cl3d_BDTeg", 		&_gentau_pT2cl3d_BDTeg);
	in_tree->SetBranchAddress("gentau_pT3cl3d_BDTeg", 		&_gentau_pT3cl3d_BDTeg);

	in_tree->SetBranchAddress("gentau_pT1cl3d_coreshowerlength",			   &_gentau_pT1cl3d_coreshowerlength);
		in_tree->SetBranchAddress("gentau_pT2cl3d_coreshowerlength",			   &_gentau_pT2cl3d_coreshowerlength);
		in_tree->SetBranchAddress("gentau_pT3cl3d_coreshowerlength",			   &_gentau_pT3cl3d_coreshowerlength);

	in_tree->SetBranchAddress("gentau_pT1cl3d_maxlayer",			   &_gentau_pT1cl3d_maxlayer);
		in_tree->SetBranchAddress("gentau_pT2cl3d_maxlayer",			   &_gentau_pT2cl3d_maxlayer);
		in_tree->SetBranchAddress("gentau_pT3cl3d_maxlayer",			   &_gentau_pT3cl3d_maxlayer);

	in_tree->SetBranchAddress("gentau_decayMode",	&_gentau_decayMode);

	for (int i=0;i<nentries;i++) {

		if(i%1000==0) cout<<"i="<<i<<endl;

		//old branches

		_gentau_n = 0;

		_gentau_isMatched = 0;

		_gentau_numberMatchedcl3d = 0;

		_gentau_pT1cl3d_pTfraction = 0;
		_gentau_pT2cl3d_pTfraction = 0;
		_gentau_pT3cl3d_pTfraction = 0;

		_gentau_pT1cl3d_BDTeg = 0;
		_gentau_pT2cl3d_BDTeg = 0;
		_gentau_pT3cl3d_BDTeg = 0;

		_gentau_pT1cl3d_coreshowerlength = 0;
		_gentau_pT2cl3d_coreshowerlength = 0;
		_gentau_pT3cl3d_coreshowerlength = 0;

		_gentau_pT1cl3d_maxlayer = 0;
		_gentau_pT2cl3d_maxlayer = 0;
		_gentau_pT3cl3d_maxlayer = 0;

		_gentau_decayMode = 0;

		int entry_ok = in_tree->GetEntry(i);	
		if(entry_ok<0) 
			continue;

		for(int i_gentau=0; i_gentau<_gentau_n; i_gentau++){

			if (!(*_gentau_isMatched)[i_gentau]) continue;

			if ((*_gentau_decayMode)[i_gentau] == 0 || (*_gentau_decayMode)[i_gentau] == 1 || (*_gentau_decayMode)[i_gentau] == 4 || (*_gentau_decayMode)[i_gentau] == 5){

				int f_gentau_numberMatchedcl3d = (*_gentau_numberMatchedcl3d)[i_gentau];

				float f_gentau_pT1cl3d_pTfraction = (*_gentau_pT1cl3d_pTfraction)[i_gentau];
				float f_gentau_pT2cl3d_pTfraction = (*_gentau_pT2cl3d_pTfraction)[i_gentau];
				float f_gentau_pT3cl3d_pTfraction = (*_gentau_pT3cl3d_pTfraction)[i_gentau];

				float f_gentau_pT1cl3d_BDTeg = (*_gentau_pT1cl3d_BDTeg)[i_gentau];
				float f_gentau_pT2cl3d_BDTeg = (*_gentau_pT2cl3d_BDTeg)[i_gentau];
				float f_gentau_pT3cl3d_BDTeg = (*_gentau_pT3cl3d_BDTeg)[i_gentau];

				int f_gentau_pT1cl3d_coreshowerlength = (*_gentau_pT1cl3d_coreshowerlength)[i_gentau];
								int f_gentau_pT2cl3d_coreshowerlength = (*_gentau_pT2cl3d_coreshowerlength)[i_gentau];
								int f_gentau_pT3cl3d_coreshowerlength = (*_gentau_pT3cl3d_coreshowerlength)[i_gentau];

				int f_gentau_pT1cl3d_maxlayer = (*_gentau_pT1cl3d_maxlayer)[i_gentau];
								int f_gentau_pT2cl3d_maxlayer = (*_gentau_pT2cl3d_maxlayer)[i_gentau];
								int f_gentau_pT3cl3d_maxlayer = (*_gentau_pT3cl3d_maxlayer)[i_gentau];

				int f_gentau_decayMode = (*_gentau_decayMode)[i_gentau];

				int f_gentau_decayMode_joined;

				if (f_gentau_decayMode == 1 || f_gentau_decayMode == 5){
					f_gentau_decayMode_joined = 6;											
								}
				else f_gentau_decayMode_joined = f_gentau_decayMode;

				TString s_gentau_numberMatchedcl3d = Form("%.0i", f_gentau_numberMatchedcl3d);

				TString s_gentau_pT1cl3d_pTfraction = Form("%.5f", f_gentau_pT1cl3d_pTfraction);
				TString s_gentau_pT2cl3d_pTfraction = Form("%.5f", f_gentau_pT2cl3d_pTfraction);
				TString s_gentau_pT3cl3d_pTfraction = Form("%.5f", f_gentau_pT3cl3d_pTfraction);

				TString s_gentau_pT1cl3d_BDTeg = Form("%.5f", f_gentau_pT1cl3d_BDTeg);
				TString s_gentau_pT2cl3d_BDTeg = Form("%.5f", f_gentau_pT2cl3d_BDTeg);
				TString s_gentau_pT3cl3d_BDTeg = Form("%.5f", f_gentau_pT3cl3d_BDTeg);

				TString s_gentau_pT1cl3d_coreshowerlength = Form("%.1i", f_gentau_pT1cl3d_coreshowerlength);
								TString s_gentau_pT2cl3d_coreshowerlength = Form("%.1i", f_gentau_pT2cl3d_coreshowerlength);
								TString s_gentau_pT3cl3d_coreshowerlength = Form("%.1i", f_gentau_pT3cl3d_coreshowerlength);

				TString s_gentau_pT1cl3d_maxlayer = Form("%.1i", f_gentau_pT1cl3d_maxlayer);
								TString s_gentau_pT2cl3d_maxlayer = Form("%.1i", f_gentau_pT2cl3d_maxlayer);
								TString s_gentau_pT3cl3d_maxlayer = Form("%.1i", f_gentau_pT3cl3d_maxlayer);
				
				TString s_gentau_decayMode = Form("%.1i", f_gentau_decayMode);

				TString s_gentau_decayMode_joined = Form("%.1i", f_gentau_decayMode_joined);
	
				//7+1
				TString dumper_v1 = s_gentau_numberMatchedcl3d+","+s_gentau_pT1cl3d_pTfraction+","+s_gentau_pT2cl3d_pTfraction+","+s_gentau_pT3cl3d_pTfraction+","+s_gentau_pT1cl3d_BDTeg+","+s_gentau_pT2cl3d_BDTeg+","+s_gentau_pT3cl3d_BDTeg+","+s_gentau_decayMode;
				out_file_v1<<dumper_v1<<'\n';
				
				//13+1
				TString dumper_v2 = s_gentau_numberMatchedcl3d+","+s_gentau_pT1cl3d_pTfraction+","+s_gentau_pT2cl3d_pTfraction+","+s_gentau_pT3cl3d_pTfraction+","+s_gentau_pT1cl3d_BDTeg+","+s_gentau_pT2cl3d_BDTeg+","+s_gentau_pT3cl3d_BDTeg+","+s_gentau_pT1cl3d_coreshowerlength+","+s_gentau_pT2cl3d_coreshowerlength+","+s_gentau_pT3cl3d_coreshowerlength+","+s_gentau_pT1cl3d_maxlayer+","+s_gentau_pT2cl3d_maxlayer+","+s_gentau_pT3cl3d_maxlayer+","+s_gentau_decayMode;
								out_file_v2<<dumper_v2<<'\n';

				//5+1
				TString dumper_v3 = s_gentau_numberMatchedcl3d+","+s_gentau_pT1cl3d_pTfraction+","+s_gentau_pT2cl3d_pTfraction+","+s_gentau_pT1cl3d_BDTeg+","+s_gentau_pT2cl3d_BDTeg+","+s_gentau_decayMode;
								out_file_v3<<dumper_v3<<'\n';

				//9+1
				TString dumper_v4 = s_gentau_numberMatchedcl3d+","+s_gentau_pT1cl3d_pTfraction+","+s_gentau_pT2cl3d_pTfraction+","+s_gentau_pT1cl3d_BDTeg+","+s_gentau_pT2cl3d_BDTeg+","+s_gentau_pT1cl3d_coreshowerlength+","+s_gentau_pT2cl3d_coreshowerlength+","+s_gentau_pT1cl3d_maxlayer+","+s_gentau_pT2cl3d_maxlayer+","+s_gentau_decayMode;
								out_file_v4<<dumper_v4<<'\n';

				//8+1
				TString dumper_v5 = s_gentau_pT1cl3d_pTfraction+","+s_gentau_pT2cl3d_pTfraction+","+s_gentau_pT1cl3d_BDTeg+","+s_gentau_pT2cl3d_BDTeg+","+s_gentau_pT1cl3d_coreshowerlength+","+s_gentau_pT2cl3d_coreshowerlength+","+s_gentau_pT1cl3d_maxlayer+","+s_gentau_pT2cl3d_maxlayer+","+s_gentau_decayMode;
								out_file_v5<<dumper_v5<<'\n';

				//12+1
				TString dumper_v6 = s_gentau_pT1cl3d_pTfraction+","+s_gentau_pT2cl3d_pTfraction+","+s_gentau_pT3cl3d_pTfraction+","+s_gentau_pT1cl3d_BDTeg+","+s_gentau_pT2cl3d_BDTeg+","+s_gentau_pT3cl3d_BDTeg+","+s_gentau_pT1cl3d_coreshowerlength+","+s_gentau_pT2cl3d_coreshowerlength+","+s_gentau_pT3cl3d_coreshowerlength+","+s_gentau_pT1cl3d_maxlayer+","+s_gentau_pT2cl3d_maxlayer+","+s_gentau_pT3cl3d_maxlayer+","+s_gentau_decayMode;
								out_file_v6<<dumper_v6<<'\n';

				//12+1 (joined DM1 and DM5)
				TString dumper_v7 = s_gentau_pT1cl3d_pTfraction+","+s_gentau_pT2cl3d_pTfraction+","+s_gentau_pT3cl3d_pTfraction+","+s_gentau_pT1cl3d_BDTeg+","+s_gentau_pT2cl3d_BDTeg+","+s_gentau_pT3cl3d_BDTeg+","+s_gentau_pT1cl3d_coreshowerlength+","+s_gentau_pT2cl3d_coreshowerlength+","+s_gentau_pT3cl3d_coreshowerlength+","+s_gentau_pT1cl3d_maxlayer+","+s_gentau_pT2cl3d_maxlayer+","+s_gentau_pT3cl3d_maxlayer+","+s_gentau_decayMode_joined;
								out_file_v7<<dumper_v7<<'\n';

			}
		}
	}

	out_file_v1.close();
	out_file_v2.close();
	out_file_v3.close();
	out_file_v4.close();
	out_file_v5.close();
	out_file_v6.close();
	out_file_v7.close();

	return;
}



void dump_variables_regressionCalibration(int nevents = -1){
	
	ofstream out_file_v1;
	ofstream out_file_v2;
	ofstream out_file_v3;
	ofstream out_file_v4;

	out_file_v1.open("/home/llr/cms/mperez/BDT_training/CMSSW_10_2_1/xgboost/HGCal_BDT/regression_calib/data/BDTv1_calib_vbles.csv");
	out_file_v2.open("/home/llr/cms/mperez/BDT_training/CMSSW_10_2_1/xgboost/HGCal_BDT/regression_calib/data/BDTv2_calib_vbles.csv");
	out_file_v3.open("/home/llr/cms/mperez/BDT_training/CMSSW_10_2_1/xgboost/HGCal_BDT/regression_calib/data/BDTv3_calib_vbles.csv");
	out_file_v4.open("/home/llr/cms/mperez/BDT_training/CMSSW_10_2_1/xgboost/HGCal_BDT/regression_calib/data/BDTv4_calib_vbles.csv");

	TChain * in_tree = new TChain("SortedTree");	
	in_tree->Add("/data_CMS/cms/mperez/HGCal_data/May19/sorted/NTuple_ZTT_PU0_sorted0p3.root");

	Long64_t nentries = in_tree->GetEntries();
	cout<<"nentries="<<in_tree->GetEntries()<<endl;
	if (nevents != -1) 
		nentries = nevents;

	// old branches used

	int _gentau_n;

	vector<bool> 	*_gentau_isMatched;

	vector<int> 	*_gentau_numberMatchedcl3d;

	vector<float>   *_gentau_pT1cl3d_eta;

	vector<float>   *_gentau_pT1cl3d_pTfraction;
	vector<float>   *_gentau_pT2cl3d_pTfraction;
	vector<float>   *_gentau_pT3cl3d_pTfraction;

	vector<float>   *_gentau_pT1cl3d_BDTeg;
	vector<float>   *_gentau_pT2cl3d_BDTeg;
	vector<float>   *_gentau_pT3cl3d_BDTeg;

	vector<int>	 *_gentau_pT1cl3d_coreshowerlength;
	vector<int>	 *_gentau_pT2cl3d_coreshowerlength;
	vector<int>	 *_gentau_pT3cl3d_coreshowerlength;

	vector<int>	 *_gentau_pT1cl3d_maxlayer;
	vector<int>	 *_gentau_pT2cl3d_maxlayer;
	vector<int>	 *_gentau_pT3cl3d_maxlayer;

	vector<float>   *_gentau_TotPtMatchedcl3d;
	vector<float> 	*_gentau_vis_pt;
	vector<float>   *_gentau_eta;

	vector<int>   *_gentau_decayMode;


	in_tree->SetBranchAddress("gentau_n",	&_gentau_n);
	in_tree->SetBranchAddress("gentau_isMatched",	&_gentau_isMatched);

	in_tree->SetBranchAddress("gentau_numberMatchedcl3d",	&_gentau_numberMatchedcl3d);

	in_tree->SetBranchAddress("gentau_pT1cl3d_eta",	&_gentau_pT1cl3d_eta);

	in_tree->SetBranchAddress("gentau_pT1cl3d_pTfraction",	&_gentau_pT1cl3d_pTfraction);
	in_tree->SetBranchAddress("gentau_pT2cl3d_pTfraction",	&_gentau_pT2cl3d_pTfraction);
	in_tree->SetBranchAddress("gentau_pT3cl3d_pTfraction",	&_gentau_pT3cl3d_pTfraction);

	in_tree->SetBranchAddress("gentau_pT1cl3d_BDTeg", 		&_gentau_pT1cl3d_BDTeg);
	in_tree->SetBranchAddress("gentau_pT2cl3d_BDTeg", 		&_gentau_pT2cl3d_BDTeg);
	in_tree->SetBranchAddress("gentau_pT3cl3d_BDTeg", 		&_gentau_pT3cl3d_BDTeg);

	in_tree->SetBranchAddress("gentau_pT1cl3d_coreshowerlength",	&_gentau_pT1cl3d_coreshowerlength);
	in_tree->SetBranchAddress("gentau_pT2cl3d_coreshowerlength",	&_gentau_pT2cl3d_coreshowerlength);
	in_tree->SetBranchAddress("gentau_pT3cl3d_coreshowerlength",	&_gentau_pT3cl3d_coreshowerlength);

	in_tree->SetBranchAddress("gentau_pT1cl3d_maxlayer",   &_gentau_pT1cl3d_maxlayer);
	in_tree->SetBranchAddress("gentau_pT2cl3d_maxlayer",   &_gentau_pT2cl3d_maxlayer);
	in_tree->SetBranchAddress("gentau_pT3cl3d_maxlayer",   &_gentau_pT3cl3d_maxlayer);

	in_tree->SetBranchAddress("gentau_decayMode",	&_gentau_decayMode);

	in_tree->SetBranchAddress("gentau_TotPtMatchedcl3d", &_gentau_TotPtMatchedcl3d);
	in_tree->SetBranchAddress("gentau_vis_pt", &_gentau_vis_pt);
	in_tree->SetBranchAddress("gentau_eta", &_gentau_eta);

	for (int i=0;i<nentries;i++) {

		if(i%1000==0) cout<<"i="<<i<<endl;

		//old branches

		_gentau_n = 0;

		_gentau_isMatched = 0;

		_gentau_numberMatchedcl3d = 0;

		_gentau_pT1cl3d_eta = 0;

		_gentau_pT1cl3d_pTfraction = 0;
		_gentau_pT2cl3d_pTfraction = 0;
		_gentau_pT3cl3d_pTfraction = 0;

		_gentau_pT1cl3d_BDTeg = 0;
		_gentau_pT2cl3d_BDTeg = 0;
		_gentau_pT3cl3d_BDTeg = 0;

		_gentau_pT1cl3d_coreshowerlength = 0;
		_gentau_pT2cl3d_coreshowerlength = 0;
		_gentau_pT3cl3d_coreshowerlength = 0;

		_gentau_pT1cl3d_maxlayer = 0;
		_gentau_pT2cl3d_maxlayer = 0;
		_gentau_pT3cl3d_maxlayer = 0;

		_gentau_decayMode = 0;

		_gentau_TotPtMatchedcl3d = 0;
		_gentau_vis_pt = 0;		
		_gentau_eta = 0;

		int entry_ok = in_tree->GetEntry(i);	
		if(entry_ok<0) 
			continue;

		for(int i_gentau=0; i_gentau<_gentau_n; i_gentau++){

			if (!(*_gentau_isMatched)[i_gentau]) continue;

			if ((*_gentau_decayMode)[i_gentau] == 0 || (*_gentau_decayMode)[i_gentau] == 1 || (*_gentau_decayMode)[i_gentau] == 4 || (*_gentau_decayMode)[i_gentau] == 5){

				int f_gentau_numberMatchedcl3d = (*_gentau_numberMatchedcl3d)[i_gentau];

				float f_gentau_pT1cl3d_eta = (*_gentau_pT1cl3d_eta)[i_gentau];
				float f_gentau_pT1cl3d_eta_abs = abs(f_gentau_pT1cl3d_eta);

				float f_gentau_pT1cl3d_pTfraction = (*_gentau_pT1cl3d_pTfraction)[i_gentau];
				float f_gentau_pT2cl3d_pTfraction = (*_gentau_pT2cl3d_pTfraction)[i_gentau];
				float f_gentau_pT3cl3d_pTfraction = (*_gentau_pT3cl3d_pTfraction)[i_gentau];

				float f_gentau_pT1cl3d_BDTeg = (*_gentau_pT1cl3d_BDTeg)[i_gentau];
				float f_gentau_pT2cl3d_BDTeg = (*_gentau_pT2cl3d_BDTeg)[i_gentau];
				float f_gentau_pT3cl3d_BDTeg = (*_gentau_pT3cl3d_BDTeg)[i_gentau];

				int f_gentau_pT1cl3d_coreshowerlength = (*_gentau_pT1cl3d_coreshowerlength)[i_gentau];
				int f_gentau_pT2cl3d_coreshowerlength = (*_gentau_pT2cl3d_coreshowerlength)[i_gentau];
				int f_gentau_pT3cl3d_coreshowerlength = (*_gentau_pT3cl3d_coreshowerlength)[i_gentau];

				int f_gentau_pT1cl3d_maxlayer = (*_gentau_pT1cl3d_maxlayer)[i_gentau];
				int f_gentau_pT2cl3d_maxlayer = (*_gentau_pT2cl3d_maxlayer)[i_gentau];
				int f_gentau_pT3cl3d_maxlayer = (*_gentau_pT3cl3d_maxlayer)[i_gentau];

				int f_gentau_decayMode = (*_gentau_decayMode)[i_gentau];

				float f_gentau_TotPtMatchedcl3d = (*_gentau_TotPtMatchedcl3d)[i_gentau];
				float f_gentau_vis_pt = (*_gentau_vis_pt)[i_gentau];
				float f_target = f_gentau_vis_pt/f_gentau_TotPtMatchedcl3d;
				
				float f_gentau_eta = (*_gentau_eta)[i_gentau];
				
				TString s_gentau_numberMatchedcl3d = Form("%.0i", f_gentau_numberMatchedcl3d);

				TString s_gentau_pT1cl3d_eta = Form("%.5f", f_gentau_pT1cl3d_eta_abs);

				TString s_gentau_pT1cl3d_pTfraction = Form("%.5f", f_gentau_pT1cl3d_pTfraction);
				TString s_gentau_pT2cl3d_pTfraction = Form("%.5f", f_gentau_pT2cl3d_pTfraction);
				TString s_gentau_pT3cl3d_pTfraction = Form("%.5f", f_gentau_pT3cl3d_pTfraction);

				TString s_gentau_pT1cl3d_BDTeg = Form("%.5f", f_gentau_pT1cl3d_BDTeg);
				TString s_gentau_pT2cl3d_BDTeg = Form("%.5f", f_gentau_pT2cl3d_BDTeg);
				TString s_gentau_pT3cl3d_BDTeg = Form("%.5f", f_gentau_pT3cl3d_BDTeg);

				TString s_gentau_pT1cl3d_coreshowerlength = Form("%.1i", f_gentau_pT1cl3d_coreshowerlength);
				TString s_gentau_pT2cl3d_coreshowerlength = Form("%.1i", f_gentau_pT2cl3d_coreshowerlength);
				TString s_gentau_pT3cl3d_coreshowerlength = Form("%.1i", f_gentau_pT3cl3d_coreshowerlength);

				TString s_gentau_pT1cl3d_maxlayer = Form("%.1i", f_gentau_pT1cl3d_maxlayer);
				TString s_gentau_pT2cl3d_maxlayer = Form("%.1i", f_gentau_pT2cl3d_maxlayer);
				TString s_gentau_pT3cl3d_maxlayer = Form("%.1i", f_gentau_pT3cl3d_maxlayer);
				
				TString s_gentau_decayMode = Form("%.1i", f_gentau_decayMode);

				TString s_target = Form("%.5f", f_target);

				TString s_gentau_TotPtMatchedcl3d = Form("%.5f", f_gentau_TotPtMatchedcl3d);
                                TString s_gentau_vis_pt = Form("%.5f", f_gentau_vis_pt);
				TString s_gentau_eta = Form("%.5f", f_gentau_eta);

				
				//13+1
				TString dumper_v1 = s_gentau_pT1cl3d_eta+","+s_gentau_pT1cl3d_pTfraction+","+s_gentau_pT2cl3d_pTfraction+","+s_gentau_pT3cl3d_pTfraction+","+s_gentau_pT1cl3d_BDTeg+","+s_gentau_pT2cl3d_BDTeg+","+s_gentau_pT3cl3d_BDTeg+","+s_gentau_pT1cl3d_coreshowerlength+","+s_gentau_pT2cl3d_coreshowerlength+","+s_gentau_pT3cl3d_coreshowerlength+","+s_gentau_pT1cl3d_maxlayer+","+s_gentau_pT2cl3d_maxlayer+","+s_gentau_pT3cl3d_maxlayer+","+s_target+","+s_gentau_vis_pt+","+s_gentau_TotPtMatchedcl3d+","+s_gentau_eta+","+s_gentau_decayMode;
				out_file_v1<<dumper_v1<<'\n';

			
				//14+1
				TString dumper_v2 = s_gentau_TotPtMatchedcl3d+","+s_gentau_pT1cl3d_eta+","+s_gentau_pT1cl3d_pTfraction+","+s_gentau_pT2cl3d_pTfraction+","+s_gentau_pT3cl3d_pTfraction+","+s_gentau_pT1cl3d_BDTeg+","+s_gentau_pT2cl3d_BDTeg+","+s_gentau_pT3cl3d_BDTeg+","+s_gentau_pT1cl3d_coreshowerlength+","+s_gentau_pT2cl3d_coreshowerlength+","+s_gentau_pT3cl3d_coreshowerlength+","+s_gentau_pT1cl3d_maxlayer+","+s_gentau_pT2cl3d_maxlayer+","+s_gentau_pT3cl3d_maxlayer+","+s_target+","+s_gentau_vis_pt+","+s_gentau_TotPtMatchedcl3d+","+s_gentau_eta+","+s_gentau_decayMode;
                                out_file_v2<<dumper_v2<<'\n';

				//11+1
				TString dumper_v3 = s_gentau_TotPtMatchedcl3d+","+s_gentau_pT1cl3d_eta+","+s_gentau_pT1cl3d_BDTeg+","+s_gentau_pT2cl3d_BDTeg+","+s_gentau_pT3cl3d_BDTeg+","+s_gentau_pT1cl3d_coreshowerlength+","+s_gentau_pT2cl3d_coreshowerlength+","+s_gentau_pT3cl3d_coreshowerlength+","+s_gentau_pT1cl3d_maxlayer+","+s_gentau_pT2cl3d_maxlayer+","+s_gentau_pT3cl3d_maxlayer+","+s_target+","+s_gentau_vis_pt+","+s_gentau_TotPtMatchedcl3d+","+s_gentau_eta+","+s_gentau_decayMode;
				out_file_v3<<dumper_v3<<'\n';

				//6+1
                                TString dumper_v4 = s_gentau_TotPtMatchedcl3d+","+s_gentau_pT1cl3d_eta+","+s_gentau_pT1cl3d_pTfraction+","+s_gentau_pT1cl3d_BDTeg+","+s_gentau_pT1cl3d_coreshowerlength+","+s_gentau_pT1cl3d_maxlayer+","+s_target+","+s_gentau_vis_pt+","+s_gentau_TotPtMatchedcl3d+","+s_gentau_eta+","+s_gentau_decayMode;
				out_file_v4<<dumper_v4<<'\n';

			}
		}
	}

	out_file_v1.close();
	out_file_v2.close();
	out_file_v3.close();
	out_file_v4.close();

	return;
}

