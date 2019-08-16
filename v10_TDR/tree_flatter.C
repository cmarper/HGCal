#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TH1F.h>
#include <iostream>

using namespace std;

void flatten_tree( TString filein, TString fileout, int nevents = -1){
	
	TFile* out_file = TFile::Open(fileout);
  	/*if(out_file!=0){
      cout<<fileout<<" already exists, please delete it before converting again"<<endl;
      return;
	}*/

	out_file = TFile::Open(fileout,"RECREATE");

	TChain * in_tree = new TChain("ClusteredTree");	
	in_tree->Add(filein);

	Long64_t nentries = in_tree->GetEntries();
	cout<<"nentries="<<in_tree->GetEntries()<<endl;
	if (nevents != -1) nentries = nevents;

	// old branches used

	int _in_event;

	int _in_gentau_n;

	vector<float>  *_in_gentau_pt; 
    vector<float>  *_in_gentau_eta;
    vector<float>  *_in_gentau_phi;
    vector<float>  *_in_gentau_energy;
    vector<float>  *_in_gentau_mass;

    vector<float>  *_in_gentau_vis_pt; 
    vector<float>  *_in_gentau_vis_eta;
    vector<float>  *_in_gentau_vis_phi;
    vector<float>  *_in_gentau_vis_energy;
    vector<float>  *_in_gentau_vis_mass;

    vector<int>    *_in_gentau_decayMode;

    vector<bool>   *_in_gentau_isMatched;

    vector<int>    *_in_gentau_matchedSC_n_cl3d;
    vector<float>  *_in_gentau_matchedSC_pt_tot;
    vector<float>  *_in_gentau_matchedSC_pt_seed;
    vector<float>  *_in_gentau_matchedSC_eta_Eweighted;
    vector<float>  *_in_gentau_matchedSC_eta_seed;
    vector<float>  *_in_gentau_matchedSC_phi_Eweighted;
    vector<float>  *_in_gentau_matchedSC_phi_seed;

    vector<int>    *_in_gentau_matchedSC_showerlength_seed;
    vector<int>    *_in_gentau_matchedSC_coreshowerlength_seed;
    vector<int>    *_in_gentau_matchedSC_firstlayer_seed;
    vector<int>    *_in_gentau_matchedSC_maxlayer_seed;
    vector<float>  *_in_gentau_matchedSC_seetot_seed;
    vector<float>  *_in_gentau_matchedSC_seemax_seed;
    vector<float>  *_in_gentau_matchedSC_spptot_seed;
    vector<float>  *_in_gentau_matchedSC_sppmax_seed;
    vector<float>  *_in_gentau_matchedSC_szz_seed;
    vector<float>  *_in_gentau_matchedSC_srrtot_seed;
    vector<float>  *_in_gentau_matchedSC_srrmax_seed;
    vector<float>  *_in_gentau_matchedSC_srrmean_seed;
    vector<float>  *_in_gentau_matchedSC_emaxe_seed;
    vector<float>  *_in_gentau_matchedSC_hoe_seed;
    vector<float>  *_in_gentau_matchedSC_meanz_seed;
    vector<float>  *_in_gentau_matchedSC_layer10_seed;
    vector<float>  *_in_gentau_matchedSC_layer50_seed;
    vector<float>  *_in_gentau_matchedSC_layer90_seed;
    vector<float>  *_in_gentau_matchedSC_ntc67_seed;
    vector<float>  *_in_gentau_matchedSC_ntc90_seed;
    vector<float>  *_in_gentau_matchedSC_bdteg_seed;
    vector<int>    *_in_gentau_matchedSC_quality_seed;


	in_tree->SetBranchAddress("event",&_in_event);

    in_tree->SetBranchAddress("gentau_n",&_in_gentau_n);

    in_tree->SetBranchAddress("gentau_pt",      &_in_gentau_pt);
    in_tree->SetBranchAddress("gentau_eta",     &_in_gentau_eta);
    in_tree->SetBranchAddress("gentau_phi",     &_in_gentau_phi);
    in_tree->SetBranchAddress("gentau_energy",  &_in_gentau_energy);
    in_tree->SetBranchAddress("gentau_mass",    &_in_gentau_mass);

    in_tree->SetBranchAddress("gentau_vis_pt",      &_in_gentau_vis_pt);
    in_tree->SetBranchAddress("gentau_vis_eta",     &_in_gentau_vis_eta);
    in_tree->SetBranchAddress("gentau_vis_phi",     &_in_gentau_vis_phi);
    in_tree->SetBranchAddress("gentau_vis_energy",  &_in_gentau_vis_energy);
    in_tree->SetBranchAddress("gentau_vis_mass",    &_in_gentau_vis_mass);

    in_tree->SetBranchAddress("gentau_decayMode",   &_in_gentau_decayMode);

    in_tree->SetBranchAddress("gentau_isMatched",	&_in_gentau_isMatched);

    in_tree->SetBranchAddress("gentau_matchedSC_n_cl3d", 		&_in_gentau_matchedSC_n_cl3d);
    in_tree->SetBranchAddress("gentau_matchedSC_pt_tot", 		&_in_gentau_matchedSC_pt_tot);
    in_tree->SetBranchAddress("gentau_matchedSC_pt_seed",		&_in_gentau_matchedSC_pt_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_eta_Eweighted", &_in_gentau_matchedSC_eta_Eweighted);
    in_tree->SetBranchAddress("gentau_matchedSC_eta_seed", 		&_in_gentau_matchedSC_eta_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_phi_Eweighted", &_in_gentau_matchedSC_phi_Eweighted);
    in_tree->SetBranchAddress("gentau_matchedSC_phi_seed", 		&_in_gentau_matchedSC_phi_seed);

    in_tree->SetBranchAddress("gentau_matchedSC_showerlength_seed", 	&_in_gentau_matchedSC_showerlength_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_coreshowerlength_seed", &_in_gentau_matchedSC_coreshowerlength_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_firstlayer_seed", 		&_in_gentau_matchedSC_firstlayer_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_maxlayer_seed", 		&_in_gentau_matchedSC_maxlayer_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_seetot_seed", 			&_in_gentau_matchedSC_seetot_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_seemax_seed", 			&_in_gentau_matchedSC_seemax_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_spptot_seed", 			&_in_gentau_matchedSC_spptot_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_sppmax_seed", 			&_in_gentau_matchedSC_sppmax_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_szz_seed", 				&_in_gentau_matchedSC_szz_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_srrtot_seed", 			&_in_gentau_matchedSC_srrtot_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_srrmax_seed", 			&_in_gentau_matchedSC_srrmax_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_srrmean_seed", 			&_in_gentau_matchedSC_srrmean_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_emaxe_seed", 			&_in_gentau_matchedSC_emaxe_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_hoe_seed", 				&_in_gentau_matchedSC_hoe_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_meanz_seed", 			&_in_gentau_matchedSC_meanz_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_layer10_seed", 			&_in_gentau_matchedSC_layer10_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_layer50_seed", 			&_in_gentau_matchedSC_layer50_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_layer90_seed", 			&_in_gentau_matchedSC_layer90_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_ntc67_seed", 			&_in_gentau_matchedSC_ntc67_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_ntc90_seed", 			&_in_gentau_matchedSC_ntc90_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_bdteg_seed", 			&_in_gentau_matchedSC_bdteg_seed);
    in_tree->SetBranchAddress("gentau_matchedSC_quality_seed", 			&_in_gentau_matchedSC_quality_seed);

	// new branches

	//TTree* out_tree=in_tree->GetTree()->CloneTree(0);
	TTree* out_tree = new TTree;
	out_tree->SetNameTitle("FlatTree","FlatTree");

	int _out_event;

	int _out_gentau_n;

	float  _out_gentau_pt; 
    float  _out_gentau_eta;
    float  _out_gentau_phi;
    float  _out_gentau_energy;
    float  _out_gentau_mass;

    float  _out_gentau_vis_pt; 
    float  _out_gentau_vis_eta;
    float  _out_gentau_vis_phi;
    float  _out_gentau_vis_energy;
    float  _out_gentau_vis_mass;

    int    _out_gentau_decayMode;

    bool   _out_gentau_isMatched;

    int     _out_gentau_matchedSC_n_cl3d;
    float   _out_gentau_matchedSC_pt_tot;
    float   _out_gentau_matchedSC_pt_seed;
    float   _out_gentau_matchedSC_eta_Eweighted;
    float   _out_gentau_matchedSC_eta_seed;
    float   _out_gentau_matchedSC_phi_Eweighted;
    float   _out_gentau_matchedSC_phi_seed;

    int    _out_gentau_matchedSC_showerlength_seed;
    int    _out_gentau_matchedSC_coreshowerlength_seed;
    int    _out_gentau_matchedSC_firstlayer_seed;
    int    _out_gentau_matchedSC_maxlayer_seed;
    float  _out_gentau_matchedSC_seetot_seed;
    float  _out_gentau_matchedSC_seemax_seed;
    float  _out_gentau_matchedSC_spptot_seed;
    float  _out_gentau_matchedSC_sppmax_seed;
    float  _out_gentau_matchedSC_szz_seed;
    float  _out_gentau_matchedSC_srrtot_seed;
    float  _out_gentau_matchedSC_srrmax_seed;
    float  _out_gentau_matchedSC_srrmean_seed;
    float  _out_gentau_matchedSC_emaxe_seed;
    float  _out_gentau_matchedSC_hoe_seed;
    float  _out_gentau_matchedSC_meanz_seed;
    float  _out_gentau_matchedSC_layer10_seed;
    float  _out_gentau_matchedSC_layer50_seed;
    float  _out_gentau_matchedSC_layer90_seed;
    float  _out_gentau_matchedSC_ntc67_seed;
    float  _out_gentau_matchedSC_ntc90_seed;
    float  _out_gentau_matchedSC_bdteg_seed;
    int    _out_gentau_matchedSC_quality_seed;

	out_tree->Branch("event",&_out_event);

    out_tree->Branch("gentau_n",&_out_gentau_n);

    out_tree->Branch("gentau_pt",      &_out_gentau_pt);
    out_tree->Branch("gentau_eta",     &_out_gentau_eta);
    out_tree->Branch("gentau_phi",     &_out_gentau_phi);
    out_tree->Branch("gentau_energy",  &_out_gentau_energy);
    out_tree->Branch("gentau_mass",    &_out_gentau_mass);

    out_tree->Branch("gentau_vis_pt",      &_out_gentau_vis_pt);
    out_tree->Branch("gentau_vis_eta",     &_out_gentau_vis_eta);
    out_tree->Branch("gentau_vis_phi",     &_out_gentau_vis_phi);
    out_tree->Branch("gentau_vis_energy",  &_out_gentau_vis_energy);
    out_tree->Branch("gentau_vis_mass",    &_out_gentau_vis_mass);

    out_tree->Branch("gentau_decayMode",   &_out_gentau_decayMode);

    out_tree->Branch("gentau_isMatched",	&_out_gentau_isMatched);

    out_tree->Branch("gentau_matchedSC_n_cl3d", 		&_out_gentau_matchedSC_n_cl3d);
    out_tree->Branch("gentau_matchedSC_pt_tot", 		&_out_gentau_matchedSC_pt_tot);
    out_tree->Branch("gentau_matchedSC_pt_seed",		&_out_gentau_matchedSC_pt_seed);
    out_tree->Branch("gentau_matchedSC_eta_Eweighted", 	&_out_gentau_matchedSC_eta_Eweighted);
    out_tree->Branch("gentau_matchedSC_eta_seed", 		&_out_gentau_matchedSC_eta_seed);
    out_tree->Branch("gentau_matchedSC_phi_Eweighted", 	&_out_gentau_matchedSC_phi_Eweighted);
    out_tree->Branch("gentau_matchedSC_phi_seed", 		&_out_gentau_matchedSC_phi_seed);

    out_tree->Branch("gentau_matchedSC_showerlength_seed", 		&_out_gentau_matchedSC_showerlength_seed);
    out_tree->Branch("gentau_matchedSC_coreshowerlength_seed", 	&_out_gentau_matchedSC_coreshowerlength_seed);
    out_tree->Branch("gentau_matchedSC_firstlayer_seed", 		&_out_gentau_matchedSC_firstlayer_seed);
    out_tree->Branch("gentau_matchedSC_maxlayer_seed", 			&_out_gentau_matchedSC_maxlayer_seed);
    out_tree->Branch("gentau_matchedSC_seetot_seed", 			&_out_gentau_matchedSC_seetot_seed);
    out_tree->Branch("gentau_matchedSC_seemax_seed", 			&_out_gentau_matchedSC_seemax_seed);
    out_tree->Branch("gentau_matchedSC_spptot_seed", 			&_out_gentau_matchedSC_spptot_seed);
    out_tree->Branch("gentau_matchedSC_sppmax_seed", 			&_out_gentau_matchedSC_sppmax_seed);
    out_tree->Branch("gentau_matchedSC_szz_seed", 				&_out_gentau_matchedSC_szz_seed);
    out_tree->Branch("gentau_matchedSC_srrtot_seed", 			&_out_gentau_matchedSC_srrtot_seed);
    out_tree->Branch("gentau_matchedSC_srrmax_seed", 			&_out_gentau_matchedSC_srrmax_seed);
    out_tree->Branch("gentau_matchedSC_srrmean_seed", 			&_out_gentau_matchedSC_srrmean_seed);
    out_tree->Branch("gentau_matchedSC_emaxe_seed", 			&_out_gentau_matchedSC_emaxe_seed);
    out_tree->Branch("gentau_matchedSC_hoe_seed", 				&_out_gentau_matchedSC_hoe_seed);
    out_tree->Branch("gentau_matchedSC_meanz_seed", 			&_out_gentau_matchedSC_meanz_seed);
    out_tree->Branch("gentau_matchedSC_layer10_seed", 			&_out_gentau_matchedSC_layer10_seed);
    out_tree->Branch("gentau_matchedSC_layer50_seed", 			&_out_gentau_matchedSC_layer50_seed);
    out_tree->Branch("gentau_matchedSC_layer90_seed", 			&_out_gentau_matchedSC_layer90_seed);
    out_tree->Branch("gentau_matchedSC_ntc67_seed", 			&_out_gentau_matchedSC_ntc67_seed);
    out_tree->Branch("gentau_matchedSC_ntc90_seed", 			&_out_gentau_matchedSC_ntc90_seed);
    out_tree->Branch("gentau_matchedSC_bdteg_seed", 			&_out_gentau_matchedSC_bdteg_seed);
    out_tree->Branch("gentau_matchedSC_quality_seed", 			&_out_gentau_matchedSC_quality_seed);

	for (int i=0;i<nentries;i++) {

		if(i%1000==0) cout<<"i="<<i<<endl;

		// old branches

		_in_event = 0;

		_in_gentau_n = 0;

		_in_gentau_pt = 0; 
        _in_gentau_eta = 0;
        _in_gentau_phi = 0;
        _in_gentau_energy = 0;
        _in_gentau_mass = 0;

        _in_gentau_vis_pt = 0; 
        _in_gentau_vis_eta = 0; 
        _in_gentau_vis_phi = 0; 
        _in_gentau_vis_energy = 0; 
        _in_gentau_vis_mass = 0; 

        _in_gentau_decayMode = 0; 

    	_in_gentau_isMatched = 0;

    	_in_gentau_matchedSC_n_cl3d = 0;
    	_in_gentau_matchedSC_pt_tot = 0;
    	_in_gentau_matchedSC_pt_seed = 0;
    	_in_gentau_matchedSC_eta_Eweighted = 0;
    	_in_gentau_matchedSC_eta_seed = 0;
    	_in_gentau_matchedSC_phi_Eweighted = 0;
    	_in_gentau_matchedSC_phi_seed = 0;

    	_in_gentau_matchedSC_showerlength_seed = 0;
    	_in_gentau_matchedSC_coreshowerlength_seed = 0;
    	_in_gentau_matchedSC_firstlayer_seed = 0;
    	_in_gentau_matchedSC_maxlayer_seed = 0;
    	_in_gentau_matchedSC_seetot_seed = 0;
    	_in_gentau_matchedSC_seemax_seed = 0;
    	_in_gentau_matchedSC_spptot_seed = 0;
    	_in_gentau_matchedSC_sppmax_seed = 0;
    	_in_gentau_matchedSC_szz_seed = 0;
    	_in_gentau_matchedSC_srrtot_seed = 0;
    	_in_gentau_matchedSC_srrmax_seed = 0;
    	_in_gentau_matchedSC_srrmean_seed = 0;
    	_in_gentau_matchedSC_emaxe_seed = 0;
    	_in_gentau_matchedSC_hoe_seed = 0;
    	_in_gentau_matchedSC_meanz_seed = 0;
    	_in_gentau_matchedSC_layer10_seed = 0;
    	_in_gentau_matchedSC_layer50_seed = 0;
    	_in_gentau_matchedSC_layer90_seed = 0;
    	_in_gentau_matchedSC_ntc67_seed = 0;
    	_in_gentau_matchedSC_ntc90_seed = 0;
    	_in_gentau_matchedSC_bdteg_seed = 0;
    	_in_gentau_matchedSC_quality_seed = 0;

		// new branches

		_out_event = 0;

		_out_gentau_n = 0;

		_out_gentau_pt = 0; 
        _out_gentau_eta = 0;
        _out_gentau_phi = 0;
        _out_gentau_energy = 0;
        _out_gentau_mass = 0;

        _out_gentau_vis_pt = 0; 
        _out_gentau_vis_eta = 0; 
        _out_gentau_vis_phi = 0; 
        _out_gentau_vis_energy = 0; 
        _out_gentau_vis_mass = 0; 

        _out_gentau_decayMode = 0; 

    	_out_gentau_isMatched = 0;

    	_out_gentau_matchedSC_n_cl3d = 0;
    	_out_gentau_matchedSC_pt_tot = 0;
    	_out_gentau_matchedSC_pt_seed = 0;
    	_out_gentau_matchedSC_eta_Eweighted = 0;
    	_out_gentau_matchedSC_eta_seed = 0;
    	_out_gentau_matchedSC_phi_Eweighted = 0;
    	_out_gentau_matchedSC_phi_seed = 0;

    	_out_gentau_matchedSC_showerlength_seed = 0;
    	_out_gentau_matchedSC_coreshowerlength_seed = 0;
    	_out_gentau_matchedSC_firstlayer_seed = 0;
    	_out_gentau_matchedSC_maxlayer_seed = 0;
    	_out_gentau_matchedSC_seetot_seed = 0;
    	_out_gentau_matchedSC_seemax_seed = 0;
    	_out_gentau_matchedSC_spptot_seed = 0;
    	_out_gentau_matchedSC_sppmax_seed = 0;
    	_out_gentau_matchedSC_szz_seed = 0;
    	_out_gentau_matchedSC_srrtot_seed = 0;
    	_out_gentau_matchedSC_srrmax_seed = 0;
    	_out_gentau_matchedSC_srrmean_seed = 0;
    	_out_gentau_matchedSC_emaxe_seed = 0;
    	_out_gentau_matchedSC_hoe_seed = 0;
    	_out_gentau_matchedSC_meanz_seed = 0;
    	_out_gentau_matchedSC_layer10_seed = 0;
    	_out_gentau_matchedSC_layer50_seed = 0;
    	_out_gentau_matchedSC_layer90_seed = 0;
    	_out_gentau_matchedSC_ntc67_seed = 0;
    	_out_gentau_matchedSC_ntc90_seed = 0;
    	_out_gentau_matchedSC_bdteg_seed = 0;
    	_out_gentau_matchedSC_quality_seed = 0;

		int entry_ok = in_tree->GetEntry(i);	
		if(entry_ok<0) continue;
		
		for(int i_gentau=0; i_gentau<_in_gentau_n; i_gentau++){

        	_out_event = _in_event;

        	_out_gentau_n = _in_gentau_n;

        	_out_gentau_pt = _in_gentau_pt->at(i_gentau);
        	_out_gentau_eta = _in_gentau_eta->at(i_gentau);
        	_out_gentau_phi = _in_gentau_phi->at(i_gentau);
        	_out_gentau_energy = _in_gentau_energy->at(i_gentau);
        	_out_gentau_mass = _in_gentau_mass->at(i_gentau);

        	_out_gentau_vis_pt = _in_gentau_vis_pt->at(i_gentau);
        	_out_gentau_vis_eta = _in_gentau_vis_eta->at(i_gentau);
        	_out_gentau_vis_phi = _in_gentau_vis_phi->at(i_gentau);
        	_out_gentau_vis_energy = _in_gentau_vis_energy->at(i_gentau);
        	_out_gentau_vis_mass = _in_gentau_vis_mass->at(i_gentau);

        	_out_gentau_decayMode = _in_gentau_decayMode->at(i_gentau);

        	_out_gentau_isMatched = _in_gentau_isMatched->at(i_gentau);

        	_out_gentau_matchedSC_n_cl3d = _in_gentau_matchedSC_n_cl3d->at(i_gentau);
        	_out_gentau_matchedSC_pt_tot = _in_gentau_matchedSC_pt_tot->at(i_gentau);
        	_out_gentau_matchedSC_pt_seed = _in_gentau_matchedSC_pt_seed->at(i_gentau);
        	_out_gentau_matchedSC_eta_Eweighted = _in_gentau_matchedSC_eta_Eweighted->at(i_gentau);
        	_out_gentau_matchedSC_eta_seed = _in_gentau_matchedSC_eta_seed->at(i_gentau);
        	_out_gentau_matchedSC_phi_Eweighted = _in_gentau_matchedSC_phi_Eweighted->at(i_gentau);
        	_out_gentau_matchedSC_phi_seed = _in_gentau_matchedSC_phi_seed->at(i_gentau);

        	_out_gentau_matchedSC_showerlength_seed = _in_gentau_matchedSC_showerlength_seed->at(i_gentau);
        	_out_gentau_matchedSC_coreshowerlength_seed = _in_gentau_matchedSC_coreshowerlength_seed->at(i_gentau);
        	_out_gentau_matchedSC_firstlayer_seed = _in_gentau_matchedSC_firstlayer_seed->at(i_gentau);
        	_out_gentau_matchedSC_maxlayer_seed = _in_gentau_matchedSC_maxlayer_seed->at(i_gentau);
        	_out_gentau_matchedSC_seetot_seed = _in_gentau_matchedSC_seetot_seed->at(i_gentau);
        	_out_gentau_matchedSC_seemax_seed = _in_gentau_matchedSC_seemax_seed->at(i_gentau);
        	_out_gentau_matchedSC_spptot_seed = _in_gentau_matchedSC_spptot_seed->at(i_gentau);
        	_out_gentau_matchedSC_sppmax_seed = _in_gentau_matchedSC_sppmax_seed->at(i_gentau);
        	_out_gentau_matchedSC_szz_seed = _in_gentau_matchedSC_szz_seed->at(i_gentau);
        	_out_gentau_matchedSC_srrtot_seed = _in_gentau_matchedSC_srrtot_seed->at(i_gentau);
        	_out_gentau_matchedSC_srrmax_seed = _in_gentau_matchedSC_srrmax_seed->at(i_gentau);
        	_out_gentau_matchedSC_srrmean_seed = _in_gentau_matchedSC_srrmean_seed->at(i_gentau);
        	_out_gentau_matchedSC_emaxe_seed = _in_gentau_matchedSC_emaxe_seed->at(i_gentau);
        	_out_gentau_matchedSC_hoe_seed = _in_gentau_matchedSC_hoe_seed->at(i_gentau);
        	_out_gentau_matchedSC_meanz_seed = _in_gentau_matchedSC_meanz_seed->at(i_gentau);
        	_out_gentau_matchedSC_layer10_seed = _in_gentau_matchedSC_layer10_seed->at(i_gentau);
        	_out_gentau_matchedSC_layer50_seed = _in_gentau_matchedSC_layer50_seed->at(i_gentau);
        	_out_gentau_matchedSC_layer90_seed = _in_gentau_matchedSC_layer90_seed->at(i_gentau);
        	_out_gentau_matchedSC_ntc67_seed = _in_gentau_matchedSC_ntc67_seed->at(i_gentau);
        	_out_gentau_matchedSC_ntc90_seed = _in_gentau_matchedSC_ntc90_seed->at(i_gentau);
        	_out_gentau_matchedSC_bdteg_seed = _in_gentau_matchedSC_bdteg_seed->at(i_gentau);
        	_out_gentau_matchedSC_quality_seed = _in_gentau_matchedSC_quality_seed->at(i_gentau);

        	out_tree->Fill();

		}

	}

	out_file->cd();
    out_tree->Write();
    out_file->Close();

	return;

}


void test(int n_events = -1){

  TString infile = "/data_CMS/cms/mperez/HGCal_data/Aug19/ntuple_RelValDiTau_Pt20To100_Eta1p6To2p9_clustered.root";
  TString outfile = "/data_CMS/cms/mperez/HGCal_data/Aug19/ntuple_RelValDiTau_Pt20To100_Eta1p6To2p9_clustered_flat.root";

  //TString infile = "/data_CMS/cms/mperez/HGCal_data/Aug19/ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_clustered.root";
  //TString outfile = "/data_CMS/cms/mperez/HGCal_data/Aug19/ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_clustered_flat.root";

  flatten_tree(infile, outfile, n_events);

}