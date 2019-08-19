/////////////////////////////////////////////////////////
///// HGCal L1 taus, C. Martin Perez, LLR, Jul 2019 /////
/////////////////////////////////////////////////////////

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <iostream>

bool pT_comparison_pairs(pair<int,TLorentzVector> pair1, pair<int,TLorentzVector> pair2){

  return (pair1.second).Pt()>(pair2.second).Pt();

}

void cluster_tree( TString filein, TString fileout, int nevents = -1, float thr_seed = 4, float thr_sec = 2, float dEta = 0.2, float dPhi = 0.4){

    bool debug = false;

    TFile* out_file = TFile::Open(fileout);
    /*if(out_file!=0){
      cout<<fileout<<" already exists, please delete it before converting again"<<endl;
      return;
    }*/

    out_file = TFile::Open(fileout,"RECREATE");

    TChain * in_tree = new TChain("SkimmedTree");    
    in_tree->Add(filein);

    Long64_t nentries = in_tree->GetEntries();
    cout<<"nentries="<<in_tree->GetEntries()<<endl;
    if (nevents != -1) 
        nentries = nevents;

    // old branches used

    int _event;

    int _gentau_n;

    vector<float> *_gentau_pt; 
    vector<float> *_gentau_eta;
    vector<float> *_gentau_phi;
    vector<float> *_gentau_energy;
    vector<float> *_gentau_mass;

    vector<float> *_gentau_vis_pt; 
    vector<float> *_gentau_vis_eta;
    vector<float> *_gentau_vis_phi;
    vector<float> *_gentau_vis_energy;
    vector<float> *_gentau_vis_mass;

    vector<int>   *_gentau_decayMode;

    int _cl3d_n;

    vector<float>  *_cl3d_pt;
    vector<float>  *_cl3d_energy;
    vector<float>  *_cl3d_eta;
    vector<float>  *_cl3d_phi;

    vector<int>    *_cl3d_showerlength;
    vector<int>    *_cl3d_coreshowerlength;
    vector<int>    *_cl3d_firstlayer;
    vector<int>    *_cl3d_maxlayer;
    vector<float>  *_cl3d_seetot;
    vector<float>  *_cl3d_seemax;
    vector<float>  *_cl3d_spptot;
    vector<float>  *_cl3d_sppmax;
    vector<float>  *_cl3d_szz;
    vector<float>  *_cl3d_srrtot;
    vector<float>  *_cl3d_srrmax;
    vector<float>  *_cl3d_srrmean;
    vector<float>  *_cl3d_emaxe;
    vector<float>  *_cl3d_hoe;
    vector<float>  *_cl3d_meanz;
    vector<float>  *_cl3d_layer10;
    vector<float>  *_cl3d_layer50;
    vector<float>  *_cl3d_layer90;
    vector<float>  *_cl3d_ntc67;
    vector<float>  *_cl3d_ntc90;
    vector<float>  *_cl3d_bdteg;
    vector<int>    *_cl3d_quality;
    vector<float>  *_cl3d_puBDT;

    in_tree->SetBranchAddress("event",&_event);

    in_tree->SetBranchAddress("gentau_n",&_gentau_n);

    in_tree->SetBranchAddress("gentau_pt",      &_gentau_pt);
    in_tree->SetBranchAddress("gentau_eta",     &_gentau_eta);
    in_tree->SetBranchAddress("gentau_phi",     &_gentau_phi);
    in_tree->SetBranchAddress("gentau_energy",  &_gentau_energy);
    in_tree->SetBranchAddress("gentau_mass",    &_gentau_mass);

    in_tree->SetBranchAddress("gentau_vis_pt",      &_gentau_vis_pt);
    in_tree->SetBranchAddress("gentau_vis_eta",     &_gentau_vis_eta);
    in_tree->SetBranchAddress("gentau_vis_phi",     &_gentau_vis_phi);
    in_tree->SetBranchAddress("gentau_vis_energy",  &_gentau_vis_energy);
    in_tree->SetBranchAddress("gentau_vis_mass",    &_gentau_vis_mass);

    in_tree->SetBranchAddress("gentau_decayMode",   &_gentau_decayMode);

    in_tree->SetBranchAddress("cl3d_n",  &_cl3d_n);

    in_tree->SetBranchAddress("cl3d_pt",               &_cl3d_pt);
    in_tree->SetBranchAddress("cl3d_energy",           &_cl3d_energy);
    in_tree->SetBranchAddress("cl3d_eta",              &_cl3d_eta);
    in_tree->SetBranchAddress("cl3d_phi",              &_cl3d_phi);

    in_tree->SetBranchAddress("cl3d_showerlength",      &_cl3d_showerlength);
    in_tree->SetBranchAddress("cl3d_coreshowerlength",  &_cl3d_coreshowerlength);
    in_tree->SetBranchAddress("cl3d_firstlayer",        &_cl3d_firstlayer);
    in_tree->SetBranchAddress("cl3d_maxlayer",          &_cl3d_maxlayer);
    in_tree->SetBranchAddress("cl3d_seetot",    &_cl3d_seetot);
    in_tree->SetBranchAddress("cl3d_seemax",    &_cl3d_seemax);
    in_tree->SetBranchAddress("cl3d_spptot",    &_cl3d_spptot);
    in_tree->SetBranchAddress("cl3d_sppmax",    &_cl3d_sppmax);
    in_tree->SetBranchAddress("cl3d_szz",       &_cl3d_szz);
    in_tree->SetBranchAddress("cl3d_srrtot",    &_cl3d_srrtot);
    in_tree->SetBranchAddress("cl3d_srrmax",    &_cl3d_srrmax);
    in_tree->SetBranchAddress("cl3d_srrmean",   &_cl3d_srrmean);
    in_tree->SetBranchAddress("cl3d_emaxe",     &_cl3d_emaxe);
    in_tree->SetBranchAddress("cl3d_hoe",       &_cl3d_hoe);
    in_tree->SetBranchAddress("cl3d_meanz",     &_cl3d_meanz);
    in_tree->SetBranchAddress("cl3d_layer10",   &_cl3d_layer10);
    in_tree->SetBranchAddress("cl3d_layer50",   &_cl3d_layer50);
    in_tree->SetBranchAddress("cl3d_layer90",   &_cl3d_layer90);
    in_tree->SetBranchAddress("cl3d_ntc67",     &_cl3d_ntc67);
    in_tree->SetBranchAddress("cl3d_ntc90",     &_cl3d_ntc90);
    in_tree->SetBranchAddress("cl3d_bdteg",     &_cl3d_bdteg);
    in_tree->SetBranchAddress("cl3d_quality",   &_cl3d_quality);
    in_tree->SetBranchAddress("cl3d_puBDT",     &_cl3d_puBDT);

    TTree* out_tree=in_tree->GetTree()->CloneTree(0);
    out_tree->SetNameTitle("ClusteredTree","ClusteredTree");

    vector<bool>   _gentau_isMatched;

    vector<int>    _gentau_matchedSC_n_cl3d;
    vector<float>  _gentau_matchedSC_pt_tot;
    vector<float>  _gentau_matchedSC_pt_seed;
    vector<float>  _gentau_matchedSC_eta_Eweighted;
    vector<float>  _gentau_matchedSC_eta_seed;
    vector<float>  _gentau_matchedSC_phi_Eweighted;
    vector<float>  _gentau_matchedSC_phi_seed;

    vector<int>    _gentau_matchedSC_showerlength_seed;
    vector<int>    _gentau_matchedSC_coreshowerlength_seed;
    vector<int>    _gentau_matchedSC_firstlayer_seed;
    vector<int>    _gentau_matchedSC_maxlayer_seed;
    vector<float>  _gentau_matchedSC_seetot_seed;
    vector<float>  _gentau_matchedSC_seemax_seed;
    vector<float>  _gentau_matchedSC_spptot_seed;
    vector<float>  _gentau_matchedSC_sppmax_seed;
    vector<float>  _gentau_matchedSC_szz_seed;
    vector<float>  _gentau_matchedSC_srrtot_seed;
    vector<float>  _gentau_matchedSC_srrmax_seed;
    vector<float>  _gentau_matchedSC_srrmean_seed;
    vector<float>  _gentau_matchedSC_emaxe_seed;
    vector<float>  _gentau_matchedSC_hoe_seed;
    vector<float>  _gentau_matchedSC_meanz_seed;
    vector<float>  _gentau_matchedSC_layer10_seed;
    vector<float>  _gentau_matchedSC_layer50_seed;
    vector<float>  _gentau_matchedSC_layer90_seed;
    vector<float>  _gentau_matchedSC_ntc67_seed;
    vector<float>  _gentau_matchedSC_ntc90_seed;
    vector<float>  _gentau_matchedSC_bdteg_seed;
    vector<int>    _gentau_matchedSC_quality_seed;
    vector<float>  _gentau_matchedSC_puBDT_seed;

    int _n_supercluster;

    vector<int>    _supercluster_n_cl3d;
    vector<float>  _supercluster_pt_tot;
    vector<float>  _supercluster_pt_seed;
    vector<float>  _supercluster_eta_seed;
    vector<float>  _supercluster_eta_Eweighted;
    vector<float>  _supercluster_phi_seed;
    vector<float>  _supercluster_phi_Eweighted;

    vector<float> _supercluster_max_eta;
    vector<int>   _supercluster_max_eta_index;
    vector<float> _supercluster_min_eta;
    vector<int>   _supercluster_min_eta_index;
    vector<float> _supercluster_max_phi;
    vector<int>   _supercluster_max_phi_index;
    vector<float> _supercluster_min_phi;
    vector<int>   _supercluster_min_phi_index;

    vector<bool>  _supercluster_isMatched;
    vector<int>   _supercluster_matchedGenTauIdx;

    vector<vector<float>>  _supercluster_cl3d_pt;
    vector<vector<float>>  _supercluster_cl3d_energy;
    vector<vector<float>>  _supercluster_cl3d_eta;
    vector<vector<float>>  _supercluster_cl3d_phi;

    vector<vector<int>>    _supercluster_cl3d_showerlength;
    vector<vector<int>>    _supercluster_cl3d_coreshowerlength;
    vector<vector<int>>    _supercluster_cl3d_firstlayer;
    vector<vector<int>>    _supercluster_cl3d_maxlayer;
    vector<vector<float>>  _supercluster_cl3d_seetot;
    vector<vector<float>>  _supercluster_cl3d_seemax;
    vector<vector<float>>  _supercluster_cl3d_spptot;
    vector<vector<float>>  _supercluster_cl3d_sppmax;
    vector<vector<float>>  _supercluster_cl3d_szz;
    vector<vector<float>>  _supercluster_cl3d_srrtot;
    vector<vector<float>>  _supercluster_cl3d_srrmax;
    vector<vector<float>>  _supercluster_cl3d_srrmean;
    vector<vector<float>>  _supercluster_cl3d_emaxe;
    vector<vector<float>>  _supercluster_cl3d_hoe;
    vector<vector<float>>  _supercluster_cl3d_meanz;
    vector<vector<float>>  _supercluster_cl3d_layer10;
    vector<vector<float>>  _supercluster_cl3d_layer50;
    vector<vector<float>>  _supercluster_cl3d_layer90;
    vector<vector<float>>  _supercluster_cl3d_ntc67;
    vector<vector<float>>  _supercluster_cl3d_ntc90;
    vector<vector<float>>  _supercluster_cl3d_bdteg;
    vector<vector<int>>    _supercluster_cl3d_quality;
    vector<vector<float>>  _supercluster_cl3d_puBDT;

    out_tree->Branch("gentau_isMatched",  &_gentau_isMatched);

    out_tree->Branch("gentau_matchedSC_n_cl3d", &_gentau_matchedSC_n_cl3d);
    out_tree->Branch("gentau_matchedSC_pt_tot", &_gentau_matchedSC_pt_tot);
    out_tree->Branch("gentau_matchedSC_pt_seed", &_gentau_matchedSC_pt_seed);
    out_tree->Branch("gentau_matchedSC_eta_Eweighted", &_gentau_matchedSC_eta_Eweighted);
    out_tree->Branch("gentau_matchedSC_eta_seed", &_gentau_matchedSC_eta_seed);
    out_tree->Branch("gentau_matchedSC_phi_Eweighted", &_gentau_matchedSC_phi_Eweighted);
    out_tree->Branch("gentau_matchedSC_phi_seed", &_gentau_matchedSC_phi_seed);

    out_tree->Branch("gentau_matchedSC_showerlength_seed", &_gentau_matchedSC_showerlength_seed);
    out_tree->Branch("gentau_matchedSC_coreshowerlength_seed", &_gentau_matchedSC_coreshowerlength_seed);
    out_tree->Branch("gentau_matchedSC_firstlayer_seed", &_gentau_matchedSC_firstlayer_seed);
    out_tree->Branch("gentau_matchedSC_maxlayer_seed", &_gentau_matchedSC_maxlayer_seed);
    out_tree->Branch("gentau_matchedSC_seetot_seed", &_gentau_matchedSC_seetot_seed);
    out_tree->Branch("gentau_matchedSC_seemax_seed", &_gentau_matchedSC_seemax_seed);
    out_tree->Branch("gentau_matchedSC_spptot_seed", &_gentau_matchedSC_spptot_seed);
    out_tree->Branch("gentau_matchedSC_sppmax_seed", &_gentau_matchedSC_sppmax_seed);
    out_tree->Branch("gentau_matchedSC_szz_seed", &_gentau_matchedSC_szz_seed);
    out_tree->Branch("gentau_matchedSC_srrtot_seed", &_gentau_matchedSC_srrtot_seed);
    out_tree->Branch("gentau_matchedSC_srrmax_seed", &_gentau_matchedSC_srrmax_seed);
    out_tree->Branch("gentau_matchedSC_srrmean_seed", &_gentau_matchedSC_srrmean_seed);
    out_tree->Branch("gentau_matchedSC_emaxe_seed", &_gentau_matchedSC_emaxe_seed);
    out_tree->Branch("gentau_matchedSC_hoe_seed", &_gentau_matchedSC_hoe_seed);
    out_tree->Branch("gentau_matchedSC_meanz_seed", &_gentau_matchedSC_meanz_seed);
    out_tree->Branch("gentau_matchedSC_layer10_seed", &_gentau_matchedSC_layer10_seed);
    out_tree->Branch("gentau_matchedSC_layer50_seed", &_gentau_matchedSC_layer50_seed);
    out_tree->Branch("gentau_matchedSC_layer90_seed", &_gentau_matchedSC_layer90_seed);
    out_tree->Branch("gentau_matchedSC_ntc67_seed", &_gentau_matchedSC_ntc67_seed);
    out_tree->Branch("gentau_matchedSC_ntc90_seed", &_gentau_matchedSC_ntc90_seed);
    out_tree->Branch("gentau_matchedSC_bdteg_seed", &_gentau_matchedSC_bdteg_seed);
    out_tree->Branch("gentau_matchedSC_quality_seed", &_gentau_matchedSC_quality_seed);
    out_tree->Branch("gentau_matchedSC_puBDT_seed", &_gentau_matchedSC_puBDT_seed);

    out_tree->Branch("n_supercluster",     &_n_supercluster);

    out_tree->Branch("supercluster_n_cl3d",         &_supercluster_n_cl3d);
    out_tree->Branch("supercluster_pt_tot",         &_supercluster_pt_tot);
    out_tree->Branch("supercluster_pt_seed",        &_supercluster_pt_seed);
    out_tree->Branch("supercluster_eta_seed",       &_supercluster_eta_seed);
    out_tree->Branch("supercluster_eta_Eweighted",  &_supercluster_eta_Eweighted);
    out_tree->Branch("supercluster_phi_seedcl3d",   &_supercluster_phi_seed);
    out_tree->Branch("supercluster_phi_Eweighted",  &_supercluster_phi_Eweighted);

    out_tree->Branch("supercluster_max_eta",        &_supercluster_max_eta);
    out_tree->Branch("supercluster_max_eta_index",  &_supercluster_max_eta_index);
    out_tree->Branch("supercluster_min_eta",        &_supercluster_min_eta);
    out_tree->Branch("supercluster_min_eta_index",  &_supercluster_min_eta_index);
    out_tree->Branch("supercluster_max_phi",        &_supercluster_max_phi);
    out_tree->Branch("supercluster_max_phi_index",  &_supercluster_max_phi_index);
    out_tree->Branch("supercluster_min_phi",        &_supercluster_min_phi);
    out_tree->Branch("supercluster_min_phi_index",  &_supercluster_min_phi_index);

    out_tree->Branch("supercluster_isMatched",        &_supercluster_isMatched);
    out_tree->Branch("supercluster_matchedGenTauIdx", &_supercluster_matchedGenTauIdx);

    out_tree->Branch("supercluster_cl3d_pt",                 &_supercluster_cl3d_pt);
    out_tree->Branch("supercluster_cl3d_energy",             &_supercluster_cl3d_energy);
    out_tree->Branch("supercluster_cl3d_eta",                &_supercluster_cl3d_eta);
    out_tree->Branch("supercluster_cl3d_phi",                &_supercluster_cl3d_phi);

    out_tree->Branch("supercluster_cl3d_showerlength",      &_supercluster_cl3d_showerlength);
    out_tree->Branch("supercluster_cl3d_coreshowerlength",  &_supercluster_cl3d_coreshowerlength);
    out_tree->Branch("supercluster_cl3d_firstlayer", &_supercluster_cl3d_firstlayer);
    out_tree->Branch("supercluster_cl3d_maxlayer",   &_supercluster_cl3d_maxlayer);
    out_tree->Branch("supercluster_cl3d_seetot",     &_supercluster_cl3d_seetot);
    out_tree->Branch("supercluster_cl3d_seemax",     &_supercluster_cl3d_seemax);
    out_tree->Branch("supercluster_cl3d_spptot",     &_supercluster_cl3d_spptot);
    out_tree->Branch("supercluster_cl3d_sppmax",     &_supercluster_cl3d_sppmax);
    out_tree->Branch("supercluster_cl3d_szz",        &_supercluster_cl3d_szz);
    out_tree->Branch("supercluster_cl3d_srrtot",     &_supercluster_cl3d_srrtot);
    out_tree->Branch("supercluster_cl3d_srrmax",     &_supercluster_cl3d_srrmax);
    out_tree->Branch("supercluster_cl3d_srrmean",    &_supercluster_cl3d_srrmean);
    out_tree->Branch("supercluster_cl3d_emaxe",      &_supercluster_cl3d_emaxe);
    out_tree->Branch("supercluster_cl3d_hoe",        &_supercluster_cl3d_hoe);
    out_tree->Branch("supercluster_cl3d_meanz",      &_supercluster_cl3d_meanz);
    out_tree->Branch("supercluster_cl3d_layer10",    &_supercluster_cl3d_layer10);
    out_tree->Branch("supercluster_cl3d_layer50",    &_supercluster_cl3d_layer50);
    out_tree->Branch("supercluster_cl3d_layer90",    &_supercluster_cl3d_layer90);
    out_tree->Branch("supercluster_cl3d_ntc67",      &_supercluster_cl3d_ntc67);
    out_tree->Branch("supercluster_cl3d_ntc90",      &_supercluster_cl3d_ntc90);
    out_tree->Branch("supercluster_cl3d_bdteg",      &_supercluster_cl3d_bdteg);
    out_tree->Branch("supercluster_cl3d_quality",    &_supercluster_cl3d_quality);
    out_tree->Branch("supercluster_cl3d_puBDT",      &_supercluster_cl3d_puBDT);

    float matched_taus = 0;
    float unmatched_taus = 0;

    for (int i=0;i<nentries;i++) {

        if(i%1000==0) cout<<"i="<<i<<endl;

        // old branches

        _event = 0;

        _gentau_n = 0;

        _gentau_pt = 0; 
        _gentau_eta = 0;
        _gentau_phi = 0;
        _gentau_energy = 0;
        _gentau_mass = 0;

        _gentau_vis_pt = 0; 
        _gentau_vis_eta = 0; 
        _gentau_vis_phi = 0; 
        _gentau_vis_energy = 0; 
        _gentau_vis_mass = 0; 

        _gentau_decayMode = 0; 

        _cl3d_n = 0; 

        _cl3d_pt = 0; 
        _cl3d_energy = 0; 
        _cl3d_eta = 0; 
        _cl3d_phi = 0; 

        _cl3d_showerlength = 0;
        _cl3d_coreshowerlength = 0;
        _cl3d_firstlayer = 0;
        _cl3d_maxlayer = 0;      
        _cl3d_seetot = 0;
        _cl3d_seemax = 0;
        _cl3d_spptot = 0;
        _cl3d_sppmax = 0;
        _cl3d_szz = 0;
        _cl3d_srrtot = 0;
        _cl3d_srrmax = 0;
        _cl3d_srrmean = 0;
        _cl3d_emaxe = 0;
        _cl3d_hoe = 0;
        _cl3d_meanz = 0;
        _cl3d_layer10 = 0;
        _cl3d_layer50 = 0;
        _cl3d_layer90 = 0;
        _cl3d_ntc67 = 0;
        _cl3d_ntc90 = 0;
        _cl3d_bdteg = 0;
        _cl3d_quality = 0;
        _cl3d_puBDT = 0;

        _gentau_isMatched.clear();

        _gentau_matchedSC_n_cl3d.clear();
        _gentau_matchedSC_pt_tot.clear();
        _gentau_matchedSC_pt_seed.clear();
        _gentau_matchedSC_eta_Eweighted.clear();
        _gentau_matchedSC_eta_seed.clear();
        _gentau_matchedSC_phi_Eweighted.clear();
        _gentau_matchedSC_phi_seed.clear();

        _gentau_matchedSC_showerlength_seed.clear();
        _gentau_matchedSC_coreshowerlength_seed.clear();
        _gentau_matchedSC_firstlayer_seed.clear();
        _gentau_matchedSC_maxlayer_seed.clear();
        _gentau_matchedSC_seetot_seed.clear();
        _gentau_matchedSC_seemax_seed.clear();
        _gentau_matchedSC_spptot_seed.clear();
        _gentau_matchedSC_sppmax_seed.clear();
        _gentau_matchedSC_szz_seed.clear();
        _gentau_matchedSC_srrtot_seed.clear();
        _gentau_matchedSC_srrmax_seed.clear();
        _gentau_matchedSC_srrmean_seed.clear();
        _gentau_matchedSC_emaxe_seed.clear();
        _gentau_matchedSC_hoe_seed.clear();
        _gentau_matchedSC_meanz_seed.clear();
        _gentau_matchedSC_layer10_seed.clear();
        _gentau_matchedSC_layer50_seed.clear();
        _gentau_matchedSC_layer90_seed.clear();
        _gentau_matchedSC_ntc67_seed.clear();
        _gentau_matchedSC_ntc90_seed.clear();
        _gentau_matchedSC_bdteg_seed.clear();
        _gentau_matchedSC_quality_seed.clear();
        _gentau_matchedSC_puBDT_seed.clear();

        _n_supercluster = 0;

        _supercluster_n_cl3d.clear(); 
        _supercluster_pt_tot.clear(); 
        _supercluster_pt_seed.clear(); 
        _supercluster_eta_seed.clear(); 
        _supercluster_eta_Eweighted.clear(); 
        _supercluster_phi_seed.clear(); 
        _supercluster_phi_Eweighted.clear(); 

        _supercluster_max_eta.clear(); 
        _supercluster_max_eta_index.clear(); 
        _supercluster_min_eta.clear(); 
        _supercluster_min_eta_index.clear(); 
        _supercluster_max_phi.clear(); 
        _supercluster_max_phi_index.clear(); 
        _supercluster_min_phi.clear(); 
        _supercluster_min_phi_index.clear(); 

        _supercluster_isMatched.clear();
        _supercluster_matchedGenTauIdx.clear();

        _supercluster_cl3d_pt.clear();
        _supercluster_cl3d_energy.clear();
        _supercluster_cl3d_eta.clear();
        _supercluster_cl3d_phi.clear();

        _supercluster_cl3d_showerlength.clear();
        _supercluster_cl3d_coreshowerlength.clear();
        _supercluster_cl3d_firstlayer.clear();
        _supercluster_cl3d_maxlayer.clear();
        _supercluster_cl3d_seetot.clear();
        _supercluster_cl3d_seemax.clear();
        _supercluster_cl3d_spptot.clear();
        _supercluster_cl3d_sppmax.clear();
        _supercluster_cl3d_szz.clear();
        _supercluster_cl3d_srrtot.clear();
        _supercluster_cl3d_srrmax.clear();
        _supercluster_cl3d_srrmean.clear();
        _supercluster_cl3d_emaxe.clear();
        _supercluster_cl3d_hoe.clear();
        _supercluster_cl3d_meanz.clear();
        _supercluster_cl3d_layer10.clear();
        _supercluster_cl3d_layer50.clear();
        _supercluster_cl3d_layer90.clear();
        _supercluster_cl3d_ntc67.clear();
        _supercluster_cl3d_ntc90.clear();
        _supercluster_cl3d_bdteg.clear();
        _supercluster_cl3d_quality.clear();
        _supercluster_cl3d_puBDT.clear();

        int entry_ok = in_tree->GetEntry(i);    
        if(entry_ok<0) 
            continue;

        //if(_cl3d_n==0) cout<<"Event "<<_event<<", with n_gentau"<<_gentau_n<<", n_cl3d "<<_cl3d_n<<endl;

        float previous_pt_seeding_supercl = -999.;


        //////////////////////////////////////////
        ////////////// CLUSTERING ////////////////
        //////////////////////////////////////////

        // SUPERCLUSTERS

        if(debug){
            cout<<" "<<endl;
            cout<<" ==================> Event "<<i<< "  <================== "<<endl;
            cout<<" "<<endl;
            cout<<" Number of reconstructed cl3d: "<<_cl3d_n<<endl;
        }

        for (int i_main=0; i_main<_cl3d_n; i_main++){

            // RAW CLUSTERING 

            vector<TLorentzVector>  candidate_supercluster;
            vector<int>   candidate_supercluster_n_cl3d;
            vector<int>   candidate_supercluster_showerlength;
            vector<int>   candidate_supercluster_coreshowerlength;
            vector<int>   candidate_supercluster_firstlayer;
            vector<int>   candidate_supercluster_maxlayer;
            vector<float> candidate_supercluster_seetot;
            vector<float> candidate_supercluster_seemax;
            vector<float> candidate_supercluster_spptot;
            vector<float> candidate_supercluster_sppmax;
            vector<float> candidate_supercluster_szz;
            vector<float> candidate_supercluster_srrtot;
            vector<float> candidate_supercluster_srrmax;
            vector<float> candidate_supercluster_srrmean;
            vector<float> candidate_supercluster_emaxe;
            vector<float> candidate_supercluster_hoe;
            vector<float> candidate_supercluster_meanz;
            vector<float> candidate_supercluster_layer10;
            vector<float> candidate_supercluster_layer50;
            vector<float> candidate_supercluster_layer90;
            vector<float> candidate_supercluster_ntc67;
            vector<float> candidate_supercluster_ntc90;
            vector<float> candidate_supercluster_bdteg;
            vector<int>   candidate_supercluster_quality;
            vector<float> candidate_supercluster_puBDT;

            candidate_supercluster.clear();
            candidate_supercluster_n_cl3d.clear();
            candidate_supercluster_showerlength.clear();
            candidate_supercluster_coreshowerlength.clear();
            candidate_supercluster_firstlayer.clear();
            candidate_supercluster_maxlayer.clear();
            candidate_supercluster_seetot.clear();
            candidate_supercluster_seemax.clear();
            candidate_supercluster_spptot.clear();
            candidate_supercluster_sppmax.clear();
            candidate_supercluster_szz.clear();
            candidate_supercluster_srrtot.clear();
            candidate_supercluster_srrmax.clear();
            candidate_supercluster_srrmean.clear();
            candidate_supercluster_emaxe.clear();
            candidate_supercluster_hoe.clear();
            candidate_supercluster_meanz.clear();
            candidate_supercluster_layer10.clear();
            candidate_supercluster_layer50.clear();
            candidate_supercluster_layer90.clear();
            candidate_supercluster_ntc67.clear();
            candidate_supercluster_ntc90.clear();
            candidate_supercluster_bdteg.clear();
            candidate_supercluster_quality.clear();
            candidate_supercluster_puBDT.clear();

            // MAIN CLUSTER CANDIDATE (SUPERCLUSTER SEED)

            TLorentzVector seed_candidate_supercluster;
            int   seed_candidate_supercluster_showerlength;
            int   seed_candidate_supercluster_coreshowerlength;
            int   seed_candidate_supercluster_firstlayer;
            int   seed_candidate_supercluster_maxlayer;
            float seed_candidate_supercluster_seetot;
            float seed_candidate_supercluster_seemax;
            float seed_candidate_supercluster_spptot;
            float seed_candidate_supercluster_sppmax;
            float seed_candidate_supercluster_szz;
            float seed_candidate_supercluster_srrtot;
            float seed_candidate_supercluster_srrmax;
            float seed_candidate_supercluster_srrmean;
            float seed_candidate_supercluster_emaxe;
            float seed_candidate_supercluster_hoe;
            float seed_candidate_supercluster_meanz;
            float seed_candidate_supercluster_layer10;
            float seed_candidate_supercluster_layer50;
            float seed_candidate_supercluster_layer90;
            float seed_candidate_supercluster_ntc67;
            float seed_candidate_supercluster_ntc90;
            float seed_candidate_supercluster_bdteg;
            int   seed_candidate_supercluster_quality;
            float seed_candidate_supercluster_puBDT;

            seed_candidate_supercluster.SetPtEtaPhiM( (*_cl3d_pt)[i_main], (*_cl3d_eta)[i_main], (*_cl3d_phi)[i_main], 0);
            seed_candidate_supercluster_showerlength = (*_cl3d_showerlength)[i_main];
            seed_candidate_supercluster_coreshowerlength = (*_cl3d_coreshowerlength)[i_main];
            seed_candidate_supercluster_firstlayer = (*_cl3d_firstlayer)[i_main];
            seed_candidate_supercluster_maxlayer = (*_cl3d_maxlayer)[i_main];
            seed_candidate_supercluster_seetot = (*_cl3d_seetot)[i_main];
            seed_candidate_supercluster_seemax = (*_cl3d_seemax)[i_main];
            seed_candidate_supercluster_spptot = (*_cl3d_spptot)[i_main];
            seed_candidate_supercluster_sppmax = (*_cl3d_sppmax)[i_main];
            seed_candidate_supercluster_szz = (*_cl3d_szz)[i_main];
            seed_candidate_supercluster_srrtot = (*_cl3d_srrtot)[i_main];
            seed_candidate_supercluster_srrmax = (*_cl3d_srrmax)[i_main];
            seed_candidate_supercluster_srrmean = (*_cl3d_srrmean)[i_main];
            seed_candidate_supercluster_emaxe = (*_cl3d_emaxe)[i_main];
            seed_candidate_supercluster_hoe = (*_cl3d_hoe)[i_main];
            seed_candidate_supercluster_meanz = (*_cl3d_meanz)[i_main];
            seed_candidate_supercluster_layer10 = (*_cl3d_layer10)[i_main];
            seed_candidate_supercluster_layer50 = (*_cl3d_layer50)[i_main];
            seed_candidate_supercluster_layer90 = (*_cl3d_layer90)[i_main];
            seed_candidate_supercluster_ntc67 = (*_cl3d_ntc67)[i_main];
            seed_candidate_supercluster_ntc90 = (*_cl3d_ntc90)[i_main];
            seed_candidate_supercluster_bdteg = (*_cl3d_bdteg)[i_main];
            seed_candidate_supercluster_quality = (*_cl3d_quality)[i_main];
            seed_candidate_supercluster_puBDT = (*_cl3d_puBDT)[i_main];

            if ( seed_candidate_supercluster_puBDT < 0 ) continue;

            if ( seed_candidate_supercluster.Pt() < thr_seed ) continue;
        
            if(debug) cout<<"   Seed cl3d candidate #"<<i_main<<": pT "<<seed_candidate_supercluster.Pt()<<", eta "<<seed_candidate_supercluster.Eta()<<", phi "<<seed_candidate_supercluster.Phi()<<endl;

            candidate_supercluster.push_back(seed_candidate_supercluster);
            candidate_supercluster_showerlength.push_back(seed_candidate_supercluster_showerlength);
            candidate_supercluster_coreshowerlength.push_back(seed_candidate_supercluster_coreshowerlength);
            candidate_supercluster_firstlayer.push_back(seed_candidate_supercluster_firstlayer);
            candidate_supercluster_maxlayer.push_back(seed_candidate_supercluster_maxlayer);
            candidate_supercluster_seetot.push_back(seed_candidate_supercluster_seetot);
            candidate_supercluster_seemax.push_back(seed_candidate_supercluster_seemax);
            candidate_supercluster_spptot.push_back(seed_candidate_supercluster_spptot);
            candidate_supercluster_sppmax.push_back(seed_candidate_supercluster_sppmax);
            candidate_supercluster_szz.push_back(seed_candidate_supercluster_szz);
            candidate_supercluster_srrtot.push_back(seed_candidate_supercluster_srrtot);
            candidate_supercluster_srrmax.push_back(seed_candidate_supercluster_srrmax);
            candidate_supercluster_srrmean.push_back(seed_candidate_supercluster_srrmean);
            candidate_supercluster_emaxe.push_back(seed_candidate_supercluster_emaxe);
            candidate_supercluster_hoe.push_back(seed_candidate_supercluster_hoe);
            candidate_supercluster_meanz.push_back(seed_candidate_supercluster_meanz);
            candidate_supercluster_layer10.push_back(seed_candidate_supercluster_layer10);
            candidate_supercluster_layer50.push_back(seed_candidate_supercluster_layer50);
            candidate_supercluster_layer90.push_back(seed_candidate_supercluster_layer90);
            candidate_supercluster_ntc67.push_back(seed_candidate_supercluster_ntc67);
            candidate_supercluster_ntc90.push_back(seed_candidate_supercluster_ntc90);
            candidate_supercluster_bdteg.push_back(seed_candidate_supercluster_bdteg);
            candidate_supercluster_quality.push_back(seed_candidate_supercluster_quality);
            candidate_supercluster_puBDT.push_back(seed_candidate_supercluster_puBDT);

            // SECONDARY CLUSTERS CANDIDATES

            for (int i_sec=0; i_sec<_cl3d_n; i_sec++){

                if( i_sec==i_main ) continue;

                TLorentzVector secondary_candidate_supercluster;
                int   secondary_candidate_supercluster_showerlength;
                int   secondary_candidate_supercluster_coreshowerlength;
                int   secondary_candidate_supercluster_firstlayer;
                int   secondary_candidate_supercluster_maxlayer;
                float secondary_candidate_supercluster_seetot;
                float secondary_candidate_supercluster_seemax;
                float secondary_candidate_supercluster_spptot;
                float secondary_candidate_supercluster_sppmax;
                float secondary_candidate_supercluster_szz;
                float secondary_candidate_supercluster_srrtot;
                float secondary_candidate_supercluster_srrmax;
                float secondary_candidate_supercluster_srrmean;
                float secondary_candidate_supercluster_emaxe;
                float secondary_candidate_supercluster_hoe;
                float secondary_candidate_supercluster_meanz;
                float secondary_candidate_supercluster_layer10;
                float secondary_candidate_supercluster_layer50;
                float secondary_candidate_supercluster_layer90;
                float secondary_candidate_supercluster_ntc67;
                float secondary_candidate_supercluster_ntc90;
                float secondary_candidate_supercluster_bdteg;
                int   secondary_candidate_supercluster_quality;
                float secondary_candidate_supercluster_puBDT;

                secondary_candidate_supercluster.SetPtEtaPhiM( (*_cl3d_pt)[i_sec], (*_cl3d_eta)[i_sec], (*_cl3d_phi)[i_sec], 0);
                secondary_candidate_supercluster_showerlength = (*_cl3d_showerlength)[i_sec];
                secondary_candidate_supercluster_coreshowerlength = (*_cl3d_coreshowerlength)[i_sec];
                secondary_candidate_supercluster_firstlayer = (*_cl3d_firstlayer)[i_sec];
                secondary_candidate_supercluster_maxlayer = (*_cl3d_maxlayer)[i_sec];
                secondary_candidate_supercluster_seetot = (*_cl3d_seetot)[i_sec];
                secondary_candidate_supercluster_seemax = (*_cl3d_seemax)[i_sec];
                secondary_candidate_supercluster_spptot = (*_cl3d_spptot)[i_sec];
                secondary_candidate_supercluster_sppmax = (*_cl3d_sppmax)[i_sec];
                secondary_candidate_supercluster_szz = (*_cl3d_szz)[i_sec];
                secondary_candidate_supercluster_srrtot = (*_cl3d_srrtot)[i_sec];
                secondary_candidate_supercluster_srrmax = (*_cl3d_srrmax)[i_sec];
                secondary_candidate_supercluster_srrmean = (*_cl3d_srrmean)[i_sec];
                secondary_candidate_supercluster_emaxe = (*_cl3d_emaxe)[i_sec];
                secondary_candidate_supercluster_hoe = (*_cl3d_hoe)[i_sec];
                secondary_candidate_supercluster_meanz = (*_cl3d_meanz)[i_sec];
                secondary_candidate_supercluster_layer10 = (*_cl3d_layer10)[i_sec];
                secondary_candidate_supercluster_layer50 = (*_cl3d_layer50)[i_sec];
                secondary_candidate_supercluster_layer90 = (*_cl3d_layer90)[i_sec];
                secondary_candidate_supercluster_ntc67 = (*_cl3d_ntc67)[i_sec];
                secondary_candidate_supercluster_ntc90 = (*_cl3d_ntc90)[i_sec];
                secondary_candidate_supercluster_bdteg = (*_cl3d_bdteg)[i_sec];
                secondary_candidate_supercluster_quality = (*_cl3d_quality)[i_sec];
                secondary_candidate_supercluster_puBDT = (*_cl3d_puBDT)[i_sec];

                if( secondary_candidate_supercluster_puBDT < 0 ) continue;

                if( secondary_candidate_supercluster.Pt() < thr_sec ) continue;
 
                if( fabs(secondary_candidate_supercluster.Eta() - seed_candidate_supercluster.Eta()) > dEta ) continue; // (TT/2)*(eta/TT) -> 3
                if( fabs(seed_candidate_supercluster.DeltaPhi(secondary_candidate_supercluster)) > dPhi ) continue; // => 9

                if(debug) cout<<"      Secondary cl3d candidate #"<<i_sec<<": pT "<<secondary_candidate_supercluster.Pt()<<", eta "<<secondary_candidate_supercluster.Eta()<<", phi "<<secondary_candidate_supercluster.Phi()<<endl;

                candidate_supercluster.push_back(secondary_candidate_supercluster);
                candidate_supercluster_showerlength.push_back(secondary_candidate_supercluster_showerlength);
                candidate_supercluster_coreshowerlength.push_back(secondary_candidate_supercluster_coreshowerlength);
                candidate_supercluster_firstlayer.push_back(secondary_candidate_supercluster_firstlayer);
                candidate_supercluster_maxlayer.push_back(secondary_candidate_supercluster_maxlayer);
                candidate_supercluster_seetot.push_back(secondary_candidate_supercluster_seetot);
                candidate_supercluster_seemax.push_back(secondary_candidate_supercluster_seemax);
                candidate_supercluster_spptot.push_back(secondary_candidate_supercluster_spptot);
                candidate_supercluster_sppmax.push_back(secondary_candidate_supercluster_sppmax);
                candidate_supercluster_szz.push_back(secondary_candidate_supercluster_szz);
                candidate_supercluster_srrtot.push_back(secondary_candidate_supercluster_srrtot);
                candidate_supercluster_srrmax.push_back(secondary_candidate_supercluster_srrmax);
                candidate_supercluster_srrmean.push_back(secondary_candidate_supercluster_srrmean);
                candidate_supercluster_emaxe.push_back(secondary_candidate_supercluster_emaxe);
                candidate_supercluster_hoe.push_back(secondary_candidate_supercluster_hoe);
                candidate_supercluster_meanz.push_back(secondary_candidate_supercluster_meanz);
                candidate_supercluster_layer10.push_back(secondary_candidate_supercluster_layer10);
                candidate_supercluster_layer50.push_back(secondary_candidate_supercluster_layer50);
                candidate_supercluster_layer90.push_back(secondary_candidate_supercluster_layer90);
                candidate_supercluster_ntc67.push_back(secondary_candidate_supercluster_ntc67);
                candidate_supercluster_ntc90.push_back(secondary_candidate_supercluster_ntc90);
                candidate_supercluster_bdteg.push_back(secondary_candidate_supercluster_bdteg);
                candidate_supercluster_quality.push_back(secondary_candidate_supercluster_quality);
                candidate_supercluster_puBDT.push_back(secondary_candidate_supercluster_puBDT);

            }

            // SORT CLUSTERS BY PT IN SUPERCLUSTER

            vector<TLorentzVector>  candidate_supercluster_pTsorted;
            vector<float> candidate_supercluster_pt_pTsorted;
            vector<float> candidate_supercluster_eta_pTsorted;
            vector<float> candidate_supercluster_phi_pTsorted;
            vector<float> candidate_supercluster_E_pTsorted;
            vector<int>   candidate_supercluster_showerlength_pTsorted;
            vector<int>   candidate_supercluster_coreshowerlength_pTsorted;
            vector<int>   candidate_supercluster_firstlayer_pTsorted;
            vector<int>   candidate_supercluster_maxlayer_pTsorted;
            vector<float> candidate_supercluster_seetot_pTsorted;
            vector<float> candidate_supercluster_seemax_pTsorted;
            vector<float> candidate_supercluster_spptot_pTsorted;
            vector<float> candidate_supercluster_sppmax_pTsorted;
            vector<float> candidate_supercluster_szz_pTsorted;
            vector<float> candidate_supercluster_srrtot_pTsorted;
            vector<float> candidate_supercluster_srrmax_pTsorted;
            vector<float> candidate_supercluster_srrmean_pTsorted;
            vector<float> candidate_supercluster_emaxe_pTsorted;
            vector<float> candidate_supercluster_hoe_pTsorted;
            vector<float> candidate_supercluster_meanz_pTsorted;
            vector<float> candidate_supercluster_layer10_pTsorted;
            vector<float> candidate_supercluster_layer50_pTsorted;
            vector<float> candidate_supercluster_layer90_pTsorted;
            vector<float> candidate_supercluster_ntc67_pTsorted;
            vector<float> candidate_supercluster_ntc90_pTsorted;
            vector<float> candidate_supercluster_bdteg_pTsorted;
            vector<int>   candidate_supercluster_quality_pTsorted;
            vector<float> candidate_supercluster_puBDT_pTsorted;

            candidate_supercluster_pTsorted.clear();
            candidate_supercluster_pt_pTsorted.clear();
            candidate_supercluster_eta_pTsorted.clear();
            candidate_supercluster_phi_pTsorted.clear();
            candidate_supercluster_E_pTsorted.clear();
            candidate_supercluster_showerlength_pTsorted.clear();
            candidate_supercluster_coreshowerlength_pTsorted.clear();
            candidate_supercluster_firstlayer_pTsorted.clear();
            candidate_supercluster_maxlayer_pTsorted.clear();
            candidate_supercluster_seetot_pTsorted.clear();
            candidate_supercluster_seemax_pTsorted.clear();
            candidate_supercluster_spptot_pTsorted.clear();
            candidate_supercluster_sppmax_pTsorted.clear();
            candidate_supercluster_szz_pTsorted.clear();
            candidate_supercluster_srrtot_pTsorted.clear();
            candidate_supercluster_srrmax_pTsorted.clear();
            candidate_supercluster_srrmean_pTsorted.clear();
            candidate_supercluster_emaxe_pTsorted.clear();
            candidate_supercluster_hoe_pTsorted.clear();
            candidate_supercluster_meanz_pTsorted.clear();
            candidate_supercluster_layer10_pTsorted.clear();
            candidate_supercluster_layer50_pTsorted.clear();
            candidate_supercluster_layer90_pTsorted.clear();
            candidate_supercluster_ntc67_pTsorted.clear();
            candidate_supercluster_ntc90_pTsorted.clear();
            candidate_supercluster_bdteg_pTsorted.clear();
            candidate_supercluster_quality_pTsorted.clear();
            candidate_supercluster_puBDT_pTsorted.clear();

            vector< pair<int,TLorentzVector> > isupercluster_supercluster_pairs;

            for (unsigned int i_cluster = 0; i_cluster<candidate_supercluster.size(); i_cluster++){

                pair<int,TLorentzVector> supercluster_pair = make_pair(i_cluster,candidate_supercluster.at(i_cluster));
                isupercluster_supercluster_pairs.push_back(supercluster_pair);

            }

            sort(isupercluster_supercluster_pairs.begin(), isupercluster_supercluster_pairs.end(), pT_comparison_pairs);

            for (unsigned int i_cluster = 0; i_cluster<candidate_supercluster.size(); i_cluster++){

                int index = isupercluster_supercluster_pairs[i_cluster].first;
                TLorentzVector superclus = isupercluster_supercluster_pairs[i_cluster].second;

                candidate_supercluster_pTsorted.push_back(superclus);
                candidate_supercluster_pt_pTsorted.push_back(superclus.Pt());
                candidate_supercluster_E_pTsorted.push_back(superclus.E());
                candidate_supercluster_eta_pTsorted.push_back(superclus.Eta());
                candidate_supercluster_phi_pTsorted.push_back(superclus.Phi());
                candidate_supercluster_showerlength_pTsorted.push_back(candidate_supercluster_showerlength.at(index));
                candidate_supercluster_coreshowerlength_pTsorted.push_back(candidate_supercluster_coreshowerlength.at(index));
                candidate_supercluster_firstlayer_pTsorted.push_back(candidate_supercluster_firstlayer.at(index));
                candidate_supercluster_maxlayer_pTsorted.push_back(candidate_supercluster_maxlayer.at(index));
                candidate_supercluster_seetot_pTsorted.push_back(candidate_supercluster_seetot.at(index));
                candidate_supercluster_seemax_pTsorted.push_back(candidate_supercluster_seemax.at(index));
                candidate_supercluster_spptot_pTsorted.push_back(candidate_supercluster_spptot.at(index));
                candidate_supercluster_sppmax_pTsorted.push_back(candidate_supercluster_sppmax.at(index));
                candidate_supercluster_szz_pTsorted.push_back(candidate_supercluster_szz.at(index));
                candidate_supercluster_srrtot_pTsorted.push_back(candidate_supercluster_srrtot.at(index));
                candidate_supercluster_srrmax_pTsorted.push_back(candidate_supercluster_srrmax.at(index));
                candidate_supercluster_srrmean_pTsorted.push_back(candidate_supercluster_srrmean.at(index));
                candidate_supercluster_emaxe_pTsorted.push_back(candidate_supercluster_emaxe.at(index));
                candidate_supercluster_hoe_pTsorted.push_back(candidate_supercluster_hoe.at(index));
                candidate_supercluster_meanz_pTsorted.push_back(candidate_supercluster_meanz.at(index));
                candidate_supercluster_layer10_pTsorted.push_back(candidate_supercluster_layer10.at(index));
                candidate_supercluster_layer50_pTsorted.push_back(candidate_supercluster_layer50.at(index));
                candidate_supercluster_layer90_pTsorted.push_back(candidate_supercluster_layer90.at(index));
                candidate_supercluster_ntc67_pTsorted.push_back(candidate_supercluster_ntc67.at(index));
                candidate_supercluster_ntc90_pTsorted.push_back(candidate_supercluster_ntc90.at(index));
                candidate_supercluster_bdteg_pTsorted.push_back(candidate_supercluster_bdteg.at(index));
                candidate_supercluster_quality_pTsorted.push_back(candidate_supercluster_quality.at(index));
                candidate_supercluster_puBDT_pTsorted.push_back(candidate_supercluster_puBDT.at(index));

            }

            // CLEANING DUPLICATED SUPERCLUSTERS

            float new_pt_seeding_clus = candidate_supercluster_pTsorted.at(0).Pt();

            if(new_pt_seeding_clus == previous_pt_seeding_supercl) 
               continue;

            previous_pt_seeding_supercl = new_pt_seeding_clus;

            // SETTING GLOBAL VARIABLES

            int n_cl3d = candidate_supercluster_pt_pTsorted.size();
            _supercluster_n_cl3d.push_back(n_cl3d);

            float pt_tot = 0;
            float pt_seed = candidate_supercluster_pt_pTsorted.at(0);
            float eta_pt_tot = 0;
            float eta_seed = candidate_supercluster_eta_pTsorted.at(0);
            float phi_pt_tot = 0;
            float phi_seed = candidate_supercluster_phi_pTsorted.at(0);

            float max_eta = -999.;
            int   index_max_eta = -1;
            float min_eta = 999.;
            int   index_min_eta = -1;
            float max_phi = -999.;
            int   index_max_phi = -1;
            float min_phi = 999.;
            int   index_min_phi = -1;

            for(unsigned int subcl3d=0;subcl3d<candidate_supercluster_pt_pTsorted.size();subcl3d++) {

                pt_tot += candidate_supercluster_pt_pTsorted.at(subcl3d);
                eta_pt_tot += ((candidate_supercluster_pt_pTsorted.at(subcl3d))*(candidate_supercluster_eta_pTsorted.at(subcl3d)));
                phi_pt_tot += ((candidate_supercluster_pt_pTsorted.at(subcl3d))*(candidate_supercluster_phi_pTsorted.at(subcl3d)));

                if( candidate_supercluster_eta_pTsorted.at(subcl3d) > max_eta ) {
                    max_eta = candidate_supercluster_eta_pTsorted.at(subcl3d);
                    index_max_eta = subcl3d;
                }

                if( candidate_supercluster_eta_pTsorted.at(subcl3d) < min_eta ) {
                    min_eta = candidate_supercluster_eta_pTsorted.at(subcl3d);
                    index_min_eta = subcl3d;
                }

                if( candidate_supercluster_phi_pTsorted.at(subcl3d) > max_phi ) {
                    max_phi = candidate_supercluster_phi_pTsorted.at(subcl3d);
                    index_max_phi = subcl3d;
                }

                if( candidate_supercluster_phi_pTsorted.at(subcl3d) < min_phi ) {
                    min_phi = candidate_supercluster_phi_pTsorted.at(subcl3d);
                    index_min_phi = subcl3d;
                }

            }

            float eta_Eweighted = eta_pt_tot/pt_tot;
            float phi_Eweighted = phi_pt_tot/pt_tot;

            _supercluster_pt_tot.push_back(pt_tot);
            _supercluster_pt_seed.push_back(pt_seed);
            _supercluster_eta_seed.push_back(eta_seed);
            _supercluster_eta_Eweighted.push_back(eta_Eweighted);
            _supercluster_phi_seed.push_back(phi_seed);
            _supercluster_phi_Eweighted.push_back(phi_Eweighted);

            _supercluster_max_eta.push_back(max_eta);
            _supercluster_max_eta_index.push_back(index_max_eta);
            _supercluster_min_eta.push_back(min_eta);
            _supercluster_min_eta_index.push_back(index_min_eta);

            _supercluster_max_phi.push_back(max_phi);
            _supercluster_max_phi_index.push_back(index_max_phi);
            _supercluster_min_phi.push_back(min_phi);
            _supercluster_min_phi_index.push_back(index_min_phi);


            // SETTING cl3d VARIABLES

            _supercluster_cl3d_pt.push_back(candidate_supercluster_pt_pTsorted);
            _supercluster_cl3d_energy.push_back(candidate_supercluster_E_pTsorted);
            _supercluster_cl3d_eta.push_back(candidate_supercluster_eta_pTsorted);
            _supercluster_cl3d_phi.push_back(candidate_supercluster_phi_pTsorted);
            _supercluster_cl3d_showerlength.push_back(candidate_supercluster_showerlength_pTsorted);
            _supercluster_cl3d_coreshowerlength.push_back(candidate_supercluster_coreshowerlength_pTsorted);
            _supercluster_cl3d_firstlayer.push_back(candidate_supercluster_firstlayer_pTsorted);
            _supercluster_cl3d_maxlayer.push_back(candidate_supercluster_maxlayer_pTsorted);
            _supercluster_cl3d_seetot.push_back(candidate_supercluster_seetot_pTsorted);
            _supercluster_cl3d_seemax.push_back(candidate_supercluster_seemax_pTsorted);
            _supercluster_cl3d_spptot.push_back(candidate_supercluster_spptot_pTsorted);
            _supercluster_cl3d_sppmax.push_back(candidate_supercluster_sppmax_pTsorted);
            _supercluster_cl3d_szz.push_back(candidate_supercluster_szz_pTsorted);
            _supercluster_cl3d_srrtot.push_back(candidate_supercluster_srrtot_pTsorted);
            _supercluster_cl3d_srrmax.push_back(candidate_supercluster_srrmax_pTsorted);
            _supercluster_cl3d_srrmean.push_back(candidate_supercluster_srrmean_pTsorted);
            _supercluster_cl3d_emaxe.push_back(candidate_supercluster_emaxe_pTsorted);
            _supercluster_cl3d_hoe.push_back(candidate_supercluster_hoe_pTsorted);
            _supercluster_cl3d_meanz.push_back(candidate_supercluster_meanz_pTsorted);
            _supercluster_cl3d_layer10.push_back(candidate_supercluster_layer10_pTsorted);
            _supercluster_cl3d_layer50.push_back(candidate_supercluster_layer50_pTsorted);
            _supercluster_cl3d_layer90.push_back(candidate_supercluster_layer90_pTsorted);
            _supercluster_cl3d_ntc67.push_back(candidate_supercluster_ntc67_pTsorted);
            _supercluster_cl3d_ntc90.push_back(candidate_supercluster_ntc90_pTsorted);
            _supercluster_cl3d_bdteg.push_back(candidate_supercluster_bdteg_pTsorted);
            _supercluster_cl3d_quality.push_back(candidate_supercluster_quality_pTsorted);
            _supercluster_cl3d_puBDT.push_back(candidate_supercluster_puBDT_pTsorted);

        }

        _n_supercluster = _supercluster_cl3d_pt.size();

        // PRINT OUT

        if(debug){
            cout<<" .................................................. "<<endl;
            cout<<" Number of reconstructed superclusters: "<<_n_supercluster<<endl;
        }

        for(int i_supercl=0; i_supercl<_n_supercluster; i_supercl++){

            TLorentzVector supercluster;
            supercluster.SetPtEtaPhiM(_supercluster_pt_tot.at(i_supercl),_supercluster_eta_Eweighted.at(i_supercl),_supercluster_phi_Eweighted.at(i_supercl),0);

            vector<float> subclusters_pt;
            subclusters_pt = _supercluster_cl3d_pt.at(i_supercl);

            vector<float> subclusters_eta;
            subclusters_eta = _supercluster_cl3d_eta.at(i_supercl);

            vector<float> subclusters_phi;
            subclusters_phi = _supercluster_cl3d_phi.at(i_supercl);

            if(debug){
                cout<<"   Supercluster #"<<i_supercl<<", with "<<subclusters_pt.size()<<" subcluster(s) cl3d:"<<endl;
                for(unsigned int i_subcl=0; i_subcl<subclusters_pt.size();i_subcl++) 
                    cout<<"     Subcluster cl3d #"<<i_subcl<<": pt "<<subclusters_pt.at(i_subcl)<<", eta "<<subclusters_eta.at(i_subcl)<<", phi "<<subclusters_phi.at(i_subcl)<<endl;
            }
        }

        if(debug) {
            cout<<" .................................................. "<<endl;
            cout<<" Number of generated taus: "<<_gentau_n<<endl;
        }

        for(int i_gentau=0; i_gentau<_gentau_n;i_gentau++){

            TLorentzVector gentau;
            gentau.SetPtEtaPhiM( (*_gentau_pt)[i_gentau], (*_gentau_eta)[i_gentau], (*_gentau_phi)[i_gentau], (*_gentau_mass)[i_gentau]);

            TLorentzVector gentauvis;
            gentauvis.SetPtEtaPhiM( (*_gentau_vis_pt)[i_gentau], (*_gentau_vis_eta)[i_gentau], (*_gentau_vis_phi)[i_gentau], (*_gentau_vis_mass)[i_gentau]);

            if(debug) cout<<"   Gentau #"<<i_gentau<<": pt "<<(*_gentau_vis_pt)[i_gentau]<<", eta "<<(*_gentau_vis_eta)[i_gentau]<<", phi "<<(*_gentau_vis_phi)[i_gentau]<<endl;

        }

        if(debug) cout<<" .................................................. "<<endl;


        // MATCHING

        for(int i_gentau=0; i_gentau<_gentau_n;i_gentau++){

            TLorentzVector gentauvis;
            gentauvis.SetPtEtaPhiM( (*_gentau_vis_pt)[i_gentau], (*_gentau_vis_eta)[i_gentau], (*_gentau_vis_phi)[i_gentau], (*_gentau_vis_mass)[i_gentau]);

            bool  tau_isMatched = false;

            int   tau_n_cl3d = 0;
            
            float tau_pt_tot = -999;
            float tau_pt_seed = -999;
            float tau_eta_Eweighted = -999;
            float tau_eta_seed = -999;
            float tau_phi_Eweighted = -999;
            float tau_phi_seed = -999;

            int   tau_showerlength_seed = -999;
            int   tau_coreshowerlength_seed = -999;
            int   tau_firstlayer_seed = -999;
            int   tau_maxlayer_seed = -999;
            float tau_seetot_seed = -999;
            float tau_seemax_seed = -999;
            float tau_spptot_seed = -999;
            float tau_sppmax_seed = -999;
            float tau_szz_seed = -999;
            float tau_srrtot_seed = -999;
            float tau_srrmax_seed = -999;
            float tau_srrmean_seed = -999;
            float tau_emaxe_seed = -999;
            float tau_hoe_seed = -999;
            float tau_meanz_seed = -999;
            float tau_layer10_seed = -999;
            float tau_layer50_seed = -999;
            float tau_layer90_seed = -999;
            float tau_ntc67_seed = -999;
            float tau_ntc90_seed = -999;
            float tau_bdteg_seed = -999;
            int   tau_quality_seed = -999;
            float tau_puBDT_seed = -999;

            for(int i_supercl=0; i_supercl<_n_supercluster; i_supercl++){

                TLorentzVector supercluster;
                supercluster.SetPtEtaPhiM(_supercluster_pt_tot.at(i_supercl),_supercluster_eta_Eweighted.at(i_supercl),_supercluster_phi_Eweighted.at(i_supercl),0);

                vector<int>   cl3ds_showerlength;
                vector<int>   cl3ds_coreshowerlength;
                vector<int>   cl3ds_firstlayer;
                vector<int>   cl3ds_maxlayer;
                vector<float> cl3ds_seetot;
                vector<float> cl3ds_seemax;
                vector<float> cl3ds_spptot;
                vector<float> cl3ds_sppmax;
                vector<float> cl3ds_szz;
                vector<float> cl3ds_srrtot;
                vector<float> cl3ds_srrmax;
                vector<float> cl3ds_srrmean;
                vector<float> cl3ds_emaxe;
                vector<float> cl3ds_hoe;
                vector<float> cl3ds_meanz;
                vector<float> cl3ds_layer10;
                vector<float> cl3ds_layer50;
                vector<float> cl3ds_layer90;
                vector<float> cl3ds_ntc67;
                vector<float> cl3ds_ntc90;
                vector<float> cl3ds_bdteg;
                vector<int>   cl3ds_quality;
                vector<float> cl3ds_puBDT;

                cl3ds_showerlength = _supercluster_cl3d_showerlength.at(i_supercl);
                cl3ds_coreshowerlength = _supercluster_cl3d_coreshowerlength.at(i_supercl);
                cl3ds_firstlayer = _supercluster_cl3d_firstlayer.at(i_supercl);
                cl3ds_maxlayer = _supercluster_cl3d_maxlayer.at(i_supercl);
                cl3ds_seetot = _supercluster_cl3d_seetot.at(i_supercl);
                cl3ds_seemax = _supercluster_cl3d_seemax.at(i_supercl);
                cl3ds_spptot = _supercluster_cl3d_spptot.at(i_supercl);
                cl3ds_sppmax = _supercluster_cl3d_sppmax.at(i_supercl);
                cl3ds_szz = _supercluster_cl3d_szz.at(i_supercl);
                cl3ds_srrtot = _supercluster_cl3d_srrtot.at(i_supercl);
                cl3ds_srrmax = _supercluster_cl3d_srrmax.at(i_supercl);
                cl3ds_srrmean = _supercluster_cl3d_srrmean.at(i_supercl);
                cl3ds_emaxe = _supercluster_cl3d_emaxe.at(i_supercl);
                cl3ds_hoe = _supercluster_cl3d_hoe.at(i_supercl);
                cl3ds_meanz = _supercluster_cl3d_meanz.at(i_supercl);
                cl3ds_layer10 = _supercluster_cl3d_layer10.at(i_supercl);
                cl3ds_layer50 = _supercluster_cl3d_layer50.at(i_supercl);
                cl3ds_layer90 = _supercluster_cl3d_layer90.at(i_supercl);
                cl3ds_ntc67 = _supercluster_cl3d_ntc67.at(i_supercl);
                cl3ds_ntc90 = _supercluster_cl3d_ntc90.at(i_supercl);
                cl3ds_bdteg = _supercluster_cl3d_bdteg.at(i_supercl);
                cl3ds_quality = _supercluster_cl3d_quality.at(i_supercl);
                cl3ds_puBDT = _supercluster_cl3d_puBDT.at(i_supercl);

                float dR = supercluster.DeltaR(gentauvis);

                tau_isMatched = (dR <= 0.3);

                if( tau_isMatched ) {                 

                    tau_n_cl3d = _supercluster_n_cl3d.at(i_supercl);

                    tau_pt_tot = _supercluster_pt_tot.at(i_supercl);
                    tau_pt_seed = _supercluster_pt_seed.at(i_supercl);
                    tau_eta_Eweighted = _supercluster_eta_Eweighted.at(i_supercl);
                    tau_eta_seed = _supercluster_eta_seed.at(i_supercl);
                    tau_phi_Eweighted = _supercluster_phi_Eweighted.at(i_supercl);
                    tau_phi_seed = _supercluster_phi_seed.at(i_supercl);

                    tau_showerlength_seed = cl3ds_showerlength.at(0);
                    tau_coreshowerlength_seed = cl3ds_coreshowerlength.at(0);
                    tau_firstlayer_seed = cl3ds_firstlayer.at(0);
                    tau_maxlayer_seed = cl3ds_maxlayer.at(0);
                    tau_seetot_seed = cl3ds_seetot.at(0);
                    tau_seemax_seed = cl3ds_seemax.at(0);
                    tau_spptot_seed = cl3ds_spptot.at(0);
                    tau_sppmax_seed = cl3ds_sppmax.at(0);
                    tau_szz_seed = cl3ds_szz.at(0);
                    tau_srrtot_seed = cl3ds_srrtot.at(0);
                    tau_srrmax_seed = cl3ds_srrmax.at(0);
                    tau_srrmean_seed = cl3ds_srrmean.at(0);
                    tau_emaxe_seed = cl3ds_emaxe.at(0);
                    tau_hoe_seed = cl3ds_hoe.at(0);
                    tau_meanz_seed = cl3ds_meanz.at(0);
                    tau_layer10_seed = cl3ds_layer10.at(0);
                    tau_layer50_seed = cl3ds_layer50.at(0);
                    tau_layer90_seed = cl3ds_layer90.at(0);
                    tau_ntc67_seed = cl3ds_ntc67.at(0);
                    tau_ntc90_seed = cl3ds_ntc90.at(0);
                    tau_bdteg_seed = cl3ds_bdteg.at(0);
                    tau_quality_seed = cl3ds_quality.at(0);
                    tau_puBDT_seed = cl3ds_puBDT.at(0);

                    if(debug){
                        cout<<" --> Matched: gen tau eta "<<gentauvis.Eta()<<" to supercluster eta "<<supercluster.Eta()<<endl;
                        cout<<"       SC n cl3d: "<<tau_n_cl3d<<endl;
                        cout<<"       SC pt tot: "<<tau_pt_tot<<endl;
                        cout<<"       SC pt seed: "<<tau_pt_seed<<endl;
                        cout<<"       SC eta Eweight: "<<tau_eta_Eweighted<<endl;
                        cout<<"       SC eta seed: "<<tau_eta_seed<<endl;
                        cout<<"       SC phi Eweight: "<<tau_phi_Eweighted<<endl;
                        cout<<"       SC phi seed: "<<tau_phi_seed<<endl;   
                    }

                    break;
                }

            }

            _gentau_isMatched.push_back(tau_isMatched);

            _gentau_matchedSC_n_cl3d.push_back(tau_n_cl3d);
            _gentau_matchedSC_pt_tot.push_back(tau_pt_tot);
            _gentau_matchedSC_pt_seed.push_back(tau_pt_seed);
            _gentau_matchedSC_eta_Eweighted.push_back(tau_eta_Eweighted);
            _gentau_matchedSC_eta_seed.push_back(tau_eta_seed);
            _gentau_matchedSC_phi_Eweighted.push_back(tau_phi_Eweighted);
            _gentau_matchedSC_phi_seed.push_back(tau_phi_seed);

            _gentau_matchedSC_showerlength_seed.push_back(tau_showerlength_seed);
            _gentau_matchedSC_coreshowerlength_seed.push_back(tau_coreshowerlength_seed);
            _gentau_matchedSC_firstlayer_seed.push_back(tau_firstlayer_seed);
            _gentau_matchedSC_maxlayer_seed.push_back(tau_maxlayer_seed);
            _gentau_matchedSC_seetot_seed.push_back(tau_seetot_seed);
            _gentau_matchedSC_seemax_seed.push_back(tau_seemax_seed);
            _gentau_matchedSC_spptot_seed.push_back(tau_spptot_seed);
            _gentau_matchedSC_sppmax_seed.push_back(tau_sppmax_seed);
            _gentau_matchedSC_szz_seed.push_back(tau_szz_seed);
            _gentau_matchedSC_srrtot_seed.push_back(tau_srrtot_seed);
            _gentau_matchedSC_srrmax_seed.push_back(tau_srrmax_seed);
            _gentau_matchedSC_srrmean_seed.push_back(tau_srrmean_seed);
            _gentau_matchedSC_emaxe_seed.push_back(tau_emaxe_seed);
            _gentau_matchedSC_hoe_seed.push_back(tau_hoe_seed);
            _gentau_matchedSC_meanz_seed.push_back(tau_meanz_seed);
            _gentau_matchedSC_layer10_seed.push_back(tau_layer10_seed);
            _gentau_matchedSC_layer50_seed.push_back(tau_layer50_seed);
            _gentau_matchedSC_layer90_seed.push_back(tau_layer90_seed);
            _gentau_matchedSC_ntc67_seed.push_back(tau_ntc67_seed);
            _gentau_matchedSC_ntc90_seed.push_back(tau_ntc90_seed);
            _gentau_matchedSC_bdteg_seed.push_back(tau_bdteg_seed);
            _gentau_matchedSC_quality_seed.push_back(tau_quality_seed);
            _gentau_matchedSC_puBDT_seed.push_back(tau_puBDT_seed);

        }

        for(int i_gentau=0; i_gentau<_gentau_n;i_gentau++){

            bool matched = _gentau_isMatched[i_gentau];
            if(matched) matched_taus += 1;
            else if(!matched) unmatched_taus += 1;

        }

        out_tree->Fill();

    }

    out_file->cd();

    out_tree->Write();
    out_file->Close();

    float efficiency = matched_taus/(matched_taus+unmatched_taus);

    cout<<" "<<endl;
    cout<<" =========================================== "<<endl;
    cout<<" Matched taus: "<<matched_taus<<"; unmatched taus: "<<unmatched_taus<<" --> matching efficiency "<<efficiency*100<<"%"<<endl;
    cout<<" =========================================== "<<endl;


    return;
}


void test(int n_events = -1){

  //TString infile = "/data_CMS/cms/mperez/HGCal_data/Aug19/skimmed/ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_skimmed.root";
  //TString outfile = "/data_CMS/cms/mperez/HGCal_data/Aug19/clustered/ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_clustered.root";

  //TString infile = "/data_CMS/cms/mperez/HGCal_data/Aug19/skimmed/ntuple_RelValDiTau_Pt20To100_Eta1p6To2p9_skimmed.root";
  //TString outfile = "/data_CMS/cms/mperez/HGCal_data/Aug19/clustered/ntuple_RelValDiTau_Pt20To100_Eta1p6To2p9_clustered.root";

  TString infile = "/data_CMS/cms/mperez/HGCal_data/Aug19/skimmed/ntuple_Nu_E10_v10_PU200_files100to150_skimmed.root";
  TString outfile = "/data_CMS/cms/mperez/HGCal_data/Aug19/clustered/ntuple_Nu_E10_v10_PU200_files100to150_clustered.root";

  cluster_tree(infile, outfile, n_events);

}
