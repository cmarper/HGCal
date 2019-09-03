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

void cluster_tree( TString filein, TString fileout, int nevents = -1, bool debug = false,
                    float thr_cl3d_seed = 4, float thr_cl3d_sec = 2, 
                    float thr_tc_seed = 0.4, float thr_tc_sec = 0.2, 
                    float dEta = 0.1, float dPhi = 0.2, 
                    float dRmax_match = 0.3){

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

    ///////////////////

    int _event;

    ///////////////////

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

    ///////////////////

    int _tc_n;

    vector<unsigned int>    *_tc_id;
    vector<int>     *_tc_subdet;
    vector<int>     *_tc_zside;
    vector<int>     *_tc_layer;
    vector<int>     *_tc_waferu;
    vector<int>     *_tc_waferv;
    vector<int>     *_tc_wafertype;
    vector<int>     *_tc_panel_number;
    vector<int>     *_tc_panel_sector;
    vector<int>     *_tc_cell;
    vector<int>     *_tc_cellu;
    vector<int>     *_tc_cellv;  
    vector<unsigned int>    *_tc_data;
    vector<unsigned int>    *_tc_uncompressedCharge;
    vector<unsigned int>    *_tc_compressedCharge;

    vector<float>   *_tc_pt;
    vector<float>   *_tc_mipPt;
    vector<float>   *_tc_energy;
    vector<float>   *_tc_eta;
    vector<float>   *_tc_phi;
    vector<float>   *_tc_x;
    vector<float>   *_tc_y;
    vector<float>   *_tc_z;

    vector<unsigned int> *_tc_cluster_id;
    vector<unsigned int> *_tc_multicluster_id;
    vector<unsigned int> *_tc_multicluster_pt;

    ///////////////////

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

    ///////////////////

    in_tree->SetBranchAddress("event",&_event);

    ///////////////////

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

    ///////////////////

    in_tree->SetBranchAddress("tc_n",   &_tc_n);

    in_tree->SetBranchAddress("tc_id",      &_tc_id);
    in_tree->SetBranchAddress("tc_subdet",  &_tc_subdet);
    in_tree->SetBranchAddress("tc_zside",   &_tc_zside);
    in_tree->SetBranchAddress("tc_layer",   &_tc_layer);
    in_tree->SetBranchAddress("tc_waferu",  &_tc_waferu);
    in_tree->SetBranchAddress("tc_waferv",  &_tc_waferv);
    in_tree->SetBranchAddress("tc_wafertype",   &_tc_wafertype);
    in_tree->SetBranchAddress("tc_panel_number",    &_tc_panel_number);
    in_tree->SetBranchAddress("tc_panel_sector",    &_tc_panel_sector);
    in_tree->SetBranchAddress("tc_cellu",           &_tc_cellu);
    in_tree->SetBranchAddress("tc_cellv",           &_tc_cellv); 
    in_tree->SetBranchAddress("tc_data",            &_tc_data);
    in_tree->SetBranchAddress("tc_uncompressedCharge",  &_tc_uncompressedCharge);
    in_tree->SetBranchAddress("tc_compressedCharge",    &_tc_compressedCharge);

    in_tree->SetBranchAddress("tc_pt",              &_tc_pt);
    in_tree->SetBranchAddress("tc_mipPt",           &_tc_mipPt);
    in_tree->SetBranchAddress("tc_energy",          &_tc_energy);
    in_tree->SetBranchAddress("tc_eta",             &_tc_eta);
    in_tree->SetBranchAddress("tc_phi",             &_tc_phi);
    in_tree->SetBranchAddress("tc_x",               &_tc_x);
    in_tree->SetBranchAddress("tc_y",               &_tc_y);
    in_tree->SetBranchAddress("tc_z",               &_tc_z);

    in_tree->SetBranchAddress("tc_cluster_id",      &_tc_cluster_id);
    in_tree->SetBranchAddress("tc_multicluster_id", &_tc_multicluster_id);
    in_tree->SetBranchAddress("tc_multicluster_pt", &_tc_multicluster_pt);

    ///////////////////

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


    ///////////////////

    TTree* out_tree=in_tree->GetTree()->CloneTree(0);
    out_tree->SetNameTitle("ClusteredTree","ClusteredTree");

    ///////////////////

    vector<bool>   _gentau_isMatchedtoSTC;

    vector<int>    _gentau_matchedSTC_n_tc;
    vector<float>  _gentau_matchedSTC_pt_tot;
    vector<float>  _gentau_matchedSTC_pt_seed;
    vector<float>  _gentau_matchedSTC_eta_Eweighted;
    vector<float>  _gentau_matchedSTC_eta_seed;
    vector<float>  _gentau_matchedSTC_phi_Eweighted;
    vector<float>  _gentau_matchedSTC_phi_seed;
    vector<float>  _gentau_matchedSTC_x_Eweighted;
    vector<float>  _gentau_matchedSTC_x_seed;
    vector<float>  _gentau_matchedSTC_y_Eweighted;
    vector<float>  _gentau_matchedSTC_y_seed;
    vector<float>  _gentau_matchedSTC_z_Eweighted;
    vector<float>  _gentau_matchedSTC_z_seed;

    ///////////////////

    vector<bool>   _gentau_isMatchedtoSCL3D;

    vector<int>    _gentau_matchedSCL3D_n_cl3d;
    vector<float>  _gentau_matchedSCL3D_pt_tot;
    vector<float>  _gentau_matchedSCL3D_pt_seed;
    vector<float>  _gentau_matchedSCL3D_eta_Eweighted;
    vector<float>  _gentau_matchedSCL3D_eta_seed;
    vector<float>  _gentau_matchedSCL3D_phi_Eweighted;
    vector<float>  _gentau_matchedSCL3D_phi_seed;

    vector<int>    _gentau_matchedSCL3D_showerlength_seed;
    vector<int>    _gentau_matchedSCL3D_coreshowerlength_seed;
    vector<int>    _gentau_matchedSCL3D_firstlayer_seed;
    vector<int>    _gentau_matchedSCL3D_maxlayer_seed;
    vector<float>  _gentau_matchedSCL3D_seetot_seed;
    vector<float>  _gentau_matchedSCL3D_seemax_seed;
    vector<float>  _gentau_matchedSCL3D_spptot_seed;
    vector<float>  _gentau_matchedSCL3D_sppmax_seed;
    vector<float>  _gentau_matchedSCL3D_szz_seed;
    vector<float>  _gentau_matchedSCL3D_srrtot_seed;
    vector<float>  _gentau_matchedSCL3D_srrmax_seed;
    vector<float>  _gentau_matchedSCL3D_srrmean_seed;
    vector<float>  _gentau_matchedSCL3D_emaxe_seed;
    vector<float>  _gentau_matchedSCL3D_hoe_seed;
    vector<float>  _gentau_matchedSCL3D_meanz_seed;
    vector<float>  _gentau_matchedSCL3D_layer10_seed;
    vector<float>  _gentau_matchedSCL3D_layer50_seed;
    vector<float>  _gentau_matchedSCL3D_layer90_seed;
    vector<float>  _gentau_matchedSCL3D_ntc67_seed;
    vector<float>  _gentau_matchedSCL3D_ntc90_seed;
    vector<float>  _gentau_matchedSCL3D_bdteg_seed;
    vector<int>    _gentau_matchedSCL3D_quality_seed;
    vector<float>  _gentau_matchedSCL3D_puBDT_seed;

    ///////////////////

    int _n_STC;

    vector<int>    _STC_n_tc;
    vector<float>  _STC_pt_tot;
    vector<float>  _STC_pt_seed;
    vector<float>  _STC_eta_seed;
    vector<float>  _STC_eta_Eweighted;
    vector<float>  _STC_phi_seed;
    vector<float>  _STC_phi_Eweighted;  
    vector<float>  _STC_x_seed;
    vector<float>  _STC_x_Eweighted;
    vector<float>  _STC_y_seed;
    vector<float>  _STC_y_Eweighted;    
    vector<float>  _STC_z_seed;
    vector<float>  _STC_z_Eweighted;  

    vector<float> _STC_max_eta;
    vector<float> _STC_min_eta;
    vector<float> _STC_max_phi;
    vector<float> _STC_min_phi; 

    vector<bool>  _STC_isMatched;
    vector<int>   _STC_matchedGenTauIdx;

    vector<vector<float>>  _STC_tc_pt;
    vector<vector<float>>  _STC_tc_energy;
    vector<vector<float>>  _STC_tc_eta;
    vector<vector<float>>  _STC_tc_phi;
    vector<vector<float>>  _STC_tc_x;
    vector<vector<float>>  _STC_tc_y;
    vector<vector<float>>  _STC_tc_z;

    ///////////////////

    int _n_SCL3D;

    vector<int>    _SCL3D_n_cl3d;
    vector<float>  _SCL3D_pt_tot;
    vector<float>  _SCL3D_pt_seed;
    vector<float>  _SCL3D_eta_seed;
    vector<float>  _SCL3D_eta_Eweighted;
    vector<float>  _SCL3D_phi_seed;
    vector<float>  _SCL3D_phi_Eweighted;

    vector<float> _SCL3D_max_eta;
    vector<float> _SCL3D_min_eta;
    vector<float> _SCL3D_max_phi;
    vector<float> _SCL3D_min_phi;

    vector<bool>  _SCL3D_isMatched;
    vector<int>   _SCL3D_matchedGenTauIdx;

    vector<vector<float>>  _SCL3D_cl3d_pt;
    vector<vector<float>>  _SCL3D_cl3d_energy;
    vector<vector<float>>  _SCL3D_cl3d_eta;
    vector<vector<float>>  _SCL3D_cl3d_phi;

    vector<vector<int>>    _SCL3D_cl3d_showerlength;
    vector<vector<int>>    _SCL3D_cl3d_coreshowerlength;
    vector<vector<int>>    _SCL3D_cl3d_firstlayer;
    vector<vector<int>>    _SCL3D_cl3d_maxlayer;
    vector<vector<float>>  _SCL3D_cl3d_seetot;
    vector<vector<float>>  _SCL3D_cl3d_seemax;
    vector<vector<float>>  _SCL3D_cl3d_spptot;
    vector<vector<float>>  _SCL3D_cl3d_sppmax;
    vector<vector<float>>  _SCL3D_cl3d_szz;
    vector<vector<float>>  _SCL3D_cl3d_srrtot;
    vector<vector<float>>  _SCL3D_cl3d_srrmax;
    vector<vector<float>>  _SCL3D_cl3d_srrmean;
    vector<vector<float>>  _SCL3D_cl3d_emaxe;
    vector<vector<float>>  _SCL3D_cl3d_hoe;
    vector<vector<float>>  _SCL3D_cl3d_meanz;
    vector<vector<float>>  _SCL3D_cl3d_layer10;
    vector<vector<float>>  _SCL3D_cl3d_layer50;
    vector<vector<float>>  _SCL3D_cl3d_layer90;
    vector<vector<float>>  _SCL3D_cl3d_ntc67;
    vector<vector<float>>  _SCL3D_cl3d_ntc90;
    vector<vector<float>>  _SCL3D_cl3d_bdteg;
    vector<vector<int>>    _SCL3D_cl3d_quality;
    vector<vector<float>>  _SCL3D_cl3d_puBDT;

    ///////////////////

    out_tree->Branch("gentau_isMatchedtoSTC", &_gentau_isMatchedtoSTC);

    out_tree->Branch("gentau_matchedSTC_n_tc", &_gentau_matchedSTC_n_tc);
    out_tree->Branch("gentau_matchedSTC_pt_tot", &_gentau_matchedSTC_pt_tot);
    out_tree->Branch("gentau_matchedSTC_pt_seed", &_gentau_matchedSTC_pt_seed);
    out_tree->Branch("gentau_matchedSTC_eta_Eweighted", &_gentau_matchedSTC_eta_Eweighted);
    out_tree->Branch("gentau_matchedSTC_eta_seed", &_gentau_matchedSTC_eta_seed);
    out_tree->Branch("gentau_matchedSTC_phi_Eweighted", &_gentau_matchedSTC_phi_Eweighted);
    out_tree->Branch("gentau_matchedSTC_phi_seed", &_gentau_matchedSTC_phi_seed);
    out_tree->Branch("gentau_matchedSTC_x_Eweighted", &_gentau_matchedSTC_x_Eweighted);
    out_tree->Branch("gentau_matchedSTC_x_seed", &_gentau_matchedSTC_x_seed);
    out_tree->Branch("gentau_matchedSTC_y_Eweighted", &_gentau_matchedSTC_y_Eweighted);
    out_tree->Branch("gentau_matchedSTC_y_seed", &_gentau_matchedSTC_y_seed);
    out_tree->Branch("gentau_matchedSTC_z_Eweighted", &_gentau_matchedSTC_z_Eweighted);
    out_tree->Branch("gentau_matchedSTC_z_seed", &_gentau_matchedSTC_z_seed);

    ///////////////////

    out_tree->Branch("gentau_isMatchedtoSCL3D",  &_gentau_isMatchedtoSCL3D);

    out_tree->Branch("gentau_matchedSCL3D_n_cl3d", &_gentau_matchedSCL3D_n_cl3d);
    out_tree->Branch("gentau_matchedSCL3D_pt_tot", &_gentau_matchedSCL3D_pt_tot);
    out_tree->Branch("gentau_matchedSCL3D_pt_seed", &_gentau_matchedSCL3D_pt_seed);
    out_tree->Branch("gentau_matchedSCL3D_eta_Eweighted", &_gentau_matchedSCL3D_eta_Eweighted);
    out_tree->Branch("gentau_matchedSCL3D_eta_seed", &_gentau_matchedSCL3D_eta_seed);
    out_tree->Branch("gentau_matchedSCL3D_phi_Eweighted", &_gentau_matchedSCL3D_phi_Eweighted);
    out_tree->Branch("gentau_matchedSCL3D_phi_seed", &_gentau_matchedSCL3D_phi_seed);

    out_tree->Branch("gentau_matchedSCL3D_showerlength_seed", &_gentau_matchedSCL3D_showerlength_seed);
    out_tree->Branch("gentau_matchedSCL3D_coreshowerlength_seed", &_gentau_matchedSCL3D_coreshowerlength_seed);
    out_tree->Branch("gentau_matchedSCL3D_firstlayer_seed", &_gentau_matchedSCL3D_firstlayer_seed);
    out_tree->Branch("gentau_matchedSCL3D_maxlayer_seed", &_gentau_matchedSCL3D_maxlayer_seed);
    out_tree->Branch("gentau_matchedSCL3D_seetot_seed", &_gentau_matchedSCL3D_seetot_seed);
    out_tree->Branch("gentau_matchedSCL3D_seemax_seed", &_gentau_matchedSCL3D_seemax_seed);
    out_tree->Branch("gentau_matchedSCL3D_spptot_seed", &_gentau_matchedSCL3D_spptot_seed);
    out_tree->Branch("gentau_matchedSCL3D_sppmax_seed", &_gentau_matchedSCL3D_sppmax_seed);
    out_tree->Branch("gentau_matchedSCL3D_szz_seed", &_gentau_matchedSCL3D_szz_seed);
    out_tree->Branch("gentau_matchedSCL3D_srrtot_seed", &_gentau_matchedSCL3D_srrtot_seed);
    out_tree->Branch("gentau_matchedSCL3D_srrmax_seed", &_gentau_matchedSCL3D_srrmax_seed);
    out_tree->Branch("gentau_matchedSCL3D_srrmean_seed", &_gentau_matchedSCL3D_srrmean_seed);
    out_tree->Branch("gentau_matchedSCL3D_emaxe_seed", &_gentau_matchedSCL3D_emaxe_seed);
    out_tree->Branch("gentau_matchedSCL3D_hoe_seed", &_gentau_matchedSCL3D_hoe_seed);
    out_tree->Branch("gentau_matchedSCL3D_meanz_seed", &_gentau_matchedSCL3D_meanz_seed);
    out_tree->Branch("gentau_matchedSCL3D_layer10_seed", &_gentau_matchedSCL3D_layer10_seed);
    out_tree->Branch("gentau_matchedSCL3D_layer50_seed", &_gentau_matchedSCL3D_layer50_seed);
    out_tree->Branch("gentau_matchedSCL3D_layer90_seed", &_gentau_matchedSCL3D_layer90_seed);
    out_tree->Branch("gentau_matchedSCL3D_ntc67_seed", &_gentau_matchedSCL3D_ntc67_seed);
    out_tree->Branch("gentau_matchedSCL3D_ntc90_seed", &_gentau_matchedSCL3D_ntc90_seed);
    out_tree->Branch("gentau_matchedSCL3D_bdteg_seed", &_gentau_matchedSCL3D_bdteg_seed);
    out_tree->Branch("gentau_matchedSCL3D_quality_seed", &_gentau_matchedSCL3D_quality_seed);
    out_tree->Branch("gentau_matchedSCL3D_puBDT_seed", &_gentau_matchedSCL3D_puBDT_seed);

    ///////////////////

    out_tree->Branch("n_STC", &_n_STC);

    out_tree->Branch("STC_n_tc",            &_STC_n_tc);
    out_tree->Branch("STC_pt_tot",          &_STC_pt_tot);
    out_tree->Branch("STC_pt_seed",         &_STC_pt_seed);
    out_tree->Branch("STC_eta_seed",        &_STC_eta_seed);
    out_tree->Branch("STC_eta_Eweighted",   &_STC_eta_Eweighted);
    out_tree->Branch("STC_phi_seed",        &_STC_phi_seed);
    out_tree->Branch("STC_phi_Eweighted",   &_STC_phi_Eweighted); 
    out_tree->Branch("STC_x_seed",          &_STC_x_seed);
    out_tree->Branch("STC_x_Eweighted",     &_STC_x_Eweighted);
    out_tree->Branch("STC_y_seed",          &_STC_y_seed);
    out_tree->Branch("STC_y_Eweighted",     &_STC_y_Eweighted);
    out_tree->Branch("STC_z_seed",          &_STC_z_seed);
    out_tree->Branch("STC_z_Eweighted",     &_STC_z_Eweighted);

    out_tree->Branch("STC_max_eta", &_STC_max_eta);
    out_tree->Branch("STC_min_eta", &_STC_min_eta);
    out_tree->Branch("STC_max_phi", &_STC_max_phi);
    out_tree->Branch("STC_min_phi", &_STC_min_phi);

    out_tree->Branch("STC_isMatched",           &_STC_isMatched);
    out_tree->Branch("STC_matchedGenTauIdx",    &_STC_matchedGenTauIdx);

    out_tree->Branch("STC_tc_pt",     &_STC_tc_pt);
    out_tree->Branch("STC_tc_energy", &_STC_tc_energy);
    out_tree->Branch("STC_tc_eta",    &_STC_tc_eta);
    out_tree->Branch("STC_tc_phi",    &_STC_tc_phi);
    out_tree->Branch("STC_tc_x",      &_STC_tc_x);
    out_tree->Branch("STC_tc_y",      &_STC_tc_y);
    out_tree->Branch("STC_tc_z",      &_STC_tc_z);

    ///////////////////

    out_tree->Branch("n_SCL3D",     &_n_SCL3D);

    out_tree->Branch("SCL3D_n_cl3d",         &_SCL3D_n_cl3d);
    out_tree->Branch("SCL3D_pt_tot",         &_SCL3D_pt_tot);
    out_tree->Branch("SCL3D_pt_seed",        &_SCL3D_pt_seed);
    out_tree->Branch("SCL3D_eta_seed",       &_SCL3D_eta_seed);
    out_tree->Branch("SCL3D_eta_Eweighted",  &_SCL3D_eta_Eweighted);
    out_tree->Branch("SCL3D_phi_seedcl3d",   &_SCL3D_phi_seed);
    out_tree->Branch("SCL3D_phi_Eweighted",  &_SCL3D_phi_Eweighted);

    out_tree->Branch("SCL3D_max_eta",        &_SCL3D_max_eta);
    out_tree->Branch("SCL3D_min_eta",        &_SCL3D_min_eta);
    out_tree->Branch("SCL3D_max_phi",        &_SCL3D_max_phi);
    out_tree->Branch("SCL3D_min_phi",        &_SCL3D_min_phi);

    out_tree->Branch("SCL3D_isMatched",        &_SCL3D_isMatched);
    out_tree->Branch("SCL3D_matchedGenTauIdx", &_SCL3D_matchedGenTauIdx);

    out_tree->Branch("SCL3D_cl3d_pt",       &_SCL3D_cl3d_pt);
    out_tree->Branch("SCL3D_cl3d_energy",   &_SCL3D_cl3d_energy);
    out_tree->Branch("SCL3D_cl3d_eta",      &_SCL3D_cl3d_eta);
    out_tree->Branch("SCL3D_cl3d_phi",      &_SCL3D_cl3d_phi);

    out_tree->Branch("SCL3D_cl3d_showerlength",      &_SCL3D_cl3d_showerlength);
    out_tree->Branch("SCL3D_cl3d_coreshowerlength",  &_SCL3D_cl3d_coreshowerlength);
    out_tree->Branch("SCL3D_cl3d_firstlayer", &_SCL3D_cl3d_firstlayer);
    out_tree->Branch("SCL3D_cl3d_maxlayer",   &_SCL3D_cl3d_maxlayer);
    out_tree->Branch("SCL3D_cl3d_seetot",     &_SCL3D_cl3d_seetot);
    out_tree->Branch("SCL3D_cl3d_seemax",     &_SCL3D_cl3d_seemax);
    out_tree->Branch("SCL3D_cl3d_spptot",     &_SCL3D_cl3d_spptot);
    out_tree->Branch("SCL3D_cl3d_sppmax",     &_SCL3D_cl3d_sppmax);
    out_tree->Branch("SCL3D_cl3d_szz",        &_SCL3D_cl3d_szz);
    out_tree->Branch("SCL3D_cl3d_srrtot",     &_SCL3D_cl3d_srrtot);
    out_tree->Branch("SCL3D_cl3d_srrmax",     &_SCL3D_cl3d_srrmax);
    out_tree->Branch("SCL3D_cl3d_srrmean",    &_SCL3D_cl3d_srrmean);
    out_tree->Branch("SCL3D_cl3d_emaxe",      &_SCL3D_cl3d_emaxe);
    out_tree->Branch("SCL3D_cl3d_hoe",        &_SCL3D_cl3d_hoe);
    out_tree->Branch("SCL3D_cl3d_meanz",      &_SCL3D_cl3d_meanz);
    out_tree->Branch("SCL3D_cl3d_layer10",    &_SCL3D_cl3d_layer10);
    out_tree->Branch("SCL3D_cl3d_layer50",    &_SCL3D_cl3d_layer50);
    out_tree->Branch("SCL3D_cl3d_layer90",    &_SCL3D_cl3d_layer90);
    out_tree->Branch("SCL3D_cl3d_ntc67",      &_SCL3D_cl3d_ntc67);
    out_tree->Branch("SCL3D_cl3d_ntc90",      &_SCL3D_cl3d_ntc90);
    out_tree->Branch("SCL3D_cl3d_bdteg",      &_SCL3D_cl3d_bdteg);
    out_tree->Branch("SCL3D_cl3d_quality",    &_SCL3D_cl3d_quality);
    out_tree->Branch("SCL3D_cl3d_puBDT",      &_SCL3D_cl3d_puBDT);

    ///////////////////

    float matchedtoSTC_taus = 0;
    float matchedtoSCL3D_taus = 0;

    float unmatchedtoSTC_taus = 0;
    float unmatchedtoSCL3D_taus = 0;

    ///////////////////

    for (int i=0;i<nentries;i++) {

        if(i%1000==0) cout<<"i="<<i<<endl;

        // old branches

        ///////////////////

        _event = 0;

        ///////////////////

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

        ///////////////////

        _tc_n = 0; 

        _tc_id = 0; 
        _tc_subdet = 0; 
        _tc_zside = 0; 
        _tc_layer = 0; 
        _tc_waferu = 0; 
        _tc_waferv = 0; 
        _tc_wafertype = 0; 
        _tc_panel_number = 0; 
        _tc_panel_sector = 0; 
        _tc_cell = 0; 
        _tc_cellu = 0; 
        _tc_cellv = 0; 
        _tc_data = 0; 
        _tc_uncompressedCharge = 0; 
        _tc_compressedCharge = 0; 

        _tc_pt = 0; 
        _tc_mipPt = 0; 
        _tc_energy = 0; 
        _tc_eta = 0; 
        _tc_phi = 0; 
        _tc_x = 0; 
        _tc_y = 0; 
        _tc_z = 0; 

        _tc_cluster_id = 0;
        _tc_multicluster_id = 0;
        _tc_multicluster_pt = 0;

        ///////////////////

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

        ///////////////////

        _gentau_isMatchedtoSTC.clear();

        _gentau_matchedSTC_n_tc.clear();
        _gentau_matchedSTC_pt_tot.clear();
        _gentau_matchedSTC_pt_seed.clear();
        _gentau_matchedSTC_eta_Eweighted.clear();
        _gentau_matchedSTC_eta_seed.clear();
        _gentau_matchedSTC_phi_Eweighted.clear();
        _gentau_matchedSTC_phi_seed.clear();
        _gentau_matchedSTC_x_Eweighted.clear();
        _gentau_matchedSTC_x_seed.clear();
        _gentau_matchedSTC_y_Eweighted.clear();
        _gentau_matchedSTC_y_seed.clear();
        _gentau_matchedSTC_z_Eweighted.clear();
        _gentau_matchedSTC_z_seed.clear();

        ///////////////////

        _gentau_isMatchedtoSCL3D.clear();

        _gentau_matchedSCL3D_n_cl3d.clear();
        _gentau_matchedSCL3D_pt_tot.clear();
        _gentau_matchedSCL3D_pt_seed.clear();
        _gentau_matchedSCL3D_eta_Eweighted.clear();
        _gentau_matchedSCL3D_eta_seed.clear();
        _gentau_matchedSCL3D_phi_Eweighted.clear();
        _gentau_matchedSCL3D_phi_seed.clear();

        _gentau_matchedSCL3D_showerlength_seed.clear();
        _gentau_matchedSCL3D_coreshowerlength_seed.clear();
        _gentau_matchedSCL3D_firstlayer_seed.clear();
        _gentau_matchedSCL3D_maxlayer_seed.clear();
        _gentau_matchedSCL3D_seetot_seed.clear();
        _gentau_matchedSCL3D_seemax_seed.clear();
        _gentau_matchedSCL3D_spptot_seed.clear();
        _gentau_matchedSCL3D_sppmax_seed.clear();
        _gentau_matchedSCL3D_szz_seed.clear();
        _gentau_matchedSCL3D_srrtot_seed.clear();
        _gentau_matchedSCL3D_srrmax_seed.clear();
        _gentau_matchedSCL3D_srrmean_seed.clear();
        _gentau_matchedSCL3D_emaxe_seed.clear();
        _gentau_matchedSCL3D_hoe_seed.clear();
        _gentau_matchedSCL3D_meanz_seed.clear();
        _gentau_matchedSCL3D_layer10_seed.clear();
        _gentau_matchedSCL3D_layer50_seed.clear();
        _gentau_matchedSCL3D_layer90_seed.clear();
        _gentau_matchedSCL3D_ntc67_seed.clear();
        _gentau_matchedSCL3D_ntc90_seed.clear();
        _gentau_matchedSCL3D_bdteg_seed.clear();
        _gentau_matchedSCL3D_quality_seed.clear();
        _gentau_matchedSCL3D_puBDT_seed.clear();

        ///////////////////

        _n_STC = 0;

        _STC_n_tc.clear(); 
        _STC_pt_tot.clear(); 
        _STC_pt_seed.clear(); 
        _STC_eta_seed.clear(); 
        _STC_eta_Eweighted.clear(); 
        _STC_phi_seed.clear(); 
        _STC_phi_Eweighted.clear(); 
        _STC_x_seed.clear(); 
        _STC_x_Eweighted.clear();
        _STC_y_seed.clear(); 
        _STC_y_Eweighted.clear();
        _STC_z_seed.clear(); 
        _STC_z_Eweighted.clear();

        _STC_max_eta.clear(); 
        _STC_min_eta.clear(); 
        _STC_max_phi.clear(); 
        _STC_min_phi.clear(); 

        _STC_isMatched.clear();
        _STC_matchedGenTauIdx.clear();

        _STC_tc_pt.clear();
        _STC_tc_energy.clear();
        _STC_tc_eta.clear();
        _STC_tc_phi.clear();

        ///////////////////

        _n_SCL3D = 0;

        _SCL3D_n_cl3d.clear(); 
        _SCL3D_pt_tot.clear(); 
        _SCL3D_pt_seed.clear(); 
        _SCL3D_eta_seed.clear(); 
        _SCL3D_eta_Eweighted.clear(); 
        _SCL3D_phi_seed.clear(); 
        _SCL3D_phi_Eweighted.clear(); 

        _SCL3D_max_eta.clear(); 
        _SCL3D_min_eta.clear(); 
        _SCL3D_max_phi.clear(); 
        _SCL3D_min_phi.clear(); 

        _SCL3D_isMatched.clear();
        _SCL3D_matchedGenTauIdx.clear();

        _SCL3D_cl3d_pt.clear();
        _SCL3D_cl3d_energy.clear();
        _SCL3D_cl3d_eta.clear();
        _SCL3D_cl3d_phi.clear();

        _SCL3D_cl3d_showerlength.clear();
        _SCL3D_cl3d_coreshowerlength.clear();
        _SCL3D_cl3d_firstlayer.clear();
        _SCL3D_cl3d_maxlayer.clear();
        _SCL3D_cl3d_seetot.clear();
        _SCL3D_cl3d_seemax.clear();
        _SCL3D_cl3d_spptot.clear();
        _SCL3D_cl3d_sppmax.clear();
        _SCL3D_cl3d_szz.clear();
        _SCL3D_cl3d_srrtot.clear();
        _SCL3D_cl3d_srrmax.clear();
        _SCL3D_cl3d_srrmean.clear();
        _SCL3D_cl3d_emaxe.clear();
        _SCL3D_cl3d_hoe.clear();
        _SCL3D_cl3d_meanz.clear();
        _SCL3D_cl3d_layer10.clear();
        _SCL3D_cl3d_layer50.clear();
        _SCL3D_cl3d_layer90.clear();
        _SCL3D_cl3d_ntc67.clear();
        _SCL3D_cl3d_ntc90.clear();
        _SCL3D_cl3d_bdteg.clear();
        _SCL3D_cl3d_quality.clear();
        _SCL3D_cl3d_puBDT.clear();

        ///////////////////

        int entry_ok = in_tree->GetEntry(i);    
        if(entry_ok<0) 
            continue;

        ///////////////////

        vector<float> previous_pt_seeding_STC;
        vector<float> previous_eta_seeding_STC;
        vector<float> previous_phi_seeding_STC;

        previous_pt_seeding_STC.clear();
        previous_eta_seeding_STC.clear();
        previous_phi_seeding_STC.clear();

        previous_pt_seeding_STC.push_back(-999);
        previous_eta_seeding_STC.push_back(-999);
        previous_phi_seeding_STC.push_back(-999);

        vector<float> previous_pt_seeding_SCL3D;
        vector<float> previous_eta_seeding_SCL3D;
        vector<float> previous_phi_seeding_SCL3D;

        previous_pt_seeding_SCL3D.clear();
        previous_eta_seeding_SCL3D.clear();
        previous_phi_seeding_SCL3D.clear();

        previous_pt_seeding_SCL3D.push_back(-999);
        previous_eta_seeding_SCL3D.push_back(-999);
        previous_phi_seeding_SCL3D.push_back(-999);


        //////////////////////////
        //// CLUSTERING OF TC ////
        //////////////////////////

        /*if(debug){
            cout<<" "<<endl;
            cout<<" ==================> Event "<<i<< "  <================== "<<endl;
            cout<<" "<<endl;
            cout<<" Number of reconstructed tc: "<<_tc_n<<endl;
        }

        for (int i_main=0; i_main<_tc_n; i_main++){

            // RAW CLUSTERING 

            vector<TLorentzVector>  candidate_STC;
            vector<int>   candidate_STC_n_tc;
            vector<float> candidate_STC_x;
            vector<float> candidate_STC_y;
            vector<float> candidate_STC_z;

            candidate_STC.clear();
            candidate_STC_n_tc.clear();
            candidate_STC_x.clear();
            candidate_STC_y.clear();
            candidate_STC_z.clear();

            // MAIN CLUSTER CANDIDATE (SUPERCLUSTER SEED)

            TLorentzVector seed_candidate_STC;
            float seed_candidate_STC_x;
            float seed_candidate_STC_y;
            float seed_candidate_STC_z;

            seed_candidate_STC.SetPtEtaPhiM( (*_tc_pt)[i_main], (*_tc_eta)[i_main], (*_tc_phi)[i_main], 0);
            seed_candidate_STC_x = (*_tc_x)[i_main];
            seed_candidate_STC_y = (*_tc_y)[i_main];
            seed_candidate_STC_z = (*_tc_z)[i_main];

            if ( seed_candidate_STC.Pt() < thr_tc_seed ) continue;
        
            //if(debug) cout<<"   Seed tc candidate #"<<i_main<<": pT "<<seed_candidate_STC.Pt()<<", eta "<<seed_candidate_STC.Eta()<<", phi "<<seed_candidate_STC.Phi()<<endl;

            candidate_STC.push_back(seed_candidate_STC);
            candidate_STC_x.push_back(seed_candidate_STC_x);
            candidate_STC_y.push_back(seed_candidate_STC_y);
            candidate_STC_z.push_back(seed_candidate_STC_z);

            // SECONDARY CLUSTERS CANDIDATES

            for (int i_sec=0; i_sec<_tc_n; i_sec++){

                if( i_sec==i_main ) continue;

                TLorentzVector secondary_candidate_STC;
                float secondary_candidate_STC_x;
                float secondary_candidate_STC_y;
                float secondary_candidate_STC_z;

                secondary_candidate_STC.SetPtEtaPhiM( (*_tc_pt)[i_sec], (*_tc_eta)[i_sec], (*_tc_phi)[i_sec], 0);
                secondary_candidate_STC_x = (*_tc_x)[i_sec];
                secondary_candidate_STC_y = (*_tc_y)[i_sec];
                secondary_candidate_STC_z = (*_tc_z)[i_sec];

                if( secondary_candidate_STC.Pt() < thr_tc_sec ) continue;
 
                if( fabs(secondary_candidate_STC.Eta() - seed_candidate_STC.Eta()) > dEta ) continue; // (TT/2)*(eta/TT) -> 3
                if( fabs(seed_candidate_STC.DeltaPhi(secondary_candidate_STC)) > dPhi ) continue; // => 9
                //if( fabs(seed_candidate_STC.DeltaR(secondary_candidate_STC)) > 0.5) continue;

                //if(debug) cout<<"      Secondary tc candidate #"<<i_sec<<": pT "<<secondary_candidate_STC.Pt()<<", eta "<<secondary_candidate_STC.Eta()<<", phi "<<secondary_candidate_STC.Phi()<<endl;

                candidate_STC.push_back(secondary_candidate_STC);
                candidate_STC_x.push_back(secondary_candidate_STC_x);
                candidate_STC_y.push_back(secondary_candidate_STC_y);
                candidate_STC_z.push_back(secondary_candidate_STC_z);

            }

            // SORT CLUSTERS BY PT IN SUPERCLUSTER

            vector<TLorentzVector>  candidate_STC_pTsorted;
            vector<float> candidate_STC_pt_pTsorted;
            vector<float> candidate_STC_eta_pTsorted;
            vector<float> candidate_STC_phi_pTsorted;
            vector<float> candidate_STC_E_pTsorted;
            vector<float> candidate_STC_x_pTsorted;
            vector<float> candidate_STC_y_pTsorted;
            vector<float> candidate_STC_z_pTsorted;

            candidate_STC_pTsorted.clear();
            candidate_STC_pt_pTsorted.clear();
            candidate_STC_eta_pTsorted.clear();
            candidate_STC_phi_pTsorted.clear();
            candidate_STC_E_pTsorted.clear();
            candidate_STC_x_pTsorted.clear();
            candidate_STC_y_pTsorted.clear();
            candidate_STC_z_pTsorted.clear();

            vector< pair<int,TLorentzVector> > iSTC_STC_pairs;

            for (unsigned int i_cluster = 0; i_cluster<candidate_STC.size(); i_cluster++){

                pair<int,TLorentzVector> STC_pair = make_pair(i_cluster,candidate_STC.at(i_cluster));
                iSTC_STC_pairs.push_back(STC_pair);

            }

            sort(iSTC_STC_pairs.begin(), iSTC_STC_pairs.end(), pT_comparison_pairs);

            for (unsigned int i_cluster = 0; i_cluster<candidate_STC.size(); i_cluster++){

                int index = iSTC_STC_pairs[i_cluster].first;
                TLorentzVector superclus = iSTC_STC_pairs[i_cluster].second;

                candidate_STC_pTsorted.push_back(superclus);
                candidate_STC_pt_pTsorted.push_back(superclus.Pt());
                candidate_STC_E_pTsorted.push_back(superclus.E());
                candidate_STC_eta_pTsorted.push_back(superclus.Eta());
                candidate_STC_phi_pTsorted.push_back(superclus.Phi());
                candidate_STC_x_pTsorted.push_back(candidate_STC_x.at(index));
                candidate_STC_y_pTsorted.push_back(candidate_STC_y.at(index));
                candidate_STC_z_pTsorted.push_back(candidate_STC_z.at(index));

            }

            // CLEANING DUPLICATED SUPERCLUSTERS

            float new_pt_seeding_STC = candidate_STC_pTsorted.at(0).Pt();
            float new_eta_seeding_STC = candidate_STC_pTsorted.at(0).Eta();
            float new_phi_seeding_STC = candidate_STC_pTsorted.at(0).Phi();

            bool foundbefore = false;
            for(unsigned int i=0; i<previous_pt_seeding_STC.size(); i++){
                if(new_pt_seeding_STC == previous_pt_seeding_STC.at(i) && new_eta_seeding_STC == previous_eta_seeding_STC.at(i) && new_phi_seeding_STC == previous_phi_seeding_STC.at(i)){
                    foundbefore = true;
                    break;
                }
            }

            if(foundbefore) continue;

            previous_pt_seeding_STC.push_back(new_pt_seeding_STC);
            previous_eta_seeding_STC.push_back(new_eta_seeding_STC);
            previous_phi_seeding_STC.push_back(new_phi_seeding_STC);

            // SETTING GLOBAL VARIABLES

            int n_tc = candidate_STC_pt_pTsorted.size();
            _STC_n_tc.push_back(n_tc);

            float pt_tot = 0;
            float pt_seed = candidate_STC_pt_pTsorted.at(0);
            float eta_pt_tot = 0;
            float eta_seed = candidate_STC_eta_pTsorted.at(0);
            float phi_pt_tot = 0;
            float phi_seed = candidate_STC_phi_pTsorted.at(0);
            float x_pt_tot = 0;
            float x_seed = candidate_STC_x_pTsorted.at(0);
            float y_pt_tot = 0;
            float y_seed = candidate_STC_y_pTsorted.at(0);
            float z_pt_tot = 0;
            float z_seed = candidate_STC_z_pTsorted.at(0);                                   

            float max_eta = -999.;
            int   index_max_eta = -1;
            float min_eta = 999.;
            int   index_min_eta = -1;
            float max_phi = -999.;
            int   index_max_phi = -1;
            float min_phi = 999.;
            int   index_min_phi = -1;

            for(unsigned int subtc=0;subtc<candidate_STC_pt_pTsorted.size();subtc++) {

                pt_tot += candidate_STC_pt_pTsorted.at(subtc);
                eta_pt_tot += ((candidate_STC_pt_pTsorted.at(subtc))*(candidate_STC_eta_pTsorted.at(subtc)));
                phi_pt_tot += ((candidate_STC_pt_pTsorted.at(subtc))*(candidate_STC_phi_pTsorted.at(subtc)));
                x_pt_tot += ((candidate_STC_pt_pTsorted.at(subtc))*(candidate_STC_x_pTsorted.at(subtc)));
                y_pt_tot += ((candidate_STC_pt_pTsorted.at(subtc))*(candidate_STC_y_pTsorted.at(subtc)));
                z_pt_tot += ((candidate_STC_pt_pTsorted.at(subtc))*(candidate_STC_z_pTsorted.at(subtc)));

                if( candidate_STC_eta_pTsorted.at(subtc) > max_eta ) {
                    max_eta = candidate_STC_eta_pTsorted.at(subtc);
                    index_max_eta = subtc;
                }

                if( candidate_STC_eta_pTsorted.at(subtc) < min_eta ) {
                    min_eta = candidate_STC_eta_pTsorted.at(subtc);
                    index_min_eta = subtc;
                }

                if( candidate_STC_phi_pTsorted.at(subtc) > max_phi ) {
                    max_phi = candidate_STC_phi_pTsorted.at(subtc);
                    index_max_phi = subtc;
                }

                if( candidate_STC_phi_pTsorted.at(subtc) < min_phi ) {
                    min_phi = candidate_STC_phi_pTsorted.at(subtc);
                    index_min_phi = subtc;
                }

            }

            float eta_Eweighted = eta_pt_tot/pt_tot;
            float phi_Eweighted = phi_pt_tot/pt_tot;
            float x_Eweighted = x_pt_tot/pt_tot;
            float y_Eweighted = y_pt_tot/pt_tot;
            float z_Eweighted = z_pt_tot/pt_tot;

            _STC_pt_tot.push_back(pt_tot);
            _STC_pt_seed.push_back(pt_seed);
            _STC_eta_seed.push_back(eta_seed);
            _STC_eta_Eweighted.push_back(eta_Eweighted);
            _STC_phi_seed.push_back(phi_seed);
            _STC_phi_Eweighted.push_back(phi_Eweighted);
            _STC_x_seed.push_back(x_seed);
            _STC_x_Eweighted.push_back(x_Eweighted);
            _STC_y_seed.push_back(y_seed);
            _STC_y_Eweighted.push_back(y_Eweighted);
            _STC_z_seed.push_back(z_seed);
            _STC_z_Eweighted.push_back(z_Eweighted);

            _STC_max_eta.push_back(max_eta);
            _STC_min_eta.push_back(min_eta);
            _STC_max_phi.push_back(max_phi);
            _STC_min_phi.push_back(min_phi);


            // SETTING tc VARIABLES

            _STC_tc_pt.push_back(candidate_STC_pt_pTsorted);
            _STC_tc_energy.push_back(candidate_STC_E_pTsorted);
            _STC_tc_eta.push_back(candidate_STC_eta_pTsorted);
            _STC_tc_phi.push_back(candidate_STC_phi_pTsorted);
            _STC_tc_x.push_back(candidate_STC_x_pTsorted);
            _STC_tc_y.push_back(candidate_STC_y_pTsorted);
            _STC_tc_z.push_back(candidate_STC_z_pTsorted);

        }

        _n_STC = _STC_tc_pt.size();

        // PRINT OUT

        if(debug){
            cout<<" .................................................. "<<endl;
            cout<<" Number of reconstructed STCs: "<<_n_STC<<endl;
        }

        for(int i_supercl=0; i_supercl<_n_STC; i_supercl++){

            TLorentzVector STC;
            STC.SetPtEtaPhiM(_STC_pt_tot.at(i_supercl),_STC_eta_Eweighted.at(i_supercl),_STC_phi_Eweighted.at(i_supercl),0);

            vector<float> subclusters_pt;
            subclusters_pt = _STC_tc_pt.at(i_supercl);

            vector<float> subclusters_eta;
            subclusters_eta = _STC_tc_eta.at(i_supercl);

            vector<float> subclusters_phi;
            subclusters_phi = _STC_tc_phi.at(i_supercl);

            vector<float> subclusters_x;
            subclusters_x = _STC_tc_x.at(i_supercl);

            vector<float> subclusters_y;
            subclusters_x = _STC_tc_y.at(i_supercl);

            vector<float> subclusters_z;
            subclusters_x = _STC_tc_z.at(i_supercl);

            if(debug){
                cout<<"   STC #"<<i_supercl<<", with "<<subclusters_pt.size()<<" subcluster(s) tc:"<<endl;
                for(unsigned int i_subcl=0; i_subcl<subclusters_pt.size();i_subcl++) 
                    cout<<"     Subcluster tc #"<<i_subcl<<": pt "<<subclusters_pt.at(i_subcl)<<", eta "<<subclusters_eta.at(i_subcl)<<", phi "<<subclusters_phi.at(i_subcl)<<endl;
            }
        }
        */

        ////////////////////////////
        //// CLUSTERING OF CL3D ////
        ////////////////////////////

        if(debug){
            cout<<" .................................................. "<<endl;
            cout<<" Number of reconstructed cl3d: "<<_cl3d_n<<endl;
        }

        for (int i_main=0; i_main<_cl3d_n; i_main++){

            // RAW CLUSTERING 

            vector<TLorentzVector>  candidate_SCL3D;
            vector<int>   candidate_SCL3D_n_cl3d;
            vector<int>   candidate_SCL3D_showerlength;
            vector<int>   candidate_SCL3D_coreshowerlength;
            vector<int>   candidate_SCL3D_firstlayer;
            vector<int>   candidate_SCL3D_maxlayer;
            vector<float> candidate_SCL3D_seetot;
            vector<float> candidate_SCL3D_seemax;
            vector<float> candidate_SCL3D_spptot;
            vector<float> candidate_SCL3D_sppmax;
            vector<float> candidate_SCL3D_szz;
            vector<float> candidate_SCL3D_srrtot;
            vector<float> candidate_SCL3D_srrmax;
            vector<float> candidate_SCL3D_srrmean;
            vector<float> candidate_SCL3D_emaxe;
            vector<float> candidate_SCL3D_hoe;
            vector<float> candidate_SCL3D_meanz;
            vector<float> candidate_SCL3D_layer10;
            vector<float> candidate_SCL3D_layer50;
            vector<float> candidate_SCL3D_layer90;
            vector<float> candidate_SCL3D_ntc67;
            vector<float> candidate_SCL3D_ntc90;
            vector<float> candidate_SCL3D_bdteg;
            vector<int>   candidate_SCL3D_quality;
            vector<float> candidate_SCL3D_puBDT;

            candidate_SCL3D.clear();
            candidate_SCL3D_n_cl3d.clear();
            candidate_SCL3D_showerlength.clear();
            candidate_SCL3D_coreshowerlength.clear();
            candidate_SCL3D_firstlayer.clear();
            candidate_SCL3D_maxlayer.clear();
            candidate_SCL3D_seetot.clear();
            candidate_SCL3D_seemax.clear();
            candidate_SCL3D_spptot.clear();
            candidate_SCL3D_sppmax.clear();
            candidate_SCL3D_szz.clear();
            candidate_SCL3D_srrtot.clear();
            candidate_SCL3D_srrmax.clear();
            candidate_SCL3D_srrmean.clear();
            candidate_SCL3D_emaxe.clear();
            candidate_SCL3D_hoe.clear();
            candidate_SCL3D_meanz.clear();
            candidate_SCL3D_layer10.clear();
            candidate_SCL3D_layer50.clear();
            candidate_SCL3D_layer90.clear();
            candidate_SCL3D_ntc67.clear();
            candidate_SCL3D_ntc90.clear();
            candidate_SCL3D_bdteg.clear();
            candidate_SCL3D_quality.clear();
            candidate_SCL3D_puBDT.clear();

            // MAIN CLUSTER CANDIDATE (SUPERCLUSTER SEED)

            TLorentzVector seed_candidate_SCL3D;
            int   seed_candidate_SCL3D_showerlength;
            int   seed_candidate_SCL3D_coreshowerlength;
            int   seed_candidate_SCL3D_firstlayer;
            int   seed_candidate_SCL3D_maxlayer;
            float seed_candidate_SCL3D_seetot;
            float seed_candidate_SCL3D_seemax;
            float seed_candidate_SCL3D_spptot;
            float seed_candidate_SCL3D_sppmax;
            float seed_candidate_SCL3D_szz;
            float seed_candidate_SCL3D_srrtot;
            float seed_candidate_SCL3D_srrmax;
            float seed_candidate_SCL3D_srrmean;
            float seed_candidate_SCL3D_emaxe;
            float seed_candidate_SCL3D_hoe;
            float seed_candidate_SCL3D_meanz;
            float seed_candidate_SCL3D_layer10;
            float seed_candidate_SCL3D_layer50;
            float seed_candidate_SCL3D_layer90;
            float seed_candidate_SCL3D_ntc67;
            float seed_candidate_SCL3D_ntc90;
            float seed_candidate_SCL3D_bdteg;
            int   seed_candidate_SCL3D_quality;
            float seed_candidate_SCL3D_puBDT;

            seed_candidate_SCL3D.SetPtEtaPhiM( (*_cl3d_pt)[i_main], (*_cl3d_eta)[i_main], (*_cl3d_phi)[i_main], 0);
            seed_candidate_SCL3D_showerlength = (*_cl3d_showerlength)[i_main];
            seed_candidate_SCL3D_coreshowerlength = (*_cl3d_coreshowerlength)[i_main];
            seed_candidate_SCL3D_firstlayer = (*_cl3d_firstlayer)[i_main];
            seed_candidate_SCL3D_maxlayer = (*_cl3d_maxlayer)[i_main];
            seed_candidate_SCL3D_seetot = (*_cl3d_seetot)[i_main];
            seed_candidate_SCL3D_seemax = (*_cl3d_seemax)[i_main];
            seed_candidate_SCL3D_spptot = (*_cl3d_spptot)[i_main];
            seed_candidate_SCL3D_sppmax = (*_cl3d_sppmax)[i_main];
            seed_candidate_SCL3D_szz = (*_cl3d_szz)[i_main];
            seed_candidate_SCL3D_srrtot = (*_cl3d_srrtot)[i_main];
            seed_candidate_SCL3D_srrmax = (*_cl3d_srrmax)[i_main];
            seed_candidate_SCL3D_srrmean = (*_cl3d_srrmean)[i_main];
            seed_candidate_SCL3D_emaxe = (*_cl3d_emaxe)[i_main];
            seed_candidate_SCL3D_hoe = (*_cl3d_hoe)[i_main];
            seed_candidate_SCL3D_meanz = (*_cl3d_meanz)[i_main];
            seed_candidate_SCL3D_layer10 = (*_cl3d_layer10)[i_main];
            seed_candidate_SCL3D_layer50 = (*_cl3d_layer50)[i_main];
            seed_candidate_SCL3D_layer90 = (*_cl3d_layer90)[i_main];
            seed_candidate_SCL3D_ntc67 = (*_cl3d_ntc67)[i_main];
            seed_candidate_SCL3D_ntc90 = (*_cl3d_ntc90)[i_main];
            seed_candidate_SCL3D_bdteg = (*_cl3d_bdteg)[i_main];
            seed_candidate_SCL3D_quality = (*_cl3d_quality)[i_main];
            seed_candidate_SCL3D_puBDT = (*_cl3d_puBDT)[i_main];

            if ( seed_candidate_SCL3D_puBDT < 0 ) continue;

            if ( seed_candidate_SCL3D.Pt() < thr_cl3d_seed ) continue;
        
            //if(debug) cout<<"   Seed cl3d candidate #"<<i_main<<": pT "<<seed_candidate_SCL3D.Pt()<<", eta "<<seed_candidate_SCL3D.Eta()<<", phi "<<seed_candidate_SCL3D.Phi()<<endl;

            candidate_SCL3D.push_back(seed_candidate_SCL3D);
            candidate_SCL3D_showerlength.push_back(seed_candidate_SCL3D_showerlength);
            candidate_SCL3D_coreshowerlength.push_back(seed_candidate_SCL3D_coreshowerlength);
            candidate_SCL3D_firstlayer.push_back(seed_candidate_SCL3D_firstlayer);
            candidate_SCL3D_maxlayer.push_back(seed_candidate_SCL3D_maxlayer);
            candidate_SCL3D_seetot.push_back(seed_candidate_SCL3D_seetot);
            candidate_SCL3D_seemax.push_back(seed_candidate_SCL3D_seemax);
            candidate_SCL3D_spptot.push_back(seed_candidate_SCL3D_spptot);
            candidate_SCL3D_sppmax.push_back(seed_candidate_SCL3D_sppmax);
            candidate_SCL3D_szz.push_back(seed_candidate_SCL3D_szz);
            candidate_SCL3D_srrtot.push_back(seed_candidate_SCL3D_srrtot);
            candidate_SCL3D_srrmax.push_back(seed_candidate_SCL3D_srrmax);
            candidate_SCL3D_srrmean.push_back(seed_candidate_SCL3D_srrmean);
            candidate_SCL3D_emaxe.push_back(seed_candidate_SCL3D_emaxe);
            candidate_SCL3D_hoe.push_back(seed_candidate_SCL3D_hoe);
            candidate_SCL3D_meanz.push_back(seed_candidate_SCL3D_meanz);
            candidate_SCL3D_layer10.push_back(seed_candidate_SCL3D_layer10);
            candidate_SCL3D_layer50.push_back(seed_candidate_SCL3D_layer50);
            candidate_SCL3D_layer90.push_back(seed_candidate_SCL3D_layer90);
            candidate_SCL3D_ntc67.push_back(seed_candidate_SCL3D_ntc67);
            candidate_SCL3D_ntc90.push_back(seed_candidate_SCL3D_ntc90);
            candidate_SCL3D_bdteg.push_back(seed_candidate_SCL3D_bdteg);
            candidate_SCL3D_quality.push_back(seed_candidate_SCL3D_quality);
            candidate_SCL3D_puBDT.push_back(seed_candidate_SCL3D_puBDT);

            // SECONDARY CLUSTERS CANDIDATES

            for (int i_sec=0; i_sec<_cl3d_n; i_sec++){

                if( i_sec==i_main ) continue;

                TLorentzVector secondary_candidate_SCL3D;
                int   secondary_candidate_SCL3D_showerlength;
                int   secondary_candidate_SCL3D_coreshowerlength;
                int   secondary_candidate_SCL3D_firstlayer;
                int   secondary_candidate_SCL3D_maxlayer;
                float secondary_candidate_SCL3D_seetot;
                float secondary_candidate_SCL3D_seemax;
                float secondary_candidate_SCL3D_spptot;
                float secondary_candidate_SCL3D_sppmax;
                float secondary_candidate_SCL3D_szz;
                float secondary_candidate_SCL3D_srrtot;
                float secondary_candidate_SCL3D_srrmax;
                float secondary_candidate_SCL3D_srrmean;
                float secondary_candidate_SCL3D_emaxe;
                float secondary_candidate_SCL3D_hoe;
                float secondary_candidate_SCL3D_meanz;
                float secondary_candidate_SCL3D_layer10;
                float secondary_candidate_SCL3D_layer50;
                float secondary_candidate_SCL3D_layer90;
                float secondary_candidate_SCL3D_ntc67;
                float secondary_candidate_SCL3D_ntc90;
                float secondary_candidate_SCL3D_bdteg;
                int   secondary_candidate_SCL3D_quality;
                float secondary_candidate_SCL3D_puBDT;

                secondary_candidate_SCL3D.SetPtEtaPhiM( (*_cl3d_pt)[i_sec], (*_cl3d_eta)[i_sec], (*_cl3d_phi)[i_sec], 0);
                secondary_candidate_SCL3D_showerlength = (*_cl3d_showerlength)[i_sec];
                secondary_candidate_SCL3D_coreshowerlength = (*_cl3d_coreshowerlength)[i_sec];
                secondary_candidate_SCL3D_firstlayer = (*_cl3d_firstlayer)[i_sec];
                secondary_candidate_SCL3D_maxlayer = (*_cl3d_maxlayer)[i_sec];
                secondary_candidate_SCL3D_seetot = (*_cl3d_seetot)[i_sec];
                secondary_candidate_SCL3D_seemax = (*_cl3d_seemax)[i_sec];
                secondary_candidate_SCL3D_spptot = (*_cl3d_spptot)[i_sec];
                secondary_candidate_SCL3D_sppmax = (*_cl3d_sppmax)[i_sec];
                secondary_candidate_SCL3D_szz = (*_cl3d_szz)[i_sec];
                secondary_candidate_SCL3D_srrtot = (*_cl3d_srrtot)[i_sec];
                secondary_candidate_SCL3D_srrmax = (*_cl3d_srrmax)[i_sec];
                secondary_candidate_SCL3D_srrmean = (*_cl3d_srrmean)[i_sec];
                secondary_candidate_SCL3D_emaxe = (*_cl3d_emaxe)[i_sec];
                secondary_candidate_SCL3D_hoe = (*_cl3d_hoe)[i_sec];
                secondary_candidate_SCL3D_meanz = (*_cl3d_meanz)[i_sec];
                secondary_candidate_SCL3D_layer10 = (*_cl3d_layer10)[i_sec];
                secondary_candidate_SCL3D_layer50 = (*_cl3d_layer50)[i_sec];
                secondary_candidate_SCL3D_layer90 = (*_cl3d_layer90)[i_sec];
                secondary_candidate_SCL3D_ntc67 = (*_cl3d_ntc67)[i_sec];
                secondary_candidate_SCL3D_ntc90 = (*_cl3d_ntc90)[i_sec];
                secondary_candidate_SCL3D_bdteg = (*_cl3d_bdteg)[i_sec];
                secondary_candidate_SCL3D_quality = (*_cl3d_quality)[i_sec];
                secondary_candidate_SCL3D_puBDT = (*_cl3d_puBDT)[i_sec];

                if( secondary_candidate_SCL3D_puBDT < 0 ) continue;

                if( secondary_candidate_SCL3D.Pt() < thr_cl3d_sec ) continue;
 
                if( fabs(secondary_candidate_SCL3D.Eta() - seed_candidate_SCL3D.Eta()) > dEta ) continue; // (TT/2)*(eta/TT) -> 3
                if( fabs(seed_candidate_SCL3D.DeltaPhi(secondary_candidate_SCL3D)) > dPhi ) continue; // => 9
                //if( fabs(seed_candidate_SCL3D.DeltaR(secondary_candidate_SCL3D)) > 0.5) continue;

                //if(debug) cout<<"      Secondary cl3d candidate #"<<i_sec<<": pT "<<secondary_candidate_SCL3D.Pt()<<", eta "<<secondary_candidate_SCL3D.Eta()<<", phi "<<secondary_candidate_SCL3D.Phi()<<endl;

                candidate_SCL3D.push_back(secondary_candidate_SCL3D);
                candidate_SCL3D_showerlength.push_back(secondary_candidate_SCL3D_showerlength);
                candidate_SCL3D_coreshowerlength.push_back(secondary_candidate_SCL3D_coreshowerlength);
                candidate_SCL3D_firstlayer.push_back(secondary_candidate_SCL3D_firstlayer);
                candidate_SCL3D_maxlayer.push_back(secondary_candidate_SCL3D_maxlayer);
                candidate_SCL3D_seetot.push_back(secondary_candidate_SCL3D_seetot);
                candidate_SCL3D_seemax.push_back(secondary_candidate_SCL3D_seemax);
                candidate_SCL3D_spptot.push_back(secondary_candidate_SCL3D_spptot);
                candidate_SCL3D_sppmax.push_back(secondary_candidate_SCL3D_sppmax);
                candidate_SCL3D_szz.push_back(secondary_candidate_SCL3D_szz);
                candidate_SCL3D_srrtot.push_back(secondary_candidate_SCL3D_srrtot);
                candidate_SCL3D_srrmax.push_back(secondary_candidate_SCL3D_srrmax);
                candidate_SCL3D_srrmean.push_back(secondary_candidate_SCL3D_srrmean);
                candidate_SCL3D_emaxe.push_back(secondary_candidate_SCL3D_emaxe);
                candidate_SCL3D_hoe.push_back(secondary_candidate_SCL3D_hoe);
                candidate_SCL3D_meanz.push_back(secondary_candidate_SCL3D_meanz);
                candidate_SCL3D_layer10.push_back(secondary_candidate_SCL3D_layer10);
                candidate_SCL3D_layer50.push_back(secondary_candidate_SCL3D_layer50);
                candidate_SCL3D_layer90.push_back(secondary_candidate_SCL3D_layer90);
                candidate_SCL3D_ntc67.push_back(secondary_candidate_SCL3D_ntc67);
                candidate_SCL3D_ntc90.push_back(secondary_candidate_SCL3D_ntc90);
                candidate_SCL3D_bdteg.push_back(secondary_candidate_SCL3D_bdteg);
                candidate_SCL3D_quality.push_back(secondary_candidate_SCL3D_quality);
                candidate_SCL3D_puBDT.push_back(secondary_candidate_SCL3D_puBDT);

            }

            // SORT CLUSTERS BY PT IN SUPERCLUSTER

            vector<TLorentzVector>  candidate_SCL3D_pTsorted;
            vector<float> candidate_SCL3D_pt_pTsorted;
            vector<float> candidate_SCL3D_eta_pTsorted;
            vector<float> candidate_SCL3D_phi_pTsorted;
            vector<float> candidate_SCL3D_E_pTsorted;
            vector<int>   candidate_SCL3D_showerlength_pTsorted;
            vector<int>   candidate_SCL3D_coreshowerlength_pTsorted;
            vector<int>   candidate_SCL3D_firstlayer_pTsorted;
            vector<int>   candidate_SCL3D_maxlayer_pTsorted;
            vector<float> candidate_SCL3D_seetot_pTsorted;
            vector<float> candidate_SCL3D_seemax_pTsorted;
            vector<float> candidate_SCL3D_spptot_pTsorted;
            vector<float> candidate_SCL3D_sppmax_pTsorted;
            vector<float> candidate_SCL3D_szz_pTsorted;
            vector<float> candidate_SCL3D_srrtot_pTsorted;
            vector<float> candidate_SCL3D_srrmax_pTsorted;
            vector<float> candidate_SCL3D_srrmean_pTsorted;
            vector<float> candidate_SCL3D_emaxe_pTsorted;
            vector<float> candidate_SCL3D_hoe_pTsorted;
            vector<float> candidate_SCL3D_meanz_pTsorted;
            vector<float> candidate_SCL3D_layer10_pTsorted;
            vector<float> candidate_SCL3D_layer50_pTsorted;
            vector<float> candidate_SCL3D_layer90_pTsorted;
            vector<float> candidate_SCL3D_ntc67_pTsorted;
            vector<float> candidate_SCL3D_ntc90_pTsorted;
            vector<float> candidate_SCL3D_bdteg_pTsorted;
            vector<int>   candidate_SCL3D_quality_pTsorted;
            vector<float> candidate_SCL3D_puBDT_pTsorted;

            candidate_SCL3D_pTsorted.clear();
            candidate_SCL3D_pt_pTsorted.clear();
            candidate_SCL3D_eta_pTsorted.clear();
            candidate_SCL3D_phi_pTsorted.clear();
            candidate_SCL3D_E_pTsorted.clear();
            candidate_SCL3D_showerlength_pTsorted.clear();
            candidate_SCL3D_coreshowerlength_pTsorted.clear();
            candidate_SCL3D_firstlayer_pTsorted.clear();
            candidate_SCL3D_maxlayer_pTsorted.clear();
            candidate_SCL3D_seetot_pTsorted.clear();
            candidate_SCL3D_seemax_pTsorted.clear();
            candidate_SCL3D_spptot_pTsorted.clear();
            candidate_SCL3D_sppmax_pTsorted.clear();
            candidate_SCL3D_szz_pTsorted.clear();
            candidate_SCL3D_srrtot_pTsorted.clear();
            candidate_SCL3D_srrmax_pTsorted.clear();
            candidate_SCL3D_srrmean_pTsorted.clear();
            candidate_SCL3D_emaxe_pTsorted.clear();
            candidate_SCL3D_hoe_pTsorted.clear();
            candidate_SCL3D_meanz_pTsorted.clear();
            candidate_SCL3D_layer10_pTsorted.clear();
            candidate_SCL3D_layer50_pTsorted.clear();
            candidate_SCL3D_layer90_pTsorted.clear();
            candidate_SCL3D_ntc67_pTsorted.clear();
            candidate_SCL3D_ntc90_pTsorted.clear();
            candidate_SCL3D_bdteg_pTsorted.clear();
            candidate_SCL3D_quality_pTsorted.clear();
            candidate_SCL3D_puBDT_pTsorted.clear();

            vector< pair<int,TLorentzVector> > iSCL3D_SCL3D_pairs;

            for (unsigned int i_cluster = 0; i_cluster<candidate_SCL3D.size(); i_cluster++){

                pair<int,TLorentzVector> SCL3D_pair = make_pair(i_cluster,candidate_SCL3D.at(i_cluster));
                iSCL3D_SCL3D_pairs.push_back(SCL3D_pair);

            }

            sort(iSCL3D_SCL3D_pairs.begin(), iSCL3D_SCL3D_pairs.end(), pT_comparison_pairs);

            for (unsigned int i_cluster = 0; i_cluster<candidate_SCL3D.size(); i_cluster++){

                int index = iSCL3D_SCL3D_pairs[i_cluster].first;
                TLorentzVector superclus = iSCL3D_SCL3D_pairs[i_cluster].second;

                candidate_SCL3D_pTsorted.push_back(superclus);
                candidate_SCL3D_pt_pTsorted.push_back(superclus.Pt());
                candidate_SCL3D_E_pTsorted.push_back(superclus.E());
                candidate_SCL3D_eta_pTsorted.push_back(superclus.Eta());
                candidate_SCL3D_phi_pTsorted.push_back(superclus.Phi());
                candidate_SCL3D_showerlength_pTsorted.push_back(candidate_SCL3D_showerlength.at(index));
                candidate_SCL3D_coreshowerlength_pTsorted.push_back(candidate_SCL3D_coreshowerlength.at(index));
                candidate_SCL3D_firstlayer_pTsorted.push_back(candidate_SCL3D_firstlayer.at(index));
                candidate_SCL3D_maxlayer_pTsorted.push_back(candidate_SCL3D_maxlayer.at(index));
                candidate_SCL3D_seetot_pTsorted.push_back(candidate_SCL3D_seetot.at(index));
                candidate_SCL3D_seemax_pTsorted.push_back(candidate_SCL3D_seemax.at(index));
                candidate_SCL3D_spptot_pTsorted.push_back(candidate_SCL3D_spptot.at(index));
                candidate_SCL3D_sppmax_pTsorted.push_back(candidate_SCL3D_sppmax.at(index));
                candidate_SCL3D_szz_pTsorted.push_back(candidate_SCL3D_szz.at(index));
                candidate_SCL3D_srrtot_pTsorted.push_back(candidate_SCL3D_srrtot.at(index));
                candidate_SCL3D_srrmax_pTsorted.push_back(candidate_SCL3D_srrmax.at(index));
                candidate_SCL3D_srrmean_pTsorted.push_back(candidate_SCL3D_srrmean.at(index));
                candidate_SCL3D_emaxe_pTsorted.push_back(candidate_SCL3D_emaxe.at(index));
                candidate_SCL3D_hoe_pTsorted.push_back(candidate_SCL3D_hoe.at(index));
                candidate_SCL3D_meanz_pTsorted.push_back(candidate_SCL3D_meanz.at(index));
                candidate_SCL3D_layer10_pTsorted.push_back(candidate_SCL3D_layer10.at(index));
                candidate_SCL3D_layer50_pTsorted.push_back(candidate_SCL3D_layer50.at(index));
                candidate_SCL3D_layer90_pTsorted.push_back(candidate_SCL3D_layer90.at(index));
                candidate_SCL3D_ntc67_pTsorted.push_back(candidate_SCL3D_ntc67.at(index));
                candidate_SCL3D_ntc90_pTsorted.push_back(candidate_SCL3D_ntc90.at(index));
                candidate_SCL3D_bdteg_pTsorted.push_back(candidate_SCL3D_bdteg.at(index));
                candidate_SCL3D_quality_pTsorted.push_back(candidate_SCL3D_quality.at(index));
                candidate_SCL3D_puBDT_pTsorted.push_back(candidate_SCL3D_puBDT.at(index));

            }

            // CLEANING DUPLICATED SUPERCLUSTERS

            float new_pt_seeding_SCL3D = candidate_SCL3D_pTsorted.at(0).Pt();
            float new_eta_seeding_SCL3D = candidate_SCL3D_pTsorted.at(0).Eta();
            float new_phi_seeding_SCL3D = candidate_SCL3D_pTsorted.at(0).Phi();

            bool foundbefore = false;
            for(unsigned int i=0; i<previous_pt_seeding_SCL3D.size(); i++){
                if(new_pt_seeding_SCL3D == previous_pt_seeding_SCL3D.at(i) && new_eta_seeding_SCL3D == previous_eta_seeding_SCL3D.at(i) && new_phi_seeding_SCL3D == previous_phi_seeding_SCL3D.at(i)){
                    foundbefore = true;
                    break;
                }
            }

            if(foundbefore) continue;

            previous_pt_seeding_SCL3D.push_back(new_pt_seeding_SCL3D);
            previous_eta_seeding_SCL3D.push_back(new_eta_seeding_SCL3D);
            previous_phi_seeding_SCL3D.push_back(new_phi_seeding_SCL3D);

            // SETTING GLOBAL VARIABLES

            int n_cl3d = candidate_SCL3D_pt_pTsorted.size();
            _SCL3D_n_cl3d.push_back(n_cl3d);

            float pt_tot = 0;
            float pt_seed = candidate_SCL3D_pt_pTsorted.at(0);
            float eta_pt_tot = 0;
            float eta_seed = candidate_SCL3D_eta_pTsorted.at(0);
            float phi_pt_tot = 0;
            float phi_seed = candidate_SCL3D_phi_pTsorted.at(0);

            float max_eta = -999.;
            int   index_max_eta = -1;
            float min_eta = 999.;
            int   index_min_eta = -1;
            float max_phi = -999.;
            int   index_max_phi = -1;
            float min_phi = 999.;
            int   index_min_phi = -1;

            for(unsigned int subcl3d=0;subcl3d<candidate_SCL3D_pt_pTsorted.size();subcl3d++) {

                pt_tot += candidate_SCL3D_pt_pTsorted.at(subcl3d);
                eta_pt_tot += ((candidate_SCL3D_pt_pTsorted.at(subcl3d))*(candidate_SCL3D_eta_pTsorted.at(subcl3d)));
                phi_pt_tot += ((candidate_SCL3D_pt_pTsorted.at(subcl3d))*(candidate_SCL3D_phi_pTsorted.at(subcl3d)));

                if( candidate_SCL3D_eta_pTsorted.at(subcl3d) > max_eta ) {
                    max_eta = candidate_SCL3D_eta_pTsorted.at(subcl3d);
                    index_max_eta = subcl3d;
                }

                if( candidate_SCL3D_eta_pTsorted.at(subcl3d) < min_eta ) {
                    min_eta = candidate_SCL3D_eta_pTsorted.at(subcl3d);
                    index_min_eta = subcl3d;
                }

                if( candidate_SCL3D_phi_pTsorted.at(subcl3d) > max_phi ) {
                    max_phi = candidate_SCL3D_phi_pTsorted.at(subcl3d);
                    index_max_phi = subcl3d;
                }

                if( candidate_SCL3D_phi_pTsorted.at(subcl3d) < min_phi ) {
                    min_phi = candidate_SCL3D_phi_pTsorted.at(subcl3d);
                    index_min_phi = subcl3d;
                }

            }

            float eta_Eweighted = eta_pt_tot/pt_tot;
            float phi_Eweighted = phi_pt_tot/pt_tot;

            _SCL3D_pt_tot.push_back(pt_tot);
            _SCL3D_pt_seed.push_back(pt_seed);
            _SCL3D_eta_seed.push_back(eta_seed);
            _SCL3D_eta_Eweighted.push_back(eta_Eweighted);
            _SCL3D_phi_seed.push_back(phi_seed);
            _SCL3D_phi_Eweighted.push_back(phi_Eweighted);

            _SCL3D_max_eta.push_back(max_eta);
            _SCL3D_min_eta.push_back(min_eta);
            _SCL3D_max_phi.push_back(max_phi);
            _SCL3D_min_phi.push_back(min_phi);


            // SETTING cl3d VARIABLES

            _SCL3D_cl3d_pt.push_back(candidate_SCL3D_pt_pTsorted);
            _SCL3D_cl3d_energy.push_back(candidate_SCL3D_E_pTsorted);
            _SCL3D_cl3d_eta.push_back(candidate_SCL3D_eta_pTsorted);
            _SCL3D_cl3d_phi.push_back(candidate_SCL3D_phi_pTsorted);
            _SCL3D_cl3d_showerlength.push_back(candidate_SCL3D_showerlength_pTsorted);
            _SCL3D_cl3d_coreshowerlength.push_back(candidate_SCL3D_coreshowerlength_pTsorted);
            _SCL3D_cl3d_firstlayer.push_back(candidate_SCL3D_firstlayer_pTsorted);
            _SCL3D_cl3d_maxlayer.push_back(candidate_SCL3D_maxlayer_pTsorted);
            _SCL3D_cl3d_seetot.push_back(candidate_SCL3D_seetot_pTsorted);
            _SCL3D_cl3d_seemax.push_back(candidate_SCL3D_seemax_pTsorted);
            _SCL3D_cl3d_spptot.push_back(candidate_SCL3D_spptot_pTsorted);
            _SCL3D_cl3d_sppmax.push_back(candidate_SCL3D_sppmax_pTsorted);
            _SCL3D_cl3d_szz.push_back(candidate_SCL3D_szz_pTsorted);
            _SCL3D_cl3d_srrtot.push_back(candidate_SCL3D_srrtot_pTsorted);
            _SCL3D_cl3d_srrmax.push_back(candidate_SCL3D_srrmax_pTsorted);
            _SCL3D_cl3d_srrmean.push_back(candidate_SCL3D_srrmean_pTsorted);
            _SCL3D_cl3d_emaxe.push_back(candidate_SCL3D_emaxe_pTsorted);
            _SCL3D_cl3d_hoe.push_back(candidate_SCL3D_hoe_pTsorted);
            _SCL3D_cl3d_meanz.push_back(candidate_SCL3D_meanz_pTsorted);
            _SCL3D_cl3d_layer10.push_back(candidate_SCL3D_layer10_pTsorted);
            _SCL3D_cl3d_layer50.push_back(candidate_SCL3D_layer50_pTsorted);
            _SCL3D_cl3d_layer90.push_back(candidate_SCL3D_layer90_pTsorted);
            _SCL3D_cl3d_ntc67.push_back(candidate_SCL3D_ntc67_pTsorted);
            _SCL3D_cl3d_ntc90.push_back(candidate_SCL3D_ntc90_pTsorted);
            _SCL3D_cl3d_bdteg.push_back(candidate_SCL3D_bdteg_pTsorted);
            _SCL3D_cl3d_quality.push_back(candidate_SCL3D_quality_pTsorted);
            _SCL3D_cl3d_puBDT.push_back(candidate_SCL3D_puBDT_pTsorted);

        }

        _n_SCL3D = _SCL3D_cl3d_pt.size();

        // PRINT OUT

        if(debug){
            cout<<" .................................................. "<<endl;
            cout<<" Number of reconstructed SCL3Ds: "<<_n_SCL3D<<endl;
        }

        for(int i_supercl=0; i_supercl<_n_SCL3D; i_supercl++){

            TLorentzVector SCL3D;
            SCL3D.SetPtEtaPhiM(_SCL3D_pt_tot.at(i_supercl),_SCL3D_eta_Eweighted.at(i_supercl),_SCL3D_phi_Eweighted.at(i_supercl),0);

            vector<float> subclusters_pt;
            subclusters_pt = _SCL3D_cl3d_pt.at(i_supercl);

            vector<float> subclusters_eta;
            subclusters_eta = _SCL3D_cl3d_eta.at(i_supercl);

            vector<float> subclusters_phi;
            subclusters_phi = _SCL3D_cl3d_phi.at(i_supercl);

            if(debug){
                cout<<"   SCL3D #"<<i_supercl<<", with "<<subclusters_pt.size()<<" subcluster(s) cl3d:"<<endl;
                for(unsigned int i_subcl=0; i_subcl<subclusters_pt.size();i_subcl++) 
                    cout<<"     Subcluster cl3d #"<<i_subcl<<": pt "<<subclusters_pt.at(i_subcl)<<", eta "<<subclusters_eta.at(i_subcl)<<", phi "<<subclusters_phi.at(i_subcl)<<endl;
            }
        }


        //////////////////
        //// GEN TAUS ////
        //////////////////

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


        ////////////////////////
        //// MATCHING TO TC ////
        ////////////////////////

        /*for(int i_gentau=0; i_gentau<_gentau_n;i_gentau++){

            TLorentzVector gentauvis;
            gentauvis.SetPtEtaPhiM( (*_gentau_vis_pt)[i_gentau], (*_gentau_vis_eta)[i_gentau], (*_gentau_vis_phi)[i_gentau], (*_gentau_vis_mass)[i_gentau]);

            bool  tau_isMatchedtoSTC = false;

            int   tau_n_tc = 0;
            
            float tau_pt_tot = -999;
            float tau_pt_seed = -999;
            float tau_eta_Eweighted = -999;
            float tau_eta_seed = -999;
            float tau_phi_Eweighted = -999;
            float tau_phi_seed = -999;
            float tau_x_seed = -999;
            float tau_x_Eweighted = -999;
            float tau_y_seed = -999;
            float tau_y_Eweighted = -999;
            float tau_z_seed = -999;
            float tau_z_Eweighted = -999;

            for(int i_supercl=0; i_supercl<_n_STC; i_supercl++){

                TLorentzVector STC;
                STC.SetPtEtaPhiM(_STC_pt_tot.at(i_supercl),_STC_eta_Eweighted.at(i_supercl),_STC_phi_Eweighted.at(i_supercl),0);

                float dR = STC.DeltaR(gentauvis);

                tau_isMatchedtoSTC = (dR <= dRmax_match);

                if( tau_isMatchedtoSTC ) {                 

                    tau_n_tc = _STC_n_tc.at(i_supercl);

                    tau_pt_tot = _STC_pt_tot.at(i_supercl);
                    tau_pt_seed = _STC_pt_seed.at(i_supercl);
                    tau_eta_Eweighted = _STC_eta_Eweighted.at(i_supercl);
                    tau_eta_seed = _STC_eta_seed.at(i_supercl);
                    tau_phi_Eweighted = _STC_phi_Eweighted.at(i_supercl);
                    tau_phi_seed = _STC_phi_seed.at(i_supercl);
                    tau_x_Eweighted = _STC_x_Eweighted.at(i_supercl);
                    tau_x_seed = _STC_x_seed.at(i_supercl);
                    tau_y_Eweighted = _STC_y_Eweighted.at(i_supercl);
                    tau_y_seed = _STC_y_seed.at(i_supercl);
                    tau_z_Eweighted = _STC_z_Eweighted.at(i_supercl);
                    tau_z_seed = _STC_z_seed.at(i_supercl);

                    if(debug){
                        cout<<" --> Matched: gen tau eta "<<gentauvis.Eta()<<" to STC eta "<<STC.Eta()<<endl;
                        cout<<"       SC n tc: "<<tau_n_tc<<endl;
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

            _gentau_isMatchedtoSTC.push_back(tau_isMatchedtoSTC);

            _gentau_matchedSTC_n_tc.push_back(tau_n_tc);
            _gentau_matchedSTC_pt_tot.push_back(tau_pt_tot);
            _gentau_matchedSTC_pt_seed.push_back(tau_pt_seed);
            _gentau_matchedSTC_eta_Eweighted.push_back(tau_eta_Eweighted);
            _gentau_matchedSTC_eta_seed.push_back(tau_eta_seed);
            _gentau_matchedSTC_phi_Eweighted.push_back(tau_phi_Eweighted);
            _gentau_matchedSTC_phi_seed.push_back(tau_phi_seed);
            _gentau_matchedSTC_x_Eweighted.push_back(tau_x_Eweighted);
            _gentau_matchedSTC_x_seed.push_back(tau_x_seed);
            _gentau_matchedSTC_y_Eweighted.push_back(tau_y_Eweighted);
            _gentau_matchedSTC_y_seed.push_back(tau_y_seed);
            _gentau_matchedSTC_z_Eweighted.push_back(tau_z_Eweighted);
            _gentau_matchedSTC_z_seed.push_back(tau_z_seed);

        }

        for(int i_gentau=0; i_gentau<_gentau_n;i_gentau++){

            //bool ishadronic = ( ((*_gentau_decayMode)[i_gentau]== 0) || ((*_gentau_decayMode)[i_gentau]== 1) || ((*_gentau_decayMode)[i_gentau]== 4) || ((*_gentau_decayMode)[i_gentau]== 5)) ;
            bool matched = _gentau_isMatchedtoSTC[i_gentau];
            if(matched) matchedtoSTC_taus += 1;
            else if(!matched) unmatchedtoSTC_taus += 1;

        }

        if(debug) cout<<" .................................................. "<<endl;
        */

        //////////////////////////
        //// MATCHING TO CL3D ////
        //////////////////////////

        for(int i_gentau=0; i_gentau<_gentau_n;i_gentau++){

            TLorentzVector gentauvis;
            gentauvis.SetPtEtaPhiM( (*_gentau_vis_pt)[i_gentau], (*_gentau_vis_eta)[i_gentau], (*_gentau_vis_phi)[i_gentau], (*_gentau_vis_mass)[i_gentau]);

            bool  tau_isMatchedtoSCL3D = false;

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

            for(int i_supercl=0; i_supercl<_n_SCL3D; i_supercl++){

                TLorentzVector SCL3D;
                SCL3D.SetPtEtaPhiM(_SCL3D_pt_tot.at(i_supercl),_SCL3D_eta_Eweighted.at(i_supercl),_SCL3D_phi_Eweighted.at(i_supercl),0);

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

                cl3ds_showerlength = _SCL3D_cl3d_showerlength.at(i_supercl);
                cl3ds_coreshowerlength = _SCL3D_cl3d_coreshowerlength.at(i_supercl);
                cl3ds_firstlayer = _SCL3D_cl3d_firstlayer.at(i_supercl);
                cl3ds_maxlayer = _SCL3D_cl3d_maxlayer.at(i_supercl);
                cl3ds_seetot = _SCL3D_cl3d_seetot.at(i_supercl);
                cl3ds_seemax = _SCL3D_cl3d_seemax.at(i_supercl);
                cl3ds_spptot = _SCL3D_cl3d_spptot.at(i_supercl);
                cl3ds_sppmax = _SCL3D_cl3d_sppmax.at(i_supercl);
                cl3ds_szz = _SCL3D_cl3d_szz.at(i_supercl);
                cl3ds_srrtot = _SCL3D_cl3d_srrtot.at(i_supercl);
                cl3ds_srrmax = _SCL3D_cl3d_srrmax.at(i_supercl);
                cl3ds_srrmean = _SCL3D_cl3d_srrmean.at(i_supercl);
                cl3ds_emaxe = _SCL3D_cl3d_emaxe.at(i_supercl);
                cl3ds_hoe = _SCL3D_cl3d_hoe.at(i_supercl);
                cl3ds_meanz = _SCL3D_cl3d_meanz.at(i_supercl);
                cl3ds_layer10 = _SCL3D_cl3d_layer10.at(i_supercl);
                cl3ds_layer50 = _SCL3D_cl3d_layer50.at(i_supercl);
                cl3ds_layer90 = _SCL3D_cl3d_layer90.at(i_supercl);
                cl3ds_ntc67 = _SCL3D_cl3d_ntc67.at(i_supercl);
                cl3ds_ntc90 = _SCL3D_cl3d_ntc90.at(i_supercl);
                cl3ds_bdteg = _SCL3D_cl3d_bdteg.at(i_supercl);
                cl3ds_quality = _SCL3D_cl3d_quality.at(i_supercl);
                cl3ds_puBDT = _SCL3D_cl3d_puBDT.at(i_supercl);

                float dR = SCL3D.DeltaR(gentauvis);

                tau_isMatchedtoSCL3D = (dR <= dRmax_match);

                if( tau_isMatchedtoSCL3D ) {                 

                    tau_n_cl3d = _SCL3D_n_cl3d.at(i_supercl);

                    tau_pt_tot = _SCL3D_pt_tot.at(i_supercl);
                    tau_pt_seed = _SCL3D_pt_seed.at(i_supercl);
                    tau_eta_Eweighted = _SCL3D_eta_Eweighted.at(i_supercl);
                    tau_eta_seed = _SCL3D_eta_seed.at(i_supercl);
                    tau_phi_Eweighted = _SCL3D_phi_Eweighted.at(i_supercl);
                    tau_phi_seed = _SCL3D_phi_seed.at(i_supercl);

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
                        cout<<" --> Matched: gen tau eta "<<gentauvis.Eta()<<" to SCL3D eta "<<SCL3D.Eta()<<endl;
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

            _gentau_isMatchedtoSCL3D.push_back(tau_isMatchedtoSCL3D);

            _gentau_matchedSCL3D_n_cl3d.push_back(tau_n_cl3d);
            _gentau_matchedSCL3D_pt_tot.push_back(tau_pt_tot);
            _gentau_matchedSCL3D_pt_seed.push_back(tau_pt_seed);
            _gentau_matchedSCL3D_eta_Eweighted.push_back(tau_eta_Eweighted);
            _gentau_matchedSCL3D_eta_seed.push_back(tau_eta_seed);
            _gentau_matchedSCL3D_phi_Eweighted.push_back(tau_phi_Eweighted);
            _gentau_matchedSCL3D_phi_seed.push_back(tau_phi_seed);

            _gentau_matchedSCL3D_showerlength_seed.push_back(tau_showerlength_seed);
            _gentau_matchedSCL3D_coreshowerlength_seed.push_back(tau_coreshowerlength_seed);
            _gentau_matchedSCL3D_firstlayer_seed.push_back(tau_firstlayer_seed);
            _gentau_matchedSCL3D_maxlayer_seed.push_back(tau_maxlayer_seed);
            _gentau_matchedSCL3D_seetot_seed.push_back(tau_seetot_seed);
            _gentau_matchedSCL3D_seemax_seed.push_back(tau_seemax_seed);
            _gentau_matchedSCL3D_spptot_seed.push_back(tau_spptot_seed);
            _gentau_matchedSCL3D_sppmax_seed.push_back(tau_sppmax_seed);
            _gentau_matchedSCL3D_szz_seed.push_back(tau_szz_seed);
            _gentau_matchedSCL3D_srrtot_seed.push_back(tau_srrtot_seed);
            _gentau_matchedSCL3D_srrmax_seed.push_back(tau_srrmax_seed);
            _gentau_matchedSCL3D_srrmean_seed.push_back(tau_srrmean_seed);
            _gentau_matchedSCL3D_emaxe_seed.push_back(tau_emaxe_seed);
            _gentau_matchedSCL3D_hoe_seed.push_back(tau_hoe_seed);
            _gentau_matchedSCL3D_meanz_seed.push_back(tau_meanz_seed);
            _gentau_matchedSCL3D_layer10_seed.push_back(tau_layer10_seed);
            _gentau_matchedSCL3D_layer50_seed.push_back(tau_layer50_seed);
            _gentau_matchedSCL3D_layer90_seed.push_back(tau_layer90_seed);
            _gentau_matchedSCL3D_ntc67_seed.push_back(tau_ntc67_seed);
            _gentau_matchedSCL3D_ntc90_seed.push_back(tau_ntc90_seed);
            _gentau_matchedSCL3D_bdteg_seed.push_back(tau_bdteg_seed);
            _gentau_matchedSCL3D_quality_seed.push_back(tau_quality_seed);
            _gentau_matchedSCL3D_puBDT_seed.push_back(tau_puBDT_seed);

        }

        for(int i_gentau=0; i_gentau<_gentau_n;i_gentau++){

            //bool ishadronic = ( ((*_gentau_decayMode)[i_gentau]== 0) || ((*_gentau_decayMode)[i_gentau]== 1) || ((*_gentau_decayMode)[i_gentau]== 4) || ((*_gentau_decayMode)[i_gentau]== 5)) ;
            bool matched = _gentau_isMatchedtoSCL3D[i_gentau];
            if(matched) matchedtoSCL3D_taus += 1;
            else if(!matched) unmatchedtoSCL3D_taus += 1;

        }

        if(debug) cout<<" .................................................. "<<endl;

        out_tree->Fill();

    }

    out_file->cd();

    out_tree->Write();
    out_file->Close();

    //Total efficiencies

    //float efficiency_STC   = matchedtoSTC_taus/(matchedtoSTC_taus+unmatchedtoSTC_taus);
    float efficiency_SCL3D = matchedtoSCL3D_taus/(matchedtoSCL3D_taus+unmatchedtoSCL3D_taus);

    cout<<" "<<endl;
    cout<<" =========================================== "<<endl;
    //cout<<" STC:    Matched taus: "<<matchedtoSTC_taus<<"; unmatched taus: "<<unmatchedtoSTC_taus<<" --> matching efficiency "<<efficiency_STC*100<<"%"<<endl;
    cout<<" SCL3DL: Matched taus: "<<matchedtoSCL3D_taus<<"; unmatched taus: "<<unmatchedtoSCL3D_taus<<" --> matching efficiency "<<efficiency_SCL3D*100<<"%"<<endl;
    cout<<" =========================================== "<<endl;


    return;
}


void test(int n_events = -1, bool debugging = false){

  //TString infile = "/data_CMS/cms/mperez/HGCal_data/Aug19/skimmed/ntuple_RelValDiTau_Pt20To100_Eta1p6To2p9_v10_PU0_skimmed.root";
  //TString outfile = "/data_CMS/cms/mperez/HGCal_data/Aug19/clustered/ntuple_RelValDiTau_Pt20To100_Eta1p6To2p9_v10_PU0_clustered_PUcut.root";

  //TString infile = "/data_CMS/cms/mperez/HGCal_data/Aug19/skimmed/ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_skimmed.root";
  //TString outfile = "/data_CMS/cms/mperez/HGCal_data/Aug19/clustered/ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_clustered_PUcut.root";

  TString infile = "/data_CMS/cms/mperez/HGCal_data/Aug19/skimmed/ntuple_RelValDiTau_Pt20To100_Eta1p6To2p9_v10_PU200_skimmed.root";
  TString outfile = "/data_CMS/cms/mperez/HGCal_data/Aug19/clustered/ntuple_RelValDiTau_Pt20To100_Eta1p6To2p9_v10_PU200_clustered_PUcut.root";

  //TString infile = "/data_CMS/cms/mperez/HGCal_data/Aug19/skimmed/ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU200_skimmed.root";
  //TString outfile = "/data_CMS/cms/mperez/HGCal_data/Aug19/clustered/ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU200_clustered_noPUcut.root";

  //TString infile = "/data_CMS/cms/mperez/HGCal_data/Aug19/skimmed/ntuple_Nu_E10_v10_PU200_files100to150_skimmed.root";
  //TString outfile = "/data_CMS/cms/mperez/HGCal_data/Aug19/clustered/ntuple_Nu_E10_v10_PU200_files100to150_clustered_noPUcut.root";

  cluster_tree(infile, outfile, n_events, debugging);

}
