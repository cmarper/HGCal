/////////////////////////////////////////////////////////
///// HGCal L1 taus, C. Martin Perez, LLR, Jul 2019 /////
/////////////////////////////////////////////////////////

#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <iostream>

bool pT_comparison_pairs(pair<int,TLorentzVector> pair1, pair<int,TLorentzVector> pair2){

  return (pair1.second).Pt()>(pair2.second).Pt();

}

void cluster_tree( TString filein, TString fileout, int nevents = -1, float pt_thr_seed = 4, float pt_thr_sec = 2, float eta_window_maxi = 0.15, float phi_window_maxi = 0.45, float eta_window_mini = 0.15, float phi_window_mini = 0.15){

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
    vector<float>  *_cl3d_bdteg;
    vector<int>    *_cl3d_showerlength;
    vector<int>    *_cl3d_coreshowerlength;
    vector<int>    *_cl3d_firstlayer;
    vector<int>    *_cl3d_maxlayer;
    vector<float>  *_cl3d_seetot;
    vector<float>  *_cl3d_spptot;
    vector<float>  *_cl3d_szz;
    vector<float>  *_cl3d_srrtot;
    vector<float>  *_cl3d_srrmean;


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
    in_tree->SetBranchAddress("cl3d_bdteg",            &_cl3d_bdteg);
    in_tree->SetBranchAddress("cl3d_showerlength",     &_cl3d_showerlength);
    in_tree->SetBranchAddress("cl3d_coreshowerlength", &_cl3d_coreshowerlength);
    in_tree->SetBranchAddress("cl3d_firstlayer",       &_cl3d_firstlayer);
    in_tree->SetBranchAddress("cl3d_maxlayer",         &_cl3d_maxlayer);
    in_tree->SetBranchAddress("cl3d_seetot",           &_cl3d_seetot);
    in_tree->SetBranchAddress("cl3d_spptot",           &_cl3d_spptot);
    in_tree->SetBranchAddress("cl3d_szz",              &_cl3d_szz);
    in_tree->SetBranchAddress("cl3d_srrtot",           &_cl3d_srrtot);
    in_tree->SetBranchAddress("cl3d_srrmean",          &_cl3d_srrmean);


    TTree* out_tree=in_tree->GetTree()->CloneTree(0);
    out_tree->SetNameTitle("ClusteredTree","ClusteredTree");

    int _maxicluster_n;

    vector<vector<float>>  _maxicluster_cl3d_pt;
    vector<vector<float>>  _maxicluster_cl3d_energy;
    vector<vector<float>>  _maxicluster_cl3d_eta;
    vector<vector<float>>  _maxicluster_cl3d_phi;
    vector<vector<float>>  _maxicluster_cl3d_bdteg;
    vector<vector<int>>    _maxicluster_cl3d_showerlength;
    vector<vector<int>>    _maxicluster_cl3d_coreshowerlength;
    vector<vector<int>>    _maxicluster_cl3d_firstlayer;
    vector<vector<int>>    _maxicluster_cl3d_maxlayer;
    vector<vector<float>>  _maxicluster_cl3d_seetot;
    vector<vector<float>>  _maxicluster_cl3d_spptot;
    vector<vector<float>>  _maxicluster_cl3d_szz;
    vector<vector<float>>  _maxicluster_cl3d_srrtot;
    vector<vector<float>>  _maxicluster_cl3d_srrmean;

    int _minicluster_n;

    vector<vector<float>>  _minicluster_cl3d_pt;
    vector<vector<float>>  _minicluster_cl3d_energy;
    vector<vector<float>>  _minicluster_cl3d_eta;
    vector<vector<float>>  _minicluster_cl3d_phi;
    vector<vector<float>>  _minicluster_cl3d_bdteg;
    vector<vector<int>>    _minicluster_cl3d_showerlength;
    vector<vector<int>>    _minicluster_cl3d_coreshowerlength;
    vector<vector<int>>    _minicluster_cl3d_firstlayer;
    vector<vector<int>>    _minicluster_cl3d_maxlayer;
    vector<vector<float>>  _minicluster_cl3d_seetot;
    vector<vector<float>>  _minicluster_cl3d_spptot;
    vector<vector<float>>  _minicluster_cl3d_szz;
    vector<vector<float>>  _minicluster_cl3d_srrtot;
    vector<vector<float>>  _minicluster_cl3d_srrmean;

    int _supercluster_n;

    vector<bool> _supercluster_isMerged;

    vector<vector<float>>  _supercluster_cl3d_pt;
    vector<vector<float>>  _supercluster_cl3d_energy;
    vector<vector<float>>  _supercluster_cl3d_eta;
    vector<vector<float>>  _supercluster_cl3d_phi;
    vector<vector<float>>  _supercluster_cl3d_bdteg;
    vector<vector<int>>    _supercluster_cl3d_showerlength;
    vector<vector<int>>    _supercluster_cl3d_coreshowerlength;
    vector<vector<int>>    _supercluster_cl3d_firstlayer;
    vector<vector<int>>    _supercluster_cl3d_maxlayer;
    vector<vector<float>>  _supercluster_cl3d_seetot;
    vector<vector<float>>  _supercluster_cl3d_spptot;
    vector<vector<float>>  _supercluster_cl3d_szz;
    vector<vector<float>>  _supercluster_cl3d_srrtot;
    vector<vector<float>>  _supercluster_cl3d_srrmean;

    out_tree->Branch("maxicluster_n",     &_maxicluster_n);

    out_tree->Branch("maxicluster_cl3d_pt",                 &_maxicluster_cl3d_pt);
    out_tree->Branch("maxicluster_cl3d_energy",             &_maxicluster_cl3d_energy);
    out_tree->Branch("maxicluster_cl3d_eta",                &_maxicluster_cl3d_eta);
    out_tree->Branch("maxicluster_cl3d_phi",                &_maxicluster_cl3d_phi);
    out_tree->Branch("maxicluster_cl3d_bdteg",              &_maxicluster_cl3d_bdteg);
    out_tree->Branch("maxicluster_cl3d_showerlength",       &_maxicluster_cl3d_showerlength);
    out_tree->Branch("maxicluster_cl3d_coreshowerlength",   &_maxicluster_cl3d_coreshowerlength);
    out_tree->Branch("maxicluster_cl3d_firstlayer",         &_maxicluster_cl3d_firstlayer);
    out_tree->Branch("maxicluster_cl3d_maxlayer",           &_maxicluster_cl3d_maxlayer);
    out_tree->Branch("maxicluster_cl3d_seetot",             &_maxicluster_cl3d_seetot);
    out_tree->Branch("maxicluster_cl3d_spptot",             &_maxicluster_cl3d_spptot);
    out_tree->Branch("maxicluster_cl3d_szz",                &_maxicluster_cl3d_szz);
    out_tree->Branch("maxicluster_cl3d_srrtot",             &_maxicluster_cl3d_srrtot);
    out_tree->Branch("maxicluster_cl3d_srrmean",            &_maxicluster_cl3d_srrmean);

    out_tree->Branch("minicluster_n",     &_minicluster_n);

    out_tree->Branch("minicluster_cl3d_pt",                 &_minicluster_cl3d_pt);
    out_tree->Branch("minicluster_cl3d_energy",             &_minicluster_cl3d_energy);
    out_tree->Branch("minicluster_cl3d_eta",                &_minicluster_cl3d_eta);
    out_tree->Branch("minicluster_cl3d_phi",                &_minicluster_cl3d_phi);
    out_tree->Branch("minicluster_cl3d_bdteg",              &_minicluster_cl3d_bdteg);
    out_tree->Branch("minicluster_cl3d_showerlength",       &_minicluster_cl3d_showerlength);
    out_tree->Branch("minicluster_cl3d_coreshowerlength",   &_minicluster_cl3d_coreshowerlength);
    out_tree->Branch("minicluster_cl3d_firstlayer",         &_minicluster_cl3d_firstlayer);
    out_tree->Branch("minicluster_cl3d_maxlayer",           &_minicluster_cl3d_maxlayer);
    out_tree->Branch("minicluster_cl3d_seetot",             &_minicluster_cl3d_seetot);
    out_tree->Branch("minicluster_cl3d_spptot",             &_minicluster_cl3d_spptot);
    out_tree->Branch("minicluster_cl3d_szz",                &_minicluster_cl3d_szz);
    out_tree->Branch("minicluster_cl3d_srrtot",             &_minicluster_cl3d_srrtot);
    out_tree->Branch("minicluster_cl3d_srrmean",            &_minicluster_cl3d_srrmean);

    out_tree->Branch("supercluster_n",     &_supercluster_n);

    out_tree->Branch("supercluster_isMerged", &_supercluster_isMerged);

    out_tree->Branch("supercluster_cl3d_pt",                 &_supercluster_cl3d_pt);
    out_tree->Branch("supercluster_cl3d_energy",             &_supercluster_cl3d_energy);
    out_tree->Branch("supercluster_cl3d_eta",                &_supercluster_cl3d_eta);
    out_tree->Branch("supercluster_cl3d_phi",                &_supercluster_cl3d_phi);
    out_tree->Branch("supercluster_cl3d_bdteg",              &_supercluster_cl3d_bdteg);
    out_tree->Branch("supercluster_cl3d_showerlength",       &_supercluster_cl3d_showerlength);
    out_tree->Branch("supercluster_cl3d_coreshowerlength",   &_supercluster_cl3d_coreshowerlength);
    out_tree->Branch("supercluster_cl3d_firstlayer",         &_supercluster_cl3d_firstlayer);
    out_tree->Branch("supercluster_cl3d_maxlayer",           &_supercluster_cl3d_maxlayer);
    out_tree->Branch("supercluster_cl3d_seetot",             &_supercluster_cl3d_seetot);
    out_tree->Branch("supercluster_cl3d_spptot",             &_supercluster_cl3d_spptot);
    out_tree->Branch("supercluster_cl3d_szz",                &_supercluster_cl3d_szz);
    out_tree->Branch("supercluster_cl3d_srrtot",             &_supercluster_cl3d_srrtot);
    out_tree->Branch("supercluster_cl3d_srrmean",            &_supercluster_cl3d_srrmean);


    for (int i=0;i<nentries;i++) {

        if(i%1000==0) cout<<"i="<<i<<endl;

        // old branches

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
        _cl3d_bdteg = 0;
        _cl3d_showerlength = 0;
        _cl3d_coreshowerlength = 0;
        _cl3d_firstlayer = 0;
        _cl3d_maxlayer = 0;
        _cl3d_seetot = 0;
        _cl3d_spptot = 0;
        _cl3d_szz = 0;
        _cl3d_srrtot = 0;
        _cl3d_srrmean = 0;

        _maxicluster_n = 0;

        _maxicluster_cl3d_pt.clear();
        _maxicluster_cl3d_energy.clear();
        _maxicluster_cl3d_eta.clear();
        _maxicluster_cl3d_phi.clear();
        _maxicluster_cl3d_bdteg.clear();
        _maxicluster_cl3d_showerlength.clear();
        _maxicluster_cl3d_coreshowerlength.clear();
        _maxicluster_cl3d_firstlayer.clear();
        _maxicluster_cl3d_maxlayer.clear();
        _maxicluster_cl3d_seetot.clear();
        _maxicluster_cl3d_spptot.clear();
        _maxicluster_cl3d_szz.clear();
        _maxicluster_cl3d_srrtot.clear();
        _maxicluster_cl3d_srrmean.clear();

        _minicluster_n = 0;

        _minicluster_cl3d_pt.clear();
        _minicluster_cl3d_energy.clear();
        _minicluster_cl3d_eta.clear();
        _minicluster_cl3d_phi.clear();
        _minicluster_cl3d_bdteg.clear();
        _minicluster_cl3d_showerlength.clear();
        _minicluster_cl3d_coreshowerlength.clear();
        _minicluster_cl3d_firstlayer.clear();
        _minicluster_cl3d_maxlayer.clear();
        _minicluster_cl3d_seetot.clear();
        _minicluster_cl3d_spptot.clear();
        _minicluster_cl3d_szz.clear();
        _minicluster_cl3d_srrtot.clear();
        _minicluster_cl3d_srrmean.clear();

        _supercluster_n = 0;

        _supercluster_isMerged.clear();

        _supercluster_cl3d_pt.clear();
        _supercluster_cl3d_energy.clear();
        _supercluster_cl3d_eta.clear();
        _supercluster_cl3d_phi.clear();
        _supercluster_cl3d_bdteg.clear();
        _supercluster_cl3d_showerlength.clear();
        _supercluster_cl3d_coreshowerlength.clear();
        _supercluster_cl3d_firstlayer.clear();
        _supercluster_cl3d_maxlayer.clear();
        _supercluster_cl3d_seetot.clear();
        _supercluster_cl3d_spptot.clear();
        _supercluster_cl3d_szz.clear();
        _supercluster_cl3d_srrtot.clear();
        _supercluster_cl3d_srrmean.clear();


        int entry_ok = in_tree->GetEntry(i);    
        if(entry_ok<0) 
            continue;

        float previous_pt_leading_maxi = -999.;
        float previous_pt_leading_mini = -999.;


        //////////////////////////////////////////
        ////////////// CLUSTERING ////////////////
        //////////////////////////////////////////

        // MAXICLUSTERS

        for (int i_main=0; i_main<_cl3d_n; i_main++){

            if ( (*_cl3d_pt)[i_main] < 0.0001 || (*_cl3d_pt)[i_main] > 10000 ) continue;

            // RAW CLUSTERING 

            vector<TLorentzVector> candidate_maxicluster;
            vector<float> candidate_maxicluster_bdteg;
            vector<int>   candidate_maxicluster_showerlength;
            vector<int>   candidate_maxicluster_coreshowerlength;
            vector<int>   candidate_maxicluster_firstlayer;
            vector<int>   candidate_maxicluster_maxlayer;
            vector<float> candidate_maxicluster_seetot;
            vector<float> candidate_maxicluster_spptot;
            vector<float> candidate_maxicluster_szz;
            vector<float> candidate_maxicluster_srrtot;
            vector<float> candidate_maxicluster_srrmean;

            candidate_maxicluster.clear();
            candidate_maxicluster_bdteg.clear();
            candidate_maxicluster_showerlength.clear();
            candidate_maxicluster_coreshowerlength.clear();
            candidate_maxicluster_firstlayer.clear();
            candidate_maxicluster_maxlayer.clear();
            candidate_maxicluster_seetot.clear();
            candidate_maxicluster_spptot.clear();
            candidate_maxicluster_szz.clear();
            candidate_maxicluster_srrtot.clear();
            candidate_maxicluster_srrmean.clear();


            // MAIN CLUSTER CANDIDATE (SUPERCLUSTER SEED)

            TLorentzVector seed_candidate_maxicluster;
            float seed_candidate_maxicluster_bdteg;
            int   seed_candidate_maxicluster_showerlength;
            int   seed_candidate_maxicluster_coreshowerlength;
            int   seed_candidate_maxicluster_firstlayer;
            int   seed_candidate_maxicluster_maxlayer;
            float seed_candidate_maxicluster_seetot;
            float seed_candidate_maxicluster_spptot;
            float seed_candidate_maxicluster_szz;
            float seed_candidate_maxicluster_srrtot;
            float seed_candidate_maxicluster_srrmean;

            seed_candidate_maxicluster.SetPtEtaPhiM( (*_cl3d_pt)[i_main], (*_cl3d_eta)[i_main], (*_cl3d_phi)[i_main], 0);
            seed_candidate_maxicluster_bdteg = (*_cl3d_bdteg)[i_main];
            seed_candidate_maxicluster_showerlength = (*_cl3d_showerlength)[i_main];
            seed_candidate_maxicluster_coreshowerlength = (*_cl3d_coreshowerlength)[i_main];
            seed_candidate_maxicluster_firstlayer = (*_cl3d_firstlayer)[i_main];
            seed_candidate_maxicluster_maxlayer = (*_cl3d_maxlayer)[i_main];
            seed_candidate_maxicluster_seetot = (*_cl3d_seetot)[i_main];
            seed_candidate_maxicluster_spptot = (*_cl3d_spptot)[i_main];
            seed_candidate_maxicluster_szz = (*_cl3d_szz)[i_main];
            seed_candidate_maxicluster_srrtot = (*_cl3d_srrtot)[i_main];
            seed_candidate_maxicluster_srrmean = (*_cl3d_srrmean)[i_main];


            if ( seed_candidate_maxicluster.Pt() < pt_thr_seed ) continue;

            candidate_maxicluster.push_back(seed_candidate_maxicluster);
            candidate_maxicluster_bdteg.push_back(seed_candidate_maxicluster_bdteg);
            candidate_maxicluster_showerlength.push_back(seed_candidate_maxicluster_showerlength);
            candidate_maxicluster_coreshowerlength.push_back(seed_candidate_maxicluster_coreshowerlength);
            candidate_maxicluster_firstlayer.push_back(seed_candidate_maxicluster_firstlayer);
            candidate_maxicluster_maxlayer.push_back(seed_candidate_maxicluster_maxlayer);
            candidate_maxicluster_seetot.push_back(seed_candidate_maxicluster_seetot);
            candidate_maxicluster_spptot.push_back(seed_candidate_maxicluster_spptot);
            candidate_maxicluster_szz.push_back(seed_candidate_maxicluster_szz);
            candidate_maxicluster_srrtot.push_back(seed_candidate_maxicluster_srrtot);
            candidate_maxicluster_srrmean.push_back(seed_candidate_maxicluster_srrmean);

            // SECONDARY CLUSTERS CANDIDATES

            for (int i_sec=0; i_sec<_cl3d_n; i_sec++){

                if( i_sec==i_main ) continue;

                if ( (*_cl3d_pt)[i_sec] < 0.0001 || (*_cl3d_pt)[i_sec] > 10000 ) continue;

                TLorentzVector secondary_candidate_maxicluster;
                float secondary_candidate_maxicluster_bdteg;
                int   secondary_candidate_maxicluster_showerlength;
                int   secondary_candidate_maxicluster_coreshowerlength;
                int   secondary_candidate_maxicluster_firstlayer;
                int   secondary_candidate_maxicluster_maxlayer;
                float secondary_candidate_maxicluster_seetot;
                float secondary_candidate_maxicluster_spptot;
                float secondary_candidate_maxicluster_szz;
                float secondary_candidate_maxicluster_srrtot;
                float secondary_candidate_maxicluster_srrmean;

                secondary_candidate_maxicluster.SetPtEtaPhiM( (*_cl3d_pt)[i_sec], (*_cl3d_eta)[i_sec], (*_cl3d_phi)[i_sec], 0);
                secondary_candidate_maxicluster_bdteg = (*_cl3d_bdteg)[i_sec];
                secondary_candidate_maxicluster_showerlength = (*_cl3d_showerlength)[i_sec];
                secondary_candidate_maxicluster_coreshowerlength = (*_cl3d_coreshowerlength)[i_sec];
                secondary_candidate_maxicluster_firstlayer = (*_cl3d_firstlayer)[i_sec];
                secondary_candidate_maxicluster_maxlayer = (*_cl3d_maxlayer)[i_sec];
                secondary_candidate_maxicluster_seetot = (*_cl3d_seetot)[i_sec];
                secondary_candidate_maxicluster_spptot = (*_cl3d_spptot)[i_sec];
                secondary_candidate_maxicluster_szz = (*_cl3d_szz)[i_sec];
                secondary_candidate_maxicluster_srrtot = (*_cl3d_srrtot)[i_sec];
                secondary_candidate_maxicluster_srrmean = (*_cl3d_srrmean)[i_sec];

                if ( secondary_candidate_maxicluster.Pt() < pt_thr_sec ) continue;

                if( abs(secondary_candidate_maxicluster.Eta() - seed_candidate_maxicluster.Eta()) > eta_window_maxi ) continue;
                if( abs(secondary_candidate_maxicluster.Phi() - seed_candidate_maxicluster.Phi()) > phi_window_maxi ) continue;

                candidate_maxicluster.push_back(secondary_candidate_maxicluster);
                candidate_maxicluster_bdteg.push_back(secondary_candidate_maxicluster_bdteg);
                candidate_maxicluster_showerlength.push_back(secondary_candidate_maxicluster_showerlength);
                candidate_maxicluster_coreshowerlength.push_back(secondary_candidate_maxicluster_coreshowerlength);
                candidate_maxicluster_firstlayer.push_back(secondary_candidate_maxicluster_firstlayer);
                candidate_maxicluster_maxlayer.push_back(secondary_candidate_maxicluster_maxlayer);
                candidate_maxicluster_seetot.push_back(secondary_candidate_maxicluster_seetot);
                candidate_maxicluster_spptot.push_back(secondary_candidate_maxicluster_spptot);
                candidate_maxicluster_szz.push_back(secondary_candidate_maxicluster_szz);
                candidate_maxicluster_srrtot.push_back(secondary_candidate_maxicluster_srrtot);
                candidate_maxicluster_srrmean.push_back(secondary_candidate_maxicluster_srrmean);

            }

            // SORT CLUSTERS BY PT IN SUPERCLUSTER

            vector<TLorentzVector> candidate_maxicluster_pTsorted;
            vector<float> candidate_maxicluster_pt_pTsorted;
            vector<float> candidate_maxicluster_E_pTsorted;
            vector<float> candidate_maxicluster_eta_pTsorted;
            vector<float> candidate_maxicluster_phi_pTsorted;
            vector<float> candidate_maxicluster_bdteg_pTsorted;
            vector<int>   candidate_maxicluster_showerlength_pTsorted;
            vector<int>   candidate_maxicluster_coreshowerlength_pTsorted;
            vector<int>   candidate_maxicluster_firstlayer_pTsorted;
            vector<int>   candidate_maxicluster_maxlayer_pTsorted;
            vector<float> candidate_maxicluster_seetot_pTsorted;
            vector<float> candidate_maxicluster_spptot_pTsorted;
            vector<float> candidate_maxicluster_szz_pTsorted;
            vector<float> candidate_maxicluster_srrtot_pTsorted;
            vector<float> candidate_maxicluster_srrmean_pTsorted;

            candidate_maxicluster_pTsorted.clear();
            candidate_maxicluster_pt_pTsorted.clear();
            candidate_maxicluster_E_pTsorted.clear();
            candidate_maxicluster_eta_pTsorted.clear();
            candidate_maxicluster_phi_pTsorted.clear();
            candidate_maxicluster_bdteg_pTsorted.clear();
            candidate_maxicluster_showerlength_pTsorted.clear();
            candidate_maxicluster_coreshowerlength_pTsorted.clear();
            candidate_maxicluster_firstlayer_pTsorted.clear();
            candidate_maxicluster_maxlayer_pTsorted.clear();
            candidate_maxicluster_seetot_pTsorted.clear();
            candidate_maxicluster_spptot_pTsorted.clear();
            candidate_maxicluster_szz_pTsorted.clear();
            candidate_maxicluster_srrtot_pTsorted.clear();
            candidate_maxicluster_srrmean_pTsorted.clear();

            vector< pair<int,TLorentzVector> > imaxicluster_maxicluster_pairs;

            for (unsigned int i_cluster = 0; i_cluster<candidate_maxicluster.size(); i_cluster++){

                pair<int,TLorentzVector> maxicluster_pair = make_pair(i_cluster,candidate_maxicluster.at(i_cluster));
                imaxicluster_maxicluster_pairs.push_back(maxicluster_pair);

            }

            sort(imaxicluster_maxicluster_pairs.begin(), imaxicluster_maxicluster_pairs.end(), pT_comparison_pairs);

            for (unsigned int i_cluster = 0; i_cluster<candidate_maxicluster.size(); i_cluster++){

                int index = imaxicluster_maxicluster_pairs[i_cluster].first;
                TLorentzVector maxiclus = imaxicluster_maxicluster_pairs[i_cluster].second;

                candidate_maxicluster_pTsorted.push_back(maxiclus);
                candidate_maxicluster_pt_pTsorted.push_back(maxiclus.Pt());
                candidate_maxicluster_E_pTsorted.push_back(maxiclus.E());
                candidate_maxicluster_eta_pTsorted.push_back(maxiclus.Eta());
                candidate_maxicluster_phi_pTsorted.push_back(maxiclus.Phi());
                candidate_maxicluster_bdteg_pTsorted.push_back(candidate_maxicluster_bdteg.at(index));
                candidate_maxicluster_showerlength_pTsorted.push_back(candidate_maxicluster_showerlength.at(index));
                candidate_maxicluster_coreshowerlength_pTsorted.push_back(candidate_maxicluster_coreshowerlength.at(index));
                candidate_maxicluster_firstlayer_pTsorted.push_back(candidate_maxicluster_firstlayer.at(index));
                candidate_maxicluster_maxlayer_pTsorted.push_back(candidate_maxicluster_maxlayer.at(index));
                candidate_maxicluster_seetot_pTsorted.push_back(candidate_maxicluster_seetot.at(index));
                candidate_maxicluster_spptot_pTsorted.push_back(candidate_maxicluster_spptot.at(index));
                candidate_maxicluster_szz_pTsorted.push_back(candidate_maxicluster_szz.at(index));
                candidate_maxicluster_srrtot_pTsorted.push_back(candidate_maxicluster_srrtot.at(index));
                candidate_maxicluster_srrmean_pTsorted.push_back(candidate_maxicluster_srrmean.at(index));

            }

            // REMOVE DUPLICATED SUPERCLUSTERS

            float new_pt_leading_maxi = candidate_maxicluster_pTsorted.at(0).Pt();

            if(new_pt_leading_maxi == previous_pt_leading_maxi) 
               continue;

            previous_pt_leading_maxi = new_pt_leading_maxi;

            _maxicluster_cl3d_pt.push_back(candidate_maxicluster_pt_pTsorted);
            _maxicluster_cl3d_energy.push_back(candidate_maxicluster_E_pTsorted);
            _maxicluster_cl3d_eta.push_back(candidate_maxicluster_eta_pTsorted);
            _maxicluster_cl3d_phi.push_back(candidate_maxicluster_phi_pTsorted);
            _maxicluster_cl3d_bdteg.push_back(candidate_maxicluster_bdteg_pTsorted);
            _maxicluster_cl3d_showerlength.push_back(candidate_maxicluster_showerlength_pTsorted);
            _maxicluster_cl3d_coreshowerlength.push_back(candidate_maxicluster_coreshowerlength_pTsorted);
            _maxicluster_cl3d_firstlayer.push_back(candidate_maxicluster_firstlayer_pTsorted);
            _maxicluster_cl3d_maxlayer.push_back(candidate_maxicluster_maxlayer_pTsorted);
            _maxicluster_cl3d_seetot.push_back(candidate_maxicluster_seetot_pTsorted);
            _maxicluster_cl3d_spptot.push_back(candidate_maxicluster_spptot_pTsorted);
            _maxicluster_cl3d_szz.push_back(candidate_maxicluster_szz_pTsorted);
            _maxicluster_cl3d_srrtot.push_back(candidate_maxicluster_srrtot_pTsorted);
            _maxicluster_cl3d_srrmean.push_back(candidate_maxicluster_srrmean_pTsorted);

        }

        _maxicluster_n = _maxicluster_cl3d_pt.size();


        // MINICLUSTERS

        for (int i_main=0; i_main<_cl3d_n; i_main++){

            if ( (*_cl3d_pt)[i_main] < 0.0001 || (*_cl3d_pt)[i_main] > 10000 ) continue;

            // RAW CLUSTERING 

            vector<TLorentzVector> candidate_minicluster;
            vector<float> candidate_minicluster_bdteg;
            vector<int>   candidate_minicluster_showerlength;
            vector<int>   candidate_minicluster_coreshowerlength;
            vector<int>   candidate_minicluster_firstlayer;
            vector<int>   candidate_minicluster_maxlayer;
            vector<float> candidate_minicluster_seetot;
            vector<float> candidate_minicluster_spptot;
            vector<float> candidate_minicluster_szz;
            vector<float> candidate_minicluster_srrtot;
            vector<float> candidate_minicluster_srrmean;

            candidate_minicluster.clear();
            candidate_minicluster_bdteg.clear();
            candidate_minicluster_showerlength.clear();
            candidate_minicluster_coreshowerlength.clear();
            candidate_minicluster_firstlayer.clear();
            candidate_minicluster_maxlayer.clear();
            candidate_minicluster_seetot.clear();
            candidate_minicluster_spptot.clear();
            candidate_minicluster_szz.clear();
            candidate_minicluster_srrtot.clear();
            candidate_minicluster_srrmean.clear();


            // MAIN CLUSTER CANDIDATE (SUPERCLUSTER SEED)

            TLorentzVector seed_candidate_minicluster;
            float seed_candidate_minicluster_bdteg;
            int   seed_candidate_minicluster_showerlength;
            int   seed_candidate_minicluster_coreshowerlength;
            int   seed_candidate_minicluster_firstlayer;
            int   seed_candidate_minicluster_maxlayer;
            float seed_candidate_minicluster_seetot;
            float seed_candidate_minicluster_spptot;
            float seed_candidate_minicluster_szz;
            float seed_candidate_minicluster_srrtot;
            float seed_candidate_minicluster_srrmean;

            seed_candidate_minicluster.SetPtEtaPhiM( (*_cl3d_pt)[i_main], (*_cl3d_eta)[i_main], (*_cl3d_phi)[i_main], 0);
            seed_candidate_minicluster_bdteg = (*_cl3d_bdteg)[i_main];
            seed_candidate_minicluster_showerlength = (*_cl3d_showerlength)[i_main];
            seed_candidate_minicluster_coreshowerlength = (*_cl3d_coreshowerlength)[i_main];
            seed_candidate_minicluster_firstlayer = (*_cl3d_firstlayer)[i_main];
            seed_candidate_minicluster_maxlayer = (*_cl3d_maxlayer)[i_main];
            seed_candidate_minicluster_seetot = (*_cl3d_seetot)[i_main];
            seed_candidate_minicluster_spptot = (*_cl3d_spptot)[i_main];
            seed_candidate_minicluster_szz = (*_cl3d_szz)[i_main];
            seed_candidate_minicluster_srrtot = (*_cl3d_srrtot)[i_main];
            seed_candidate_minicluster_srrmean = (*_cl3d_srrmean)[i_main];


            if ( seed_candidate_minicluster.Pt() < pt_thr_seed ) continue;

            candidate_minicluster.push_back(seed_candidate_minicluster);
            candidate_minicluster_bdteg.push_back(seed_candidate_minicluster_bdteg);
            candidate_minicluster_showerlength.push_back(seed_candidate_minicluster_showerlength);
            candidate_minicluster_coreshowerlength.push_back(seed_candidate_minicluster_coreshowerlength);
            candidate_minicluster_firstlayer.push_back(seed_candidate_minicluster_firstlayer);
            candidate_minicluster_maxlayer.push_back(seed_candidate_minicluster_maxlayer);
            candidate_minicluster_seetot.push_back(seed_candidate_minicluster_seetot);
            candidate_minicluster_spptot.push_back(seed_candidate_minicluster_spptot);
            candidate_minicluster_szz.push_back(seed_candidate_minicluster_szz);
            candidate_minicluster_srrtot.push_back(seed_candidate_minicluster_srrtot);
            candidate_minicluster_srrmean.push_back(seed_candidate_minicluster_srrmean);

            // SECONDARY CLUSTERS CANDIDATES

            for (int i_sec=0; i_sec<_cl3d_n; i_sec++){

                if( i_sec==i_main ) continue;

                if ( (*_cl3d_pt)[i_sec] < 0.0001 || (*_cl3d_pt)[i_sec] > 10000 ) continue;

                TLorentzVector secondary_candidate_minicluster;
                float secondary_candidate_minicluster_bdteg;
                int   secondary_candidate_minicluster_showerlength;
                int   secondary_candidate_minicluster_coreshowerlength;
                int   secondary_candidate_minicluster_firstlayer;
                int   secondary_candidate_minicluster_maxlayer;
                float secondary_candidate_minicluster_seetot;
                float secondary_candidate_minicluster_spptot;
                float secondary_candidate_minicluster_szz;
                float secondary_candidate_minicluster_srrtot;
                float secondary_candidate_minicluster_srrmean;

                secondary_candidate_minicluster.SetPtEtaPhiM( (*_cl3d_pt)[i_sec], (*_cl3d_eta)[i_sec], (*_cl3d_phi)[i_sec], 0);
                secondary_candidate_minicluster_bdteg = (*_cl3d_bdteg)[i_sec];
                secondary_candidate_minicluster_showerlength = (*_cl3d_showerlength)[i_sec];
                secondary_candidate_minicluster_coreshowerlength = (*_cl3d_coreshowerlength)[i_sec];
                secondary_candidate_minicluster_firstlayer = (*_cl3d_firstlayer)[i_sec];
                secondary_candidate_minicluster_maxlayer = (*_cl3d_maxlayer)[i_sec];
                secondary_candidate_minicluster_seetot = (*_cl3d_seetot)[i_sec];
                secondary_candidate_minicluster_spptot = (*_cl3d_spptot)[i_sec];
                secondary_candidate_minicluster_szz = (*_cl3d_szz)[i_sec];
                secondary_candidate_minicluster_srrtot = (*_cl3d_srrtot)[i_sec];
                secondary_candidate_minicluster_srrmean = (*_cl3d_srrmean)[i_sec];

                if ( secondary_candidate_minicluster.Pt() < pt_thr_sec ) continue;

                if( abs(secondary_candidate_minicluster.Eta() - seed_candidate_minicluster.Eta()) > eta_window_mini ) continue;
                if( abs(secondary_candidate_minicluster.Phi() - seed_candidate_minicluster.Phi()) > phi_window_mini ) continue;

                candidate_minicluster.push_back(secondary_candidate_minicluster);
                candidate_minicluster_bdteg.push_back(secondary_candidate_minicluster_bdteg);
                candidate_minicluster_showerlength.push_back(secondary_candidate_minicluster_showerlength);
                candidate_minicluster_coreshowerlength.push_back(secondary_candidate_minicluster_coreshowerlength);
                candidate_minicluster_firstlayer.push_back(secondary_candidate_minicluster_firstlayer);
                candidate_minicluster_maxlayer.push_back(secondary_candidate_minicluster_maxlayer);
                candidate_minicluster_seetot.push_back(secondary_candidate_minicluster_seetot);
                candidate_minicluster_spptot.push_back(secondary_candidate_minicluster_spptot);
                candidate_minicluster_szz.push_back(secondary_candidate_minicluster_szz);
                candidate_minicluster_srrtot.push_back(secondary_candidate_minicluster_srrtot);
                candidate_minicluster_srrmean.push_back(secondary_candidate_minicluster_srrmean);

            }

            // SORT CLUSTERS BY PT IN SUPERCLUSTER

            vector<TLorentzVector> candidate_minicluster_pTsorted;
            vector<float> candidate_minicluster_pt_pTsorted;
            vector<float> candidate_minicluster_E_pTsorted;
            vector<float> candidate_minicluster_eta_pTsorted;
            vector<float> candidate_minicluster_phi_pTsorted;
            vector<float> candidate_minicluster_bdteg_pTsorted;
            vector<int>   candidate_minicluster_showerlength_pTsorted;
            vector<int>   candidate_minicluster_coreshowerlength_pTsorted;
            vector<int>   candidate_minicluster_firstlayer_pTsorted;
            vector<int>   candidate_minicluster_maxlayer_pTsorted;
            vector<float> candidate_minicluster_seetot_pTsorted;
            vector<float> candidate_minicluster_spptot_pTsorted;
            vector<float> candidate_minicluster_szz_pTsorted;
            vector<float> candidate_minicluster_srrtot_pTsorted;
            vector<float> candidate_minicluster_srrmean_pTsorted;

            candidate_minicluster_pTsorted.clear();
            candidate_minicluster_pt_pTsorted.clear();
            candidate_minicluster_E_pTsorted.clear();
            candidate_minicluster_eta_pTsorted.clear();
            candidate_minicluster_phi_pTsorted.clear();
            candidate_minicluster_bdteg_pTsorted.clear();
            candidate_minicluster_showerlength_pTsorted.clear();
            candidate_minicluster_coreshowerlength_pTsorted.clear();
            candidate_minicluster_firstlayer_pTsorted.clear();
            candidate_minicluster_maxlayer_pTsorted.clear();
            candidate_minicluster_seetot_pTsorted.clear();
            candidate_minicluster_spptot_pTsorted.clear();
            candidate_minicluster_szz_pTsorted.clear();
            candidate_minicluster_srrtot_pTsorted.clear();
            candidate_minicluster_srrmean_pTsorted.clear();

            vector< pair<int,TLorentzVector> > iminicluster_minicluster_pairs;

            for (unsigned int i_cluster = 0; i_cluster<candidate_minicluster.size(); i_cluster++){

                pair<int,TLorentzVector> minicluster_pair = make_pair(i_cluster,candidate_minicluster.at(i_cluster));
                iminicluster_minicluster_pairs.push_back(minicluster_pair);

            }

            sort(iminicluster_minicluster_pairs.begin(), iminicluster_minicluster_pairs.end(), pT_comparison_pairs);

            for (unsigned int i_cluster = 0; i_cluster<candidate_minicluster.size(); i_cluster++){

                int index = iminicluster_minicluster_pairs[i_cluster].first;
                TLorentzVector miniclus = iminicluster_minicluster_pairs[i_cluster].second;

                candidate_minicluster_pTsorted.push_back(miniclus);
                candidate_minicluster_pt_pTsorted.push_back(miniclus.Pt());
                candidate_minicluster_E_pTsorted.push_back(miniclus.E());
                candidate_minicluster_eta_pTsorted.push_back(miniclus.Eta());
                candidate_minicluster_phi_pTsorted.push_back(miniclus.Phi());
                candidate_minicluster_bdteg_pTsorted.push_back(candidate_minicluster_bdteg.at(index));
                candidate_minicluster_showerlength_pTsorted.push_back(candidate_minicluster_showerlength.at(index));
                candidate_minicluster_coreshowerlength_pTsorted.push_back(candidate_minicluster_coreshowerlength.at(index));
                candidate_minicluster_firstlayer_pTsorted.push_back(candidate_minicluster_firstlayer.at(index));
                candidate_minicluster_maxlayer_pTsorted.push_back(candidate_minicluster_maxlayer.at(index));
                candidate_minicluster_seetot_pTsorted.push_back(candidate_minicluster_seetot.at(index));
                candidate_minicluster_spptot_pTsorted.push_back(candidate_minicluster_spptot.at(index));
                candidate_minicluster_szz_pTsorted.push_back(candidate_minicluster_szz.at(index));
                candidate_minicluster_srrtot_pTsorted.push_back(candidate_minicluster_srrtot.at(index));
                candidate_minicluster_srrmean_pTsorted.push_back(candidate_minicluster_srrmean.at(index));

            }

            // REMOVE DUPLICATED SUPERCLUSTERS

            float new_pt_leading_mini = candidate_minicluster_pTsorted.at(0).Pt();

            if(new_pt_leading_mini == previous_pt_leading_mini) 
               continue;

            previous_pt_leading_mini = new_pt_leading_mini;

            _minicluster_cl3d_pt.push_back(candidate_minicluster_pt_pTsorted);
            _minicluster_cl3d_energy.push_back(candidate_minicluster_E_pTsorted);
            _minicluster_cl3d_eta.push_back(candidate_minicluster_eta_pTsorted);
            _minicluster_cl3d_phi.push_back(candidate_minicluster_phi_pTsorted);
            _minicluster_cl3d_bdteg.push_back(candidate_minicluster_bdteg_pTsorted);
            _minicluster_cl3d_showerlength.push_back(candidate_minicluster_showerlength_pTsorted);
            _minicluster_cl3d_coreshowerlength.push_back(candidate_minicluster_coreshowerlength_pTsorted);
            _minicluster_cl3d_firstlayer.push_back(candidate_minicluster_firstlayer_pTsorted);
            _minicluster_cl3d_maxlayer.push_back(candidate_minicluster_maxlayer_pTsorted);
            _minicluster_cl3d_seetot.push_back(candidate_minicluster_seetot_pTsorted);
            _minicluster_cl3d_spptot.push_back(candidate_minicluster_spptot_pTsorted);
            _minicluster_cl3d_szz.push_back(candidate_minicluster_szz_pTsorted);
            _minicluster_cl3d_srrtot.push_back(candidate_minicluster_srrtot_pTsorted);
            _minicluster_cl3d_srrmean.push_back(candidate_minicluster_srrmean_pTsorted);

        }

        _minicluster_n = _minicluster_cl3d_pt.size();


        // MERGING AND OVERLAP REMOVAL

        vector<float> maxiclust_cl3d_pt;
        vector<float> maxiclust_cl3d_energy;
        vector<float> maxiclust_cl3d_eta;
        vector<float> maxiclust_cl3d_phi;
        /*vector<float> maxiclust_cl3d_bdteg;
        vector<int>   maxiclust_cl3d_showerlength;
        vector<int>   maxiclust_cl3d_coreshowerlength;
        vector<int>   maxiclust_cl3d_firstlayer;
        vector<int>   maxiclust_cl3d_maxlayer;
        vector<float> maxiclust_cl3d_seetot;
        vector<float> maxiclust_cl3d_spptot;
        vector<float> maxiclust_cl3d_szz;
        vector<float> maxiclust_cl3d_srrtot;
        vector<float> maxiclust_cl3d_srrmean;*/

        maxiclust_cl3d_pt.clear();
        maxiclust_cl3d_energy.clear();
        maxiclust_cl3d_eta.clear();
        maxiclust_cl3d_phi.clear();
        /*maxiclust_cl3d_bdteg.clear();
        maxiclust_cl3d_showerlength.clear();
        maxiclust_cl3d_coreshowerlength.clear();
        maxiclust_cl3d_firstlayer.clear();
        maxiclust_cl3d_maxlayer.clear();
        maxiclust_cl3d_seetot.clear();
        maxiclust_cl3d_spptot.clear();
        maxiclust_cl3d_szz.clear();
        maxiclust_cl3d_srrtot.clear();
        maxiclust_cl3d_srrmean.clear();*/

        //cout<<"Number of maxiclusters: "<<_maxicluster_n<<endl;

        for (int i_maxicluster=0; i_maxicluster<_maxicluster_n; i_maxicluster++){

            //cout<<"With subclusters: "<<_maxicluster_cl3d_pt.at(i_maxicluster).size()<<endl;

            maxiclust_cl3d_pt  = _maxicluster_cl3d_pt.at(i_maxicluster);
            maxiclust_cl3d_energy = _maxicluster_cl3d_energy.at(i_maxicluster);
            maxiclust_cl3d_eta = _maxicluster_cl3d_eta.at(i_maxicluster);
            maxiclust_cl3d_phi = _maxicluster_cl3d_phi.at(i_maxicluster);

            cout<<"Maxicluster #"<<i_maxicluster<<",with leading cl3d (pT,eta,phi): ("<<maxiclust_cl3d_pt.at(0)<<","<<maxiclust_cl3d_eta.at(0)<<","<<maxiclust_cl3d_phi.at(0)<<")"<<endl;

            /*maxiclust_cl3d_bdteg = _maxicluster_cl3d_bdteg.at(i_maxicluster);
            maxiclust_cl3d_showerlength = _maxicluster_cl3d_showerlength.at(i_maxicluster);
            maxiclust_cl3d_coreshowerlength = _maxicluster_cl3d_coreshowerlength.at(i_maxicluster);
            maxiclust_cl3d_firstlayer = _maxicluster_cl3d_firstlayer.at(i_maxicluster);
            maxiclust_cl3d_maxlayer = _maxicluster_cl3d_maxlayer.at(i_maxicluster);
            maxiclust_cl3d_seetot = _maxicluster_cl3d_seetot.at(i_maxicluster);
            maxiclust_cl3d_spptot = _maxicluster_cl3d_spptot.at(i_maxicluster);
            maxiclust_cl3d_szz = _maxicluster_cl3d_szz.at(i_maxicluster);
            maxiclust_cl3d_srrtot = _maxicluster_cl3d_srrtot.at(i_maxicluster);
            maxiclust_cl3d_srrmean = _maxicluster_cl3d_srrmean.at(i_maxicluster);*/

            /*for (int i_cl3d=0; i_cl3d<_maxicluster_cl3d_pt.size(); i_cl3d++)

                _supercluster_cl3d_pt.push_back(maxiclust_cl3d_pt.at(i_cl3d));
                _supercluster_cl3d_energy.push_back(maxiclust_cl3d_energy.at(i_cl3d));
                _supercluster_cl3d_eta.push_back(maxiclust_cl3d_eta.at(i_cl3d));
                _supercluster_cl3d_phi.push_back(maxiclust_cl3d_phi.at(i_cl3d));
                _supercluster_cl3d_bdteg.push_back(maxiclust_cl3d_bdteg.at(i_cl3d));
                _supercluster_cl3d_showerlength.push_back(maxiclust_cl3d_showerlength.at(i_cl3d));
                _supercluster_cl3d_coreshowerlength.push_back(maxiclust_cl3d_coreshowerlength.at(i_cl3d));
                _supercluster_cl3d_firstlayer.push_back(maxiclust_cl3d_firstlayer.at(i_cl3d));
                _supercluster_cl3d_maxlayer.push_back(maxiclust_cl3d_maxlayer.at(i_cl3d));
                _supercluster_cl3d_seetot.push_back(maxiclust_cl3d_seetot.at(i_cl3d));
                _supercluster_cl3d_spptot.push_back(maxiclust_cl3d_spptot.at(i_cl3d));
                _supercluster_cl3d_szz.push_back(maxiclust_cl3d_szz.at(i_cl3d));
                _supercluster_cl3d_srrtot.push_back(maxiclust_cl3d_srrtot.at(i_cl3d));
                _supercluster_cl3d_srrmean.push_back(maxiclust_cl3d_srrmean.at(i_cl3d));

            }*/

            int isMerged = false;

            vector<float> miniclust_cl3d_pt;
            vector<float> miniclust_cl3d_energy;
            vector<float> miniclust_cl3d_eta;
            vector<float> miniclust_cl3d_phi;
            /*vector<float> miniclust_cl3d_bdteg;
            vector<int>   miniclust_cl3d_showerlength;
            vector<int>   miniclust_cl3d_coreshowerlength;
            vector<int>   miniclust_cl3d_firstlayer;
            vector<int>   miniclust_cl3d_maxlayer;
            vector<float> miniclust_cl3d_seetot;
            vector<float> miniclust_cl3d_spptot;
            vector<float> miniclust_cl3d_szz;
            vector<float> miniclust_cl3d_srrtot;
            vector<float> miniclust_cl3d_srrmean;*/

            miniclust_cl3d_pt.clear();
            miniclust_cl3d_energy.clear();
            miniclust_cl3d_eta.clear();
            miniclust_cl3d_phi.clear();
            /*miniclust_cl3d_bdteg.clear();
            miniclust_cl3d_showerlength.clear();
            miniclust_cl3d_coreshowerlength.clear();
            miniclust_cl3d_firstlayer.clear();
            miniclust_cl3d_maxlayer.clear();
            miniclust_cl3d_seetot.clear();
            miniclust_cl3d_spptot.clear();
            miniclust_cl3d_szz.clear();
            miniclust_cl3d_srrtot.clear();
            miniclust_cl3d_srrmean.clear();*/

            for (int i_minicluster=0; i_minicluster<_minicluster_n; i_minicluster++){

                miniclust_cl3d_pt  = _minicluster_cl3d_pt.at(i_minicluster);
                miniclust_cl3d_energy = _minicluster_cl3d_energy.at(i_minicluster);
                miniclust_cl3d_eta = _minicluster_cl3d_eta.at(i_minicluster);
                miniclust_cl3d_phi = _minicluster_cl3d_phi.at(i_minicluster);

                /*miniclust_cl3d_bdteg = _minicluster_cl3d_bdteg.at(i_minicluster);
                miniclust_cl3d_showerlength = _minicluster_cl3d_showerlength.at(i_minicluster);
                miniclust_cl3d_coreshowerlength = _minicluster_cl3d_coreshowerlength.at(i_minicluster);
                miniclust_cl3d_firstlayer = _minicluster_cl3d_firstlayer.at(i_minicluster);
                miniclust_cl3d_maxlayer = _minicluster_cl3d_maxlayer.at(i_minicluster);
                miniclust_cl3d_seetot = _minicluster_cl3d_seetot.at(i_minicluster);
                miniclust_cl3d_spptot = _minicluster_cl3d_spptot.at(i_minicluster);
                miniclust_cl3d_szz = _minicluster_cl3d_szz.at(i_minicluster);
                miniclust_cl3d_srrtot = _minicluster_cl3d_srrtot.at(i_minicluster);
                miniclust_cl3d_srrmean = _minicluster_cl3d_srrmean.at(i_minicluster);*/

                if(maxiclust_cl3d_pt.at(0) == miniclust_cl3d_pt.at(0)) continue;

                if(fabs(maxiclust_cl3d_eta.at(0)-miniclust_cl3d_eta.at(0))>=0.15) continue;
                if(fabs(maxiclust_cl3d_phi.at(0)-miniclust_cl3d_phi.at(0))<=0.45 || fabs(maxiclust_cl3d_phi.at(0)-miniclust_cl3d_phi.at(0))>=0.60) continue;


                cout<<"Minicluster #"<<i_minicluster<<",with leading cl3d (pT,eta,phi): ("<<miniclust_cl3d_pt.at(0)<<","<<miniclust_cl3d_eta.at(0)<<","<<miniclust_cl3d_phi.at(0)<<")"<<endl;

                /*_supercluster_cl3d_pt.push_back(miniclust_cl3d_pt);
                _supercluster_cl3d_energy.push_back(miniclust_cl3d_energy);
                _supercluster_cl3d_eta.push_back(miniclust_cl3d_eta);
                _supercluster_cl3d_phi.push_back(miniclust_cl3d_phi);
                _supercluster_cl3d_bdteg.push_back(miniclust_cl3d_bdteg);
                _supercluster_cl3d_showerlength.push_back(miniclust_cl3d_showerlength);
                _supercluster_cl3d_coreshowerlength.push_back(miniclust_cl3d_coreshowerlength);
                _supercluster_cl3d_firstlayer.push_back(miniclust_cl3d_firstlayer);
                _supercluster_cl3d_maxlayer.push_back(miniclust_cl3d_maxlayer);
                _supercluster_cl3d_seetot.push_back(miniclust_cl3d_seetot);
                _supercluster_cl3d_spptot.push_back(miniclust_cl3d_spptot);
                _supercluster_cl3d_szz.push_back(miniclust_cl3d_szz);
                _supercluster_cl3d_srrtot.push_back(miniclust_cl3d_srrtot);
                _supercluster_cl3d_srrmean.push_back(miniclust_cl3d_srrmean);*/

                isMerged = true;

            }

            _supercluster_isMerged.push_back(isMerged);

        }

        _supercluster_n = _supercluster_cl3d_pt.size();


        out_tree->Fill();

    }

    out_file->cd();

    out_tree->Write();
    out_file->Close();

    return;

}

void test(int n_events = -1, TString pu = "0", float ptthr_seed = 4, float ptthr_sec = 2, float etawindowmax = 0.15, float phiwindowmax = 0.45, float etawindowmin = 0.15, float phiwindowmin = 0.15){

  TString dir = "/data_CMS/cms/mperez/HGCal_data/May19/";

  //cout<<"pu"<<pu<<endl;

  int i_ptthr_seed = ptthr_seed;
  int i_ptthr_sec = ptthr_sec;
  int i_etawindowmax = 100*etawindowmax;
  int i_phiwindowmax = 100*phiwindowmax;
  int i_etawindowmin = 100*etawindowmin;
  int i_phiwindowmin = 100*phiwindowmin;

  TString s_ptthr_seed = to_string(i_ptthr_seed);
  TString s_ptthr_sec  = to_string(i_ptthr_sec);
  TString s_etawindowmax  = to_string(i_etawindowmax);
  TString s_phiwindowmax  = to_string(i_phiwindowmax);
  TString s_etawindowmin  = to_string(i_etawindowmin);
  TString s_phiwindowmin  = to_string(i_phiwindowmin);

  TString infile = dir+"skimmed/NTuple_ZTT_PU"+pu+"_skimmed.root";
  TString outfile = dir+"clustered/NTuple_ZTT_PU"+pu+"_clustered_pT"+s_ptthr_seed+"-"+s_ptthr_sec+"_maxi"+s_etawindowmax+"-"+s_phiwindowmax+"_mini"+s_etawindowmin+"-"+s_phiwindowmin+".root";

  cluster_tree(infile, outfile, n_events, ptthr_seed, ptthr_sec, etawindowmax, phiwindowmax, etawindowmin, phiwindowmin);

}
