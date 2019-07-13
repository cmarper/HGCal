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

void cluster_tree( TString filein, TString fileout, int nevents = -1, float pt_thr_seed = 4, float pt_thr_sec = 2, float eta_window = 0.3, float phi_window = 0.5){

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

    TH1F* h_pT_cl1=new TH1F("h_pT_cl1","h_pT_cl1",100,0,50);
    TH1F* h_pT_cl2=new TH1F("h_pT_cl2","h_pT_cl2",100,0,50);
    TH1F* h_pT_cl3=new TH1F("h_pT_cl3","h_pT_cl3",100,0,50);
    TH1F* h_pT_cl4=new TH1F("h_pT_cl4","h_pT_cl4",100,0,50);

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
    vector<int>    *_cl3d_coreshowerlength;
    vector<int>    *_cl3d_maxlayer;

    in_tree->SetBranchAddress("gentau_n",&_gentau_n);

    in_tree->SetBranchAddress("gentau_pt",        &_gentau_pt);
    in_tree->SetBranchAddress("gentau_eta",        &_gentau_eta);
    in_tree->SetBranchAddress("gentau_phi",        &_gentau_phi);
    in_tree->SetBranchAddress("gentau_energy",    &_gentau_energy);
    in_tree->SetBranchAddress("gentau_mass",    &_gentau_mass);

    in_tree->SetBranchAddress("gentau_vis_pt",        &_gentau_vis_pt);
    in_tree->SetBranchAddress("gentau_vis_eta",        &_gentau_vis_eta);
    in_tree->SetBranchAddress("gentau_vis_phi",        &_gentau_vis_phi);
    in_tree->SetBranchAddress("gentau_vis_energy",    &_gentau_vis_energy);
    in_tree->SetBranchAddress("gentau_vis_mass",    &_gentau_vis_mass);

    in_tree->SetBranchAddress("gentau_decayMode",    &_gentau_decayMode);

    in_tree->SetBranchAddress("cl3d_n",    &_cl3d_n);

    in_tree->SetBranchAddress("cl3d_pt",                 &_cl3d_pt);
    in_tree->SetBranchAddress("cl3d_energy",             &_cl3d_energy);
    in_tree->SetBranchAddress("cl3d_eta",                 &_cl3d_eta);
    in_tree->SetBranchAddress("cl3d_phi",                 &_cl3d_phi);
    in_tree->SetBranchAddress("cl3d_bdteg",                 &_cl3d_bdteg);
    in_tree->SetBranchAddress("cl3d_coreshowerlength",     &_cl3d_coreshowerlength);
    in_tree->SetBranchAddress("cl3d_maxlayer",             &_cl3d_maxlayer);

    TTree* out_tree=in_tree->GetTree()->CloneTree(0);
    out_tree->SetNameTitle("ClusteredTree","ClusteredTree");

    float _n_supercl3ds;

    vector<int>    _supercl3d_n_cl3d; //number of sub cl3d
    vector<float>  _supercl3d_pt; //sum of cl3d pT
    vector<float>  _supercl3d_energy; //sum of cl3d E
    vector<float>  _supercl3d_eta; //eta of main cl3d
    vector<float>  _supercl3d_phi; //phi of main cl3d

    vector<vector<float>>  _supercl3d_cl3d_pt;
    vector<vector<float>>  _supercl3d_cl3d_energy;
    vector<vector<float>>  _supercl3d_cl3d_eta;
    vector<vector<float>>  _supercl3d_cl3d_phi;
    vector<vector<float>>  _supercl3d_cl3d_bdteg;
    vector<vector<int>>    _supercl3d_cl3d_coreshowerlength;
    vector<vector<int>>    _supercl3d_cl3d_maxlayer;

    vector<bool>               _gentau_isMatched;
    vector<int>               _gentau_numberMatchedcl3d;
    vector<vector<float> >    _gentau_PtMatchedcl3d;
    vector<vector<float> >    _gentau_EtaMatchedcl3d;
    vector<vector<float> >    _gentau_PhiMatchedcl3d;
    vector<vector<float> >    _gentau_BDTegMatchedcl3d;
    vector<vector<int> >    _gentau_CoreShowerLengthMatchedcl3d;
    vector<vector<int> >    _gentau_MaxLayerMatchedcl3d;

    out_tree->Branch("n_supercl3d",                 &_n_supercl3ds);

    out_tree->Branch("supercl3d_n_cl3d",            &_supercl3d_n_cl3d);
    out_tree->Branch("supercl3d_pt",                &_supercl3d_pt);
    out_tree->Branch("supercl3d_energy",            &_supercl3d_energy);
    out_tree->Branch("supercl3d_eta",                &_supercl3d_eta);
    out_tree->Branch("supercl3d_phi",                &_supercl3d_phi);

    out_tree->Branch("supercl3d_cl3d_pt",                    &_supercl3d_cl3d_pt);
    out_tree->Branch("supercl3d_cl3d_energy",                &_supercl3d_cl3d_energy);
    out_tree->Branch("supercl3d_cl3d_eta",                    &_supercl3d_cl3d_eta);
    out_tree->Branch("supercl3d_cl3d_phi",                    &_supercl3d_cl3d_phi);
    out_tree->Branch("supercl3d_cl3d_bdteg",                &_supercl3d_cl3d_bdteg);
    out_tree->Branch("supercl3d_cl3d_coreshowerlength",        &_supercl3d_cl3d_coreshowerlength);
    out_tree->Branch("supercl3d_cl3d_maxlayer",                &_supercl3d_cl3d_maxlayer);

    out_tree->Branch("gentau_isMatched", &_gentau_isMatched);
    out_tree->Branch("gentau_numberMatchedcl3d", &_gentau_numberMatchedcl3d);
    out_tree->Branch("gentau_PtMatchedcl3d", &_gentau_PtMatchedcl3d);
    out_tree->Branch("gentau_EtaMatchedcl3d", &_gentau_EtaMatchedcl3d);
    out_tree->Branch("gentau_PhiMatchedcl3d", &_gentau_PhiMatchedcl3d);
    out_tree->Branch("gentau_BDTegMatchedcl3d", &_gentau_BDTegMatchedcl3d);
    out_tree->Branch("gentau_CoreShowerLengthMatchedcl3d", &_gentau_CoreShowerLengthMatchedcl3d);
    out_tree->Branch("gentau_MaxLayerMatchedcl3d", &_gentau_MaxLayerMatchedcl3d);


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
        _cl3d_coreshowerlength = 0;
        _cl3d_maxlayer = 0;

        _n_supercl3ds = 0;

        _supercl3d_n_cl3d.clear();
        _supercl3d_pt.clear();
        _supercl3d_energy.clear();
        _supercl3d_eta.clear();
        _supercl3d_phi.clear();

        _supercl3d_cl3d_pt.clear();
        _supercl3d_cl3d_energy.clear();
        _supercl3d_cl3d_eta.clear();
        _supercl3d_cl3d_phi.clear();
        _supercl3d_cl3d_bdteg.clear();
        _supercl3d_cl3d_coreshowerlength.clear();
        _supercl3d_cl3d_maxlayer.clear();

        _gentau_isMatched.clear();
        _gentau_numberMatchedcl3d.clear();
        _gentau_PtMatchedcl3d.clear();
        _gentau_EtaMatchedcl3d.clear();
        _gentau_PhiMatchedcl3d.clear();
        _gentau_BDTegMatchedcl3d.clear();
        _gentau_CoreShowerLengthMatchedcl3d.clear();
        _gentau_MaxLayerMatchedcl3d.clear();


        int entry_ok = in_tree->GetEntry(i);    
        if(entry_ok<0) 
            continue;

        float previous_pt_leading = -999.;


        //////////////////////////////////////////
        ////////////// CLUSTERING ////////////////
        //////////////////////////////////////////

        for (int i_main=0; i_main<_cl3d_n; i_main++){

            if ( (*_cl3d_pt)[i_main] < 0.0001 || (*_cl3d_pt)[i_main] > 10000 ) continue;

            // RAW CLUSTERING 

            vector<TLorentzVector> clusters;
            vector<float> clusters_bdteg;
            vector<int> clusters_coreshowerlength;
            vector<int> clusters_maxlayer;

            clusters.clear();
            clusters_bdteg.clear();
            clusters_coreshowerlength.clear();
            clusters_maxlayer.clear();

            // MAIN CLUSTER CANDIDATE (SUPERCLUSTER SEED)

            TLorentzVector main;
            float main_bdteg;
            int main_coreshowerlength;
            int main_maxlayer;

            main.SetPtEtaPhiM( (*_cl3d_pt)[i_main], (*_cl3d_eta)[i_main], (*_cl3d_phi)[i_main], 0);
            main_bdteg = (*_cl3d_bdteg)[i_main];
            main_coreshowerlength = (*_cl3d_coreshowerlength)[i_main];
            main_maxlayer = (*_cl3d_maxlayer)[i_main];


            if ( main.Pt() < pt_thr_seed ) continue;

            clusters.push_back(main);
            clusters_bdteg.push_back(main_bdteg);
            clusters_coreshowerlength.push_back(main_coreshowerlength);
            clusters_maxlayer.push_back(main_maxlayer);

            // SECONDARY CLUSTERS CANDIDATES

            for (int i_sec=0; i_sec<_cl3d_n; i_sec++){

                if( i_sec==i_main ) continue;

                if ( (*_cl3d_pt)[i_sec] < 0.0001 || (*_cl3d_pt)[i_sec] > 10000 ) continue;

                TLorentzVector secondary;
                float secondary_bdteg;
                int secondary_coreshowerlength;
                int secondary_maxlayer;

                secondary.SetPtEtaPhiM( (*_cl3d_pt)[i_sec], (*_cl3d_eta)[i_sec], (*_cl3d_phi)[i_sec], 0);
                secondary_bdteg = (*_cl3d_bdteg)[i_sec];
                secondary_coreshowerlength = (*_cl3d_coreshowerlength)[i_sec];
                secondary_maxlayer = (*_cl3d_maxlayer)[i_sec];

                if ( secondary.Pt() < pt_thr_sec ) continue;

                if( abs(secondary.Eta() - main.Eta()) > eta_window ) continue;
                if( abs(secondary.Phi() - main.Phi()) > phi_window ) continue;

                clusters.push_back(secondary);
                clusters_bdteg.push_back(secondary_bdteg);
                clusters_coreshowerlength.push_back(secondary_coreshowerlength);
                clusters_maxlayer.push_back(secondary_maxlayer);

            }

            // SORT CLUSTERS BY PT IN SUPERCLUSTER

            vector<TLorentzVector> clusters_pTsorted;
            vector<float> clusters_pt_pTsorted;
            vector<float> clusters_E_pTsorted;
            vector<float> clusters_eta_pTsorted;
            vector<float> clusters_phi_pTsorted;
            vector<float> clusters_bdteg_pTsorted;
            vector<int> clusters_coreshowerlength_pTsorted;
            vector<int> clusters_maxlayer_pTsorted;

            clusters_pTsorted.clear();
            clusters_bdteg_pTsorted.clear();
            clusters_coreshowerlength_pTsorted.clear();
            clusters_maxlayer_pTsorted.clear();

            vector< pair<int,TLorentzVector> > icluster_cluster_pairs;

            for (unsigned int i_cluster = 0; i_cluster<clusters.size(); i_cluster++){

                pair<int,TLorentzVector> cluster_pair = make_pair(i_cluster,clusters.at(i_cluster));
                icluster_cluster_pairs.push_back(cluster_pair);

            }

            sort(icluster_cluster_pairs.begin(), icluster_cluster_pairs.end(), pT_comparison_pairs);

            for (unsigned int i_cluster = 0; i_cluster<clusters.size(); i_cluster++){

                int index = icluster_cluster_pairs[i_cluster].first;
                TLorentzVector tlv = icluster_cluster_pairs[i_cluster].second;

                clusters_pTsorted.push_back(tlv);
                clusters_pt_pTsorted.push_back(tlv.Pt());
                clusters_E_pTsorted.push_back(tlv.E());
                clusters_eta_pTsorted.push_back(tlv.Eta());
                clusters_phi_pTsorted.push_back(tlv.Phi());
                clusters_bdteg_pTsorted.push_back(clusters_bdteg.at(index));
                clusters_coreshowerlength_pTsorted.push_back(clusters_coreshowerlength.at(index));
                clusters_maxlayer_pTsorted.push_back(clusters_maxlayer.at(index));

            }

            // REMOVE DUPLICATED SUPERCLUSTERS

            float new_pt_leading = clusters_pTsorted.at(0).Pt();

            if(new_pt_leading == previous_pt_leading) 
               continue;

            previous_pt_leading = new_pt_leading;

            _supercl3d_cl3d_pt.push_back(clusters_pt_pTsorted);
            _supercl3d_cl3d_energy.push_back(clusters_E_pTsorted);
            _supercl3d_cl3d_eta.push_back(clusters_eta_pTsorted);
            _supercl3d_cl3d_phi.push_back(clusters_phi_pTsorted);
            _supercl3d_cl3d_bdteg.push_back(clusters_bdteg_pTsorted);
            _supercl3d_cl3d_coreshowerlength.push_back(clusters_coreshowerlength_pTsorted);
            _supercl3d_cl3d_maxlayer.push_back(clusters_maxlayer_pTsorted);

        }

        _n_supercl3ds = _supercl3d_cl3d_pt.size();

        // BUILD SUPERCLUSTER INFO

        for (int i_supercl3d=0; i_supercl3d<_n_supercl3ds; i_supercl3d++){

            vector<float> cl3d_pt = _supercl3d_cl3d_pt.at(i_supercl3d);
            vector<float> cl3d_energy = _supercl3d_cl3d_energy.at(i_supercl3d);
            vector<float> cl3d_eta = _supercl3d_cl3d_eta.at(i_supercl3d);
            vector<float> cl3d_phi = _supercl3d_cl3d_phi.at(i_supercl3d);

            int n_cl3d = cl3d_pt.size();
            _supercl3d_n_cl3d.push_back(n_cl3d);

            float pt_supercluster = 0;
            float e_supercluster = 0;

            _supercl3d_eta.push_back(cl3d_eta.at(0));
            _supercl3d_phi.push_back(cl3d_phi.at(0));
            
            for (int i_subcl3d=0; i_subcl3d<n_cl3d; i_subcl3d++){
                pt_supercluster += cl3d_pt.at(i_subcl3d);
                e_supercluster += cl3d_energy.at(i_subcl3d);
            }

            _supercl3d_pt.push_back(pt_supercluster);
            _supercl3d_energy.push_back(e_supercluster);

            if(n_cl3d>0) h_pT_cl1->Fill(cl3d_pt.at(0));
            if(n_cl3d>1) h_pT_cl2->Fill(cl3d_pt.at(1));
            if(n_cl3d>2) h_pT_cl3->Fill(cl3d_pt.at(2));
            if(n_cl3d>3) h_pT_cl4->Fill(cl3d_pt.at(3));

        }

        // MATCHING TO GEN TAUS

        for (int i_gentau=0; i_gentau<_gentau_n; i_gentau++){

            TLorentzVector gentau;
            gentau.SetPtEtaPhiM( (*_gentau_pt)[i_gentau], (*_gentau_eta)[i_gentau], (*_gentau_phi)[i_gentau], (*_gentau_mass)[i_gentau]);

            TLorentzVector gentauvis;
            gentauvis.SetPtEtaPhiM( (*_gentau_vis_pt)[i_gentau], (*_gentau_vis_eta)[i_gentau], (*_gentau_vis_phi)[i_gentau], (*_gentau_vis_mass)[i_gentau]);

            float matched_cl3d_n = 0;
            vector<float> matched_cl3d_pt;
            vector<float> matched_cl3d_energy;
            vector<float> matched_cl3d_eta;
            vector<float> matched_cl3d_phi;
            vector<float> matched_cl3d_bdteg;
            vector<int> matched_cl3d_coreshowerlength;
            vector<int> matched_cl3d_maxlayer;

            matched_cl3d_pt.clear();
            matched_cl3d_energy.clear();
            matched_cl3d_eta.clear();
            matched_cl3d_phi.clear();
            matched_cl3d_bdteg.clear();
            matched_cl3d_coreshowerlength.clear();
            matched_cl3d_maxlayer.clear();

            float dRmin = 1.0;

            for (int i_supercl3d=0; i_supercl3d<_n_supercl3ds; i_supercl3d++){

                TLorentzVector supercluster;
                supercluster.SetPtEtaPhiE( _supercl3d_pt[i_supercl3d], _supercl3d_eta[i_supercl3d], _supercl3d_phi[i_supercl3d], _supercl3d_energy[i_supercl3d]);    

                bool isMatch = false;

                float dR_gentauvis = supercluster.DeltaR(gentauvis);

                isMatch = (dR_gentauvis <= 0.3);

                if(dR_gentauvis<dRmin){

                    dRmin = dR_gentauvis;

                    matched_cl3d_n      = _supercl3d_n_cl3d.at(i_supercl3d);
                    matched_cl3d_pt     = _supercl3d_cl3d_pt.at(i_supercl3d);
                    matched_cl3d_energy = _supercl3d_cl3d_energy.at(i_supercl3d);
                    matched_cl3d_eta    = _supercl3d_cl3d_eta.at(i_supercl3d);
                    matched_cl3d_phi    = _supercl3d_cl3d_phi.at(i_supercl3d);
                    matched_cl3d_bdteg  = _supercl3d_cl3d_bdteg.at(i_supercl3d);
                    matched_cl3d_coreshowerlength = _supercl3d_cl3d_coreshowerlength.at(i_supercl3d);
                    matched_cl3d_maxlayer = _supercl3d_cl3d_maxlayer.at(i_supercl3d);

                }

            }

            _gentau_isMatched.push_back(dRmin<=0.3);

            if(dRmin<=0.3) {

                _gentau_numberMatchedcl3d.push_back(matched_cl3d_n);
                _gentau_PtMatchedcl3d.push_back(matched_cl3d_pt);
                _gentau_EtaMatchedcl3d.push_back(matched_cl3d_eta);
                _gentau_PhiMatchedcl3d.push_back(matched_cl3d_phi);
                _gentau_BDTegMatchedcl3d.push_back(matched_cl3d_bdteg);
                _gentau_CoreShowerLengthMatchedcl3d.push_back(matched_cl3d_coreshowerlength);
                _gentau_MaxLayerMatchedcl3d.push_back(matched_cl3d_maxlayer);

            }

        }

        out_tree->Fill();

    }

    out_file->cd();

    h_pT_cl1->Write();
    h_pT_cl2->Write();
    h_pT_cl3->Write();
    h_pT_cl4->Write();

    out_tree->Write();
    out_file->Close();

    return;

}

void test(int n_events = -1, TString pu = "0", float ptthr_seed = 4, float ptthr_sec = 2, float etawindow = 0.3, float phiwindow = 0.5){

  TString dir = "/data_CMS/cms/mperez/HGCal_data/May19/";

  //cout<<"pu"<<pu<<endl;

  int i_ptthr_seed = ptthr_seed;
  int i_ptthr_sec = ptthr_sec;
  int i_etawindow = 10*etawindow;
  int i_phiwindow = 10*phiwindow;

  TString s_ptthr_seed = to_string(i_ptthr_seed);
  TString s_ptthr_sec  = to_string(i_ptthr_sec);
  TString s_etawindow  = to_string(i_etawindow);
  TString s_phiwindow  = to_string(i_phiwindow);

  TString infile = dir+"skimmed/NTuple_ZTT_PU"+pu+"_skimmed.root";
  TString outfile = dir+"clustered/NTuple_ZTT_PU"+pu+"_clustered_pTseed"+s_ptthr_seed+"_pTsec"+s_ptthr_sec+"_eta"+s_etawindow+"_phi"+s_phiwindow+".root";

  cluster_tree(infile, outfile, n_events, ptthr_seed, ptthr_sec, etawindow, phiwindow);

}
