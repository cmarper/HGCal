#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>


bool pT_comparison_pairs(pair<int,TLorentzVector> pair1, pair<int,TLorentzVector> pair2){

  return (pair1.second).Pt()>(pair2.second).Pt();

}

bool BDT_comparison_pairs(pair<int,float> pair1, pair<int,float> pair2){

  return pair1.second>pair2.second;

}



void sort_tree( TString filein, TString fileout, int nevents = -1){

	TFile* out_file = TFile::Open(fileout);
  	/*if(out_file!=0){
      cout<<fileout<<" already exists, please delete it before converting again"<<endl;
      return;
	}*/

	out_file = TFile::Open(fileout,"RECREATE");

	TChain * in_tree = new TChain("MatchedTree");	
	in_tree->Add(filein);

	TH2F* h2_DMall_icl3d_pticl3d=new TH2F("h2_DMall_icl3d_pticl3d","h2_DMall_icl3d_pticl3d",10,1,11,40,0,200);
	TH2F* h2_DM0_icl3d_pticl3d=new TH2F("h2_DM0_icl3d_pticl3d","h2_DM0_icl3d_pticl3d",10,1,11,40,0,200);
	TH2F* h2_DM1_icl3d_pticl3d=new TH2F("h2_DM1_icl3d_pticl3d","h2_DM1_icl3d_pticl3d",10,1,11,40,0,200);
	TH2F* h2_DM4_icl3d_pticl3d=new TH2F("h2_DM4_icl3d_pticl3d","h2_DM4_icl3d_pticl3d",10,1,11,40,0,200);
	TH2F* h2_DM5_icl3d_pticl3d=new TH2F("h2_DM5_icl3d_pticl3d","h2_DM5_icl3d_pticl3d",10,1,11,40,0,200);

	TH2F* h2_DMall_icl3d_pticl3d_o_ptMaxcl3d=new TH2F("h2_DMall_icl3d_pticl3d_o_ptMaxcl3d","h2_DMall_icl3d_pticl3d_o_ptMaxcl3d",9,2,11,40,0,2);
	TH2F* h2_DM0_icl3d_pticl3d_o_ptMaxcl3d=new TH2F("h2_DM0_icl3d_pticl3d_o_ptMaxcl3d","h2_DM0_icl3d_pticl3d_o_ptMaxcl3d",9,2,11,40,0,2);
	TH2F* h2_DM1_icl3d_pticl3d_o_ptMaxcl3d=new TH2F("h2_DM1_icl3d_pticl3d_o_ptMaxcl3d","h2_DM1_icl3d_pticl3d_o_ptMaxcl3d",9,2,11,40,0,2);
	TH2F* h2_DM4_icl3d_pticl3d_o_ptMaxcl3d=new TH2F("h2_DM4_icl3d_pticl3d_o_ptMaxcl3d","h2_DM4_icl3d_pticl3d_o_ptMaxcl3d",9,2,11,40,0,2);
	TH2F* h2_DM5_icl3d_pticl3d_o_ptMaxcl3d=new TH2F("h2_DM5_icl3d_pticl3d_o_ptMaxcl3d","h2_DM5_icl3d_pticl3d_o_ptMaxcl3d",9,2,11,40,0,2);

	TH2F* h2_DMall_bdtcl3d1_bdtcl3d2=new TH2F("h2_DMall_bdtcl3d1_bdtcl3d2","h2_DMall_bdtcl3d1_bdtcl3d2",15,-1,0.5,15,-1,0.5);
	TH2F* h2_DM0_bdtcl3d1_bdtcl3d2=new TH2F("h2_DM0_bdtcl3d1_bdtcl3d2","h2_DM0_bdtcl3d1_bdtcl3d2",15,-1,0.5,15,-1,0.5);
	TH2F* h2_DM1_bdtcl3d1_bdtcl3d2=new TH2F("h2_DM1_bdtcl3d1_bdtcl3d2","h2_DM1_bdtcl3d1_bdtcl3d2",15,-1,0.5,15,-1,0.5);
	TH2F* h2_DM4_bdtcl3d1_bdtcl3d2=new TH2F("h2_DM4_bdtcl3d1_bdtcl3d2","h2_DM4_bdtcl3d1_bdtcl3d2",15,-1,0.5,15,-1,0.5);
	TH2F* h2_DM5_bdtcl3d1_bdtcl3d2=new TH2F("h2_DM5_bdtcl3d1_bdtcl3d2","h2_DM5_bdtcl3d1_bdtcl3d2",15,-1,0.5,15,-1,0.5);

	TH2F* h2_DMall_icl3d_bdtegicl3d=new TH2F("h2_DMall_icl3d_bdtegicl3d","h2_DMall_icl3d_bdtegicl3d",10,1,11,15,-1,0.5);
	TH2F* h2_DM0_icl3d_bdtegicl3d=new TH2F("h2_DM0_icl3d_bdtegicl3d","h2_DM0_icl3d_bdtegicl3d",10,1,11,15,-1,0.5);
	TH2F* h2_DM1_icl3d_bdtegicl3d=new TH2F("h2_DM1_icl3d_bdtegicl3d","h2_DM1_icl3d_bdtegicl3d",10,1,11,15,-1,0.5);
	TH2F* h2_DM4_icl3d_bdtegicl3d=new TH2F("h2_DM4_icl3d_bdtegicl3d","h2_DM4_icl3d_bdtegicl3d",10,1,11,15,-1,0.5);
	TH2F* h2_DM5_icl3d_bdtegicl3d=new TH2F("h2_DM5_icl3d_bdtegicl3d","h2_DM5_icl3d_bdtegicl3d",10,1,11,15,-1,0.5);

	Long64_t nentries = in_tree->GetEntries();
	cout<<"nentries="<<in_tree->GetEntries()<<endl;
	if (nevents != -1) nentries = nevents;

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

	vector<bool>  *_gentau_isMatched;
	vector<int>   *_gentau_nMatchedcl3d;
	vector<float> *_gentau_PtMaxPtMatchedcl3d;
	vector<float> *_gentau_EtaMaxPtMatchedcl3d;
	vector<float> *_gentau_PhiMaxPtMatchedcl3d;
	vector<float> *_gentau_iMaxPtMatchedcl3d;
	vector<float> *_gentau_PtTotMatchedcl3d;

	int _cl3d_n;

	vector<float>  *_cl3d_pt;
	vector<float>  *_cl3d_energy;
	vector<float>  *_cl3d_eta;
	vector<float>  *_cl3d_phi;

	vector<float>  *_cl3d_bdteg;

	vector<int>    *_cl3d_decayMode;

	vector<bool>   *_cl3d_isMatched;
	vector<int>    *_cl3d_iMatchedgentau;
	vector<bool>   *_cl3d_isMaxPtcl3d;
	vector<float>  *_cl3d_dRtoMaxPtcl3d;

	in_tree->SetBranchAddress("gentau_n",&_gentau_n);

	in_tree->SetBranchAddress("gentau_pt",&_gentau_pt);
	in_tree->SetBranchAddress("gentau_eta",&_gentau_eta);
	in_tree->SetBranchAddress("gentau_phi",&_gentau_phi);
	in_tree->SetBranchAddress("gentau_energy",&_gentau_energy);
	in_tree->SetBranchAddress("gentau_mass",&_gentau_mass);

	in_tree->SetBranchAddress("gentau_vis_pt",&_gentau_vis_pt);
	in_tree->SetBranchAddress("gentau_vis_eta",&_gentau_vis_eta);
	in_tree->SetBranchAddress("gentau_vis_phi",&_gentau_vis_phi);
	in_tree->SetBranchAddress("gentau_vis_energy",&_gentau_vis_energy);
	in_tree->SetBranchAddress("gentau_vis_mass",&_gentau_vis_mass);

	in_tree->SetBranchAddress("gentau_decayMode",&_gentau_decayMode);

	in_tree->SetBranchAddress("gentau_isMatched",&_gentau_isMatched);
	in_tree->SetBranchAddress("gentau_nMatchedcl3d",&_gentau_nMatchedcl3d);
	in_tree->SetBranchAddress("gentau_PtMaxPtMatchedcl3d",&_gentau_PtMaxPtMatchedcl3d);
	in_tree->SetBranchAddress("gentau_EtaMaxPtMatchedcl3d",&_gentau_EtaMaxPtMatchedcl3d);
	in_tree->SetBranchAddress("gentau_PhiMaxPtMatchedcl3d",&_gentau_PhiMaxPtMatchedcl3d);
	in_tree->SetBranchAddress("gentau_iMaxPtMatchedcl3d",&_gentau_iMaxPtMatchedcl3d);
	in_tree->SetBranchAddress("gentau_PtTotMatchedcl3d",&_gentau_PtTotMatchedcl3d);

	in_tree->SetBranchAddress("cl3d_n",&_cl3d_n);

	in_tree->SetBranchAddress("cl3d_pt",&_cl3d_pt);
	in_tree->SetBranchAddress("cl3d_energy",&_cl3d_energy);
	in_tree->SetBranchAddress("cl3d_eta",&_cl3d_eta);
	in_tree->SetBranchAddress("cl3d_phi",&_cl3d_phi);

	in_tree->SetBranchAddress("cl3d_bdteg",&_cl3d_bdteg);

	in_tree->SetBranchAddress("cl3d_decayMode",&_cl3d_decayMode);

	in_tree->SetBranchAddress("cl3d_isMatched",&_cl3d_isMatched);
	in_tree->SetBranchAddress("cl3d_iMatchedgentau",&_cl3d_iMatchedgentau);
	in_tree->SetBranchAddress("cl3d_isMaxPtcl3d",&_cl3d_isMaxPtcl3d);
	in_tree->SetBranchAddress("cl3d_dRtoMaxPtcl3d",&_cl3d_dRtoMaxPtcl3d);


	// new branches

	TTree* out_tree=in_tree->GetTree()->CloneTree(0);
	out_tree->SetNameTitle("SortedTree","SortedTree");

	int _cl3d_matched_n;

	vector<float>  _cl3d_matched_pt;
	vector<float>  _cl3d_matched_energy;
	vector<float>  _cl3d_matched_eta;
	vector<float>  _cl3d_matched_phi;
	vector<float>  _cl3d_matched_bdteg;
	vector<int>    _cl3d_matched_decayMode;
	vector<bool>   _cl3d_matched_isMaxPtcl3d;
	vector<int>    _cl3d_matched_iMatchedgentau;

	vector<float>  _cl3d_matched_ptsort_pt;
	vector<float>  _cl3d_matched_ptsort_energy;
	vector<float>  _cl3d_matched_ptsort_eta;
	vector<float>  _cl3d_matched_ptsort_phi;
	vector<float>  _cl3d_matched_ptsort_bdteg;
	vector<int>    _cl3d_matched_ptsort_decayMode;
	vector<bool>   _cl3d_matched_ptsort_isMaxPtcl3d;
	vector<int>    _cl3d_matched_ptsort_iMatchedgentau;

	vector<float>  _cl3d_matched_bdtegsort_pt;
	vector<float>  _cl3d_matched_bdtegsort_energy;
	vector<float>  _cl3d_matched_bdtegsort_eta;
	vector<float>  _cl3d_matched_bdtegsort_phi;
	vector<float>  _cl3d_matched_bdtegsort_bdteg;
	vector<int>    _cl3d_matched_bdtegsort_decayMode;
	vector<bool>   _cl3d_matched_bdtegsort_isMaxPtcl3d;
	vector<int>    _cl3d_matched_bdtegsort_iMatchedgentau;

	out_tree->Branch("cl3d_matched_n",&_cl3d_matched_n);

	out_tree->Branch("cl3d_matched_pt",&_cl3d_matched_pt);
	out_tree->Branch("cl3d_matched_energy",&_cl3d_matched_energy);
	out_tree->Branch("cl3d_matched_eta",&_cl3d_matched_eta);
	out_tree->Branch("cl3d_matched_phi",&_cl3d_matched_phi);
	out_tree->Branch("cl3d_matched_bdteg",&_cl3d_matched_bdteg);
	out_tree->Branch("cl3d_matched_decayMode",&_cl3d_matched_decayMode);
	out_tree->Branch("cl3d_matched_isMaxPtcl3d",&_cl3d_matched_isMaxPtcl3d);
	out_tree->Branch("cl3d_matched_iMatchedgentau",&_cl3d_matched_iMatchedgentau);

	out_tree->Branch("cl3d_matched_ptsort_pt",&_cl3d_matched_ptsort_pt);
	out_tree->Branch("cl3d_matched_ptsort_energy",&_cl3d_matched_ptsort_energy);
	out_tree->Branch("cl3d_matched_ptsort_eta",&_cl3d_matched_ptsort_eta);
	out_tree->Branch("cl3d_matched_ptsort_phi",&_cl3d_matched_ptsort_phi);
	out_tree->Branch("cl3d_matched_ptsort_bdteg",&_cl3d_matched_ptsort_bdteg);
	out_tree->Branch("cl3d_matched_ptsort_decayMode",&_cl3d_matched_ptsort_decayMode);
	out_tree->Branch("cl3d_matched_ptsort_isMaxPtcl3d",&_cl3d_matched_ptsort_isMaxPtcl3d);
	out_tree->Branch("cl3d_matched_ptsort_iMatchedgentau",&_cl3d_matched_ptsort_iMatchedgentau);

	out_tree->Branch("cl3d_matched_bdtegsort_pt",&_cl3d_matched_bdtegsort_pt);
	out_tree->Branch("cl3d_matched_bdtegsort_energy",&_cl3d_matched_bdtegsort_energy);
	out_tree->Branch("cl3d_matched_bdtegsort_eta",&_cl3d_matched_bdtegsort_eta);
	out_tree->Branch("cl3d_matched_bdtegsort_phi",&_cl3d_matched_bdtegsort_phi);
	out_tree->Branch("cl3d_matched_bdtegsort_bdteg",&_cl3d_matched_bdtegsort_bdteg);
	out_tree->Branch("cl3d_matched_bdtegsort_decayMode",&_cl3d_matched_bdtegsort_decayMode);
	out_tree->Branch("cl3d_matched_bdtegsort_isMaxPtcl3d",&_cl3d_matched_bdtegsort_isMaxPtcl3d);
	out_tree->Branch("cl3d_matched_bdtegsort_iMatchedgentau",&_cl3d_matched_bdtegsort_iMatchedgentau);


	for (int i=0;i<nentries;i++) {

		if(i%1000==0) cout<<"i="<<i<<endl;

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
		_gentau_isMatched = 0;
		_gentau_nMatchedcl3d = 0;
		_gentau_PtMaxPtMatchedcl3d = 0;
		_gentau_EtaMaxPtMatchedcl3d = 0;
		_gentau_PhiMaxPtMatchedcl3d = 0;
		_gentau_iMaxPtMatchedcl3d = 0;
		_gentau_PtTotMatchedcl3d = 0;

		_cl3d_n = 0;
		_cl3d_pt = 0;
		_cl3d_energy = 0;
		_cl3d_eta = 0;
		_cl3d_phi = 0;
		_cl3d_bdteg = 0;
		_cl3d_decayMode = 0;
		_cl3d_isMatched = 0;
		_cl3d_iMatchedgentau = 0;
		_cl3d_isMaxPtcl3d = 0;
		_cl3d_dRtoMaxPtcl3d = 0;

		_cl3d_matched_n = 0;


		int entry_ok = in_tree->GetEntry(i);	
		if(entry_ok<0) continue;


		for(int i_gentau=0; i_gentau<_gentau_n; i_gentau++){

			_cl3d_matched_pt.clear();
			_cl3d_matched_energy.clear();
			_cl3d_matched_eta.clear();
			_cl3d_matched_phi.clear();
			_cl3d_matched_bdteg.clear();
			_cl3d_matched_decayMode.clear();
			_cl3d_matched_isMaxPtcl3d.clear();
			_cl3d_matched_iMatchedgentau.clear();

			_cl3d_matched_ptsort_pt.clear();
			_cl3d_matched_ptsort_energy.clear();
			_cl3d_matched_ptsort_eta.clear();
			_cl3d_matched_ptsort_phi.clear();
			_cl3d_matched_ptsort_bdteg.clear();
			_cl3d_matched_ptsort_decayMode.clear();
			_cl3d_matched_ptsort_isMaxPtcl3d.clear();
			_cl3d_matched_ptsort_iMatchedgentau.clear();

			_cl3d_matched_bdtegsort_pt.clear();
			_cl3d_matched_bdtegsort_energy.clear();
			_cl3d_matched_bdtegsort_eta.clear();
			_cl3d_matched_bdtegsort_phi.clear();
			_cl3d_matched_bdtegsort_bdteg.clear();
			_cl3d_matched_bdtegsort_decayMode.clear();
			_cl3d_matched_bdtegsort_isMaxPtcl3d.clear();
			_cl3d_matched_bdtegsort_iMatchedgentau.clear();


			/////////////////////////////////
			//////// matched cl3d ///////////
			/////////////////////////////////

			for(unsigned int i_cl3d=0; i_cl3d<(*_cl3d_iMatchedgentau).size(); i_cl3d++){

				if (!(*_cl3d_isMatched)[i_cl3d]) continue;
				if ((*_cl3d_iMatchedgentau)[i_cl3d] != i_gentau) continue;
				if ((*_cl3d_pt)[i_cl3d]<0.00001 || (*_cl3d_pt)[i_cl3d]>1e+10) continue;

				_cl3d_matched_pt.push_back((*_cl3d_pt)[i_cl3d]);
				_cl3d_matched_energy.push_back((*_cl3d_energy)[i_cl3d]);
				_cl3d_matched_eta.push_back((*_cl3d_eta)[i_cl3d]);
				_cl3d_matched_phi.push_back((*_cl3d_phi)[i_cl3d]);
				_cl3d_matched_bdteg.push_back((*_cl3d_bdteg)[i_cl3d]);
				_cl3d_matched_decayMode.push_back((*_cl3d_decayMode)[i_cl3d]);
				_cl3d_matched_isMaxPtcl3d.push_back((*_cl3d_isMaxPtcl3d)[i_cl3d]);
				_cl3d_matched_iMatchedgentau.push_back((*_cl3d_iMatchedgentau)[i_cl3d]);

			}

			_cl3d_matched_n = _cl3d_matched_pt.size();

			/////////////////////////////////
			///////// pT sorting ////////////
			/////////////////////////////////

			vector< pair<int,TLorentzVector> > matchedcl3d_pairs;
			TLorentzVector matchedcl3d;

			for(int i_matched_cl3d=0; i_matched_cl3d<_cl3d_matched_n; i_matched_cl3d++){

				matchedcl3d.SetPtEtaPhiM( _cl3d_matched_pt[i_matched_cl3d], _cl3d_matched_eta[i_matched_cl3d], _cl3d_matched_phi[i_matched_cl3d], 0 );
				pair<int,TLorentzVector> matchedcl3d_pair = make_pair(i_matched_cl3d,matchedcl3d);
				matchedcl3d_pairs.push_back(matchedcl3d_pair);

			}

			sort(matchedcl3d_pairs.begin(), matchedcl3d_pairs.end(), pT_comparison_pairs);


			for(int i_matched_cl3d=0; i_matched_cl3d<_cl3d_matched_n; i_matched_cl3d++){

				int i_mat_cl3d = matchedcl3d_pairs[i_matched_cl3d].first;
				TLorentzVector mat_cl3d = matchedcl3d_pairs[i_matched_cl3d].second;

				_cl3d_matched_ptsort_pt.push_back(mat_cl3d.Pt());
				_cl3d_matched_ptsort_eta.push_back(mat_cl3d.Eta());
				_cl3d_matched_ptsort_phi.push_back(mat_cl3d.Phi());
				_cl3d_matched_ptsort_energy.push_back(_cl3d_matched_energy[i_mat_cl3d]);
				_cl3d_matched_ptsort_bdteg.push_back(_cl3d_matched_bdteg[i_mat_cl3d]);
				_cl3d_matched_ptsort_decayMode.push_back(_cl3d_matched_decayMode[i_mat_cl3d]);
				_cl3d_matched_ptsort_isMaxPtcl3d.push_back(_cl3d_matched_isMaxPtcl3d[i_mat_cl3d]);
				_cl3d_matched_ptsort_iMatchedgentau.push_back(_cl3d_matched_iMatchedgentau[i_mat_cl3d]);

			}


			/////////////////////////////////
			/////// BDT EG sorting //////////
			/////////////////////////////////

			vector< pair<int,float> > i_cl3d_BDT_pairs;

			for(int i_matched_cl3d=0; i_matched_cl3d<_cl3d_matched_n; i_matched_cl3d++){

				float BDT = _cl3d_matched_bdteg[i_matched_cl3d];
				pair<int,float> BDT_pair = make_pair(i_matched_cl3d,BDT);
				i_cl3d_BDT_pairs.push_back(BDT_pair);

			}

			sort(i_cl3d_BDT_pairs.begin(), i_cl3d_BDT_pairs.end(), BDT_comparison_pairs);

			for(int i_matched_cl3d=0; i_matched_cl3d<_cl3d_matched_n; i_matched_cl3d++){

				int i_mat_cl3d = i_cl3d_BDT_pairs[i_matched_cl3d].first;
				float bdt_cl3d = i_cl3d_BDT_pairs[i_matched_cl3d].second;

				_cl3d_matched_bdtegsort_pt.push_back(_cl3d_matched_pt[i_mat_cl3d]);
				_cl3d_matched_bdtegsort_eta.push_back(_cl3d_matched_eta[i_mat_cl3d]);
				_cl3d_matched_bdtegsort_phi.push_back(_cl3d_matched_phi[i_mat_cl3d]);
				_cl3d_matched_bdtegsort_energy.push_back(_cl3d_matched_energy[i_mat_cl3d]);
				_cl3d_matched_bdtegsort_bdteg.push_back(_cl3d_matched_bdteg[i_mat_cl3d]);
				_cl3d_matched_bdtegsort_decayMode.push_back(_cl3d_matched_decayMode[i_mat_cl3d]);
				_cl3d_matched_bdtegsort_isMaxPtcl3d.push_back(_cl3d_matched_isMaxPtcl3d[i_mat_cl3d]);
				_cl3d_matched_bdtegsort_iMatchedgentau.push_back(_cl3d_matched_iMatchedgentau[i_mat_cl3d]);

			}


			/////////////////////////////////
			//////// pT histograms //////////
			/////////////////////////////////


			for(unsigned int i_matched_cl3d=0; i_matched_cl3d<_cl3d_matched_ptsort_pt.size(); i_matched_cl3d++){

				h2_DMall_icl3d_pticl3d->Fill(i_matched_cl3d+1,_cl3d_matched_ptsort_pt[i_matched_cl3d]);

				float pt_rel = _cl3d_matched_ptsort_pt[i_matched_cl3d]/_cl3d_matched_ptsort_pt[0];
				h2_DMall_icl3d_pticl3d_o_ptMaxcl3d->Fill(i_matched_cl3d+1,pt_rel);

				if(_cl3d_matched_ptsort_decayMode[i_matched_cl3d]==0) {
					h2_DM0_icl3d_pticl3d->Fill(i_matched_cl3d+1,_cl3d_matched_ptsort_pt[i_matched_cl3d]);
					h2_DM0_icl3d_pticl3d_o_ptMaxcl3d->Fill(i_matched_cl3d+1,pt_rel);
				}

				else if(_cl3d_matched_ptsort_decayMode[i_matched_cl3d]==1) {
					h2_DM1_icl3d_pticl3d->Fill(i_matched_cl3d+1,_cl3d_matched_ptsort_pt[i_matched_cl3d]);
					h2_DM1_icl3d_pticl3d_o_ptMaxcl3d->Fill(i_matched_cl3d+1,pt_rel);
				}

				else if(_cl3d_matched_ptsort_decayMode[i_matched_cl3d]==4) {
					h2_DM4_icl3d_pticl3d->Fill(i_matched_cl3d+1,_cl3d_matched_ptsort_pt[i_matched_cl3d]);
					h2_DM4_icl3d_pticl3d_o_ptMaxcl3d->Fill(i_matched_cl3d+1,pt_rel);
				}

				else if(_cl3d_matched_ptsort_decayMode[i_matched_cl3d]==5) {
					h2_DM5_icl3d_pticl3d->Fill(i_matched_cl3d+1,_cl3d_matched_ptsort_pt[i_matched_cl3d]);
					h2_DM5_icl3d_pticl3d_o_ptMaxcl3d->Fill(i_matched_cl3d+1,pt_rel);
				}

			}

			if(_cl3d_matched_ptsort_pt.size()>1){

				h2_DMall_bdtcl3d1_bdtcl3d2->Fill(_cl3d_matched_ptsort_bdteg[0],_cl3d_matched_ptsort_bdteg[1]);

				if (_cl3d_matched_ptsort_decayMode[0] == 0 && _cl3d_matched_ptsort_decayMode[1] == 0 ) h2_DM0_bdtcl3d1_bdtcl3d2->Fill(_cl3d_matched_ptsort_bdteg[0],_cl3d_matched_ptsort_bdteg[1]);
				else if (_cl3d_matched_ptsort_decayMode[0] == 1 && _cl3d_matched_ptsort_decayMode[1] == 1 ) h2_DM1_bdtcl3d1_bdtcl3d2->Fill(_cl3d_matched_ptsort_bdteg[0],_cl3d_matched_ptsort_bdteg[1]);
				else if (_cl3d_matched_ptsort_decayMode[0] == 4 && _cl3d_matched_ptsort_decayMode[1] == 4 ) h2_DM4_bdtcl3d1_bdtcl3d2->Fill(_cl3d_matched_ptsort_bdteg[0],_cl3d_matched_ptsort_bdteg[1]);
				else if (_cl3d_matched_ptsort_decayMode[0] == 5 && _cl3d_matched_ptsort_decayMode[1] == 5 ) h2_DM5_bdtcl3d1_bdtcl3d2->Fill(_cl3d_matched_ptsort_bdteg[0],_cl3d_matched_ptsort_bdteg[1]);

			}


			/////////////////////////////////
			//////// BDT histograms /////////
			/////////////////////////////////

			for(unsigned int i_matched_cl3d=0; i_matched_cl3d<_cl3d_matched_bdtegsort_pt.size(); i_matched_cl3d++){

				h2_DMall_icl3d_bdtegicl3d->Fill(i_matched_cl3d+1,_cl3d_matched_bdtegsort_bdteg[i_matched_cl3d]);

				if(_cl3d_matched_bdtegsort_decayMode[i_matched_cl3d]==0) h2_DM0_icl3d_bdtegicl3d->Fill(i_matched_cl3d+1,_cl3d_matched_bdtegsort_bdteg[i_matched_cl3d]);
				else if (_cl3d_matched_bdtegsort_decayMode[i_matched_cl3d]==1) h2_DM1_icl3d_bdtegicl3d->Fill(i_matched_cl3d+1,_cl3d_matched_bdtegsort_bdteg[i_matched_cl3d]);
				else if (_cl3d_matched_bdtegsort_decayMode[i_matched_cl3d]==4) h2_DM4_icl3d_bdtegicl3d->Fill(i_matched_cl3d+1,_cl3d_matched_bdtegsort_bdteg[i_matched_cl3d]);
				else if (_cl3d_matched_bdtegsort_decayMode[i_matched_cl3d]==5) h2_DM5_icl3d_bdtegicl3d->Fill(i_matched_cl3d+1,_cl3d_matched_bdtegsort_bdteg[i_matched_cl3d]);

			}

		}

		out_tree->Fill();

	}

	out_file->cd();
    out_tree->Write();

    h2_DMall_icl3d_pticl3d->Write();
    h2_DM0_icl3d_pticl3d->Write();
    h2_DM1_icl3d_pticl3d->Write();
    h2_DM4_icl3d_pticl3d->Write();
    h2_DM5_icl3d_pticl3d->Write();
    
    h2_DMall_icl3d_pticl3d_o_ptMaxcl3d->Write();
    h2_DM0_icl3d_pticl3d_o_ptMaxcl3d->Write();
    h2_DM1_icl3d_pticl3d_o_ptMaxcl3d->Write();
    h2_DM4_icl3d_pticl3d_o_ptMaxcl3d->Write();
    h2_DM5_icl3d_pticl3d_o_ptMaxcl3d->Write();

    h2_DMall_bdtcl3d1_bdtcl3d2->Write();
    h2_DM0_bdtcl3d1_bdtcl3d2->Write();
    h2_DM1_bdtcl3d1_bdtcl3d2->Write();
    h2_DM4_bdtcl3d1_bdtcl3d2->Write();
    h2_DM5_bdtcl3d1_bdtcl3d2->Write();

    h2_DMall_icl3d_bdtegicl3d->Write();
	h2_DM0_icl3d_bdtegicl3d->Write();
	h2_DM1_icl3d_bdtegicl3d->Write();
	h2_DM4_icl3d_bdtegicl3d->Write();
	h2_DM5_icl3d_bdtegicl3d->Write();

    out_file->Close();

	return;

}


void test(int n_events){

  TString infile = "/data_CMS/cms/mperez/HGCal_data/L1taus/NTuple_ZTT_PU0_matched_0p3_BDTeg.root";
  TString outfile = "/data_CMS/cms/mperez/HGCal_data/L1taus/NTuple_ZTT_PU0_sorted_0p3_BDTeg.root";

  sort_tree(infile, outfile, n_events);

}