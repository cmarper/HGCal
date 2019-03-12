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

void match_tree( TString filein, TString fileout, int nevents = -1, float dR_max = 0.3){
	
	TFile* out_file = TFile::Open(fileout);
  	if(out_file!=0){
      cout<<fileout<<" already exists, please delete it before converting again"<<endl;
      return;
	}

	out_file = TFile::Open(fileout,"RECREATE");

	TChain * in_tree = new TChain("SkimmedTree");	
	in_tree->Add(filein);

	TH1F* h_dR_match=new TH1F("h_dR_match","h_dR_match",10,0,1);
	h_dR_match->Fill(dR_max);

	TH1F* h_n_total_tau=new TH1F("h_n_tot_tau","h_n_tot_tau",1000000,0,1000000);
	TH1F* h_n_matched_tau=new TH1F("h_n_matched_tau","h_n_matched_tau",1000000,0,1000000);
	TH1F* h_match_eff=new TH1F("h_match_eff","h_match_eff",100,0,1);

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

	int _cl3d_n;

	vector<float>  *_cl3d_pt;
	vector<float>  *_cl3d_energy;
	vector<float>  *_cl3d_eta;
	vector<float>  *_cl3d_phi;


	in_tree->SetBranchAddress("gentau_n",&_gentau_n);

	in_tree->SetBranchAddress("gentau_pt",		&_gentau_pt);
	in_tree->SetBranchAddress("gentau_eta",		&_gentau_eta);
	in_tree->SetBranchAddress("gentau_phi",		&_gentau_phi);
	in_tree->SetBranchAddress("gentau_energy",	&_gentau_energy);
	in_tree->SetBranchAddress("gentau_mass",	&_gentau_mass);

	in_tree->SetBranchAddress("gentau_vis_pt",		&_gentau_vis_pt);
	in_tree->SetBranchAddress("gentau_vis_eta",		&_gentau_vis_eta);
	in_tree->SetBranchAddress("gentau_vis_phi",		&_gentau_vis_phi);
	in_tree->SetBranchAddress("gentau_vis_energy",	&_gentau_vis_energy);
	in_tree->SetBranchAddress("gentau_vis_mass",	&_gentau_vis_mass);

	in_tree->SetBranchAddress("gentau_decayMode",	&_gentau_decayMode);

	in_tree->SetBranchAddress("cl3d_n",&_cl3d_n);

	in_tree->SetBranchAddress("cl3d_pt",	 &_cl3d_pt);
	in_tree->SetBranchAddress("cl3d_energy", &_cl3d_energy);
	in_tree->SetBranchAddress("cl3d_eta",	 &_cl3d_eta);
	in_tree->SetBranchAddress("cl3d_phi",	 &_cl3d_phi);


	// new branches

	TTree* out_tree=in_tree->GetTree()->CloneTree(0);
	out_tree->SetNameTitle("MatchedTree","MatchedTree");

	vector<float>  _dR_cl3d_gentau;

	vector<bool>   			_gentau_isMatched;
	vector<vector<float> >	_gentau_iMatchedcl3d;
	vector<int>   			_gentau_nMatchedcl3d;
	vector<float>   		_gentau_PtMaxMatchedcl3d;
	vector<float>   		_gentau_iPtMaxMatchedcl3d;
	vector<float>   		_gentau_PtTotMatchedcl3d;

	vector<bool>   _cl3d_isMatched;
	vector<int>    _cl3d_iMatchedgentau;

	out_tree->Branch("dR_cl3d_gentau",&_dR_cl3d_gentau);

	out_tree->Branch("gentau_isMatched",			&_gentau_isMatched);
	out_tree->Branch("gentau_iMatchedcl3d",			&_gentau_iMatchedcl3d);
	out_tree->Branch("gentau_nMatchedcl3d",			&_gentau_nMatchedcl3d);
	out_tree->Branch("gentau_PtMaxMatchedcl3d",		&_gentau_PtMaxMatchedcl3d);
	out_tree->Branch("gentau_iPtMaxMatchedcl3d",	&_gentau_iPtMaxMatchedcl3d);
	out_tree->Branch("gentau_PtTotMatchedcl3d",		&_gentau_PtTotMatchedcl3d);

	out_tree->Branch("cl3d_isMatched", &_cl3d_isMatched);
	out_tree->Branch("cl3d_iMatchedgentau", &_cl3d_iMatchedgentau);

	float n_total_gentau = 0;
	float n_total_gentau_matched = 0;

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

		// new branches

		_dR_cl3d_gentau.clear();

		_gentau_isMatched.clear();
		_gentau_iMatchedcl3d.clear();
		_gentau_nMatchedcl3d.clear();
		_gentau_PtMaxMatchedcl3d.clear();
		_gentau_iPtMaxMatchedcl3d.clear();
		_gentau_PtTotMatchedcl3d.clear();

		_cl3d_isMatched.clear();
		_cl3d_iMatchedgentau.clear();


		int entry_ok = in_tree->GetEntry(i);	
		if(entry_ok<0) continue;
		
	
		for(int i_gentau=0; i_gentau<_gentau_n; i_gentau++){

			n_total_gentau++;

			TLorentzVector matching_tau;
			matching_tau.SetPtEtaPhiM( (*_gentau_vis_pt)[i_gentau], (*_gentau_vis_eta)[i_gentau], (*_gentau_vis_phi)[i_gentau], (*_gentau_vis_mass)[i_gentau]);

			int n_match = 0;

			vector<float> matched_cl3ds;
			matched_cl3ds.clear();

			float maxPt = -1;
			int i_maxPt = -1;
			float totPt = 0;

			for (int i_cl3d=0; i_cl3d<_cl3d_n; i_cl3d++){

				TLorentzVector matching_cl3d;
				matching_cl3d.SetPtEtaPhiM( (*_cl3d_pt)[i_cl3d], (*_cl3d_eta)[i_cl3d], (*_cl3d_phi)[i_cl3d], 0);

				float _myidr = sqrt( ((*_gentau_vis_eta)[i_gentau]-(*_cl3d_eta)[i_cl3d])*((*_gentau_vis_eta)[i_gentau]-(*_cl3d_eta)[i_cl3d]) + ((*_gentau_vis_phi)[i_gentau]-(*_cl3d_phi)[i_cl3d])*((*_gentau_vis_phi)[i_gentau]-(*_cl3d_phi)[i_cl3d]) );

				float _i_dR = matching_tau.DeltaR(matching_cl3d);
				_dR_cl3d_gentau.push_back(_i_dR);

				bool _isMatch = _i_dR <= dR_max;
				_cl3d_isMatched.push_back(_isMatch);

				if (_isMatch){	

					_cl3d_iMatchedgentau.push_back(i_gentau);
					matched_cl3ds.push_back(i_cl3d);
					n_match++;

					totPt = totPt + matching_cl3d.Pt();

					if (matching_cl3d.Pt()>maxPt){
						maxPt = matching_cl3d.Pt();
						i_maxPt = i_cl3d;
					}

				}

				else if (!_isMatch){
					_cl3d_iMatchedgentau.push_back(-1);	
				}

			}


			if (n_match>0) {

				_gentau_isMatched.push_back(true);
				_gentau_iMatchedcl3d.push_back(matched_cl3ds);
				_gentau_PtMaxMatchedcl3d.push_back(maxPt);
				_gentau_iPtMaxMatchedcl3d.push_back(i_maxPt);
				_gentau_PtTotMatchedcl3d.push_back(totPt);
				n_total_gentau_matched++;
			
			}

			else if (n_match==0) {
				_gentau_isMatched.push_back(false);
			}

			_gentau_nMatchedcl3d.push_back(n_match);

		}

		out_tree->Fill();

	}

	h_n_total_tau->Fill(n_total_gentau);
	h_n_matched_tau->Fill(n_total_gentau_matched);
	h_match_eff->Fill(n_total_gentau_matched/n_total_gentau);

	out_file->cd();

	h_n_total_tau->Write();
	h_n_matched_tau->Write();
	h_dR_match->Write();

    out_tree->Write();
    out_file->Close();

	return;

}


void test(int n_events){

  TString infile = "/data_CMS/cms/mperez/HGCal_data/L1taus/NTuple_ZTT_PU0_skimmed_BDTeg.root";
  TString outfile = "/data_CMS/cms/mperez/HGCal_data/L1taus/NTuple_ZTT_PU0_matched_0p3_BDTeg.root";

  match_tree(infile, outfile, n_events, 0.3);

}