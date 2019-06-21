/////////////////////////////////////////////////////////
///// HGCal L1 taus, C. Martin Perez, LLR, May 2019 /////
/////////////////////////////////////////////////////////

#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <iostream>

using namespace std;

void match_tree( TString filein, TString fileout, int nevents = -1, float dR_max = 0.3){
	
	TFile* out_file = TFile::Open(fileout);
  	/*if(out_file!=0){
      cout<<fileout<<" already exists, please delete it before converting again"<<endl;
      return;
	}*/

	out_file = TFile::Open(fileout,"RECREATE");

	TChain * in_tree = new TChain("SkimmedTree");	
	in_tree->Add(filein);

	TH1F* h_dR_match=new TH1F("h_dR_match","h_dR_match",100,0,1);
	h_dR_match->Fill(dR_max);


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
	vector<int>  *_cl3d_coreshowerlength;
        vector<int>  *_cl3d_maxlayer;

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
	in_tree->SetBranchAddress("cl3d_bdteg",	 &_cl3d_bdteg);
	in_tree->SetBranchAddress("cl3d_coreshowerlength", &_cl3d_coreshowerlength);
	in_tree->SetBranchAddress("cl3d_maxlayer", &_cl3d_maxlayer);

	// new branches

	TTree* out_tree=in_tree->GetTree()->CloneTree(0);
	out_tree->SetNameTitle("MatchedTree","MatchedTree");

	vector<bool>   		_gentau_isMatched;
	vector<vector<float> >	_gentau_indicesMatchedcl3d;
	vector<int>   		_gentau_numberMatchedcl3d;

	vector<vector<float> >	_gentau_dRMatchedcl3d;
	vector<vector<float> >	_gentau_PtMatchedcl3d;
	vector<vector<float> >	_gentau_EtaMatchedcl3d;
	vector<vector<float> >	_gentau_PhiMatchedcl3d;
	vector<vector<float> >	_gentau_BDTegMatchedcl3d;
	vector<vector<int> >  _gentau_CoreshowerlengthMatchedcl3d;
	vector<vector<int> >  _gentau_MaxlayerMatchedcl3d;

	vector<bool>   _cl3d_isMatched;
	vector<int>    _cl3d_indexMatchedgentau;

	vector<float>  _cl3d_dRMatchedgentau;
	vector<float>  _cl3d_dRMatchedgentauvis;

	vector<float>  _cl3d_PtMatchedgentau;
	vector<float>  _cl3d_EtaMatchedgentau;
	vector<float>  _cl3d_PhiMatchedgentau;
	vector<float>  _cl3d_EnergyMatchedgentau;
	vector<float>  _cl3d_MassMatchedgentau;
	vector<float>  _cl3d_decayModeMatchedgentau;

	vector<float>  _cl3d_PtMatchedgentauvis;
	vector<float>  _cl3d_EtaMatchedgentauvis;
	vector<float>  _cl3d_PhiMatchedgentauvis;
	vector<float>  _cl3d_EnergyMatchedgentauvis;
	vector<float>  _cl3d_MassMatchedgentauvis;


	out_tree->Branch("gentau_isMatched",		&_gentau_isMatched);
	out_tree->Branch("gentau_indicesMatchedcl3d",	&_gentau_indicesMatchedcl3d);
	out_tree->Branch("gentau_numberMatchedcl3d",	&_gentau_numberMatchedcl3d);

	out_tree->Branch("gentau_dRMatchedcl3d",	&_gentau_dRMatchedcl3d);
	out_tree->Branch("gentau_PtMatchedcl3d",	&_gentau_PtMatchedcl3d);
	out_tree->Branch("gentau_EtaMatchedcl3d",	&_gentau_EtaMatchedcl3d);
	out_tree->Branch("gentau_PhiMatchedcl3d",	&_gentau_PhiMatchedcl3d);
	out_tree->Branch("gentau_BDTegMatchedcl3d",	&_gentau_BDTegMatchedcl3d);
	out_tree->Branch("gentau_CoreshowerlengthMatchedcl3d",	&_gentau_CoreshowerlengthMatchedcl3d);
	out_tree->Branch("gentau_MaxlayerMatchedcl3d",	&_gentau_MaxlayerMatchedcl3d);

	out_tree->Branch("cl3d_isMatched", 		&_cl3d_isMatched);
	out_tree->Branch("cl3d_indexMatchedgentau", 	&_cl3d_indexMatchedgentau);

	out_tree->Branch("cl3d_dRMatchedgentau",	&_cl3d_dRMatchedgentau);
	out_tree->Branch("cl3d_dRMatchedgentauvis",	&_cl3d_dRMatchedgentauvis);
	
	out_tree->Branch("cl3d_PtMatchedgentau", 	&_cl3d_PtMatchedgentau);
	out_tree->Branch("cl3d_EtaMatchedgentau", 	&_cl3d_EtaMatchedgentau);
	out_tree->Branch("cl3d_PhiMatchedgentau",	&_cl3d_PhiMatchedgentau);
	out_tree->Branch("cl3d_EnergyMatchedgentau",	&_cl3d_EnergyMatchedgentau);
	out_tree->Branch("cl3d_MassMatchedgentau",	&_cl3d_MassMatchedgentau);
	out_tree->Branch("cl3d_decayModeMatchedgentau",	&_cl3d_decayModeMatchedgentau);

	out_tree->Branch("cl3d_PtMatchedgentauvis", 	&_cl3d_PtMatchedgentauvis);
	out_tree->Branch("cl3d_EtaMatchedgentauvis", 	&_cl3d_EtaMatchedgentauvis);
	out_tree->Branch("cl3d_PhiMatchedgentauvis",	&_cl3d_PhiMatchedgentauvis);
	out_tree->Branch("cl3d_EnergyMatchedgentauvis",	&_cl3d_EnergyMatchedgentauvis);
	out_tree->Branch("cl3d_MassMatchedgentauvis",	&_cl3d_MassMatchedgentauvis);


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

		// new branches

		_gentau_isMatched.clear();
		_gentau_indicesMatchedcl3d.clear();
		_gentau_numberMatchedcl3d.clear();

		_gentau_dRMatchedcl3d.clear();
		_gentau_PtMatchedcl3d.clear();
		_gentau_EtaMatchedcl3d.clear();
		_gentau_PhiMatchedcl3d.clear();
		_gentau_BDTegMatchedcl3d.clear();
		_gentau_CoreshowerlengthMatchedcl3d.clear();
		_gentau_MaxlayerMatchedcl3d.clear();
		
		_cl3d_isMatched.clear();
		_cl3d_indexMatchedgentau.clear();

		_cl3d_dRMatchedgentau.clear();
		_cl3d_dRMatchedgentauvis.clear();

		_cl3d_PtMatchedgentau.clear();
		_cl3d_EtaMatchedgentau.clear();
		_cl3d_PhiMatchedgentau.clear();
		_cl3d_EnergyMatchedgentau.clear();
		_cl3d_MassMatchedgentau.clear();
		_cl3d_decayModeMatchedgentau.clear();

		_cl3d_PtMatchedgentauvis.clear();
		_cl3d_EtaMatchedgentauvis.clear();
		_cl3d_PhiMatchedgentauvis.clear();
		_cl3d_EnergyMatchedgentauvis.clear();
		_cl3d_MassMatchedgentauvis.clear();


		int entry_ok = in_tree->GetEntry(i);	
		if(entry_ok<0) 
			continue;
	     
		//loop over gentaus

		for(int i_gentau=0; i_gentau<_gentau_n; i_gentau++){

			TLorentzVector matching_gentau;
			matching_gentau.SetPtEtaPhiM( (*_gentau_pt)[i_gentau], (*_gentau_eta)[i_gentau], (*_gentau_phi)[i_gentau], (*_gentau_mass)[i_gentau]);

			TLorentzVector matching_gentauvis;
			matching_gentauvis.SetPtEtaPhiM( (*_gentau_vis_pt)[i_gentau], (*_gentau_vis_eta)[i_gentau], (*_gentau_vis_phi)[i_gentau], (*_gentau_vis_mass)[i_gentau]);

			bool isMatch = false;
			int n_match = 0;

			vector<float> indicesMatchedcl3d;
			vector<float> dRMatchedcl3d;
			vector<float> PtMatchedcl3d;
			vector<float> EtaMatchedcl3d;
			vector<float> PhiMatchedcl3d;
			vector<float> BDTegMatchedcl3d;
			vector<int> CoreshowerlengthMatchedcl3d;
			vector<int> MaxlayerMatchedcl3d;

			indicesMatchedcl3d.clear();
			dRMatchedcl3d.clear();
			PtMatchedcl3d.clear();
			EtaMatchedcl3d.clear();
			PhiMatchedcl3d.clear();
			BDTegMatchedcl3d.clear();
			CoreshowerlengthMatchedcl3d.clear();
			MaxlayerMatchedcl3d.clear();
		
			// loop over cl3d

			for (int i_cl3d=0; i_cl3d<_cl3d_n; i_cl3d++){

				if ((*_cl3d_maxlayer)[i_cl3d]<6) continue;
				
				TLorentzVector matching_cl3d;
				matching_cl3d.SetPtEtaPhiM( (*_cl3d_pt)[i_cl3d], (*_cl3d_eta)[i_cl3d], (*_cl3d_phi)[i_cl3d], 0);

				float dR_gentau = matching_cl3d.DeltaR(matching_gentau);
				_cl3d_dRMatchedgentau.push_back(dR_gentau);

				float dR_gentauvis = matching_cl3d.DeltaR(matching_gentauvis);
				_cl3d_dRMatchedgentauvis.push_back(dR_gentauvis);

				isMatch = (dR_gentauvis <= dR_max);
				_cl3d_isMatched.push_back(isMatch);

				
				if (isMatch){

					n_match++;

					indicesMatchedcl3d.push_back(i_cl3d);
					dRMatchedcl3d.push_back(dR_gentauvis);
					PtMatchedcl3d.push_back((*_cl3d_pt)[i_cl3d]);
					EtaMatchedcl3d.push_back((*_cl3d_eta)[i_cl3d]);
					PhiMatchedcl3d.push_back((*_cl3d_phi)[i_cl3d]);
					BDTegMatchedcl3d.push_back((*_cl3d_bdteg)[i_cl3d]);
					CoreshowerlengthMatchedcl3d.push_back((*_cl3d_coreshowerlength)[i_cl3d]);
					MaxlayerMatchedcl3d.push_back((*_cl3d_maxlayer)[i_cl3d]);
					
					_cl3d_indexMatchedgentau.push_back(i_gentau);
					_cl3d_dRMatchedgentau.push_back(dR_gentau);
					_cl3d_dRMatchedgentauvis.push_back(dR_gentauvis);

	    				_cl3d_PtMatchedgentau.push_back((*_gentau_pt)[i_gentau]);
					_cl3d_EtaMatchedgentau.push_back((*_gentau_eta)[i_gentau]);
					_cl3d_PhiMatchedgentau.push_back((*_gentau_phi)[i_gentau]);
					_cl3d_EnergyMatchedgentau.push_back((*_gentau_energy)[i_gentau]);
					_cl3d_MassMatchedgentau.push_back((*_gentau_mass)[i_gentau]);
					_cl3d_decayModeMatchedgentau.push_back((*_gentau_decayMode)[i_gentau]);

					_cl3d_PtMatchedgentauvis.push_back((*_gentau_vis_pt)[i_gentau]);
					_cl3d_EtaMatchedgentauvis.push_back((*_gentau_vis_eta)[i_gentau]);
					_cl3d_PhiMatchedgentauvis.push_back((*_gentau_vis_phi)[i_gentau]);
					_cl3d_EnergyMatchedgentauvis.push_back((*_gentau_vis_energy)[i_gentau]);
					_cl3d_MassMatchedgentauvis.push_back((*_gentau_vis_mass)[i_gentau]);


				}

				else if (!isMatch){

					_cl3d_indexMatchedgentau.push_back(-999.);
					_cl3d_dRMatchedgentau.push_back(-999.);
					_cl3d_dRMatchedgentauvis.push_back(-999.);

	    				_cl3d_PtMatchedgentau.push_back(-999.);
					_cl3d_EtaMatchedgentau.push_back(-999.);
					_cl3d_PhiMatchedgentau.push_back(-999.);
					_cl3d_EnergyMatchedgentau.push_back(-999.);
					_cl3d_MassMatchedgentau.push_back(-999.);
					_cl3d_decayModeMatchedgentau.push_back(-999.);

					_cl3d_PtMatchedgentauvis.push_back(-999.);
					_cl3d_EtaMatchedgentauvis.push_back(-999.);
					_cl3d_PhiMatchedgentauvis.push_back(-999.);
					_cl3d_EnergyMatchedgentauvis.push_back(-999.);
					_cl3d_MassMatchedgentauvis.push_back(-999.);

				}

			} //end of loop over cl3d

			_gentau_numberMatchedcl3d.push_back(n_match);
			_gentau_isMatched.push_back(n_match>0);

			_gentau_indicesMatchedcl3d.push_back(indicesMatchedcl3d);
			_gentau_dRMatchedcl3d.push_back(dRMatchedcl3d);
			_gentau_PtMatchedcl3d.push_back(PtMatchedcl3d);
			_gentau_EtaMatchedcl3d.push_back(EtaMatchedcl3d);
			_gentau_PhiMatchedcl3d.push_back(PhiMatchedcl3d);
			_gentau_BDTegMatchedcl3d.push_back(BDTegMatchedcl3d);
			_gentau_CoreshowerlengthMatchedcl3d.push_back(CoreshowerlengthMatchedcl3d);
			_gentau_MaxlayerMatchedcl3d.push_back(MaxlayerMatchedcl3d);

		} //end of loop over gentaus

		out_tree->Fill();

	}

	out_file->cd();

	h_dR_match->Write();

	out_tree->Write();
	out_file->Close();

	return;

}


void test(int n_events = -1, TString pu = "0", TString dr = "0p3"){

  TString dir = "/data_CMS/cms/mperez/HGCal_data/May19/";
  
  float dR;
  if(dr=="0p1") dR = 0.1;
  else if(dr=="0p2") dR = 0.2;
  else if(dr=="0p3") dR = 0.3;
  else if(dr=="0p4") dR = 0.4;
  else if(dr=="0p5") dR = 0.5;
  else if(dr=="0p6") dR = 0.6;

  cout<<"pu"<<pu<<endl;
  cout<<"dR"<<dR<<endl;

  TString infile = dir+"skimmed/NTuple_ZTT_PU"+pu+"_skimmed.root";
  TString outfile = dir+"matched/NTuple_ZTT_PU"+pu+"_matched"+dr+"_cutMaxLayer6.root";

  match_tree(infile, outfile, n_events, dR);

}
