/////////////////////////////////////////////////////////
///// HGCal L1 taus, C. Martin Perez, LLR, May 2019 /////
/////////////////////////////////////////////////////////

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TH1F.h>
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

	//cout<<"here0b"<<endl;

	TChain * in_tree = new TChain("MatchedTree");	
	in_tree->Add(filein);

	Long64_t nentries = in_tree->GetEntries();
	cout<<"nentries="<<in_tree->GetEntries()<<endl;
	if (nevents != -1) 
		nentries = nevents;

	// old branches used

	int _gentau_n;

	vector<bool>   			*_gentau_isMatched;
	vector<vector<float> >	*_gentau_indicesMatchedcl3d;
	vector<int>   			*_gentau_numberMatchedcl3d;

	vector<vector<float> >	*_gentau_dRMatchedcl3d;
	vector<vector<float> >	*_gentau_PtMatchedcl3d;
	vector<vector<float> >	*_gentau_EtaMatchedcl3d;
	vector<vector<float> >	*_gentau_PhiMatchedcl3d;
	vector<vector<float> >	*_gentau_BDTegMatchedcl3d;

	in_tree->SetBranchAddress("gentau_n",	&_gentau_n);

	in_tree->SetBranchAddress("gentau_isMatched",			&_gentau_isMatched);
	in_tree->SetBranchAddress("gentau_indicesMatchedcl3d",	&_gentau_indicesMatchedcl3d);
	in_tree->SetBranchAddress("gentau_numberMatchedcl3d",	&_gentau_numberMatchedcl3d);

	in_tree->SetBranchAddress("gentau_dRMatchedcl3d",		&_gentau_dRMatchedcl3d);
	in_tree->SetBranchAddress("gentau_PtMatchedcl3d",		&_gentau_PtMatchedcl3d);
	in_tree->SetBranchAddress("gentau_EtaMatchedcl3d",		&_gentau_EtaMatchedcl3d);
	in_tree->SetBranchAddress("gentau_PhiMatchedcl3d",		&_gentau_PhiMatchedcl3d);
	in_tree->SetBranchAddress("gentau_BDTegMatchedcl3d",	&_gentau_BDTegMatchedcl3d);

	// new branches

	TTree* out_tree=in_tree->GetTree()->CloneTree(0);
	out_tree->SetNameTitle("SortedTree","SortedTree"); 

	vector<float>   _gentau_TotPtMatchedcl3d;
	vector<float>   _gentau_TotBDTMatchedcl3d;
	vector<float>   _gentau_MaxPtMatchedcl3d;
	vector<float>   _gentau_MaxBDTMatchedcl3d;

	vector<vector<float> >	_gentau_dRMatchedcl3d_pTsorted;
	vector<vector<float> >	_gentau_PtMatchedcl3d_pTsorted;
	vector<vector<float> >	_gentau_EtaMatchedcl3d_pTsorted;
	vector<vector<float> >	_gentau_PhiMatchedcl3d_pTsorted;
	vector<vector<float> >	_gentau_BDTegMatchedcl3d_pTsorted;

	vector<vector<float> >	_gentau_dRMatchedcl3d_BDTsorted;
	vector<vector<float> >	_gentau_PtMatchedcl3d_BDTsorted;
	vector<vector<float> >	_gentau_EtaMatchedcl3d_BDTsorted;
	vector<vector<float> >	_gentau_PhiMatchedcl3d_BDTsorted;
	vector<vector<float> >	_gentau_BDTegMatchedcl3d_BDTsorted;

	out_tree->Branch("gentau_TotPtMatchedcl3d",&_gentau_TotPtMatchedcl3d);
	out_tree->Branch("gentau_TotBDTMatchedcl3d",&_gentau_TotBDTMatchedcl3d);
	out_tree->Branch("gentau_MaxPtMatchedcl3d",&_gentau_MaxPtMatchedcl3d);
	out_tree->Branch("gentau_MaxBDTMatchedcl3d",&_gentau_MaxBDTMatchedcl3d);

	out_tree->Branch("gentau_dRMatchedcl3d_pTsorted",		&_gentau_dRMatchedcl3d_pTsorted);
	out_tree->Branch("gentau_PtMatchedcl3d_pTsorted",		&_gentau_PtMatchedcl3d_pTsorted);
	out_tree->Branch("gentau_EtaMatchedcl3d_pTsorted",		&_gentau_EtaMatchedcl3d_pTsorted);
	out_tree->Branch("gentau_PhiMatchedcl3d_pTsorted",		&_gentau_PhiMatchedcl3d_pTsorted);
	out_tree->Branch("gentau_BDTegMatchedcl3d_pTsorted",	&_gentau_BDTegMatchedcl3d_pTsorted);

	out_tree->Branch("gentau_dRMatchedcl3d_BDTsorted",		&_gentau_dRMatchedcl3d_BDTsorted);
	out_tree->Branch("gentau_PtMatchedcl3d_BDTsorted",		&_gentau_PtMatchedcl3d_BDTsorted);
	out_tree->Branch("gentau_EtaMatchedcl3d_BDTsorted",		&_gentau_EtaMatchedcl3d_BDTsorted);
	out_tree->Branch("gentau_PhiMatchedcl3d_BDTsorted",		&_gentau_PhiMatchedcl3d_BDTsorted);
	out_tree->Branch("gentau_BDTegMatchedcl3d_BDTsorted",	&_gentau_BDTegMatchedcl3d_BDTsorted);


	// variables for BDT

	int _numberMatchedcl3d;

	float _pT1cl3d_pTfraction;
	float _pT2cl3d_pTfraction;
	float _pT3cl3d_pTfraction;

	float _pT1cl3d_BDTeg;
	float _pT2cl3d_BDTeg;
	float _pT3cl3d_BDTeg;

	float _BDTeg1cl3d_pTfraction;
	float _BDTeg2cl3d_pTfraction;
	float _BDTeg3cl3d_pTfraction;

	float _BDTeg1cl3d_BDTeg;
	float _BDTeg2cl3d_BDTeg;
	float _BDTeg3cl3d_BDTeg;

	out_tree->Branch("numberMatchedcl3d",&_numberMatchedcl3d);

	out_tree->Branch("pT1cl3d_pTfraction",&_pT1cl3d_pTfraction);
	out_tree->Branch("pT2cl3d_pTfraction",&_pT2cl3d_pTfraction);
	out_tree->Branch("pT3cl3d_pTfraction",&_pT3cl3d_pTfraction);

	out_tree->Branch("pT1cl3d_BDTeg",&_pT1cl3d_BDTeg);
	out_tree->Branch("pT2cl3d_BDTeg",&_pT2cl3d_BDTeg);
	out_tree->Branch("pT3cl3d_BDTeg",&_pT3cl3d_BDTeg);

	out_tree->Branch("BDTeg1cl3d_pTfraction",&_BDTeg1cl3d_pTfraction);
	out_tree->Branch("BDTeg2cl3d_pTfraction",&_BDTeg2cl3d_pTfraction);
	out_tree->Branch("BDTeg3cl3d_pTfraction",&_BDTeg3cl3d_pTfraction);

	out_tree->Branch("BDTeg1cl3d_BDTeg",&_BDTeg1cl3d_BDTeg);
	out_tree->Branch("BDTeg1cl3d_BDTeg",&_BDTeg1cl3d_BDTeg);
	out_tree->Branch("BDTeg1cl3d_BDTeg",&_BDTeg1cl3d_BDTeg);


	for (int i=0;i<nentries;i++) {

		if(i%1000==0) cout<<"i="<<i<<endl;

		//old branches

		_gentau_n = 0;

		_gentau_isMatched = 0;
		_gentau_indicesMatchedcl3d = 0;
		_gentau_numberMatchedcl3d = 0;

		_gentau_dRMatchedcl3d = 0;
		_gentau_PtMatchedcl3d = 0;
		_gentau_EtaMatchedcl3d = 0;
		_gentau_PhiMatchedcl3d = 0;
		_gentau_BDTegMatchedcl3d = 0;

		//new branches

		_gentau_TotPtMatchedcl3d.clear();
		_gentau_TotBDTMatchedcl3d.clear();
		_gentau_MaxPtMatchedcl3d.clear();
		_gentau_MaxBDTMatchedcl3d.clear();

		_gentau_dRMatchedcl3d_pTsorted.clear();
		_gentau_PtMatchedcl3d_pTsorted.clear();
		_gentau_EtaMatchedcl3d_pTsorted.clear();
		_gentau_PhiMatchedcl3d_pTsorted.clear();
		_gentau_BDTegMatchedcl3d_pTsorted.clear();

		_gentau_dRMatchedcl3d_BDTsorted.clear();
		_gentau_PtMatchedcl3d_BDTsorted.clear();
		_gentau_EtaMatchedcl3d_BDTsorted.clear();
		_gentau_PhiMatchedcl3d_BDTsorted.clear();
		_gentau_BDTegMatchedcl3d_BDTsorted.clear();

		_numberMatchedcl3d = 0;

		_pT1cl3d_pTfraction = 0;
		_pT2cl3d_pTfraction = 0;
		_pT3cl3d_pTfraction = 0;

		_pT1cl3d_BDTeg = 0;
		_pT2cl3d_BDTeg = 0;
		_pT3cl3d_BDTeg = 0;

		_BDTeg1cl3d_pTfraction = 0;
		_BDTeg2cl3d_pTfraction = 0;
		_BDTeg3cl3d_pTfraction = 0;

		_BDTeg1cl3d_BDTeg = 0;
		_BDTeg2cl3d_BDTeg = 0;
		_BDTeg3cl3d_BDTeg = 0;


		int entry_ok = in_tree->GetEntry(i);	
		if(entry_ok<0) 
			continue;


		for(int i_gentau=0; i_gentau<_gentau_n; i_gentau++){

			if (!(*_gentau_isMatched)[i_gentau]) continue;

			vector<float> dRMatchedcl3d = (*_gentau_dRMatchedcl3d)[i_gentau];
			vector<float> PtMatchedcl3d = (*_gentau_PtMatchedcl3d)[i_gentau];
			vector<float> EtaMatchedcl3d = (*_gentau_EtaMatchedcl3d)[i_gentau];
			vector<float> PhiMatchedcl3d = (*_gentau_PhiMatchedcl3d)[i_gentau];
			vector<float> BDTegMatchedcl3d = (*_gentau_BDTegMatchedcl3d)[i_gentau];

			/////////////////////////////////
			///////// pT sorting ////////////
			/////////////////////////////////

			vector< pair<int,TLorentzVector> > icl3d_cl3dtlv_pairs;
			TLorentzVector matchedcl3d;

			for (unsigned int i_cl3d = 0; i_cl3d<PtMatchedcl3d.size(); i_cl3d++){

				matchedcl3d.SetPtEtaPhiM( PtMatchedcl3d[i_cl3d], EtaMatchedcl3d[i_cl3d], PhiMatchedcl3d[i_cl3d], 0 );
				pair<int,TLorentzVector> matchedcl3d_pair = make_pair(i_cl3d,matchedcl3d);
				icl3d_cl3dtlv_pairs.push_back(matchedcl3d_pair);

			}

			sort(icl3d_cl3dtlv_pairs.begin(), icl3d_cl3dtlv_pairs.end(), pT_comparison_pairs);

			vector<float> dRMatchedcl3d_pTsorted;
			vector<float> PtMatchedcl3d_pTsorted;
			vector<float> EtaMatchedcl3d_pTsorted;
			vector<float> PhiMatchedcl3d_pTsorted;
			vector<float> BDTegMatchedcl3d_pTsorted;

			dRMatchedcl3d_pTsorted.clear();
			PtMatchedcl3d_pTsorted.clear();
			EtaMatchedcl3d_pTsorted.clear();
			PhiMatchedcl3d_pTsorted.clear();
			BDTegMatchedcl3d_pTsorted.clear();

			for (unsigned int i_cl3d = 0; i_cl3d<PtMatchedcl3d.size(); i_cl3d++){

				int i_mat_cl3d = icl3d_cl3dtlv_pairs[i_cl3d].first;
				TLorentzVector mat_cl3d = icl3d_cl3dtlv_pairs[i_cl3d].second;

				dRMatchedcl3d_pTsorted.push_back(dRMatchedcl3d[i_mat_cl3d]);
				PtMatchedcl3d_pTsorted.push_back(mat_cl3d.Pt());
				EtaMatchedcl3d_pTsorted.push_back(mat_cl3d.Eta());
				PhiMatchedcl3d_pTsorted.push_back(mat_cl3d.Phi());
				BDTegMatchedcl3d_pTsorted.push_back(BDTegMatchedcl3d[i_mat_cl3d]);

			}

			_gentau_dRMatchedcl3d_pTsorted.push_back(dRMatchedcl3d_pTsorted);
			_gentau_PtMatchedcl3d_pTsorted.push_back(PtMatchedcl3d_pTsorted);
			_gentau_EtaMatchedcl3d_pTsorted.push_back(EtaMatchedcl3d_pTsorted);
			_gentau_PhiMatchedcl3d_pTsorted.push_back(PhiMatchedcl3d_pTsorted);
			_gentau_BDTegMatchedcl3d_pTsorted.push_back(BDTegMatchedcl3d_pTsorted);


			/////////////////////////////////
			/////// BDT EG sorting //////////
			/////////////////////////////////

			vector< pair<int,float> > icl3d_BDT_pairs;

			for(unsigned int i_cl3d=0; i_cl3d<PtMatchedcl3d.size(); i_cl3d++){

				float BDT = BDTegMatchedcl3d[i_cl3d];
				pair<int,float> BDT_pair = make_pair(i_cl3d,BDT);
				icl3d_BDT_pairs.push_back(BDT_pair);

			}

			sort(icl3d_BDT_pairs.begin(), icl3d_BDT_pairs.end(), BDT_comparison_pairs);

			vector<float> dRMatchedcl3d_BDTsorted;
			vector<float> PtMatchedcl3d_BDTsorted;
			vector<float> EtaMatchedcl3d_BDTsorted;
			vector<float> PhiMatchedcl3d_BDTsorted;
			vector<float> BDTegMatchedcl3d_BDTsorted;

			dRMatchedcl3d_BDTsorted.clear();
			PtMatchedcl3d_BDTsorted.clear();
			EtaMatchedcl3d_BDTsorted.clear();
			PhiMatchedcl3d_BDTsorted.clear();
			BDTegMatchedcl3d_BDTsorted.clear();

			for (unsigned int i_cl3d = 0; i_cl3d<PtMatchedcl3d.size(); i_cl3d++){

				int i_mat_cl3d = icl3d_BDT_pairs[i_cl3d].first;
				float bdt_cl3d = icl3d_BDT_pairs[i_cl3d].second;

				dRMatchedcl3d_BDTsorted.push_back(dRMatchedcl3d[i_mat_cl3d]);
				PtMatchedcl3d_BDTsorted.push_back(PtMatchedcl3d[i_mat_cl3d]);
				EtaMatchedcl3d_BDTsorted.push_back(EtaMatchedcl3d[i_mat_cl3d]);
				PhiMatchedcl3d_BDTsorted.push_back(PhiMatchedcl3d[i_mat_cl3d]);
				BDTegMatchedcl3d_BDTsorted.push_back(bdt_cl3d);

			}

			_gentau_dRMatchedcl3d_BDTsorted.push_back(dRMatchedcl3d_BDTsorted);
			_gentau_PtMatchedcl3d_BDTsorted.push_back(PtMatchedcl3d_BDTsorted);
			_gentau_EtaMatchedcl3d_BDTsorted.push_back(EtaMatchedcl3d_BDTsorted);
			_gentau_PhiMatchedcl3d_BDTsorted.push_back(PhiMatchedcl3d_BDTsorted);
			_gentau_BDTegMatchedcl3d_BDTsorted.push_back(BDTegMatchedcl3d_BDTsorted);


			/*cout<<" --- "<<endl;
			cout<<" unsorted, size: "<<PtMatchedcl3d.size()<<endl;
			for (unsigned int i_cl3d = 0; i_cl3d<PtMatchedcl3d.size(); i_cl3d++){	
				cout<<"  i_cl3d: pT "<<PtMatchedcl3d[i_cl3d]<<", eta "<<EtaMatchedcl3d[i_cl3d]<<",phi "<<PhiMatchedcl3d[i_cl3d]<<",bdt "<<BDTegMatchedcl3d[i_cl3d]<<endl;
			}*/

			cout<<" sorted by pT, size: "<<PtMatchedcl3d_pTsorted.size()<<endl;
			for (unsigned int i_cl3d = 0; i_cl3d<PtMatchedcl3d_pTsorted.size(); i_cl3d++){				
				cout<<"  i_cl3d: pT "<<PtMatchedcl3d_pTsorted[i_cl3d]<<", eta "<<EtaMatchedcl3d_pTsorted[i_cl3d]<<",phi "<<PhiMatchedcl3d_pTsorted[i_cl3d]<<",bdt "<<BDTegMatchedcl3d_pTsorted[i_cl3d]<<endl;
			}

			/*cout<<" sorted by BDT, size: "<<PtMatchedcl3d_BDTsorted.size()<<endl;
			for (unsigned int i_cl3d = 0; i_cl3d<PtMatchedcl3d_BDTsorted.size(); i_cl3d++){				
				cout<<"  i_cl3d: pT "<<PtMatchedcl3d_BDTsorted[i_cl3d]<<", eta "<<EtaMatchedcl3d_BDTsorted[i_cl3d]<<",phi "<<PhiMatchedcl3d_BDTsorted[i_cl3d]<<",bdt "<<BDTegMatchedcl3d_BDTsorted[i_cl3d]<<endl;
			}*/

			float TotPt = 0.;
			float TotBDT = 0.;

			for (unsigned int i_cl3d = 0; i_cl3d<PtMatchedcl3d.size(); i_cl3d++){
				TotPt += PtMatchedcl3d[i_cl3d];
				TotBDT += BDTegMatchedcl3d[i_cl3d];
			}

			float MaxPt = PtMatchedcl3d_pTsorted[0];
			float MaxBDT = BDTegMatchedcl3d_BDTsorted[0];

			//cout<<"MaxPt "<<MaxPt<<",MaxBDT "<<MaxBDT<<",TotPt "<<TotPt<<",TotBDT "<<TotBDT<<endl;

			_gentau_TotPtMatchedcl3d.push_back(TotPt);
			_gentau_TotBDTMatchedcl3d.push_back(TotBDT);
			_gentau_MaxPtMatchedcl3d.push_back(MaxPt);
			_gentau_MaxBDTMatchedcl3d.push_back(MaxBDT);


			// Variables for BDT

			_numberMatchedcl3d = (*_gentau_numberMatchedcl3d)[i_gentau];

			_pT1cl3d_pTfraction = -1.;
			_pT2cl3d_pTfraction = -1.;
			_pT3cl3d_pTfraction = -1.;

			_pT1cl3d_BDTeg = -1.;
			_pT2cl3d_BDTeg = -1.;
			_pT3cl3d_BDTeg = -1.;

			_BDTeg1cl3d_pTfraction = -1;
			_BDTeg2cl3d_pTfraction = -1;
			_BDTeg3cl3d_pTfraction = -1;

			_BDTeg1cl3d_BDTeg = -1;
			_BDTeg2cl3d_BDTeg = -1;
			_BDTeg3cl3d_BDTeg = -1;

			if(_numberMatchedcl3d>=1){

				cout<<"1 matched"<<endl;
				cout<<"PtMatchedcl3d_pTsorted[0] "<<PtMatchedcl3d_pTsorted[0]<<endl;

				_pT1cl3d_pTfraction = PtMatchedcl3d_pTsorted[0]/TotPt;
				_pT1cl3d_BDTeg = BDTegMatchedcl3d_pTsorted[0];
				_BDTeg1cl3d_pTfraction = PtMatchedcl3d_BDTsorted[0]/TotPt;
				_BDTeg1cl3d_BDTeg = BDTegMatchedcl3d_BDTsorted[0];

			}

			if(_numberMatchedcl3d>=2){

				cout<<"2 matched"<<endl;
				cout<<"PtMatchedcl3d_pTsorted[1] "<<PtMatchedcl3d_pTsorted[1]<<endl;

				_pT2cl3d_pTfraction = PtMatchedcl3d_pTsorted[1]/TotPt;
				_pT2cl3d_BDTeg = BDTegMatchedcl3d_pTsorted[1];
				_BDTeg2cl3d_pTfraction = PtMatchedcl3d_BDTsorted[1]/TotPt;
				_BDTeg2cl3d_BDTeg = BDTegMatchedcl3d_BDTsorted[1];

			}

			if(_numberMatchedcl3d>=3){

				cout<<"3 matched"<<endl;
				cout<<"PtMatchedcl3d_pTsorted[2] "<<PtMatchedcl3d_pTsorted[2]<<endl;

				_pT3cl3d_pTfraction = PtMatchedcl3d_pTsorted[2]/TotPt;
				_pT3cl3d_BDTeg = BDTegMatchedcl3d_pTsorted[2];
				_BDTeg3cl3d_pTfraction = PtMatchedcl3d_BDTsorted[2]/TotPt;
				_BDTeg3cl3d_BDTeg = BDTegMatchedcl3d_BDTsorted[2];

			}

		}

		out_tree->Fill();

	}


	out_file->cd();

	out_tree->Write();
	out_file->Close();

	return;

}


void test(int n_events = -1){

  TString dir = "/data_CMS/cms/mperez/HGCal_data/May19/";

  TString infile = dir+"NTuple_ZTT_PU0_matched0p3.root";
  TString outfile = dir+"NTuple_ZTT_PU0_sorted.root";

  sort_tree(infile, outfile, n_events);

}

