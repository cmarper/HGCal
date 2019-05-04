#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <iostream>

using namespace std;


void skim_tree( vector<TString> filein, TString fileout, int nevents ){
	
	TFile* out_file = TFile::Open(fileout);
  	if(out_file!=0){
      cout<<fileout<<" already exists, please delete it before converting again"<<endl;
      return;
	}

	out_file = TFile::Open(fileout,"RECREATE");

	TChain * in_tree = new TChain("hgcalTriggerNtuplizer/HGCalTriggerNtuple");

	for (unsigned int ifile = 0; ifile<filein.size(); ifile++){
		in_tree->Add(filein[ifile]);
	}

	Long64_t nentries = in_tree->GetEntries();
	cout<<"nentries="<<in_tree->GetEntries()<<endl;
	if (nevents != -1) nentries = nevents;

	// old branches used

	int _in_run;
	int _in_event;
	int _in_lumi;

	vector<float> *_in_gentau_pt; 
	vector<float> *_in_gentau_eta;
	vector<float> *_in_gentau_phi;
	vector<float> *_in_gentau_energy;
	vector<float> *_in_gentau_mass;

	vector<float> *_in_gentau_vis_pt; 
	vector<float> *_in_gentau_vis_eta;
	vector<float> *_in_gentau_vis_phi;
	vector<float> *_in_gentau_vis_energy;
	vector<float> *_in_gentau_vis_mass;

	vector<vector<float> > *_in_gentau_products_pt;
	vector<vector<float> > *_in_gentau_products_eta;
	vector<vector<float> > *_in_gentau_products_phi;
	vector<vector<float> > *_in_gentau_products_energy;
	vector<vector<float> > *_in_gentau_products_mass;
	vector<vector<float> > *_in_gentau_products_id;

	vector<int>   *_in_gentau_decayMode;
	vector<int>   *_in_gentau_totNproducts;
	vector<int>   *_in_gentau_totNgamma;
	vector<int>   *_in_gentau_totNpiZero;
	vector<int>   *_in_gentau_totNcharged;

	int _in_cl3d_n;

	vector<unsigned int> *_in_cl3d_id;
	vector<float> 		 *_in_cl3d_pt;
	vector<float> 		 *_in_cl3d_energy;
	vector<float> 		 *_in_cl3d_eta;
	vector<float> 		 *_in_cl3d_phi;

	vector<int>   					*_in_cl3d_clusters_n;
	vector<vector<unsigned int> >   *_in_cl3d_clusters_id;

	vector<int>   *_in_cl3d_showerlength;
	vector<int>   *_in_cl3d_coreshowerlength;
	vector<int>   *_in_cl3d_firstlayer;
	vector<int>   *_in_cl3d_maxlayer;

	vector<float> *_in_cl3d_seetot;
	vector<float> *_in_cl3d_seemax;
	vector<float> *_in_cl3d_spptot;
	vector<float> *_in_cl3d_sppmax;
	vector<float> *_in_cl3d_szz;
	vector<float> *_in_cl3d_srrtot;
	vector<float> *_in_cl3d_srrmax;
	vector<float> *_in_cl3d_srrmean;
	vector<float> *_in_cl3d_emaxe;
	vector<float> *_in_cl3d_bdteg;
	vector<int>   *_in_cl3d_quality;

	in_tree->SetBranchAddress("run",	&_in_run);
	in_tree->SetBranchAddress("event",	&_in_event);
	in_tree->SetBranchAddress("lumi",	&_in_lumi);

	in_tree->SetBranchAddress("gentau_pt",		&_in_gentau_pt);
	in_tree->SetBranchAddress("gentau_eta",		&_in_gentau_eta);
	in_tree->SetBranchAddress("gentau_phi",		&_in_gentau_phi);
	in_tree->SetBranchAddress("gentau_energy",	&_in_gentau_energy);
	in_tree->SetBranchAddress("gentau_mass",	&_in_gentau_mass);

	in_tree->SetBranchAddress("gentau_vis_pt",		&_in_gentau_vis_pt);
	in_tree->SetBranchAddress("gentau_vis_eta",		&_in_gentau_vis_eta);
	in_tree->SetBranchAddress("gentau_vis_phi",		&_in_gentau_vis_phi);
	in_tree->SetBranchAddress("gentau_vis_energy",	&_in_gentau_vis_energy);
	in_tree->SetBranchAddress("gentau_vis_mass",	&_in_gentau_vis_mass);

	in_tree->SetBranchAddress("gentau_products_pt",		&_in_gentau_products_pt);
	in_tree->SetBranchAddress("gentau_products_eta",	&_in_gentau_products_eta);
	in_tree->SetBranchAddress("gentau_products_phi",	&_in_gentau_products_phi);
	in_tree->SetBranchAddress("gentau_products_energy",	&_in_gentau_products_energy);
	in_tree->SetBranchAddress("gentau_products_mass",	&_in_gentau_products_mass);
	in_tree->SetBranchAddress("gentau_products_id",		&_in_gentau_products_id);

	in_tree->SetBranchAddress("gentau_decayMode",		&_in_gentau_decayMode);
	in_tree->SetBranchAddress("gentau_totNproducts",	&_in_gentau_totNproducts);
	in_tree->SetBranchAddress("gentau_totNgamma",		&_in_gentau_totNgamma);
	in_tree->SetBranchAddress("gentau_totNpiZero",		&_in_gentau_totNpiZero);
	in_tree->SetBranchAddress("gentau_totNcharged",		&_in_gentau_totNcharged);

	in_tree->SetBranchAddress("cl3d_n",&_in_cl3d_n);

	in_tree->SetBranchAddress("cl3d_id",		&_in_cl3d_id);
	in_tree->SetBranchAddress("cl3d_pt",		&_in_cl3d_pt);
	in_tree->SetBranchAddress("cl3d_energy",	&_in_cl3d_energy);
	in_tree->SetBranchAddress("cl3d_eta",		&_in_cl3d_eta);
	in_tree->SetBranchAddress("cl3d_phi",		&_in_cl3d_phi);

	in_tree->SetBranchAddress("cl3d_clusters_n",	&_in_cl3d_clusters_n);
	in_tree->SetBranchAddress("cl3d_clusters_id",	&_in_cl3d_clusters_id);

	in_tree->SetBranchAddress("cl3d_showerlength",		&_in_cl3d_showerlength);
	in_tree->SetBranchAddress("cl3d_coreshowerlength",	&_in_cl3d_coreshowerlength);
	in_tree->SetBranchAddress("cl3d_firstlayer",		&_in_cl3d_firstlayer);
	in_tree->SetBranchAddress("cl3d_maxlayer",			&_in_cl3d_maxlayer);

	in_tree->SetBranchAddress("cl3d_seetot",	&_in_cl3d_seetot);
	in_tree->SetBranchAddress("cl3d_seemax",	&_in_cl3d_seemax);
	in_tree->SetBranchAddress("cl3d_spptot",	&_in_cl3d_spptot);
	in_tree->SetBranchAddress("cl3d_sppmax",	&_in_cl3d_sppmax);
	in_tree->SetBranchAddress("cl3d_szz",		&_in_cl3d_szz);
	in_tree->SetBranchAddress("cl3d_srrtot",	&_in_cl3d_srrtot);
	in_tree->SetBranchAddress("cl3d_srrmax",	&_in_cl3d_srrmax);
	in_tree->SetBranchAddress("cl3d_srrmean",	&_in_cl3d_srrmean);
	in_tree->SetBranchAddress("cl3d_emaxe",		&_in_cl3d_emaxe);
	in_tree->SetBranchAddress("cl3d_bdteg",		&_in_cl3d_bdteg);
	in_tree->SetBranchAddress("cl3d_quality",	&_in_cl3d_quality);


	TTree* out_tree = new TTree("SkimmedTree","SkimmedTree");

	int _out_run;
	int _out_event;
	int _out_lumi;

	int _out_gentau_n;

	vector<float> _out_gentau_pt; 
	vector<float> _out_gentau_eta;
	vector<float> _out_gentau_phi;
	vector<float> _out_gentau_energy;
	vector<float> _out_gentau_mass;

	vector<float> _out_gentau_vis_pt; 
	vector<float> _out_gentau_vis_eta;
	vector<float> _out_gentau_vis_phi;
	vector<float> _out_gentau_vis_energy;
	vector<float> _out_gentau_vis_mass;

	vector<vector<float> > _out_gentau_products_pt;
	vector<vector<float> > _out_gentau_products_eta;
	vector<vector<float> > _out_gentau_products_phi;
	vector<vector<float> > _out_gentau_products_energy;
	vector<vector<float> > _out_gentau_products_mass;
	vector<vector<float> > _out_gentau_products_id;

	vector<int>   _out_gentau_decayMode;
	vector<int>   _out_gentau_totNproducts;
	vector<int>   _out_gentau_totNgamma;
	vector<int>   _out_gentau_totNpiZero;
	vector<int>   _out_gentau_totNcharged;

	int _out_cl3d_n;

	vector<unsigned int> _out_cl3d_id;
	vector<float> 		 _out_cl3d_pt;
	vector<float> 		 _out_cl3d_energy;
	vector<float> 		 _out_cl3d_eta;
	vector<float> 		 _out_cl3d_phi;

	vector<int>   					_out_cl3d_clusters_n;
	vector<vector<unsigned int> >   _out_cl3d_clusters_id;

	vector<int>   _out_cl3d_showerlength;
	vector<int>   _out_cl3d_coreshowerlength;
	vector<int>   _out_cl3d_firstlayer;
	vector<int>   _out_cl3d_maxlayer;

	vector<float> _out_cl3d_seetot;
	vector<float> _out_cl3d_seemax;
	vector<float> _out_cl3d_spptot;
	vector<float> _out_cl3d_sppmax;
	vector<float> _out_cl3d_szz;
	vector<float> _out_cl3d_srrtot;
	vector<float> _out_cl3d_srrmax;
	vector<float> _out_cl3d_srrmean;
	vector<float> _out_cl3d_emaxe;
	vector<float> _out_cl3d_bdteg;
	vector<int>   _out_cl3d_quality;


	out_tree->Branch("run",		&_in_run);
	out_tree->Branch("event",	&_in_event);
	out_tree->Branch("lumi",	&_in_lumi);

	out_tree->Branch("gentau_n",&_out_gentau_n);

	out_tree->Branch("gentau_pt",		&_out_gentau_pt);
	out_tree->Branch("gentau_eta",		&_out_gentau_eta);
	out_tree->Branch("gentau_phi",		&_out_gentau_phi);
	out_tree->Branch("gentau_energy",	&_out_gentau_energy);
	out_tree->Branch("gentau_mass",		&_out_gentau_mass);

	out_tree->Branch("gentau_vis_pt",		&_out_gentau_vis_pt);
	out_tree->Branch("gentau_vis_eta",		&_out_gentau_vis_eta);
	out_tree->Branch("gentau_vis_phi",		&_out_gentau_vis_phi);
	out_tree->Branch("gentau_vis_energy",	&_out_gentau_vis_energy);
	out_tree->Branch("gentau_vis_mass",		&_out_gentau_vis_mass);

	out_tree->Branch("gentau_products_pt",		&_out_gentau_products_pt);
	out_tree->Branch("gentau_products_eta",		&_out_gentau_products_eta);
	out_tree->Branch("gentau_products_phi",		&_out_gentau_products_phi);
	out_tree->Branch("gentau_products_energy",	&_out_gentau_products_energy);
	out_tree->Branch("gentau_products_mass",	&_out_gentau_products_mass);
	out_tree->Branch("gentau_products_id",		&_out_gentau_products_id);

	out_tree->Branch("gentau_decayMode",		&_out_gentau_decayMode);
	out_tree->Branch("gentau_totNproducts",		&_out_gentau_totNproducts);
	out_tree->Branch("gentau_totNgamma",		&_out_gentau_totNgamma);
	out_tree->Branch("gentau_totNpiZero",		&_out_gentau_totNpiZero);
	out_tree->Branch("gentau_totNcharged",		&_out_gentau_totNcharged);

	out_tree->Branch("cl3d_n",&_out_cl3d_n);

	out_tree->Branch("cl3d_id",		&_out_cl3d_id);
	out_tree->Branch("cl3d_pt",		&_out_cl3d_pt);
	out_tree->Branch("cl3d_energy",	&_out_cl3d_energy);
	out_tree->Branch("cl3d_eta",		&_out_cl3d_eta);
	out_tree->Branch("cl3d_phi",		&_out_cl3d_phi);

	out_tree->Branch("cl3d_clusters_n",	&_out_cl3d_clusters_n);
	out_tree->Branch("cl3d_clusters_id",	&_out_cl3d_clusters_id);

	out_tree->Branch("cl3d_showerlength",		&_out_cl3d_showerlength);
	out_tree->Branch("cl3d_coreshowerlength",	&_out_cl3d_coreshowerlength);
	out_tree->Branch("cl3d_firstlayer",			&_out_cl3d_firstlayer);
	out_tree->Branch("cl3d_maxlayer",			&_out_cl3d_maxlayer);

	out_tree->Branch("cl3d_seetot",	&_out_cl3d_seetot);
	out_tree->Branch("cl3d_seemax",	&_out_cl3d_seemax);
	out_tree->Branch("cl3d_spptot",	&_out_cl3d_spptot);
	out_tree->Branch("cl3d_sppmax",	&_out_cl3d_sppmax);
	out_tree->Branch("cl3d_szz",		&_out_cl3d_szz);
	out_tree->Branch("cl3d_srrtot",	&_out_cl3d_srrtot);
	out_tree->Branch("cl3d_srrmax",	&_out_cl3d_srrmax);
	out_tree->Branch("cl3d_srrmean",	&_out_cl3d_srrmean);
	out_tree->Branch("cl3d_emaxe",	&_out_cl3d_emaxe);
	out_tree->Branch("cl3d_bdteg",	&_out_cl3d_bdteg);
	out_tree->Branch("cl3d_quality",	&_out_cl3d_quality);


	for (int i=0;i<nentries;i++) {

		if(i%1000==0) cout<<"i="<<i<<endl;

		//old branches

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
		
		_in_gentau_products_pt = 0;
		_in_gentau_products_eta = 0;
		_in_gentau_products_phi = 0;
		_in_gentau_products_energy = 0;
		_in_gentau_products_mass = 0;
		_in_gentau_products_id = 0;
		
		_in_gentau_decayMode = 0;
		_in_gentau_totNproducts = 0;
		_in_gentau_totNgamma = 0;
		_in_gentau_totNpiZero = 0;
		_in_gentau_totNcharged = 0;

		_in_cl3d_n = 0;
		
		_in_cl3d_id = 0;
		_in_cl3d_pt = 0;
		_in_cl3d_energy = 0;
		_in_cl3d_eta = 0;
		_in_cl3d_phi = 0;
		
		_in_cl3d_clusters_n = 0;
		_in_cl3d_clusters_id = 0;
		
		_in_cl3d_showerlength = 0;
		_in_cl3d_coreshowerlength = 0;
		_in_cl3d_firstlayer = 0;
		_in_cl3d_maxlayer = 0;
		_in_cl3d_quality = 0;
		
		_in_cl3d_seetot = 0;
		_in_cl3d_seemax = 0;
		_in_cl3d_spptot = 0;
		_in_cl3d_sppmax = 0;
		_in_cl3d_szz = 0;
		_in_cl3d_srrtot = 0;
		_in_cl3d_srrmax = 0;
		_in_cl3d_srrmean = 0;
		_in_cl3d_emaxe = 0;
		_in_cl3d_bdteg = 0;
		_in_cl3d_quality = 0;


		//new branches

		_out_gentau_n = 0;

		_out_gentau_pt.clear(); 
		_out_gentau_eta.clear();
		_out_gentau_phi.clear();
		_out_gentau_energy.clear();
		_out_gentau_mass.clear();
		
		_out_gentau_vis_pt.clear(); 
		_out_gentau_vis_eta.clear();
		_out_gentau_vis_phi.clear();
		_out_gentau_vis_energy.clear();
		_out_gentau_vis_mass.clear();
		
		_out_gentau_products_pt.clear();
		_out_gentau_products_eta.clear();
		_out_gentau_products_phi.clear();
		_out_gentau_products_energy.clear();
		_out_gentau_products_mass.clear();
		_out_gentau_products_id.clear();
		
		_out_gentau_decayMode.clear();
		_out_gentau_totNproducts.clear();
		_out_gentau_totNgamma.clear();
		_out_gentau_totNpiZero.clear();
		_out_gentau_totNcharged.clear();


		_out_cl3d_n = 0;
		
		_out_cl3d_id.clear();
		_out_cl3d_pt.clear();
		_out_cl3d_energy.clear();
		_out_cl3d_eta.clear();
		_out_cl3d_phi.clear();
		
		_out_cl3d_clusters_n.clear();
		_out_cl3d_clusters_id.clear();
		
		_out_cl3d_showerlength.clear();
		_out_cl3d_coreshowerlength.clear();
		_out_cl3d_firstlayer.clear();
		_out_cl3d_maxlayer.clear();
		_out_cl3d_quality.clear();
		
		_out_cl3d_seetot.clear();
		_out_cl3d_seemax.clear();
		_out_cl3d_spptot.clear();
		_out_cl3d_sppmax.clear();
		_out_cl3d_szz.clear();
		_out_cl3d_srrtot.clear();
		_out_cl3d_srrmax.clear();
		_out_cl3d_srrmean.clear();
		_out_cl3d_emaxe.clear();
		_out_cl3d_bdteg.clear();
		_out_cl3d_quality.clear();


		//loop through entries

		int entry_ok = in_tree->GetEntry(i);
		if(entry_ok<0) continue;

		_out_run   = _in_run;
		_out_event = _in_event;
		_out_lumi  = _in_lumi;

		// loop over gentaus to select the ones in the endcaps

		int n_gentaus = (*_in_gentau_pt).size();

		for (int i_gentau=0; i_gentau<n_gentaus; i_gentau++){

			if ( abs( (*_in_gentau_eta)[i_gentau] ) <= 1.7 || abs( (*_in_gentau_eta)[i_gentau] ) >= 2.8 ) continue;

			_out_gentau_pt.push_back((*_in_gentau_pt)[i_gentau]);
			_out_gentau_eta.push_back((*_in_gentau_eta)[i_gentau]);
			_out_gentau_phi.push_back((*_in_gentau_phi)[i_gentau]);
			_out_gentau_energy.push_back((*_in_gentau_energy)[i_gentau]);
			_out_gentau_mass.push_back((*_in_gentau_mass)[i_gentau]);
		
			_out_gentau_vis_pt.push_back((*_in_gentau_vis_pt)[i_gentau]);
			_out_gentau_vis_eta.push_back((*_in_gentau_vis_eta)[i_gentau]);
			_out_gentau_vis_phi.push_back((*_in_gentau_vis_phi)[i_gentau]);
			_out_gentau_vis_energy.push_back((*_in_gentau_vis_energy)[i_gentau]);
			_out_gentau_vis_mass.push_back((*_in_gentau_vis_mass)[i_gentau]);
		
			_out_gentau_products_pt.push_back((*_in_gentau_products_pt)[i_gentau]);
			_out_gentau_products_eta.push_back((*_in_gentau_products_eta)[i_gentau]);
			_out_gentau_products_phi.push_back((*_in_gentau_products_phi)[i_gentau]);
			_out_gentau_products_energy.push_back((*_in_gentau_products_energy)[i_gentau]);
			_out_gentau_products_mass.push_back((*_in_gentau_products_mass)[i_gentau]);
			_out_gentau_products_id.push_back((*_in_gentau_products_id)[i_gentau]);
		
			_out_gentau_decayMode.push_back((*_in_gentau_decayMode)[i_gentau]);
			_out_gentau_totNproducts.push_back((*_in_gentau_totNproducts)[i_gentau]);
			_out_gentau_totNgamma.push_back((*_in_gentau_totNgamma)[i_gentau]);
			_out_gentau_totNpiZero.push_back((*_in_gentau_totNpiZero)[i_gentau]);
			_out_gentau_totNcharged.push_back((*_in_gentau_totNcharged)[i_gentau]);			

		}

		_out_gentau_n = _out_gentau_pt.size();


		// loop over 3d clusters

		_out_cl3d_n = _in_cl3d_n;

		for (int i_cl3d=0; i_cl3d<_in_cl3d_n; i_cl3d++){	
		
			_out_cl3d_id.push_back((*_in_cl3d_id)[i_cl3d]);
			_out_cl3d_pt.push_back((*_in_cl3d_pt)[i_cl3d]);
			_out_cl3d_energy.push_back((*_in_cl3d_energy)[i_cl3d]);
			_out_cl3d_eta.push_back((*_in_cl3d_eta)[i_cl3d]);
			_out_cl3d_phi.push_back((*_in_cl3d_phi)[i_cl3d]);
		
			_out_cl3d_clusters_n.push_back((*_in_cl3d_clusters_n)[i_cl3d]);
			_out_cl3d_clusters_id.push_back((*_in_cl3d_clusters_id)[i_cl3d]);
		
			_out_cl3d_showerlength.push_back((*_in_cl3d_showerlength)[i_cl3d]);
			_out_cl3d_coreshowerlength.push_back((*_in_cl3d_coreshowerlength)[i_cl3d]);
			_out_cl3d_firstlayer.push_back((*_in_cl3d_firstlayer)[i_cl3d]);
			_out_cl3d_maxlayer.push_back((*_in_cl3d_maxlayer)[i_cl3d]);
			_out_cl3d_quality.push_back((*_in_cl3d_quality)[i_cl3d]);
		
			_out_cl3d_seetot.push_back((*_in_cl3d_seetot)[i_cl3d]);
			_out_cl3d_seemax.push_back((*_in_cl3d_seemax)[i_cl3d]);
			_out_cl3d_spptot.push_back((*_in_cl3d_spptot)[i_cl3d]);
			_out_cl3d_sppmax.push_back((*_in_cl3d_sppmax)[i_cl3d]);
			_out_cl3d_szz.push_back((*_in_cl3d_szz)[i_cl3d]);
			_out_cl3d_srrtot.push_back((*_in_cl3d_srrtot)[i_cl3d]);
			_out_cl3d_srrmax.push_back((*_in_cl3d_srrmax)[i_cl3d]);
			_out_cl3d_srrmean.push_back((*_in_cl3d_srrmean)[i_cl3d]);
			_out_cl3d_emaxe.push_back((*_in_cl3d_emaxe)[i_cl3d]);
			_out_cl3d_bdteg.push_back((*_in_cl3d_bdteg)[i_cl3d]);
			_out_cl3d_quality.push_back((*_in_cl3d_quality)[i_cl3d]);

		}

		out_tree->Fill();

	}

	out_file->cd();
    out_tree->Write();
    out_file->Close();

	return;

}


void test(int n_events){

  /*TString indir = "root://polgrid4.in2p3.fr//store/user/cmartinp/HGCAL/ZTT_Tauola_All_hadronic_14TeV_TuneCUETP8M1_Pythia8/ZTauTau_PU0_L1TFall17/190228_160553/0000/";

  vector<TString> infiles;
  infiles.push_back(indir+Form("ntuple_Ztt_5k_1.root"));
  infiles.push_back(indir+Form("ntuple_Ztt_5k_2.root"));
  infiles.push_back(indir+Form("ntuple_Ztt_5k_3.root"));
  infiles.push_back(indir+Form("ntuple_Ztt_5k_4.root"));
  infiles.push_back(indir+Form("ntuple_Ztt_5k_5.root"));
  infiles.push_back(indir+Form("ntuple_Ztt_5k_6.root"));
  infiles.push_back(indir+Form("ntuple_Ztt_5k_7.root"));

  TString outfile = "/data_CMS/cms/mperez/HGCal_data/New_TPG/NTuple_ZTT_PU0_skimmed.root";*/

  //NEW BDTeg 12.03.19
  TString indir = "root://polgrid4.in2p3.fr//store/user/cmartinp/HGCAL/ZTT_Tauola_All_hadronic_14TeV_TuneCUETP8M1_Pythia8/ZTauTau_PU0_L1TFall17_new-3dcl-BDTeg/190308_150157/0000/";

  vector<TString> infiles;
  infiles.push_back(indir+Form("ntuple_Ztt_5k_1.root"));
  infiles.push_back(indir+Form("ntuple_Ztt_5k_2.root"));
  infiles.push_back(indir+Form("ntuple_Ztt_5k_3.root"));
  infiles.push_back(indir+Form("ntuple_Ztt_5k_4.root"));
  infiles.push_back(indir+Form("ntuple_Ztt_5k_5.root"));
  infiles.push_back(indir+Form("ntuple_Ztt_5k_6.root"));
  infiles.push_back(indir+Form("ntuple_Ztt_5k_7.root"));

  TString outfile = "/data_CMS/cms/mperez/HGCal_data/New_TPG/NTuple_ZTT_PU0_skimmed_BDTeg.root";
  

  skim_tree(infiles, outfile, n_events);

}


