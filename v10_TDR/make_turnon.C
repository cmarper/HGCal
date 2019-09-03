#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <iostream>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TROOT.h>
#include <sstream>
#include <TBranchElement.h>
#include <fstream>
#include <TGraphAsymmErrors.h>
#include <stdio.h>
#include <math.h>

using namespace std;

void make_turnon(){

	TString InputFileName = "/data_CMS/cms/mperez/HGCal_data/Aug19/calibrated/ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU200_flat_PUcut_calibsPU200.root";
	TFile f(InputFileName.Data(),"READ");
	TTree* inTree = (TTree*)f.Get("FlatTree");

	TString OutputFileName = "/data_CMS/cms/mperez/HGCal_data/Aug19/turnons/turnon_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU200_flat_PUcut_calibsPU200.root";
	TFile turnOns(OutputFileName.Data(),"RECREATE");

	// Input variables
	float in_gentau_pt;
	float in_gentau_eta;
	float in_gentau_phi;
	float in_gentau_energy;
	float in_gentau_vis_pt;
	float in_gentau_vis_eta;
	float in_gentau_vis_phi;
	float in_gentau_vis_energy;
	int   in_gentau_decayMode;
	bool  in_gentau_isMatchedtoSCL3D;

	float in_gentau_matchedSCL3D_pt_tot;
	float in_gentau_matchedSCL3D_calib_factor;
	float in_gentau_matchedSCL3D_pt_tot_calib;
	float in_gentau_matchedSCL3D_eta_Eweighted;
	float in_gentau_matchedSCL3D_phi_Eweighted;

	inTree->SetBranchAddress("gentau_pt",			&in_gentau_pt);
	inTree->SetBranchAddress("gentau_eta",			&in_gentau_eta);
	inTree->SetBranchAddress("gentau_phi",			&in_gentau_phi);
	inTree->SetBranchAddress("gentau_energy",		&in_gentau_energy);
	inTree->SetBranchAddress("gentau_vis_pt",		&in_gentau_vis_pt);
	inTree->SetBranchAddress("gentau_vis_eta",		&in_gentau_vis_eta);
	inTree->SetBranchAddress("gentau_vis_phi",		&in_gentau_vis_phi);
	inTree->SetBranchAddress("gentau_vis_energy",	&in_gentau_vis_energy);
	inTree->SetBranchAddress("gentau_decayMode",	&in_gentau_decayMode);
	inTree->SetBranchAddress("gentau_isMatchedtoSCL3D",	&in_gentau_isMatchedtoSCL3D);

	inTree->SetBranchAddress("gentau_matchedSCL3D_pt_tot", 		&in_gentau_matchedSCL3D_pt_tot);
	inTree->SetBranchAddress("gentau_matchedSCL3D_calib_factor", 	&in_gentau_matchedSCL3D_calib_factor);
	inTree->SetBranchAddress("gentau_matchedSCL3D_pt_tot_calib", 	&in_gentau_matchedSCL3D_pt_tot_calib);
	inTree->SetBranchAddress("gentau_matchedSCL3D_eta_Eweighted", 	&in_gentau_matchedSCL3D_eta_Eweighted);
	inTree->SetBranchAddress("gentau_matchedSCL3D_phi_Eweighted", 	&in_gentau_matchedSCL3D_phi_Eweighted);

	double binning[16] = {18,20,22,24,26,28,30,32,35,40,45,50,60,70,80,100};

	TH1F* pt = new TH1F("pt","pt",15,binning);
	TH1F* pt_pass = new TH1F("pt_pass","pt_pass",15,binning);

	TH1F* pt_DM0 = new TH1F("pt_DM0","pt_DM0",15,binning);
	TH1F* pt_pass_DM0 = new TH1F("pt_pass_DM0","pt_pass_DM0",15,binning);

	TH1F* pt_DM1 = new TH1F("pt_DM1","pt_DM1",15,binning);
	TH1F* pt_pass_DM1 = new TH1F("pt_pass_DM1","pt_pass_DM1",15,binning);

	TH1F* pt_DM4 = new TH1F("pt_DM4","pt_DM4",15,binning);
	TH1F* pt_pass_DM4 = new TH1F("pt_pass_DM4","pt_pass_DM4",15,binning);

	TH1F* pt_DM5 = new TH1F("pt_DM5","pt_DM5",15,binning);
	TH1F* pt_pass_DM5 = new TH1F("pt_pass_DM5","pt_pass_DM5",15,binning);

	TString s_Threshold = "30";
	float f_Threshold = 30;

	for(int i = 0 ; i < inTree->GetEntries() ; ++i){

		inTree->GetEntry(i);

		pt->Fill(in_gentau_vis_pt);

		if(in_gentau_decayMode==0) pt_DM0->Fill(in_gentau_vis_pt);
		if(in_gentau_decayMode==1) pt_DM1->Fill(in_gentau_vis_pt);
		if(in_gentau_decayMode==4) pt_DM4->Fill(in_gentau_vis_pt);
		if(in_gentau_decayMode==5) pt_DM5->Fill(in_gentau_vis_pt);

		if(in_gentau_matchedSCL3D_pt_tot_calib>=f_Threshold) pt_pass->Fill(in_gentau_vis_pt);

		if(in_gentau_matchedSCL3D_pt_tot_calib>=f_Threshold && in_gentau_decayMode==0) pt_pass_DM0->Fill(in_gentau_vis_pt);
		if(in_gentau_matchedSCL3D_pt_tot_calib>=f_Threshold && in_gentau_decayMode==1) pt_pass_DM1->Fill(in_gentau_vis_pt);
		if(in_gentau_matchedSCL3D_pt_tot_calib>=f_Threshold && in_gentau_decayMode==4) pt_pass_DM4->Fill(in_gentau_vis_pt);
		if(in_gentau_matchedSCL3D_pt_tot_calib>=f_Threshold && in_gentau_decayMode==5) pt_pass_DM5->Fill(in_gentau_vis_pt);

	}

	TGraphAsymmErrors* turnOn = new TGraphAsymmErrors(pt_pass,pt,"cp");
	TGraphAsymmErrors* turnOn_DM0 = new TGraphAsymmErrors(pt_pass_DM0,pt_DM0,"cp");
	TGraphAsymmErrors* turnOn_DM1 = new TGraphAsymmErrors(pt_pass_DM1,pt_DM1,"cp");
	TGraphAsymmErrors* turnOn_DM4 = new TGraphAsymmErrors(pt_pass_DM4,pt_DM4,"cp");
	TGraphAsymmErrors* turnOn_DM5 = new TGraphAsymmErrors(pt_pass_DM5,pt_DM5,"cp");
  	
  	turnOn->Write();
  	turnOn_DM0->Write();
  	turnOn_DM1->Write();
  	turnOn_DM4->Write();
  	turnOn_DM5->Write();

}
