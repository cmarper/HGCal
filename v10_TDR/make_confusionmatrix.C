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

void make_confmatrix(){

	TString InputFileName = "/data_CMS/cms/mperez/HGCal_data/Aug19/DM/ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_flat_PUcut_calibsPU0_DM.root";
	TFile f(InputFileName.Data(),"READ");
	TTree* inTree = (TTree*)f.Get("FlatTree");

	TString OutputFileName = "/data_CMS/cms/mperez/HGCal_data/Aug19/confmatrices/confmatrices_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_flat_PUcut_calibsPU0_DM.root";
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

	int in_gentau_matchedSCL3D_predDM;
	float in_gentau_matchedSCL3D_probDM0;
	float in_gentau_matchedSCL3D_probDM1;
	float in_gentau_matchedSCL3D_probDM4;
	float in_gentau_matchedSCL3D_probDM5;

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

	inTree->SetBranchAddress("gentau_matchedSCL3D_pt_tot", 			&in_gentau_matchedSCL3D_pt_tot);
	inTree->SetBranchAddress("gentau_matchedSCL3D_calib_factor", 	&in_gentau_matchedSCL3D_calib_factor);
	inTree->SetBranchAddress("gentau_matchedSCL3D_pt_tot_calib", 	&in_gentau_matchedSCL3D_pt_tot_calib);
	inTree->SetBranchAddress("gentau_matchedSCL3D_eta_Eweighted", 	&in_gentau_matchedSCL3D_eta_Eweighted);
	inTree->SetBranchAddress("gentau_matchedSCL3D_phi_Eweighted", 	&in_gentau_matchedSCL3D_phi_Eweighted);

	inTree->SetBranchAddress("gentau_matchedSCL3D_predDM",  &in_gentau_matchedSCL3D_predDM);
	inTree->SetBranchAddress("gentau_matchedSCL3D_probDM0", &in_gentau_matchedSCL3D_probDM0);
	inTree->SetBranchAddress("gentau_matchedSCL3D_probDM1", &in_gentau_matchedSCL3D_probDM1);
	inTree->SetBranchAddress("gentau_matchedSCL3D_probDM4", &in_gentau_matchedSCL3D_probDM4);
	inTree->SetBranchAddress("gentau_matchedSCL3D_probDM5", &in_gentau_matchedSCL3D_probDM5);

	float pred0_true0 = 0;
	float pred0_true1 = 0;
	float pred0_true4 = 0;
	float pred0_true5 = 0;

	float pred1_true0 = 0;
	float pred1_true1 = 0;
	float pred1_true4 = 0;
	float pred1_true5 = 0;

	float pred4_true0 = 0;
	float pred4_true1 = 0;
	float pred4_true4 = 0;
	float pred4_true5 = 0;

	float pred5_true0 = 0;
	float pred5_true1 = 0;
	float pred5_true4 = 0;
	float pred5_true5 = 0;

	TH2F* conf_matrix = new TH2F("conf_matrix","conf_matrix",4,0,4,4,0,4);

	TH1F* prob_pred0_true0 = new TH1F("prob_pred0_true0","prob_pred0_true0",10,0,1);
	TH1F* prob_pred0_true1 = new TH1F("prob_pred0_true1","prob_pred0_true1",10,0,1);
	TH1F* prob_pred0_true4 = new TH1F("prob_pred0_true4","prob_pred0_true4",10,0,1);
	TH1F* prob_pred0_true5 = new TH1F("prob_pred0_true5","prob_pred0_true5",10,0,1);

	TH1F* prob_pred1_true0 = new TH1F("prob_pred1_true0","prob_pred1_true0",10,0,1);
	TH1F* prob_pred1_true1 = new TH1F("prob_pred1_true1","prob_pred1_true1",10,0,1);
	TH1F* prob_pred1_true4 = new TH1F("prob_pred1_true4","prob_pred1_true4",10,0,1);
	TH1F* prob_pred1_true5 = new TH1F("prob_pred1_true5","prob_pred1_true5",10,0,1);

	TH1F* prob_pred4_true0 = new TH1F("prob_pred4_true0","prob_pred4_true0",10,0,1);
	TH1F* prob_pred4_true1 = new TH1F("prob_pred4_true1","prob_pred4_true1",10,0,1);
	TH1F* prob_pred4_true4 = new TH1F("prob_pred4_true4","prob_pred4_true4",10,0,1);
	TH1F* prob_pred4_true5 = new TH1F("prob_pred4_true5","prob_pred4_true5",10,0,1);

	TH1F* prob_pred5_true0 = new TH1F("prob_pred5_true0","prob_pred5_true0",10,0,1);
	TH1F* prob_pred5_true1 = new TH1F("prob_pred5_true1","prob_pred5_true1",10,0,1);
	TH1F* prob_pred5_true4 = new TH1F("prob_pred5_true4","prob_pred5_true4",10,0,1);
	TH1F* prob_pred5_true5 = new TH1F("prob_pred5_true5","prob_pred5_true5",10,0,1);

	for(int i = 0 ; i < inTree->GetEntries() ; ++i){

		inTree->GetEntry(i);

		if(in_gentau_decayMode==0 && in_gentau_matchedSCL3D_predDM == 0) pred0_true0 += 1;
		if(in_gentau_decayMode==1 && in_gentau_matchedSCL3D_predDM == 0) pred0_true1 += 1;
		if(in_gentau_decayMode==4 && in_gentau_matchedSCL3D_predDM == 0) pred0_true4 += 1;
		if(in_gentau_decayMode==5 && in_gentau_matchedSCL3D_predDM == 0) pred0_true5 += 1;

		if(in_gentau_decayMode==0 && in_gentau_matchedSCL3D_predDM == 1) pred1_true0 += 1;
		if(in_gentau_decayMode==1 && in_gentau_matchedSCL3D_predDM == 1) pred1_true1 += 1;
		if(in_gentau_decayMode==4 && in_gentau_matchedSCL3D_predDM == 1) pred1_true4 += 1;
		if(in_gentau_decayMode==5 && in_gentau_matchedSCL3D_predDM == 1) pred1_true5 += 1;

		if(in_gentau_decayMode==0 && in_gentau_matchedSCL3D_predDM == 4) pred4_true0 += 1;
		if(in_gentau_decayMode==1 && in_gentau_matchedSCL3D_predDM == 4) pred4_true1 += 1;
		if(in_gentau_decayMode==4 && in_gentau_matchedSCL3D_predDM == 4) pred4_true4 += 1;
		if(in_gentau_decayMode==5 && in_gentau_matchedSCL3D_predDM == 4) pred4_true5 += 1;

		if(in_gentau_decayMode==0 && in_gentau_matchedSCL3D_predDM == 5) pred5_true0 += 1;
		if(in_gentau_decayMode==1 && in_gentau_matchedSCL3D_predDM == 5) pred5_true1 += 1;
		if(in_gentau_decayMode==4 && in_gentau_matchedSCL3D_predDM == 5) pred5_true4 += 1;
		if(in_gentau_decayMode==5 && in_gentau_matchedSCL3D_predDM == 5) pred5_true5 += 1;

		if(in_gentau_decayMode==0){
			prob_pred0_true0->Fill(in_gentau_matchedSCL3D_probDM0);
			prob_pred1_true0->Fill(in_gentau_matchedSCL3D_probDM1);
			prob_pred4_true0->Fill(in_gentau_matchedSCL3D_probDM4);
			prob_pred5_true0->Fill(in_gentau_matchedSCL3D_probDM5);
		}

		if(in_gentau_decayMode==1){
			prob_pred0_true1->Fill(in_gentau_matchedSCL3D_probDM0);
			prob_pred1_true1->Fill(in_gentau_matchedSCL3D_probDM1);
			prob_pred4_true1->Fill(in_gentau_matchedSCL3D_probDM4);
			prob_pred5_true1->Fill(in_gentau_matchedSCL3D_probDM5);
		}

		if(in_gentau_decayMode==4){
			prob_pred0_true4->Fill(in_gentau_matchedSCL3D_probDM0);
			prob_pred1_true4->Fill(in_gentau_matchedSCL3D_probDM1);
			prob_pred4_true4->Fill(in_gentau_matchedSCL3D_probDM4);
			prob_pred5_true4->Fill(in_gentau_matchedSCL3D_probDM5);
		}

		if(in_gentau_decayMode==5){
			prob_pred0_true5->Fill(in_gentau_matchedSCL3D_probDM0);
			prob_pred1_true5->Fill(in_gentau_matchedSCL3D_probDM1);
			prob_pred4_true5->Fill(in_gentau_matchedSCL3D_probDM4);
			prob_pred5_true5->Fill(in_gentau_matchedSCL3D_probDM5);
		}

	}

	float total = 
		pred0_true0 + pred0_true1 + pred0_true4 + pred0_true5 + 
		pred1_true0 + pred1_true1 + pred1_true4 + pred1_true5 +
		pred4_true0 + pred4_true1 + pred4_true4 + pred4_true5 +
		pred5_true0 + pred5_true1 + pred5_true4 + pred5_true5;

	float eff = (pred0_true0 + pred1_true1 + pred4_true4 + pred5_true5) / total;
	cout<<"eff "<<eff<<endl;

	conf_matrix->Fill(0.5,0.5,pred0_true0/total);
	conf_matrix->Fill(0.5,1.5,pred0_true1/total);
	conf_matrix->Fill(0.5,2.5,pred0_true4/total);
	conf_matrix->Fill(0.5,3.5,pred0_true5/total);

	conf_matrix->Fill(1.5,0.5,pred1_true0/total);
	conf_matrix->Fill(1.5,1.5,pred1_true1/total);
	conf_matrix->Fill(1.5,2.5,pred1_true4/total);
	conf_matrix->Fill(1.5,3.5,pred1_true5/total);

	conf_matrix->Fill(2.5,0.5,pred4_true0/total);
	conf_matrix->Fill(2.5,1.5,pred4_true1/total);
	conf_matrix->Fill(2.5,2.5,pred4_true4/total);
	conf_matrix->Fill(2.5,3.5,pred4_true5/total);

	conf_matrix->Fill(3.5,0.5,pred5_true0/total);
	conf_matrix->Fill(3.5,1.5,pred5_true1/total);
	conf_matrix->Fill(3.5,2.5,pred5_true4/total);
	conf_matrix->Fill(3.5,3.5,pred5_true5/total);
  	
  	conf_matrix->Write();

  	prob_pred0_true0->Write();
	prob_pred0_true1->Write();
	prob_pred0_true4->Write();
	prob_pred0_true5->Write();

	prob_pred1_true0->Write();
	prob_pred1_true1->Write();
	prob_pred1_true4->Write();
	prob_pred1_true5->Write();

	prob_pred4_true0->Write();
	prob_pred4_true1->Write();
	prob_pred4_true4->Write();
	prob_pred4_true5->Write();

	prob_pred5_true0->Write();
	prob_pred5_true1->Write();
	prob_pred5_true4->Write();
	prob_pred5_true5->Write();

}