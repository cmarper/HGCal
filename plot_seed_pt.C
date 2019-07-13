#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <cmath>
#include <typeinfo>

#include <TF1.h>
#include <TF2.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TEventList.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TGaxis.h>
#include <TMarker.h>
#include <TText.h>
#include <TMath.h>
#include <TBox.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <THStack.h>
#include <TSpectrum2.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TNtuple.h>
#include <TProfile.h>
#include <TFrame.h>
#include <TBenchmark.h>
#include <TInterpreter.h>
#include <TMultiGraph.h>
#include <TMatrixD.h>
#include <TPaletteAxis.h>
#include <TPaveText.h>
#include "TSystem.h"
#include "TPaletteAxis.h"

void plot_cl1_pt(){

	TFile file_pu0("/data_CMS/cms/mperez/HGCal_data/May19/clustered/NTuple_ZTT_PU0_clustered_pTseed0_pTsec0_eta3_phi5.root","READ");
	TFile file_pu200("/data_CMS/cms/mperez/HGCal_data/May19/clustered/NTuple_ZTT_PU200_clustered_pTseed0_pTsec0_eta3_phi5.root","READ");

	TH1F* cl1_pu0  = (TH1F*)file_pu0.Get("h_pT_cl1");
	TH1F* cl2_pu0  = (TH1F*)file_pu0.Get("h_pT_cl2");
	TH1F* cl3_pu0  = (TH1F*)file_pu0.Get("h_pT_cl3");
	TH1F* cl4_pu0  = (TH1F*)file_pu0.Get("h_pT_cl4");

	TH1F* cl1_pu200  = (TH1F*)file_pu200.Get("h_pT_cl1");
	TH1F* cl2_pu200  = (TH1F*)file_pu200.Get("h_pT_cl2");
	TH1F* cl3_pu200  = (TH1F*)file_pu200.Get("h_pT_cl3");
	TH1F* cl4_pu200  = (TH1F*)file_pu200.Get("h_pT_cl4");

	TLatex tex;
  	tex.SetTextSize(0.03);
  	tex.DrawLatexNDC(0.13,0.91,"#scale[1.2]{CMS #bf{#it{HGCal Simulation}}}");
  	tex.DrawLatexNDC(0.65,0.91,"#scale[1.2]{#font[40]{Fall17 Z#rightarrow#tau#tau}}");

	TCanvas* c1=new TCanvas("c1","c1",650,600);
  	c1->SetLeftMargin(0.12);
  	c1->SetGridx();
  	c1->SetGridy();

  	gStyle->SetOptStat(0);  

  	cl1_pu0->SetTitle(" ");
  	cl1_pu0->GetXaxis()->SetTitle("1st cl3d p{T}");
  	cl1_pu0->GetYaxis()->SetTitle("Normalized events");	
  	cl1_pu0->GetYaxis()->SetRangeUser(0,20);
  	cl1_pu0->GetYaxis()->SetTitleOffset(1.7);
  	cl1_pu0->GetXaxis()->SetTitleOffset(1.2);

  	TLegend* leg=new TLegend(0.48,0.65,0.65,0.85);
  	leg->SetBorderSize(0);
  	leg->SetTextSize(0.03);
  	leg->SetFillColor(0);
  	leg->AddEntry(cl1_pu0,"PU=0");
  	leg->AddEntry(cl1_pu200,"PU=200");

  	cl1_pu0->SetLineWidth(2);
  	cl1_pu0->SetFillColor(0);
  	cl1_pu0->SetLineColor(kRed);

  	cl1_pu200->SetLineWidth(2);
  	cl1_pu200->SetFillColor(0);
  	cl1_pu200->SetLineColor(kBlue);

  	cl1_pu0->DrawNormalized();
  	cl1_pu200->DrawNormalized("same");
  
  	leg->Draw("same");

  	c1->SaveAs("../plots/pdf/cl1_pt.pdf");
  	c1->SaveAs("../plots/png/cl1_pt.png");


  	TCanvas* c2=new TCanvas("c2","c2",650,600);
  	c2->SetLeftMargin(0.12);
  	c2->SetGridx();
  	c2->SetGridy();

  	gStyle->SetOptStat(0);  

  	cl2_pu0->SetTitle(" ");
  	cl2_pu0->GetXaxis()->SetTitle("2nd cl3d p{T}");
  	cl2_pu0->GetYaxis()->SetTitle("Normalized events");	
  	cl2_pu0->GetYaxis()->SetRangeUser(0,20);
  	cl2_pu0->GetYaxis()->SetTitleOffset(1.7);
  	cl2_pu0->GetXaxis()->SetTitleOffset(1.2);

  	TLegend* leg=new TLegend(0.48,0.65,0.65,0.85);
  	leg->SetBorderSize(0);
  	leg->SetTextSize(0.03);
  	leg->SetFillColor(0);
  	leg->AddEntry(cl2_pu0,"PU=0");
  	leg->AddEntry(cl2_pu200,"PU=200");

  	cl2_pu0->SetLineWidth(2);
  	cl2_pu0->SetFillColor(0);
  	cl2_pu0->SetLineColor(kRed);

  	cl2_pu200->SetLineWidth(2);
  	cl2_pu200->SetFillColor(0);
  	cl2_pu200->SetLineColor(kBlue);

  	cl2_pu0->DrawNormalized();
  	cl2_pu200->DrawNormalized("same");
  
  	leg->Draw("same");

  	c2->SaveAs("../plots/pdf/cl2_pt.pdf");
  	c2->SaveAs("../plots/png/cl2_pt.png");

  	TCanvas* c3=new TCanvas("c3","c3",650,600);
  	c3->SetLeftMargin(0.12);
  	c3->SetGridx();
  	c3->SetGridy();

  	gStyle->SetOptStat(0);  

  	cl3_pu0->SetTitle(" ");
  	cl3_pu0->GetXaxis()->SetTitle("3rd cl3d p{T}");
  	cl3_pu0->GetYaxis()->SetTitle("Normalized events");	
  	cl3_pu0->GetYaxis()->SetRangeUser(0,20);
  	cl3_pu0->GetYaxis()->SetTitleOffset(1.7);
  	cl3_pu0->GetXaxis()->SetTitleOffset(1.2);

  	TLegend* leg=new TLegend(0.48,0.65,0.65,0.85);
  	leg->SetBorderSize(0);
  	leg->SetTextSize(0.03);
  	leg->SetFillColor(0);
  	leg->AddEntry(cl3_pu0,"PU=0");
  	leg->AddEntry(cl3_pu200,"PU=200");

  	cl3_pu0->SetLineWidth(2);
  	cl3_pu0->SetFillColor(0);
  	cl3_pu0->SetLineColor(kRed);

  	cl3_pu200->SetLineWidth(2);
  	cl3_pu200->SetFillColor(0);
  	cl3_pu200->SetLineColor(kBlue);

  	cl3_pu0->DrawNormalized();
  	cl3_pu200->DrawNormalized("same");
  
  	leg->Draw("same");

  	c3->SaveAs("../plots/pdf/cl3_pt.pdf");
  	c3->SaveAs("../plots/png/cl3_pt.png");

  	return;

}