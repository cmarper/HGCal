#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <cmath>
#include <typeinfo>

#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TLegend.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <TChain.h>
#include <TMath.h>
#include <TF1.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>

#include "Helpers.C"

using namespace std;


void plot_turnon_inclusive_pu0(){

  TString indir = "/data_CMS/cms/mperez/HGCal_data/Aug19/turnons/";

  TString infile_tau_pu0 = indir+Form("turnon_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_flat_PUcut_calibs_DM.root");

  TFile* file_tau_pu0 = TFile::Open(infile_tau_pu0,"READ");

  TGraphAsymmErrors* turnon_pu0_inclusive = (TGraphAsymmErrors*)file_tau_pu0->Get("divide_pt_pass_by_pt");

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();

  gStyle->SetOptStat(0);  

  turnon_pu0_inclusive->SetTitle(" ");
  turnon_pu0_inclusive->GetXaxis()->SetTitle("p_{T} (gen. vis. #tau) [GeV]");
  turnon_pu0_inclusive->GetYaxis()->SetTitle("Efficiency");
  turnon_pu0_inclusive->GetYaxis()->SetTitleOffset(1.4);
  turnon_pu0_inclusive->GetXaxis()->SetTitleOffset(1.2);
  turnon_pu0_inclusive->GetXaxis()->SetRangeUser(0,100);

  TLegend* leg=new TLegend(0.68,0.55,0.85,0.65);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->AddEntry(turnon_pu0_inclusive,"PU = 0");

  turnon_pu0_inclusive->SetLineWidth(2);
  turnon_pu0_inclusive->SetFillColor(0);
  turnon_pu0_inclusive->SetFillColor(0);

  turnon_pu0_inclusive->SetLineColor(kRed);

  turnon_pu0_inclusive->Draw();
  
  leg->Draw("same");

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.13,0.91,"#scale[1.3]{CMS #bf{#it{Preliminary}}}");
  tex.DrawLatexNDC(0.67,0.91,"#scale[1.3]{MC Simulation}");
  tex.DrawLatexNDC(0.52,0.64,"#scale[1.2]{HGCal calorimeter-only}");

  tex.DrawLatexNDC(0.57,0.48,"#font[40]{E_{T} (L1 #tau) > 30 GeV}");

  c->SaveAs("../plots/pdf/turnon_pu0_TDR.pdf");
  c->SaveAs("../plots/png/turnon_pu0_TDR.png");


  return;

}



void plot_turnon_inclusive_allpu(){

  TString indir = "/data_CMS/cms/mperez/HGCal_data/Aug19/turnons/";

  TString infile_tau_pu0 = indir+Form("turnon_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_flat_PUcut_calibsPU0.root");
  TString infile_tau_pu200 = indir+Form("turnon_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU200_flat_PUcut_calibsPU200.root");

  TFile* file_tau_pu0 = TFile::Open(infile_tau_pu0,"READ");
  TFile* file_tau_pu200 = TFile::Open(infile_tau_pu200,"READ");

  TGraphAsymmErrors* turnon_pu0_inclusive = (TGraphAsymmErrors*)file_tau_pu0->Get("divide_pt_pass_by_pt");
  TGraphAsymmErrors* turnon_pu200_inclusive = (TGraphAsymmErrors*)file_tau_pu200->Get("divide_pt_pass_by_pt");

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();

  gStyle->SetOptStat(0);  

  turnon_pu0_inclusive->SetTitle(" ");
  turnon_pu0_inclusive->GetXaxis()->SetTitle("p_{T} (gen. vis. #tau) [GeV]");
  turnon_pu0_inclusive->GetYaxis()->SetTitle("Efficiency");
  turnon_pu0_inclusive->GetYaxis()->SetTitleOffset(1.4);
  turnon_pu0_inclusive->GetXaxis()->SetTitleOffset(1.2);
  turnon_pu0_inclusive->GetXaxis()->SetRangeUser(20,100);

  TLegend* leg=new TLegend(0.58,0.45,0.75,0.55);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->AddEntry(turnon_pu0_inclusive,"PU = 0");
  leg->AddEntry(turnon_pu200_inclusive,"PU = 200");

  turnon_pu0_inclusive->SetLineWidth(2);
  turnon_pu0_inclusive->SetFillColor(0);
  turnon_pu0_inclusive->SetFillColor(0);

  turnon_pu200_inclusive->SetLineWidth(2);
  turnon_pu200_inclusive->SetFillColor(0);
  turnon_pu200_inclusive->SetFillColor(0);

  turnon_pu0_inclusive->SetLineColor(kRed);
  turnon_pu200_inclusive->SetLineColor(kBlue);

  turnon_pu0_inclusive->Draw();
  turnon_pu200_inclusive->Draw("same");
  
  leg->Draw("same");

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.13,0.91,"#scale[1.3]{CMS #bf{#it{Preliminary}}}");
  tex.DrawLatexNDC(0.67,0.91,"#scale[1.3]{MC Simulation}");
  tex.DrawLatexNDC(0.50,0.6,"#scale[1.2]{HGCal calorimeter-only}");

  tex.DrawLatexNDC(0.57,0.32,"#font[40]{E_{T} (L1 #tau) > 30 GeV}");

  c->SaveAs("../plots/pdf/turnon_allpu_TDR_2.pdf");
  c->SaveAs("../plots/png/turnon_allpu_TDR_2.png");

  return;

}


void plot_turnon_DM_pu0(){

  TString indir = "/data_CMS/cms/mperez/HGCal_data/Aug19/turnons/";

  TString infile_tau_pu0 = indir+Form("turnon_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_flat_PUcut_calibs_DM.root");

  TFile* file_tau_pu0 = TFile::Open(infile_tau_pu0,"READ");

  TGraphAsymmErrors* turnon_pu0_DM0 = (TGraphAsymmErrors*)file_tau_pu0->Get("divide_pt_pass_DM0_by_pt_DM0");
  TGraphAsymmErrors* turnon_pu0_DM1 = (TGraphAsymmErrors*)file_tau_pu0->Get("divide_pt_pass_DM1_by_pt_DM1");
  TGraphAsymmErrors* turnon_pu0_DM4 = (TGraphAsymmErrors*)file_tau_pu0->Get("divide_pt_pass_DM4_by_pt_DM4");
  TGraphAsymmErrors* turnon_pu0_DM5 = (TGraphAsymmErrors*)file_tau_pu0->Get("divide_pt_pass_DM5_by_pt_DM5");

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();

  gStyle->SetOptStat(0);  

  turnon_pu0_DM0->SetTitle(" ");
  turnon_pu0_DM0->GetXaxis()->SetTitle("p_{T} (gen. vis. #tau) [GeV]");
  turnon_pu0_DM0->GetYaxis()->SetTitle("Efficiency");
  turnon_pu0_DM0->GetYaxis()->SetTitleOffset(1.4);
  turnon_pu0_DM0->GetXaxis()->SetTitleOffset(1.2);
  turnon_pu0_DM0->GetXaxis()->SetRangeUser(20,100);

  TLegend* leg=new TLegend(0.56,0.4,0.73,0.58);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->AddEntry(turnon_pu0_DM0,"1-prong");
  leg->AddEntry(turnon_pu0_DM1,"1-prong+#pi^{0}'s");
  leg->AddEntry(turnon_pu0_DM4,"3-prongs");
  leg->AddEntry(turnon_pu0_DM5,"3-prongs+#pi^{0}'s");

  turnon_pu0_DM0->SetLineWidth(2);
  turnon_pu0_DM0->SetFillColor(0);
  turnon_pu0_DM0->SetFillColor(0);

  turnon_pu0_DM1->SetLineWidth(2);
  turnon_pu0_DM1->SetFillColor(0);
  turnon_pu0_DM1->SetFillColor(0);

  turnon_pu0_DM4->SetLineWidth(2);
  turnon_pu0_DM4->SetFillColor(0);
  turnon_pu0_DM4->SetFillColor(0);

  turnon_pu0_DM5->SetLineWidth(2);
  turnon_pu0_DM5->SetFillColor(0);
  turnon_pu0_DM5->SetFillColor(0);

  turnon_pu0_DM0->SetLineColor(kRed);
  turnon_pu0_DM1->SetLineColor(kBlue);
  turnon_pu0_DM4->SetLineColor(kGreen);
  turnon_pu0_DM5->SetLineColor(kMagenta);

  turnon_pu0_DM0->Draw();
  turnon_pu0_DM1->Draw("same");
  turnon_pu0_DM4->Draw("same");
  turnon_pu0_DM5->Draw("same");
  
  leg->Draw("same");

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.13,0.91,"#scale[1.3]{CMS #bf{#it{Preliminary}}}");
  tex.DrawLatexNDC(0.67,0.91,"#scale[1.3]{MC Simulation}");
  tex.DrawLatexNDC(0.52,0.64,"#scale[1.2]{HGCal calorimeter-only}");

  tex.DrawLatexNDC(0.57,0.28,"#font[40]{E_{T} (L1 #tau) > 30 GeV}");

  c->SaveAs("../plots/pdf/turnon_pu0_DM_TDR.pdf");
  c->SaveAs("../plots/png/turnon_pu0_DM_TDR.png");


  return;

}

