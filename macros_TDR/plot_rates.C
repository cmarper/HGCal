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

//#include "Helpers.C"

using namespace std;


void plot_rate_inclusive(){

  TString indir = "/data_CMS_upgrade/mperez/HGCal_data/Sep19/rates/";

  TString infile = indir+Form("rates_Nu_E10_v10_PU200.root");

  TFile* file = TFile::Open(infile,"READ");

  TH1F* rate = (TH1F*)file->Get("RateTau");

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();
  c->SetLogy();

  gStyle->SetOptStat(0);  

  rate->SetTitle(" ");
  rate->GetXaxis()->SetTitle("E_{T}^{L1, #tau} [GeV]");
  rate->GetYaxis()->SetTitle("Rate [kHz]");
  rate->GetYaxis()->SetTitleOffset(1.2);
  rate->GetXaxis()->SetTitleOffset(1.1);
  rate->GetXaxis()->SetTitleSize(0.04);
  rate->GetYaxis()->SetTitleSize(0.04);
  rate->GetXaxis()->SetRangeUser(10,99);

  rate->SetLineWidth(3);
  rate->SetFillColor(0);
  rate->SetLineColor(kRed);

  rate->Draw();

  TLatex tex;
  tex.SetTextSize(0.03);
  //tex.DrawLatexNDC(0.13,0.91,"#scale[1.3]{CMS #bf{#it{Preliminary}}}");
  tex.DrawLatexNDC(0.13,0.92,"#scale[1.3]{CMS #bf{#it{Phase-2 Simulation}}}"); 
  tex.DrawLatexNDC(0.54,0.83,"#scale[1.3]{#bf{HGCal calorimeter-only}}");
  tex.DrawLatexNDC(0.79,0.92,"#scale[1.3]{#bf{200 PU}}");

  tex.DrawLatexNDC(0.48,0.68,"#scale[1.2]{#bf{Inst. luminosity 7.5e34 cm^{-2}s^{-1}}}");

  c->SaveAs("/home/llr/cms/mperez/HGCal/v10_geometry/tau_algorithm/plots/rate_TDR3.pdf");
  c->SaveAs("/home/llr/cms/mperez/HGCal/v10_geometry/tau_algorithm/plots/rate_TDR3.png");

  return;

}

void plot_rate_byDM(){

  TString indir = "/data_CMS_upgrade/mperez/HGCal_data/Sep19/rates/";

  TString infile_dm0 = indir+Form("rates_Nu_E10_v10_PU200_DM0.root");
  TString infile_dm1 = indir+Form("rates_Nu_E10_v10_PU200_DM1.root");
  TString infile_dm2 = indir+Form("rates_Nu_E10_v10_PU200_DM2.root");

  TFile* file_dm0 = TFile::Open(infile_dm0,"READ");
  TFile* file_dm1 = TFile::Open(infile_dm1,"READ");
  TFile* file_dm2 = TFile::Open(infile_dm2,"READ");

  TH1F* rate_dm0 = (TH1F*)file_dm0->Get("RateTau");
  TH1F* rate_dm1 = (TH1F*)file_dm1->Get("RateTau");
  TH1F* rate_dm2 = (TH1F*)file_dm2->Get("RateTau");

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();
  c->SetLogy();

  gStyle->SetOptStat(0);  

  rate_dm0->SetTitle(" ");
  rate_dm0->GetXaxis()->SetTitle("E_{T}^{L1, #tau} [GeV]");
  rate_dm0->GetYaxis()->SetTitle("Rate [kHz]");
  rate_dm0->GetYaxis()->SetTitleOffset(1.2);
  rate_dm0->GetXaxis()->SetTitleOffset(1.1);
  rate_dm0->GetXaxis()->SetTitleSize(0.04);
  rate_dm0->GetYaxis()->SetTitleSize(0.04);
  rate_dm0->GetXaxis()->SetRangeUser(10,99);

  rate_dm0->SetLineWidth(3);
  rate_dm0->SetFillColor(0);
  rate_dm0->SetLineColor(kRed);

  rate_dm1->SetLineWidth(3);
  rate_dm1->SetFillColor(0);
  rate_dm1->SetLineColor(kBlue);

  rate_dm2->SetLineWidth(3);
  rate_dm2->SetFillColor(0);
  rate_dm2->SetLineColor(kOrange);

  rate_dm0->Draw();
  rate_dm1->Draw("same");
  rate_dm2->Draw("same");

  TLegend* leg=new TLegend(0.50,0.45,0.75,0.64);
  leg->SetHeader("Predicted #tau decay mode:");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);

  leg->AddEntry(rate_dm0, "1-prong", "apl");
  leg->AddEntry(rate_dm1, "1-prong + #pi^{0}'s", "apl");
  leg->AddEntry(rate_dm2, "3-prongs (+ #pi^{0}'s)", "apl");

  leg->Draw("same");

  TLatex tex;
  tex.SetTextSize(0.03);
  //tex.DrawLatexNDC(0.13,0.91,"#scale[1.3]{CMS #bf{#it{Preliminary}}}");
  tex.DrawLatexNDC(0.13,0.92,"#scale[1.3]{CMS #bf{#it{Phase-2 Simulation}}}"); 
  tex.DrawLatexNDC(0.54,0.83,"#scale[1.3]{#bf{HGCal calorimeter-only}}");
  tex.DrawLatexNDC(0.79,0.92,"#scale[1.3]{#bf{200 PU}}");

  tex.DrawLatexNDC(0.48,0.68,"#scale[1.2]{#bf{Inst. luminosity 7.5e34 cm^{-2}s^{-1}}}");

  c->SaveAs("/home/llr/cms/mperez/HGCal/v10_geometry/tau_algorithm/plots/rate_DM_TDR3.pdf");
  c->SaveAs("/home/llr/cms/mperez/HGCal/v10_geometry/tau_algorithm/plots/rate_DM_TDR3.png");

  return;

}

void plot_rate_byDM_norm_inclusive(){

  TString indir = "/data_CMS_upgrade/mperez/HGCal_data/Sep19/rates/";

  TString infile_all = indir+Form("rates_Nu_E10_v10_PU200.root");
  TString infile_dm0 = indir+Form("rates_Nu_E10_v10_PU200_DM0_norm_inclusive.root");
  TString infile_dm1 = indir+Form("rates_Nu_E10_v10_PU200_DM1_norm_inclusive.root");
  TString infile_dm2 = indir+Form("rates_Nu_E10_v10_PU200_DM2_norm_inclusive.root");

  TFile* file_all = TFile::Open(infile_all,"READ");
  TFile* file_dm0 = TFile::Open(infile_dm0,"READ");
  TFile* file_dm1 = TFile::Open(infile_dm1,"READ");
  TFile* file_dm2 = TFile::Open(infile_dm2,"READ");

  TH1F* rate_all = (TH1F*)file_all->Get("RateTau");
  TH1F* rate_dm0 = (TH1F*)file_dm0->Get("RateTau");
  TH1F* rate_dm1 = (TH1F*)file_dm1->Get("RateTau");
  TH1F* rate_dm2 = (TH1F*)file_dm2->Get("RateTau");

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();
  c->SetLogy();

  gStyle->SetOptStat(0);  

  rate_all->SetTitle(" ");
  rate_all->GetXaxis()->SetTitle("E_{T}^{L1, #tau} [GeV]");
  rate_all->GetYaxis()->SetTitle("Single-#tau rate [kHz]");
  rate_all->GetYaxis()->SetTitleOffset(1.2);
  rate_all->GetXaxis()->SetTitleOffset(1.1);
  rate_all->GetXaxis()->SetTitleSize(0.04);
  rate_all->GetYaxis()->SetTitleSize(0.04);
  rate_all->GetXaxis()->SetRangeUser(10,99);

  rate_all->SetLineWidth(3);
  rate_all->SetFillColor(0);
  rate_all->SetLineColor(kBlack);

  rate_dm0->SetLineWidth(3);
  rate_dm0->SetFillColor(0);
  rate_dm0->SetLineColor(kRed);

  rate_dm1->SetLineWidth(3);
  rate_dm1->SetFillColor(0);
  rate_dm1->SetLineColor(kBlue);

  rate_dm2->SetLineWidth(3);
  rate_dm2->SetFillColor(0);
  rate_dm2->SetLineColor(kOrange);

  rate_all->Draw();
  rate_dm0->Draw("same");
  rate_dm1->Draw("same");
  rate_dm2->Draw("same");

  TLegend* leg=new TLegend(0.46,0.45,0.71,0.64);
  leg->SetHeader("L1 #tau reconstructed decay mode:");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);

  leg->AddEntry(rate_all, "All", "apl");
  leg->AddEntry(rate_dm0, "1-prong", "apl");
  leg->AddEntry(rate_dm1, "1-prong + #pi^{0}'s", "apl");
  leg->AddEntry(rate_dm2, "3-prong (+ #pi^{0}'s)", "apl");

  leg->Draw("same");

  TLatex tex;
  tex.SetTextSize(0.03);
  //tex.DrawLatexNDC(0.13,0.91,"#scale[1.3]{CMS #bf{#it{Preliminary}}}");
  tex.DrawLatexNDC(0.13,0.92,"#scale[1.3]{CMS #bf{#it{Phase-2 Simulation}}}"); 
  tex.DrawLatexNDC(0.54,0.85,"#scale[1.3]{#bf{HGCal calorimeter-only}}");
  tex.DrawLatexNDC(0.79,0.92,"#scale[1.3]{#bf{200 PU}}");
 
  tex.DrawLatexNDC(0.48,0.77,"#scale[1.2]{#bf{Inst. luminosity 7.5e34 cm^{-2}s^{-1}}}");
  tex.DrawLatexNDC(0.55,0.69,"#scale[1.2]{#bf{1.6 < | #eta_{gen,#tau} | < 2.9}}");

  c->SaveAs("/home/llr/cms/mperez/HGCal/v10_geometry/tau_algorithm/plots/rate_DM_norm_inclusive_TDR4.pdf");
  c->SaveAs("/home/llr/cms/mperez/HGCal/v10_geometry/tau_algorithm/plots/rate_DM_norm_inclusive_TDR4.png");

  return;

}
