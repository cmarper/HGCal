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


void plot_puBDT ( ){

  TString indir = "/data_CMS/cms/mperez/HGCal_data/Aug19/skimmed/";

  vector<TString> infile_tau;
  infile_tau.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_skimmed.root"));

  vector<TString> infile_pu;
  infile_pu.push_back(indir+Form("ntuple_Nu_E10_v10_PU200_files100to150_skimmed.root"));
 
  vector<TH1F*> histograms;
  histograms.push_back(single_plot(infile_pu,"SkimmedTree","cl3d_puBDT","1",30,-1.5,1.5));
  histograms.push_back(single_plot(infile_tau,"SkimmedTree","cl3d_puBDT","1",30,-1.5,1.5));

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();

  gStyle->SetOptStat(0);  

  histograms[0]->SetTitle(" ");
  histograms[0]->GetXaxis()->SetTitle("3D-cluster BDT score");
  histograms[0]->GetYaxis()->SetTitle("a. u.");
  //histograms[0]->GetXaxis()->SetRangeUser(xmin,xmax);
  histograms[0]->GetYaxis()->SetTitleOffset(1.4);
  histograms[0]->GetXaxis()->SetTitleOffset(1.1);

  TLegend* leg=new TLegend(0.4,0.62,0.50,0.74);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->AddEntry(histograms[1],"#tau_{h}");
  leg->AddEntry(histograms[0],"PU");

  for(unsigned int i=0; i<histograms.size(); i++){
    histograms[i]->SetLineWidth(2);
    histograms[i]->SetFillColor(0);
    histograms[i]->SetFillColor(0);
  }

  histograms[0]->SetLineColor(kGreen+1);
  histograms[1]->SetLineColor(kMagenta+1);

  histograms[0]->DrawNormalized();
  
  for(unsigned int i=1; i<histograms.size(); i++){    
    histograms[i]->DrawNormalized("same");
  }
  
  leg->Draw("same");

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.13,0.91,"#scale[1.3]{CMS #bf{#it{Preliminary}}}");
  tex.DrawLatexNDC(0.67,0.91,"#scale[1.3]{MC Simulation}");
  tex.DrawLatexNDC(0.53,0.84,"#scale[1.2]{HGCal calorimeter-only}");

  c->SaveAs("../plots/pdf/PUvsTauBDT_TDR.pdf");
  c->SaveAs("../plots/png/PUvsTauBDT_TDR.png");


  return;

}


