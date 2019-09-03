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


void plot_eta_res(){

  TString indir = "/data_CMS/cms/mperez/HGCal_data/Aug19/flat/";

  vector<TString> infile_tau;
  //infile_tau.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_flat_calibs_PUcut.root"));
  infile_tau.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_flat_PUcut.root"));
  infile_tau.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU200_flat_PUcut.root"));

  vector<TH1F*> histograms;
  histograms.push_back(single_plot(infile_tau[0],"FlatTree","gentau_matchedSCL3D_eta_Eweighted-gentau_vis_eta","gentau_isMatchedtoSCL3D && gentau_vis_pt>20",50,-0.1,0.1));
  histograms.push_back(single_plot(infile_tau[1],"FlatTree","gentau_matchedSCL3D_eta_Eweighted-gentau_vis_eta","gentau_isMatchedtoSCL3D && gentau_vis_pt>20",50,-0.1,0.1));

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();

  gStyle->SetOptStat(0);  

  histograms[0]->SetTitle(" ");
  histograms[0]->GetXaxis()->SetTitle("#eta (L1 #tau) - #eta (gen. vis. #tau)");
  histograms[0]->GetYaxis()->SetTitle("a. u.");
  //histograms[0]->GetXaxis()->SetRangeUser(xmin,xmax);
  histograms[0]->GetYaxis()->SetTitleOffset(1.5);
  histograms[0]->GetXaxis()->SetTitleOffset(1.2);
  histograms[0]->SetLineColor(kRed);
  histograms[0]->SetLineWidth(2);
  histograms[0]->SetFillColor(0);
  histograms[0]->SetFillColor(0);

  histograms[1]->SetLineColor(kBlue);
  histograms[1]->SetLineWidth(2);
  histograms[1]->SetFillColor(0);
  histograms[1]->SetFillColor(0);

  histograms[0]->DrawNormalized();
  histograms[1]->DrawNormalized("same");

    cout<<"PU=0 "<<histograms[0]->GetMean()<<","<<histograms[0]->GetRMS()<<endl;
  cout<<"PU=200 "<<histograms[1]->GetMean()<<","<<histograms[1]->GetRMS()<<endl;

  TLegend* leg=new TLegend(0.2,0.62,0.3,0.74);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->AddEntry(histograms[0],"PU=0");
  leg->AddEntry(histograms[1],"PU=200");
  
  leg->Draw("same");

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.13,0.91,"#scale[1.3]{CMS #bf{#it{Preliminary}}}");
  tex.DrawLatexNDC(0.67,0.91,"#scale[1.3]{MC Simulation}");
  tex.DrawLatexNDC(0.56,0.84,"#scale[1.2]{HGCal calorimeter-only}");

  //c->SaveAs("../plots/pdf/eta_resol_PUcut.pdf");
  //c->SaveAs("../plots/png/eta_resol_PUcut.png");

  c->SaveAs("../plots/pdf/eta_resol_PUcut_TDR_PU200.pdf");
  c->SaveAs("../plots/png/eta_resol_PUcut_TDR_PU200.png");


  return;

}

void plot_phi_res(){

  TString indir = "/data_CMS/cms/mperez/HGCal_data/Aug19/flat/";

  vector<TString> infile_tau;
  //infile_tau.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_flat_calibs_PUcut.root"));
  infile_tau.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_flat_PUcut.root"));
  infile_tau.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU200_flat_PUcut.root"));

  vector<TH1F*> histograms;
  histograms.push_back(single_plot(infile_tau[0],"FlatTree","gentau_matchedSCL3D_phi_Eweighted-gentau_vis_phi","gentau_isMatchedtoSCL3D && gentau_vis_pt>20",50,-0.1,0.1));
  histograms.push_back(single_plot(infile_tau[1],"FlatTree","gentau_matchedSCL3D_phi_Eweighted-gentau_vis_phi","gentau_isMatchedtoSCL3D && gentau_vis_pt>20",50,-0.1,0.1));

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();

  gStyle->SetOptStat(0);  

  histograms[0]->SetTitle(" ");
  histograms[0]->GetXaxis()->SetTitle("#phi (L1 #tau) - #phi (gen. vis. #tau)");
  histograms[0]->GetYaxis()->SetTitle("a. u.");
  //histograms[0]->GetXaxis()->SetRangeUser(xmin,xmax);
  histograms[0]->GetYaxis()->SetTitleOffset(1.5);
  histograms[0]->GetXaxis()->SetTitleOffset(1.2);
  histograms[0]->SetLineColor(kRed);
  histograms[0]->SetLineWidth(2);
  histograms[0]->SetFillColor(0);
  histograms[0]->SetFillColor(0);
  histograms[1]->SetLineColor(kBlue);
  histograms[1]->SetLineWidth(2);
  histograms[1]->SetFillColor(0);
  histograms[1]->SetFillColor(0);

  histograms[0]->DrawNormalized();
  histograms[1]->DrawNormalized("same");

  cout<<"PU=0 "<<histograms[0]->GetMean()<<","<<histograms[0]->GetRMS()<<endl;
  cout<<"PU=200 "<<histograms[1]->GetMean()<<","<<histograms[1]->GetRMS()<<endl;

  TLegend* leg=new TLegend(0.2,0.62,0.3,0.74);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->AddEntry(histograms[0],"PU=0");
  leg->AddEntry(histograms[1],"PU=200");
  
  leg->Draw("same");

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.13,0.91,"#scale[1.3]{CMS #bf{#it{Preliminary}}}");
  tex.DrawLatexNDC(0.67,0.91,"#scale[1.3]{MC Simulation}");
  tex.DrawLatexNDC(0.56,0.84,"#scale[1.2]{HGCal calorimeter-only}");

  //c->SaveAs("../plots/pdf/phi_resol_PUcut.pdf");
  //c->SaveAs("../plots/png/phi_resol_PUcut.png");

  c->SaveAs("../plots/pdf/phi_resol_PUcut_TDR_PU200.pdf");
  c->SaveAs("../plots/png/phi_resol_PUcut_TDR_PU200.png");


  return;

}



