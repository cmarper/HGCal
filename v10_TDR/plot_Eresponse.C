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


void plot_Eresponse(){

  TString indir = "/data_CMS/cms/mperez/HGCal_data/Aug19/calibrated/";

  vector<TString> infile_tau;
  //infile_tau.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_flat_calibs_PUcut.root"));
  infile_tau.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_flat_PUcut_calibsPU0.root"));
  //infile_tau.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU200_flat_PUcut_calibsPU200.root"));

  vector<TH1F*> histograms;
  histograms.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot_calib/gentau_vis_pt","gentau_isMatchedtoSCL3D",40,0,2));
  histograms.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot/gentau_vis_pt","gentau_isMatchedtoSCL3D",40,0,2));

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();

  gStyle->SetOptStat(0);  

  histograms[0]->SetTitle(" ");
  histograms[0]->GetXaxis()->SetTitle("E_{T} (L1 #tau) / p_{T} (gen. vis. #tau)");
  histograms[0]->GetYaxis()->SetTitle("a. u.");
  //histograms[0]->GetXaxis()->SetRangeUser(xmin,xmax);
  histograms[0]->GetYaxis()->SetTitleOffset(1.5);
  histograms[0]->GetXaxis()->SetTitleOffset(1.2);

  cout<<"Calibrated "<<histograms[0]->GetMean()<<","<<histograms[0]->GetRMS()<<endl;
  cout<<"Raw "<<histograms[1]->GetMean()<<","<<histograms[1]->GetRMS()<<endl;

  TLegend* leg=new TLegend(0.68,0.70,0.85,0.80);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->AddEntry(histograms[1],"Raw");
  leg->AddEntry(histograms[0],"Calibrated");

  for(unsigned int i=0; i<histograms.size(); i++){
    histograms[i]->SetLineWidth(2);
    histograms[i]->SetFillColor(0);
    histograms[i]->SetFillColor(0);
  }

  histograms[0]->SetLineColor(kRed);
  histograms[1]->SetLineColor(kBlue);

  histograms[0]->DrawNormalized();
  histograms[1]->DrawNormalized("same");
  
  leg->Draw("same");

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.13,0.91,"#scale[1.3]{CMS #bf{#it{Preliminary}}}");
  tex.DrawLatexNDC(0.67,0.91,"#scale[1.3]{MC Simulation}");
  tex.DrawLatexNDC(0.56,0.84,"#scale[1.2]{HGCal calorimeter-only}");

  //c->SaveAs("../plots/pdf/response_inclusive_PUcut.pdf");
  //c->SaveAs("../plots/png/response_inclusive_PUcut.png");

  c->SaveAs("../plots/pdf/response_inclusive_PUcut_PU200_TDR.pdf");
  c->SaveAs("../plots/png/response_inclusive_PUcut_PU200_TDR.png");


  return;

}


void plot_Eresponse_perDecayMode(){

  TString indir = "/data_CMS/cms/mperez/HGCal_data/Aug19/calibrated/";

  vector<TString> infile_tau;
  infile_tau.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_flat_PUcut_calibsPU0.root"));
  //infile_tau.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU200_flat_PUcut_calibsPU200.root"));

  vector<TH1F*> histograms;
  histograms.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot_calib/gentau_vis_pt","gentau_decayMode==0 && gentau_isMatchedtoSCL3D",40,0,2));
  histograms.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot_calib/gentau_vis_pt","gentau_decayMode==1 && gentau_isMatchedtoSCL3D",40,0,2));
  histograms.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot_calib/gentau_vis_pt","gentau_decayMode==4 && gentau_isMatchedtoSCL3D",40,0,2));
  histograms.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot_calib/gentau_vis_pt","gentau_decayMode==5 && gentau_isMatchedtoSCL3D",40,0,2));

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();

  gStyle->SetOptStat(0);  

  histograms[3]->SetTitle(" ");
  histograms[3]->GetXaxis()->SetTitle("E_{T} (L1 #tau) / p_{T} (gen. vis. #tau)");
  histograms[3]->GetYaxis()->SetTitle("a. u.");
  //histograms[0]->GetXaxis()->SetRangeUser(xmin,xmax);
  histograms[3]->GetYaxis()->SetTitleOffset(1.5);
  histograms[3]->GetXaxis()->SetTitleOffset(1.2);

  cout<<"0 "<<histograms[0]->GetMean()<<","<<histograms[0]->GetRMS()<<endl;
  cout<<"1 "<<histograms[1]->GetMean()<<","<<histograms[1]->GetRMS()<<endl;
  cout<<"4 "<<histograms[2]->GetMean()<<","<<histograms[2]->GetRMS()<<endl;
  cout<<"5 "<<histograms[3]->GetMean()<<","<<histograms[3]->GetRMS()<<endl;

  TLegend* leg=new TLegend(0.68,0.60,0.85,0.80);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->AddEntry(histograms[0],"1prong");
  leg->AddEntry(histograms[1],"1prong+#pi^{0}");
  leg->AddEntry(histograms[2],"3prongs");
  leg->AddEntry(histograms[3],"3prongs+#pi^{0}");

  for(unsigned int i=0; i<histograms.size(); i++){
    histograms[i]->SetLineWidth(2);
    histograms[i]->SetFillColor(0);
    histograms[i]->SetLineColor(i+2);
    if(i>2)
      histograms[i]->SetLineColor(i+3);
    if(i>6)
      histograms[i]->SetLineColor(i+4);
    histograms[i]->SetFillColor(0);
  }

  histograms[3]->DrawNormalized();
  histograms[1]->DrawNormalized("same");
  histograms[2]->DrawNormalized("same");
  histograms[0]->DrawNormalized("same");
  
  leg->Draw("same");

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.13,0.91,"#scale[1.3]{CMS #bf{#it{Preliminary}}}");
  tex.DrawLatexNDC(0.67,0.91,"#scale[1.3]{MC Simulation}");
  tex.DrawLatexNDC(0.56,0.84,"#scale[1.2]{HGCal calorimeter-only}");

  //c->SaveAs("../plots/pdf/response_inclusive_PUcut.pdf");
  //c->SaveAs("../plots/png/response_inclusive_PUcut.png");

  c->SaveAs("../plots/pdf/response_perDM_PUcut_PU0_TDR.pdf");
  c->SaveAs("../plots/png/response_perDM_PUcut_PU0_TDR.png");


  return;

}

void plot_Eresponse_seed_tot(){

  TString indir = "/data_CMS/cms/mperez/HGCal_data/Aug19/flat/";

  vector<TString> infile_tau;
  //infile_tau.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_flat_calibs_PUcut.root"));
  infile_tau.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_flat_noPUcut_calibs.root"));

  vector<TH1F*> histograms;
  histograms.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_seed/gentau_vis_pt","gentau_isMatchedtoSCL3D",50,0,1.5));
  histograms.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot/gentau_vis_pt","gentau_isMatchedtoSCL3D",50,0,1.5));

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();

  gStyle->SetOptStat(0);  

  histograms[1]->SetTitle(" ");
  histograms[1]->GetXaxis()->SetTitle("E_{T}^{#tau,L1} / p_{T}^{#tau,gen.}");
  histograms[1]->GetYaxis()->SetTitle("a. u.");
  //histograms[1]->GetXaxis()->SetRangeUser(xmin,xmax);
  histograms[1]->GetYaxis()->SetTitleOffset(1.5);
  histograms[1]->GetXaxis()->SetTitleOffset(1.2);

  TLegend* leg=new TLegend(0.68,0.75,0.85,0.85);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->AddEntry(histograms[1],"Total raw E");
  leg->AddEntry(histograms[0],"Seed raw E");

  for(unsigned int i=0; i<histograms.size(); i++){
    histograms[i]->SetLineWidth(2);
    histograms[i]->SetFillColor(0);
    histograms[i]->SetFillColor(0);
  }

  histograms[0]->SetLineColor(kGreen);
  histograms[1]->SetLineColor(kBlue);

  histograms[1]->DrawNormalized();
  histograms[0]->DrawNormalized("same");
  
  leg->Draw("same");

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.13,0.91,"#scale[1.2]{CMS #bf{#it{Preliminary}}}");
  tex.DrawLatexNDC(0.67,0.91,"#scale[1.2]{HGCal calo-only}");

  //c->SaveAs("../plots/pdf/response_seed_tot_PUcut.pdf");
  //c->SaveAs("../plots/png/response_seed_tot_PUcut.png");

  c->SaveAs("../plots/pdf/response_seed_tot_noPUcut.pdf");
  c->SaveAs("../plots/png/response_seed_tot_noPUcut.png");


  return;

}


void plot_Eresponse_bins_pT(){

  TString indir = "/data_CMS/cms/mperez/HGCal_data/Aug19/calibrated/";

  vector<TString> infile_tau;
  //infile_tau.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_flat_calibs_PUcut.root"));
  //infile_tau.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_flat_PUcut_calibs_DM.root"));
  infile_tau.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU200_flat_PUcut_calibsPU200.root"));

  vector<TH1F*> histograms_calib;
  histograms_calib.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot_calib/gentau_vis_pt","gentau_isMatchedtoSCL3D && gentau_vis_pt>=20 && gentau_vis_pt<25",100,0,10));
  histograms_calib.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot_calib/gentau_vis_pt","gentau_isMatchedtoSCL3D && gentau_vis_pt>=25 && gentau_vis_pt<30",100,0,10));
  histograms_calib.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot_calib/gentau_vis_pt","gentau_isMatchedtoSCL3D && gentau_vis_pt>=30 && gentau_vis_pt<35",100,0,10));
  histograms_calib.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot_calib/gentau_vis_pt","gentau_isMatchedtoSCL3D && gentau_vis_pt>=35 && gentau_vis_pt<40",100,0,10));
  histograms_calib.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot_calib/gentau_vis_pt","gentau_isMatchedtoSCL3D && gentau_vis_pt>=40 && gentau_vis_pt<45",100,0,10));
  histograms_calib.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot_calib/gentau_vis_pt","gentau_isMatchedtoSCL3D && gentau_vis_pt>=45 && gentau_vis_pt<50",100,0,10));
  histograms_calib.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot_calib/gentau_vis_pt","gentau_isMatchedtoSCL3D && gentau_vis_pt>=50 && gentau_vis_pt<60",100,0,10));
  histograms_calib.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot_calib/gentau_vis_pt","gentau_isMatchedtoSCL3D && gentau_vis_pt>=60 && gentau_vis_pt<70",100,0,10));
  histograms_calib.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot_calib/gentau_vis_pt","gentau_isMatchedtoSCL3D && gentau_vis_pt>=70 && gentau_vis_pt<=100",100,0,10));

  vector<TH1F*> histograms_raw;
  histograms_raw.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot/gentau_vis_pt","gentau_isMatchedtoSCL3D && gentau_vis_pt>=20 && gentau_vis_pt<25",100,0,10));
  histograms_raw.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot/gentau_vis_pt","gentau_isMatchedtoSCL3D && gentau_vis_pt>=25 && gentau_vis_pt<30",100,0,10));
  histograms_raw.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot/gentau_vis_pt","gentau_isMatchedtoSCL3D && gentau_vis_pt>=30 && gentau_vis_pt<35",100,0,10));
  histograms_raw.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot/gentau_vis_pt","gentau_isMatchedtoSCL3D && gentau_vis_pt>=35 && gentau_vis_pt<40",100,0,10));
  histograms_raw.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot/gentau_vis_pt","gentau_isMatchedtoSCL3D && gentau_vis_pt>=40 && gentau_vis_pt<45",100,0,10));
  histograms_raw.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot/gentau_vis_pt","gentau_isMatchedtoSCL3D && gentau_vis_pt>=45 && gentau_vis_pt<50",100,0,10));
  histograms_raw.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot/gentau_vis_pt","gentau_isMatchedtoSCL3D && gentau_vis_pt>=50 && gentau_vis_pt<60",100,0,10));
  histograms_raw.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot/gentau_vis_pt","gentau_isMatchedtoSCL3D && gentau_vis_pt>=60 && gentau_vis_pt<70",100,0,10));
  histograms_raw.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot/gentau_vis_pt","gentau_isMatchedtoSCL3D && gentau_vis_pt>=70 && gentau_vis_pt<=100",100,0,10));
 
  vector<float> means_calib;
  means_calib.push_back(histograms_calib[0]->GetMean());
  means_calib.push_back(histograms_calib[1]->GetMean());
  means_calib.push_back(histograms_calib[2]->GetMean());
  means_calib.push_back(histograms_calib[3]->GetMean());
  means_calib.push_back(histograms_calib[4]->GetMean());
  means_calib.push_back(histograms_calib[5]->GetMean());
  means_calib.push_back(histograms_calib[6]->GetMean());
  means_calib.push_back(histograms_calib[7]->GetMean());
  means_calib.push_back(histograms_calib[8]->GetMean());

  vector<float> rms_calib;
  rms_calib.push_back(histograms_calib[0]->GetRMS());
  rms_calib.push_back(histograms_calib[1]->GetRMS());
  rms_calib.push_back(histograms_calib[2]->GetRMS());
  rms_calib.push_back(histograms_calib[3]->GetRMS());
  rms_calib.push_back(histograms_calib[4]->GetRMS());
  rms_calib.push_back(histograms_calib[5]->GetRMS());
  rms_calib.push_back(histograms_calib[6]->GetRMS());
  rms_calib.push_back(histograms_calib[7]->GetRMS());
  rms_calib.push_back(histograms_calib[8]->GetRMS());

  vector<float> res_calib;
  res_calib.push_back(rms_calib[0]/means_calib[0]);
  res_calib.push_back(rms_calib[1]/means_calib[1]);
  res_calib.push_back(rms_calib[2]/means_calib[2]);
  res_calib.push_back(rms_calib[3]/means_calib[3]);
  res_calib.push_back(rms_calib[4]/means_calib[4]);
  res_calib.push_back(rms_calib[5]/means_calib[5]);
  res_calib.push_back(rms_calib[6]/means_calib[6]);
  res_calib.push_back(rms_calib[7]/means_calib[7]);
  res_calib.push_back(rms_calib[8]/means_calib[8]);

  vector<float> means_raw;
  means_raw.push_back(histograms_raw[0]->GetMean());
  means_raw.push_back(histograms_raw[1]->GetMean());
  means_raw.push_back(histograms_raw[2]->GetMean());
  means_raw.push_back(histograms_raw[3]->GetMean());
  means_raw.push_back(histograms_raw[4]->GetMean());
  means_raw.push_back(histograms_raw[5]->GetMean());
  means_raw.push_back(histograms_raw[6]->GetMean());
  means_raw.push_back(histograms_raw[7]->GetMean());
  means_raw.push_back(histograms_raw[8]->GetMean());

  vector<float> rms_raw;
  rms_raw.push_back(histograms_raw[0]->GetRMS());
  rms_raw.push_back(histograms_raw[1]->GetRMS());
  rms_raw.push_back(histograms_raw[2]->GetRMS());
  rms_raw.push_back(histograms_raw[3]->GetRMS());
  rms_raw.push_back(histograms_raw[4]->GetRMS());
  rms_raw.push_back(histograms_raw[5]->GetRMS());
  rms_raw.push_back(histograms_raw[6]->GetRMS());
  rms_raw.push_back(histograms_raw[7]->GetRMS());
  rms_raw.push_back(histograms_raw[8]->GetRMS());

  vector<float> res_raw;
  res_raw.push_back(rms_raw[0]/means_raw[0]);
  res_raw.push_back(rms_raw[1]/means_raw[1]);
  res_raw.push_back(rms_raw[2]/means_raw[2]);
  res_raw.push_back(rms_raw[3]/means_raw[3]);
  res_raw.push_back(rms_raw[4]/means_raw[4]);
  res_raw.push_back(rms_raw[5]/means_raw[5]);
  res_raw.push_back(rms_raw[6]/means_raw[6]);
  res_raw.push_back(rms_raw[7]/means_raw[7]);
  res_raw.push_back(rms_raw[8]/means_raw[8]);

  const Int_t nbins = 9;
  float pT_edges[nbins+1] = { 20, 25, 30, 35, 40, 45, 50, 60, 70, 100}; 

  TH1F *h_means_calib = new TH1F("h_means_calib","h_means_calib",nbins,pT_edges);
  for (unsigned int i=0; i<means_calib.size(); i++){ 
    h_means_calib->AddBinContent(i+1,means_calib[i]);
    cout<<"Mean: "<<means_calib[i]<<endl;
    h_means_calib->SetBinError(i+1,0.000001);
  }

  TH1F *h_means_raw = new TH1F("h_means_raw","h_means_raw",nbins,pT_edges);
  for (unsigned int i=0; i<means_raw.size(); i++){ 
    h_means_raw->AddBinContent(i+1,means_raw[i]);
    h_means_raw->SetBinError(i+1,0.000001);
  }

  TH1F *h_rms_calib = new TH1F("h_rms_calib","h_rms_calib",nbins,pT_edges);
  for (unsigned int i=0; i<rms_calib.size(); i++){ 
    h_rms_calib->AddBinContent(i+1,rms_calib[i]);
    cout<<"RMS: "<<rms_calib[i]<<endl;
    h_rms_calib->SetBinError(i+1,0.000001);
  }

  TH1F *h_rms_raw = new TH1F("h_rms_raw","h_rms_raw",nbins,pT_edges);
  for (unsigned int i=0; i<rms_raw.size(); i++){ 
    h_rms_raw->AddBinContent(i+1,rms_raw[i]);
    h_rms_raw->SetBinError(i+1,0.000001);
  }

  TH1F *h_res_calib = new TH1F("h_res_calib","h_res_calib",nbins,pT_edges);
  for (unsigned int i=0; i<res_calib.size(); i++){ 
    h_res_calib->AddBinContent(i+1,res_calib[i]);
    cout<<"Res: "<<res_calib[i]<<endl;
    h_res_calib->SetBinError(i+1,0.000001);
  }

  TH1F *h_res_raw = new TH1F("h_res_raw","h_res_raw",nbins,pT_edges);
  for (unsigned int i=0; i<res_raw.size(); i++){ 
    h_res_raw->AddBinContent(i+1,res_raw[i]);
    h_res_raw->SetBinError(i+1,0.000001);
  }

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();

  gStyle->SetOptStat(0);  

  h_res_calib->SetLineColor(kRed);
  h_res_calib->SetMarkerStyle(20);
  h_res_calib->SetMarkerSize(1.3);
  h_res_calib->SetMarkerColor(kRed);
  h_res_calib->SetLineWidth(2);
  h_res_calib->SetMaximum(0.43);
  h_res_calib->SetMinimum(0.11);
  h_res_calib->GetXaxis()->SetTitleOffset(1.2);
  h_res_calib->GetYaxis()->SetTitleOffset(1.5);
  h_res_calib->GetXaxis()->SetTitle("p_{T} (gen. vis. #tau) [GeV]");
  h_res_calib->GetYaxis()->SetTitle("RMS / < E_{T} (L1 #tau) / p_{T} (gen. vis. #tau) >");
  h_res_calib->SetTitle(" ");
  h_res_calib->Draw("E");

  h_res_raw->SetLineColor(kBlue);
  h_res_raw->SetMarkerStyle(20);
  h_res_raw->SetMarkerSize(1.3);
  h_res_raw->SetMarkerColor(kBlue);
  h_res_raw->SetLineWidth(2);
  h_res_raw->SetTitle(" ");
  h_res_raw->Draw("E same");

  TLegend* leg=new TLegend(0.68,0.70,0.85,0.80);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->AddEntry(h_res_raw,"Raw","apl");
  leg->AddEntry(h_res_calib,"Calibrated","apl");
  
  leg->Draw("same");

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.13,0.91,"#scale[1.3]{CMS #bf{#it{Preliminary}}}");
  tex.DrawLatexNDC(0.67,0.91,"#scale[1.3]{MC Simulation}");
  tex.DrawLatexNDC(0.56,0.84,"#scale[1.2]{HGCal calorimeter-only}");

  //c->SaveAs("../plots/pdf/response_ptbins_PUcut.pdf");
  //c->SaveAs("../plots/png/response_ptbins_PUcut.png");

  c->SaveAs("../plots/pdf/response_ptbins_PUcut_TDR_PU200.pdf");
  c->SaveAs("../plots/png/response_ptbins_PUcut_TDR_PU200.png");

  return;

}



void plot_Eresponse_bins_eta(){

  TString indir = "/data_CMS/cms/mperez/HGCal_data/Aug19/flat/";

  vector<TString> infile_tau;
  //infile_tau.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_flat_calibs_PUcut.root"));
  infile_tau.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_flat_PUcut_calibs_DM.root"));

  vector<TH1F*> histograms_calib;
  histograms_calib.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot_calib/gentau_vis_pt","gentau_isMatchedtoSCL3D && abs(gentau_vis_eta)>=1.6 && abs(gentau_vis_eta)<2.0",100,0,10));
  histograms_calib.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot_calib/gentau_vis_pt","gentau_isMatchedtoSCL3D && abs(gentau_vis_eta)>=2.0 && abs(gentau_vis_eta)<2.4",100,0,10));
  histograms_calib.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot_calib/gentau_vis_pt","gentau_isMatchedtoSCL3D && abs(gentau_vis_eta)>=2.4 && abs(gentau_vis_eta)<2.9",100,0,10));
  
  vector<TH1F*> histograms_raw;
  histograms_raw.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot/gentau_vis_pt","gentau_isMatchedtoSCL3D && abs(gentau_vis_eta)>=1.6 && abs(gentau_vis_eta)<2.0",100,0,10));
  histograms_raw.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot/gentau_vis_pt","gentau_isMatchedtoSCL3D && abs(gentau_vis_eta)>=2.0 && abs(gentau_vis_eta)<2.4",100,0,10));
  histograms_raw.push_back(single_plot(infile_tau,"FlatTree","gentau_matchedSCL3D_pt_tot/gentau_vis_pt","gentau_isMatchedtoSCL3D && abs(gentau_vis_eta)>=2.4 && abs(gentau_vis_eta)<2.9",100,0,10));
 
  vector<float> means_calib;
  means_calib.push_back(histograms_calib[0]->GetMean());
  means_calib.push_back(histograms_calib[1]->GetMean());
  means_calib.push_back(histograms_calib[2]->GetMean());

  vector<float> rms_calib;
  rms_calib.push_back(histograms_calib[0]->GetRMS());
  rms_calib.push_back(histograms_calib[1]->GetRMS());
  rms_calib.push_back(histograms_calib[2]->GetRMS());

  vector<float> res_calib;
  res_calib.push_back(rms_calib[0]/means_calib[0]);
  res_calib.push_back(rms_calib[1]/means_calib[1]);
  res_calib.push_back(rms_calib[2]/means_calib[2]);

  vector<float> means_raw;
  means_raw.push_back(histograms_raw[0]->GetMean());
  means_raw.push_back(histograms_raw[1]->GetMean());
  means_raw.push_back(histograms_raw[2]->GetMean());

  vector<float> rms_raw;
  rms_raw.push_back(histograms_raw[0]->GetRMS());
  rms_raw.push_back(histograms_raw[1]->GetRMS());
  rms_raw.push_back(histograms_raw[2]->GetRMS());

  vector<float> res_raw;
  res_raw.push_back(rms_raw[0]/means_raw[0]);
  res_raw.push_back(rms_raw[1]/means_raw[1]);
  res_raw.push_back(rms_raw[2]/means_raw[2]);

  const Int_t nbins = 3;
  float pT_edges[nbins+1] = { 1.6, 2.0, 2.4, 2.9}; 

  TH1F *h_means_calib = new TH1F("h_means_calib","h_means_calib",nbins,pT_edges);
  for (unsigned int i=0; i<means_calib.size(); i++){ 
    h_means_calib->AddBinContent(i+1,means_calib[i]);
    cout<<"Mean: "<<means_calib[i]<<endl;
    h_means_calib->SetBinError(i+1,0.000001);
  }

  TH1F *h_means_raw = new TH1F("h_means_raw","h_means_raw",nbins,pT_edges);
  for (unsigned int i=0; i<means_raw.size(); i++){ 
    h_means_raw->AddBinContent(i+1,means_raw[i]);
    h_means_raw->SetBinError(i+1,0.000001);
  }

  TH1F *h_rms_calib = new TH1F("h_rms_calib","h_rms_calib",nbins,pT_edges);
  for (unsigned int i=0; i<rms_calib.size(); i++){ 
    h_rms_calib->AddBinContent(i+1,rms_calib[i]);
    cout<<"RMS: "<<rms_calib[i]<<endl;
    h_rms_calib->SetBinError(i+1,0.000001);
  }

  TH1F *h_rms_raw = new TH1F("h_rms_raw","h_rms_raw",nbins,pT_edges);
  for (unsigned int i=0; i<rms_raw.size(); i++){ 
    h_rms_raw->AddBinContent(i+1,rms_raw[i]);
    h_rms_raw->SetBinError(i+1,0.000001);
  }

  TH1F *h_res_calib = new TH1F("h_res_calib","h_res_calib",nbins,pT_edges);
  for (unsigned int i=0; i<res_calib.size(); i++){ 
    h_res_calib->AddBinContent(i+1,res_calib[i]);
    cout<<"Res: "<<res_calib[i]<<endl;
    h_res_calib->SetBinError(i+1,0.000001);
  }

  TH1F *h_res_raw = new TH1F("h_res_raw","h_res_raw",nbins,pT_edges);
  for (unsigned int i=0; i<res_raw.size(); i++){ 
    h_res_raw->AddBinContent(i+1,res_raw[i]);
    h_res_raw->SetBinError(i+1,0.000001);
  }

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();

  gStyle->SetOptStat(0);  

  h_res_calib->SetLineColor(kRed);
  h_res_calib->SetMarkerStyle(20);
  h_res_calib->SetMarkerSize(1.3);
  h_res_calib->SetMarkerColor(kRed);
  h_res_calib->SetLineWidth(2);
  h_res_calib->SetMaximum(0.58);
  h_res_calib->SetMinimum(0.24);
  h_res_calib->GetXaxis()->SetTitleOffset(1.2);
  h_res_calib->GetYaxis()->SetTitleOffset(1.5);
  h_res_calib->GetXaxis()->SetTitle("| #eta (gen. vis. #tau) |");
  h_res_calib->GetYaxis()->SetTitle("RMS / < E_{T} (L1 #tau) / p_{T} (gen. vis. #tau) >");
  h_res_calib->SetTitle(" ");
  h_res_calib->Draw("E");

  h_res_raw->SetLineColor(kBlue);
  h_res_raw->SetMarkerStyle(20);
  h_res_raw->SetMarkerSize(1.3);
  h_res_raw->SetMarkerColor(kBlue);
  h_res_raw->SetLineWidth(2);
  h_res_raw->SetTitle(" ");
  h_res_raw->Draw("E same");

  TLegend* leg=new TLegend(0.68,0.7,0.85,0.8);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->AddEntry(h_res_raw,"Raw","apl");
  leg->AddEntry(h_res_calib,"Calibrated","apl");
  
  leg->Draw("same");

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.13,0.91,"#scale[1.3]{CMS #bf{#it{Preliminary}}}");
  tex.DrawLatexNDC(0.67,0.91,"#scale[1.3]{MC Simulation}");
  tex.DrawLatexNDC(0.56,0.84,"#scale[1.2]{HGCal calorimeter-only}");

  //c->SaveAs("../plots/pdf/response_etabins_PUcut.pdf");
  //c->SaveAs("../plots/png/response_etabins_PUcut.png");

  c->SaveAs("../plots/pdf/response_etabins_PUcut_TDR.pdf");
  c->SaveAs("../plots/png/response_etabins_PUcut_TDR.png");

  return;

}

