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

void plot_conf_matrix(){

  TString indir = "/data_CMS/cms/mperez/HGCal_data/Aug19/confmatrices/";

  TString infile = indir+Form("confmatrices_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_flat_PUcut_calibsPU0_DM.root");

  TFile* file = TFile::Open(infile,"READ");

  TH2F* conf_matrix = (TH2F*)file->Get("conf_matrix");

  TString DMs[4] = { "1-prong", "1-prong+#pi^{0}'s", "3-prongs", "3-prongs+#pi^{0}'s"};

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.17);
  c->SetRightMargin(0.12);
  c->SetGridx();
  c->SetGridy();

  gStyle->SetOptStat(0);  

  conf_matrix->SetTitle(" ");
  conf_matrix->GetXaxis()->SetTitle("Predicted");
  conf_matrix->GetYaxis()->SetTitle("True");
  conf_matrix->GetYaxis()->SetTitleOffset(1.6);
  conf_matrix->GetXaxis()->SetTitleOffset(1.3);
  conf_matrix->GetXaxis()->SetRangeUser(0,4);
  conf_matrix->GetYaxis()->SetRangeUser(0,4);
  //conf_matrix->SetMaximum(0.1);

  for(int i=0; i<4; i++){
    conf_matrix->GetXaxis()->SetBinLabel(i+1,DMs[i]);
    conf_matrix->GetYaxis()->SetBinLabel(i+1,DMs[i]);
  }

  conf_matrix->GetXaxis()->SetLabelSize(0.045);
  conf_matrix->GetYaxis()->SetLabelSize(0.045);

  conf_matrix->Draw("colz");

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.18,0.91,"#scale[1.3]{CMS #bf{#it{Preliminary}}}");
  tex.DrawLatexNDC(0.65,0.91,"#scale[1.3]{MC Simulation}");
  tex.DrawLatexNDC(0.54,0.84,"#scale[1.2]#color[White]{HGCal calorimeter-only}");

  c->SaveAs("../plots/pdf/conf_matrix.pdf");
  c->SaveAs("../plots/png/conf_matrix.png");

  return;

}


void plot_pred0(){

  TString indir = "/data_CMS/cms/mperez/HGCal_data/Aug19/confmatrices/";

  TString infile = indir+Form("confmatrices_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_flat_PUcut_calibsPU0_DM.root");

  TFile* file = TFile::Open(infile,"READ");

  TH1F* h_pred0_true0 = (TH1F*)file->Get("prob_pred0_true0");
  TH1F* h_pred0_true1 = (TH1F*)file->Get("prob_pred0_true1");
  TH1F* h_pred0_true4 = (TH1F*)file->Get("prob_pred0_true4");
  TH1F* h_pred0_true5 = (TH1F*)file->Get("prob_pred0_true5");

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();

  gStyle->SetOptStat(0);  

  h_pred0_true1->SetTitle(" ");
  h_pred0_true1->GetXaxis()->SetTitle("Predicted 1-prong probability");
  h_pred0_true1->GetYaxis()->SetTitle("a. u.");
  h_pred0_true1->GetYaxis()->SetTitleOffset(1.4);
  h_pred0_true1->GetXaxis()->SetTitleOffset(1.2);
  h_pred0_true1->GetXaxis()->SetRangeUser(0,1);

  TLegend* leg=new TLegend(0.58,0.56,0.75,0.75);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->AddEntry(h_pred0_true0,"True 1-prong");
  leg->AddEntry(h_pred0_true1,"True 1-prong+#pi^{0}'s");
  leg->AddEntry(h_pred0_true4,"True 3-prongs");
  leg->AddEntry(h_pred0_true5,"True 3-prongs+#pi^{0}'s");

  h_pred0_true0->SetLineWidth(2);
  h_pred0_true0->SetFillColor(0);
  h_pred0_true0->SetFillColor(0);

  h_pred0_true1->SetLineWidth(2);
  h_pred0_true1->SetFillColor(0);
  h_pred0_true1->SetFillColor(0);

  h_pred0_true4->SetLineWidth(2);
  h_pred0_true4->SetFillColor(0);
  h_pred0_true4->SetFillColor(0);

  h_pred0_true5->SetLineWidth(2);
  h_pred0_true5->SetFillColor(0);
  h_pred0_true5->SetFillColor(0);

  h_pred0_true0->SetLineColor(kRed);
  h_pred0_true1->SetLineColor(kBlue);
  h_pred0_true4->SetLineColor(kGreen);
  h_pred0_true5->SetLineColor(kMagenta);

  h_pred0_true1->DrawNormalized();
  h_pred0_true0->DrawNormalized("same");
  h_pred0_true4->DrawNormalized("same");
  h_pred0_true5->DrawNormalized("same");
  
  leg->Draw("same");

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.13,0.91,"#scale[1.3]{CMS #bf{#it{Preliminary}}}");
  tex.DrawLatexNDC(0.67,0.91,"#scale[1.3]{MC Simulation}");
  tex.DrawLatexNDC(0.52,0.84,"#scale[1.2]{HGCal calorimeter-only}");

  c->SaveAs("../plots/pdf/probs_pred0.pdf");
  c->SaveAs("../plots/png/probs_pred0.png");

  return;

}


void plot_pred1(){

  TString indir = "/data_CMS/cms/mperez/HGCal_data/Aug19/confmatrices/";

  TString infile = indir+Form("confmatrices_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_flat_PUcut_calibsPU0_DM.root");

  TFile* file = TFile::Open(infile,"READ");

  TH1F* h_pred1_true0 = (TH1F*)file->Get("prob_pred1_true0");
  TH1F* h_pred1_true1 = (TH1F*)file->Get("prob_pred1_true1");
  TH1F* h_pred1_true4 = (TH1F*)file->Get("prob_pred1_true4");
  TH1F* h_pred1_true5 = (TH1F*)file->Get("prob_pred1_true5");

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();

  gStyle->SetOptStat(0);  

  h_pred1_true1->SetTitle(" ");
  h_pred1_true1->GetXaxis()->SetTitle("Predicted 1-prong+#pi^{0}'s probability");
  h_pred1_true1->GetYaxis()->SetTitle("a. u.");
  h_pred1_true1->GetYaxis()->SetTitleOffset(1.4);
  h_pred1_true1->GetXaxis()->SetTitleOffset(1.2);
  h_pred1_true1->GetXaxis()->SetRangeUser(0,1);

  TLegend* leg=new TLegend(0.34,0.61,0.51,0.80);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->AddEntry(h_pred1_true0,"True 1-prong");
  leg->AddEntry(h_pred1_true1,"True 1-prong+#pi^{0}'s");
  leg->AddEntry(h_pred1_true4,"True 3-prongs");
  leg->AddEntry(h_pred1_true5,"True 3-prongs+#pi^{0}'s");

  h_pred1_true0->SetLineWidth(2);
  h_pred1_true0->SetFillColor(0);
  h_pred1_true0->SetFillColor(0);

  h_pred1_true1->SetLineWidth(2);
  h_pred1_true1->SetFillColor(0);
  h_pred1_true1->SetFillColor(0);

  h_pred1_true4->SetLineWidth(2);
  h_pred1_true4->SetFillColor(0);
  h_pred1_true4->SetFillColor(0);

  h_pred1_true5->SetLineWidth(2);
  h_pred1_true5->SetFillColor(0);
  h_pred1_true5->SetFillColor(0);

  h_pred1_true0->SetLineColor(kRed);
  h_pred1_true1->SetLineColor(kBlue);
  h_pred1_true4->SetLineColor(kGreen);
  h_pred1_true5->SetLineColor(kMagenta);

  h_pred1_true1->DrawNormalized();
  h_pred1_true0->DrawNormalized("same");
  h_pred1_true4->DrawNormalized("same");
  h_pred1_true5->DrawNormalized("same");
  
  leg->Draw("same");

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.13,0.91,"#scale[1.3]{CMS #bf{#it{Preliminary}}}");
  tex.DrawLatexNDC(0.67,0.91,"#scale[1.3]{MC Simulation}");
  tex.DrawLatexNDC(0.42,0.84,"#scale[1.2]{HGCal calorimeter-only}");

  c->SaveAs("../plots/pdf/probs_pred1.pdf");
  c->SaveAs("../plots/png/probs_pred1.png");

  return;

}


void plot_pred4(){

  TString indir = "/data_CMS/cms/mperez/HGCal_data/Aug19/confmatrices/";

  TString infile = indir+Form("confmatrices_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_flat_PUcut_calibsPU0_DM.root");

  TFile* file = TFile::Open(infile,"READ");

  TH1F* h_pred4_true0 = (TH1F*)file->Get("prob_pred4_true0");
  TH1F* h_pred4_true1 = (TH1F*)file->Get("prob_pred4_true1");
  TH1F* h_pred4_true4 = (TH1F*)file->Get("prob_pred4_true4");
  TH1F* h_pred4_true5 = (TH1F*)file->Get("prob_pred4_true5");

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();

  gStyle->SetOptStat(0);  

  h_pred4_true1->SetTitle(" ");
  h_pred4_true1->GetXaxis()->SetTitle("Predicted 3-prongs probability");
  h_pred4_true1->GetYaxis()->SetTitle("a. u.");
  h_pred4_true1->GetYaxis()->SetTitleOffset(1.4);
  h_pred4_true1->GetXaxis()->SetTitleOffset(1.2);
  h_pred4_true1->GetXaxis()->SetRangeUser(0,1);

  TLegend* leg=new TLegend(0.58,0.56,0.75,0.75);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->AddEntry(h_pred4_true0,"True 1-prong");
  leg->AddEntry(h_pred4_true1,"True 1-prong+#pi^{0}'s");
  leg->AddEntry(h_pred4_true4,"True 3-prongs");
  leg->AddEntry(h_pred4_true5,"True 3-prongs+#pi^{0}'s");

  h_pred4_true0->SetLineWidth(2);
  h_pred4_true0->SetFillColor(0);
  h_pred4_true0->SetFillColor(0);

  h_pred4_true1->SetLineWidth(2);
  h_pred4_true1->SetFillColor(0);
  h_pred4_true1->SetFillColor(0);

  h_pred4_true4->SetLineWidth(2);
  h_pred4_true4->SetFillColor(0);
  h_pred4_true4->SetFillColor(0);

  h_pred4_true5->SetLineWidth(2);
  h_pred4_true5->SetFillColor(0);
  h_pred4_true5->SetFillColor(0);

  h_pred4_true0->SetLineColor(kRed);
  h_pred4_true1->SetLineColor(kBlue);
  h_pred4_true4->SetLineColor(kGreen);
  h_pred4_true5->SetLineColor(kMagenta);

  h_pred4_true1->DrawNormalized();
  h_pred4_true0->DrawNormalized("same");
  h_pred4_true4->DrawNormalized("same");
  h_pred4_true5->DrawNormalized("same");
  
  leg->Draw("same");

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.13,0.91,"#scale[1.3]{CMS #bf{#it{Preliminary}}}");
  tex.DrawLatexNDC(0.67,0.91,"#scale[1.3]{MC Simulation}");
  tex.DrawLatexNDC(0.52,0.84,"#scale[1.2]{HGCal calorimeter-only}");

  c->SaveAs("../plots/pdf/probs_pred4.pdf");
  c->SaveAs("../plots/png/probs_pred4.png");

  return;

}


void plot_pred5(){

  TString indir = "/data_CMS/cms/mperez/HGCal_data/Aug19/confmatrices/";

  TString infile = indir+Form("confmatrices_RelValDiTau_Pt20To100_Etam1p6Tom2p9_v10_PU0_flat_PUcut_calibsPU0_DM.root");

  TFile* file = TFile::Open(infile,"READ");

  TH1F* h_pred5_true0 = (TH1F*)file->Get("prob_pred5_true0");
  TH1F* h_pred5_true1 = (TH1F*)file->Get("prob_pred5_true1");
  TH1F* h_pred5_true4 = (TH1F*)file->Get("prob_pred5_true4");
  TH1F* h_pred5_true5 = (TH1F*)file->Get("prob_pred5_true5");

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();

  gStyle->SetOptStat(0);  

  h_pred5_true0->SetTitle(" ");
  h_pred5_true0->GetXaxis()->SetTitle("Predicted 3-prongs+#pi^{0}'s probability");
  h_pred5_true0->GetYaxis()->SetTitle("a. u.");
  h_pred5_true0->GetYaxis()->SetTitleOffset(1.4);
  h_pred5_true0->GetXaxis()->SetTitleOffset(1.2);
  h_pred5_true0->GetXaxis()->SetRangeUser(0,1);

  TLegend* leg=new TLegend(0.58,0.56,0.75,0.75);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->AddEntry(h_pred5_true0,"True 1-prong");
  leg->AddEntry(h_pred5_true1,"True 1-prong+#pi^{0}'s");
  leg->AddEntry(h_pred5_true4,"True 3-prongs");
  leg->AddEntry(h_pred5_true5,"True 3-prongs+#pi^{0}'s");

  h_pred5_true0->SetLineWidth(2);
  h_pred5_true0->SetFillColor(0);
  h_pred5_true0->SetFillColor(0);

  h_pred5_true1->SetLineWidth(2);
  h_pred5_true1->SetFillColor(0);
  h_pred5_true1->SetFillColor(0);

  h_pred5_true4->SetLineWidth(2);
  h_pred5_true4->SetFillColor(0);
  h_pred5_true4->SetFillColor(0);

  h_pred5_true5->SetLineWidth(2);
  h_pred5_true5->SetFillColor(0);
  h_pred5_true5->SetFillColor(0);

  h_pred5_true0->SetLineColor(kRed);
  h_pred5_true1->SetLineColor(kBlue);
  h_pred5_true4->SetLineColor(kGreen);
  h_pred5_true5->SetLineColor(kMagenta);

  h_pred5_true0->DrawNormalized();
  h_pred5_true1->DrawNormalized("same");
  h_pred5_true4->DrawNormalized("same");
  h_pred5_true5->DrawNormalized("same");
  
  leg->Draw("same");

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.13,0.91,"#scale[1.3]{CMS #bf{#it{Preliminary}}}");
  tex.DrawLatexNDC(0.67,0.91,"#scale[1.3]{MC Simulation}");
  tex.DrawLatexNDC(0.52,0.84,"#scale[1.2]{HGCal calorimeter-only}");

  c->SaveAs("../plots/pdf/probs_pred5.pdf");
  c->SaveAs("../plots/png/probs_pred5.png");

  return;

}