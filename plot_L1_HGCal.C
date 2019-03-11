#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <iostream>
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


void plot_inclusive (  TString var, double bins, double xmin, double xmax, TString axislabel, TString filename, bool drawnorm ){

  TString indir = "/data_CMS/cms/mperez/HGCal_data/L1taus/";

  vector<TString> infiles;
  infiles.push_back(indir+Form("NTuple_ZTT_PU0_matched_0p3.root"));
 
  vector<TH1F*> histograms;
  histograms.push_back(single_plot(infiles,"MatchedTree",var,"1",bins,xmin,xmax));

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();

  gStyle->SetOptStat(0);  

  histograms[0]->SetTitle(" ");
  histograms[0]->GetXaxis()->SetTitle(axislabel);
  if(drawnorm) histograms[0]->GetYaxis()->SetTitle("Normalized events");
  else if(!drawnorm) histograms[0]->GetYaxis()->SetTitle("Events");
  histograms[0]->GetXaxis()->SetRangeUser(xmin,xmax);
  histograms[0]->GetYaxis()->SetTitleOffset(1.7);
  histograms[0]->GetXaxis()->SetTitleOffset(1.2);

  TLegend* leg=new TLegend(0.68,0.75,0.83,0.85);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(histograms[0],"#tau");

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

  if(drawnorm) histograms[0]->DrawNormalized();
  else if(!drawnorm) histograms[0]->Draw();
  
  for(unsigned int i=1; i<histograms.size(); i++){    
    if(drawnorm) histograms[i]->DrawNormalized("same");
    else if(!drawnorm) histograms[i]->Draw("same");
  }
  
  //leg->Draw("same");

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.13,0.91,"#scale[1.2]{HGCal #bf{#it{Simulation}}}");
  //tex.DrawLatexNDC(0.6,0.91,"#scale[1.2]{#font[40]{Fall17 Single-#tau, PU = 0}}");
  tex.DrawLatexNDC(0.62,0.91,"#scale[1.2]{#font[40]{Fall17 Z#rightarrow#tau#tau, PU = 0}}");

  if(drawnorm) {
    c->SaveAs("../plots/pdf/"+filename+"_inclusive_norm.pdf");
    c->SaveAs("../plots/png/"+filename+"_inclusive_norm.png");
  }

  if(!drawnorm) {
    c->SaveAs("../plots/pdf/"+filename+"_inclusive.pdf");
    c->SaveAs("../plots/png/"+filename+"_inclusive.png");
  }

  return;

}


void plot_perdecaymode (  TString var, double bins, double xmin, double xmax, TString axislabel, TString filename, bool drawnorm, float mymax ){

  //github.com/PFCal-dev/cmssw/blob/hgc-tpg-devel-CMSSW_10_5_0_pre1/L1Trigger/L1THGCalUtilities/plugins/ntuples/HGCalTriggerNtupleGenTau.cc#L315-L327
  /* 
  0->1prong
  1->1prong+pi0
  4->3prongs
  5->3prongs+pi0
  */

  TString basecut = "1";

  TString indir = "/data_CMS/cms/mperez/HGCal_data/L1taus/";

  vector<TString> infiles;
  infiles.push_back(indir+Form("NTuple_ZTT_PU0_matched_0p3.root"));
 
  vector<TH1F*> histograms;
  histograms.push_back(single_plot(infiles,"MatchedTree",var,basecut+" && gentau_decayMode==0", bins,xmin,xmax));
  histograms.push_back(single_plot(infiles,"MatchedTree",var,basecut+" && gentau_decayMode==1", bins,xmin,xmax));
  histograms.push_back(single_plot(infiles,"MatchedTree",var,basecut+" && gentau_decayMode==4", bins,xmin,xmax));
  histograms.push_back(single_plot(infiles,"MatchedTree",var,basecut+" && gentau_decayMode==5", bins,xmin,xmax));

  TCanvas* c=new TCanvas("c","c",650,600);
  c->SetLeftMargin(0.12);
  c->SetGridx();
  c->SetGridy();

  gStyle->SetOptStat(0);  

  double ymax = mymax;

  histograms[0]->SetTitle(" ");
  histograms[0]->GetXaxis()->SetTitle(axislabel);
  if(drawnorm) histograms[0]->GetYaxis()->SetTitle("Normalized events");
  else if(!drawnorm) histograms[0]->GetYaxis()->SetTitle("Events");
  histograms[0]->GetXaxis()->SetRangeUser(xmin,xmax);
  //if(!drawnorm) histograms[0]->GetYaxis()->SetRangeUser(0,ymax);
  histograms[0]->GetYaxis()->SetRangeUser(0,ymax);
  histograms[0]->GetYaxis()->SetTitleOffset(1.7);
  histograms[0]->GetXaxis()->SetTitleOffset(1.2);

  TLegend* leg=new TLegend(0.68,0.65,0.85,0.85);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
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

  if(drawnorm) histograms[0]->DrawNormalized();
  else if(!drawnorm) histograms[0]->Draw();
  
  for(unsigned int i=1; i<histograms.size(); i++){    
    if(drawnorm) histograms[i]->DrawNormalized("same");
    else if(!drawnorm) histograms[i]->Draw("same");
  }
  
  leg->Draw("same");

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.13,0.91,"#scale[1.2]{HGCal #bf{#it{Simulation}}}");
  //tex.DrawLatexNDC(0.6,0.91,"#scale[1.2]{#font[40]{Fall17 Single-#tau, PU = 0}}");
  tex.DrawLatexNDC(0.62,0.91,"#scale[1.2]{#font[40]{Fall17 Z#rightarrow#tau#tau, PU = 0}}");

  if(drawnorm) {
    c->SaveAs("../plots/pdf/"+filename+"_perdecaymode_norm.pdf");
    c->SaveAs("../plots/png/"+filename+"_perdecaymode_norm.png");
  }

  if(!drawnorm) {
    c->SaveAs("../plots/pdf/"+filename+"_perdecaymode.pdf");
    c->SaveAs("../plots/png/"+filename+"_perdecaymode.png");
  }

  return;

}

/*void plot_inclusive_all(){

  plot_perdecaymode("cl3d_n",40,0,40,"3D cluster N","cl3d_n",false);
  plot_perdecaymode("cl3d_pt",40,0,20,"3D cluster p_{T} [GeV]","cl3d_pt",false);
  plot_perdecaymode("cl3d_energy",40,0,80,"3D cluster E [GeV]","cl3d_energy",false);
  plot_perdecaymode("cl3d_eta",45,-4,4,"3D cluster #eta","cl3d_eta",false);
  plot_perdecaymode("cl3d_phi",45,-4,4,"3D cluster #phi","cl3d_phi",false);
  plot_perdecaymode("cl3d_clusters_n",50,0,150,"3D cluster N clusters","cl3d_clusters_n",false);
  plot_perdecaymode("cl3d_showerlength",30,0,60,"3D cluster shower length","cl3d_showerlength",false);
  plot_perdecaymode("cl3d_coreshowerlength",40,0,40,"3D cluster core shower length","cl3d_coreshowerlength",false);
  plot_perdecaymode("cl3d_firstlayer",60,0,60,"3D cluster first layer","cl3d_firstlayer",false);
  plot_perdecaymode("cl3d_maxlayer",60,0,60,"3D cluster max. layer","cl3d_maxlayer",false);

}*/

void plot_perdecaymode_all(){

  plot_perdecaymode("dR_cl3d_gentau",100,0,1,"dR(gentau,cl3d)","dR_cl3d_gentau",false, 10);
  plot_perdecaymode("gentau_isMatched",2,0,2,"is Matched gentau (dR=0.3)","gentau_isMatched",false, 10);
  plot_perdecaymode("gentau_nMatchedcl3d",10,0,10,"number of matched cl3d per gentau","gentau_nMatchedcl3d",false, 10);
  plot_perdecaymode("gentau_PtMaxMatchedcl3d/gentau_vis_pt",50,0,2,"p^{T}_{max}(3dcl) / p^{T}(gentau_vis)","res_ptmax",true, 500);
  plot_perdecaymode("gentau_PtTotMatchedcl3d/gentau_vis_pt",10,0,10,"p^{T}_{tot}(3dcl)/p^{T}(gentau)","res_pttot",false, 10);

}