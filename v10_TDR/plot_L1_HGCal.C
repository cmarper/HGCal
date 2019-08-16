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


void plot_inclusive (  TString var, double bins, double xmin, double xmax, TString axislabel, TString filename, bool drawnorm ){

  TString indir = "/data_CMS/cms/mperez/HGCal_data/L1taus/";

  vector<TString> infiles;
  //infiles.push_back(indir+Form("NTuple_ZTT_PU0_matched_0p3.root"));
  infiles.push_back(indir+Form("NTuple_ZTT_PU0_matched_0p3_BDTeg.root"));
 
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
  //histograms[0]->GetXaxis()->SetRangeUser(xmin,xmax);
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


void plot_perdecaymode (  TString var, TString mybasecut, double bins, double xmin, double xmax, TString axislabel, TString filename, bool drawnorm, float mymax){

  //github.com/PFCal-dev/cmssw/blob/hgc-tpg-devel-CMSSW_10_5_0_pre1/L1Trigger/L1THGCalUtilities/plugins/ntuples/HGCalTriggerNtupleGenTau.cc#L315-L327
  /* 
  0->1prong
  1->1prong+pi0
  4->3prongs
  5->3prongs+pi0
  */

  //TString basecut = "1";
  TString basecut = mybasecut;

  //TString indir = "/data_CMS/cms/mperez/HGCal_data/L1taus/";
  TString indir = "/data_CMS/cms/mperez/HGCal_data/Aug19/";

  vector<TString> infiles;
  //infiles.push_back(indir+Form("NTuple_ZTT_PU0_matched_0p3_BDTeg.root"));
  infiles.push_back(indir+Form("ntuple_RelValDiTau_Pt20To100_Etam1p6Tom2p9_8k_clustered.root"));
 
  vector<TH1F*> histograms;
  /*histograms.push_back(single_plot(infiles,"MatchedTree",var,basecut+" && gentau_decayMode==0", bins,xmin,xmax));
  histograms.push_back(single_plot(infiles,"MatchedTree",var,basecut+" && gentau_decayMode==1", bins,xmin,xmax));
  histograms.push_back(single_plot(infiles,"MatchedTree",var,basecut+" && gentau_decayMode==4", bins,xmin,xmax));
  histograms.push_back(single_plot(infiles,"MatchedTree",var,basecut+" && gentau_decayMode==5", bins,xmin,xmax));*/

  /*histograms.push_back(single_plot(infiles,"MatchedTree",var,basecut+" && cl3d_decayMode==0", bins,xmin,xmax));
  histograms.push_back(single_plot(infiles,"MatchedTree",var,basecut+" && cl3d_decayMode==1", bins,xmin,xmax));
  histograms.push_back(single_plot(infiles,"MatchedTree",var,basecut+" && cl3d_decayMode==4", bins,xmin,xmax));
  histograms.push_back(single_plot(infiles,"MatchedTree",var,basecut+" && cl3d_decayMode==5", bins,xmin,xmax));*/

  histograms.push_back(single_plot(infiles,"ClusteredTree",var,basecut+" && gentau_decayMode==0", bins,xmin,xmax));
  histograms.push_back(single_plot(infiles,"ClusteredTree",var,basecut+" && gentau_decayMode==1", bins,xmin,xmax));
  histograms.push_back(single_plot(infiles,"ClusteredTree",var,basecut+" && gentau_decayMode==4", bins,xmin,xmax));
  histograms.push_back(single_plot(infiles,"ClusteredTree",var,basecut+" && gentau_decayMode==5", bins,xmin,xmax));

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
  //histograms[0]->GetXaxis()->SetRangeUser(xmin,xmax);
  histograms[0]->GetYaxis()->SetRangeUser(0,ymax);
  histograms[0]->GetYaxis()->SetTitleOffset(1.7);
  histograms[0]->GetXaxis()->SetTitleOffset(1.2);

  TLegend* leg=new TLegend(0.68,0.65,0.85,0.85);
  //TLegend* leg=new TLegend(0.18,0.65,0.35,0.85);
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
  //tex.DrawLatexNDC(0.62,0.91,"#scale[1.2]{#font[40]{Fall17 Z#rightarrow#tau#tau, PU = 0}}");
  tex.DrawLatexNDC(0.55,0.91,"#scale[1.2]{#font[40]{Di-#tau, v10 geometry, PU = 0}}");

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



void plot_2D (  TString histo, int DM, TString xlabel, TString ylabel, TString filename, float xmin, float xmax, float ymin, float ymax){

  //github.com/PFCal-dev/cmssw/blob/hgc-tpg-devel-CMSSW_10_5_0_pre1/L1Trigger/L1THGCalUtilities/plugins/ntuples/HGCalTriggerNtupleGenTau.cc#L315-L327
  /* 
  0->1prong
  1->1prong+pi0
  4->3prongs
  5->3prongs+pi0
  */

  TFile f("/data_CMS/cms/mperez/HGCal_data/L1taus/NTuple_ZTT_PU0_sorted_0p3_BDTeg.root","READ");
 
  TH2F* histogram = (TH2F*)f.Get(histo);

  TCanvas* c=new TCanvas("c","c",600,600);
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.13);
  c->SetGridx();
  c->SetGridy();

  gStyle->SetOptStat(0);  

  histogram->SetTitle(" ");
  histogram->GetXaxis()->SetTitle(xlabel);
  histogram->GetYaxis()->SetTitle(ylabel);
  histogram->GetXaxis()->SetRangeUser(xmin,xmax);
  histogram->GetYaxis()->SetRangeUser(ymin,ymax);
  histogram->GetXaxis()->SetTitleOffset(1.2);
  histogram->GetYaxis()->SetTitleOffset(1.5);
  //histogram->GetXaxis()->CenterLabels(true);

  histogram->Draw("colz");

  TString myDM = "1";
  if(DM==0) myDM = "#scale[1.4]{#font[40]{1prong with #geq 2 cl3d's}}";
  else if(DM==1) myDM = "#scale[1.4]{#font[40]{1prong+#pi^{0} with #geq 2 cl3d's}}";
  else if(DM==4) myDM = "#scale[1.4]{#font[40]{3prongs with  #geq 2 cl3d's}}";
  else if(DM==5) myDM = "#scale[1.4]{#font[40]{3prongs+#pi^{0} with #geq 2 cl3d's}}";

  /*if(DM==0) myDM = "#scale[1.4]{#font[40]{1prong}}";
  else if(DM==1) myDM = "#scale[1.4]{#font[40]{1prong+#pi^{0}}}";
  else if(DM==4) myDM = "#scale[1.4]{#font[40]{3prongs}}";
  else if(DM==5) myDM = "#scale[1.4]{#font[40]{3prongs+#pi^{0}}}";*/

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.13,0.91,"#scale[1.2]{HGCal #bf{#it{Simulation}}}");
  tex.DrawLatexNDC(0.59,0.91,"#scale[1.2]{#font[40]{Fall17 Z#rightarrow#tau#tau, PU = 0}}");
  tex.DrawLatexNDC(0.4,0.83,myDM);
  //tex.DrawLatexNDC(0.65,0.78,myDM);


  c->SaveAs("../plots/pdf/"+filename+".pdf");
  c->SaveAs("../plots/png/"+filename+".png");

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

  /*plot_perdecaymode("dR_cl3d_gentau",100,0,1,"dR(gentau,cl3d)","dR_cl3d_gentau",false, 10);
  plot_perdecaymode("gentau_isMatched",2,0,2,"is Matched gentau (dR=0.3)","gentau_isMatched",false, 10);
  plot_perdecaymode("gentau_nMatchedcl3d",10,0,10,"number of matched cl3d per gentau","gentau_nMatchedcl3d",false, 10);
  plot_perdecaymode("gentau_PtMaxMatchedcl3d/gentau_vis_pt",50,0,2,"p^{T}_{max}(3dcl) / p^{T}(gentau_vis)","res_ptmax",true, 500);
  plot_perdecaymode("gentau_PtTotMatchedcl3d/gentau_vis_pt",10,0,10,"p^{T}_{tot}(3dcl)/p^{T}(gentau)","res_pttot",false, 10);*/

  /*plot_perdecaymode("cl3d_bdteg","cl3d_isMatched"50,-1,0.5,"3D cluster EG BDT","cl3d_bdteg",true);
  plot_perdecaymode("cl3d_bdteg","cl3d_isMatched && cl3d_isMaxPtcl3d", 50,-1,0.5,"max. p_{T} 3D cluster EG BDT","cl3d_bdteg_maxPt",true,300);
  plot_perdecaymode("cl3d_dRtoMaxPtcl3d","cl3d_isMatched",30,0,0.6,"dR(cl3d,max p_{T} cl3d)","cl3d_dRmaxPt3dcl",true, 300);
  plot_perdecaymode("cl3d_showerlength","cl3d_isMatched",30,0,30,"3D cluster shower length","cl3d_showerlength",true, 300);*/

  //plot_perdecaymode("gentau_isMatched","1",2,0,2,"Gen. #tau is matched", "gentau_ismatched",true,2000);
  //plot_perdecaymode("gentau_matchedSC_n_cl3d","1",5,0,5,"Number of cl3d in matched Scl3d", "gentau_ncl3d",true,2000);
  //plot_perdecaymode("gentau_matchedSC_pt_tot/gentau_vis_pt","gentau_isMatched",40,0,1.6,"p_{T}(Scl3d) / p_{T}(gen. #tau)","res_pttot",true,250);
  //plot_perdecaymode("gentau_matchedSC_pt_seed/gentau_vis_pt","gentau_isMatched",40,0,1.6,"p_{T}(seed Scl3d) / p_{T}(gen. #tau)","res_ptseed",true,200);

  //plot_perdecaymode("gentau_matchedSC_showerlength_seed","gentau_isMatched", 60,0,60, "Shower length of Scl3d seed","showerlength_seed",true,150);
  //plot_perdecaymode("gentau_matchedSC_coreshowerlength_seed","gentau_isMatched", 40,0,40, "Core shower length of Scl3d seed","coreshowerlength_seed",true,250);
  //plot_perdecaymode("gentau_matchedSC_firstlayer_seed","gentau_isMatched", 50,0,50, "First layer of Scl3d seed","firstlayer_seed",true,1200);
  //plot_perdecaymode("gentau_matchedSC_maxlayer_seed","gentau_isMatched", 50,0,50, "Max. layer of Scl3d seed","maxlayer_seed",true,500);

  plot_perdecaymode("gentau_matchedSC_seetot_seed","gentau_isMatched", 50,0,0.1, "Scl3d seed #sigma_{#eta#eta}^{tot}","seetot_seed",true,200);
  plot_perdecaymode("gentau_matchedSC_seemax_seed","gentau_isMatched", 50,0,0.15,"Scl3d seed #sigma_{#eta#eta}^{max}","seemax_seed",true,150);
  plot_perdecaymode("gentau_matchedSC_spptot_seed","gentau_isMatched", 50,0,0.1, "Scl3d seed #sigma_{#phi#phi}^{tot}","spptot_seed",true,200);
  plot_perdecaymode("gentau_matchedSC_sppmax_seed","gentau_isMatched", 50,0,0.15,"Scl3d seed #sigma_{#phi#phi}^{max}","sppmax_seed",true,150);
  /*plot_perdecaymode("gentau_matchedSC_szz_seed","gentau_isMatched",
  plot_perdecaymode("gentau_matchedSC_srrtot_seed","gentau_isMatched",
  plot_perdecaymode("gentau_matchedSC_srrmax_seed","gentau_isMatched",
  plot_perdecaymode("gentau_matchedSC_srrmean_seed","gentau_isMatched",
  plot_perdecaymode("gentau_matchedSC_emaxe_seed","gentau_isMatched",
  plot_perdecaymode("gentau_matchedSC_hoe_seed","gentau_isMatched",
  plot_perdecaymode("gentau_matchedSC_meanz_seed","gentau_isMatched",
  plot_perdecaymode("gentau_matchedSC_layer10_seed","gentau_isMatched",
  plot_perdecaymode("gentau_matchedSC_layer50_seed","gentau_isMatched",
  plot_perdecaymode("gentau_matchedSC_layer90_seed","gentau_isMatched",
  plot_perdecaymode("gentau_matchedSC_ntc67_seed","gentau_isMatched",
  plot_perdecaymode("gentau_matchedSC_ntc90_seed","gentau_isMatched",
  plot_perdecaymode("gentau_matchedSC_bdteg_seed","gentau_isMatched",
  plot_perdecaymode("gentau_matchedSC_quality_seed","gentau_isMatched",*/

}

void plot_2D_all(){

  /*plot_2D ("h2_DM0_icl3d_pticl3d_o_ptMaxcl3d",0,"i^{th} cl3d (p_{T}-sorted)","p_{T}(i^{th} cl3d) / p_{T}(max. p_{T} cl3d)","h2_ptrel_DM0");
  plot_2D ("h2_DM1_icl3d_pticl3d_o_ptMaxcl3d",1,"i^{th} cl3d (p_{T}-sorted)","p_{T}(i^{th} cl3d) / p_{T}(max. p_{T} cl3d)","h2_ptrel_DM1");
  plot_2D ("h2_DM4_icl3d_pticl3d_o_ptMaxcl3d",4,"i^{th} cl3d (p_{T}-sorted)","p_{T}(i^{th} cl3d) / p_{T}(max. p_{T} cl3d)","h2_ptrel_DM4");
  plot_2D ("h2_DM5_icl3d_pticl3d_o_ptMaxcl3d",5,"i^{th} cl3d (p_{T}-sorted)","p_{T}(i^{th} cl3d) / p_{T}(max. p_{T} cl3d)","h2_ptrel_DM5");*/

  /*plot_2D ("h2_DM0_icl3d_pticl3d",0,"i^{th} cl3d (p_{T}-sorted)","p_{T}(i^{th} cl3d) [GeV]","h2_pt_DM0",1,11,0,110);
  plot_2D ("h2_DM1_icl3d_pticl3d",1,"i^{th} cl3d (p_{T}-sorted)","p_{T}(i^{th} cl3d) [GeV]","h2_pt_DM1",1,11,0,110);
  plot_2D ("h2_DM4_icl3d_pticl3d",4,"i^{th} cl3d (p_{T}-sorted)","p_{T}(i^{th} cl3d) [GeV]","h2_pt_DM4",1,11,0,110);
  plot_2D ("h2_DM5_icl3d_pticl3d",5,"i^{th} cl3d (p_{T}-sorted)","p_{T}(i^{th} cl3d) [GeV]","h2_pt_DM5",1,11,0,110);*/

  /*plot_2D ("h2_DM0_bdtcl3d1_bdtcl3d2",0,"leading cl3d BDT EG","subleading cl3d BDT EG","h2_bdt1_bdt2_DM0",-1,0.5,-1,0.5);
  plot_2D ("h2_DM1_bdtcl3d1_bdtcl3d2",1,"leading cl3d BDT EG","subleading cl3d BDT EG","h2_bdt1_bdt2_DM1",-1,0.5,-1,0.5);
  plot_2D ("h2_DM4_bdtcl3d1_bdtcl3d2",4,"leading cl3d BDT EG","subleading cl3d BDT EG","h2_bdt1_bdt2_DM4",-1,0.5,-1,0.5);
  plot_2D ("h2_DM5_bdtcl3d1_bdtcl3d2",5,"leading cl3d BDT EG","subleading cl3d BDT EG","h2_bdt1_bdt2_DM5",-1,0.5,-1,0.5);*/

  /*plot_2D ("h2_DM0_icl3d_bdtegicl3d",0,"i^{th} cl3d (BDT EG score-sorted)","BDT EG score (i^{th} cl3d)","h2_bdtscore_DM0",1,11,-1,0.5);
  plot_2D ("h2_DM1_icl3d_bdtegicl3d",1,"i^{th} cl3d (BDT EG score-sorted)","BDT EG score (i^{th} cl3d)","h2_bdtscore_DM1",1,11,-1,0.5);
  plot_2D ("h2_DM4_icl3d_bdtegicl3d",4,"i^{th} cl3d (BDT EG score-sorted)","BDT EG score (i^{th} cl3d)","h2_bdtscore_DM4",1,11,-1,0.5);
  plot_2D ("h2_DM5_icl3d_bdtegicl3d",5,"i^{th} cl3d (BDT EG score-sorted)","BDT EG score (i^{th} cl3d)","h2_bdtscore_DM5",1,11,-1,0.5);*/



}