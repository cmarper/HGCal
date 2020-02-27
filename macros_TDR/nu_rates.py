#!/usr/bin/env python 
from ROOT import *

indir = '/data_CMS_upgrade/mperez/HGCal_data/Sep19/rates/'
outdir = '/data_CMS_upgrade/mperez/HGCal_data/Sep19/rates/'

f = TFile( indir+'ntuple_Nu_E10_v10_PU200_rates.root ')
fileout=outdir+'rates_Nu_E10_v10_PU200.root'

#f = TFile( indir+'ntuple_Nu_E10_v10_PU200_rates_DM0.root')
#fileout=outdir+'rates_Nu_E10_v10_PU200_DM0.root'

#f = TFile( indir+'ntuple_Nu_E10_v10_PU200_rates_DM1.root')
#fileout=outdir+'rates_Nu_E10_v10_PU200_DM1.root'

#f = TFile( indir+'ntuple_Nu_E10_v10_PU200_rates_DM2.root')
#fileout=outdir+'rates_Nu_E10_v10_PU200_DM2.root'

#f = TFile( indir+'ntuple_Nu_E10_v10_PU200_rates.root ')
#fileout=outdir+'rates_Nu_E10_v10_PU200_DM0_norm_inclusive.root'

#f = TFile( indir+'ntuple_Nu_E10_v10_PU200_rates.root ')
#fileout=outdir+'rates_Nu_E10_v10_PU200_DM1_norm_inclusive.root'

#f = TFile( indir+'ntuple_Nu_E10_v10_PU200_rates.root ')
#fileout=outdir+'rates_Nu_E10_v10_PU200_DM2_norm_inclusive.root'

print fileout

tree = f.Get("SkimmedTree")

totalrate=2760*11246./1000  # bunches * frequency , /1000 to pass to kHz

entries=tree.GetEntries()


def fillLead(namePt,nameEta,bins,start,end,binsEta,etaMin,etaMax):
	histoPt=TH1F(namePt,namePt,bins,start,end)
	histoEta=TH1F(nameEta,nameEta,binsEta,-5,5)

	for event in tree:

		pt = getattr(event, "cl3d_pt_c3")
		eta = getattr(event, "cl3d_eta")
                dm = getattr(event, "cl3d_predDM")                
       
                #if dm != 2:
		#	continue

		if  abs(eta)<etaMin or abs(eta)>etaMax:
			continue

		histoPt.Fill(pt)

		if pt>20:  # just an example to cut out the low pt stuff 
			histoEta.Fill(eta)

	return histoPt,histoEta 


def doRate(histo,name):
	historate=histo.Clone()
	historate.SetName(name)
	for i in range(0,historate.GetNbinsX()):
		integral=historate.Integral(i,-1)
		historate.SetBinContent(i,integral)
	historate.Scale(totalrate/entries)
	return historate


# make a bunch of rates and pt/eta plots:
TauPTLead,TauETALead=fillLead("TauPTLead","TauETALead",100,0,100,100,1.6,2.9)
rateTau=doRate(TauPTLead,"RateTau")

# save all the rates...
out=TFile(fileout,"RECREATE")
out.cd()
rateTau.Write()
TauETALead.Write()
TauPTLead.Write()

