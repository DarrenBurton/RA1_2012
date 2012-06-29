#!/usr/bin/env python
from ROOT import *
import ROOT as r
import logging,itertools
import os,fnmatch,sys
import glob, errno
from time import strftime
from optparse import OptionParser
import array, ast
import math as m
#from plottingUtils import *
from Btag_8TeV_Plots import *


'''
Setting Intrustions

dirs - HT bins
Plots - Plots to produce Histograms for
Lumo - Scale Factor for MC
----- Webpage arguments
Webpage - Webpage mode, just keep as btag
Category - "Photon_","OneMuon","DiMuon" Used in string search when producing webpage
WebBinning - Bins that you want to show information for on website.
'''

settings = {
  "dirs":["275_325","325_375","375_475","475_575","575_675","675_775","775_875","875",],
  #"AlphaTSlices":["0.52_0.53","0.53_0.55"],
  "Plots":["Muon_Iso_PU_Interactions_all"],
  #"Plots":["AlphaT_all","HT_after_alphaT_55_all","Btag_Post_AlphaT_5_55_all","EffectiveMass_after_alphaT_55_all", "JetMultiplicityAfterAlphaT_55_all","MHT_after_alphaT_55_all","Number_Primary_verticies_after_alphaT_55_all","AlphaT_Zoomed_all","MHTOvMET_PF_TypeI_afterAlphaT_55_all" ,"MHT_after_alphaT_55_all","MET_PF_TypeI_afterAlphaT_55_all","PhotonIso_all","Photonrho25_all","PhotonPt_all"],
  #"Plots":[ "AlphaT_all","HT_after_alphaT_55_all","CaloMET_TypeI_afterAlphaT_55_all","CaloMET_Uncorrected_afterAlphaT_55_all", "MET_PF_afterAlphaT_55_all", "MET_PF_TypeI_afterAlphaT_55_all","MHTOvCaloMET_TypeI_afterAlphaT_55_all", "MHTOvCaloMET_Uncorrected_afterAlphaT_55_all", "MHTOvMET_PF_afterAlphaT_55_all" , "MHTOvMET_PF_TypeI_afterAlphaT_55_all"     ,"MHT_after_alphaT_55_all"  ],
  "Lumo" : 1.0,
  "Webpage":"btag",
  "Category":"OneMuon",
  "WebBinning":["375_upwards"],
  #"Misc":["Normalise"]
  "Misc":[],
  "Trigger":{"275":1.0,"325":1.0,"375":1.0,"475":1.0,"575":1.0,"675":1.0,"775":1.0,"875":1.0}
  #"WebBinning":["275_325","325_375","375_475"]
      }


'''
Sample Instructions

1st argument - Path to root File
2nd argument - Prefix to ht bin, "Photon_","OneMuon_","DiMuon_"
3rd argument - MC Type, ( WJets,TTbar,Zinv,DY,Di-Boson,QCD,Single_Top)
4th argument - Sample Type, "Photon","Muon","DiMuon"
5th argument - Btag type, "Inclusive"(Baseline),"One","Two" etc

'''

muon_plots = {
     "nMuon":("./NoIso_Muon/Muon_Data.root","OneMuon_","Data","Muon","Inclusive"), 
     "mc9":("./NoIso_Muon/Muon_MC.root","OneMuon_","TTbar","Muon","Inclusive"),
    
    }


'''

Plotter Instructions

Imported from Btag_plots.py
Plotter will produce all plots specified in settings["Plots"]. Option file for rebinning etc found in Btag_Plots. 

Additionally
-----------
Currently has a hack to correct for Zinv NLO scaling factor. We have it as k-factor 1.28 where as WJets etc has 1.12.
Rob suggested for plots to scale Zinv by 1.12/1.28. This is done in MC_Plot function in Btag_Plots.
Also to account for this in the Combined MC. I have written a function to add together all the MC contributions within
the code rather than using the MC_Combined root file. I think this will only work for TH1D's but they are the only plots
within my plotting op. 

'''



if __name__=="__main__":
  a = Plotter(settings,muon_plots,jet_multiplicity = "False")
  #b = Plotter(settings,muon_morethanzero_btag_plots,jet_multiplicity = "False")
  #c = Plotter(settings,muon_two_btag_plots,jet_multiplicity = "False")
  #d = Plotter(settings,muon_zero_btag_plots,jet_multiplicity = "False")
  #e = Plotter(settings,muon_one_btag_plots,jet_multiplicity = "False")


  finish = Webpage_Maker(settings["Plots"],settings["WebBinning"],settings["Category"],option=settings["Webpage"])
