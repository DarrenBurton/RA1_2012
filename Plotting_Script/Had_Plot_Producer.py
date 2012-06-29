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
Category - "","OneMuon","DiMuon" Used in string search when producing webpage
WebBinning - Bins that you want to show information for on website.
'''

settings = {
  "dirs":["275_325","325_375","375_475","475_575","575_675","675_775","775_875","875",],
  "Plots":["AlphaT_all","HT_after_alphaT_55_all","Btag_Post_AlphaT_5_55_all","EffectiveMass_after_alphaT_55_all", "JetMultiplicityAfterAlphaT_55_all","MHT_after_alphaT_55_all","Number_Primary_verticies_after_alphaT_55_all","AlphaT_Zoomed_all","MHTOvMET_PF_TypeI_afterAlphaT_55_all" ,"MHT_after_alphaT_55_all","MET_PF_TypeI_afterAlphaT_55_all"],
  #"Plots":[ "AlphaT_all","HT_after_alphaT_55_all","CaloMET_TypeI_afterAlphaT_55_all","CaloMET_Uncorrected_afterAlphaT_55_all", "MET_PF_afterAlphaT_55_all", "MET_PF_TypeI_afterAlphaT_55_all","MHTOvCaloMET_TypeI_afterAlphaT_55_all", "MHTOvCaloMET_Uncorrected_afterAlphaT_55_all", "MHTOvMET_PF_afterAlphaT_55_all" , "MHTOvMET_PF_TypeI_afterAlphaT_55_all"     ,"MHT_after_alphaT_55_all"  ],
  "Lumo" : 49.80,
  "Webpage":"btag",
  "Category":"Had",
  "WebBinning":["275_325","325_375","375_upwards"],
  #"Misc":["Normalise"]
  "Misc":[],
  "Trigger":{"275":0.91,"325":0.99,"375":0.99,"475":0.99,"575":1.0,"675":1.0,"775":1.0,"875":1.0}
      }


'''
Sample Instructions

1st argument - Path to root File
2nd argument - Prefix to ht bin, "","OneMuon_","DiMuon_"
3rd argument - MC Type, ( WJets,TTbar,Zinv,DY,Di-Boson,QCD,Single_Top)
4th argument - Sample Type, "Had","Muon","DiMuon"
5th argument - Btag type, "Inclusive"(Baseline),"One","Two" etc

'''

muon_plots = {
     "nMuon":("./June26_5fb_Capped5/Had_Data.root","","Data","Had","Inclusive"), 
     "mc1":("./June26_5fb_Capped5/Had_MC.root","","MC Combined","Had","Inclusive"),
     "mc2":("./June26_5fb_Capped5/Had_WJets.root","","WJets","Had","Inclusive"),
     "mc3":("./June26_5fb_Capped5/Had_TTbar.root","","TTbar","Had","Inclusive"),
     "mc4":("./June26_5fb_Capped5/Had_Zinv.root","","Zinv","Had","Inclusive"),
     "mc5":("./June26_5fb_Capped5/Had_DY.root","","DY","Had","Inclusive"),
     "mc7":("./June26_5fb_Capped5/Had_DiBoson.root","","Di-Boson","Had","Inclusive"),
     #"mc8":("./June26_5fb_Capped5/Had_QCD.root","","QCD","Had","Inclusive"), 
     "mc9":("./June26_5fb_Capped5/Had_SingleTop.root","","Single_Top","Had","Inclusive"),
    
    }

muon_one_btag_plots = {
     "nbMuon":("./June26_5fb_Capped5/Had_Data.root","btag_one_","Data","Had","One"), 
     #"mcb1":("./June26_5fb_Capped5/Had_MC.root","btag_one_","MC Combined","Had","One"),
     "mcb2":("./June26_5fb_Capped5/Had_WJets.root","btag_one_","WJets","Had","One"),
     "mcb3":("./June26_5fb_Capped5/Had_TTbar.root","btag_one_","TTbar","Had","One"),
     "mcb4":("./June26_5fb_Capped5/Had_Zinv.root","btag_one_","Zinv","Had","One"),
     "mcb5":("./June26_5fb_Capped5/Had_DY.root","btag_one_","DY","Had","One"),
     "mcb6":("./June26_5fb_Capped5/Had_SingleTop.root","btag_one_","Single_Top","Had","One"),
     "mcb7":("./June26_5fb_Capped5/Had_DiBoson.root","btag_one_","Di-Boson","Had","One"),
     #"mcb8":("./June26_5fb_Capped5/Had_QCD.root","btag_one_","QCD","Had","One"),
        
    }


muon_two_btag_plots = {
     "nbMuon":("./June26_5fb_Capped5/Had_Data.root","btag_two_","Data","Had","Two"), 
     "mcb2":("./June26_5fb_Capped5/Had_WJets.root","btag_two_","WJets","Had","Two"),
     "mcb3":("./June26_5fb_Capped5/Had_TTbar.root","btag_two_","TTbar","Had","Two"),
     "mcb4":("./June26_5fb_Capped5/Had_Zinv.root","btag_two_","Zinv","Had","Two"),
     "mcb5":("./June26_5fb_Capped5/Had_DY.root","btag_two_","DY","Had","Two"),
     "mcb6":("./June26_5fb_Capped5/Had_SingleTop.root","btag_two_","Single_Top","Had","Two"),
     "mcb7":("./June26_5fb_Capped5/Had_DiBoson.root","btag_two_","Di-Boson","Had","Two"), 
     #"mcb8":("./June26_5fb_Capped5/Had_QCD.root","btag_two_","QCD","Had","Two"),   
    }


muon_zero_btag_plots = {
     "nbMuon":("./June26_5fb_Capped5/Had_Data.root","btag_zero_","Data","Had","Zero"), 
     "mcb2":("./June26_5fb_Capped5/Had_WJets.root","btag_zero_","WJets","Had","Zero"),
     "mcb3":("./June26_5fb_Capped5/Had_TTbar.root","btag_zero_","TTbar","Had","Zero"),
     "mcb4":("./June26_5fb_Capped5/Had_Zinv.root","btag_zero_","Zinv","Had","Zero"),
     "mcb5":("./June26_5fb_Capped5/Had_DY.root","btag_zero_","DY","Had","Zero"),
     "mcb6":("./June26_5fb_Capped5/Had_SingleTop.root","btag_zero_","Single_Top","Had","Zero"),
     "mcb7":("./June26_5fb_Capped5/Had_DiBoson.root","btag_zero_","Di-Boson","Had","Zero"),
     #"mcb8":("./June26_5fb_Capped5/Had_QCD.root","btag_zero_","QCD","Had","Zero"),
        
    }

muon_morethanzero_btag_plots = {
     "nbMuon":("./June26_5fb_Capped5/Had_Data.root","btag_morethanzero_","Data","Had","Zero"), 
     "mcb2":("./June26_5fb_Capped5/Had_WJets.root","btag_morethanzero_","WJets","Had","Zero"),
     "mcb3":("./June26_5fb_Capped5/Had_TTbar.root","btag_morethanzero_","TTbar","Had","Zero"),
     "mcb4":("./June26_5fb_Capped5/Had_Zinv.root","btag_morethanzero_","Zinv","Had","Zero"),
     "mcb5":("./June26_5fb_Capped5/Had_DY.root","btag_morethanzero_","DY","Had","Zero"),
     "mcb7":("./June26_5fb_Capped5/Had_DiBoson.root","btag_morethanzero_","Di-Boson","Had","Zero"),
     #"mcb8":("./June26_5fb_Capped5/Had_QCD.root","btag_morethanzero_","QCD","Had","Zero"), 
     "mcb9":("./June26_5fb_Capped5/Had_SingleTop.root","btag_morethanzero_","Single_Top","Had","Zero"), 
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
  b = Plotter(settings,muon_morethanzero_btag_plots,jet_multiplicity = "False")
  c = Plotter(settings,muon_two_btag_plots,jet_multiplicity = "False")
  d = Plotter(settings,muon_zero_btag_plots,jet_multiplicity = "False")
  e = Plotter(settings,muon_one_btag_plots,jet_multiplicity = "False")


  finish = Webpage_Maker(settings["Plots"],settings["WebBinning"],settings["Category"],option=settings["Webpage"])
