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


settings = {
  "dirs":["275_325","325_375","375_475","475_575","575_675","675_775","775_875","875",],
  "Plots":["MuPFIso_all","MuPFIso_Less15_all","MuPFIso_15_20_all","MuPFIso_More_20_all" ], 
  #"Plots":["MuPt_all","MuPFIso_all","MT_all","EffectiveMass_all","MHT_all","AlphaT_all","AlphaT_Zoomed_all","JetMultiplicity_all","Number_Primary_verticies_all","HT_all", "Btag_Pre_AlphaT_5__all","MHTOvMET_PF_TypeI_all","MHTOvMET_PF_TypeI_all" ],
  #"Plots":[ "MHT_all","AlphaT_all","AlphaT_Zoomed_all","JetMultiplicity_all","HT_all","Btag_Pre_AlphaT_5__all","EffectiveMass_all", "Number_Primary_verticies_all","MHTOvCaloMET_Lepton_all","MHTOvCaloMET_all" ],
  #"Plots":["MHT_all","AlphaT_all","AlphaT_Zoomed_all","JetMultiplicity_all","HT_all","Btag_Pre_AlphaT_5__all","EffectiveMass_all","MuPt_all","MuEIso_all","MuTrIso_all","MuHIso_all","MuCso_all","Number_Primary_verticies_all","MHTOvCaloMET_all","MHTOvCaloMET_Lepton_all"],
  "Lumo" : 49.80,
  "Webpage":"btag",
  "Category":"OneMuon",
  "WebBinning":["275_325","325_375","375_upwards"],
  "Misc":[],
  "Trigger":{"275":0.88,"325":0.88,"375":0.88,"475":0.88,"575":0.88,"675":0.88,"775":0.88,"875":0.88}
  }
      

muon_plots = {
     "nMuon":("./NoIso_Muon/Muon_Data.root","OneMuon_","Data","Muon","Inclusive"), 
     "mc2":("./NoIso_Muon/Muon_WJets.root","OneMuon_","WJets","Muon","Inclusive"),
     "mc3":("./NoIso_Muon/Muon_TTbar.root","OneMuon_","TTbar","Muon","Inclusive"),
     "mc4":("./NoIso_Muon/Muon_Zinv.root","OneMuon_","Zinv","Muon","Inclusive"),
     "mc5":("./NoIso_Muon/Muon_DY.root","OneMuon_","DY","Muon","Inclusive"),
     "mc7":("./NoIso_Muon/Muon_DiBoson.root","OneMuon_","Di-Boson","Muon","Inclusive"),
     #"mc8":("./NoIso_Muon/Muon_QCD.root","OneMuon_","QCD","Muon","Inclusive"), 
     "mc9":("./NoIso_Muon/Muon_SingleTop.root","OneMuon_","Single_Top","Muon","Inclusive"),
    
    }

muon_one_btag_plots = {
     "nbMuon":("./NoIso_Muon/Muon_Data.root","btag_one_OneMuon_","Data","Muon","One"), 
     "mcb2":("./NoIso_Muon/Muon_WJets.root","btag_one_OneMuon_","WJets","Muon","One"),
     "mcb3":("./NoIso_Muon/Muon_TTbar.root","btag_one_OneMuon_","TTbar","Muon","One"),
     "mcb4":("./NoIso_Muon/Muon_Zinv.root","btag_one_OneMuon_","Zinv","Muon","One"),
     "mcb5":("./NoIso_Muon/Muon_DY.root","btag_one_OneMuon_","DY","Muon","One"),
     "mcb6":("./NoIso_Muon/Muon_SingleTop.root","btag_one_OneMuon_","Single_Top","Muon","One"),
     "mcb7":("./NoIso_Muon/Muon_DiBoson.root","btag_one_OneMuon_","Di-Boson","Muon","One"),
    # "mcb8":("./NoIso_Muon/Muon_QCD.root","btag_one_OneMuon_","QCD","Muon","One"),
        
    }


muon_two_btag_plots = {
     "nbMuon":("./NoIso_Muon/Muon_Data.root","btag_two_OneMuon_","Data","Muon","Two"), 
     "mcb2":("./NoIso_Muon/Muon_WJets.root","btag_two_OneMuon_","WJets","Muon","Two"),
     "mcb3":("./NoIso_Muon/Muon_TTbar.root","btag_two_OneMuon_","TTbar","Muon","Two"),
     "mcb4":("./NoIso_Muon/Muon_Zinv.root","btag_two_OneMuon_","Zinv","Muon","Two"),
     "mcb5":("./NoIso_Muon/Muon_DY.root","btag_two_OneMuon_","DY","Muon","Two"),
     "mcb6":("./NoIso_Muon/Muon_SingleTop.root","btag_two_OneMuon_","Single_Top","Muon","Two"),
     "mcb7":("./NoIso_Muon/Muon_DiBoson.root","btag_two_OneMuon_","Di-Boson","Muon","Two"), 
     #"mcb8":("./NoIso_Muon/Muon_QCD.root","btag_two_OneMuon_","QCD","Muon","Two"),   
    }


muon_zero_btag_plots = {
     "nbMuon":("./NoIso_Muon/Muon_Data.root","btag_zero_OneMuon_","Data","Muon","Zero"), 
     "mcb2":("./NoIso_Muon/Muon_WJets.root","btag_zero_OneMuon_","WJets","Muon","Zero"),
     "mcb3":("./NoIso_Muon/Muon_TTbar.root","btag_zero_OneMuon_","TTbar","Muon","Zero"),
     "mcb4":("./NoIso_Muon/Muon_Zinv.root","btag_zero_OneMuon_","Zinv","Muon","Zero"),
     "mcb5":("./NoIso_Muon/Muon_DY.root","btag_zero_OneMuon_","DY","Muon","Zero"),
     "mcb6":("./NoIso_Muon/Muon_SingleTop.root","btag_zero_OneMuon_","Single_Top","Muon","Zero"),
     "mcb7":("./NoIso_Muon/Muon_DiBoson.root","btag_zero_OneMuon_","Di-Boson","Muon","Zero"),
     #"mcb8":("./NoIso_Muon/Muon_QCD.root","btag_zero_OneMuon_","QCD","Muon","Zero"),
        
    }


muon_morethanzero_btag_plots = {
     "nbMuon":("./NoIso_Muon/Muon_Data.root","btag_morethanzero_OneMuon_","Data","Muon","Zero"), 
     "mcb2":("./NoIso_Muon/Muon_WJets.root","btag_morethanzero_OneMuon_","WJets","Muon","Zero"),
     "mcb3":("./NoIso_Muon/Muon_TTbar.root","btag_morethanzero_OneMuon_","TTbar","Muon","Zero"),
     "mcb4":("./NoIso_Muon/Muon_Zinv.root","btag_morethanzero_OneMuon_","Zinv","Muon","Zero"),
     "mcb5":("./NoIso_Muon/Muon_DY.root","btag_morethanzero_OneMuon_","DY","Muon","Zero"),
     "mcb7":("./NoIso_Muon/Muon_DiBoson.root","btag_morethanzero_OneMuon_","Di-Boson","Muon","Zero"),
     #"mcb8":("./NoIso_Muon/Muon_QCD.root","btag_morethanzero_OneMuon_","QCD","Muon","Zero"), 
     "mcb9":("./NoIso_Muon/Muon_SingleTop.root","btag_morethanzero_OneMuon_","Single_Top","Muon","Zero"),
    }

if __name__=="__main__":
  a = Plotter(settings,muon_plots,jet_multiplicity = "False")
  b = Plotter(settings,muon_morethanzero_btag_plots,jet_multiplicity = "False")
  c = Plotter(settings,muon_two_btag_plots,jet_multiplicity = "False")
  d = Plotter(settings,muon_zero_btag_plots,jet_multiplicity = "False")
  e = Plotter(settings,muon_one_btag_plots,jet_multiplicity = "False")

  finish = Webpage_Maker(settings["Plots"],settings["WebBinning"],settings["Category"],option=settings["Webpage"])
