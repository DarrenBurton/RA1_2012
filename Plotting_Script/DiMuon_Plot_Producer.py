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
  "Plots":["EffectiveMass_all","MHT_all","AlphaT_all","AlphaT_Zoomed_all","HT_all","JetMultiplicity_all","Btag_Pre_AlphaT_5__all","Number_Primary_verticies_all", "MHTOvMET_PF_TypeI_all", "MuPt_all","MuPFIso_all"],
  "Lumo" : 49.80,
  "Webpage":"btag",
  "Category":"DiMuon",
  "WebBinning":["275_325","325_375","375_upwards"],
  "Misc":[],
  "Trigger":{"275":0.95,"325":0.96,"375":0.96,"475":0.97,"575":0.97,"675":0.97,"775":0.98,"875":0.98}
      }

muon_plots = {
     "nMuon":("./June26_5fb_Capped5/Muon_Data.root","DiMuon_","Data","DiMuon","Inclusive"), 
     "mc2":("./June26_5fb_Capped5/Muon_WJets.root","DiMuon_","WJets","DiMuon","Inclusive"),
     "mc3":("./June26_5fb_Capped5/Muon_TTbar.root","DiMuon_","TTbar","DiMuon","Inclusive"),
     "mc4":("./June26_5fb_Capped5/Muon_Zinv.root","DiMuon_","Zinv","DiMuon","Inclusive"),
     "mc5":("./June26_5fb_Capped5/Muon_DY.root","DiMuon_","DY","DiMuon","Inclusive"),
     "mc6":("./June26_5fb_Capped5/Muon_SingleTop.root","DiMuon_","Single_Top","DiMuon","Inclusive"),
     "mc7":("./June26_5fb_Capped5/Muon_DiBoson.root","DiMuon_","Di-Boson","DiMuon","Inclusive"),
     #"mc8":("./June26_5fb_Capped5/Muon_QCD.root","DiMuon_","QCD","DiMuon","Inclusive"),
    
    
    }

muon_one_btag_plots = {
     "nbMuon":("./June26_5fb_Capped5/Muon_Data.root","btag_one_DiMuon_","Data","DiMuon","One"), 
     "mcb1":("./June26_5fb_Capped5/Muon_MC.root","btag_one_DiMuon_","MC Combined","DiMuon","One"),
     "mcb2":("./June26_5fb_Capped5/Muon_WJets.root","btag_one_DiMuon_","WJets","DiMuon","One"),
     "mcb3":("./June26_5fb_Capped5/Muon_TTbar.root","btag_one_DiMuon_","TTbar","DiMuon","One"),
     "mcb4":("./June26_5fb_Capped5/Muon_Zinv.root","btag_one_DiMuon_","Zinv","DiMuon","One"),
     "mcb5":("./June26_5fb_Capped5/Muon_DY.root","btag_one_DiMuon_","DY","DiMuon","One"),
     "mcb6":("./June26_5fb_Capped5/Muon_SingleTop.root","btag_one_DiMuon_","Single_Top","DiMuon","One"),
     "mcb7":("./June26_5fb_Capped5/Muon_DiBoson.root","btag_one_DiMuon_","Di-Boson","DiMuon","One"),
     #"mcb8":("./June26_5fb_Capped5/Muon_QCD.root","btag_one_DiMuon_","QCD","DiMuon","One"),
        
    }


muon_two_btag_plots = {
     "nbMuon":("./June26_5fb_Capped5/Muon_Data.root","btag_two_DiMuon_","Data","DiMuon","Two"), 
     "mcb2":("./June26_5fb_Capped5/Muon_WJets.root","btag_two_DiMuon_","WJets","DiMuon","Two"),
     "mcb3":("./June26_5fb_Capped5/Muon_TTbar.root","btag_two_DiMuon_","TTbar","DiMuon","Two"),
     "mcb4":("./June26_5fb_Capped5/Muon_Zinv.root","btag_two_DiMuon_","Zinv","DiMuon","Two"),
     "mcb5":("./June26_5fb_Capped5/Muon_DY.root","btag_two_DiMuon_","DY","DiMuon","Two"),
     "mcb6":("./June26_5fb_Capped5/Muon_SingleTop.root","btag_two_DiMuon_","Single_Top","DiMuon","Two"),
     "mcb7":("./June26_5fb_Capped5/Muon_DiBoson.root","btag_two_DiMuon_","Di-Boson","DiMuon","Two"), 
     #"mcb8":("./June26_5fb_Capped5/Muon_QCD.root","btag_two_DiMuon_","QCD","DiMuon","Two"),   
    }


muon_zero_btag_plots = {
     "nbMuon":("./June26_5fb_Capped5/Muon_Data.root","btag_zero_DiMuon_","Data","DiMuon","Zero"), 
     "mcb2":("./June26_5fb_Capped5/Muon_WJets.root","btag_zero_DiMuon_","WJets","DiMuon","Zero"),
     "mcb3":("./June26_5fb_Capped5/Muon_TTbar.root","btag_zero_DiMuon_","TTbar","DiMuon","Zero"),
     "mcb4":("./June26_5fb_Capped5/Muon_Zinv.root","btag_zero_DiMuon_","Zinv","DiMuon","Zero"),
     "mcb5":("./June26_5fb_Capped5/Muon_DY.root","btag_zero_DiMuon_","DY","DiMuon","Zero"),
     "mcb6":("./June26_5fb_Capped5/Muon_SingleTop.root","btag_zero_DiMuon_","Single_Top","DiMuon","Zero"),
     "mcb7":("./June26_5fb_Capped5/Muon_DiBoson.root","btag_zero_DiMuon_","Di-Boson","DiMuon","Zero"),
     #"mcb8":("./June26_5fb_Capped5/Muon_QCD.root","btag_zero_DiMuon_","QCD","DiMuon","Zero"),
        
    }


muon_morethanzero_btag_plots = {

     "nbMuon":("./June26_5fb_Capped5/Muon_Data.root","btag_morethanzero_DiMuon_","Data","DiMuon","Zero"), 
     "mcb2":("./June26_5fb_Capped5/Muon_WJets.root","btag_morethanzero_DiMuon_","WJets","DiMuon","Zero"),
     "mcb3":("./June26_5fb_Capped5/Muon_TTbar.root","btag_morethanzero_DiMuon_","TTbar","DiMuon","Zero"),
     "mcb4":("./June26_5fb_Capped5/Muon_Zinv.root","btag_morethanzero_DiMuon_","Zinv","DiMuon","Zero"),
     "mcb5":("./June26_5fb_Capped5/Muon_DY.root","btag_morethanzero_DiMuon_","DY","DiMuon","Zero"),
     "mcb6":("./June26_5fb_Capped5/Muon_SingleTop.root","btag_morethanzero_DiMuon_","Single_Top","DiMuon","Zero"),
     "mcb7":("./June26_5fb_Capped5/Muon_DiBoson.root","btag_morethanzero_DiMuon_","Di-Boson","DiMuon","Zero"),
     #"mcb8":("./June26_5fb_Capped5/Muon_QCD.root","btag_morethanzero_DiMuon_","QCD","DiMuon","Zero"),
        
    }


if __name__=="__main__":
  a = Plotter(settings,muon_plots,jet_multiplicity = "False")
  b = Plotter(settings,muon_one_btag_plots,jet_multiplicity = "False")
  c = Plotter(settings,muon_two_btag_plots,jet_multiplicity = "False")
  d = Plotter(settings,muon_morethanzero_btag_plots,jet_multiplicity = "False")
  e = Plotter(settings,muon_zero_btag_plots,jet_multiplicity = "False")

  finish = Webpage_Maker(settings["Plots"],settings["WebBinning"],settings["Category"],option=settings["Webpage"])
