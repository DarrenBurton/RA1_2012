#!/usr/bin/env python
from ROOT import *
import ROOT as r
import logging,itertools
import os,fnmatch,sys
import glob, errno
from time import strftime
from optparse import OptionParser
import array  #, ast
import math as m
from plottingUtils import *
from Bryn_Numbers import *


settings = {
  #"dirs":["275_325"],
  "dirs":["275_325","325_375","375_475","475_575","575_675","675_775","775_875","875",],  #HT Bins
  "plots":["AlphaT_all",],  # Histogram that Yields are taken from
  "AlphaTSlices":["0.55_20"], # AlphaT Slices
  "Lumo":4.98, # Luminosity in fb
  "Multi_Lumi":{'Had':4.98,'Muon':4.963,'DiMuon':4.963,'Photon':4.988},
  "Analysis":"8TeV"
      }


'''
Sample Dictionary Instructions

eg "nMuon":("../June27_5fb_Reweighting/Muon_Data","btag_two_OneMuon_","Data","Muon"),
if n at start of name entry then the file is data and will no be scaled to luminosity.
first argument is path to root file
second argument is prefix to ht bin. i.e OneMuon_275_325
third argument is data/mc type, i.e. Data, WJets250 - MC relating to the binned WJets 250-300 HT sample
fourth argument is sample Type, Had/DiMuon/Muon. 

the only thing that will have to be changed is the second argument depending on wether you are running btag multiplicity/baseline
'''
btag_two_samples = {

    "nHad":("../June27_5fb_Reweighting/Had_Data","btag_two_","Data","Had"),
    
     "mcHadW1":("../June27_5fb_Reweighting/Had_WJets","","WJetsInc","Had"),
     "mcHadttbar":("../June27_5fb_Reweighting/Had_TTbar","","TTbar","Had"),
     "mcHadzinv":("../June27_5fb_Reweighting/Had_Zinv","","Zinv50","Had"),
     "mcHadsingt":("../June27_5fb_Reweighting/Had_SingleTop","","Single_Tbar_t","Had"),
     "mcHaddiboson":("../June27_5fb_Reweighting/Had_DiBoson","","ZZ","Had"),
     "mcHadDY":("../June27_5fb_Reweighting/Had_DY","","DY","Had"),

    "nMuon":("../June27_5fb_Reweighting/Muon_Data","btag_two_OneMuon_","Data","Muon"),
    
     "mcMuonW1":("../June27_5fb_Reweighting/Muon_WJets","OneMuon_","WJetsInc","Muon"),
     "mcMuonttbar":("../June27_5fb_Reweighting/Muon_TTbar","OneMuon_","TTbar","Muon"),
     "mcMuonzinv":("../June27_5fb_Reweighting/Muon_Zinv","OneMuon_","Zinv50","Muon"),
     "mcMuonsingt":("../June27_5fb_Reweighting/Muon_SingleTop","OneMuon_","Single_Tbar_t","Muon"),
     "mcMuondiboson":("../June27_5fb_Reweighting/Muon_DiBoson","OneMuon_","ZZ","Muon"),
     "mcMuonDY":("../June27_5fb_Reweighting/Muon_DY","OneMuon_","DY","Muon"),


    "nDiMuon":("../June27_5fb_Reweighting/Muon_Data","btag_two_DiMuon_","Data","DiMuon"),
    
     "mcDiMuonW1":("../June27_5fb_Reweighting/Muon_WJets","DiMuon_","WJetsInc","DiMuon"),
     "mcDiMuonttbar":("../June27_5fb_Reweighting/Muon_TTbar","DiMuon_","TTbar","DiMuon"),
     "mcDiMuonzinv":("../June27_5fb_Reweighting/Muon_Zinv","DiMuon_","Zinv50","DiMuon"),
     "mcDiMuonsingt":("../June27_5fb_Reweighting/Muon_SingleTop","DiMuon_","Single_Tbar_t","DiMuon"),
     "mcDiMuondiboson":("../June27_5fb_Reweighting/Muon_DiBoson","DiMuon_","ZZ","DiMuon"),
     "mcDiMuonDY":("../June27_5fb_Reweighting/Muon_DY","DiMuon_","DY","DiMuon"),

     "mcPhoton":("../June27_5fb_Reweighting/Photon_MC","Photon_","Photon","Photon"),
     "ncPhoton":("../June27_5fb_Reweighting/Photon_Data","btag_two_Photon_","Data","Photon"),


    }


btag_one_samples = {

    "nHad":("../June27_5fb_Reweighting/Had_Data","btag_one_","Data","Had"),
    
     "mcHadW1":("../June27_5fb_Reweighting/Had_WJets","","WJetsInc","Had"),
     "mcHadttbar":("../June27_5fb_Reweighting/Had_TTbar","","TTbar","Had"),
     "mcHadzinv":("../June27_5fb_Reweighting/Had_Zinv","","Zinv50","Had"),
     "mcHadsingt":("../June27_5fb_Reweighting/Had_SingleTop","","Single_Tbar_t","Had"),
     "mcHaddiboson":("../June27_5fb_Reweighting/Had_DiBoson","","ZZ","Had"),
     "mcHadDY":("../June27_5fb_Reweighting/Had_DY","","DY","Had"),


    "nMuon":("../June27_5fb_Reweighting/Muon_Data","btag_one_OneMuon_","Data","Muon"),
    
     "mcMuonW1":("../June27_5fb_Reweighting/Muon_WJets","OneMuon_","WJetsInc","Muon"),
     "mcMuonttbar":("../June27_5fb_Reweighting/Muon_TTbar","OneMuon_","TTbar","Muon"),
     "mcMuonzinv":("../June27_5fb_Reweighting/Muon_Zinv","OneMuon_","Zinv50","Muon"),
     "mcMuonsingt":("../June27_5fb_Reweighting/Muon_SingleTop","OneMuon_","Single_Tbar_t","Muon"),
     "mcMuondiboson":("../June27_5fb_Reweighting/Muon_DiBoson","OneMuon_","ZZ","Muon"),
     "mcMuonDY":("../June27_5fb_Reweighting/Muon_DY","OneMuon_","DY","Muon"),


    "nDiMuon":("../June27_5fb_Reweighting/Muon_Data","btag_one_DiMuon_","Data","DiMuon"),
    
     "mcDiMuonW1":("../June27_5fb_Reweighting/Muon_WJets","DiMuon_","WJetsInc","DiMuon"),
     "mcDiMuonttbar":("../June27_5fb_Reweighting/Muon_TTbar","DiMuon_","TTbar","DiMuon"),
     "mcDiMuonzinv":("../June27_5fb_Reweighting/Muon_Zinv","DiMuon_","Zinv50","DiMuon"),
     "mcDiMuonsingt":("../June27_5fb_Reweighting/Muon_SingleTop","DiMuon_","Single_Tbar_t","DiMuon"),
     "mcDiMuondiboson":("../June27_5fb_Reweighting/Muon_DiBoson","DiMuon_","ZZ","DiMuon"),
     "mcDiMuonDY":("../June27_5fb_Reweighting/Muon_DY","DiMuon_","DY","DiMuon"),

     "mcPhoton":("../June27_5fb_Reweighting/Photon_MC","Photon_","Photon","Photon"),
     "ncPhoton":("../June27_5fb_Reweighting/Photon_Data","btag_one_Photon_","Data","Photon"),


    }



btag_zero_samples = {

    "nHad":("../June27_5fb_Reweighting/Had_Data","btag_zero_","Data","Had"),
    
     "mcHadW1":("../June27_5fb_Reweighting/Had_WJets","","WJetsInc","Had"),
     "mcHadttbar":("../June27_5fb_Reweighting/Had_TTbar","","TTbar","Had"),
     "mcHadzinv":("../June27_5fb_Reweighting/Had_Zinv","","Zinv50","Had"),
     "mcHadsingt":("../June27_5fb_Reweighting/Had_SingleTop","","Single_Tbar_t","Had"),
     "mcHaddiboson":("../June27_5fb_Reweighting/Had_DiBoson","","ZZ","Had"),
     "mcHadDY":("../June27_5fb_Reweighting/Had_DY","","DY","Had"),


    "nMuon":("../June27_5fb_Reweighting/Muon_Data","btag_zero_OneMuon_","Data","Muon"),
    
     "mcMuonW1":("../June27_5fb_Reweighting/Muon_WJets","OneMuon_","WJetsInc","Muon"),
     "mcMuonttbar":("../June27_5fb_Reweighting/Muon_TTbar","OneMuon_","TTbar","Muon"),
     "mcMuonzinv":("../June27_5fb_Reweighting/Muon_Zinv","OneMuon_","Zinv50","Muon"),
     "mcMuonsingt":("../June27_5fb_Reweighting/Muon_SingleTop","OneMuon_","Single_Tbar_t","Muon"),
     "mcMuondiboson":("../June27_5fb_Reweighting/Muon_DiBoson","OneMuon_","ZZ","Muon"),
     "mcMuonDY":("../June27_5fb_Reweighting/Muon_DY","OneMuon_","DY","Muon"),


    "nDiMuon":("../June27_5fb_Reweighting/Muon_Data","btag_zero_DiMuon_","Data","DiMuon"),
    
     "mcDiMuonW1":("../June27_5fb_Reweighting/Muon_WJets","DiMuon_","WJetsInc","DiMuon"),
     "mcDiMuonttbar":("../June27_5fb_Reweighting/Muon_TTbar","DiMuon_","TTbar","DiMuon"),
     "mcDiMuonzinv":("../June27_5fb_Reweighting/Muon_Zinv","DiMuon_","Zinv50","DiMuon"),
     "mcDiMuonsingt":("../June27_5fb_Reweighting/Muon_SingleTop","DiMuon_","Single_Tbar_t","DiMuon"),
     "mcDiMuondiboson":("../June27_5fb_Reweighting/Muon_DiBoson","DiMuon_","ZZ","DiMuon"),
     "mcDiMuonDY":("../June27_5fb_Reweighting/Muon_DY","DiMuon_","DY","DiMuon"),

     "mcPhoton":("../June27_5fb_Reweighting/Photon_MC","Photon_","Photon","Photon"),
     "ncPhoton":("../June27_5fb_Reweighting/Photon_Data","btag_zero_Photon_","Data","Photon"),

    }


btag_more_than_two_samples = {

    "nHad":("../June27_5fb_Reweighting/Had_Data","btag_morethantwo_","Data","Had"),
    
     "mcHadW1":("../June27_5fb_Reweighting/Had_WJets","","WJetsInc","Had"),
     "mcHadttbar":("../June27_5fb_Reweighting/Had_TTbar","","TTbar","Had"),
     "mcHadzinv":("../June27_5fb_Reweighting/Had_Zinv","","Zinv50","Had"),
     "mcHadsingt":("../June27_5fb_Reweighting/Had_SingleTop","","Single_Tbar_t","Had"),
     "mcHaddiboson":("../June27_5fb_Reweighting/Had_DiBoson","","ZZ","Had"),
     "mcHadDY":("../June27_5fb_Reweighting/Had_DY","","DY","Had"),


    "nMuon":("../June27_5fb_Reweighting/Muon_Data","btag_morethantwo_OneMuon_","Data","Muon"),
    
     "mcMuonW1":("../June27_5fb_Reweighting/Muon_WJets","OneMuon_","WJetsInc","Muon"),
     "mcMuonttbar":("../June27_5fb_Reweighting/Muon_TTbar","OneMuon_","TTbar","Muon"),
     "mcMuonzinv":("../June27_5fb_Reweighting/Muon_Zinv","OneMuon_","Zinv50","Muon"),
     "mcMuonsingt":("../June27_5fb_Reweighting/Muon_SingleTop","OneMuon_","Single_Tbar_t","Muon"),
     "mcMuondiboson":("../June27_5fb_Reweighting/Muon_DiBoson","OneMuon_","ZZ","Muon"),
     "mcMuonDY":("../June27_5fb_Reweighting/Muon_DY","OneMuon_","DY","Muon"),


    "nDiMuon":("../June27_5fb_Reweighting/Muon_Data","btag_morethantwo_DiMuon_","Data","DiMuon"),
    
     "mcDiMuonW1":("../June27_5fb_Reweighting/Muon_WJets","DiMuon_","WJetsInc","DiMuon"),
     "mcDiMuonttbar":("../June27_5fb_Reweighting/Muon_TTbar","DiMuon_","TTbar","DiMuon"),
     "mcDiMuonzinv":("../June27_5fb_Reweighting/Muon_Zinv","DiMuon_","Zinv50","DiMuon"),
     "mcDiMuonsingt":("../June27_5fb_Reweighting/Muon_SingleTop","DiMuon_","Single_Tbar_t","DiMuon"),
     "mcDiMuondiboson":("../June27_5fb_Reweighting/Muon_DiBoson","DiMuon_","ZZ","DiMuon"),
     "mcDiMuonDY":("../June27_5fb_Reweighting/Muon_DY","DiMuon_","DY","DiMuon"),

     "mcPhoton":("../June27_5fb_Reweighting/Photon_MC","Photon_","Photon","Photon"),
     "ncPhoton":("../June27_5fb_Reweighting/Photon_Data","btag_morethantwo_Photon_","Data","Photon"),


    }


btag_more_than_zero_samples = {

    "nHad":("../June27_5fb_Reweighting/Had_Data","btag_morethanzero_","Data","Had"),
    
     "mcHadW1":("../June27_5fb_Reweighting/Had_WJets","","WJetsInc","Had"),
     "mcHadttbar":("../June27_5fb_Reweighting/Had_TTbar","","TTbar","Had"),
     "mcHadzinv":("../June27_5fb_Reweighting/Had_Zinv","","Zinv50","Had"),
     "mcHadsingt":("../June27_5fb_Reweighting/Had_SingleTop","","Single_Tbar_t","Had"),
     "mcHaddiboson":("../June27_5fb_Reweighting/Had_DiBoson","","ZZ","Had"),
     "mcHadDY":("../June27_5fb_Reweighting/Had_DY","","DY","Had"),


    "nMuon":("../June27_5fb_Reweighting/Muon_Data","btag_morethanzero_OneMuon_","Data","Muon"),
    
     "mcMuonW1":("../June27_5fb_Reweighting/Muon_WJets","OneMuon_","WJetsInc","Muon"),
     "mcMuonttbar":("../June27_5fb_Reweighting/Muon_TTbar","OneMuon_","TTbar","Muon"),
     "mcMuonzinv":("../June27_5fb_Reweighting/Muon_Zinv","OneMuon_","Zinv50","Muon"),
     "mcMuonsingt":("../June27_5fb_Reweighting/Muon_SingleTop","OneMuon_","Single_Tbar_t","Muon"),
     "mcMuondiboson":("../June27_5fb_Reweighting/Muon_DiBoson","OneMuon_","ZZ","Muon"),
     "mcMuonDY":("../June27_5fb_Reweighting/Muon_DY","OneMuon_","DY","Muon"),


    "nDiMuon":("../June27_5fb_Reweighting/Muon_Data","btag_morethanzero_DiMuon_","Data","DiMuon"),
    
     "mcDiMuonW1":("../June27_5fb_Reweighting/Muon_WJets","DiMuon_","WJetsInc","DiMuon"),
     "mcDiMuonttbar":("../June27_5fb_Reweighting/Muon_TTbar","DiMuon_","TTbar","DiMuon"),
     "mcDiMuonzinv":("../June27_5fb_Reweighting/Muon_Zinv","DiMuon_","Zinv50","DiMuon"),
     "mcDiMuonsingt":("../June27_5fb_Reweighting/Muon_SingleTop","DiMuon_","Single_Tbar_t","DiMuon"),
     "mcDiMuondiboson":("../June27_5fb_Reweighting/Muon_DiBoson","DiMuon_","ZZ","DiMuon"),
     "mcDiMuonDY":("../June27_5fb_Reweighting/Muon_DY","DiMuon_","DY","DiMuon"),

     "mcPhoton":("../June27_5fb_Reweighting/Photon_MC","Photon_","Photon","Photon"),

    }


btag_more_than_one_samples = {

    "nHad":("../June27_5fb_Reweighting/Had_Data","btag_morethanone_","Data","Had"),
    
     "mcHadW1":("../June27_5fb_Reweighting/Had_WJets","","WJetsInc","Had"),
     "mcHadttbar":("../June27_5fb_Reweighting/Had_TTbar","","TTbar","Had"),
     "mcHadzinv":("../June27_5fb_Reweighting/Had_Zinv","","Zinv50","Had"),
     "mcHadsingt":("../June27_5fb_Reweighting/Had_SingleTop","","Single_Tbar_t","Had"),
     "mcHaddiboson":("../June27_5fb_Reweighting/Had_DiBoson","","ZZ","Had"),
     "mcHadDY":("../June27_5fb_Reweighting/Had_DY","","DY","Had"),


    "nMuon":("../June27_5fb_Reweighting/Muon_Data","btag_morethanone_OneMuon_","Data","Muon"),
    
     "mcMuonW1":("../June27_5fb_Reweighting/Muon_WJets","OneMuon_","WJetsInc","Muon"),
     "mcMuonttbar":("../June27_5fb_Reweighting/Muon_TTbar","OneMuon_","TTbar","Muon"),
     "mcMuonzinv":("../June27_5fb_Reweighting/Muon_Zinv","OneMuon_","Zinv50","Muon"),
     "mcMuonsingt":("../June27_5fb_Reweighting/Muon_SingleTop","OneMuon_","Single_Tbar_t","Muon"),
     "mcMuondiboson":("../June27_5fb_Reweighting/Muon_DiBoson","OneMuon_","ZZ","Muon"),
     "mcMuonDY":("../June27_5fb_Reweighting/Muon_DY","OneMuon_","DY","Muon"),


    "nDiMuon":("../June27_5fb_Reweighting/Muon_Data","btag_morethanone_DiMuon_","Data","DiMuon"),
    
     "mcDiMuonW1":("../June27_5fb_Reweighting/Muon_WJets","DiMuon_","WJetsInc","DiMuon"),
     "mcDiMuonttbar":("../June27_5fb_Reweighting/Muon_TTbar","DiMuon_","TTbar","DiMuon"),
     "mcDiMuonzinv":("../June27_5fb_Reweighting/Muon_Zinv","DiMuon_","Zinv50","DiMuon"),
     "mcDiMuonsingt":("../June27_5fb_Reweighting/Muon_SingleTop","DiMuon_","Single_Tbar_t","DiMuon"),
     "mcDiMuondiboson":("../June27_5fb_Reweighting/Muon_DiBoson","DiMuon_","ZZ","DiMuon"),
     "mcDiMuonDY":("../June27_5fb_Reweighting/Muon_DY","DiMuon_","DY","DiMuon"),

     "mcPhoton":("../June27_5fb_Reweighting/Photon_MC","Photon_","Photon","Photon"),

    }




inclusive_samples = {

    "nHad":("../June27_5fb_Reweighting/Had_Data","","Data","Had"),
    
     "mcHadW1":("../June27_5fb_Reweighting/Had_WJets","","WJetsInc","Had"),
     "mcHadttbar":("../June27_5fb_Reweighting/Had_TTbar","","TTbar","Had"),
     "mcHadzinv":("../June27_5fb_Reweighting/Had_Zinv","","Zinv50","Had"),
     "mcHadsingt":("../June27_5fb_Reweighting/Had_SingleTop","","Single_Tbar_t","Had"),
     "mcHaddiboson":("../June27_5fb_Reweighting/Had_DiBoson","","ZZ","Had"),
     "mcHadDY":("../June27_5fb_Reweighting/Had_DY","","DY","Had"),


    "nMuon":("../June27_5fb_Reweighting/Muon_Data","OneMuon_","Data","Muon"),
    
     "mcMuonW1":("../June27_5fb_Reweighting/Muon_WJets","OneMuon_","WJetsInc","Muon"),
     "mcMuonttbar":("../June27_5fb_Reweighting/Muon_TTbar","OneMuon_","TTbar","Muon"),
     "mcMuonzinv":("../June27_5fb_Reweighting/Muon_Zinv","OneMuon_","Zinv50","Muon"),
     "mcMuonsingt":("../June27_5fb_Reweighting/Muon_SingleTop","OneMuon_","Single_Tbar_t","Muon"),
     "mcMuondiboson":("../June27_5fb_Reweighting/Muon_DiBoson","OneMuon_","ZZ","Muon"),
     "mcMuonDY":("../June27_5fb_Reweighting/Muon_DY","OneMuon_","DY","Muon"),


    "nDiMuon":("../June27_5fb_Reweighting/Muon_Data","DiMuon_","Data","DiMuon"),
    
     "mcDiMuonW1":("../June27_5fb_Reweighting/Muon_WJets","DiMuon_","WJetsInc","DiMuon"),
     "mcDiMuonttbar":("../June27_5fb_Reweighting/Muon_TTbar","DiMuon_","TTbar","DiMuon"),
     "mcDiMuonzinv":("../June27_5fb_Reweighting/Muon_Zinv","DiMuon_","Zinv50","DiMuon"),
     "mcDiMuonsingt":("../June27_5fb_Reweighting/Muon_SingleTop","DiMuon_","Single_Tbar_t","DiMuon"),
     "mcDiMuondiboson":("../June27_5fb_Reweighting/Muon_DiBoson","DiMuon_","ZZ","DiMuon"),
     "mcDiMuonDY":("../June27_5fb_Reweighting/Muon_DY","DiMuon_","DY","DiMuon"),

     "mcPhoton":("../June27_5fb_Reweighting/Photon_MC","Photon_","Photon","Photon"),

    }


calc_file = {
     "mchad":("../June26_5fb_NoReweighting/Had_MC.root","Had",""),
     "mchadzinv":("../June26_5fb_NoReweighting/Had_Zinv.root","Had_Zinv",""),
     "mcmuon":("../June26_5fb_NoReweighting/Muon_MC.root","Muon","OneMuon_"),
     "mcdimuon":("../June26_5fb_NoReweighting/Muon_MC.root","DiMuon","DiMuon_"),
     "mcphoton":("../June26_5fb_NoReweighting/Had_Zinv.root","Photon",""),

}

if __name__=="__main__":
  a = Number_Extractor(settings,btag_two_samples,"Two_btags",Triggers = "False",AlphaT="False",Calculation=calc_file,Stats = "True",Split_Lumi = "True")
  b = Number_Extractor(settings,btag_one_samples,"One_btag",Triggers = "False",AlphaT="False",Calculation=calc_file,Stats = "True",Split_Lumi = "True")
  c = Number_Extractor(settings,btag_zero_samples,"Zero_btags",Triggers = "False",AlphaT="False",Calculation=calc_file,Stats = "True",Split_Lumi = "True")
  d = Number_Extractor(settings,btag_more_than_two_samples,"More_Than_Two_btag",Triggers = "False",AlphaT="False",Calculation=calc_file,Stats = "True",Split_Lumi = "True")
  

  # We dont use these in the fit
  #e  = Number_Extractor(settings,btag_more_than_zero_samples,"More_Than_Zero_btag",Triggers = "False",AlphaT="False",Calculation=calc_file,Stats = "True",Split_Lumi = "True")
  #f  = Number_Extractor(settings,btag_more_than_one_samples,"More_Than_One_btag",Triggers = "False",AlphaT="False",Calculation=calc_file,Stats = "True",Split_Lumi = "True")
  #g  = Number_Extractor(settings,inclusive_samples,"Inclusive",Triggers = "False",AlphaT="False",Calculation=calc_file,Stats = "True",Split_Lumi = "True")
