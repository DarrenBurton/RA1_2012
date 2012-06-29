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
from plottingUtils import *
from Bryn_Numbers import *
from Jad_Compute import *



'''
Jads Closure Tests are used to show how we can relax the alphaT cut in the control samples
Therefore the AlphaTSlices are given as 0.55-10 and 0.01-10
'''

settings = {
  #"dirs":["275_325"],
  "dirs":["275_325","325_375","375_475","475_575","575_675","675_775","775_875","875",],
  "plots":["AlphaT_all",],
  "AlphaTSlices":["0.55_10","0.01_10"],
  "Lumo":4.98, 
  "Multi_Lumi":{'Had':4.98,'Muon':4.963,'DiMuon':4.963,'Photon':4.988},
  "Analysis":"8TeV"
      }

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

     #"mcPhoton":("../June27_5fb_Reweighting/Photon_MC_7TeV","Photon_","Photon","Photon"),
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

     #"mcPhoton":("../June27_5fb_Reweighting/Photon_MC_7TeV","Photon_","Photon","Photon"),
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

     #"mcPhoton":("../June27_5fb_Reweighting/Photon_MC_7TeV","Photon_","Photon","Photon"),
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
     
     #"mcPhoton":("../June27_5fb_Reweighting/Photon_MC_7TeV","Photon_","Photon","Photon"),
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

     #"mcPhoton":("../June27_5fb_Reweighting/Photon_MC_7TeV","Photon_","Photon","Photon"),
     "mcPhoton":("../June27_5fb_Reweighting/Photon_MC","Photon_","Photon","Photon"),
     "ncPhoton":("../June27_5fb_Reweighting/Photon_Data","btag_morethanzero_Photon_","Data","Photon"),


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

     #"mcPhoton":("../June27_5fb_Reweighting/Photon_MC_7TeV","Photon_","Photon","Photon"),
     "mcPhoton":("../June27_5fb_Reweighting/Photon_MC","Photon_","Photon","Photon"),
     "ncPhoton":("../June27_5fb_Reweighting/Photon_Data","btag_morethanone_Photon_","Data","Photon"),


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

     #"mcPhoton":("../June27_5fb_Reweighting/Photon_MC_7TeV","Photon_","Photon","Photon"),
     "mcPhoton":("../June27_5fb_Reweighting/Photon_MC","Photon_","Photon","Photon"),
     "ncPhoton":("../June27_5fb_Reweighting/Photon_Data","Photon_","Data","Photon"),


    }

calc_file = {
     "mchad":("../June26_5fb_NoReweighting/Had_MC.root","Had",""),
     "mchadzinv":("../June26_5fb_NoReweighting/Had_Zinv.root","Had_Zinv",""),
     "mcmuon":("../June26_5fb_NoReweighting/Muon_MC.root","Muon","OneMuon_"),
     "mcdimuon":("../June26_5fb_NoReweighting/Muon_MC.root","DiMuon","DiMuon_"),
     "mcphoton":("../June26_5fb_NoReweighting/Had_Zinv.root","Photon",""),

}

'''
All Sample are passed through Number_Extractor and each of the relevant dictionaries are output to c_file.
Jad_Compute which is then called and closure test png's are produced.
3rd argument here is entered into dictionary to identify btag multiplicity. Dont change these strings

Addtional options
Classic - If true then produces 'classic baseline' closure tests. i.e. Photon -> dimuon, muon -> dimuon. Use with basline files only, no btags

Look in Jad_Compute.py for any further comments
'''


if __name__=="__main__":
  LIST_FOR_JAD = []
  a = Number_Extractor(settings,btag_two_samples,"Two_btags",c_file = LIST_FOR_JAD,Closure = "True",Triggers = "True",AlphaT="True",Calculation=calc_file,Split_Lumi = "True")
  b = Number_Extractor(settings,btag_one_samples,"One_btag",c_file = LIST_FOR_JAD,Closure = "True",Triggers = "True",AlphaT="True",Calculation=calc_file,Split_Lumi = "True")
  c = Number_Extractor(settings,btag_zero_samples,"Zero_btags",c_file = LIST_FOR_JAD,Closure = "True",Triggers = "True",AlphaT="True",Calculation=calc_file,Split_Lumi = "True")
  d = Number_Extractor(settings,btag_more_than_two_samples,"More_Than_Two_btag",c_file = LIST_FOR_JAD,Closure = "True",Triggers = "True",AlphaT="True",Calculation=calc_file,Split_Lumi = "True")
  e  = Number_Extractor(settings,btag_more_than_zero_samples,"More_Than_Zero_btag",c_file = LIST_FOR_JAD,Closure = "True",Triggers = "True",AlphaT="True",Calculation=calc_file,Split_Lumi = "True")
  f  = Number_Extractor(settings,btag_more_than_one_samples,"More_Than_One_btag",c_file = LIST_FOR_JAD,Closure = "True",Triggers = "True",AlphaT="True",Calculation=calc_file,Split_Lumi = "True")
  g  = Number_Extractor(settings,inclusive_samples,"Inclusive",c_file = LIST_FOR_JAD,Closure = "True",Triggers = "True",AlphaT="True",Calculation=calc_file,Split_Lumi = "True")
  h = Jad_Compute(LIST_FOR_JAD,classic ="False",Lumo = settings["Lumo"])
