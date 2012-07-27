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
  "dirs":["275_325","325_375","375_475","475_575","575_675","675_775","775_875","875",],  #HT Bins
  "plots":["AlphaT_all",],  # Histogram that Yields are taken from
  "AlphaTSlices":["0.55_20"], # AlphaT Slices
  "Lumo":4.98, # Luminosity in fb
  "Multi_Lumi":{'Had':4.98,'Muon':4.963,'DiMuon':4.963,'Photon':4.988},  # Different Luminosity per sample, used when SplitLumi = True
  "Analysis":"8TeV" # Differentiate between 7 and 8 TeV analysis i.e. uses alphaT cut in lowest two bins if 7TeV is selected
      }

'''
Sample Dictionary Instructions
eg "nMuon":("./May_11_root_Files/Muon_Data","btag_two_OneMuon_","Data","Muon"),
if n at start of name entry then the file is data and will no be scaled to luminosity.
first argument is path to root file
second argument is prefix to ht bin. i.e OneMuon_275_325
third argument is data/mc type, i.e. Data, WJets250 - MC relating to the binned WJets 250-300 HT sample
fourth argument is sample Type, Had/DiMuon/Muon. 

the only thing that will have to be changed is the second argument depending on wether you are running btag multiplicity/baseline
'''
btag_two_samples = {

    "nHad":("../BtagAlgo7_5fb/Had_Data","btag_two_","Data","Had"),
    
     "mcHadW1":("../BtagAlgo7_5fb/Had_WJets","","WJetsInc","Had"),
     "mcHadttbar":("../BtagAlgo7_5fb/Had_TTbar","","TTbar","Had"),
     "mcHadzinv":("../BtagAlgo7_5fb/Had_Zinv","","Zinv50","Had"),
     "mcHadsingt":("../BtagAlgo7_5fb/Had_SingleTop","","Single_Tbar_t","Had"),
     "mcHaddiboson":("../BtagAlgo7_5fb/Had_DiBoson","","ZZ","Had"),
     "mcHadDY":("../BtagAlgo7_5fb/Had_DY","","DY","Had"),

    "nMuon":("../BtagAlgo7_5fb/Muon_Data","btag_two_OneMuon_","Data","Muon"),
    
     "mcMuonW1":("../BtagAlgo7_5fb/Muon_WJets","OneMuon_","WJetsInc","Muon"),
     "mcMuonttbar":("../BtagAlgo7_5fb/Muon_TTbar","OneMuon_","TTbar","Muon"),
     "mcMuonzinv":("../BtagAlgo7_5fb/Muon_Zinv","OneMuon_","Zinv50","Muon"),
     "mcMuonsingt":("../BtagAlgo7_5fb/Muon_SingleTop","OneMuon_","Single_Tbar_t","Muon"),
     "mcMuondiboson":("../BtagAlgo7_5fb/Muon_DiBoson","OneMuon_","ZZ","Muon"),
     "mcMuonDY":("../BtagAlgo7_5fb/Muon_DY","OneMuon_","DY","Muon"),



    "nDiMuon":("../BtagAlgo7_5fb/Muon_Data","btag_two_DiMuon_","Data","DiMuon"),
    
     "mcDiMuonW1":("../BtagAlgo7_5fb/Muon_WJets","DiMuon_","WJetsInc","DiMuon"),
     "mcDiMuonttbar":("../BtagAlgo7_5fb/Muon_TTbar","DiMuon_","TTbar","DiMuon"),
     "mcDiMuonzinv":("../BtagAlgo7_5fb/Muon_Zinv","DiMuon_","Zinv50","DiMuon"),
     "mcDiMuonsingt":("../BtagAlgo7_5fb/Muon_SingleTop","DiMuon_","Single_Tbar_t","DiMuon"),
     "mcDiMuondiboson":("../BtagAlgo7_5fb/Muon_DiBoson","DiMuon_","ZZ","DiMuon"),
     "mcDiMuonDY":("../BtagAlgo7_5fb/Muon_DY","DiMuon_","DY","DiMuon"),

     "mcPhoton":("../BtagAlgo7_5fb/Photon_MC","Photon_","Photon","Photon"),
     "nPhoton":("../BtagAlgo7_5fb/Photon_Data","btag_two_Photon_","Data","Photon"),



    }

btag_two_uncorrected_samples = {

    "nHad":("../BtagAlgo7_5fb/Had_Data","btag_two_","Data","Had"),
    
     "mcHadW1":("../BtagAlgo7_5fb/Had_WJets","btag_two_","WJetsInc","Had"),
     "mcHadttbar":("../BtagAlgo7_5fb/Had_TTbar","btag_two_","TTbar","Had"),
     "mcHadzinv":("../BtagAlgo7_5fb/Had_Zinv","btag_two_","Zinv50","Had"),
     "mcHadsingt":("../BtagAlgo7_5fb/Had_SingleTop","btag_two_","Single_Tbar_t","Had"),
     "mcHaddiboson":("../BtagAlgo7_5fb/Had_DiBoson","btag_two_","ZZ","Had"),
     "mcHadDY":("../BtagAlgo7_5fb/Had_DY","btag_two_","DY","Had"),

    "nMuon":("../BtagAlgo7_5fb/Muon_Data","btag_two_OneMuon_","Data","Muon"),
    
     "mcMuonW1":("../BtagAlgo7_5fb/Muon_WJets","btag_two_OneMuon_","WJetsInc","Muon"),
     "mcMuonttbar":("../BtagAlgo7_5fb/Muon_TTbar","btag_two_OneMuon_","TTbar","Muon"),
     "mcMuonzinv":("../BtagAlgo7_5fb/Muon_Zinv","btag_two_OneMuon_","Zinv50","Muon"),
     "mcMuonsingt":("../BtagAlgo7_5fb/Muon_SingleTop","btag_two_OneMuon_","Single_Tbar_t","Muon"),
     "mcMuondiboson":("../BtagAlgo7_5fb/Muon_DiBoson","btag_two_OneMuon_","ZZ","Muon"),
     "mcMuonDY":("../BtagAlgo7_5fb/Muon_DY","btag_two_OneMuon_","DY","Muon"),



    "nDiMuon":("../BtagAlgo7_5fb/Muon_Data","btag_two_DiMuon_","Data","DiMuon"),
    
     "mcDiMuonW1":("../BtagAlgo7_5fb/Muon_WJets","btag_two_DiMuon_","WJetsInc","DiMuon"),
     "mcDiMuonttbar":("../BtagAlgo7_5fb/Muon_TTbar","btag_two_DiMuon_","TTbar","DiMuon"),
     "mcDiMuonzinv":("../BtagAlgo7_5fb/Muon_Zinv","btag_two_DiMuon_","Zinv50","DiMuon"),
     "mcDiMuonsingt":("../BtagAlgo7_5fb/Muon_SingleTop","btag_two_DiMuon_","Single_Tbar_t","DiMuon"),
     "mcDiMuondiboson":("../BtagAlgo7_5fb/Muon_DiBoson","btag_two_DiMuon_","ZZ","DiMuon"),
     "mcDiMuonDY":("../BtagAlgo7_5fb/Muon_DY","btag_two_DiMuon_","DY","DiMuon"),

     #"mcPhoton":("../BtagAlgo7_5fb/Photon_MC_7TeV","Photon_","Photon","Photon"),



    }


btag_one_samples = {

    "nHad":("../BtagAlgo7_5fb/Had_Data","btag_one_","Data","Had"),
    
     "mcHadW1":("../BtagAlgo7_5fb/Had_WJets","","WJetsInc","Had"),
     "mcHadttbar":("../BtagAlgo7_5fb/Had_TTbar","","TTbar","Had"),
     "mcHadzinv":("../BtagAlgo7_5fb/Had_Zinv","","Zinv50","Had"),
     "mcHadsingt":("../BtagAlgo7_5fb/Had_SingleTop","","Single_Tbar_t","Had"),
     "mcHaddiboson":("../BtagAlgo7_5fb/Had_DiBoson","","ZZ","Had"),
     "mcHadDY":("../BtagAlgo7_5fb/Had_DY","","DY","Had"),


    "nMuon":("../BtagAlgo7_5fb/Muon_Data","btag_one_OneMuon_","Data","Muon"),
    
    "mcMuonW1":("../BtagAlgo7_5fb/Muon_WJets","OneMuon_","WJetsInc","Muon"),
     "mcMuonttbar":("../BtagAlgo7_5fb/Muon_TTbar","OneMuon_","TTbar","Muon"),
     "mcMuonzinv":("../BtagAlgo7_5fb/Muon_Zinv","OneMuon_","Zinv50","Muon"),
     "mcMuonsingt":("../BtagAlgo7_5fb/Muon_SingleTop","OneMuon_","Single_Tbar_t","Muon"),
     "mcMuondiboson":("../BtagAlgo7_5fb/Muon_DiBoson","OneMuon_","ZZ","Muon"),
     "mcMuonDY":("../BtagAlgo7_5fb/Muon_DY","OneMuon_","DY","Muon"),


    "nDiMuon":("../BtagAlgo7_5fb/Muon_Data","btag_one_DiMuon_","Data","DiMuon"),
    
     "mcDiMuonW1":("../BtagAlgo7_5fb/Muon_WJets","DiMuon_","WJetsInc","DiMuon"),
     "mcDiMuonttbar":("../BtagAlgo7_5fb/Muon_TTbar","DiMuon_","TTbar","DiMuon"),
     "mcDiMuonzinv":("../BtagAlgo7_5fb/Muon_Zinv","DiMuon_","Zinv50","DiMuon"),
     "mcDiMuonsingt":("../BtagAlgo7_5fb/Muon_SingleTop","DiMuon_","Single_Tbar_t","DiMuon"),
     "mcDiMuondiboson":("../BtagAlgo7_5fb/Muon_DiBoson","DiMuon_","ZZ","DiMuon"),
     "mcDiMuonDY":("../BtagAlgo7_5fb/Muon_DY","DiMuon_","DY","DiMuon"),
     
     "mcPhoton":("../BtagAlgo7_5fb/Photon_MC","Photon_","Photon","Photon"),
     "nPhoton":("../BtagAlgo7_5fb/Photon_Data","btag_one_Photon_","Data","Photon"),



    }

btag_zero_samples = {
     
     "nHad":("../BtagAlgo7_5fb/Had_Data","btag_zero_","Data","Had"),
    
     "mcHadW1":("../BtagAlgo7_5fb/Had_WJets","","WJetsInc","Had"),
     "mcHadttbar":("../BtagAlgo7_5fb/Had_TTbar","","TTbar","Had"),
     "mcHadzinv":("../BtagAlgo7_5fb/Had_Zinv","","Zinv50","Had"),
     "mcHadsingt":("../BtagAlgo7_5fb/Had_SingleTop","","Single_Tbar_t","Had"),
     "mcHaddiboson":("../BtagAlgo7_5fb/Had_DiBoson","","ZZ","Had"),
     "mcHadDY":("../BtagAlgo7_5fb/Had_DY","","DY","Had"),


    "nMuon":("../BtagAlgo7_5fb/Muon_Data","btag_zero_OneMuon_","Data","Muon"),
    
     "mcMuonW1":("../BtagAlgo7_5fb/Muon_WJets","OneMuon_","WJetsInc","Muon"),
     "mcMuonttbar":("../BtagAlgo7_5fb/Muon_TTbar","OneMuon_","TTbar","Muon"),
     "mcMuonzinv":("../BtagAlgo7_5fb/Muon_Zinv","OneMuon_","Zinv50","Muon"),
     "mcMuonsingt":("../BtagAlgo7_5fb/Muon_SingleTop","OneMuon_","Single_Tbar_t","Muon"),
     "mcMuondiboson":("../BtagAlgo7_5fb/Muon_DiBoson","OneMuon_","ZZ","Muon"),
     "mcMuonDY":("../BtagAlgo7_5fb/Muon_DY","OneMuon_","DY","Muon"),


    "nDiMuon":("../BtagAlgo7_5fb/Muon_Data","btag_zero_DiMuon_","Data","DiMuon"),
    
     "mcDiMuonW1":("../BtagAlgo7_5fb/Muon_WJets","DiMuon_","WJetsInc","DiMuon"),
     "mcDiMuonttbar":("../BtagAlgo7_5fb/Muon_TTbar","DiMuon_","TTbar","DiMuon"),
     "mcDiMuonzinv":("../BtagAlgo7_5fb/Muon_Zinv","DiMuon_","Zinv50","DiMuon"),
     "mcDiMuonsingt":("../BtagAlgo7_5fb/Muon_SingleTop","DiMuon_","Single_Tbar_t","DiMuon"),
     "mcDiMuondiboson":("../BtagAlgo7_5fb/Muon_DiBoson","DiMuon_","ZZ","DiMuon"),
     "mcDiMuonDY":("../BtagAlgo7_5fb/Muon_DY","DiMuon_","DY","DiMuon"),
     
     "mcPhoton":("../BtagAlgo7_5fb/Photon_MC","Photon_","Photon","Photon"),
     "ncPhoton":("../BtagAlgo7_5fb/Photon_Data","btag_zero_Photon_","Data","Photon"),



    }


btag_one_uncorrected_samples = {

     "nHad":("../BtagAlgo7_5fb/Had_Data","btag_one_","Data","Had"),
    
     "mcHadW1":("../BtagAlgo7_5fb/Had_WJets","btag_one_","WJetsInc","Had"),
     "mcHadttbar":("../BtagAlgo7_5fb/Had_TTbar","btag_one_","TTbar","Had"),
     "mcHadzinv":("../BtagAlgo7_5fb/Had_Zinv","btag_one_","Zinv50","Had"),
     "mcHadsingt":("../BtagAlgo7_5fb/Had_SingleTop","btag_one_","Single_Tbar_t","Had"),
     "mcHaddiboson":("../BtagAlgo7_5fb/Had_DiBoson","btag_one_","ZZ","Had"),
     "mcHadDY":("../BtagAlgo7_5fb/Had_DY","btag_one_","DY","Had"),


    "nMuon":("../BtagAlgo7_5fb/Muon_Data","btag_one_OneMuon_","Data","Muon"),
    
     "mcMuonW1":("../BtagAlgo7_5fb/Muon_WJets","btag_one_OneMuon_","WJetsInc","Muon"),
     "mcMuonttbar":("../BtagAlgo7_5fb/Muon_TTbar","btag_one_OneMuon_","TTbar","Muon"),
     "mcMuonzinv":("../BtagAlgo7_5fb/Muon_Zinv","btag_one_OneMuon_","Zinv50","Muon"),
     "mcMuonsingt":("../BtagAlgo7_5fb/Muon_SingleTop","btag_one_OneMuon_","Single_Tbar_t","Muon"),
     "mcMuondiboson":("../BtagAlgo7_5fb/Muon_DiBoson","btag_one_OneMuon_","ZZ","Muon"),
     "mcMuonDY":("../BtagAlgo7_5fb/Muon_DY","btag_one_OneMuon_","DY","Muon"),


    "nDiMuon":("../BtagAlgo7_5fb/Muon_Data","btag_one_DiMuon_","Data","DiMuon"),
    
     "mcDiMuonW1":("../BtagAlgo7_5fb/Muon_WJets","btag_one_DiMuon_","WJetsInc","DiMuon"),
     "mcDiMuonttbar":("../BtagAlgo7_5fb/Muon_TTbar","btag_one_DiMuon_","TTbar","DiMuon"),
     "mcDiMuonzinv":("../BtagAlgo7_5fb/Muon_Zinv","btag_one_DiMuon_","Zinv50","DiMuon"),
     "mcDiMuonsingt":("../BtagAlgo7_5fb/Muon_SingleTop","btag_one_DiMuon_","Single_Tbar_t","DiMuon"),
     "mcDiMuondiboson":("../BtagAlgo7_5fb/Muon_DiBoson","btag_one_DiMuon_","ZZ","DiMuon"),
     "mcDiMuonDY":("../BtagAlgo7_5fb/Muon_DY","btag_one_DiMuon_","DY","DiMuon"),

    }



btag_zero_uncorrected_samples = {

    "nHad":("../BtagAlgo7_5fb/Had_Data","btag_zero_","Data","Had"),
    
     "mcHadW1":("../BtagAlgo7_5fb/Had_WJets","btag_zero_","WJetsInc","Had"),
     "mcHadttbar":("../BtagAlgo7_5fb/Had_TTbar","btag_zero_","TTbar","Had"),
     "mcHadzinv":("../BtagAlgo7_5fb/Had_Zinv","btag_zero_","Zinv50","Had"),
     "mcHadsingt":("../BtagAlgo7_5fb/Had_SingleTop","btag_zero_","Single_Tbar_t","Had"),
     "mcHaddiboson":("../BtagAlgo7_5fb/Had_DiBoson","btag_zero_","ZZ","Had"),
     "mcHadDY":("../BtagAlgo7_5fb/Had_DY","btag_zero_","DY","Had"),


    "nMuon":("../BtagAlgo7_5fb/Muon_Data","btag_zero_OneMuon_","Data","Muon"),
    
     "mcMuonW1":("../BtagAlgo7_5fb/Muon_WJets","btag_zero_OneMuon_","WJetsInc","Muon"),
     "mcMuonttbar":("../BtagAlgo7_5fb/Muon_TTbar","btag_zero_OneMuon_","TTbar","Muon"),
     "mcMuonzinv":("../BtagAlgo7_5fb/Muon_Zinv","btag_zero_OneMuon_","Zinv50","Muon"),
     "mcMuonsingt":("../BtagAlgo7_5fb/Muon_SingleTop","btag_zero_OneMuon_","Single_Tbar_t","Muon"),
     "mcMuondiboson":("../BtagAlgo7_5fb/Muon_DiBoson","btag_zero_OneMuon_","ZZ","Muon"),
     "mcMuonDY":("../BtagAlgo7_5fb/Muon_DY","btag_zero_OneMuon_","DY","Muon"),

    "nDiMuon":("../BtagAlgo7_5fb/Muon_Data","btag_zero_DiMuon_","Data","DiMuon"),
    
     "mcDiMuonW1":("../BtagAlgo7_5fb/Muon_WJets","btag_zero_DiMuon_","WJetsInc","DiMuon"),
     "mcDiMuonttbar":("../BtagAlgo7_5fb/Muon_TTbar","btag_zero_DiMuon_","TTbar","DiMuon"),
     "mcDiMuonzinv":("../BtagAlgo7_5fb/Muon_Zinv","btag_zero_DiMuon_","Zinv50","DiMuon"),
     "mcDiMuonsingt":("../BtagAlgo7_5fb/Muon_SingleTop","btag_zero_DiMuon_","Single_Tbar_t","DiMuon"),
     "mcDiMuondiboson":("../BtagAlgo7_5fb/Muon_DiBoson","btag_zero_DiMuon_","ZZ","DiMuon"),
     "mcDiMuonDY":("../BtagAlgo7_5fb/Muon_DY","btag_zero_DiMuon_","DY","DiMuon"),

    }


btag_more_than_two_samples = {

    "nHad":("../BtagAlgo7_5fb/Had_Data","btag_morethantwo_","Data","Had"),
    
     "mcHadW1":("../BtagAlgo7_5fb/Had_WJets","","WJetsInc","Had"),
     "mcHadttbar":("../BtagAlgo7_5fb/Had_TTbar","","TTbar","Had"),
     "mcHadzinv":("../BtagAlgo7_5fb/Had_Zinv","","Zinv50","Had"),
     "mcHadsingt":("../BtagAlgo7_5fb/Had_SingleTop","","Single_Tbar_t","Had"),
     "mcHaddiboson":("../BtagAlgo7_5fb/Had_DiBoson","","ZZ","Had"),
     "mcHadDY":("../BtagAlgo7_5fb/Had_DY","","DY","Had"),


    "nMuon":("../BtagAlgo7_5fb/Muon_Data","btag_morethantwo_OneMuon_","Data","Muon"),
    
     "mcMuonW1":("../BtagAlgo7_5fb/Muon_WJets","OneMuon_","WJetsInc","Muon"),
     "mcMuonttbar":("../BtagAlgo7_5fb/Muon_TTbar","OneMuon_","TTbar","Muon"),
     "mcMuonzinv":("../BtagAlgo7_5fb/Muon_Zinv","OneMuon_","Zinv50","Muon"),
     "mcMuonsingt":("../BtagAlgo7_5fb/Muon_SingleTop","OneMuon_","Single_Tbar_t","Muon"),
     "mcMuondiboson":("../BtagAlgo7_5fb/Muon_DiBoson","OneMuon_","ZZ","Muon"),
     "mcMuonDY":("../BtagAlgo7_5fb/Muon_DY","OneMuon_","DY","Muon"),


    "nDiMuon":("../BtagAlgo7_5fb/Muon_Data","btag_morethantwo_DiMuon_","Data","DiMuon"),
    
     "mcDiMuonW1":("../BtagAlgo7_5fb/Muon_WJets","DiMuon_","WJetsInc","DiMuon"),
     "mcDiMuonttbar":("../BtagAlgo7_5fb/Muon_TTbar","DiMuon_","TTbar","DiMuon"),
     "mcDiMuonzinv":("../BtagAlgo7_5fb/Muon_Zinv","DiMuon_","Zinv50","DiMuon"),
     "mcDiMuonsingt":("../BtagAlgo7_5fb/Muon_SingleTop","DiMuon_","Single_Tbar_t","DiMuon"),
     "mcDiMuondiboson":("../BtagAlgo7_5fb/Muon_DiBoson","DiMuon_","ZZ","DiMuon"),
     "mcDiMuonDY":("../BtagAlgo7_5fb/Muon_DY","DiMuon_","DY","DiMuon"),
     
     "mcPhoton":("../BtagAlgo7_5fb/Photon_MC","Photon_","Photon","Photon"),
     "ncPhoton":("../BtagAlgo7_5fb/Photon_Data","btag_morethantwo_Photon_","Data","Photon"),



    }

btag_more_than_two_uncorrected_samples = {

    "nHad":("../BtagAlgo7_5fb/Had_Data","btag_morethantwo_","Data","Had"),
    
     "mcHadW1":("../BtagAlgo7_5fb/Had_WJets","btag_morethantwo_","WJetsInc","Had"),
     "mcHadttbar":("../BtagAlgo7_5fb/Had_TTbar","btag_morethantwo_","TTbar","Had"),
     "mcHadzinv":("../BtagAlgo7_5fb/Had_Zinv","btag_morethantwo_","Zinv50","Had"),
     "mcHadsingt":("../BtagAlgo7_5fb/Had_SingleTop","btag_morethantwo_","Single_Tbar_t","Had"),
     "mcHaddiboson":("../BtagAlgo7_5fb/Had_DiBoson","btag_morethantwo_","ZZ","Had"),
     "mcHadDY":("../BtagAlgo7_5fb/Had_DY","btag_morethantwo_","DY","Had"),


    "nMuon":("../BtagAlgo7_5fb/Muon_Data","btag_morethantwo_OneMuon_","Data","Muon"),
    
     "mcMuonW1":("../BtagAlgo7_5fb/Muon_WJets","btag_morethantwo_OneMuon_","WJetsInc","Muon"),
     "mcMuonttbar":("../BtagAlgo7_5fb/Muon_TTbar","btag_morethantwo_OneMuon_","TTbar","Muon"),
     "mcMuonzinv":("../BtagAlgo7_5fb/Muon_Zinv","btag_morethantwo_OneMuon_","Zinv50","Muon"),
     "mcMuonsingt":("../BtagAlgo7_5fb/Muon_SingleTop","btag_morethantwo_OneMuon_","Single_Tbar_t","Muon"),
     "mcMuondiboson":("../BtagAlgo7_5fb/Muon_DiBoson","btag_morethantwo_OneMuon_","ZZ","Muon"),
     "mcMuonDY":("../BtagAlgo7_5fb/Muon_DY","btag_morethantwo_OneMuon_","DY","Muon"),

    "nDiMuon":("../BtagAlgo7_5fb/Muon_Data","btag_morethantwo_DiMuon_","Data","DiMuon"),
    
     "mcDiMuonW1":("../BtagAlgo7_5fb/Muon_WJets","btag_morethantwo_DiMuon_","WJetsInc","DiMuon"),
     "mcDiMuonttbar":("../BtagAlgo7_5fb/Muon_TTbar","btag_morethantwo_DiMuon_","TTbar","DiMuon"),
     "mcDiMuonzinv":("../BtagAlgo7_5fb/Muon_Zinv","btag_morethantwo_DiMuon_","Zinv50","DiMuon"),
     "mcDiMuonsingt":("../BtagAlgo7_5fb/Muon_SingleTop","btag_morethantwo_DiMuon_","Single_Tbar_t","DiMuon"),
     "mcDiMuondiboson":("../BtagAlgo7_5fb/Muon_DiBoson","btag_morethantwo_DiMuon_","ZZ","DiMuon"),
     "mcDiMuonDY":("../BtagAlgo7_5fb/Muon_DY","btag_morethantwo_DiMuon_","DY","DiMuon"),

    }

inclusive_samples = {

    "nHad":("../BtagAlgo7_5fb/Had_Data","","Data","Had"),
    
     "mcHadW1":("../BtagAlgo7_5fb/Had_WJets","","WJetsInc","Had"),
     "mcHadttbar":("../BtagAlgo7_5fb/Had_TTbar","","TTbar","Had"),
     "mcHadzinv":("../BtagAlgo7_5fb/Had_Zinv","","Zinv50","Had"),
     "mcHadsingt":("../BtagAlgo7_5fb/Had_SingleTop","","Single_Tbar_t","Had"),
     "mcHaddiboson":("../BtagAlgo7_5fb/Had_DiBoson","","ZZ","Had"),
     "mcHadDY":("../BtagAlgo7_5fb/Had_DY","","DY","Had"),


    "nMuon":("../BtagAlgo7_5fb/Muon_Data","OneMuon_","Data","Muon"),
    
     "mcMuonW1":("../BtagAlgo7_5fb/Muon_WJets","OneMuon_","WJetsInc","Muon"),
     "mcMuonttbar":("../BtagAlgo7_5fb/Muon_TTbar","OneMuon_","TTbar","Muon"),
     "mcMuonzinv":("../BtagAlgo7_5fb/Muon_Zinv","OneMuon_","Zinv50","Muon"),
     "mcMuonsingt":("../BtagAlgo7_5fb/Muon_SingleTop","OneMuon_","Single_Tbar_t","Muon"),
     "mcMuondiboson":("../BtagAlgo7_5fb/Muon_DiBoson","OneMuon_","ZZ","Muon"),
     "mcMuonDY":("../BtagAlgo7_5fb/Muon_DY","OneMuon_","DY","Muon"),


    "nDiMuon":("../BtagAlgo7_5fb/Muon_Data","DiMuon_","Data","DiMuon"),
    
     "mcDiMuonW1":("../BtagAlgo7_5fb/Muon_WJets","DiMuon_","WJetsInc","DiMuon"),
     "mcDiMuonttbar":("../BtagAlgo7_5fb/Muon_TTbar","DiMuon_","TTbar","DiMuon"),
     "mcDiMuonzinv":("../BtagAlgo7_5fb/Muon_Zinv","DiMuon_","Zinv50","DiMuon"),
     "mcDiMuonsingt":("../BtagAlgo7_5fb/Muon_SingleTop","DiMuon_","Single_Tbar_t","DiMuon"),
     "mcDiMuondiboson":("../BtagAlgo7_5fb/Muon_DiBoson","DiMuon_","ZZ","DiMuon"),
     "mcDiMuonDY":("../BtagAlgo7_5fb/Muon_DY","DiMuon_","DY","DiMuon"),

     "mcPhoton":("../BtagAlgo7_5fb/Photon_MC","Photon_","Photon","Photon"),
     "ncPhoton":("../BtagAlgo7_5fb/Photon_Data","Photon_","Data","Photon"),


    }


calc_file = {
     "mchad":("../BtagAlgo7_5fb/Had_MC.root","Had",""),
     "mchadzinv":("../BtagAlgo7_5fb/Had_Zinv.root","Had_Zinv",""),
     "mcmuon":("../BtagAlgo7_5fb/Muon_MC.root","Muon","OneMuon_"),
     "mcdimuon":("../BtagAlgo7_5fb/Muon_MC.root","DiMuon","DiMuon_"),
     "mcphoton":("../BtagAlgo7_5fb/Had_Zinv.root","Photon",""),

}



if __name__=="__main__":

  # Formula Method
  a = Number_Extractor(settings,btag_two_samples,"Two_btags",Triggers = "True",AlphaT="False",Calculation=calc_file,Split_Lumi = "True")
  #b = Number_Extractor(settings,btag_one_samples,"One_btag",Triggers = "True",AlphaT="False",Calculation=calc_file,Split_Lumi = "True")
  #c = Number_Extractor(settings,btag_zero_samples,"Zero_btags",Triggers = "True",AlphaT="False",Calculation=calc_file,Split_Lumi = "True")
  #d = Number_Extractor(settings,btag_more_than_two_samples,"More_Than_Two_btag",Triggers = "True",AlphaT="False",Calculation=calc_file,Split_Lumi = "True")

  # Vanilla Yields
  #a = Number_Extractor(settings,btag_two_uncorrected_samples,"Two_btags",Triggers = "True",AlphaT="False",Split_Lumi = "True")
  #b = Number_Extractor(settings,btag_one_uncorrected_samples,"One_btag",Triggers = "True",AlphaT="False",Split_Lumi = "True")
  #c = Number_Extractor(settings,btag_zero_uncorrected_samples,"Zero_btags",Triggers = "True",AlphaT="False",Split_Lumi = "True")
  #d = Number_Extractor(settings,btag_more_than_two_uncorrected_samples,"More_Than_Two_btag",Triggers = "True",AlphaT="False",Split_Lumi = "True")
  
  #Inclusive sample
  #g  = Number_Extractor(settings,inclusive_samples,"Inclusive",Triggers = "True",AlphaT="False",Split_Lumi = "True",Calculation=calc_file)
  #g  = Number_Extractor(settings,inclusive_samples,"Inclusive",Triggers = "True",AlphaT="False",Split_Lumi = "True")


