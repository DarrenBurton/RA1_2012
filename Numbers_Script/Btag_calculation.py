#!/usr/bin/env python
from ROOT import *
import ROOT as r
import logging,itertools
import os,fnmatch,sys
import glob, errno
from time import strftime
from optparse import OptionParser
import array, ast
from math import *
from plottingUtils import *
from Bryn_Numbers import *

class Btag_Calc(object):
        
   def __init__(self,settings,samples,btag_measure,number,alphaT_check,lumi_dict,analysis):
        print "In btag Calc taking samples"

        print "\n\n Making Predictions for %s \n\n" %number
        self.Keep_AlphaT = alphaT_check
        self.settings = settings
        self.Analysis = analysis
        self.btag_multiplicity = btag_measure
        self.Btag_Rate(btag_measure)
        self.Rate_Correct(self.Btag_Efficiencies) # hack to fix wierd fluctuation in mist rate in muon HT 475 sample

        #Use to produce Btag - Mistage plots
        #self.DiMuon_Fit(self.Btag_Efficiencies,"Had","Mistag")
        #self.DiMuon_Fit(self.Btag_Efficiencies,"Muon","Mistag")
        #self.DiMuon_Fit(self.Btag_Efficiencies,"Had","Btag")
        #self.DiMuon_Fit(self.Btag_Efficiencies,"Muon","Btag")
        self.lumi_dict = lumi_dict

   def Rate_Correct(self,dictionary):
  
      dictionary['Muon']['475']["Mistag_Efficiency"] =  dictionary['Muon']['375']["Mistag_Efficiency"]
      dictionary['Muon']['475']["Mistag_Error"] = dictionary['Muon']['375']["Mistag_Error"]

      return dictionary

   def DiMuon_Fit(self,dictionary,sample,tag):
        
        c1= r.TCanvas("Yields", "Yields",0,0,900,600)
        c1.cd()
        fit = r.TF1("fit","pol0",275,975)
        data = r.TGraphAsymmErrors(8)
        i = 0
        for num,entry in sorted(dictionary[sample].iteritems()):
                i+=1
                data.SetPoint(i,float(num),dictionary[sample][num]["%s_Efficiency"%tag ])
                data.SetPointEYhigh(i,dictionary[sample][num]["%s_Error"%tag])
                data.SetPointEYlow(i,dictionary[sample][num]["%s_Error"%tag])

        data.GetXaxis().SetRangeUser(275,975.)
        if tag == "Mistag" :data.GetYaxis().SetRangeUser(0,0.1)
        if tag == "Btag" : data.GetYaxis().SetRangeUser(0.5,1.0)
        data.SetTitle("")
        data.GetXaxis().SetTitle("H_{T} (GeV)")
        data.GetYaxis().SetTitle("%s Efficiency"%tag)
        data.GetYaxis().SetTitleOffset(1.4)
        data.GetXaxis().SetTitleOffset(0.8)
        data.SetLineWidth(3)
        data.SetMarkerStyle(20)
        data.SetMarkerSize(1.5)
        gStyle.SetOptStat(0)
        data.Draw("AP")
        #data.Fit(fit)
        c1.SaveAs("Uncorrected_%s_%s_Fit.png"%(sample,tag) )
        #print fit.GetParameter(0)
        #for num,entry in sorted(dictionary[sample].iteritems()):
        #        self.Btag_Efficiencies[sample][num]["Mistag_Efficiency"] = fit.GetParameter(0)

   def Make_Dict(self,settings,samples,number):
        htbins = ["275_325","325_375","375_475","475_575","575_675","675_775","775_875","875"]
        table_entries = "{"
        for key,fi in sorted(samples.iteritems()):
              file = r.TFile("%s.root"%fi[0])
              DirKeys = file.GetListOfKeys()  
              i = 0
              for dir in settings['dirs']:
                fixed_dir = dir
                for alphat in settings['AlphaTSlices']:
                   dir = fixed_dir
                   lower = alphat.split('_')[0]
                   higher = alphat.split('_')[1]
                   table_entries += "\t\"%s_%d\"  : "%(key,i)
                   i+=1
                   table_entries += "{\"HT\":\"%s\","%(dir[0:3])
                   checkht = dir
                   for entry in DirKeys:
                        subdirect = file.FindObjectAny(entry.GetName())
                        sample_dir = fi[1]+dir 
                        subdirect.GetName()
                        if sample_dir == subdirect.GetName():  
                            for subkey in subdirect.GetListOfKeys():
                                if fi[3] =="Had" or fi[3] == "Photon" or (str(dir[0:3]) in ["275","325"] and self.Analysis == "7TeV") or (str(lower) == "0.55" and self.Keep_AlphaT == "True"):
                                    if subkey.GetName() == "Matched_vs_Matched_noB_alphaT_all": 
                                      plot = file.Get(sample_dir+"/"+subkey.GetName())
                                else: 
                                    if subkey.GetName() == "Matched_vs_Matched_noB_all": 
                                      plot = file.Get(sample_dir+"/"+subkey.GetName())
                   if fi[2] != "Data":
                      table_entries += self.Make_Prediction(plot,fi[3],fi[2],number,dir[0:3],lower)
                      
                   else: table_entries += self.Data_Yield(fi[0],fi[1],dir,lower,higher,fi[2],fi[3])
        #if number == "Inclusive" or number == "BaseLine": table_entries += self.photon_baseline_dict()
        #if number == "Zero_btags": table_entries += self.photon_zero_dict()
        #if number == "One_btag": table_entries += self.photon_one_dict()
        #if number == "Two_btags": table_entries += self.photon_two_dict()
        #if number == "More_Than_Two_btag": table_entries += self.photon_more_two_dict()
        table_entries += "}"
        return_dict = ast.literal_eval(table_entries)
        
        return return_dict
 
   
   def photon_one_dict(self):


            if self.Analysis == "7TeV":
                    s = "\"nphot1\":{\"HT\":\"275\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot2\":{\"HT\":\"325\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot3\":{\"HT\":\"375\",\"Yield\":126.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot4\":{\"HT\":\"475\",\"Yield\":43.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot5\":{\"HT\":\"575\",\"Yield\":19.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot6\":{\"HT\":\"675\",\"Yield\":5.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot7\":{\"HT\":\"775\",\"Yield\":5.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot8\":{\"HT\":\"875\",\"Yield\":2.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"

            if self.Analysis == "8TeV":
                    s= "\"phot1\":{\"HT\":\"275\",\"Yield\":0,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"phot2\":{\"HT\":\"325\",\"Yield\":0,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"phot3\":{\"HT\":\"375\",\"Yield\":0,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"phot4\":{\"HT\":\"475\",\"Yield\":16.75,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":2.28},\n"
                    s+= "\"phot5\":{\"HT\":\"575\",\"Yield\":5.77,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":1.31},\n"
                    s+= "\"phot6\":{\"HT\":\"675\",\"Yield\":1.98,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0.59},\n"
                    s+= "\"phot7\":{\"HT\":\"775\",\"Yield\":0.54,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0.27},\n"
                    s+= "\"phot8\":{\"HT\":\"875\",\"Yield\":1.33,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0.7},\n"
                    s+= "\"nphot1\":{\"HT\":\"275\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot2\":{\"HT\":\"325\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot3\":{\"HT\":\"375\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot4\":{\"HT\":\"475\",\"Yield\":29,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot5\":{\"HT\":\"575\",\"Yield\":6.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot6\":{\"HT\":\"675\",\"Yield\":5.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot7\":{\"HT\":\"775\",\"Yield\":1.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot8\":{\"HT\":\"875\",\"Yield\":0.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"

            return s

   def photon_zero_dict(self):
            
            if self.Analysis == "7TeV":
                    s = "\"nphot1\":{\"HT\":\"275\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot2\":{\"HT\":\"325\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot3\":{\"HT\":\"375\",\"Yield\":909.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot4\":{\"HT\":\"475\",\"Yield\":328.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot5\":{\"HT\":\"575\",\"Yield\":109.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot6\":{\"HT\":\"675\",\"Yield\":50.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot7\":{\"HT\":\"775\",\"Yield\":13.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot8\":{\"HT\":\"875\",\"Yield\":12.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
           
            if self.Analysis == "8TeV":
                    s= "\"phot1\":{\"HT\":\"275\",\"Yield\":0,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"phot2\":{\"HT\":\"325\",\"Yield\":0,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"phot3\":{\"HT\":\"375\",\"Yield\":0,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"phot4\":{\"HT\":\"475\",\"Yield\":153.52,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":11.62},\n"
                    s+= "\"phot5\":{\"HT\":\"575\",\"Yield\":55.42,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":8.},\n"
                    s+= "\"phot6\":{\"HT\":\"675\",\"Yield\":20.65,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":4.0},\n"
                    s+= "\"phot7\":{\"HT\":\"775\",\"Yield\":8.09,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":2.5},\n"
                    s+= "\"phot8\":{\"HT\":\"875\",\"Yield\":5.75,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":1.73},\n"
                    s+= "\"nphot1\":{\"HT\":\"275\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot2\":{\"HT\":\"325\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot3\":{\"HT\":\"375\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot4\":{\"HT\":\"475\",\"Yield\":135,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot5\":{\"HT\":\"575\",\"Yield\":64.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot6\":{\"HT\":\"675\",\"Yield\":16.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot7\":{\"HT\":\"775\",\"Yield\":6.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot8\":{\"HT\":\"875\",\"Yield\":3.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
 
            return s

   def photon_two_dict(self):
            if self.Analysis == "7TeV":
                    s = "\"nphot1\":{\"HT\":\"275\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot2\":{\"HT\":\"325\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot3\":{\"HT\":\"375\",\"Yield\":10.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot4\":{\"HT\":\"475\",\"Yield\":4.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot5\":{\"HT\":\"575\",\"Yield\":2.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot6\":{\"HT\":\"675\",\"Yield\":4.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot7\":{\"HT\":\"775\",\"Yield\":0.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot8\":{\"HT\":\"875\",\"Yield\":0.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"

            elif self.Analysis == "8TeV":

                    s= "\"phot1\":{\"HT\":\"275\",\"Yield\":0,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"phot2\":{\"HT\":\"325\",\"Yield\":0,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"phot3\":{\"HT\":\"375\",\"Yield\":0,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"phot4\":{\"HT\":\"475\",\"Yield\":0.84,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0.34},\n"
                    s+= "\"phot5\":{\"HT\":\"575\",\"Yield\":1.76,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0.8},\n"
                    s+= "\"phot6\":{\"HT\":\"675\",\"Yield\":0.03,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0.03},\n"
                    s+= "\"phot7\":{\"HT\":\"775\",\"Yield\":0.003,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0.003},\n"
                    s+= "\"phot8\":{\"HT\":\"875\",\"Yield\":0.3,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0.3},\n"
                    s+= "\"nphot1\":{\"HT\":\"275\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot2\":{\"HT\":\"325\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot3\":{\"HT\":\"375\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot4\":{\"HT\":\"475\",\"Yield\":4,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot5\":{\"HT\":\"575\",\"Yield\":2.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot6\":{\"HT\":\"675\",\"Yield\":0.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot7\":{\"HT\":\"775\",\"Yield\":0.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot8\":{\"HT\":\"875\",\"Yield\":0.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"

            return s

   def photon_more_two_dict(self):
            s = "\"nphot1\":{\"HT\":\"275\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
            s+= "\"nphot2\":{\"HT\":\"325\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
            s+= "\"nphot3\":{\"HT\":\"375\",\"Yield\":1.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
            s+= "\"nphot4\":{\"HT\":\"475\",\"Yield\":0.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
            s+= "\"nphot5\":{\"HT\":\"575\",\"Yield\":0.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
            s+= "\"nphot6\":{\"HT\":\"675\",\"Yield\":0.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
            s+= "\"nphot7\":{\"HT\":\"775\",\"Yield\":0.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
            s+= "\"nphot8\":{\"HT\":\"875\",\"Yield\":0.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"

            return s

   def photon_baseline_dict(self):

            if self.Analysis == "7TeV":
                    s= "\"nphot1\":{\"HT\":\"275\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot2\":{\"HT\":\"325\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot3\":{\"HT\":\"375\",\"Yield\":1046.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot4\":{\"HT\":\"475\",\"Yield\":375.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot5\":{\"HT\":\"575\",\"Yield\":130.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot6\":{\"HT\":\"675\",\"Yield\":59.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot7\":{\"HT\":\"775\",\"Yield\":18.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot8\":{\"HT\":\"875\",\"Yield\":14.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"

            if self.Analysis == "8TeV":
                    s= "\"phot1\":{\"HT\":\"275\",\"Yield\":0,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"phot2\":{\"HT\":\"325\",\"Yield\":0,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"phot3\":{\"HT\":\"375\",\"Yield\":1316.,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":40.},\n"
                    s+= "\"phot4\":{\"HT\":\"475\",\"Yield\":444.,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":20.},\n"
                    s+= "\"phot5\":{\"HT\":\"575\",\"Yield\":180.,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":10.},\n"
                    s+= "\"phot6\":{\"HT\":\"675\",\"Yield\":58.,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":8.0},\n"
                    s+= "\"phot7\":{\"HT\":\"775\",\"Yield\":20.,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":5.0},\n"
                    s+= "\"phot8\":{\"HT\":\"875\",\"Yield\":14.,\"SampleType\":\"Photon\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":3.0},\n"
                    s+= "\"nphot1\":{\"HT\":\"275\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot2\":{\"HT\":\"325\",\"Yield\":0,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot3\":{\"HT\":\"375\",\"Yield\":150.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot4\":{\"HT\":\"475\",\"Yield\":53.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot5\":{\"HT\":\"575\",\"Yield\":18.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot6\":{\"HT\":\"675\",\"Yield\":6.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot7\":{\"HT\":\"775\",\"Yield\":1.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"
                    s+= "\"nphot8\":{\"HT\":\"875\",\"Yield\":0.,\"SampleType\":\"Data\",\"Category\":\"Photon\",\"AlphaT\":0.55,\"Error\":0},\n"

            return s
                    
   def Data_Yield(self,file,dir_path,dir,lower,higher,sample_type,category):
      for histName in self.settings['plots']:
          checkht = dir
          dir = dir_path+dir
          normal =  GetSumHist(File = ["%s.root"%file], Directories = [dir], Hist = histName, Col = r.kBlack, Norm = None, LegendText = "nBtag")  
          normal.HideOverFlow()
          if self.Keep_AlphaT == "True" and str(lower)=="0.55":
            err = r.Double(0.0)
            if category != "Had": normal.hObj.IntegralAndError(int(float(lower)/0.01)+1,int(float(higher)/0.01),err)
            else: normal.hObj.IntegralAndError(int(float(0.55)/0.01)+1,int(float(higher)/0.01),err)
            table_string =" \"Yield\": %.3e ,\"Error\":\"%s\",\"SampleType\":\"%s\",\"Category\":\"%s\",\"AlphaT\":%s},\n"%((normal.hObj.Integral(int(float(0.55)/0.01)+1,int(float(10)/0.01)) if category =="Had" else (normal.hObj.Integral(int(float(lower)/0.01)+1,int(float(higher)/0.01)))),err,sample_type,category,lower)
          else:
            err = r.Double(0.0)
            if category == "Had" or (str(checkht[0:3]) in ["275","325"] and self.Analysis == "7TeV") :normal.hObj.IntegralAndError(int(float(lower)/0.01)+1,int(float(higher)/0.01),err)
            else: normal.hObj.IntegralAndError(1,2000,err)
            table_string =" \"Yield\": %.3e ,\"Error\":\"%s\",\"SampleType\":\"%s\",\"Category\":\"%s\",\"AlphaT\":%s},\n"%((normal.hObj.Integral(int(float(0.55)/0.01)+1,int(float(higher)/0.01)) if category =="Had" or (str(checkht[0:3]) in ["275","325"] and self.Analysis == "7TeV") else (normal.hObj.Integral())),err,sample_type,category,lower)
          normal.a.Close()
      return table_string    

   def Make_Prediction(self,plot,sample,category,btag_number,htbin,alphaT):  
      
        def Pred_Zero(plot,htbin,sample,e,m):
             
          sum_zero = 0
          pred_zero = 0

          for b in range(0,plot.GetNbinsX()):
            for noB in range(0,plot.GetNbinsY()):
                      sum_zero += pow((plot.GetBinError(b+1,noB+1)*pow((1-e),b)*pow((1-m),noB)),2)
          error_zero = sqrt(sum_zero)       

          for b in range (0,plot.GetNbinsX()):
            for noB in range(0,plot.GetNbinsY()):
              if b == 3: pred_zero += pow((1-e),b)*pow((1-m),noB)*plot.GetBinContent(b+1,noB+1)
              else:pred_zero += pow((1-e),b)*pow((1-m),noB)*plot.GetBinContent(b+1,noB+1)


          return (pred_zero,error_zero)    
                          
        def Pred_One(plot,htbin,sample,e,m):

          pred_one = 0
          sum_one = 0
          
          for b in range(1,plot.GetNbinsX()):
            for noB in range(1,plot.GetNbinsY()):
               sum_one += pow((pow((1-e),b-1)*pow((1-m),noB-1)*((noB*m*(1-e))+(noB*e*(1-m)))),2)*pow(plot.GetBinError(b+1,noB+1),2)         
          error_one = sqrt( (m*m*pow(plot.GetBinError(1,2),2))+(e*e*pow(plot.GetBinError(2,1),2))+sum_one)

          for b in range (0,plot.GetNbinsX()):
            for noB in range(1,plot.GetNbinsY()):
              if b == 3:pred_one += pow((1-e),b)*noB*m*pow((1-m),noB-1)*plot.GetBinContent(b+1,noB+1)
              else:pred_one += pow((1-e),b)*noB*m*pow((1-m),noB-1)*plot.GetBinContent(b+1,noB+1)
          
          for b in range(1,plot.GetNbinsX()):
            for noB in range(0,plot.GetNbinsY()):
              if b ==3:pred_one += pow((1-m),noB)*b*e*pow((1-e),b-1)*plot.GetBinContent(b+1,noB+1)
              else:pred_one += pow((1-m),noB)*b*e*pow((1-e),b-1)*plot.GetBinContent(b+1,noB+1)
          
          return (pred_one,error_one)
              
        def Pred_Two(plot,htbin,sample,e,m) :

          pred_two = 0
          sum_two = 0
          for b in range (2,plot.GetNbinsX()):
            for noB in range(2,plot.GetNbinsY()):
              sum_two += pow(plot.GetBinError(b+1,noB+1),2)*((pow((1-e),b)*noB*0.5*(noB-1)*m*m*pow((1-m),noB-2)+pow((1-m),noB)*0.5*b*(b-1)*e*e*pow((1-e),b-2))+(b*e*pow((1-e),b-1)*noB*m*pow((1-m),m-1)))

          error_two = sqrt(m*m*pow(plot.GetBinError(1,3),2)+(pow(((1-e)*m*m)+(2*e*m*(1-m)) ,2)*pow(plot.GetBinError(2,3),2))+pow(e*plot.GetBinError(3,1),2)+pow(e*m*plot.GetBinError(2,2),2)+pow( (((1-m)*e*e)+(2*e*m*(1-m)))  ,2)*pow(plot.GetBinError(3,2),2)+sum_two)


          for b in range (0,plot.GetNbinsX()):
            for noB in range(2,plot.GetNbinsY()):
              #print "Pred_Two = %s" %pred_two
              if b == 3: pred_two += pow((1-e),b)*noB*0.5*(noB-1)*m*m*pow((1-m),noB-2)*plot.GetBinContent(b+1,noB+1)
              else:pred_two += pow((1-e),b)*noB*0.5*(noB-1)*m*m*pow((1-m),noB-2)*plot.GetBinContent(b+1,noB+1)
          
          for b in range (2,plot.GetNbinsX()):
            for noB in range(0,plot.GetNbinsY()):
              if b == 3:pred_two += pow((1-m),noB)*0.5*b*(b-1)*e*e*pow((1-e),b-2)*plot.GetBinContent(b+1,noB+1)
              else: pred_two += pow((1-m),noB)*0.5*b*(b-1)*e*e*pow((1-e),b-2)*plot.GetBinContent(b+1,noB+1)
          
          for b in range (1,plot.GetNbinsX()):
            for noB in range(1,plot.GetNbinsY()):
              if b ==3:pred_two += b*e*pow((1-e),b-1)*noB*m*pow((1-m),m-1)*plot.GetBinContent(b+1,noB+1)
              else:pred_two += b*e*pow((1-e),b-1)*noB*m*pow((1-m),m-1)*plot.GetBinContent(b+1,noB+1)

          return(pred_two,error_two)  
 
        def Pred_Three(plot,htbin,sample,e,m) :

          pred_three = 0

          error_three = (e*e*m)*sqrt(pow(plot.GetBinError(3,2),2 )+(4*(1-m)*(1-m)*pow(plot.GetBinError(3,3),2))+(9*pow((1-m),4)*pow(plot.GetBinError(3,4),2))+(16*pow((1-m),6)*pow(plot.GetBinError(3,5),2))+(25*pow((1-m),8)*pow(plot.GetBinError(3,6),2))     )
  
          for b in range(2,plot.GetNbinsX()):
            for i in range(1,plot.GetNbinsY()):
              if b == 3: pred_three += i*m*pow((1-m),i-1)*b*(b-1)*0.5*e*e*pow((1-e),b-2)*plot.GetBinContent(b+1,i+1)
              else: pred_three += i*m*pow((1-m),i-1)*b*(b-1)*0.5*e*e*pow((1-e),b-2)*plot.GetBinContent(b+1,i+1)
          
          for b in range(3,plot.GetNbinsX()):
            for noB in range (1,plot.GetNbinsY()):
              if b == 3:pred_three += pow((1-m),noB)*b*(b-1)*(b-2)*e*e*e*pow((1-e),b-3)*(1/6)*plot.GetBinContent(b+1,noB+1)
              else:pred_three += pow((1-m),noB)*b*(b-1)*(b-2)*e*e*e*pow((1-e),b-3)*(1/6)*plot.GetBinContent(b+1,noB+1)
          
          for b in range (0,plot.GetNbinsX()):
            for noB in range(3,plot.GetNbinsY()):
              if b ==3:pred_three += pow((1-e),b)*noB*(noB-1)*(noB-2)*(1/6)*m*m*m*pow((1-m),noB-3)*plot.GetBinContent(b+1,noB+1)
              else:pred_three += pow((1-e),b)*noB*(noB-1)*(noB-2)*(1/6)*m*m*m*pow((1-m),noB-3)*plot.GetBinContent(b+1,noB+1)
          
          for b in range (1,plot.GetNbinsX()):
            for noB in range(2,plot.GetNbinsY()):
              if b == 3:pred_three += noB*e*pow((1-e),b-1)*noB*(noB-1)*0.5*m*m*pow((1-m),noB-2)*plot.GetBinContent(b+1,noB+1)
              else:pred_three += noB*e*pow((1-e),b-1)*noB*(noB-1)*0.5*m*m*pow((1-m),noB-2)*plot.GetBinContent(b+1,noB+1)


          return (pred_three,error_three)


        def Pred_Four(plot,htbin,sample,e,m) : 
        
          pred_four = 0
          sum_four = 0

          for s in range (4,plot.GetNbinsY()):
            sum_four += pow(s*(s-1)*(s-2)*(s-3)*(1/24)*pow(m,4)*pow((1-m),s-4),2)*pow(plot.GetBinError(1,s+1),2)
            sum_four += pow((( (1-e)*(s-1)*(s-2)*(s-3)*(1/24)*pow(m,4)*pow((1-m),s-4) ) + (e*s*(s-1)*(s-2)*(1/6)*pow(m,3)*pow((1-m),s-3))  ),2)*pow(plot.GetBinError(2,s+1),2)

            sum_four += pow(((pow((1-e),2)*s*(s-1)*(s-2)*(s-3)*(1/24)*pow(m,4)*pow((1-m),s-4))   + (2*e*(1-e)*s*(s-1)*(s-2)*(1/6)*pow(m,3)*pow((1-m),s-3)) +  (e*e*s*(s-1)*0.5*m*m*pow((1-m),s-2)) ),2)*pow(plot.GetBinError(3,s+1),2)

            sum_four += pow(( (pow((1-e),3)*s*(s-1)*(s-2)*(s-3)*(1/24)*pow(m,4)*pow((1-m),s-4))   + (3*e*pow((1-e),2)*s*(s-1)*(s-2)*(1/6)*pow(m,3)*pow((1-m),s-3))    + (3*e*e*(1-e)*s*(s-1)*0.5*m*m*pow((1-m),s-2))  + e*e*e*s*m*pow((1-m),s-1)    ),2)*pow(plot.GetBinError(4,s+1),2)
         

          for b in range (4,plot.GetNbinsX()):
            sum_four += pow(b*plot.GetBinError(2,s+1)*(b-1)*(b-2)*(b-3)*(1/24)*pow(e,4)*pow((1-e),b-4),2)*pow(plot.GetBinError(b+1,1),2)
            sum_four +=pow(( ( b*(b-1)*(b-2)*(b-3)*(1/24)*pow(e,4)*pow((1-e),b-4)*(1-m)) + ( b*(b-1)*(b-2)*(b-3)*(1/6)*e*e*e*pow((1-e),b-3)*m ) ),2)*pow(plot.GetBinError(b+1,2),2)

          T = pow((plot.GetBinError(2,4)*e*m*m*m),2)+pow((plot.GetBinError(3,3)*e*e*m*m),2)+pow((plot.GetBinError(4,2)*e*e*e*m),2)+pow((plot.GetBinError(3,4)*((2*e*(1-e)*pow(m,3))+(3*e*e*m*m*(1-m)))),2)+pow((plot.GetBinError(4,3)*((3*e*e*(1-e)*m*m)+(2*e*e*e*m*(1-m)) )),2)+pow((plot.GetBinError(4,4)*((3*e*pow((1-e),2)*pow(m,3))+(9*e*e*(1-e)*m*m*(1-m))*(3*e*e*e*m*pow((1-m),2)))),2)
          error_four = sqrt(T+sum_four)
          
          
          for b in range(0,plot.GetNbinsX()):
             for s in range (4,plot.GetNbinsY()):
                pred_four += pow((1-e),b)*(s*(s-1)*(s-2)*(s-3))*(1/24)*pow(m,4)*pow((1-m),s-4)*plot.GetBinContent(b+1,s+1)
          
          for b in range(1,plot.GetNbinsX()):
             for s in range (3,plot.GetNbinsY()):
                pred_four += b*e*pow((1-e),b-1)*(s*(s-1)*(s-2))*(1/6)*pow(m,3)*pow((1-m),s-3)*plot.GetBinContent(b+1,s+1)
    
          for b in range(2,plot.GetNbinsX()):
             for s in range (2,plot.GetNbinsY()):
                pred_four += b*(b-1)*0.5*pow(e,2)*pow((1-e),b-2)*s*(s-1)*0.5*m*m*pow((1-m),s-2)*plot.GetBinContent(b+1,s+1)

          for b in range(3,plot.GetNbinsX()):
             for s in range (1,plot.GetNbinsY()):
                pred_four += b*(b-1)*(b-2)*(1/6)*e*e*e*pow((1-e),b-3)*s*m*pow((1-m),s-1)*plot.GetBinContent(b+1,s+1)

          for b in range(4,plot.GetNbinsX()):
             for s in range (0,plot.GetNbinsY()):
                pred_four += b*(b-1)*(b-2)*(b-3)*(1/24)*pow(e,4)*pow((1-e),b-4)*pow((1-m),s)*plot.GetBinContent(b+1,s+1)
          

          return (pred_four,error_four)


        prediction_dictionary = {"Zero_btags":['Pred_Zero'],"One_btag":['Pred_One'],"Two_btags":["Pred_Two"],"More_Than_Two_btag":['Pred_Three'],"More_Than_Zero_btag":['Pred_One','Pred_Two','Pred_Three'],"More_Than_One_btag":['Pred_Two','Pred_Three'],"Inclusive":['Pred_Zero','Pred_One','Pred_Two','Pred_Three'],"More_Than_Three_btag":['Pred_Four']}
        yield_pred = 0
        error_pred = 0
        for entry in  prediction_dictionary[btag_number]:
           if (category == "Zinv50" or category == "WJetsInc") and sample == "Had":
              btag_eff = self.Btag_Efficiencies['Had_Zinv'][htbin]['Btag_Efficiency'] # +  self.Btag_Efficiencies['Had_Zinv'][htbin]['Btag_Error']
              mistag_eff = self.Btag_Efficiencies['Had_Zinv'][htbin]['Mistag_Efficiency']  #+  self.Btag_Efficiencies['Had_Zinv'][htbin]['Mistag_Error']

           elif category == "Photon":
              btag_eff = self.Btag_Efficiencies['Photon'][htbin]['Btag_Efficiency']  #+  self.Btag_Efficiencies['Photon'][htbin]['Btag_Error']
              mistag_eff = self.Btag_Efficiencies['Photon'][htbin]['Mistag_Efficiency']  #+  self.Btag_Efficiencies['Photon'][htbin]['Mistag_Error']

           elif category == "WJetsInc" and sample == "Muon":
              btag_eff = self.Btag_Efficiencies['Had_Zinv'][htbin]['Btag_Efficiency'] #+  self.Btag_Efficiencies['Had_Zinv'][htbin]['Btag_Error'] 
              mistag_eff = self.Btag_Efficiencies['Had_Zinv'][htbin]['Mistag_Efficiency']  #+  self.Btag_Efficiencies['Had_Zinv'][htbin]['Mistag_Error']
           else: 
              if sample == "DiMuon": btag_eff = self.Btag_Efficiencies['Muon'][htbin]['Btag_Efficiency']  #+  self.Btag_Efficiencies['Muon'][htbin]['Btag_Error']
              else: btag_eff = self.Btag_Efficiencies[sample][htbin]['Btag_Efficiency'] #+ self.Btag_Efficiencies[sample][htbin]['Btag_Error']
              mistag_eff = self.Btag_Efficiencies[sample][htbin]['Mistag_Efficiency'] #+ self.Btag_Efficiencies[sample][htbin]['Mistag_Error']

           calc = eval("%s(plot,htbin,sample,btag_eff,mistag_eff)"%entry)
           yield_pred +=  calc[0]
           error_pred +=  pow(calc[1],2)
        error_pred = sqrt(error_pred)
        Luminosity = self.lumi_dict[sample]
        table_string =" \"Yield\": %.3e ,\"Error\":\"%s\",\"SampleType\":\"%s\",\"Category\":\"%s\",\"AlphaT\":%s},\n"%(yield_pred*(10*Luminosity),error_pred*(Luminosity*10),category,sample,alphaT)
        return table_string
  

   def Btag_Rate(self,btag_measurement):
        self.Btag_Efficiencies = {'Had':{},'Muon':{},'DiMuon':{},'Had_Zinv':{},'Photon':{}}
        self.bins = ('275','325','375','475','575','675','775','875')
        dict_entries = ('Btag_Efficiency','Mistag_Efficiency','Btag_Error','Mistag_Error')
        for key in self.Btag_Efficiencies:
             self.Btag_Efficiencies[key] = dict.fromkeys(self.bins)
             for a in self.bins: self.Btag_Efficiencies[key][a] = dict.fromkeys(dict_entries,0)


        htbins = ["275_325","325_375","375_475","475_575","575_675","675_775","775_875","875"]
        for key,fi in sorted(btag_measurement.iteritems()):
           file = r.TFile.Open(fi[0])
           DirKeys = file.GetListOfKeys()
           
           for num,bin in enumerate(htbins):
             if fi[2] == "Photon_" and bin in ["275_325","325_375"] : continue
             for entry in DirKeys:
                subdirect = file.FindObjectAny(entry.GetName())
                dir = fi[2]+bin
                subdirect.GetName()
                if dir == subdirect.GetName():
                   for subkey in subdirect.GetListOfKeys():
                       #=========================================#
                       if subkey.GetName() == "GenJetPt_all":
                           err = r.Double(0.0)     
                           plot = file.Get(dir+"/"+subkey.GetName())
                           self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Btag_Efficiency"] = plot.Integral()
                           plot.IntegralAndError(1,10000,err)
                           self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Btag_Error"] = err/plot.Integral()


                       if subkey.GetName() == "Btagged_GenJetPt_SFb_all":
                           err = r.Double(0.0)
                           bplot = file.Get(dir+"/"+subkey.GetName())
                           self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]['Btag_Efficiency'] = bplot.Integral()/(self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]['Btag_Efficiency'])
                           bplot.IntegralAndError(1,10000,err)
                           self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Btag_Error"] =  self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]['Btag_Efficiency']*sqrt(pow(self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Btag_Error"],2)+pow(err/bplot.Integral(),2))


                       #========================================#
                       if subkey.GetName() == "GenJetPt_noB_all":
                           mistag_plot = file.Get(dir+"/"+subkey.GetName())
                           err = r.Double(0.0)
                           self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]['Mistag_Efficiency'] =mistag_plot.Integral()
                           plot.IntegralAndError(1,10000,err)
                           self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Mistag_Error"] = err/mistag_plot.Integral()

                       if subkey.GetName() == "Btagged_GenJetPt_noB_SFlight_all":
                           aplot = file.Get(dir+"/"+subkey.GetName())
                           err = r.Double(0.0)
                           self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]['Mistag_Efficiency'] = aplot.Integral()/(self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]['Mistag_Efficiency'])
                           aplot.IntegralAndError(1,10000,err)
                           self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Mistag_Error"] =  self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]['Mistag_Efficiency']*sqrt(pow(self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Mistag_Error"],2)+pow(err/aplot.Integral(),2))


                        #=======================================#
 
