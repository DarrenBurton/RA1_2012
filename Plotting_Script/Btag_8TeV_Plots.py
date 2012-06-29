#!/usr/bin/env python
from ROOT import *
import ROOT as r
import logging,itertools
import os,fnmatch,sys
import glob, errno
from time import strftime
from optparse import OptionParser
import array
import math as m

class Plotter(object):

  def __init__(self,settings,sample_list,jet_multiplicity = ""):
    #================================== Preamble
    print " Selecting tdr style"
    r.gROOT.ProcessLine(".L tdrstyle.C")
    r.setstyle()
    r.gROOT.SetBatch(True)
    r.gStyle.SetOptStat(0)
    #==================================
    self.settings = settings
    self.jet_multi = jet_multiplicity
    self.sample_list = sample_list
    # Apply options
    self.splash_screen()
    self.Hist_Getter(settings,sample_list)
    self.Plotting_Option(settings,sample_list)
    

  def splash_screen(self):
    print "\n|============================================================================|"
    print "| XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX |"
    print "|==================  STARTING  BINNED ALPHA T PLOTTING ======================|"
    print "| XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX |"
    print "|============================================================================|"

  def ensure_dir(self,dir):
        try: os.makedirs(dir)
        except OSError as exc: pass

  def Hist_Getter(self,settings,sample_list):
       for key,sample in sorted(sample_list.iteritems()):
          if "n" == key[0]:
             self.DATA_FILE = sample[0]
             self.Directory_Name = sample[1]
             temp =  r.TFile.Open(self.DATA_FILE)
             DirKeys = temp.GetListOfKeys()
             self.Path_List =[]
             self.Hist_List = []
             for key in DirKeys:
               subdirect = temp.FindObjectAny(key.GetName())
               for bin in settings["dirs"]:
                dir = sample[1]+bin
                if dir ==  subdirect.GetName():
                  for subkey in subdirect.GetListOfKeys() :
                    if subkey.GetName() in settings["Plots"]:
                      if self.jet_multi == "True":
                       base = subkey.GetName().strip('all')
                       for entry in ['all','2','3']:
                        self.Path_List.append("%s/%s" % (subdirect.GetName(),base+entry))
                        self.Hist_List.append(base+entry)
                      else: 
                        self.Path_List.append("%s/%s" % (subdirect.GetName(),subkey.GetName()))
                        self.Hist_List.append(subkey.GetName())
             temp.Close()
  
  def Directory_Maker(self):
    htBins = self.settings["dirs"]
    self.Dir_Binning = []
    print "\n Making Directory ::: Plots :::" 
    self.ensure_dir("Plots")
    self.base = os.getcwd()
    os.chdir("Plots")
    owd = os.getcwd()
    for bin in htBins: self.Dir_Binning.append(bin)
    for path in self.Dir_Binning: 
      self.ensure_dir(path)
      os.chdir(path)
      for h in self.Hist_List: self.ensure_dir(str(h))
      os.chdir(owd)

  def Plotting_Option(self,settings,sample_list):
        self.Directory_Maker()
        for num,histpath in enumerate(self.Path_List):
          name = self.Hist_List[num]
          ht = (histpath.split('/')[0]).strip(self.Directory_Name)
          if ht =="275_325" or ht =="325_375" or ht =="375_475":self.Make_Plots(ht,name,histpath)
          if ht == "375_475": self.Make_Plots(ht,name,histpath,combine = "True",histpath= (histpath.split('/')[0]).rstrip("375_475")) 
        os.chdir(self.base)
 
  def MakePad(self,plot):
      if type(plot) != type(TH2D()):
        # Make 2 pads for the Data/MC ratio plot fit.
        self.up = r.TPad("u","u",0.01,0.25,0.99,0.99)
        self.dp = r.TPad("d","d",0.01,0.01,0.99,0.25)
        self.up.SetNumber(1)
        self.dp.SetNumber(2)
        self.up.Draw()
        self.dp.Draw()

  def MakeMCRatio(self,histname,data,mc):
        if type(data) != type(TH2D()):
          self.c1.cd(2)
          self.ratio = data.Clone()
          #self.ratio.Add(mc,-1)
          self.ratio.Divide(mc)
          self.ratio.GetYaxis().SetRangeUser(0,2)
          self.ratio.GetYaxis().SetTitle("Data/MC")
          self.ratio.GetYaxis().SetTitleOffset(0.5)
          self.ratio.GetYaxis().SetLabelSize(0.11)
          self.ratio.GetYaxis().SetNdivisions(505)
          self.ratio.GetXaxis().SetTitle("")
          self.ratio.GetXaxis().SetLabelSize(0.12)
          self.ratio.SetTitleSize(0.1,"Y")
          for i in range(0,self.ratio.GetNbinsX()):
              try:self.ratio.SetBinError(i+1,self.ratio.GetBinContent(i+1)*(data.GetBinError(i+1)/data.GetBinContent(i+1)))
              except ZeroDivisionError: self.ratio.SetBinError(i+1,0)
          self.Hist_Options(histname,self.ratio,norebin=true)
          if self.max_min[1] != 0: 
            self.lv = r.TLine(self.max_min[0],1,self.max_min[1],1)
            fit = r.TF1("fit","pol0", self.max_min[0], self.max_min[1])
          else:
            fit = r.TF1("fit","pol0", data.GetXaxis().GetBinLowEdge(1), data.GetXaxis().GetBinUpEdge(data.GetNbinsX()))
            self.lv = r.TLine(data.GetXaxis().GetBinLowEdge(1),1,data.GetXaxis().GetBinUpEdge(data.GetNbinsX()),1)
          error_changer = self.ratio.Clone()
          for i in range(0,error_changer.GetNbinsX()):
            error_changer.SetBinContent(i,1)
            try: error_changer.SetBinError(i+1,mc.GetBinError(i+1)/mc.GetBinContent(i+1))
            except ZeroDivisionError: error_changer.SetBinError(i+1,0)
          self.ratio_error = TGraphErrors(error_changer)
          #self.ratio.Fit(fit)
          self.ratio_error.SetFillStyle(3244)
          self.ratio_error.SetFillColor(17)
          #self.bv.SetFillColor(39)
          #self.bv.SetFillStyle(3002)
          self.lv.SetLineWidth(3)
          self.ratio.Draw("p")
          self.ratio_error.Draw("SAME2")
          #self.bv.Draw("SAME")
          self.lv.Draw("SAME")
          self.ratio.Draw("PSAME")

         
  def Poission(self,plot):
    poisson_eh = [ 1.15, 1.36, 1.53, 1.73, 1.98, 2.21, 2.42, 2.61, 2.80, 3.00, 3.16 ]
    poisson_el = [ 0.00, 1.00, 2.00, 2.14, 2.30, 2.49, 2.68, 2.86, 3.03, 3.19, 3.16 ]

    for i in range(1,plot.GetNbinsX()):
      x = plot.GetBinContent(i)
      if x<10:
        # Apply Poission errors
        n = int(x)
        plot.SetBinError(i,poisson_el[n])
    return plot 

  def Make_Plots(self,htbin,histname,rootpath,combine = "",histpath=""):
      print "At Histogram %s %s" %(rootpath,histname)
      self.DataFile =  r.TFile.Open("../%s" %self.DATA_FILE)
      self.Plot_Closer = [self.DataFile]
      plot  = self.DataFile.Get("%s" %rootpath)
      self.c1 = r.TCanvas("canvas"+str(rootpath),"canname"+str(histname),1200,1200)
      #if histname.split('_')[-1] == '3':
      #  for jet_mult in range (4,16):
      #    add_jet_mult_plot = self.DataFile.Get(("%s/%s" %(rootpath.split('/')[0],histname.rstrip('3')+str(jet_mult))))
      #    plot.Add(add_jet_mult_plot,1)
      if combine == "True":
        for bin in self.settings["dirs"][3:]:
          #if histname.split('_')[-1] != '3':
          add_plot = self.DataFile.Get(("%s_%s/%s" %(histpath,bin,histname) if histpath !="" else "%s%s/%s" %(histpath,bin,histname) ))
          #else:
          #  for jet_mult in range(4,16):
          #    if len(histpath) == 0: add_sub_plot = self.DataFile.Get("%s/%s" %(bin,histname.rstrip('3')+str(jet_mult)))
          #    else: add_sub_plot = self.DataFile.Get("%s_%s/%s" %(histpath,bin,histname.rstrip('3')+str(jet_mult)))
          #    if jet_mult == 4: add_plot = add_sub_plot.Clone()
          #    else: add_plot.Add(add_sub_plot)
          plot.Add(add_plot,1)
      else: plot  = self.DataFile.Get("%s" %rootpath)
      self.MakePad(plot)
      self.max_min = [0,0]
      self.c1.cd(1)
      plot.GetSumw2()
      self.iflog = 0
      self.reversed = 0
      self.histclone_keeper = []
      self.total_mc_maker = []
      self.ewk_mc_maker = []
      self.qcd_only = []
      self.max_maker = [plot]
      self.draw_error = []
      plot.SetMarkerStyle(20)
      plot.SetMarkerSize(1.5)
      plot.SetLineWidth(2)
      self.Hist_Options(histname,plot,canvas = self.up,word=True)
      self.Legend_Maker()
      self.leg.AddEntry(plot,"Data","P")
      self.Poission(plot)
      self.MC_Draw(rootpath,histname,htbin, histpath = histpath,combine = combine)
      self.stackleg = self.leg.Clone()
      self.Draw_Total_MC(self.total_mc_maker)
      self.Draw_Total_MC(self.ewk_mc_maker, NoDraw = "True")
      self.yaxis_maximiser(plot)
      self.Drawer(self.max_maker)
      self.Drawer(self.draw_error,error="True")
      self.leg.Draw("SAME")
      self.TextBox(plot,(rootpath.split("/")[0] if not combine else (rootpath.split("/")[0]).rstrip("475")),histname)
      self.TagBox()
      plot.Draw("PSAME")
      self.MakeMCRatio(histname,plot,self.totmcplot)
      if combine == "True":
        for ext in ['png','pdf','C']: self.c1.SaveAs("%s/%s/%s_%s_upwards.%s" %(htbin,histname,histname,(rootpath.split("/")[0]).rsplit('_',1)[0],ext))
        self.Stack_Draw(plot,self.total_mc_maker,htbin,histname,rootpath,combine = combine)
        self.Simple_Draw(plot,self.ewk_background,self.qcd_only,htbin,histname,rootpath,combine = combine)
        if "Normalise" in self.settings["Misc"]: self.Normalise_Plots(plot,self.total_mc_maker,htbin,histname,rootpath,combine = combine)

      else:
        for ext in ['png','pdf','C']: self.c1.SaveAs("%s/%s/%s_%s.%s" %(htbin,histname,histname,rootpath.split("/")[0],ext))
        self.Stack_Draw(plot,self.total_mc_maker,htbin,histname,rootpath)
        self.Simple_Draw(plot,self.ewk_background,self.qcd_only,htbin,histname,rootpath)
        if "Normalise" in self.settings["Misc"]: self.Normalise_Plots(plot,self.total_mc_maker,htbin,histname,rootpath)
      if self.reversed ==1:
        self.ReversedDrawer(histname,self.c1,self.histclone_keeper,rootpath,combined = combine)
        if combine == "True": self.c1.SaveAs("%s/%s/%s_%s_upwards_reversed.png" %(htbin,histname,histname,(rootpath.split("/")[0]).rsplit('_',1)[0]))
        else: self.c1.SaveAs("%s/%s/%s_%s_reversed.png" %(htbin,histname,histname,rootpath.split("/")[0]))
      
      for file in self.Plot_Closer: file.Close()

  def Normalise_Plots(self,plot,mc_array,htbin,histname,rootpath,combine=""):

    self.c1.cd(1) 
    for num,mc_hist in enumerate(mc_array):
      if num == 0: mc_plot = mc_hist
      else: mc_plot.Add(mc_hist)
      print num
    mc_plot.GetSumw2()  
    mc_int = mc_plot.Integral()
    plot_int = plot.Integral()
    
    for i in range (plot.GetNbinsX()): 
      try:plot.SetBinContent(i+1,plot.GetBinContent(i+1)/plot_int)
      except ZeroDivisionError: plot.SetBinContent(i+1,0)
      try:mc_plot.SetBinContent(i+1,mc_plot.GetBinContent(i+1)/mc_int)
      except ZeroDivisionError: mc_plot.SetBinContent(i+1,0)


    leg = r.TLegend(0.68,0.53,0.90,0.75) 
    leg.SetTextSize(0.02)
    leg.SetShadowColor(0)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    leg.AddEntry(plot,"8TeV Data","L")
    leg.AddEntry(mc_plot,"7TeV Data","L")


    mc_plot.SetLineColor(3)
    mc_plot.SetFillColor(0)
    mc_plot.GetYaxis().SetRangeUser(0.01,1)
    mc_plot.Draw("HIST")
    plot.SetLineColor(2)
    plot.Draw("HISTSAME")
    leg.Draw("SAME")
    self.TextBox(plot,(rootpath.split("/")[0] if not combine else (rootpath.split("/")[0]).rstrip("475")),histname)
    self.TagBox()
    self.MakeMCRatio(histname,plot,mc_plot)

    if combine == "True":
      for ext in ['png','pdf','C']: self.c1.SaveAs("%s/%s/Simplified_%s_%s_upwards.%s" %(htbin,histname,histname,(rootpath.split("/")[0]).rsplit('_',1)[0],ext))
    else:
      for ext in ['png','pdf','C']: self.c1.SaveAs("%s/%s/Simplified_%s_%s.%s" %(htbin,histname,histname,rootpath.split("/")[0],ext))  


  def Stack_Draw(self,data,mc_array,htbin,histname,rootpath,combine=""):
    self.c1.cd(1)
    mc_stack = THStack()
    self.Stack_Sorter(mc_array)
    for plot in mc_array:
      plot.SetFillColor(plot.GetLineColor())
      mc_stack.Add(plot)
    data.Draw()
    mc_stack.Draw("HISTSAME")
    self.stackleg.Draw("SAME")
    data.Draw("PSAME")
    self.TextBox(data,(rootpath.split("/")[0] if not combine else (rootpath.split("/")[0]).rstrip("475")),histname)
    self.TagBox()
    self.MakeMCRatio(histname,data,self.totmcplot)
    if combine == "True":
      for ext in ['png','pdf','C']: self.c1.SaveAs("%s/%s/Stacked_%s_%s_upwards.%s" %(htbin,histname,histname,(rootpath.split("/")[0]).rsplit('_',1)[0],ext))
    else:
      for ext in ['png','pdf','C']: self.c1.SaveAs("%s/%s/Stacked_%s_%s.%s" %(htbin,histname,histname,rootpath.split("/")[0],ext))  

  def Stack_Sorter(self,mc_array):
    temp = []
    for num,plot in enumerate(mc_array):
        temp.append(plot.Integral())
    swapped = True  
    while swapped:  
       swapped = False  
       for i in range(len(temp)-1):  
          if temp[i] > temp[i+1]:  
             temp[i], temp[i+1] = temp[i+1], temp[i]
             mc_array[i+1], mc_array[i] = mc_array[i], mc_array[i+1]
             swapped = True 
    return mc_array
  
  def Simple_Draw(self,data,mc_combined,qcd_only,htbin,histname,rootpath,combine = ''):
    self.c1.cd(1)
    self.Legend_Maker()
    self.leg.AddEntry(data,"Data","P")
    #self.leg.AddEntry(self.simple_draw_zinv,"Zinv","L")
    self.leg.AddEntry(mc_combined,"Combined EWK","L")
    if qcd_only: self.leg.AddEntry(qcd_only[0],"QCD","L")
    data.Draw()
    mc_combined.SetMarkerColor(3)
    mc_combined.SetMarkerStyle(20)
    #self.simple_draw_zinv.SetMarkerStyle(20)
    #self.simple_draw_zinv.Draw("SAMEHIST")
    mc_stacker = THStack()
    mc_stacker.Add(mc_combined)
    if qcd_only: 
      qcd_only[0].SetFillColor(0)
      mc_stacker.Add(qcd_only[0])
    mc_stacker.Draw("SAMEHIST")
    self.leg.Draw("SAME")
    self.TextBox(data,(rootpath.split("/")[0] if not combine else (rootpath.split("/")[0]).rstrip("475")),histname)
    self.TagBox()
    data.Draw("PSAME")
    self.MakeMCRatio(histname,data,mc_combined)
    if combine == "True": 
        for ext in ['png','pdf','C']: self.c1.SaveAs("%s/%s/Simplified_%s_%s_upwards.%s" %(htbin,histname,histname,(rootpath.split("/")[0]).rsplit('_',1)[0],ext))
    else:
        for ext in ['png','pdf','C']: self.c1.SaveAs("%s/%s/Simplified_%s_%s.%s" %(htbin,histname,histname,rootpath.split("/")[0],ext))
    
  
  def Draw_Total_MC(self,mc_list,NoDraw = ""):
    if NoDraw == "True":
      for num,hist in enumerate(mc_list):
        if num == 0:self.ewk_background = hist.Clone()
        else: self.ewk_background.Add(hist)
      self.ewk_background.SetLineColor(3)
      self.ewk_background.SetLineWidth(3)

    else:
      for num,hist in enumerate(mc_list):
        if num == 0:self.totmcplot = hist.Clone()
        else: self.totmcplot.Add(hist)
      errorbarplot = TGraphErrors(self.totmcplot)
      errorbarplot.SetFillColor(3)
      errorbarplot.SetFillStyle(3008)
      self.draw_error.append(errorbarplot)
      self.totmcplot.SetLineColor(3)
      self.leg.AddEntry(self.totmcplot,"Combined SM","L")
      self.totmcplot.SetLineWidth(3)
      self.max_maker.append(self.totmcplot)

  def Drawer(self,hist_collection,error=''):
    for num,plot in enumerate(hist_collection):
      if error:
        plot.Draw("SAME2")
      else: 
        if num == 0:
          if type(plot) == type(TH2D()): plot.Draw("P")
          else : plot.Draw("PE0") 
        else: plot.Draw("HISTSAME")

  def ReversedDrawer(self,histname,canvas,hist_collection,rootpath):
    for num,plot in enumerate(hist_collection):
      if num == 0 :self.Hist_Options(histname,plot,canvas = canvas,norebin=true)
      else: self.Hist_Options(histname,plot,norebin=true)
      self.yaxis_maximiser(hist_collection[0])
      self.Drawer(hist_collection)
      self.leg.Draw("SAME")
      self.TextBox(hist_collection[0],rootpath.split("/")[0])
      self.TagBox()
  
  def yaxis_maximiser(self,plot):
      highest = 0
      if type(plot) != type(TH2D()):
        for max in self.max_maker:
          temp = max.GetMaximum()
          if temp > highest: highest = temp
        if self.iflog and "Normalise" not in self.settings["Misc"]: self.ymax = highest * 10
        else: self.ymax = highest * 1.5
        self.max_maker[0].GetYaxis().SetRangeUser(self.iflog,self.ymax)

  def TextBox(self,plot,htbin,histname):
    Textbox = TLatex()
    Textbox.SetNDC()
    Textbox.SetTextAlign(12)
    Textbox.SetTextSize(0.04)
    Textbox.DrawLatex(0.1,0.95, htbin+'    Jet Multiplicity'+('>=' if histname.split('_')[-1] == '3' else '=')+histname.split('_')[-1])

  def TagBox(self):
    Textbox = TLatex()
    title = "CMS Preliminary 2012, 8TeV"
    lumi= "\int L dt = %s fb^{-1}" % (self.settings['Lumo']/10)
    Textbox.SetNDC()
    Textbox.SetTextAlign(12)
    Textbox.SetTextSize(0.03)
    Textbox.DrawLatex(0.6,0.85,title )
    Textbox.DrawLatex(0.6,0.79,lumi)

  def Plot_Combiner(self,passed_plot,htbins,histpath,histname,File,add_jet_mult = ''):
    
    for bin in htbins:
      #if add_jet_mult == "False": 
      another_plot = File.Get(("%s_%s/%s" %(histpath,bin,histname) if histpath !="" else "%s%s/%s" %(histpath,bin,histname) ) )        
      #else:
      #    for jet_mult in range(4,16):
      #        add_sub_plot = File.Get(("%s_%s/%s" %(histpath,bin,histname.rstrip('3')+str(jet_mult)) if histpath !="" else "%s%s/%s" %(histpath,bin,histname.rstrip('3')+str(jet_mult))))
      #        if jet_mult == 4: another_plot = add_sub_plot.Clone()
      #        else:another_plot.Add(add_sub_plot)
      passed_plot.Add(another_plot,1)

  def MC_Draw(self,rootpath,histname,htbin,histpath,combine = ""):
      for key,sample in sorted(self.sample_list.iteritems()):
        if "n" != key[0]:
          if sample[2] == "WJets": 
            self.WJetsFile = r.TFile.Open("../%s" %sample[0])
            self.Add_MCplot("wjet_plot",rootpath,histname,htbin,4,"W + Jets",self.WJetsFile,histpath,combine=combine) 
          if sample[2] == "Zinv":
            self.ZinvFile = r.TFile.Open("../%s" %sample[0]) 
            self.Add_MCplot("zinv_plot",rootpath,histname,htbin,5,"Z\\rightarrow \\nu\\bar{\\nu}",self.ZinvFile,histpath,combine=combine) 
          if sample[2] == "TTbar": 
            self.TTbarFile = r.TFile.Open("../%s" %sample[0])
            self.Add_MCplot("ttbar_plot",rootpath,histname,htbin,2,"t\\bar{t}",self.TTbarFile,histpath,combine=combine)           
          if sample[2] == "Top": 
            self.TopFile = r.TFile.Open("../%s" %sample[0])
            self.Add_MCplot("top_plot",rootpath,histname,htbin,2,"t\\bar{t}/ Single t",self.TopFile,histpath,combine=combine)  
          if sample[2] == "Di-Boson": 
            self.Di_Boson = r.TFile.Open("../%s" %sample[0])
            self.Add_MCplot("di_boson_plot",rootpath,histname,htbin,44,"WW/ZZ/WZ",self.Di_Boson,histpath,combine=combine)#,style="9")
          if sample[2] == "DY": 
            self.DiMuonFile = r.TFile.Open("../%s" %sample[0])
            self.Add_MCplot("dimuon_plot",rootpath,histname,htbin,6,"DY + Jets",self.DiMuonFile,histpath,combine=combine)
          if sample[2] == "Single_Top": 
            self.SingleTopFile = r.TFile.Open("../%s" %sample[0]) 
            self.Add_MCplot("sing_top_plot",rootpath,histname,htbin,40,"Single Top",self.SingleTopFile,histpath,combine=combine)#,style="9")
          if sample[2] == "Photon": 
            self.PhotonFile = r.TFile.Open("../%s" %sample[0])
            self.Add_MCplot("photon_plot",rootpath,histname,htbin,8,"\gamma + Jets",self.PhotonFile,histpath,combine=combine)
          if sample[2] == "QCD": 
            self.QCDFile = r.TFile.Open("../%s" %sample[0])
            self.Add_MCplot("qcd_plot",rootpath,histname,htbin,7,"QCD",self.QCDFile,histpath,combine=combine,style="10")
          
  def Add_MCplot(self,mcplot,rootpath,histname,htbin,color,leg_entry,File,histpath,combine = "",style=""):
      mcplot  = File.Get("%s" %rootpath)
      #if histname.split('_')[-1] == '3':
        #for jet_mult in range (4,16):
          #add_jet_mult_mcplot = File.Get(("%s/%s" %(rootpath.split('/')[0],histname.rstrip('3')+str(jet_mult))))
          #mcplot.Add(add_jet_mult_mcplot,1)
      if combine: self.Plot_Combiner(mcplot,self.settings["dirs"][3:],histpath,histname,File,add_jet_mult = ("True" if histname.split('_')[-1] == '3' else "False")) 
      mcplot.GetSumw2()
      mcplot.Scale(float(self.settings["Lumo"]))
      mcplot.SetLineColor(int(color))
      self.leg.AddEntry(mcplot,str(leg_entry),"L")
      mcplot.SetLineWidth(3)
      #Temp Scale Factor for Trigger Efficiency
      mcplot.Scale(self.settings["Trigger"][htbin.split('_')[0]])
      if style: mcplot.SetLineStyle(int(style))
      self.Hist_Options(histname,mcplot)
      if leg_entry =="Z\\rightarrow \\nu\\bar{\\nu}":
        self.simple_draw_zinv = mcplot.Clone()
      if leg_entry == "QCD": 
        errorbarplot = TGraphErrors(mcplot)
        errorbarplot.SetFillColor(int(color))
        errorbarplot.SetFillStyle(3008)
        self.qcd_only.append(mcplot)
        if style: 
          errorbarplot.SetLineStyle(int(style))
          self.draw_error.append(errorbarplot)
      self.max_maker.append(mcplot)
      self.total_mc_maker.append(mcplot)
      if leg_entry != "QCD": self.ewk_mc_maker.append(mcplot)
      self.Plot_Closer.append(File)

  def Legend_Maker(self): 
      self.leg = r.TLegend(0.68,0.53,0.90,0.75) 
      self.leg.SetTextSize(0.02)
      self.leg.SetShadowColor(0)
      self.leg.SetBorderSize(0)
      self.leg.SetFillColor(0)
      self.leg.SetLineColor(0)

  def Hist_Options(self,histogram,plot,canvas="",word = "",norebin=""):
 
        if word: print "Applying %s Options" % histogram
       
        if "EffectiveMass_after_alphaT_55" in histogram:
          if canvas:self.Log_Setter(plot,canvas,0.5)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events / 50 GeV")
          if not norebin:
            plot.Rebin(50)
            self.OverFlow_Bin(plot,0,3000,1500)
            #self.Reversed_Integrator(plot) 
         
        if histogram == "EffectiveMass_all" or histogram == "EffectiveMass_2" or histogram == "EffectiveMass_3":
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events / 50 GeV")
          if not norebin:
            plot.Rebin(50)
            self.OverFlow_Bin(plot,0,2500,1600)
            #self.Reversed_Integrator(plot)  

        if histogram == "MHT_all" or histogram == "MHT_2" or histogram == "MHT_3":
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events / 50 GeV")
          if not norebin:
            plot.Rebin(50)
            self.OverFlow_Bin(plot,0,600,500)
        
        if histogram == "MT__all" or histogram == "MT__2" or histogram == "MT__3":
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events / 50 GeV")
          if not norebin:
            plot.Rebin(50)
            self.OverFlow_Bin(plot,0,2000,800)

        if "MuCso_" in histogram:
          if canvas: self.Log_Setter(plot,canvas,0.1)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events / 0.02")
          if not norebin:
            plot.Rebin(20)
            self.OverFlow_Bin(plot,0,2.,0.2)
            #self.Reversed_Integrator(plot) 

        if "MuPFIso_" in histogram:
          if canvas: self.Log_Setter(plot,canvas,0.1)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events / 0.02")
          if not norebin:
            plot.Rebin(20)
            self.OverFlow_Bin(plot,0,2.,0.16)
            #self.Reversed_Integrator(plot)  



        if "MuHIso_" in histogram:
          if canvas: self.Log_Setter(plot,canvas,0.1)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events / 0.02")
          if not norebin:
            plot.Rebin(20)
            self.OverFlow_Bin(plot,0,2.,0.14)
            #self.Reversed_Integrator(plot)  

        if "MuEIso_" in histogram:
          if canvas: self.Log_Setter(plot,canvas,0.1)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events / 0.02")
          if not norebin:
            plot.Rebin(20)
            self.OverFlow_Bin(plot,0,2.,0.14)
            #self.Reversed_Integrator(plot)  

        if "MuTrIso_" in histogram:
          if canvas: self.Log_Setter(plot,canvas,0.1)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events / 0.02")
          if not norebin:
            plot.Rebin(20)
            self.OverFlow_Bin(plot,0,2.,0.2)
            #self.Reversed_Integrator(plot)  

        if "MuPt_" in histogram:
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events / 25 GeV")
          if not norebin:
            plot.Rebin(20)
          self.OverFlow_Bin(plot,10,2010,1000)

        if histogram == "HT_all" or histogram == "HT_2" or histogram == "HT_3":
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events / 50 GeV")
          if not norebin:
            plot.Rebin(50)  
          self.OverFlow_Bin(plot,0,2000,1250)

        if "HT_after_alphaT_55_" in histogram:
          if canvas:self.Log_Setter(plot,canvas,0.5)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events / 50 GeV")
          if not norebin:
            plot.Rebin(50)
            self.OverFlow_Bin(plot,0,2000,1250)
        
        if histogram == "AlphaT_all" or histogram == "AlphaT_2" or histogram == "AlphaT_3":  
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events / 0.1")
          if not norebin: 
            plot.Rebin(10)
            self.OverFlow_Bin(plot,0.1,10,3.0)
            #self.Reversed_Integrator(plot)

        if "AlphaT_Zero_" in histogram:  
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events / 0.1")
          if not norebin: 
            plot.Rebin(10)
            self.OverFlow_Bin(plot,0.0,10,3.0)
            #self.Reversed_Integrator(plot)

        if "AlphaT_Seven_" in histogram:  
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events / 0.1")
          if not norebin: 
            plot.Rebin(10)
            self.OverFlow_Bin(plot,0.0,10,3.0)
            #self.Reversed_Integrator(plot)

        if "AlphaT_Twelve_" in histogram:  
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events / 0.1")
          if not norebin: 
            plot.Rebin(10)
            self.OverFlow_Bin(plot,0.0,10,3.0)
            #self.Reversed_Integrator(plot)

        if histogram == "AlphaT_Zoomed_all" or histogram == "AlphaT_Zoomed_2" or histogram == "AlphaT_Zoomed_3":
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events / 0.08")
          if not norebin:
            plot.Rebin(5)
            #self.OverFlow_Bin(plot,0.45,0.6,0.6)
            #self.Reversed_Integrator(plot) 

        if histogram == "JetMultiplicity_all" or histogram == "JetMultiplicity_2" or histogram == "JetMultiplicity_3": 
          if canvas: self.Log_Setter(plot,canvas,0.5)
        
        if "JetMultiplicityAfterAlphaT_55_" in histogram: 
          if canvas: self.Log_Setter(plot,canvas,0.5)

        if "BiasedDeltaPhi_after_alphaT_55_" in histogram:
           if canvas: self.Log_Setter(plot,canvas,0.5)
           if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events / 0.25 rad")
           plot.Rebin(5)

        if histogram == "Btag_Pre_AlphaT_4_all" or histogram == "Btag_Pre_AlphaT_4_2" or histogram == "Btag_Pre_AlphaT_4_3": 
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if word: 
              plot.GetXaxis().SetTitleOffset(1.3)
              plot.GetXaxis().SetTitle("N_{b-tag}")
 
        if "Btag_Post_AlphaT_4_55_" in histogram: 
          if canvas: self.Log_Setter(plot,canvas,0.5) 
          if word: 
            plot.GetXaxis().SetTitleOffset(1.3)
            plot.GetXaxis().SetTitle("N_{b-tag}")

        if histogram == "Btag_Pre_AlphaT_5__all" or  histogram == "Btag_Pre_AlphaT_5__2" or histogram == "Btag_Pre_AlphaT_5__3": 
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if word: 
            plot.GetXaxis().SetTitleOffset(1.3)
            plot.GetXaxis().SetTitle("N_{b-tag}")
 
        if "Btag_Post_AlphaT_5_55_" in histogram: 
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if word: 
            plot.GetXaxis().SetTitleOffset(1.3)
            plot.GetXaxis().SetTitle("N_{b-tag}")
        
        if "Number_Primary_verticies_after_alphaT_55_" in histogram:
          if canvas: self.Log_Setter(plot,canvas,0.5)

        if "Number_Primary_verticies_" in histogram:
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if not norebin:
            pass
            #print "Rebinning Primary Verticies"
            #plot.Rebin(5)

        if "Total_Number_Primary_verticies_" in histogram:
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if not norebin:
            print "Rebinning Total Primary Verticies"
            #plot.Rebin(5)

        if "MHTOvCaloMET_TypeI_" in histogram or "MHTOvCaloMET_Uncorrected_" in histogram or "MHTOvMET_PF_" in histogram:
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if not norebin:
             plot.Rebin(2)
             self.OverFlow_Bin(plot,0.50,50.0,5.0)

        if histogram == "CaloMET_TypeI_all" or histogram == "CaloMET_Uncorrected_all" or histogram == "MET_PF_all" or histogram == "CaloMET_TypeI_afterAlphaT_55_all" or histogram == "CaloMET_Uncorrected_afterAlphaT_55_all" or histogram == "MET_PF_afterAlphaT_55_all" or histogram == "MET_PF_TypeI_afterAlphaT_55_all" or histogram == "MET_PF_TypeI_all":
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if not norebin:
             plot.Rebin(50)
             self.OverFlow_Bin(plot,0.0,1500,800)

        if "PhotonPt_" in histogram:
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events / 25 GeV")
          if not norebin:
            plot.Rebin(25)
          self.OverFlow_Bin(plot,0,1500,650)

        if "PhotonIso_" in histogram:
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events")
          if not norebin:
            plot.Rebin(5)
          self.OverFlow_Bin(plot,0,10,8)

        if "Photonrho25_" in histogram:
          if canvas: self.Log_Setter(plot,canvas,0.5)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events")
          if not norebin:
            plot.Rebin(4)
          self.OverFlow_Bin(plot,0,40,32)

        if "Muon_Iso_PU_Interactions_" in histogram:
          #if canvas: self.Log_Setter(plot,canvas,0.5)
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)
            plot.GetYaxis().SetTitle("Events")
          #if not norebin:
          #  plot.Rebin(2)
          self.OverFlow_Bin(plot,0,80,45)




  def Log_Setter(self,plot,canvas,min):
      self.iflog = min
      plot.SetMinimum(float(min))
      canvas.SetLogy(1)

  def OverFlow_Bin(self,hist,xmin,xmax,overflow):
      overflow_bin = int(float(overflow)/float(xmax) * hist.GetNbinsX())
      set_overflow = float(hist.Integral(overflow_bin,hist.GetNbinsX())) 
      hist.SetBinContent(overflow_bin + 1,set_overflow)
      hist.SetBinError(overflow_bin + 1,sqrt(set_overflow))
      hist.SetAxisRange(xmin,overflow,"X")
      additional_bin = (xmax/hist.GetNbinsX())
      self.max_min = [xmin,overflow+additional_bin]

  def Reversed_Integrator(self,hist):
      self.histclone = hist.Clone()
      integral_keeper = []
      self.reversed = 1
      for i in range(self.histclone.GetNbinsX()): integral_keeper.append(self.histclone.Integral(self.histclone.GetNbinsX()-i,self.histclone.GetNbinsX()))
      for k in range(self.histclone.GetNbinsX()): 
        self.histclone.SetBinContent(self.histclone.GetNbinsX()-k,integral_keeper[k])  
        self.histclone.SetBinError(self.histclone.GetNbinsX()-k,sqrt(integral_keeper[k]))
      self.histclone_keeper.append(self.histclone) 

class Webpage_Maker(object):

      def __init__(self,plotnames,foldername,category,option=""):
        if category == "Had":
          self.category = ""
          self.title = category
        else: 
          self.category = category
          self.title = category
        self.binning = foldername
        self.Make_Page(plotnames,self.binning,option)

      def ensure_dir(self,dir):
        try: os.makedirs(dir)
        except OSError as exc: pass
      
      def Make_Page(self,plotnames,foldername,option=""): 
        print "\n       ================================" 
        print "       ======== Making Webpage ========"
        print "       ********************************\n\n"
        self.webdir = self.title+"_plots_"+strftime("%d_%b_%H")
        self.ensure_dir("/home/hep/db1110/public_html/Website_Plots/"+self.webdir)
        for root,dirs,files in os.walk('./Plots'):
          for filename in fnmatch.filter(files,'*'):
              name = os.path.join(root,filename)
              os.system('cp ' +name+ ' ~/public_html/Website_Plots/'+self.webdir+'/')
        
        if option == "Normal":
          for i in plotnames:
              counter = 0
              htF = open('home/hep/db1110/public_html/Website_Plots/'+self.webdir+'/'+i+'.html','w')
              htF.write('Author: Darren Burton <br> \n')
              htF.write('<script language="Javascript"> \n document.write("Last Modified: " + document.lastModified + ""); \n </script> <br> \n ')
              htF.write('<center>\n <p> \n <font size="5"> Binned Muon Control Sample </font>\n </p>\n') 
              htF.write('<font size="3">Results for :  '+i+' </font><br> \n')
              htF.write('Hist Name: ') 
              for k in plotnames:
                counter += 1
                htF.write('<a href=\"'+k+'.html\">'+k+'</a>')
                htF.write('   |    ')
                if counter == 5:
                  htF.write('<br> \n')
                  counter = 0
              htF.write('<br> \n')

              for root,dirs,files in os.walk('./'+self.webdir):
                sorter = []
                for filenames in fnmatch.filter(files,i+'_*.png'): 
                  sorter.append(filenames)
                  sorter.sort()
                for plot in sorter: htF.write('<a href='+plot+'><img height=\"400\" src=\"'+plot+'\"></a> \n') 

        if option == "btag":
          self.btag_slices = {'Zero':"0-btag",'One':"1-btag",'Two':"2-btag","Inclusive":"Inclusive",'More_Than_Zero':"A btag",'More_Than_Two':"More Than Two"}
          self.btag_names = {'More_Than_Two':"_btag_morethantwo_"+self.category,'More_Than_Zero':"_btag_morethanzero_"+self.category,'Zero':"_btag_zero_"+self.category,'One':"_btag_one_"+self.category,'Two':"_btag_two_"+self.category,"Inclusive":'_'+self.category }

          self.Alpha_Webpage(foldername,plotnames,link="Zero",outertitle="HT Bins:  ")
          self.Alpha_Webpage(self.btag_slices,plotnames,link=foldername[0],outertitle="Btag Multiplicities:  ",slice="True")
          
          #Simplified plots
          self.Alpha_Webpage(foldername,plotnames,link="Zero",outertitle="HT Bins:  ",simplified = "True")
          self.Alpha_Webpage(self.btag_slices,plotnames,link=foldername[0],outertitle="Btag Multiplicities:  ",slice="True",simplified = "True")

          #Stacked plots
          self.Alpha_Webpage(foldername,plotnames,link="Zero",outertitle="HT Bins:  ",stacked = "True")
          self.Alpha_Webpage(self.btag_slices,plotnames,link=foldername[0],outertitle="Btag Multiplicities:  ",slice="True",stacked = "True")

      def Alpha_Webpage(self,outer,inner,link="",outertitle="",slice="",simplified="",stacked=""):
 
        for i in outer:
          for j in inner:
            counter = 0
            if simplified:htF = open('/home/hep/db1110/public_html/Website_Plots/'+self.webdir+'/Simplified_'+j+'_'+i+'.html','w')
            elif stacked: htF = open('/home/hep/db1110/public_html/Website_Plots/'+self.webdir+'/Stacked_'+j+'_'+i+'.html','w')
            else: htF = open('/home/hep/db1110/public_html/Website_Plots/'+self.webdir+'/'+j+'_'+i+'.html','w')
            htF.write('Author: Darren Burton <br> \n')
            htF.write('<script language="Javascript"> \n document.write("Last Modified: " + document.lastModified + ""); \n </script> <br> \n ')
            htF.write('<center>\n <p> \n <font size="5"> '+self.title+' Plots </font>\n </p>\n') 
            htF.write('<font size="3">Results for '+j+'_'+i+' </font><br> \n')
            htF.write('Hist Name: ')
            for k in inner:
              counter += 1
              if simplified: htF.write('<a href=\"Simplified_'+k+'_'+i+'.html\">'+k+'</a>   ')
              elif stacked: htF.write('<a href=\"Stacked_'+k+'_'+i+'.html\">'+k+'</a>   ')
              else: htF.write('<a href=\"'+k+'_'+i+'.html\">'+k+'</a>   ')
              htF.write('    |     ')
              if counter == 4:
                htF.write('<br> \n')
                counter = 0
            htF.write('<br> \n')
            htF.write('<br> \n')
            htF.write(outertitle)
            for k in outer: 
              if simplified:htF.write('<a href=\"Simplified_'+j+'_'+k+'.html\">'+(self.btag_slices[k] if slice else k)+'</a>     /    ')
              elif stacked:htF.write('<a href=\"Stacked_'+j+'_'+k+'.html\">'+(self.btag_slices[k] if slice else k)+'</a>     /    ')
              else:htF.write('<a href=\"'+j+'_'+k+'.html\">'+(self.btag_slices[k] if slice else k)+'</a>     /    ')
            htF.write('<br> \n')
            htF.write('<br> \n')
            if simplified:htF.write('Change Evolution Type: <a href=\"Simplified_'+j+'_'+link+'.html\">'+ ('HT Evolution' if slice else 'Btag Evolution')+'</a>')
            elif stacked:htF.write('Change Evolution Type: <a href=\"Stacked_'+j+'_'+link+'.html\">'+ ('HT Evolution' if slice else 'Btag Evolution')+'</a>')
            else:htF.write('Change Evolution Type: <a href=\"'+j+'_'+link+'.html\">'+ ('HT Evolution' if slice else 'Btag Evolution')+'</a>')
            htF.write('<br> \n')
            htF.write('<br> \n')
            htF.write(' Toggle Full/Basic/Stacked Plots:')
            if simplified:
              htF.write('<a href=\"'+j+'_'+i+'.html\">' + '  Full </a>')
              htF.write('     |     ')
              htF.write('<a href=\"Stacked_'+j+'_'+i+'.html\">' + '   Stacked </a>')
            elif stacked:
              htF.write('<a href=\"'+j+'_'+i+'.html\">' + '   Full </a>')
              htF.write('     |     ')
              htF.write('<a href=\"Simplified_'+j+'_'+i+'.html\">' + '   Simplified </a>')
            else:
              htF.write('<a href=\"Stacked_'+j+'_'+i+'.html\">' + '   Stacked </a>')
              htF.write('     |     ')
              htF.write('<a href=\"Simplified_'+j+'_'+i+'.html\">' + '   Simplified </a>')
            
            htF.write('<br> \n')
            htF.write('<br> \n')
            htF.write(' Return to Home Page:')
            htF.write('<a href=\"../RA1_Website_Plots.html\"> Go </a>')
            htF.write('<br><br>')
            
            jet_array = ['2','3','all']
            btag_array = ['_','_btag_morethanzero_','_btag_morethanone_','_btag_morethantwo_','_btag_zero_','_btag_one_','_btag_two_']
            if self.title == "Had": 
              for num,entry in enumerate(btag_array): btag_array[num] = (entry.rstrip('_'))
              for num,entry in enumerate(self.btag_names): 
                self.btag_names[entry] = self.btag_names[entry].rstrip('_')
                if num == 0: self.btag_names[entry] = ""
            for root,dirs,files in os.walk('/home/hep/db1110/public_html/Website_Plots/'+self.webdir):
              sorter = []
              test_sorter = []
              if not slice:
                for multi in btag_array:
                  for label in jet_array:
                    if not stacked:
                      for filenames in fnmatch.filter(files,('Simplified_'+j.strip('all')+label+multi+self.category+'_'+i+'*.png' if simplified == "True" else j.strip('all')+label+multi+self.category+'_'+i+'*.png')):
                        sorter.append(filenames)
                    else:
                      for filenames in fnmatch.filter(files,'Stacked_'+j.strip('all')+label+multi+self.category+'_'+i+'*.png'):
                        sorter.append(filenames)
              else:
                for bin in self.binning:
                  for label in jet_array:
                    if not stacked:
                      for filenames in fnmatch.filter(files,('Simplified_'+j.strip('all')+label+self.btag_names[i]+'_'+bin+'.png' if simplified == "True" else j.strip('all')+label+self.btag_names[i]+'_'+bin+'.png')):
                        sorter.append(filenames)
                    else: 
                      for filenames in fnmatch.filter(files,'Stacked_'+j.strip('all')+label+self.btag_names[i]+'_'+bin+'.png'):
                        sorter.append(filenames)
              for plot in sorter: htF.write('<a href='+plot+'><img width=\"30%\" src=\"'+plot+'\"></a> \n') 


