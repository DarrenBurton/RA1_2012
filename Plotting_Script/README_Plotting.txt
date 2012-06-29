To Use Plotting Script. Run Muon/Had/DiMuon/Photon_Plot_Producer.py

options are specified in settings:

  dirs : Directories in the root files ( Dont change)
  Plots : plots to be drawn, give the full string of the plot you wish to be plotted
  Lumo : Luminosity of the sample
  Webpage : btag (Dont change) just sets style of webpage produced
  Category : Had/OneMuon/DiMuon/Photon necessary for plot titles and sort of webpages so no need to change
  WebBinning : Which HT bins to display on the website
  Misc: misc options, only current option is "Normalise" this normalises the data and MC to unit area 1. used just for shape comparison plots
  Trigger: Trigger efficiencies per HT bin


At the bottom Plotter is called which is withing Btag_8TeV_Plots.py  

  options passed:
    . settings
    . dictionary of files to be run on
    . jet_multiplicity : if True will split plots into 2, >=3 and inclusive. i.e 3x the usual amount of plots will be made.


Plots are output to a folder in the directory the code is run called Plots and then transfered to the webpage output
Webpage output is directed to my public space at IC, change this path in Btag_8TeV_Plots line 805,814,809



#============================================

Changes to the way plots are drawn can be specified within the Hist_Options function in Btag_8TeV_Plots.py ~ line 485

eg. 

if "EffectiveMass_after_alphaT_55" in histogram:         <======= String name of plot
          if canvas:self.Log_Setter(plot,canvas,0.5)     <============ set log scale, delete if you dont want log scale, last argument specifies starting point of y-scale
          if word: 
            plot.GetYaxis().SetTitleOffset(1.3)           
            plot.GetYaxis().SetTitle("Events / 50 GeV") 
          if not norebin:
            plot.Rebin(50)                               <=========== Rebinning amount, always include within if norebin or else plots are continually rebinned.
            self.OverFlow_Bin(plot,0,3000,1500)          <============ Setoverflow bin, (plot,xmin,xmas,overflow value)
            #self.Reversed_Integrator(plot) 



#==============================================

Homepage site can be run by running ./Website_Maker.py within the directory containing the other plot folders. Its produces a link within the webpage based upon the first string before '_'  of the output plot folder i.e HadPlots_June18 would have a hyper link of HadPlots. 
