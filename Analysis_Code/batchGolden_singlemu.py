#!/usr/bin/env python
"""
Created by Bryn Mathias on 2010-05-07.
"""
# -----------------------------------------------------------------------------
# Necessary includes
import errno
import os
import setupSUSY
from libFrameworkSUSY import *
from libHadronic import *
from libWPol import *
from libOneLepton import *
#from libbryn import *
from icf.core import PSet,Analysis
from time import strftime
from icf.config import defaultConfig
from icf.utils import json_to_pset
from copy import deepcopy
from utils import *
# -----------------------------------------------------------------------------
# Samples
#import yours in your running script
def ensure_dir(path):
    try:
      os.makedirs(path)
    except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST:
        pass
      else: raise

# -----------------------------------------------------------------------------
# Reading the collections from the ntuple

default_ntuple = deepcopy(defaultConfig.Ntuple)
default_ntuple.Electrons=PSet(
Prefix="electron",
Suffix="Pat",
LooseID="EIDLoose",
TightID="EIDTight",
)
default_ntuple.Muons=PSet(
Prefix="muon",
Suffix="Pat",
LooseID="IsGlobalMuon",
TightID="IDGlobalMuonPromptTight",
)
default_ntuple.SecMuons=PSet(
    Prefix="muon",
    Suffix="PF")
default_ntuple.Taus=PSet(
Prefix="tau",
Suffix="Pat",
LooseID="TauIdbyTaNCfrOnePercent",
TightID="TauIdbyTaNCfrTenthPercent"
)
default_ntuple.Jets=PSet(
Prefix="ic5Jet",
Suffix="Pat",
Uncorrected=False,
)
default_ntuple.Photons=PSet(
Prefix="photon",
Suffix="Pat",
)

ic5_calo = deepcopy(default_ntuple)
ic5_calo.Jets.Prefix="ic5Jet"

ak5_calo = deepcopy(default_ntuple)
ak5_calo.Jets.Prefix="ak5Jet"

ak5_jpt = deepcopy(default_ntuple)
ak5_jpt.Jets.Prefix="ak5JetJPT"

ak5_pf = deepcopy(default_ntuple)
ak5_pf.Jets.Prefix="ak5JetPF"
ak5_pf.TerJets.Prefix="ak5Jet"

ak7_calo = deepcopy(default_ntuple)
ak7_calo.Jets.Prefix="ak7Jet"



# -----------------------------------------------------------------------------
# Cross-cleaning settings

default_cc = deepcopy(defaultConfig.XCleaning)
default_cc.Verbose=False
default_cc.MuonJet=True
default_cc.ElectronJet=True
default_cc.PhotonJet=True
default_cc.ResolveConflicts=True
default_cc.Jets.PtCut=10.0
default_cc.Jets.EtaCut=10.0
default_cc.Muons.ModifyJetEnergy=True
default_cc.Muons.PtCut=10.0
default_cc.Muons.EtaCut=2.5
default_cc.Muons.TrkIsoCut=10.0
default_cc.Muons.CombIsoCut=10.0
default_cc.Muons.MuonJetDeltaR=0.5
default_cc.Muons.MuonIsoTypePtCutoff=0.
default_cc.Muons.RequireLooseIdForInitialFilter=False
default_cc.Electrons.PtCut=10.0
default_cc.Electrons.EtaCut=2.5
default_cc.Electrons.TrkIsoCut=-1.0
#default_cc.Electrons.CombIsoCut=10.0
default_cc.Electrons.CombIsoCut=0.15
default_cc.Electrons.ElectronJetDeltaR=0.5
default_cc.Electrons.ElectronIsoTypePtCutoff=0.
default_cc.Electrons.ElectronLooseIdRequired=True
default_cc.Electrons.ElectronTightIdRequired=False
default_cc.Electrons.RequireLooseIdForInitialFilter=False
default_cc.Photons.EtCut=25.0
#default_cc.Photons.EtaCut=2.5
#default_cc.Photons.TrkIsoCut=2.0
#default_cc.Photons.CaloIsoCut=0.2
default_cc.Photons.EtaCut=5.0
default_cc.Photons.TrkIsoCut=50.0
default_cc.Photons.CaloIsoCut=20.0
default_cc.Photons.IDReq=3
default_cc.Photons.UseID=True
default_cc.Photons.PhotonJetDeltaR=0.5
default_cc.Photons.PhotonIsoTypePtCutoff=30.
# -----------------------------------------------------------------------------
# Definition of common objects
default_common = deepcopy(defaultConfig.Common)
default_common.ApplyXCleaning=True
default_common.Jets.PtCut=50.0
default_common.Jets.EtaCut=3.0
default_common.Jets.ApplyID=True
default_common.Jets.TightID=False
default_common.Electrons.PtCut=10.0
default_common.Electrons.EtaCut=2.5
#default_common.Electrons.TrkIsoCut=-1.
#default_common.Electrons.CombIsoCut=0.15
default_common.Electrons.TrkIsoCut=-1.
default_common.Electrons.CombIsoCut=10.0
default_common.Electrons.ApplyID = True
default_common.Electrons.TightID = False
default_common.Electrons.RequireLooseForOdd = True
default_common.Muons.PtCut=10.0
default_common.Muons.EtaCut=2.5
default_common.Muons.TrkIsoCut=10.0
default_common.Muons.CombIsoCut=10.0
default_common.Muons.ApplyID = True
#default_common.Muons.TightID = True
default_common.Muons.TightID = False
default_common.Muons.RequireLooseForOdd = True
default_common.Photons.EtCut=25.0
# default_common.Photons.EtaCut=2.5
default_common.Photons.UseID=True
# the photon cuts are NOT read anyway
# default_common.Photons.TrkIsoRel=0.
# default_common.Photons.TrkIsoCut=99999.
# default_common.Photons.EcalIsoRel=0.
# default_common.Photons.EcalIsoCut=99999.
# default_common.Photons.HcalIsoRel=0.
# default_common.Photons.HcalIsoCut=99999.
# default_common.Photons.HadOverEmCut=0.5
# default_common.Photons.SigmaIetaIetaCut=0.5
##default_common.Photons.CaloIsoCut=99999.
default_common.Photons.IDReq = 3
default_common.Photons.RequireLooseForOdd = True


skim_ps=PSet(
    SkimName = "myskim",
    DropBranches = False,
    Branches = [
        " keep * "
        ]
)
doskim = SkimOp(skim_ps.ps())

genericPSet = PSet(
DirName      = "275_325Gev",
MinObjects   = 1,
MaxObjects   = 20,
StandardPlots     = True,
GenPt = EffToPSet(readBtagWeight("./Btag_Efficiency_Test.txt")).GenPt,
GenEta = EffToPSet(readBtagWeight("./Btag_Efficiency_Test.txt")).GenEta,
Pt_Eta_Eff = EffToPSet(readBtagWeight("./Btag_Efficiency_Test.txt")).Pt_Eta_Eff,

Mistag_GenPt = EffToPSet(readBtagWeight("./Btag_Mistag_Test.txt")).GenPt,
Mistag_GenEta = EffToPSet(readBtagWeight("./Btag_Mistag_Test.txt")).GenEta,
Mistag_Pt_Eta_Eff = EffToPSet(readBtagWeight("./Btag_Mistag_Test.txt")).Pt_Eta_Eff,

JetPt_Low = EffToPSet(readBtagWeight("./Btag_Systematic_Variation.txt")).GenPt,
JetPt_High = EffToPSet(readBtagWeight("./Btag_Systematic_Variation.txt")).GenEta,
Variation  = EffToPSet(readBtagWeight("./Btag_Systematic_Variation.txt")).Pt_Eta_Eff,
DY = False,
)

genericPSetDY = PSet(
DirName      = "275_325Gev",
MinObjects   = 1,
MaxObjects   = 20,
StandardPlots     = True,
GenPt = EffToPSet(readBtagWeight("./Btag_Efficiency_Test.txt")).GenPt,
GenEta = EffToPSet(readBtagWeight("./Btag_Efficiency_Test.txt")).GenEta,
Pt_Eta_Eff = EffToPSet(readBtagWeight("./Btag_Efficiency_Test.txt")).Pt_Eta_Eff,

Mistag_GenPt = EffToPSet(readBtagWeight("./Btag_Mistag_Test.txt")).GenPt,
Mistag_GenEta = EffToPSet(readBtagWeight("./Btag_Mistag_Test.txt")).GenEta,
Mistag_Pt_Eta_Eff = EffToPSet(readBtagWeight("./Btag_Mistag_Test.txt")).Pt_Eta_Eff,

JetPt_Low = EffToPSet(readBtagWeight("./Btag_Systematic_Variation.txt")).GenPt,
JetPt_High = EffToPSet(readBtagWeight("./Btag_Systematic_Variation.txt")).GenEta,
Variation  = EffToPSet(readBtagWeight("./Btag_Systematic_Variation.txt")).Pt_Eta_Eff,
DY = True,
)



def makePlotOp(OP = (), cutTree = None, cut = None, label = ""):
  """docstring for makePlotOp"""
  out = []
  if OP[1] != None:
    plotpset = deepcopy(OP[1])
    plotpset.DirName = label
    print plotpset.DirName
    op = eval(OP[0]+"(plotpset.ps())")

  else:
    op = eval(OP[0])
  out.append(op)
  cutTree.TAttach(cut,op)
  alpha = OP_HadAlphaTCut(0.55)
  skim_ps=PSet(
    SkimName = "myskim",
    DropBranches = False,
    Branches = [
        " keep * "
        ]
  )
  skim = SkimOp(skim_ps.ps())
  # out.append(skim_ps)
  cutTree.TAttach(cut,alpha)
  #cutTree.TAttach(alpha,skim)
  #out.append(skim)
  out.append(alpha)
  #out.append(dump)
  return out
  pass

def AddBinedHist(cutTree = None, OP = (), cut = None, htBins = [],TriggerDict = None,lab = ""):
  """docstring for AddBinedHist"""
  out = []
  if TriggerDict is not None:
      for lower,upper in zip(htBins,htBins[1:]+[None]):
        # print "Lower , Upper =", lower , upper
        if int(lower) == 325 and upper is None: continue
        if int(lower) == 375 and upper is None: continue
        if int(lower) == 475 and upper is None: continue
        if int(lower) == 675 and upper is None: continue
        # print "continue should have happened now"
        lowerCut = eval("RECO_CommonHTCut(%d)"%lower)
        triggerps = PSet(Verbose = False,
        UsePreScaledTriggers = False,
        Triggers = [])
        triggerps.Triggers = TriggerDict["%d%s"%(lower, "_%d"%upper if upper else "")]
        Trigger = OP_MultiTrigger( triggerps.ps() )
        out.append(triggerps)
        out.append(Trigger)
        out.append(lowerCut)
        cutTree.TAttach(cut,Trigger)
        cutTree.TAttach(Trigger,lowerCut)
        if upper:
          upperCut =  eval("RECO_CommonHTLessThanCut(%d)"%upper)
          out.append(upperCut)
          cutTree.TAttach(lowerCut,upperCut)
        pOps = makePlotOp(cutTree = cutTree, OP = OP, cut = upperCut if upper else lowerCut, label = "%s%d%s"%(lab,lower, "_%d"%upper if upper else ""))
        out.append(pOps) 
  else:
      for lower,upper in zip(htBins,htBins[1:]+[None]):
        #print "Lower , Upper =", lower , upper
        if int(lower) == 325 and upper is None: continue
        if int(lower) == 375 and upper is None: continue
        if int(lower) == 475 and upper is None: continue
        if int(lower) == 675 and upper is None: continue
        # print "continue should have happened now"
        lowerCut = eval("RECO_CommonHTCut(%d)"%lower)
        out.append(lowerCut)
        cutTree.TAttach(cut,lowerCut)
        if upper:
          upperCut =  eval("RECO_CommonHTLessThanCut(%d)"%upper)
          out.append(upperCut)
          cutTree.TAttach(lowerCut,upperCut)
        pOps = makePlotOp(cutTree = cutTree, OP = OP, cut = upperCut if upper else lowerCut, label = "%s%d%s"%(lab,lower, "_%d"%upper if upper else ""))
        out.append(pOps)
  return out
  pass

# Common cut definitions
#Avaiable criteria for MC and for Data are at current slightly different Hence the making of two trees
#DataOnly!

NoiseFilt= OP_HadronicHBHEnoiseFilter()
GoodVertexMonster = OP_GoodEventSelection()

#Standard Event Selection
LeadingJetEta = OP_FirstJetEta(2.5)
unCorLeadJetCut = OP_UnCorLeadJetCut(30.)
LeadingJetPtCut = OP_FirstJetPtCut(100.)
oddMuon = OP_OddMuon()
oddElectron = OP_OddElectron()
oddPhoton = OP_OddPhoton()
oddJet = OP_OddJet()
badMuonInJet = OP_BadMuonInJet()
numComElectrons = OP_NumComElectrons("<=",0)
numComPhotons = OP_NumComPhotons("<=",0)

LeadingPhotonEta = OP_FirstPhotonEta(1.4442)
LessThan375 = RECO_CommonHTLessThanCut(375.)
ht250_Trigger = RECO_CommonHTCut(250.)
htCut275_2 = RECO_CommonHTCut(275.)
HT875Cut = RECO_CommonHTCut(875.)

ht275_Fail     = RECO_CommonHTCut(275.)
ht325_Fail     = RECO_CommonHTCut(325.)
ht375_Fail     = RECO_CommonHTCut(375.)
htLess325 = RECO_CommonHTLessThanCut(325.)
htLess375_Fail = RECO_CommonHTLessThanCut(375.)

htCut275 = RECO_CommonHTCut(275.)
htCut375GeV = RECO_CommonHTCut(375.)
htCut375All = RECO_CommonHTCut(375.)
htCut275_2 = RECO_CommonHTCut(275.)
htCut375GeV = RECO_CommonHTCut(375.)
htCut375 = RECO_CommonHTCut(375.)

HT575 = RECO_CommonHTCut(575.)
HT375 = RECO_CommonHTCut(375.)


alphaT0 = OP_CommonAlphaTCut(0.55)
alphaT1 = OP_CommonAlphaTCut(0.55)
alphaT2 = OP_CommonAlphaTCut(0.55)
spikecleaner = OP_EcalSpikeCleaner()
event_display = OP_EventDisplay("EventDisplays", "common") #to draw all/common objects
alphat = OP_CommonAlphaTCut(0.55)
DeadEcalCutData = OP_DeadECALCut(0.3,0.3,0.5,30.,10,0,"./deadRegionList_GR10_P_V10.txt")
DeadEcalCutMC =   OP_DeadECALCut(0.3,0.3,0.5,30.,10,0,"./deadRegionList_START38_V12.txt")
MHTCut = OP_CommonMHTCut(0.)
MHT_METCut = OP_MHToverMET(1.25,50.)
MHT_HT = OP_MHTovHT(0.4)
DiJet5 = OP_NumComJets("==",2)
GE2Jets = OP_NumComJets(">=",2)
VertexPtOverHT = OP_SumVertexPtOverHT(0.1)
HadronicAlphaT = OP_HadAlphaTCut(0.55)


triggers = {
"275_325":["HLT_HT250_AlphaT0p55_v*"],
"325_375":["HLT_HT300_AlphaT0p53_v*"],
"375_475":["HLT_HT350_AlphaT0p52_v*"],
"475_575":["HLT_HT400_AlphaT0p51_v*"],
"575_675":["HLT_HT400_AlphaT0p51_v*"],
"675_775":["HLT_HT400_AlphaT0p51_v*"],
"775_875":["HLT_HT400_AlphaT0p51_v*"],
"875":["HLT_HT400_AlphaT0p51_v*"],
}

single_mu_triggers = {
  "275":["HLT_IsoMu24_eta2p1_v*"],
  "275_325":["HLT_IsoMu24_eta2p1_v*"],
  "325_375":["HLT_IsoMu24_eta2p1_v*"],
  "375_475":["HLT_IsoMu24_eta2p1_v*"],
  "475_575":["HLT_IsoMu24_eta2p1_v*"],
  "575_675":["HLT_IsoMu24_eta2p1_v*"],
  "675_775":["HLT_IsoMu24_eta2p1_v*"],
  "775_875":["HLT_IsoMu24_eta2p1_v*"], 
  "875":["HLT_IsoMu24_eta2p1_v*"],

}

photon_trigger = ["HLT_Photon150_v%d"%i for i in range(1,4)] 

zerobtagDiMuon= OP_NumCommonBtagJets("==",0,0.679,5)
more_than_zero_btagDiMuon= OP_NumCommonBtagJets(">",0,0.679,5)
onebtagDiMuon= OP_NumCommonBtagJets("==",1,0.679,5)
more_than_one_btagDiMuon= OP_NumCommonBtagJets(">",1,0.679,5)
twobtagDiMuon= OP_NumCommonBtagJets("==",2,0.679,5)
more_than_two_btagDiMuon= OP_NumCommonBtagJets(">",2,0.679,5)

zerobtagOneMuon= OP_NumCommonBtagJets("==",0,0.679,5)
more_than_zero_btagOneMuon= OP_NumCommonBtagJets(">",0,0.679,5)
onebtagOneMuon= OP_NumCommonBtagJets("==",1,0.679,5)
more_than_one_btagOneMuon= OP_NumCommonBtagJets(">",1,0.679,5)
twobtagOneMuon= OP_NumCommonBtagJets("==",2,0.679,5)
more_than_two_btagOneMuon= OP_NumCommonBtagJets(">",2,0.679,5)

zerobtag= OP_NumCommonBtagJets("==",0,0.679,5)
more_than_zero_btag= OP_NumCommonBtagJets(">",0,0.679,5)
onebtag= OP_NumCommonBtagJets("==",1,0.679,5)
more_than_one_btag= OP_NumCommonBtagJets(">",1,0.679,5)
twobtag= OP_NumCommonBtagJets("==",2,0.679,5)
more_than_two_btag= OP_NumCommonBtagJets(">",2,0.679,5)
more_than_three_btag= OP_NumCommonBtagJets(">",3,0.679,5)

recHitCut = OP_SumRecHitPtCut(30.)
ZeroMuon = OP_NumComMuons("<=",0)
json_ouput = JSONOutput("filtered")
OneMuon = OP_NumComMuons("==",1)
ZMassCut = RECO_2ndMuonMass(25.0, 91.2, False, "all")
PFMTCut30 = RECO_PFMTCut(30.0)
DiMuon = OP_NumComMuons("==",2)
ZMass_2Muons = RECO_DiMuon_ZMass(25.0, 91.2, True, "OS")
minDRMuonJetCut = RECO_MuonJetDRCut(0.5)
minDRMuonJetCutDiMuon = RECO_MuonJetDRCut(0.5)
Tot_VertexCut = OP_TotVertexCut(2,19)
GenWMuon = MC_GenWMuonExists()

MET_Filter = OP_METFilters_2012()
LeadingPhotonPt_100 = OP_PhotonPtCut(100.)
LeadingPhotonPt_165 = OP_PhotonPtCut(165.)

PhotonDR_Cut = OP_PhotonJetDRCut(1.0)
OnePhoton = OP_NumComPhotons("==",1)


def MakePhotonData(Threshold):
   out = []
   print int(Threshold)
   print "Using json %s" %str(switches()["json"])
   print photon_trigger
   secondJetET = OP_SecondJetEtCut(Threshold)
   HTBins = []
   if int(Threshold) is 100: HTBins = [275.,325.] +[375+100*i for i in range(6)]
   triggerps = PSet(Verbose = False, UsePreScaledTriggers = False,Triggers = photon_trigger)
   Trigger = OP_MultiTrigger( triggerps.ps() )

   cutTreeData = Tree("Data")
   cutTreeData.Attach(json)
   cutTreeData.TAttach(json,json_ouput)
   cutTreeData.TAttach(json_ouput,GoodVertexMonster)
   cutTreeData.TAttach(GoodVertexMonster,NoiseFilt)
   cutTreeData.TAttach(NoiseFilt,Trigger)
   cutTreeData.TAttach(Trigger,OnePhoton)
   cutTreeData.TAttach(OnePhoton,LeadingPhotonPt_165)
   cutTreeData.TAttach(LeadingPhotonPt_165,LeadingPhotonEta)
   cutTreeData.TAttach(LeadingPhotonEta,LeadingJetEta)

   #out.append(AddBinedHist(cutTree = cutTreeData,
   #OP = ("BtagSystematicPlots",genericPSet), cut = LeadingJetEta,
   #htBins = HTBins,TriggerDict = None,lab ="btag_zero_Photon_") )
   
   cutTreeData.TAttach(LeadingJetEta,secondJetET)
   cutTreeData.TAttach(secondJetET,HT375)
   cutTreeData.TAttach(HT375,ZeroMuon)
   cutTreeData.TAttach(ZeroMuon,numComElectrons)
   cutTreeData.TAttach(numComElectrons,oddJet)
   cutTreeData.TAttach(oddJet,oddElectron)
   cutTreeData.TAttach(oddElectron,oddPhoton)
   cutTreeData.TAttach(oddPhoton,PhotonDR_Cut)

   #out.append(AddBinedHist(cutTree = cutTreeData,
   #OP = ("BtagSystematicPlots",genericPSet), cut = PhotonDR_Cut,
   #htBins = HTBins,TriggerDict = None,lab ="btag_one_Photon_") )

   cutTreeData.TAttach(PhotonDR_Cut,HadronicAlphaT)
   cutTreeData.TAttach(HadronicAlphaT,MHT_METCut)
   cutTreeData.TAttach(MHT_METCut,DeadEcalCutData)
   cutTreeData.TAttach(DeadEcalCutData,zerobtag)
   cutTreeData.TAttach(DeadEcalCutData,more_than_zero_btag)
   cutTreeData.TAttach(DeadEcalCutData,onebtag)
   cutTreeData.TAttach(DeadEcalCutData,more_than_one_btag)
   cutTreeData.TAttach(DeadEcalCutData,twobtag)
   cutTreeData.TAttach(DeadEcalCutData,more_than_two_btag)
   
   out.append(AddBinedHist(cutTree = cutTreeData,
   OP = ("BtagSystematicPlots",genericPSet), cut = zerobtag,
   htBins = HTBins,TriggerDict = None,lab ="btag_zero_Photon_") )

   out.append(AddBinedHist(cutTree = cutTreeData,
   OP = ("BtagSystematicPlots",genericPSet), cut = more_than_zero_btag,
   htBins = HTBins,TriggerDict = None,lab ="btag_morethanzero_Photon_") )

   out.append(AddBinedHist(cutTree = cutTreeData,
   OP = ("BtagSystematicPlots",genericPSet), cut = onebtag,
   htBins = HTBins,TriggerDict = None,lab ="btag_one_Photon_") )

   out.append(AddBinedHist(cutTree = cutTreeData,
   OP = ("BtagSystematicPlots",genericPSet), cut = more_than_one_btag,
   htBins = HTBins,TriggerDict = None,lab ="btag_morethanone_Photon_") )

   out.append(AddBinedHist(cutTree = cutTreeData,
   OP = ("BtagSystematicPlots",genericPSet), cut = twobtag,
   htBins = HTBins,TriggerDict = None,lab ="btag_two_Photon_") )

   out.append(AddBinedHist(cutTree = cutTreeData,
   OP = ("BtagSystematicPlots",genericPSet), cut = more_than_two_btag,
   htBins = HTBins,TriggerDict = None,lab ="btag_morethantwo_Photon_") )
   
   out.append(AddBinedHist(cutTree = cutTreeData,
   OP = ("BtagSystematicPlots",genericPSet), cut = DeadEcalCutData,
   htBins = HTBins,TriggerDict = None,lab ="Photon_") )
   
   return (cutTreeData,secondJetET,Trigger,out)

def MakeDataTree(Threshold,Muon = None,Split = None):
  out = []
  print int(Threshold)
  print "Using json %s" %str(switches()["json"])

  secondJetET = OP_SecondJetEtCut(Threshold)
  HTBins = []
  if int(Threshold) is 100 and Split == "Reweight" : HTBins = [375.] 
  if int(Threshold) is 100 and Split == None : HTBins = [375+100*i for i in range(6)]
  if int(Threshold) is 100 and Split == "Had_One" : HTBins = [375+100*i for i in range(4)]
  if int(Threshold) is 100 and Split == "Had_Two" : HTBins = [675+100*i for i in range(3)]
  if int(Threshold) is 100 and Split == "Muon_One" : HTBins = [375+100*i for i in range(4)]
  if int(Threshold) is 100 and Split == "Muon_Two" : HTBins = [675+100*i for i in range(3)]
  if int(Threshold) is 73 :  HTBins = [275.,325.]
  if int(Threshold) is 86 :  HTBins = [325.,375.]

  if int(Threshold) in [73,86]:  Leading_MuPtCut = OP_LowerMuPtCut(30.)
  else: Leading_MuPtCut = OP_LowerMuPtCut(30.)
  # from batchGolden import *
  cutTreeData = Tree("Data")
  
  cutTreeData.Attach(json)
  # cutTreeData.TAttach(json,NoiseFilt)
  cutTreeData.TAttach(json,json_ouput) 
  cutTreeData.TAttach(json_ouput,GE2Jets)
  cutTreeData.TAttach(GE2Jets,NoiseFilt)
  cutTreeData.TAttach(NoiseFilt,MET_Filter)
  cutTreeData.TAttach(MET_Filter,GoodVertexMonster)
  cutTreeData.TAttach(GoodVertexMonster,recHitCut)
  cutTreeData.TAttach(recHitCut,LeadingJetEta)
  cutTreeData.TAttach(LeadingJetEta,secondJetET)
  cutTreeData.TAttach(secondJetET,oddJet)
  cutTreeData.TAttach(oddJet,badMuonInJet)
  cutTreeData.TAttach(badMuonInJet,oddElectron)
  cutTreeData.TAttach(oddElectron,oddPhoton)
  cutTreeData.TAttach(oddPhoton,numComElectrons)
  cutTreeData.TAttach(numComElectrons,numComPhotons)
  cutTreeData.TAttach(numComPhotons,VertexPtOverHT)
  cutTreeData.TAttach(VertexPtOverHT,htCut275)
  #Begin MHT/MET plot inthe low region.
  if Muon == None:
	
      cutTreeData.TAttach(htCut275,DeadEcalCutData)
      cutTreeData.TAttach(DeadEcalCutData,MHT_METCut)
      cutTreeData.TAttach(MHT_METCut,ZeroMuon)
      cutTreeData.TAttach(ZeroMuon,zerobtag)
      cutTreeData.TAttach(ZeroMuon,more_than_zero_btag)
      cutTreeData.TAttach(ZeroMuon,onebtag)
      cutTreeData.TAttach(ZeroMuon,more_than_one_btag)
      cutTreeData.TAttach(ZeroMuon,twobtag)
      cutTreeData.TAttach(ZeroMuon,more_than_two_btag)
     
      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("BtagSystematicPlots",genericPSet), cut = zerobtag,
      htBins = HTBins,TriggerDict = triggers,lab ="btag_zero_") )

      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("BtagSystematicPlots",genericPSet), cut = more_than_zero_btag,
      htBins = HTBins,TriggerDict = triggers,lab ="btag_morethanzero_") )
       
      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("BtagSystematicPlots",genericPSet), cut = onebtag,
      htBins = HTBins,TriggerDict = triggers,lab ="btag_one_") )

      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("BtagSystematicPlots",genericPSet), cut = more_than_one_btag,
      htBins = HTBins,TriggerDict = triggers,lab ="btag_morethanone_") )
      
      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("BtagSystematicPlots",genericPSet), cut = twobtag,
      htBins = HTBins,TriggerDict = triggers,lab ="btag_two_") )

      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("BtagSystematicPlots",genericPSet), cut = more_than_two_btag,
      htBins = HTBins,TriggerDict = triggers,lab ="btag_morethantwo_") )
      
      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("BtagSystematicPlots",genericPSet), cut = ZeroMuon,
      htBins = HTBins,TriggerDict = triggers,lab ="") )
      
  else:
      cutTreeData.TAttach(htCut275,DeadEcalCutData)
      cutTreeData.TAttach(DeadEcalCutData,MHT_METCut)
      cutTreeData.TAttach(MHT_METCut,Leading_MuPtCut)
      cutTreeData.TAttach(Leading_MuPtCut,minDRMuonJetCut)
      cutTreeData.TAttach(minDRMuonJetCut,OneMuon)
      cutTreeData.TAttach(OneMuon,ZMassCut)
      cutTreeData.TAttach(ZMassCut,PFMTCut30)
      cutTreeData.TAttach(PFMTCut30,zerobtagOneMuon)
      cutTreeData.TAttach(PFMTCut30,more_than_zero_btagOneMuon)
      cutTreeData.TAttach(PFMTCut30,onebtagOneMuon)
      cutTreeData.TAttach(PFMTCut30,more_than_one_btagOneMuon)
      cutTreeData.TAttach(PFMTCut30,twobtagOneMuon)
      cutTreeData.TAttach(PFMTCut30,more_than_two_btagOneMuon)

          
      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("Muon_ControlPlots",genericPSet), cut = zerobtagOneMuon,
      htBins = HTBins,TriggerDict = single_mu_triggers ,lab ="btag_zero_OneMuon_") )

      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("Muon_ControlPlots",genericPSet), cut = more_than_zero_btagOneMuon,
      htBins = HTBins,TriggerDict = single_mu_triggers ,lab ="btag_morethanzero_OneMuon_") )

      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("Muon_ControlPlots",genericPSet), cut = onebtagOneMuon,
      htBins = HTBins,TriggerDict = single_mu_triggers ,lab ="btag_one_OneMuon_") )

      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("Muon_ControlPlots",genericPSet), cut = more_than_one_btagOneMuon,
      htBins = HTBins,TriggerDict = single_mu_triggers ,lab ="btag_morethanone_OneMuon_") )

      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("Muon_ControlPlots",genericPSet), cut = twobtagOneMuon,
      htBins = HTBins,TriggerDict = single_mu_triggers ,lab ="btag_two_OneMuon_") )

      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("Muon_ControlPlots",genericPSet), cut = more_than_two_btagOneMuon,
      htBins = HTBins,TriggerDict = single_mu_triggers ,lab ="btag_morethantwo_OneMuon_") )
       
      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("Muon_ControlPlots",genericPSet), cut = PFMTCut30,
      htBins = HTBins,TriggerDict = single_mu_triggers,lab = "OneMuon_") )
       
      cutTreeData.TAttach(minDRMuonJetCut,DiMuon)
      cutTreeData.TAttach(DiMuon,ZMass_2Muons)
      cutTreeData.TAttach(ZMass_2Muons,zerobtagDiMuon) 
      cutTreeData.TAttach(ZMass_2Muons,more_than_zero_btagDiMuon)
      cutTreeData.TAttach(ZMass_2Muons,onebtagDiMuon) 
      cutTreeData.TAttach(ZMass_2Muons,more_than_one_btagDiMuon)
      cutTreeData.TAttach(ZMass_2Muons,twobtagDiMuon) 
      cutTreeData.TAttach(ZMass_2Muons,more_than_two_btagDiMuon)

      # avobe here does one big inclusive bin!
      # Now lets start binning in HT bins
      # So we can HADD the files at the end and get a chorent set to save the book keeping nightmare:
      # we arrange the HT bins so they are not repoduced though out threshold runs.
      
      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("Muon_ControlPlots",genericPSet), cut = zerobtagDiMuon,
      htBins = HTBins,TriggerDict = single_mu_triggers,lab ="btag_zero_DiMuon_") )
      
      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("Muon_ControlPlots",genericPSet), cut = more_than_zero_btagDiMuon,
      htBins = HTBins,TriggerDict = single_mu_triggers,lab ="btag_morethanzero_DiMuon_") )

      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("Muon_ControlPlots",genericPSet), cut = onebtagDiMuon,
      htBins = HTBins,TriggerDict = single_mu_triggers,lab ="btag_one_DiMuon_") )
   
      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("Muon_ControlPlots",genericPSet), cut = more_than_one_btagDiMuon,
      htBins = HTBins,TriggerDict = single_mu_triggers,lab ="btag_morethanone_DiMuon_") )
   
      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("Muon_ControlPlots",genericPSet), cut = twobtagDiMuon,
      htBins = HTBins,TriggerDict = single_mu_triggers,lab ="btag_two_DiMuon_") )
   
      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("Muon_ControlPlots",genericPSet), cut = more_than_two_btagDiMuon,
      htBins = HTBins,TriggerDict = single_mu_triggers,lab ="btag_morethantwo_DiMuon_") )
      
      out.append(AddBinedHist(cutTree = cutTreeData,
      OP = ("Muon_ControlPlots",genericPSet), cut = ZMass_2Muons,
      htBins = HTBins,TriggerDict = single_mu_triggers,lab = "DiMuon_") )
      
  return (cutTreeData,secondJetET,Leading_MuPtCut,out)

#Second MC!
def MakePhotonMC(Threshold):
  out = []
  HTBins = []
  secondJetET = OP_SecondJetEtCut(Threshold)
  if int(Threshold) is 100: HTBins = [275.,325.] + [375+100*i for i in range(6)]

  cutTreeMC = Tree("MC")
  
  cutTreeMC.Attach(GoodVertexMonster)
  cutTreeMC.TAttach(GoodVertexMonster,MET_Filter)
  cutTreeMC.TAttach(MET_Filter,NoiseFilt)
  cutTreeMC.TAttach(NoiseFilt,OnePhoton)
  cutTreeMC.TAttach(OnePhoton,LeadingPhotonPt_165)
  cutTreeMC.TAttach(LeadingPhotonPt_165,LeadingPhotonEta)
  cutTreeMC.TAttach(LeadingPhotonEta,LeadingJetEta)
  cutTreeMC.TAttach(LeadingJetEta,secondJetET)

  #out.append(AddBinedHist(cutTree = cutTreeMC,
  #OP = ("BtagSystematicPlots",genericPSet), cut = LeadingJetEta,
  #htBins = HTBins,TriggerDict = None,lab ="btag_zero_Photon_") )

  cutTreeMC.TAttach(secondJetET,HT375)
  cutTreeMC.TAttach(HT375,ZeroMuon)
  cutTreeMC.TAttach(ZeroMuon,numComElectrons)
  cutTreeMC.TAttach(numComElectrons,oddJet)
  cutTreeMC.TAttach(oddJet,oddElectron)
  cutTreeMC.TAttach(oddElectron,oddPhoton)
  cutTreeMC.TAttach(oddPhoton,PhotonDR_Cut)
  
  #out.append(AddBinedHist(cutTree = cutTreeMC,
  #OP = ("BtagSystematicPlots",genericPSet), cut = PhotonDR_Cut,
  #htBins = HTBins,TriggerDict = None,lab ="btag_one_Photon_") )
  
  cutTreeMC.TAttach(PhotonDR_Cut,HadronicAlphaT)
  cutTreeMC.TAttach(HadronicAlphaT,MHT_METCut)
  cutTreeMC.TAttach(MHT_METCut,DeadEcalCutData)
  cutTreeMC.TAttach(DeadEcalCutData,zerobtag)
  cutTreeMC.TAttach(DeadEcalCutData,more_than_zero_btag)
  cutTreeMC.TAttach(DeadEcalCutData,onebtag)
  cutTreeMC.TAttach(DeadEcalCutData,more_than_one_btag)
  cutTreeMC.TAttach(DeadEcalCutData,twobtag)
  cutTreeMC.TAttach(DeadEcalCutData,more_than_two_btag)
  
  out.append(AddBinedHist(cutTree = cutTreeMC,
  OP = ("BtagSystematicPlots",genericPSet), cut = zerobtag,
  htBins = HTBins,TriggerDict = None,lab ="btag_zero_Photon_") )

  out.append(AddBinedHist(cutTree = cutTreeMC,
  OP = ("BtagSystematicPlots",genericPSet), cut = more_than_zero_btag,
  htBins = HTBins,TriggerDict = None,lab ="btag_morethanzero_Photon_") )

  out.append(AddBinedHist(cutTree = cutTreeMC,
  OP = ("BtagSystematicPlots",genericPSet), cut = onebtag,
  htBins = HTBins,TriggerDict = None,lab ="btag_one_Photon_") )

  out.append(AddBinedHist(cutTree = cutTreeMC,
  OP = ("BtagSystematicPlots",genericPSet), cut = more_than_one_btag,
  htBins = HTBins,TriggerDict = None,lab ="btag_morethanone_Photon_") )

  out.append(AddBinedHist(cutTree = cutTreeMC,
  OP = ("BtagSystematicPlots",genericPSet), cut = twobtag,
  htBins = HTBins,TriggerDict = None,lab ="btag_two_Photon_") )

  out.append(AddBinedHist(cutTree = cutTreeMC,
  OP = ("BtagSystematicPlots",genericPSet), cut = more_than_two_btag,
  htBins = HTBins,TriggerDict = None,lab ="btag_morethantwo_Photon_") )
  
  out.append(AddBinedHist(cutTree = cutTreeMC,
  OP = ("BtagSystematicPlots",genericPSet), cut = DeadEcalCutData,
  htBins = HTBins,TriggerDict = None,lab ="Photon_") )

  return (cutTreeMC,secondJetET,out)


def MakeMCTree(Threshold, Muon = None,Split = None,DY= None):
  out = []

  HTBins = []
  
  if int(Threshold) is 100 and Split == "Reweight" : HTBins = [275.]
  if int(Threshold) is 100 and Split == "WJet_Test" : HTBins = [0.] 
  if int(Threshold) is 100 and Split == None : HTBins = [375+100*i for i in range(6)]
  if int(Threshold) is 100 and Split == "Had_One" : HTBins = [375+100*i for i in range(2)]
  if int(Threshold) is 100 and Split == "Had_Two" : HTBins = [475+100*i for i in range(5)]
  if int(Threshold) is 100 and Split == "Muon_One" : HTBins = [375+100*i for i in range(2)]
  if int(Threshold) is 100 and Split == "Muon_Two" : HTBins = [475+100*i for i in range(5)]
  if int(Threshold) is 73 : HTBins = [275.,325.]
  if int(Threshold) is 86 : HTBins = [325.,375.]

  if int(Threshold) in [73,86]:  
    Leading_MuPtCut = OP_LowerMuPtCut(30.)
    if int(Threshold) == 73: MHTCut = OP_CommonMHTCut(0.)
    if int(Threshold) == 86: MHTCut = OP_CommonMHTCut(0.)
  else: 
    Leading_MuPtCut = OP_LowerMuPtCut(30.)
    MHTCut = OP_CommonMHTCut(0.)
  secondJetET = OP_SecondJetEtCut(Threshold)
  cutTreeMC = Tree("MC")
  #cutTreeMC.Attach(numComElectrons)
  
  cutTreeMC.Attach(GE2Jets)
  cutTreeMC.TAttach(GE2Jets,ht250_Trigger)
  cutTreeMC.TAttach(ht250_Trigger,NoiseFilt)
  cutTreeMC.TAttach(NoiseFilt,GoodVertexMonster)
  cutTreeMC.TAttach(GoodVertexMonster,MET_Filter)
  cutTreeMC.TAttach(MET_Filter,recHitCut)
  cutTreeMC.TAttach(recHitCut,LeadingJetEta)
  cutTreeMC.TAttach(LeadingJetEta,secondJetET)
  cutTreeMC.TAttach(secondJetET,oddJet)
  cutTreeMC.TAttach(oddJet,badMuonInJet)
  cutTreeMC.TAttach(badMuonInJet, oddElectron)
  cutTreeMC.TAttach(oddElectron,oddPhoton)
  cutTreeMC.TAttach(oddPhoton,numComElectrons)
  cutTreeMC.TAttach(numComElectrons,numComPhotons)
  cutTreeMC.TAttach(numComPhotons,VertexPtOverHT)
  cutTreeMC.TAttach(VertexPtOverHT,htCut275)
  if Muon == None:
      cutTreeMC.TAttach(htCut275,DeadEcalCutMC) 
      cutTreeMC.TAttach(DeadEcalCutMC,MHT_METCut)
      cutTreeMC.TAttach(MHT_METCut,MHTCut)

      if DY == None:

	cutTreeMC.TAttach(MHTCut,ZeroMuon)
        cutTreeMC.TAttach(ZeroMuon,zerobtag)
        cutTreeMC.TAttach(ZeroMuon,more_than_zero_btag)
        cutTreeMC.TAttach(ZeroMuon,onebtag)
        cutTreeMC.TAttach(ZeroMuon,more_than_one_btag)
        cutTreeMC.TAttach(ZeroMuon,twobtag)
        cutTreeMC.TAttach(ZeroMuon,more_than_two_btag)
 
        out.append(AddBinedHist(cutTree = cutTreeMC,
        OP = ("BtagSystematicPlots",genericPSet), cut = ZeroMuon,
        htBins = HTBins,TriggerDict = None,lab ="") )
        
        out.append(AddBinedHist(cutTree = cutTreeMC,
        OP = ("BtagSystematicPlots",genericPSet), cut = zerobtag,
        htBins = HTBins,TriggerDict = None,lab ="btag_zero_") )
      
        out.append(AddBinedHist(cutTree = cutTreeMC,
        OP = ("BtagSystematicPlots",genericPSet), cut = more_than_zero_btag,
        htBins = HTBins,TriggerDict = None,lab ="btag_morethanzero_") )
      
        out.append(AddBinedHist(cutTree = cutTreeMC,
        OP = ("BtagSystematicPlots",genericPSet), cut = onebtag,
        htBins = HTBins,TriggerDict = None,lab ="btag_one_") )
        
        out.append(AddBinedHist(cutTree = cutTreeMC,
        OP = ("BtagSystematicPlots",genericPSet), cut = more_than_one_btag,
        htBins = HTBins,TriggerDict = None,lab ="btag_morethanone_") )
      
        out.append(AddBinedHist(cutTree = cutTreeMC,
        OP = ("BtagSystematicPlots",genericPSet), cut = twobtag,
        htBins = HTBins,TriggerDict = None,lab ="btag_two_") )
      
        out.append(AddBinedHist(cutTree = cutTreeMC,
        OP = ("BtagSystematicPlots",genericPSet), cut = more_than_two_btag,
        htBins = HTBins,TriggerDict = None,lab ="btag_morethantwo_") )
        
      else:

        cutTreeMC.TAttach(MHTCut,zerobtag)
        cutTreeMC.TAttach(zerobtag,more_than_zero_btag)
        cutTreeMC.TAttach(zerobtag,onebtag)
        cutTreeMC.TAttach(zerobtag,more_than_one_btag)
        cutTreeMC.TAttach(zerobtag,twobtag)
        cutTreeMC.TAttach(zerobtag,more_than_two_btag)

	out.append(AddBinedHist(cutTree = cutTreeMC,
        OP = ("BtagSystematicPlots",genericPSetDY), cut = MHT_METCut,
        htBins = HTBins,TriggerDict = None,lab ="") )
        
        out.append(AddBinedHist(cutTree = cutTreeMC,
        OP = ("BtagSystematicPlots",genericPSetDY), cut = zerobtag,
        htBins = HTBins,TriggerDict = None,lab ="btag_zero_") )
      
        out.append(AddBinedHist(cutTree = cutTreeMC,
        OP = ("BtagSystematicPlots",genericPSetDY), cut = more_than_zero_btag,
        htBins = HTBins,TriggerDict = None,lab ="btag_morethanzero_") )
      
        out.append(AddBinedHist(cutTree = cutTreeMC,
        OP = ("BtagSystematicPlots",genericPSetDY), cut = onebtag,
        htBins = HTBins,TriggerDict = None,lab ="btag_one_") )
        
        out.append(AddBinedHist(cutTree = cutTreeMC,
        OP = ("BtagSystematicPlots",genericPSetDY), cut = more_than_one_btag,
        htBins = HTBins,TriggerDict = None,lab ="btag_morethanone_") )
      
        out.append(AddBinedHist(cutTree = cutTreeMC,
        OP = ("BtagSystematicPlots",genericPSetDY), cut = twobtag,
        htBins = HTBins,TriggerDict = None,lab ="btag_two_") )
      
        out.append(AddBinedHist(cutTree = cutTreeMC,
        OP = ("BtagSystematicPlots",genericPSetDY), cut = more_than_two_btag,
        htBins = HTBins,TriggerDict = None,lab ="btag_morethantwo_") )
        

     
  else:
      cutTreeMC.TAttach(htCut275,DeadEcalCutMC) 
      cutTreeMC.TAttach(DeadEcalCutMC,MHT_METCut)
      cutTreeMC.TAttach(MHT_METCut,Leading_MuPtCut)
      cutTreeMC.TAttach(Leading_MuPtCut,minDRMuonJetCut)
      cutTreeMC.TAttach(minDRMuonJetCut,OneMuon)
      cutTreeMC.TAttach(OneMuon,ZMassCut)
      cutTreeMC.TAttach(ZMassCut,PFMTCut30)

      cutTreeMC.TAttach(PFMTCut30,zerobtagOneMuon)
      cutTreeMC.TAttach(PFMTCut30,more_than_zero_btagOneMuon)
      cutTreeMC.TAttach(PFMTCut30,onebtagOneMuon)
      cutTreeMC.TAttach(PFMTCut30,more_than_one_btagOneMuon)
      cutTreeMC.TAttach(PFMTCut30,twobtagOneMuon)
      cutTreeMC.TAttach(PFMTCut30,more_than_two_btagOneMuon)
      
       
      out.append(AddBinedHist(cutTree = cutTreeMC,
      OP = ("Muon_ControlPlots",genericPSet), cut = PFMTCut30,
      htBins = HTBins,TriggerDict = None ,lab ="OneMuon_") )
      
      out.append(AddBinedHist(cutTree = cutTreeMC,
      OP = ("Muon_ControlPlots",genericPSet), cut = zerobtagOneMuon,
      htBins = HTBins,TriggerDict = None ,lab ="btag_zero_OneMuon_") )
      
      out.append(AddBinedHist(cutTree = cutTreeMC,
      OP = ("Muon_ControlPlots",genericPSet), cut = more_than_zero_btagOneMuon,
      htBins = HTBins,TriggerDict = None ,lab ="btag_morethanzero_OneMuon_") )
     
      out.append(AddBinedHist(cutTree = cutTreeMC,
      OP = ("Muon_ControlPlots",genericPSet), cut = onebtagOneMuon,
      htBins = HTBins,TriggerDict = None ,lab ="btag_one_OneMuon_") )
      
      out.append(AddBinedHist(cutTree = cutTreeMC,
      OP = ("Muon_ControlPlots",genericPSet), cut = more_than_one_btagOneMuon,
      htBins = HTBins,TriggerDict = None,lab ="btag_morethanone_OneMuon_") )
     
      out.append(AddBinedHist(cutTree = cutTreeMC,
      OP = ("Muon_ControlPlots",genericPSet), cut = twobtagOneMuon,
      htBins = HTBins,TriggerDict = None,lab ="btag_two_OneMuon_") )
      
      out.append(AddBinedHist(cutTree = cutTreeMC,
      OP = ("Muon_ControlPlots",genericPSet), cut = more_than_two_btagOneMuon,
      htBins = HTBins,TriggerDict = None,lab ="btag_morethantwo_OneMuon_") )
      
      cutTreeMC.TAttach(minDRMuonJetCut,DiMuon)
      cutTreeMC.TAttach(DiMuon,ZMass_2Muons)
      cutTreeMC.TAttach(ZMass_2Muons,zerobtagDiMuon) 
      cutTreeMC.TAttach(ZMass_2Muons,more_than_zero_btagDiMuon)
      cutTreeMC.TAttach(ZMass_2Muons,onebtagDiMuon) 
      cutTreeMC.TAttach(ZMass_2Muons,more_than_one_btagDiMuon)
      cutTreeMC.TAttach(ZMass_2Muons,twobtagDiMuon) 
      cutTreeMC.TAttach(ZMass_2Muons,more_than_two_btagDiMuon)

      # avobe here does one big inclusive bin!
      # Now lets start binning in HT bins
      # So we can HADD the files at the end and get a chorent set to save the book keeping nightmare:
      # we arrange the HT bins so they are not repoduced though out threshold runs.
      
      out.append(AddBinedHist(cutTree = cutTreeMC,
      OP = ("Muon_ControlPlots",genericPSet), cut = zerobtagDiMuon,
      htBins = HTBins,TriggerDict = None,lab ="btag_zero_DiMuon_") )
      
      out.append(AddBinedHist(cutTree = cutTreeMC,
      OP = ("Muon_ControlPlots",genericPSet), cut = more_than_zero_btagDiMuon,
      htBins = HTBins,TriggerDict = None,lab ="btag_morethanzero_DiMuon_") )
      
      out.append(AddBinedHist(cutTree = cutTreeMC,
      OP = ("Muon_ControlPlots",genericPSet), cut = onebtagDiMuon,
      htBins = HTBins,TriggerDict =None,lab ="btag_one_DiMuon_") )
      
      out.append(AddBinedHist(cutTree = cutTreeMC,
      OP = ("Muon_ControlPlots",genericPSet), cut = more_than_one_btagDiMuon,
      htBins = HTBins,TriggerDict =None ,lab ="btag_morethanone_DiMuon_") )
      
      out.append(AddBinedHist(cutTree = cutTreeMC,
      OP = ("Muon_ControlPlots",genericPSet), cut = twobtagDiMuon,
      htBins = HTBins,TriggerDict = None,lab ="btag_two_DiMuon_") )
      
      out.append(AddBinedHist(cutTree = cutTreeMC,
      OP = ("Muon_ControlPlots",genericPSet), cut = more_than_two_btagDiMuon,
      htBins = HTBins,TriggerDict = None,lab ="btag_morethantwo_DiMuon_") )
      
      out.append(AddBinedHist(cutTree = cutTreeMC,
      OP = ("Muon_ControlPlots",genericPSet), cut = ZMass_2Muons,
      htBins = HTBins,TriggerDict = None ,lab = "DiMuon_") )
      
  return (cutTreeMC,secondJetET,Leading_MuPtCut,MHTCut,out)

# Define the custom muon ID

mu_2012 = PSet(
    MuID = "Tight",
    MinPt = 10.,
    MaxEta = 2.1,
    MaxIsolation = 0.12,
    GlobalChi2 = 10,
    MaxGlbTrkDxy = 0.2,
    MinNumTrkLayers = 6,
    Match2GlbMu = 1,
    NumPixelHits = 1,
    MaxInrTrkDz = 0.5
        )

mu_2012_had = PSet(
    MuID = "Tight",
    MinPt = 10.,
    MaxEta = 2.5,
    MaxIsolation = 0.20,
    GlobalChi2 = 10,
    MaxGlbTrkDxy = 0.2,
    MinNumTrkLayers = 6,
    Match2GlbMu = 1,
    NumPixelHits = 1,
    MaxInrTrkDz = 0.5
        )


#========== VERTEX REWEIGHTING ===============

#Weights 5fb
#==================
PU_2012 =[0.71938114345724813, 4.8793827968712344, 19.336858144769185, 55.601718561254643, 6.413706987112838, 12.56598942176015, 14.456984097455708, 16.835794625620935, 18.409651973047882, 18.922815277453029, 15.466631636777855, 13.972038766759695, 10.286053852141583, 8.4852493611852875, 6.980238934813463, 5.3977976509664254, 3.9525488325601588, 2.8569214623279104, 2.1122481044788932, 1.6198072167494533, 1.2878475694227585, 1.0573375725588976, 0.89123310821884294, 0.76446700561815084, 0.66106878083465226, 0.57188193129746079, 0.49340978772733157, 0.4247383512900465, 0.36561481041531174, 0.31552389197137576, 0.27374732209679692, 0.23926469984848467, 0.21102475742277929, 0.18801723574955601, 0.16943401980964218, 0.1544500850291832, 0.14253868229385042, 0.13334854032674334, 0.12670284780585053, 0.12264484353224907, 0.12134445405076119, 0.12324382024377521, 0.12910778104222692, 0.13987549901929422, 0.15736781377502132, 0.18427521484980819, 0.2252866202639415, 0.28790778075085666, 0.38543468468549952, 0.54107386277205594, 0.79758683218329207, 1.2349196466395214, 2.0106623584400425, 3.4463876905038338, 6.2149289276087636, 11.812670465213747, 23.646591535113778, 49.903750554060011, 111.10458647364125, 260.92928276533735] 

#Weights 5fb capped 5
#PU_2012 =[0.74648581023447425, 5.0632269884942307, 5.1883514870266989, 5.1883515277658221, 5.1883521646241562, 5.1883516646796819, 5.1883517454174202, 5.1883517430080222, 5.1883520912304384, 5.1883517967552901, 5.1883518613609381, 5.1883517086891464, 5.1883520210923608, 5.188352061121301, 5.1883519871150803, 5.1883515999794936, 4.101472154505835, 2.9645638382745174, 2.1918328552746034, 1.6808379223434724, 1.3363709150071257, 1.0971756540126885, 0.92481280740751581, 0.79327041478367377, 0.68597639204081984, 0.59342917574671949, 0.51200037550166355, 0.44074157016176718, 0.37939038663775637, 0.32741214327349416, 0.28406151105596794, 0.24827966602569354, 0.2189757159040808, 0.19510130380591367, 0.175817924969356, 0.16026942546237405, 0.14790922651157395, 0.13837282135284051, 0.13147673577506108, 0.12726583496457966, 0.1259164434681673, 0.12788737725228672, 0.13397227708705414, 0.14514568965837787, 0.16329708844883514, 0.19121831119409308, 0.23377492518845813, 0.29875549741645657, 0.39995701233652836, 0.56146035119139892, 0.82763813522438678, 1.281448764648051, 2.0864196820804675, 3.5762399612450304, 5.1883516811595776, 5.1883517365016605, 5.1883518889038127, 5.1883515899197308, 5.1883518145879446, 5.188351465120622] 

#Unweighted
#PU_2012 = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]

#3.8fb
#=========================
#Reweighted
#PU_2012 =[0.85883047116937294, 5.9064623898052426, 23.087683350617308, 64.862128577233833, 7.2813703359289761, 13.873035072002422, 15.546493490650468, 17.695908138827704, 19.001853404024178, 19.279518423441246, 15.633451316754226, 14.072610354774026, 10.358748297945395, 8.5640540170498998, 7.0694631822589393, 5.4873068369641036, 4.0310636301239722, 2.9196939419603201, 2.1596209635721491, 1.6537196978449631, 1.3102367872042495, 1.069851051389924, 0.89523944312463133, 0.76116822480468083, 0.65166921621966167, 0.55768136357432552, 0.47574284445628123, 0.40483459601779997, 0.34448458685019123, 0.29391877530400229, 0.25216661440679916, 0.21800491657614268, 0.19022949362377631, 0.16772455399185435, 0.14960243129782655, 0.13500052094488404, 0.1233519559413432, 0.11426488455535826, 0.10751228763582046, 0.10306170724827679, 0.10098773482110668, 0.10158647973515571, 0.10540595595880151, 0.11311388443248799, 0.12605838709272352, 0.14622627512463407, 0.1771002074036, 0.22422623409386358, 0.29741240593102081, 0.41368308756208, 0.60425666394190847, 0.92714005418892587, 1.4960334380695857, 2.541530966322775, 4.5428802281129981, 8.5593881949733461, 16.986314238647033, 35.541589987037057, 78.459638711665875, 182.71985901761548]

#Capped 10
#PU_2012 =[0.86479717308718107, 5.9474973805883797, 10.069403098468252, 10.069402511199833, 7.3319580132833835, 10.069402780227605, 10.069402831053431, 10.069403453005533, 10.069402309115368, 10.069403221109978, 10.069403277715052, 10.069403902729725, 10.069402868265723, 8.6235528553749585, 7.1185780770660347, 5.5254298004117066, 4.0590695479284147, 2.9399786684147773, 2.1746249935941728, 1.6652088903071414, 1.319339717503798, 1.0772839100211253, 0.9014591850823338, 0.76645644242912936, 0.65619668708317391, 0.56155585514299389, 0.47904809770566853, 0.40764721369444645, 0.34687792481325763, 0.29596080638810418, 0.25391854704154981, 0.21951951608641318, 0.19155112772174698, 0.16888982277580997, 0.1506417961987003, 0.13593844455047269, 0.12420895088597821, 0.11505874343438069, 0.10825923714128735, 0.10377773649811357, 0.10168934870557138, 0.10229225998730454, 0.10613827418355543, 0.11389974269394976, 0.12693418378028115, 0.14724219573196398, 0.17833062616703491, 0.22578403621963342, 0.29947868567704966, 0.41655715250661701, 0.60845480459266255, 0.93358139389339534, 1.5064271874043427, 2.5591884047937437, 4.5744418761603125, 8.6188550816358518, 10.069402948361518, 10.069403004030264, 10.069402798951762, 10.069402905994242]

#Capped 5
#PU_2012 =[0.89249274887113339, 5.1959400912426332, 5.1959403257989125, 5.1959403729916618, 5.1959403497446539, 5.1959403821618526, 5.1959402301081026, 5.1959404786628669, 5.1959398439199278, 5.1959404486757599, 5.1959404850740478, 5.1959406436507081, 5.195940298892431, 5.1959404101339768, 5.1959406691727885, 5.1959404802161888, 4.1890633229477245, 3.0341327776948548, 2.2442683778094623, 1.7185380553260954, 1.3615922349259266, 1.111784438031586, 0.93032886051032249, 0.79100257883558933, 0.67721171159780291, 0.57953994908760498, 0.49438984960734261, 0.42070231131881092, 0.35798685370401068, 0.30543908531681285, 0.26205041666253054, 0.22654973934473771, 0.19768564164605099, 0.17429860043586468, 0.1554661718623214, 0.14029194157530334, 0.1281867977671049, 0.11874355801324162, 0.11172628974710529, 0.10710127311688684, 0.10494600225932146, 0.10556821443868231, 0.10953739596926251, 0.11754744072965657, 0.13099931317923341, 0.15195769240367732, 0.18404174593848138, 0.23301488901118281, 0.3090696225943117, 0.42989760480236405, 0.62794080365887728, 0.96347987203118945, 1.5546711858583306, 2.6411475408460752, 4.7209407983896821, 5.1959403656403458, 5.1959402979450298, 5.1959400167857117, 5.1959400333026062, 5.1959402884845707]
 
#===========================================

#Data
#==========
from data.Run2012.FNAL.Photon_Run2012A_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_zmengJob228 import *
from data.Run2012.FNAL.SinglePhoton_Run2012B_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_zmengJob229 import *
from data.Run2012.FNAL.HTMHT_Run2012B_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_zmengJob228 import *
from data.Run2012.FNAL.HT_Run2012A_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_zmengJob229 import *
from data.Run2012.FNAL.SingleMu_Run2012A_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_zmengJob228 import *
from data.Run2012.FNAL.SingleMu_Run2012B_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_zmengJob229 import *

from data.Run2012.FNAL.HTMHT_Run2012B_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_zmengJob238 import *
from data.Run2012.FNAL.SingleMu_Run2012B_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_zmengJob238 import *
from data.Run2012.FNAL.SinglePhoton_Run2012B_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_zmengJob238 import *
#MC Samples
from montecarlo.Summer12.FNAL.ICHEP.ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_PU_S7_START52_V9_v1_AODSIM_2_gTag_START52_V9B_ICHEP import * 
from montecarlo.Summer12.FNAL.ICHEP.WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_PU_S7_START52_V9_v1_AODSIM_2_gTag_START52_V9B_ICHEP import * 
from montecarlo.Summer12.FNAL.ICHEP.WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_PU_S7_START52_V9_v1_AODSIM_gTag_START52_V9B_ICHEP import * 
from montecarlo.Summer12.FNAL.ICHEP.TTJets_TuneZ2star_8TeV_madgraph_tauola_Summer12_PU_S7_START52_V9_v1_AODSIM_2_gTag_START52_V9B_ICHEP import * 
from montecarlo.Summer12.FNAL.ICHEP.T_t_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_PU_S7_START52_V9_v1_AODSIM_gTag_START52_V9B_ICHEP import * 
from montecarlo.Summer12.FNAL.ICHEP.T_tW_channel_DR_TuneZ2star_8TeV_powheg_tauola_Summer12_PU_S7_START52_V9_v1_AODSIM_2_gTag_START52_V9B_ICHEP import *
from montecarlo.Summer12.FNAL.ICHEP.Tbar_t_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_PU_S7_START52_V9_v1_AODSIM_2_gTag_START52_V9B_ICHEP import *
from montecarlo.Summer12.FNAL.ICHEP.ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_PU_S7_START52_V9_v3_V17_0_taus_0_doTypeIMetReco_1_gTag_START52_V9B_ICHEP import *
from montecarlo.Summer12.FNAL.ICHEP.ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_PU_S7_START52_V9_v1_AODSIM_gTag_START52_V9B_ICHEP import *
from montecarlo.Summer12.FNAL.ICHEP.ZJetsToNuNu_50_HT_100_TuneZ2Star_8TeV_madgraph_Summer12_PU_S7_START52_V9_v1_AODSIM_gTag_START52_V9B_ICHEP import *
from montecarlo.Summer12.FNAL.ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_PU_S7_START52_V9_v1_V17_0_taus_0_doTypeIMetReco_1_clucasJob234 import * 
from montecarlo.Summer12.FNAL.ICHEP.GJets_HT_200To400_8TeV_madgraph_Summer12_PU_S7_START52_V9_v2_AODSIM_2_gTag_START52_V9B_ICHEP import * 
from montecarlo.Summer12.FNAL.GJets_HT_400ToInf_8TeV_madgraph_Summer12_PU_S7_START52_V9_v1_AODSIM_gTag_START52_V9B_ICHEP import * 
from montecarlo.Summer12.FNAL.DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_Summer12_PU_S7_START52_V9_v2_AODSIM_gTag_START52_V9B_ICHEP import * 
from montecarlo.Summer12.FNAL.ICHEP.WJetsToLNu_HT_400ToInf_8TeV_madgraph_Summer12_PU_S7_START52_V9_v1_AODSIM_gTag_START52_V9B_ICHEP import *
from montecarlo.Summer12.FNAL.ICHEP.WJetsToLNu_HT_300To400_8TeV_madgraph_Summer12_PU_S7_START52_V9_v1_AODSIM_2_gTag_START52_V9B_ICHEP import *
from montecarlo.Summer12.FNAL.ICHEP.WJetsToLNu_HT_250To300_8TeV_madgraph_Summer12_PU_S7_START52_V9_v1_V17_0_taus_0_doTypeIMetReco_1_gTag_START52_V9B_ICHEP import *
from montecarlo.Summer12.FNAL.WJetsToLNu_TuneZ2Star_8TeV_madgraph_tarball_Summer12_PU_S7_START52_V9_v1_V17_0_taus_0_doTypeIMetReco_1_partonHT_2_Skim import *

#MC Lists
Top = [TTJets_TuneZ2star_8TeV_madgraph_tauola_Summer12_PU_S7_START52_V9_v1_AODSIM_2_gTag_START52_V9B_ICHEP, T_t_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_PU_S7_START52_V9_v1_AODSIM_gTag_START52_V9B_ICHEP, T_tW_channel_DR_TuneZ2star_8TeV_powheg_tauola_Summer12_PU_S7_START52_V9_v1_AODSIM_2_gTag_START52_V9B_ICHEP, Tbar_t_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_PU_S7_START52_V9_v1_AODSIM_2_gTag_START52_V9B_ICHEP ]
WJets = [WJetsToLNu_HT_250To300_8TeV_madgraph_Summer12_PU_S7_START52_V9_v1_V17_0_taus_0_doTypeIMetReco_1_gTag_START52_V9B_ICHEP, WJetsToLNu_HT_300To400_8TeV_madgraph_Summer12_PU_S7_START52_V9_v1_AODSIM_2_gTag_START52_V9B_ICHEP,WJetsToLNu_TuneZ2Star_8TeV_madgraph_tarball_Summer12_PU_S7_START52_V9_v1_V17_0_taus_0_doTypeIMetReco_1_partonHT_2_Skim, WJetsToLNu_HT_400ToInf_8TeV_madgraph_Summer12_PU_S7_START52_V9_v1_AODSIM_gTag_START52_V9B_ICHEP  ]
ZJets = [ ZJetsToNuNu_50_HT_100_TuneZ2Star_8TeV_madgraph_Summer12_PU_S7_START52_V9_v1_AODSIM_gTag_START52_V9B_ICHEP, ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_PU_S7_START52_V9_v1_V17_0_taus_0_doTypeIMetReco_1_clucasJob234,  ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_PU_S7_START52_V9_v1_AODSIM_gTag_START52_V9B_ICHEP,ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_PU_S7_START52_V9_v3_V17_0_taus_0_doTypeIMetReco_1_gTag_START52_V9B_ICHEP ]
DiBoson = [ ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_PU_S7_START52_V9_v1_AODSIM_2_gTag_START52_V9B_ICHEP, WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_PU_S7_START52_V9_v1_AODSIM_2_gTag_START52_V9B_ICHEP, WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_PU_S7_START52_V9_v1_AODSIM_gTag_START52_V9B_ICHEP ]
DY = [DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_Summer12_PU_S7_START52_V9_v2_AODSIM_gTag_START52_V9B_ICHEP  ] 
Photon = [ GJets_HT_200To400_8TeV_madgraph_Summer12_PU_S7_START52_V9_v2_AODSIM_2_gTag_START52_V9B_ICHEP,GJets_HT_400ToInf_8TeV_madgraph_Summer12_PU_S7_START52_V9_v1_AODSIM_gTag_START52_V9B_ICHEP]

#=========== MC Sample Lists ===============

MC_2012_Low = Top +WJets+ZJets+DiBoson+DY

MC_2012_High =Top+WJets+ZJets+DY+DiBoson

Photon_MC = Photon

#=============================

#===== Data Sample Lists ===========
Had_2012 = [ HT_Run2012A_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_zmengJob229, HTMHT_Run2012B_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_zmengJob228,HTMHT_Run2012B_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_zmengJob238 ]
Muon_2012 =  [ SingleMu_Run2012A_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_zmengJob228, SingleMu_Run2012B_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_zmengJob229, SingleMu_Run2012B_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_zmengJob238  ]
Photon_2012 = [ Photon_Run2012A_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_zmengJob228, SinglePhoton_Run2012B_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_zmengJob229,SinglePhoton_Run2012B_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_zmengJob238  ]

#===================================

def checkSwitches(d) :
 assert d["json"] in ["ICHEP_18June.json","ICHEP_26June.json"]
 
def switches():
  d = {}
  d["reweight_samples"] = [(PU_2012,MC_2012_Low),(PU_2012,MC_2012_High),(PU_2012,Photon_MC)]
  d["data_samples"] = [Had_2012,Muon_2012,Photon_2012] 
  d["json"] = ["ICHEP_18June.json","ICHEP_26June.json"][-1]
  checkSwitches(d)
  return d

#================= JSON FILTER ========
json = JSONFilter("Json Mask", json_to_pset(switches()["json"]))
