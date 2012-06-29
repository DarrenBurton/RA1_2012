#!/usr/bin/env python
import setupSUSY
from libFrameworkSUSY import *
#from libbryn import *
from libHadronic import *
from libOneLepton import *
from icf.core import PSet,Analysis
from time import strftime
import icf.utils as Utils
from batchGolden_singlemu import *
from ra1objectid.vbtfElectronId_cff import *
from ra1objectid.vbtfMuonId_cff import *
#from ra1objectid.ra3PhotonId_cff import *
from ra1objectid.ra3PhotonId2012_cff import *

#JetSmear = JetSmear(0.1,30)
vbtfMuonId_cff = Muon_IDFilter( vbtfmuonidps.ps()  )
CustomEleID = Electron_Egamma_Veto()
cutTreeMC,junkVar,junkVar2,junkVar3,l = MakeMCTree(100.,Muon = None)
CustomMuID = OL_TightMuID(mu_2012_had.ps())
ra3PhotonIdFilter    = Photon_IDFilter2012( ra3photonid2012ps.ps() )

#CHOOSE SAMPLE NUMBER
#0 = MC_2012_Low #1 = MC_2012_High #2=Photon #3=WJet Binned #4=WJet_Inclusive, #5=DY_Zmimic
number = 1
#==============

#vertex_reweight = GoodVertexReweighting(PSet(GoodVertexWeights = switches()["reweight_samples"][number][0]).ps())
pileup_reweight = PileUpReweighting(PSet(PileUpWeights =switches()["reweight_samples"][number][0] ).ps())

def addCutFlowMC(b) :
  #b.AddWeightFilter("Weight", vertex_reweight)
  b.AddWeightFilter("Weight", pileup_reweight)
  b.AddMuonFilter("PreCC",CustomMuID)
  b.AddPhotonFilter("PreCC",ra3PhotonIdFilter)
  b.AddElectronFilter("PreCC",CustomEleID)
  b+=cutTreeMC

#AK5 Calo
conf_ak5_caloMC = deepcopy(defaultConfig)
conf_ak5_caloMC.Ntuple = deepcopy(ak5_calo)
conf_ak5_caloMC.XCleaning = deepcopy(default_cc)
conf_ak5_caloMC.Common = deepcopy(default_common)
conf_ak5_caloMC.Common.print_out()
anal_ak5_caloMC=Analysis("AK5Calo")
addCutFlowMC(anal_ak5_caloMC)

outDir = "../../results_"+strftime("%d_%b")+"//NoSmear/"
ensure_dir(outDir)

#from CMSSM_Skim import *
anal_ak5_caloMC.Run(outDir,conf_ak5_caloMC,switches()["reweight_samples"][number][1])

