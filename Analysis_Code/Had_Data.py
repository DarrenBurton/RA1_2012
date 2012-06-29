#!/usr/bin/env python

import setupSUSY
from libFrameworkSUSY import *
#from libbryn import *
from libHadronic import *
from libOneLepton import *
from icf.core import PSet,Analysis
from time import strftime
from time import strftime
from batchGolden_singlemu import *

from ra1objectid.vbtfElectronId_cff import *
from ra1objectid.vbtfMuonId_cff import *
from ra1objectid.ra3PhotonId_cff import *
from ra1objectid.ra3PhotonId2012_cff import *

vbtfMuonId_cff = Muon_IDFilter( vbtfmuonidps.ps()  )
vbtfElectronIdFilter = Electron_IDFilter( vbtfelectronidWP95ps.ps() )
ra3PhotonIdFilter    = Photon_IDFilter2012( ra3photonid2012ps.ps() )
CustomEleID = Electron_Egamma_Veto()
CustomMuID = OL_TightMuID(mu_2012_had.ps())

cutTree,blah,blah2,l = MakeDataTree(100., Muon = None)

def addCutFlowData(a) :
  a.AddMuonFilter("PreCC",CustomMuID)
  a.AddPhotonFilter("PreCC",ra3PhotonIdFilter)
  a.AddElectronFilter("PreCC",CustomEleID)
  a+=cutTree

# AK5 Calo

conf_ak5_caloData = deepcopy(defaultConfig)
conf_ak5_caloData.Ntuple = deepcopy(ak5_calo)
conf_ak5_caloData.XCleaning = deepcopy(default_cc)
conf_ak5_caloData.Common = deepcopy(default_common)
# conf_ak5_calo.Common.print_out()
anal_ak5_caloData=Analysis("AK5Calo")
addCutFlowData(anal_ak5_caloData)


#from data.Run2011.HT_Run2011AB import *
#from HT_L1OffSet import *
#from CMSSM_Skim import *

# 8 TeV Data
#==================

#===================

#HT_Run2011AB.File = HT_Run2011AB.File[:1]

outdir = "../../results_"+strftime("%d_%b")+"/Data/"


ensure_dir(outdir)

anal_ak5_caloData.Run(outdir,conf_ak5_caloData,switches()["data_samples"][0])
