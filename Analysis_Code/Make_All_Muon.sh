#!/bin/sh

hadd Muon_Data.root MuonData*/*.root
hadd Muon_WJets.root MuonNoSmear*/*WJet*.root
hadd Muon_TTbar.root MuonNoSmear*/*TT*.root
hadd Muon_DY.root MuonNoSmear*/*DY*.root
hadd Muon_SingleTop.root MuonNoSmear*/*Tbar*.root MuonNoSmear*/*_T_*.root
hadd Muon_DiBoson.root MuonNoSmear*/*WW*.root MuonNoSmear*/*WZ*.root MuonNoSmear*/*ZZ*.root
hadd Muon_Zinv.root MuonNoSmear*/*ZJet*.root
