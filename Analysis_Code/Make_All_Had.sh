#!/bin/sh

hadd Had_Data.root Data*/*.root
hadd Had_TTbar.root NoSmear*/*TT*.root
hadd Had_WJets.root NoSmear*/*WJet*.root
hadd Had_DY.root NoSmear*/*DY*.root
hadd Had_SingleTop.root NoSmear*/*Tbar*.root NoSmear*/*_T_*.root
hadd Had_DiBoson.root NoSmear*/*WW*.root NoSmear*/*WZ*.root NoSmear*/*ZZ*.root 
hadd Had_Zinv.root NoSmear*/*ZJet*.root


