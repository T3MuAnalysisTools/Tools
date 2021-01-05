import ROOT
from ROOT import TMVA
from array import array
import numpy as np

basedir = '../CommonFiles/'

def main():

   Muon1_cLP = array('f', [-1.0]);
   Muon1_cLM = array('f', [-1.0]);
   Muon1_staRelChi2 = array('f', [-1.0]);
   Muon1_trkRelChi2 = array('f', [-1.0]);
   Muon1_glbdEP = array('f', [-1.0]);
   Muon1_trkKink = array('f', [-1.0]);
   Muon1_glbKink = array('f', [-1.0]);
   Muon1_glbTrkP = array('f', [-1.0]);
   Muon1_nTVH = array('f', [-1.0]);
   Muon1_nVPH = array('f', [-1.0]);
   Muon1_vMHC = array('f', [-1.0]);
   Muon1_nMS = array('f', [-1.0]);
   Muon1_segComp = array('f', [-1.0]);
   Muon1_tIpOnOut = array('f', [-1.0]);
   Muon1_glbNChi2 = array('f', [-1.0]);
   Muon1_inner_nChi2 = array('f', [-1.0]);
   Muon1_outer_nChi2 = array('f', [-1.0]);
   Muon1_innner_VF = array('f', [-1.0]);

   mu1_eta = array('f', [-1.0])
   mu1_pt = array('f', [-1.0])
   mu1_phi = array('f', [-1.0])
   mu1_SoftMVA = array('f', [-1.0])

   Muon2_cLP = array('f', [-1.0]);
   Muon2_cLM = array('f', [-1.0]);
   Muon2_staRelChi2 = array('f', [-1.0]);
   Muon2_trkRelChi2 = array('f', [-1.0]);
   Muon2_glbdEP = array('f', [-1.0]);
   Muon2_trkKink = array('f', [-1.0]);
   Muon2_glbKink = array('f', [-1.0]);
   Muon2_glbTrkP = array('f', [-1.0]);
   Muon2_nTVH = array('f', [-1.0]);
   Muon2_nVPH = array('f', [-1.0]);
   Muon2_vMHC = array('f', [-1.0]);
   Muon2_nMS = array('f', [-1.0]);
   Muon2_segComp = array('f', [-1.0]);
   Muon2_tIpOnOut = array('f', [-1.0]);
   Muon2_glbNChi2 = array('f', [-1.0]);
   Muon2_inner_nChi2 = array('f', [-1.0]);
   Muon2_outer_nChi2 = array('f', [-1.0]);
   Muon2_innner_VF = array('f', [-1.0]);

   mu2_eta = array('f', [-1.0])
   mu2_pt = array('f', [-1.0])
   mu2_phi = array('f', [-1.0])
   mu2_SoftMVA = array('f', [-1.0])

   Muon3_cLP = array('f', [-1.0]);
   Muon3_cLM = array('f', [-1.0]);
   Muon3_staRelChi2 = array('f', [-1.0]);
   Muon3_trkRelChi2 = array('f', [-1.0]);
   Muon3_glbdEP = array('f', [-1.0]);
   Muon3_trkKink = array('f', [-1.0]);
   Muon3_glbKink = array('f', [-1.0]);
   Muon3_glbTrkP = array('f', [-1.0]);
   Muon3_nTVH = array('f', [-1.0]);
   Muon3_nVPH = array('f', [-1.0]);
   Muon3_vMHC = array('f', [-1.0]);
   Muon3_nMS = array('f', [-1.0]);
   Muon3_segComp = array('f', [-1.0]);
   Muon3_tIpOnOut = array('f', [-1.0]);
   Muon3_glbNChi2 = array('f', [-1.0]);
   Muon3_inner_nChi2 = array('f', [-1.0]);
   Muon3_outer_nChi2 = array('f', [-1.0]);
   Muon3_innner_VF = array('f', [-1.0]);

   mu3_eta = array('f', [-1.0])
   mu3_pt = array('f', [-1.0])
   mu3_phi = array('f', [-1.0])
   mu3_SoftMVA = array('f', [-1.0])

   var_dnnSegCompMuMin = array('f', [-1.0])

   # Muon Id 1
   reader_Muon1Id_barrel = TMVA.Reader("!Color:!Silent");
   reader_Muon1Id_barrel.AddVariable("mu_combinedQuality_chi2LocalMomentum>250?250:mu_combinedQuality_chi2LocalMomentum", Muon1_cLM);
   reader_Muon1Id_barrel.AddVariable("mu_combinedQuality_chi2LocalPosition>50?50:mu_combinedQuality_chi2LocalPosition", Muon1_cLP);
   reader_Muon1Id_barrel.AddVariable("mu_combinedQuality_staRelChi2>50?50:mu_combinedQuality_staRelChi2", Muon1_staRelChi2);
   reader_Muon1Id_barrel.AddVariable("mu_combinedQuality_trkRelChi2>50?50:mu_combinedQuality_trkRelChi2", Muon1_trkRelChi2);
   reader_Muon1Id_barrel.AddVariable("mu_combinedQuality_globalDeltaEtaPhi", Muon1_glbdEP);
   reader_Muon1Id_barrel.AddVariable("log(mu_combinedQuality_trkKink)", Muon1_trkKink);
   reader_Muon1Id_barrel.AddVariable("log(mu_combinedQuality_glbKink)", Muon1_glbKink);
   reader_Muon1Id_barrel.AddVariable("mu_combinedQuality_glbTrackProbability>150?150:mu_combinedQuality_glbTrackProbability", Muon1_glbTrkP);
   reader_Muon1Id_barrel.AddVariable("mu_Numberofvalidpixelhits", Muon1_nVPH);
   reader_Muon1Id_barrel.AddVariable("mu_trackerLayersWithMeasurement", Muon1_nTVH);
   reader_Muon1Id_barrel.AddVariable("mu_validMuonHitComb", Muon1_vMHC);
   reader_Muon1Id_barrel.AddVariable("mu_numberOfMatchedStations", Muon1_nMS);
   reader_Muon1Id_barrel.AddVariable("mu_segmentCompatibility", Muon1_segComp);
   reader_Muon1Id_barrel.AddVariable("mu_timeAtIpInOutErr", Muon1_tIpOnOut);
   reader_Muon1Id_barrel.AddVariable("mu_GLnormChi2>50?50:mu_GLnormChi2", Muon1_glbNChi2);
   reader_Muon1Id_barrel.AddVariable("mu_innerTrack_normalizedChi2>50?50:mu_innerTrack_normalizedChi2", Muon1_inner_nChi2);
   reader_Muon1Id_barrel.AddVariable("mu_outerTrack_normalizedChi2>80?80:mu_outerTrack_normalizedChi2", Muon1_outer_nChi2);
   reader_Muon1Id_barrel.AddVariable("mu_innerTrack_validFraction", Muon1_innner_VF);
   reader_Muon1Id_barrel.AddSpectator("mu_eta",mu1_eta);
   reader_Muon1Id_barrel.AddSpectator("mu_pt",mu1_pt);
   reader_Muon1Id_barrel.AddSpectator("mu_phi",mu1_phi);
   reader_Muon1Id_barrel.AddSpectator("mu_SoftMVA",mu1_SoftMVA);
   reader_Muon1Id_barrel.BookMVA( "BDT", basedir+"weights/weights_barrel/TMVA_new_BDT.weights.xml" ); # weights weights.xml file after training, place it to CommonFiles

   reader_Muon1Id_endcap = TMVA.Reader("!Color:!Silent");
   reader_Muon1Id_endcap.AddVariable("mu_combinedQuality_chi2LocalMomentum>250?250:mu_combinedQuality_chi2LocalMomentum", Muon1_cLM);
   reader_Muon1Id_endcap.AddVariable("mu_combinedQuality_chi2LocalPosition>50?50:mu_combinedQuality_chi2LocalPosition", Muon1_cLP);
   reader_Muon1Id_endcap.AddVariable("mu_combinedQuality_staRelChi2>50?50:mu_combinedQuality_staRelChi2", Muon1_staRelChi2);
   reader_Muon1Id_endcap.AddVariable("mu_combinedQuality_trkRelChi2>50?50:mu_combinedQuality_trkRelChi2", Muon1_trkRelChi2);
   reader_Muon1Id_endcap.AddVariable("mu_combinedQuality_globalDeltaEtaPhi", Muon1_glbdEP);
   reader_Muon1Id_endcap.AddVariable("log(mu_combinedQuality_trkKink)", Muon1_trkKink);
   reader_Muon1Id_endcap.AddVariable("log(mu_combinedQuality_glbKink)", Muon1_glbKink);
   reader_Muon1Id_endcap.AddVariable("mu_combinedQuality_glbTrackProbability>150?150:mu_combinedQuality_glbTrackProbability", Muon1_glbTrkP);
   reader_Muon1Id_endcap.AddVariable("mu_Numberofvalidpixelhits", Muon1_nVPH);
   reader_Muon1Id_endcap.AddVariable("mu_trackerLayersWithMeasurement", Muon1_nTVH);
   reader_Muon1Id_endcap.AddVariable("mu_validMuonHitComb", Muon1_vMHC);
   reader_Muon1Id_endcap.AddVariable("mu_numberOfMatchedStations", Muon1_nMS);
   reader_Muon1Id_endcap.AddVariable("mu_segmentCompatibility", Muon1_segComp);
   reader_Muon1Id_endcap.AddVariable("mu_timeAtIpInOutErr", Muon1_tIpOnOut);
   reader_Muon1Id_endcap.AddVariable("mu_GLnormChi2>50?50:mu_GLnormChi2", Muon1_glbNChi2);
   reader_Muon1Id_endcap.AddVariable("mu_innerTrack_normalizedChi2>50?50:mu_innerTrack_normalizedChi2", Muon1_inner_nChi2);
   reader_Muon1Id_endcap.AddVariable("mu_outerTrack_normalizedChi2>80?80:mu_outerTrack_normalizedChi2", Muon1_outer_nChi2);
   reader_Muon1Id_endcap.AddVariable("mu_innerTrack_validFraction", Muon1_innner_VF);
   reader_Muon1Id_endcap.AddSpectator("mu_eta",mu1_eta);
   reader_Muon1Id_endcap.AddSpectator("mu_pt",mu1_pt);
   reader_Muon1Id_endcap.AddSpectator("mu_phi",mu1_phi);
   reader_Muon1Id_endcap.AddSpectator("mu_SoftMVA",mu1_SoftMVA);
   reader_Muon1Id_endcap.BookMVA( "BDT", basedir+"weights/weights_endcap/TMVA_new_BDT.weights.xml" ); # weights weights.xml file after training, place it to CommonFiles


   # (MuonId 2)
   reader_Muon2Id_barrel = TMVA.Reader("!Color:!Silent");
   reader_Muon2Id_barrel.AddVariable("mu_combinedQuality_chi2LocalMomentum>250?250:mu_combinedQuality_chi2LocalMomentum", Muon2_cLM);
   reader_Muon2Id_barrel.AddVariable("mu_combinedQuality_chi2LocalPosition>50?50:mu_combinedQuality_chi2LocalPosition", Muon2_cLP);
   reader_Muon2Id_barrel.AddVariable("mu_combinedQuality_staRelChi2>50?50:mu_combinedQuality_staRelChi2", Muon2_staRelChi2);
   reader_Muon2Id_barrel.AddVariable("mu_combinedQuality_trkRelChi2>50?50:mu_combinedQuality_trkRelChi2", Muon2_trkRelChi2);
   reader_Muon2Id_barrel.AddVariable("mu_combinedQuality_globalDeltaEtaPhi", Muon2_glbdEP);
   reader_Muon2Id_barrel.AddVariable("log(mu_combinedQuality_trkKink)", Muon2_trkKink);
   reader_Muon2Id_barrel.AddVariable("log(mu_combinedQuality_glbKink)", Muon2_glbKink);
   reader_Muon2Id_barrel.AddVariable("mu_combinedQuality_glbTrackProbability>150?150:mu_combinedQuality_glbTrackProbability", Muon2_glbTrkP);
   reader_Muon2Id_barrel.AddVariable("mu_Numberofvalidpixelhits", Muon2_nVPH);
   reader_Muon2Id_barrel.AddVariable("mu_trackerLayersWithMeasurement", Muon2_nTVH);
   reader_Muon2Id_barrel.AddVariable("mu_validMuonHitComb", Muon2_vMHC);
   reader_Muon2Id_barrel.AddVariable("mu_numberOfMatchedStations", Muon2_nMS);
   reader_Muon2Id_barrel.AddVariable("mu_segmentCompatibility", Muon2_segComp);
   reader_Muon2Id_barrel.AddVariable("mu_timeAtIpInOutErr", Muon2_tIpOnOut);
   reader_Muon2Id_barrel.AddVariable("mu_GLnormChi2>50?50:mu_GLnormChi2", Muon2_glbNChi2);
   reader_Muon2Id_barrel.AddVariable("mu_innerTrack_normalizedChi2>50?50:mu_innerTrack_normalizedChi2", Muon2_inner_nChi2);
   reader_Muon2Id_barrel.AddVariable("mu_outerTrack_normalizedChi2>80?80:mu_outerTrack_normalizedChi2", Muon2_outer_nChi2);
   reader_Muon2Id_barrel.AddVariable("mu_innerTrack_validFraction", Muon2_innner_VF);
   reader_Muon2Id_barrel.AddSpectator("mu_eta",mu2_eta);
   reader_Muon2Id_barrel.AddSpectator("mu_pt",mu2_pt);
   reader_Muon2Id_barrel.AddSpectator("mu_phi",mu2_phi);
   reader_Muon2Id_barrel.AddSpectator("mu_SoftMVA",mu2_SoftMVA);
   reader_Muon2Id_barrel.BookMVA( "BDT", basedir+"/weights/weights_barrel/TMVA_new_BDT.weights.xml" ); # weights weights.xml file after training, place it to CommonFiles


   reader_Muon2Id_endcap = TMVA.Reader("!Color:!Silent");
   reader_Muon2Id_endcap.AddVariable("mu_combinedQuality_chi2LocalMomentum>250?250:mu_combinedQuality_chi2LocalMomentum", Muon2_cLM);
   reader_Muon2Id_endcap.AddVariable("mu_combinedQuality_chi2LocalPosition>50?50:mu_combinedQuality_chi2LocalPosition", Muon2_cLP);
   reader_Muon2Id_endcap.AddVariable("mu_combinedQuality_staRelChi2>50?50:mu_combinedQuality_staRelChi2", Muon2_staRelChi2);
   reader_Muon2Id_endcap.AddVariable("mu_combinedQuality_trkRelChi2>50?50:mu_combinedQuality_trkRelChi2", Muon2_trkRelChi2);
   reader_Muon2Id_endcap.AddVariable("mu_combinedQuality_globalDeltaEtaPhi", Muon2_glbdEP);
   reader_Muon2Id_endcap.AddVariable("log(mu_combinedQuality_trkKink)", Muon2_trkKink);
   reader_Muon2Id_endcap.AddVariable("log(mu_combinedQuality_glbKink)", Muon2_glbKink);
   reader_Muon2Id_endcap.AddVariable("mu_combinedQuality_glbTrackProbability>150?150:mu_combinedQuality_glbTrackProbability", Muon2_glbTrkP);
   reader_Muon2Id_endcap.AddVariable("mu_Numberofvalidpixelhits", Muon2_nVPH);
   reader_Muon2Id_endcap.AddVariable("mu_trackerLayersWithMeasurement", Muon2_nTVH);
   reader_Muon2Id_endcap.AddVariable("mu_validMuonHitComb", Muon2_vMHC);
   reader_Muon2Id_endcap.AddVariable("mu_numberOfMatchedStations", Muon2_nMS);
   reader_Muon2Id_endcap.AddVariable("mu_segmentCompatibility", Muon2_segComp);
   reader_Muon2Id_endcap.AddVariable("mu_timeAtIpInOutErr", Muon2_tIpOnOut);
   reader_Muon2Id_endcap.AddVariable("mu_GLnormChi2>50?50:mu_GLnormChi2", Muon2_glbNChi2);
   reader_Muon2Id_endcap.AddVariable("mu_innerTrack_normalizedChi2>50?50:mu_innerTrack_normalizedChi2", Muon2_inner_nChi2);
   reader_Muon2Id_endcap.AddVariable("mu_outerTrack_normalizedChi2>80?80:mu_outerTrack_normalizedChi2", Muon2_outer_nChi2);
   reader_Muon2Id_endcap.AddVariable("mu_innerTrack_validFraction", Muon2_innner_VF);
   reader_Muon2Id_endcap.AddSpectator("mu_eta",mu2_eta);
   reader_Muon2Id_endcap.AddSpectator("mu_pt",mu2_pt);
   reader_Muon2Id_endcap.AddSpectator("mu_phi",mu2_phi);
   reader_Muon2Id_endcap.AddSpectator("mu_SoftMVA",mu2_SoftMVA);
   reader_Muon2Id_endcap.BookMVA( "BDT", basedir+"/weights/weights_endcap/TMVA_new_BDT.weights.xml" ); # weights weights.xml file after training, place it to CommonFiles

   # (MuonId 3)
   reader_Muon3Id_barrel = TMVA.Reader("!Color:!Silent");
   reader_Muon3Id_barrel.AddVariable("mu_combinedQuality_chi2LocalMomentum>250?250:mu_combinedQuality_chi2LocalMomentum", Muon3_cLM);
   reader_Muon3Id_barrel.AddVariable("mu_combinedQuality_chi2LocalPosition>50?50:mu_combinedQuality_chi2LocalPosition", Muon3_cLP);
   reader_Muon3Id_barrel.AddVariable("mu_combinedQuality_staRelChi2>50?50:mu_combinedQuality_staRelChi2", Muon3_staRelChi2);
   reader_Muon3Id_barrel.AddVariable("mu_combinedQuality_trkRelChi2>50?50:mu_combinedQuality_trkRelChi2", Muon3_trkRelChi2);
   reader_Muon3Id_barrel.AddVariable("mu_combinedQuality_globalDeltaEtaPhi", Muon3_glbdEP);
   reader_Muon3Id_barrel.AddVariable("log(mu_combinedQuality_trkKink)", Muon3_trkKink);
   reader_Muon3Id_barrel.AddVariable("log(mu_combinedQuality_glbKink)", Muon3_glbKink);
   reader_Muon3Id_barrel.AddVariable("mu_combinedQuality_glbTrackProbability>150?150:mu_combinedQuality_glbTrackProbability", Muon3_glbTrkP);
   reader_Muon3Id_barrel.AddVariable("mu_Numberofvalidpixelhits", Muon3_nVPH);
   reader_Muon3Id_barrel.AddVariable("mu_trackerLayersWithMeasurement", Muon3_nTVH);
   reader_Muon3Id_barrel.AddVariable("mu_validMuonHitComb", Muon3_vMHC);
   reader_Muon3Id_barrel.AddVariable("mu_numberOfMatchedStations", Muon3_nMS);
   reader_Muon3Id_barrel.AddVariable("mu_segmentCompatibility", Muon3_segComp);
   reader_Muon3Id_barrel.AddVariable("mu_timeAtIpInOutErr", Muon3_tIpOnOut);
   reader_Muon3Id_barrel.AddVariable("mu_GLnormChi2>50?50:mu_GLnormChi2", Muon3_glbNChi2);
   reader_Muon3Id_barrel.AddVariable("mu_innerTrack_normalizedChi2>50?50:mu_innerTrack_normalizedChi2", Muon3_inner_nChi2);
   reader_Muon3Id_barrel.AddVariable("mu_outerTrack_normalizedChi2>80?80:mu_outerTrack_normalizedChi2", Muon3_outer_nChi2);
   reader_Muon3Id_barrel.AddVariable("mu_innerTrack_validFraction", Muon3_innner_VF);
   reader_Muon3Id_barrel.AddSpectator("mu_eta",mu2_eta);
   reader_Muon3Id_barrel.AddSpectator("mu_pt",mu2_pt);
   reader_Muon3Id_barrel.AddSpectator("mu_phi",mu2_phi);
   reader_Muon3Id_barrel.AddSpectator("mu_SoftMVA",mu2_SoftMVA);
   reader_Muon3Id_barrel.BookMVA( "BDT", basedir+"/weights/weights_barrel/TMVA_new_BDT.weights.xml" ); # weights weights.xml file after training, place it to CommonFiles

   reader_Muon3Id_endcap = TMVA.Reader("!Color:!Silent");
   reader_Muon3Id_endcap.AddVariable("mu_combinedQuality_chi2LocalMomentum>250?250:mu_combinedQuality_chi2LocalMomentum", Muon3_cLM);
   reader_Muon3Id_endcap.AddVariable("mu_combinedQuality_chi2LocalPosition>50?50:mu_combinedQuality_chi2LocalPosition", Muon3_cLP);
   reader_Muon3Id_endcap.AddVariable("mu_combinedQuality_staRelChi2>50?50:mu_combinedQuality_staRelChi2", Muon3_staRelChi2);
   reader_Muon3Id_endcap.AddVariable("mu_combinedQuality_trkRelChi2>50?50:mu_combinedQuality_trkRelChi2", Muon3_trkRelChi2);
   reader_Muon3Id_endcap.AddVariable("mu_combinedQuality_globalDeltaEtaPhi", Muon3_glbdEP);
   reader_Muon3Id_endcap.AddVariable("log(mu_combinedQuality_trkKink)", Muon3_trkKink);
   reader_Muon3Id_endcap.AddVariable("log(mu_combinedQuality_glbKink)", Muon3_glbKink);
   reader_Muon3Id_endcap.AddVariable("mu_combinedQuality_glbTrackProbability>150?150:mu_combinedQuality_glbTrackProbability", Muon3_glbTrkP);
   reader_Muon3Id_endcap.AddVariable("mu_Numberofvalidpixelhits", Muon3_nVPH);
   reader_Muon3Id_endcap.AddVariable("mu_trackerLayersWithMeasurement", Muon3_nTVH);
   reader_Muon3Id_endcap.AddVariable("mu_validMuonHitComb", Muon3_vMHC);
   reader_Muon3Id_endcap.AddVariable("mu_numberOfMatchedStations", Muon3_nMS);
   reader_Muon3Id_endcap.AddVariable("mu_segmentCompatibility", Muon3_segComp);
   reader_Muon3Id_endcap.AddVariable("mu_timeAtIpInOutErr", Muon3_tIpOnOut);
   reader_Muon3Id_endcap.AddVariable("mu_GLnormChi2>50?50:mu_GLnormChi2", Muon3_glbNChi2);
   reader_Muon3Id_endcap.AddVariable("mu_innerTrack_normalizedChi2>50?50:mu_innerTrack_normalizedChi2", Muon3_inner_nChi2);
   reader_Muon3Id_endcap.AddVariable("mu_outerTrack_normalizedChi2>80?80:mu_outerTrack_normalizedChi2", Muon3_outer_nChi2);
   reader_Muon3Id_endcap.AddVariable("mu_innerTrack_validFraction", Muon3_innner_VF);
   reader_Muon3Id_endcap.AddSpectator("mu_eta",mu3_eta);
   reader_Muon3Id_endcap.AddSpectator("mu_pt",mu3_pt);
   reader_Muon3Id_endcap.AddSpectator("mu_phi",mu3_phi);
   reader_Muon3Id_endcap.AddSpectator("mu_SoftMVA",mu3_SoftMVA);
   reader_Muon3Id_endcap.BookMVA( "BDT", basedir+"/weights/weights_endcap/TMVA_new_BDT.weights.xml" ); # weights weights.xml file after training, place it to CommonFiles
   file_ = ROOT.TFile("/eos/user/b/bjoshi/RunIITau23Mu/AnalysisTrees/T3MSelectionTreeInput_combined_with_keras_scores.root","READ")

   old_tree_bkg = file_.Get('TreeB')
   old_tree_ds = file_.Get('TreeS_Ds')
   old_tree_bu = file_.Get('TreeS_Bu')
   old_tree_bd = file_.Get('TreeS_Bd')

   new_file = ROOT.TFile("/eos/user/b/bjoshi/RunIITau23Mu/AnalysisTrees/T3MSelectionTreeInput_combined_globalMuonId.root","RECREATE")

   new_tree_bkg = old_tree_bkg.CloneTree(0)
   new_tree_ds = old_tree_ds.CloneTree(0)
   new_tree_bu = old_tree_bu.CloneTree(0)
   new_tree_bd = old_tree_bd.CloneTree(0)

   old_tree_list = [ old_tree_bkg, old_tree_ds, old_tree_bu, old_tree_bd ]
   new_tree_list = [ new_tree_bkg, new_tree_ds, new_tree_bu, new_tree_bd ]
   globalMuon1Id = [ array('f', [-99.]), array('f', [-99.]), array('f', [-99.]), array('f', [-99.]) ]
   globalMuon2Id = [ array('f', [-99.]), array('f', [-99.]), array('f', [-99.]), array('f', [-99.]) ]
   globalMuon3Id = [ array('f', [-99.]), array('f', [-99.]), array('f', [-99.]), array('f', [-99.]) ]
   var_dnnSegCompMuMin = [ array('f', [-99.]), array('f', [-99.]), array('f', [-99.]), array('f', [-99.]) ]

   for i in xrange(4):
       new_tree_list[i].Branch("globalMuon1Id", globalMuon1Id[i], "globalMuon1Id/F")
       new_tree_list[i].Branch("globalMuon2Id", globalMuon2Id[i], "globalMuon2Id/F")
       new_tree_list[i].Branch("globalMuon3Id", globalMuon3Id[i], "globalMuon3Id/F")
       new_tree_list[i].Branch("var_dnnSegCompMuMin", var_dnnSegCompMuMin[i], "var_dnnSegCompMuMin/F")

   for i, tree in enumerate(old_tree_list):
      nentries = tree.GetEntriesFast()
      print tree.GetName()
      for j in xrange(nentries):
         if (j%1000==0): print "Processing ",j,"/",nentries," ..."
         tree.GetEntry(j)
         muon1_score = -99.0
         muon2_score = -99.0
         muon3_score = -99.0
         Muon1_cLP[0] = tree.Muon1_combinedQuality_chi2LocalPosition;
         Muon1_cLM[0] = tree.Muon1_combinedQuality_chi2LocalMomentum;
         Muon1_staRelChi2[0] = tree.Muon1_combinedQuality_staRelChi2;
         Muon1_trkRelChi2[0] = tree.Muon1_combinedQuality_trkRelChi2;
         Muon1_glbdEP[0] = tree.Muon1_combinedQuality_globalDeltaEtaPhi;
         Muon1_trkKink[0] = np.log(0.01+tree.Muon1_combinedQuality_trkKink);
         Muon1_glbKink[0] = np.log(0.01+tree.Muon1_combinedQuality_glbKink);
         Muon1_glbTrkP[0] = tree.Muon1_combinedQuality_glbTrackProbability;
         Muon1_nTVH[0] = tree.Muon1_trackerLayersWithMeasurement;
         Muon1_nVPH[0] = tree.Muon1_numberofValidPixelHits;
         Muon1_vMHC[0] = tree.Muon1_vmuonhitcomb_reco;
         Muon1_nMS[0] = tree.Muon1_numberOfMatchedStations;
         Muon1_segComp[0] = tree.Muon1_segmentCompatibility;
         Muon1_tIpOnOut[0] = tree.var_Muon1_timeAtIpInOutErr;
         Muon1_glbNChi2[0] = tree.Muon1_normChi2;
         Muon1_inner_nChi2[0] = tree.Muon1_innerTrack_normalizedChi2;
         Muon1_outer_nChi2[0] = tree.Muon1_outerTrack_normalizedChi2;
         Muon1_innner_VF[0] = tree.Muon1_innerTrack_validFraction;

         Muon2_cLM[0] = tree.Muon2_combinedQuality_chi2LocalMomentum;
         Muon2_cLP[0] = tree.Muon2_combinedQuality_chi2LocalPosition;
         Muon2_staRelChi2[0] = tree.Muon2_combinedQuality_staRelChi2;
         Muon2_trkRelChi2[0] = tree.Muon2_combinedQuality_trkRelChi2;
         Muon2_glbdEP[0] = tree.Muon2_combinedQuality_globalDeltaEtaPhi;
         Muon2_trkKink[0] = np.log(0.01+tree.Muon2_combinedQuality_trkKink);
         Muon2_glbKink[0] = np.log(0.01+tree.Muon2_combinedQuality_glbKink);
         Muon2_glbTrkP[0] = tree.Muon2_combinedQuality_glbTrackProbability;
         Muon2_nTVH[0] = tree.Muon2_trackerLayersWithMeasurement;
         Muon2_nVPH[0] = tree.Muon2_numberofValidPixelHits;
         Muon2_vMHC[0] = tree.Muon2_vmuonhitcomb_reco;
         Muon2_nMS[0] = tree.Muon2_numberOfMatchedStations;
         Muon2_segComp[0] = tree.Muon2_segmentCompatibility;
         Muon2_tIpOnOut[0] = tree.var_Muon2_timeAtIpInOutErr;
         Muon2_glbNChi2[0] = tree.Muon2_normChi2;
         Muon2_inner_nChi2[0] = tree.Muon2_innerTrack_normalizedChi2;
         Muon2_outer_nChi2[0] = tree.Muon2_outerTrack_normalizedChi2;
         Muon2_innner_VF[0] = tree.Muon2_innerTrack_validFraction;

         Muon3_cLM[0] = tree.Muon3_combinedQuality_chi2LocalMomentum;
         Muon3_cLP[0] = tree.Muon3_combinedQuality_chi2LocalPosition;
         Muon3_staRelChi2[0] = tree.Muon3_combinedQuality_staRelChi2;
         Muon3_trkRelChi2[0] = tree.Muon3_combinedQuality_trkRelChi2;
         Muon3_glbdEP[0] = tree.Muon3_combinedQuality_globalDeltaEtaPhi;
         Muon3_trkKink[0] = np.log(0.01+tree.Muon3_combinedQuality_trkKink);
         Muon3_glbKink[0] = np.log(0.01+tree.Muon3_combinedQuality_glbKink);
         Muon3_glbTrkP[0] = tree.Muon3_combinedQuality_glbTrackProbability;
         Muon3_nTVH[0] = tree.Muon3_trackerLayersWithMeasurement;
         Muon3_nVPH[0] = tree.Muon3_numberofValidPixelHits;
         Muon3_vMHC[0] = tree.Muon3_vmuonhitcomb_reco;
         Muon3_nMS[0] = tree.Muon3_numberOfMatchedStations;
         Muon3_segComp[0] = tree.Muon3_segmentCompatibility;
         Muon3_tIpOnOut[0] = tree.var_Muon3_timeAtIpInOutErr;
         Muon3_glbNChi2[0] = tree.Muon3_normChi2;
         Muon3_inner_nChi2[0] = tree.Muon3_innerTrack_normalizedChi2;
         Muon3_outer_nChi2[0] = tree.Muon3_outerTrack_normalizedChi2;
         Muon3_innner_VF[0] = tree.Muon3_innerTrack_validFraction;

         if (abs(tree.var_Muon1_Eta)<1.2): muon1_score = reader_Muon1Id_barrel.EvaluateMVA("BDT")
         else: muon1_score = reader_Muon1Id_endcap.EvaluateMVA("BDT")

         if (abs(tree.var_Muon2_Eta)<1.2): muon2_score = reader_Muon2Id_barrel.EvaluateMVA("BDT")
         else: muon2_score = reader_Muon2Id_endcap.EvaluateMVA("BDT")

         if (abs(tree.var_Muon3_Eta)<1.2): muon3_score = reader_Muon3Id_barrel.EvaluateMVA("BDT")
         else: muon3_score = reader_Muon3Id_endcap.EvaluateMVA("BDT")

         globalMuon1Id[i][0] = muon1_score
         globalMuon2Id[i][0] = muon2_score
         if (tree.threeGlobal==1): globalMuon3Id[i][0] = muon3_score
         else: globalMuon3Id[i][0] = -99.0
         var_dnnSegCompMuMin[i][0] = min(tree.muon1_seg_comp_dnn, tree.muon2_seg_comp_dnn, tree.muon3_seg_comp_dnn)
         new_tree_list[i].Fill()

   new_file.cd()
   new_tree_bkg.Write()
   new_tree_ds.Write()
   new_tree_bu.Write()
   new_tree_bd.Write()
   new_file.Close()
   file_.Close()

if __name__=='__main__':
    main()
