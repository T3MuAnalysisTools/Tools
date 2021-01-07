
// Fill DsPhiPiTree
#include "DsPhiPiTree.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace std;

DsPhiPiTree::DsPhiPiTree(TString Name_, TString id_):
   Selection(Name_,id_),
   PhiMassLow(1.00),
   PhiMassHigh(1.04),
   sidebandDsMin(1.70),
   sidebandDsMax(1.80),
   peakDsMin(1.92),
   peakDsMax(2.01),
   dsMassMin(1.68),
   dsMassMax(2.1)
{

   TString basedir = "";
   basedir = (TString)std::getenv("workdir")+"/Code/CommonFiles/PileUp/Collisions2018";
   PUWeightFile = new TFile(basedir+"/PUWeights_Run2018.root");
   puWeights = (TH1D*)PUWeightFile->Get("h1_weights");
}

DsPhiPiTree::~DsPhiPiTree(){
   for(unsigned int j=0; j<Npassed.size(); j++){
      Logger(Logger::Info) << "Selection Summary before: "
         << Npassed.at(j).GetBinContent(1)  << " +/- " << Npassed.at(j).GetBinError(1)  << " after: "
         << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
   }


   Logger(Logger::Info) << "complete." << std::endl;
}

void  DsPhiPiTree::Configure(){


   // Create ntuple branches
   DsPhiPi_Tree= new TTree("tree","tree");
   
   DsPhiPi_Tree->Branch("eventWeight",&eventWeight);
   DsPhiPi_Tree->Branch("genMatcheddR", &genMatcheddR);
   DsPhiPi_Tree->Branch("era",&era);
   DsPhiPi_Tree->Branch("dataMCType",dataMCType);
   DsPhiPi_Tree->Branch("var_relPt_iso05",&var_relPt_iso05);
   DsPhiPi_Tree->Branch("var_nTracks_iso05",&var_nTracks_iso05);
   DsPhiPi_Tree->Branch("var_mindca_iso",&var_mindca_iso);
   DsPhiPi_Tree->Branch("Muon1_TriggerMatchdR",&Muon1_TriggerMatchdR);
   DsPhiPi_Tree->Branch("Muon2_TriggerMatchdR",&Muon2_TriggerMatchdR);
   DsPhiPi_Tree->Branch("Track_TriggerMatchdR",&Track_TriggerMatchdR);
   DsPhiPi_Tree->Branch("Muon1_isGlobal",&Muon1_isGlobal);
   DsPhiPi_Tree->Branch("Muon2_isGlobal",&Muon2_isGlobal);
   DsPhiPi_Tree->Branch("Muon1_isStandAlone",&Muon1_isStandAlone);
   DsPhiPi_Tree->Branch("Muon2_isStandAlone",&Muon2_isStandAlone);
   DsPhiPi_Tree->Branch("Muon1_isTracker",&Muon1_isTracker);
   DsPhiPi_Tree->Branch("Muon2_isTracker",&Muon2_isTracker);
   DsPhiPi_Tree->Branch("Muon1_isCalo",&Muon1_isCalo);
   DsPhiPi_Tree->Branch("Muon2_isCalo",&Muon2_isCalo);
   DsPhiPi_Tree->Branch("Muon1_isIsolationValid",&Muon1_isIsolationValid);
   DsPhiPi_Tree->Branch("Muon2_isIsolationValid",&Muon2_isIsolationValid);
   DsPhiPi_Tree->Branch("NVtx_woPUWeights",&NVtx_woPUWeights);
   DsPhiPi_Tree->Branch("NVtx",&NVtx);
   DsPhiPi_Tree->Branch("PhiMass",&PhiMass);
   DsPhiPi_Tree->Branch("TripleMass",&TripleMass);
   DsPhiPi_Tree->Branch("Track_Pt",&Track_Pt);
   DsPhiPi_Tree->Branch("Track_Eta",&Track_Eta);
   DsPhiPi_Tree->Branch("Track_Phi",&Track_Phi);
   DsPhiPi_Tree->Branch("TrackP",&TrackP);
   DsPhiPi_Tree->Branch("Track_vx",&Track_vx);
   DsPhiPi_Tree->Branch("Track_vy",&Track_vy);
   DsPhiPi_Tree->Branch("Track_vz",&Track_vz);
   DsPhiPi_Tree->Branch("Track_normalizedChi2",&Track_normalizedChi2);
   DsPhiPi_Tree->Branch("Track_numberOfValidHits",&Track_numberOfValidHits);
   DsPhiPi_Tree->Branch("Track_charge",&Track_charge);
   DsPhiPi_Tree->Branch("Track_dxy",&Track_dxy);
   DsPhiPi_Tree->Branch("Track_dz",&Track_dz);
   DsPhiPi_Tree->Branch("Muon1_Pt",&Muon1_Pt);
   DsPhiPi_Tree->Branch("Muon1_Eta",&Muon1_Eta);
   DsPhiPi_Tree->Branch("Muon1_Phi",&Muon1_Phi);
   DsPhiPi_Tree->Branch("Muon1_P",&Muon1_P);
   DsPhiPi_Tree->Branch("Muon1_vx",&Muon1_vx);
   DsPhiPi_Tree->Branch("Muon1_vy",&Muon1_vy);
   DsPhiPi_Tree->Branch("Muon1_vz",&Muon1_vz);
   DsPhiPi_Tree->Branch("Muon1_cLP",&Muon1_cLP);
   DsPhiPi_Tree->Branch("Muon1_SegmentCompatibility",&Muon1_SegmentCompatibility);
   DsPhiPi_Tree->Branch("Muon1_NumberOfMatches",&Muon1_NumberOfMatches);
   DsPhiPi_Tree->Branch("Muon1_NumberOfMatchedStations",&Muon1_NumberOfMatchedStations);
   DsPhiPi_Tree->Branch("Muon1_kink",&Muon1_kink);
   DsPhiPi_Tree->Branch("Muon2_Pt",&Muon2_Pt);
   DsPhiPi_Tree->Branch("Muon2_Eta",&Muon2_Eta);
   DsPhiPi_Tree->Branch("Muon2_Phi",&Muon2_Phi);
   DsPhiPi_Tree->Branch("Muon2_P",&Muon2_P);
   DsPhiPi_Tree->Branch("Muon2_vx",&Muon2_vx);
   DsPhiPi_Tree->Branch("Muon2_vy",&Muon2_vy);
   DsPhiPi_Tree->Branch("Muon2_vz",&Muon2_vz);
   DsPhiPi_Tree->Branch("Muon2_cLP",&Muon2_cLP);
   DsPhiPi_Tree->Branch("Muon2_SegmentCompatibility",&Muon2_SegmentCompatibility);
   DsPhiPi_Tree->Branch("Muon2_NumberOfMatches",&Muon2_NumberOfMatches);
   DsPhiPi_Tree->Branch("Muon2_NumberOfMatchedStations",&Muon2_NumberOfMatchedStations);
   DsPhiPi_Tree->Branch("Muon2_kink",&Muon2_kink);
   DsPhiPi_Tree->Branch("DimuondR",&DimuondR);
   DsPhiPi_Tree->Branch("Muon1TrkdR",&Muon1TrkdR);
   DsPhiPi_Tree->Branch("Muon2TrkdR",&Muon2TrkdR);
   DsPhiPi_Tree->Branch("VertexKFChi2",&VertexKFChi2);
   DsPhiPi_Tree->Branch("SVPVDsDirAngle",&SVPVDsDirAngle);
   DsPhiPi_Tree->Branch("NtracksClose",&NtracksClose);
   DsPhiPi_Tree->Branch("NSV",&NSV);
   DsPhiPi_Tree->Branch("MinMuon_chi2LocalPosition",&MinMuon_chi2LocalPosition);
   DsPhiPi_Tree->Branch("MinDca",&MinDca);
   DsPhiPi_Tree->Branch("MinD0SigSV",&MinD0SigSV);
   DsPhiPi_Tree->Branch("MinD0SigPV",&MinD0SigPV);
   DsPhiPi_Tree->Branch("MaxVertexPairQuality",&MaxVertexPairQuality);
   DsPhiPi_Tree->Branch("MaxdeltaMuZ",&MaxdeltaMuZ);
   DsPhiPi_Tree->Branch("MaxDca",&MaxDca);
   DsPhiPi_Tree->Branch("MaxD0SigSV",&MaxD0SigSV);
   DsPhiPi_Tree->Branch("MaxD0SigPV",&MaxD0SigPV);
   DsPhiPi_Tree->Branch("FLSignificance",&FLSignificance);
   DsPhiPi_Tree->Branch("TransverseFLSignificance", &TransverseFLSignificance);
   DsPhiPi_Tree->Branch("DecayLength",&DecayLength);
   DsPhiPi_Tree->Branch("var_isPrompt",&var_isPrompt);
   DsPhiPi_Tree->Branch("BSSV_dxy", &bssv_dxy);
   DsPhiPi_Tree->Branch("BSSV_significance", &bssv_significance);
   DsPhiPi_Tree->Branch("PVSV_dxy", &pssv_dxy);
   DsPhiPi_Tree->Branch("PVSV_dz", &pvsv_dz);

   // TrackerMuonId quantities
   DsPhiPi_Tree->Branch("Muon1_InnerNC2", &Muon1_InnerNC2);
   DsPhiPi_Tree->Branch("Muon1_ValidFraction",&Muon1_ValidFraction);
   DsPhiPi_Tree->Branch("Muon1_InnerNValidHits",&Muon1_InnerNValidHits);
   DsPhiPi_Tree->Branch("Muon1_NValidPixelHits",&Muon1_NValidPixelHits);
   DsPhiPi_Tree->Branch("Muon1_NValidTrackerHits",&Muon1_NValidTrackerHits);
   DsPhiPi_Tree->Branch("Muon1_NLostTrackerHits",&Muon1_NLostTrackerHits);
   DsPhiPi_Tree->Branch("Muon1_NLostTrackerHitsInner",&Muon1_NLostTrackerHitsInner);
   DsPhiPi_Tree->Branch("Muon1_NLostTrackerHitsOuter",&Muon1_NLostTrackerHitsOuter);
   DsPhiPi_Tree->Branch("Muon1_PixelLayers",&Muon1_PixelLayers);
   DsPhiPi_Tree->Branch("Muon1_TrackerLayers",&Muon1_TrackerLayers);
   DsPhiPi_Tree->Branch("Muon1_PtErrPt",&Muon1_PtErrPt);
   DsPhiPi_Tree->Branch("Muon1_CaloComp",&Muon1_CaloComp);
   DsPhiPi_Tree->Branch("Muon1_HadS9",&Muon1_HadS9);
   DsPhiPi_Tree->Branch("Muon1_Had",&Muon1_Had);
   DsPhiPi_Tree->Branch("Muon1_EM",&Muon1_EM);
   DsPhiPi_Tree->Branch("Muon1_EMS9",&Muon1_EMS9);
   DsPhiPi_Tree->Branch("Muon1_EMS25",&Muon1_EMS25);
   DsPhiPi_Tree->Branch("Muon1_Kink",&Muon1_Kink);

   DsPhiPi_Tree->Branch("Muon2_InnerNC2", &Muon2_InnerNC2);
   DsPhiPi_Tree->Branch("Muon2_ValidFraction",&Muon2_ValidFraction);
   DsPhiPi_Tree->Branch("Muon2_InnerNValidHits",&Muon2_InnerNValidHits);
   DsPhiPi_Tree->Branch("Muon2_NValidPixelHits",&Muon2_NValidPixelHits);
   DsPhiPi_Tree->Branch("Muon2_NValidTrackerHits",&Muon2_NValidTrackerHits);
   DsPhiPi_Tree->Branch("Muon2_NLostTrackerHits",&Muon2_NLostTrackerHits);
   DsPhiPi_Tree->Branch("Muon2_NLostTrackerHitsInner",&Muon2_NLostTrackerHitsInner);
   DsPhiPi_Tree->Branch("Muon2_NLostTrackerHitsOuter",&Muon2_NLostTrackerHitsOuter);
   DsPhiPi_Tree->Branch("Muon2_PixelLayers",&Muon2_PixelLayers);
   DsPhiPi_Tree->Branch("Muon2_TrackerLayers",&Muon2_TrackerLayers);
   DsPhiPi_Tree->Branch("Muon2_PtErrPt",&Muon2_PtErrPt);
   DsPhiPi_Tree->Branch("Muon2_CaloComp",&Muon2_CaloComp);
   DsPhiPi_Tree->Branch("Muon2_HadS9",&Muon2_HadS9);
   DsPhiPi_Tree->Branch("Muon2_Had",&Muon2_Had);
   DsPhiPi_Tree->Branch("Muon2_EM",&Muon2_EM);
   DsPhiPi_Tree->Branch("Muon2_EMS9",&Muon2_EMS9);
   DsPhiPi_Tree->Branch("Muon2_EMS25",&Muon2_EMS25);
   DsPhiPi_Tree->Branch("Muon2_Kink",&Muon2_Kink);

   DsPhiPi_Tree->Branch("Track_InnerNC2", &Track_InnerNC2);
   DsPhiPi_Tree->Branch("Track_ValidFraction",&Track_ValidFraction);
   DsPhiPi_Tree->Branch("Track_InnerNValidHits",&Track_InnerNValidHits);
   DsPhiPi_Tree->Branch("Track_NValidPixelHits",&Track_NValidPixelHits);
   DsPhiPi_Tree->Branch("Track_NValidTrackerHits",&Track_NValidTrackerHits);
   DsPhiPi_Tree->Branch("Track_NLostTrackerHits",&Track_NLostTrackerHits);
   DsPhiPi_Tree->Branch("Track_NLostTrackerHitsInner",&Track_NLostTrackerHitsInner);
   DsPhiPi_Tree->Branch("Track_NLostTrackerHitsOuter",&Track_NLostTrackerHitsOuter);
   DsPhiPi_Tree->Branch("Track_PixelLayers",&Track_PixelLayers);
   DsPhiPi_Tree->Branch("Track_TrackerLayers",&Track_TrackerLayers);
   DsPhiPi_Tree->Branch("Track_PtErrPt",&Track_PtErrPt);
   DsPhiPi_Tree->Branch("Track_CaloComp",&Track_CaloComp);
   DsPhiPi_Tree->Branch("Track_HadS9",&Track_HadS9);
   DsPhiPi_Tree->Branch("Track_Had",&Track_Had);
   DsPhiPi_Tree->Branch("Track_EM",&Track_EM);
   DsPhiPi_Tree->Branch("Track_EMS9",&Track_EMS9);
   DsPhiPi_Tree->Branch("Track_EMS25",&Track_EMS25);
   DsPhiPi_Tree->Branch("Track_Kink",&Track_Kink);

   // bit to check the status of the track
   DsPhiPi_Tree->Branch("Track_isTrackerMuon", &Track_isTrackerMuon);
   DsPhiPi_Tree->Branch("Track_isGlobalMuon", &Track_isGlobalMuon);

   DsPhiPi_Tree->Branch("Muon1_timeAtIpInOut", &Muon1_timeAtIpInOut);
   DsPhiPi_Tree->Branch("Muon1_timeAtIpInOutErr", &Muon1_timeAtIpInOutErr);
   DsPhiPi_Tree->Branch("Muon2_timeAtIpInOut", &Muon2_timeAtIpInOut);
   DsPhiPi_Tree->Branch("Muon2_timeAtIpInOutErr", &Muon2_timeAtIpInOutErr);
   DsPhiPi_Tree->Branch("Track_timeAtIpInOut", &Muon3_timeAtIpInOut);
   DsPhiPi_Tree->Branch("Track_timeAtIpInOutErr", &Muon3_timeAtIpInOutErr);

   DsPhiPi_Tree->Branch("ds_pt", &ds_pt);
   DsPhiPi_Tree->Branch("ds_eta", &ds_eta);
   DsPhiPi_Tree->Branch("ds_phi", &ds_phi);
   DsPhiPi_Tree->Branch("ds_e", &ds_e);
   DsPhiPi_Tree->Branch("ds_motherPdgId", &ds_motherPdgId);

   // Muon station branches
   // Muon station information (Muon1)
   DsPhiPi_Tree->Branch("Muon1_station1_status", &Muon1_station_vars[0][0]);
   DsPhiPi_Tree->Branch("Muon1_station1_TrackX", &Muon1_station_vars[0][1]);
   DsPhiPi_Tree->Branch("Muon1_station1_TrackY", &Muon1_station_vars[0][2]);
   DsPhiPi_Tree->Branch("Muon1_station1_pullX", &Muon1_station_vars[0][3]);
   DsPhiPi_Tree->Branch("Muon1_station1_pullY", &Muon1_station_vars[0][4]);
   DsPhiPi_Tree->Branch("Muon1_station1_pullDxDz", &Muon1_station_vars[0][5]);
   DsPhiPi_Tree->Branch("Muon1_station1_pullDyDz", &Muon1_station_vars[0][6]);

   DsPhiPi_Tree->Branch("Muon1_station2_status", &Muon1_station_vars[1][0]);
   DsPhiPi_Tree->Branch("Muon1_station2_TrackX", &Muon1_station_vars[1][1]);
   DsPhiPi_Tree->Branch("Muon1_station2_TrackY", &Muon1_station_vars[1][2]);
   DsPhiPi_Tree->Branch("Muon1_station2_pullX", &Muon1_station_vars[1][3]);
   DsPhiPi_Tree->Branch("Muon1_station2_pullY", &Muon1_station_vars[1][4]);
   DsPhiPi_Tree->Branch("Muon1_station2_pullDxDz", &Muon1_station_vars[1][5]);
   DsPhiPi_Tree->Branch("Muon1_station2_pullDyDz", &Muon1_station_vars[1][6]);

   DsPhiPi_Tree->Branch("Muon1_station3_status", &Muon1_station_vars[2][0]);
   DsPhiPi_Tree->Branch("Muon1_station3_TrackX", &Muon1_station_vars[2][1]);
   DsPhiPi_Tree->Branch("Muon1_station3_TrackY", &Muon1_station_vars[2][2]);
   DsPhiPi_Tree->Branch("Muon1_station3_pullX", &Muon1_station_vars[2][3]);
   DsPhiPi_Tree->Branch("Muon1_station3_pullY", &Muon1_station_vars[2][4]);
   DsPhiPi_Tree->Branch("Muon1_station3_pullDxDz", &Muon1_station_vars[2][5]);
   DsPhiPi_Tree->Branch("Muon1_station3_pullDyDz", &Muon1_station_vars[2][6]);

   DsPhiPi_Tree->Branch("Muon1_station4_status", &Muon1_station_vars[3][0]);
   DsPhiPi_Tree->Branch("Muon1_station4_TrackX", &Muon1_station_vars[3][1]);
   DsPhiPi_Tree->Branch("Muon1_station4_TrackY", &Muon1_station_vars[3][2]);
   DsPhiPi_Tree->Branch("Muon1_station4_pullX", &Muon1_station_vars[3][3]);
   DsPhiPi_Tree->Branch("Muon1_station4_pullY", &Muon1_station_vars[3][4]);
   DsPhiPi_Tree->Branch("Muon1_station4_pullDxDz", &Muon1_station_vars[3][5]);
   DsPhiPi_Tree->Branch("Muon1_station4_pullDyDz", &Muon1_station_vars[3][6]);

   DsPhiPi_Tree->Branch("Muon1_station5_status", &Muon1_station_vars[4][0]);
   DsPhiPi_Tree->Branch("Muon1_station5_TrackX", &Muon1_station_vars[4][1]);
   DsPhiPi_Tree->Branch("Muon1_station5_TrackY", &Muon1_station_vars[4][2]);
   DsPhiPi_Tree->Branch("Muon1_station5_pullX", &Muon1_station_vars[4][3]);
   DsPhiPi_Tree->Branch("Muon1_station5_pullY", &Muon1_station_vars[4][4]);
   DsPhiPi_Tree->Branch("Muon1_station5_pullDxDz", &Muon1_station_vars[4][5]);
   DsPhiPi_Tree->Branch("Muon1_station5_pullDyDz", &Muon1_station_vars[4][6]);

   DsPhiPi_Tree->Branch("Muon1_station6_status", &Muon1_station_vars[5][0]);
   DsPhiPi_Tree->Branch("Muon1_station6_TrackX", &Muon1_station_vars[5][1]);
   DsPhiPi_Tree->Branch("Muon1_station6_TrackY", &Muon1_station_vars[5][2]);
   DsPhiPi_Tree->Branch("Muon1_station6_pullX", &Muon1_station_vars[5][3]);
   DsPhiPi_Tree->Branch("Muon1_station6_pullY", &Muon1_station_vars[5][4]);
   DsPhiPi_Tree->Branch("Muon1_station6_pullDxDz", &Muon1_station_vars[5][5]);
   DsPhiPi_Tree->Branch("Muon1_station6_pullDyDz", &Muon1_station_vars[5][6]);

   DsPhiPi_Tree->Branch("Muon1_station7_status", &Muon1_station_vars[6][0]);
   DsPhiPi_Tree->Branch("Muon1_station7_TrackX", &Muon1_station_vars[6][1]);
   DsPhiPi_Tree->Branch("Muon1_station7_TrackY", &Muon1_station_vars[6][2]);
   DsPhiPi_Tree->Branch("Muon1_station7_pullX", &Muon1_station_vars[6][3]);
   DsPhiPi_Tree->Branch("Muon1_station7_pullY", &Muon1_station_vars[6][4]);
   DsPhiPi_Tree->Branch("Muon1_station7_pullDxDz", &Muon1_station_vars[6][5]);
   DsPhiPi_Tree->Branch("Muon1_station7_pullDyDz", &Muon1_station_vars[6][6]);

   DsPhiPi_Tree->Branch("Muon1_station8_status", &Muon1_station_vars[7][0]);
   DsPhiPi_Tree->Branch("Muon1_station8_TrackX", &Muon1_station_vars[7][1]);
   DsPhiPi_Tree->Branch("Muon1_station8_TrackY", &Muon1_station_vars[7][2]);
   DsPhiPi_Tree->Branch("Muon1_station8_pullX", &Muon1_station_vars[7][3]);
   DsPhiPi_Tree->Branch("Muon1_station8_pullY", &Muon1_station_vars[7][4]);
   DsPhiPi_Tree->Branch("Muon1_station8_pullDxDz", &Muon1_station_vars[7][5]);
   DsPhiPi_Tree->Branch("Muon1_station8_pullDyDz", &Muon1_station_vars[7][6]);

   // Muon station information (Muon2)
   DsPhiPi_Tree->Branch("Muon2_station1_status", &Muon2_station_vars[0][0]);
   DsPhiPi_Tree->Branch("Muon2_station1_TrackX", &Muon2_station_vars[0][1]);
   DsPhiPi_Tree->Branch("Muon2_station1_TrackY", &Muon2_station_vars[0][2]);
   DsPhiPi_Tree->Branch("Muon2_station1_pullX", &Muon2_station_vars[0][3]);
   DsPhiPi_Tree->Branch("Muon2_station1_pullY", &Muon2_station_vars[0][4]);
   DsPhiPi_Tree->Branch("Muon2_station1_pullDxDz", &Muon2_station_vars[0][5]);
   DsPhiPi_Tree->Branch("Muon2_station1_pullDyDz", &Muon2_station_vars[0][6]);

   DsPhiPi_Tree->Branch("Muon2_station2_status", &Muon2_station_vars[1][0]);
   DsPhiPi_Tree->Branch("Muon2_station2_TrackX", &Muon2_station_vars[1][1]);
   DsPhiPi_Tree->Branch("Muon2_station2_TrackY", &Muon2_station_vars[1][2]);
   DsPhiPi_Tree->Branch("Muon2_station2_pullX", &Muon2_station_vars[1][3]);
   DsPhiPi_Tree->Branch("Muon2_station2_pullY", &Muon2_station_vars[1][4]);
   DsPhiPi_Tree->Branch("Muon2_station2_pullDxDz", &Muon2_station_vars[1][5]);
   DsPhiPi_Tree->Branch("Muon2_station2_pullDyDz", &Muon2_station_vars[1][6]);

   DsPhiPi_Tree->Branch("Muon2_station3_status", &Muon2_station_vars[2][0]);
   DsPhiPi_Tree->Branch("Muon2_station3_TrackX", &Muon2_station_vars[2][1]);
   DsPhiPi_Tree->Branch("Muon2_station3_TrackY", &Muon2_station_vars[2][2]);
   DsPhiPi_Tree->Branch("Muon2_station3_pullX", &Muon2_station_vars[2][3]);
   DsPhiPi_Tree->Branch("Muon2_station3_pullY", &Muon2_station_vars[2][4]);
   DsPhiPi_Tree->Branch("Muon2_station3_pullDxDz", &Muon2_station_vars[2][5]);
   DsPhiPi_Tree->Branch("Muon2_station3_pullDyDz", &Muon2_station_vars[2][6]);

   DsPhiPi_Tree->Branch("Muon2_station4_status", &Muon2_station_vars[3][0]);
   DsPhiPi_Tree->Branch("Muon2_station4_TrackX", &Muon2_station_vars[3][1]);
   DsPhiPi_Tree->Branch("Muon2_station4_TrackY", &Muon2_station_vars[3][2]);
   DsPhiPi_Tree->Branch("Muon2_station4_pullX", &Muon2_station_vars[3][3]);
   DsPhiPi_Tree->Branch("Muon2_station4_pullY", &Muon2_station_vars[3][4]);
   DsPhiPi_Tree->Branch("Muon2_station4_pullDxDz", &Muon2_station_vars[3][5]);
   DsPhiPi_Tree->Branch("Muon2_station4_pullDyDz", &Muon2_station_vars[3][6]);

   DsPhiPi_Tree->Branch("Muon2_station5_status", &Muon2_station_vars[4][0]);
   DsPhiPi_Tree->Branch("Muon2_station5_TrackX", &Muon2_station_vars[4][1]);
   DsPhiPi_Tree->Branch("Muon2_station5_TrackY", &Muon2_station_vars[4][2]);
   DsPhiPi_Tree->Branch("Muon2_station5_pullX", &Muon2_station_vars[4][3]);
   DsPhiPi_Tree->Branch("Muon2_station5_pullY", &Muon2_station_vars[4][4]);
   DsPhiPi_Tree->Branch("Muon2_station5_pullDxDz", &Muon2_station_vars[4][5]);
   DsPhiPi_Tree->Branch("Muon2_station5_pullDyDz", &Muon2_station_vars[4][6]);

   DsPhiPi_Tree->Branch("Muon2_station6_status", &Muon2_station_vars[5][0]);
   DsPhiPi_Tree->Branch("Muon2_station6_TrackX", &Muon2_station_vars[5][1]);
   DsPhiPi_Tree->Branch("Muon2_station6_TrackY", &Muon2_station_vars[5][2]);
   DsPhiPi_Tree->Branch("Muon2_station6_pullX", &Muon2_station_vars[5][3]);
   DsPhiPi_Tree->Branch("Muon2_station6_pullY", &Muon2_station_vars[5][4]);
   DsPhiPi_Tree->Branch("Muon2_station6_pullDxDz", &Muon2_station_vars[5][5]);
   DsPhiPi_Tree->Branch("Muon2_station6_pullDyDz", &Muon2_station_vars[5][6]);

   DsPhiPi_Tree->Branch("Muon2_station7_status", &Muon2_station_vars[6][0]);
   DsPhiPi_Tree->Branch("Muon2_station7_TrackX", &Muon2_station_vars[6][1]);
   DsPhiPi_Tree->Branch("Muon2_station7_TrackY", &Muon2_station_vars[6][2]);
   DsPhiPi_Tree->Branch("Muon2_station7_pullX", &Muon2_station_vars[6][3]);
   DsPhiPi_Tree->Branch("Muon2_station7_pullY", &Muon2_station_vars[6][4]);
   DsPhiPi_Tree->Branch("Muon2_station7_pullDxDz", &Muon2_station_vars[6][5]);
   DsPhiPi_Tree->Branch("Muon2_station7_pullDyDz", &Muon2_station_vars[6][6]);

   DsPhiPi_Tree->Branch("Muon2_station8_status", &Muon2_station_vars[7][0]);
   DsPhiPi_Tree->Branch("Muon2_station8_TrackX", &Muon2_station_vars[7][1]);
   DsPhiPi_Tree->Branch("Muon2_station8_TrackY", &Muon2_station_vars[7][2]);
   DsPhiPi_Tree->Branch("Muon2_station8_pullX", &Muon2_station_vars[7][3]);
   DsPhiPi_Tree->Branch("Muon2_station8_pullY", &Muon2_station_vars[7][4]);
   DsPhiPi_Tree->Branch("Muon2_station8_pullDxDz", &Muon2_station_vars[7][5]);
   DsPhiPi_Tree->Branch("Muon2_station8_pullDyDz", &Muon2_station_vars[7][6]);

   // TrackerMuonId
   //DsPhiPi_Tree->Branch("Muon1_TrackerMuonId", &Muon1_TrackerMuonId);
   //DsPhiPi_Tree->Branch("Muon2_TrackerMuonId", &Muon2_TrackerMuonId);
   //DsPhiPi_Tree->Branch("Track_TrackerMuonId", &Track_TrackerMuonId);

   for(int i=0; i<NCuts;i++){
      cut.push_back(0);
      value.push_back(0);
      pass.push_back(false);     
      if(i==L1Fired)             cut.at(L1Fired)=1;
      if(i==HLTFired)            cut.at(HLTFired)=1;
      if(i==TwoMuTrkCandidate)  cut.at(TwoMuTrkCandidate)=1;
      if(i==OSMuons)            cut.at(OSMuons)=0.1;
      if(i==Mu1PtCut)           cut.at(Mu1PtCut)=3.0;
      if(i==Mu2PtCut)           cut.at(Mu2PtCut)=3.0;
      if(i==TriggerMatch)    cut.at(TriggerMatch)=1;
      if(i==MuonID)             cut.at(MuonID)=1;
      if(i==MuMuMassCut)        cut.at(MuMuMassCut)=1;
      if(i==TrackPtCut)         cut.at(TrackPtCut)=2.0;
      if(i==NTrackHits)         cut.at(NTrackHits)=6;
      //if(i==ChiSqCut)         cut.at(ChiSqCut)=15; // Remove chi sqaure cut
      if(i==DsMassCut)          cut.at(DsMassCut)=1;
   }
   TString hlabel;
   TString htitle;

   for(int i=0; i<NCuts; i++){
      title.push_back("");
      distindx.push_back(false);
      dist.push_back(std::vector<float>());
      TString c="_Cut_";c+=i;
      if(i==L1Fired){
         title.at(i)="L1 Fired";
         hlabel="DoubleMu/TripleMu fired";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_L1Fired_",htitle,2,-0.5,1.5,hlabel,"Entries"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_L1Fired_",htitle,2,-0.5,1.5,hlabel,"Entries"));
      }
      else if(i==HLTFired){
         title.at(i)="HLT Fired";
         hlabel="DoubleMu3_Trk_Tau3mu";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLTFired_",htitle,2,-0.5,1.5,hlabel,"Entries"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLTFired_",htitle,2,-0.5,1.5,hlabel,"Entries"));
      }
      else if(i==TwoMuTrkCandidate){
         title.at(i)="is dimu+trk candidate";
         hlabel="is 2mu+trk candidate";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TwoMuTrkCandidate_",htitle,19,1.0,20.0,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TwoMuTrkCandidate_",htitle,19,1.0,20.0,hlabel,"Events"));
      }
      else if(i==Mu1PtCut){
         title.at(i)="Mu1 Pt $>$ ";
         title.at(i)+=cut.at(Mu1PtCut);
         title.at(i)+=" GeV";
         hlabel="Muon1 PT, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1PtCut_",htitle,50,0,25,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1PtCut_",htitle,50,0,25,hlabel,"Events"));
      }
      else if(i==OSMuons){
         title.at(i)="abs(Mu1+Mu2 Charge) < 0.1 ";
         hlabel="Sum of muon charges";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OSMuons_",htitle,5,-2.5,2.5,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OSMuons_",htitle,5,-2.5,2.5,hlabel,"Events"));
      }
      else if(i==Mu2PtCut){
         title.at(i)="Mu2 Pt $>$ ";
         title.at(i)+=cut.at(Mu2PtCut);
         title.at(i)+=" GeV";
         hlabel="Muon2 PT, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2PtCut_",htitle,40,0,20,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2PtCut_",htitle,40,0,20,hlabel,"Events"));
      }
      else if(i==MuonID){
         title.at(i)="All mu pass ID";
         hlabel="pass MuonID";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
      else if(i==MuMuMassCut){
         title.at(i)+="MuMuMass";
         title.at(i)+=" GeV";
         hlabel="InvMass(MuMu), GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuMuMassCut_",htitle,40,0.8,1.2,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuMuMassCut_",htitle,40,0.8,1.2,hlabel,"Events"));
      }
      else if(i==TrackPtCut){
         title.at(i)="Trk Pt $>$ ";
         title.at(i)+=cut.at(TrackPtCut);
         title.at(i)+=" GeV";
         hlabel="Track PT, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TrackPtCut_",htitle,40,0,20,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TrackPtCut_",htitle,40,0,20,hlabel,"Events"));
      }
      else if(i==NTrackHits){
         title.at(i)="Number of hits in tracker $>$";
         title.at(i)+=cut.at(NTrackHits);
         hlabel="N Tracker Hits";
         Nminus1.push_back(HConfig.GetTH1D(Name+"_Nminus1_NTrackHits_",htitle,50,0,50,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+"_Nminus0_NTrackHits_",htitle,50,0,50,hlabel,"Events"));
      }
      /*
         else if(i==ChiSqCut){
         title.at(i)="Chi Sq of vertex fitting $<$ ";
         title.at(i)+=cut.at(ChiSqCut);
         hlabel="Chi Sq of vertex fitting";
         Nminus1.push_back(HConfig.GetTH1D(Name+"_Nminus1_ChiSq_",htitle,50,0,50,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+"_Nminus0_ChiSq_",htitle,50,0,50,hlabel,"Events"));
         }
         */
      else if(i==TriggerMatch){
         title.at(i)="Trigger Matching";
         hlabel="TriggerMatching";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
      else if(i==DsMassCut){
         title.at(i)="Ds Mass Cut [1.68, 2.02]";
         hlabel="Ds Mass Cut";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_DsMassCut_",htitle,85,1.4,2.25,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_DsMassCut_",htitle,85,1.4,2.25,hlabel,"Events"));
      }
   }

   // Setup NPassed Histogams
   Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove

   // Setup Extra Histograms

   Selection::ConfigureHistograms(); //do not remove
   HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}

void  DsPhiPiTree::Store_ExtraDist(){

   //////////////////////////////////////////////////////////////////////////////////////////////////////
   // Here you must push back all analysis histograms, otherwise they wont be propagated to the output //
   //////////////////////////////////////////////////////////////////////////////////////////////////////

}

void  DsPhiPiTree::doEvent(){
   unsigned int t;
   int id(Ntp->GetMCID());
   if (id==1){
      if(Ntp->WhichEra(2018).Contains("RunA")){ era = 0; }
      if(Ntp->WhichEra(2018).Contains("RunB")){ era = 1; }
      if(Ntp->WhichEra(2018).Contains("RunC")){ era = 2; }
      if(Ntp->WhichEra(2018).Contains("RunD")){ era = 3; }
   }
   else era=-1;

   bool DEBUG = false;
   bool HLTOk(false);
   bool L1Ok(false);
   bool DoubleMuFired(false);
   bool TripleMuFired(false);

   double wobs=1;
   double w;

   // Isolation variables
   double mindca_iso05 = 99.0;
   double mindca_iso = 99.0;
   double mindca_ds = 99.0;

   double sumPtTracks_mu1 = 0;
   double sumPtTracks_mu3 = 0;
   double sumPtTracks_mu2 = 0;

   double sumPtTracks_ds = 0.;
   double sumPtTracks_iso05 = 0.;

   int nTracks_iso05 = 0;
   int nTracks_ds = 0;

   // Initialize all cut values
   value.at(L1Fired)=0;
   value.at(HLTFired)=0;
   value.at(TwoMuTrkCandidate)=0;
   value.at(Mu1PtCut)=0;
   value.at(Mu2PtCut)=0;
   value.at(OSMuons)=1.0;
   value.at(TrackPtCut)=0;
   value.at(NTrackHits)=0;
   value.at(MuonID)=0;
   value.at(MuMuMassCut)=0;
   value.at(TriggerMatch)=0;
   value.at(DsMassCut)=0;

   random_num = rndm.Rndm();

   if(!Ntp->isData()){
      w = (puWeights->GetBinContent(Ntp->TruthNumberOfInteraction())); // Weight MC according to truth number of vertices
   }
   //  No weights to data
   else{w=1;}

   eventWeight = w;

   if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}


   for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
      TString HLT = Ntp->HLTName(iTrigger);
      if(HLT.Contains("HLT_DoubleMu3_Trk_Tau3mu_v") && Ntp->HLTDecision(iTrigger)  ) { HLTOk = true;}
   }

   for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
      TString L1TriggerName = Ntp->L1Name(il1);
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4") && Ntp->L1Decision(il1)) { DoubleMuFired = true; }
      if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9") && Ntp->L1Decision(il1)) { TripleMuFired = true; }
      //if( id!=1 && random_num>0.3516 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMuFired = true; }
      //if( id==1 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMuFired = true; }
   }

   if (DoubleMuFired || TripleMuFired) L1Ok = true;
   value.at(L1Fired) = L1Ok;
   value.at(HLTFired) = HLTOk;

   pass.at(L1Fired) = ( value.at(L1Fired) == cut.at(L1Fired) );
   pass.at(HLTFired) = ( value.at(HLTFired) == cut.at(HLTFired) );

   value.at(TwoMuTrkCandidate) = Ntp->NTwoMuonsTrack(); // Number of two muons + track candidates
   pass.at(TwoMuTrkCandidate) = ( value.at(TwoMuTrkCandidate) >= cut.at(TwoMuTrkCandidate) );
   // Selection of the best candidate

   vector<unsigned int> selectedIndices;
   vector<unsigned int> candidateRank;
   unsigned int final_idx = 0;
   double minChiSq = 999.0;
   bool status;

   if (Ntp->NTwoMuonsTrack()>0){
      for (size_t j=0; j<Ntp->NTwoMuonsTrack(); ++j){

         unsigned int Muon1_idx = Ntp->TwoMuonsTrackMuonIndices(j).at(0);
         unsigned int Muon2_idx = Ntp->TwoMuonsTrackMuonIndices(j).at(1);
         unsigned int Track_idx = Ntp->TwoMuonsTrackTrackIndex(j).at(0);

         unsigned int Mu1_pt_idx = (Ntp->Muon_P4(Muon1_idx).Pt() > Ntp->Muon_P4(Muon2_idx).Pt() ? Muon1_idx:Muon2_idx);
         unsigned int Mu2_pt_idx = (Ntp->Muon_P4(Muon1_idx).Pt() > Ntp->Muon_P4(Muon2_idx).Pt() ? Muon2_idx:Muon1_idx);

         TLorentzVector mu1_p4 = Ntp->Muon_P4(Mu1_pt_idx);
         TLorentzVector mu2_p4 = Ntp->Muon_P4(Mu2_pt_idx);
         TLorentzVector track_p4 = Ntp->Track_P4(Track_idx);

         double mumutrkmass = (mu1_p4+mu2_p4+track_p4).M();

         value.at(Mu1PtCut) = mu1_p4.Pt();
         value.at(Mu2PtCut) = mu2_p4.Pt();
         value.at(MuonID) = ( Ntp->Muon_isPFMuon(Mu1_pt_idx) && Ntp->Muon_isPFMuon(Mu2_pt_idx) &&
               Ntp->Muon_isGlobalMuon(Mu1_pt_idx) && Ntp->Muon_isGlobalMuon(Mu2_pt_idx));
        /*
         value.at(MuonID) = ( ( (Ntp->Muon_isGlobalMuon(Mu1_pt_idx) && !Ntp->Muon_isGlobalMuon(Mu2_pt_idx) && Ntp->Muon_isTrackerMuon(Mu2_pt_idx) ) || 
           (!Ntp->Muon_isGlobalMuon(Mu1_pt_idx) && Ntp->Muon_isTrackerMuon(Mu1_pt_idx) && Ntp->Muon_isGlobalMuon(Mu2_pt_idx) ) || 
           (!Ntp->Muon_isGlobalMuon(Mu1_pt_idx) && Ntp->Muon_isTrackerMuon(Mu1_pt_idx) &&
            !Ntp->Muon_isGlobalMuon(Mu1_pt_idx) && Ntp->Muon_isTrackerMuon(Mu2_pt_idx) ))); 
            */
         value.at(OSMuons) = abs(Ntp->Muon_charge(Mu1_pt_idx)+Ntp->Muon_charge(Mu2_pt_idx));

         float PhiMass = (mu1_p4+mu2_p4).M();
         value.at(MuMuMassCut) = PhiMass;
         value.at(TrackPtCut) = track_p4.Pt();
         value.at(NTrackHits) = Ntp->Track_numberOfValidHits(Track_idx);
         //value.at(ChiSqCut) = Ntp->TwoMuonsTrack_SV_Chi2(j);

         // Trigger matching
         vector<TLorentzVector> trigobjTriplet;

         for (int i=0; i<Ntp->NTriggerObjects(); i++){
            TString name = Ntp->TriggerObject_name(i);
            if (!(name.Contains("tau3muDisplaced3muFltr"))) continue;
            TLorentzVector tmp;
            tmp.SetPtEtaPhiM(Ntp->TriggerObject_pt(i), Ntp->TriggerObject_eta(i), Ntp->TriggerObject_phi(i), PDG_Var::Muon_mass());
            trigobjTriplet.push_back(tmp);
         }

         std::vector<TLorentzVector> muonTriplet;
         muonTriplet.push_back(Ntp->Muon_P4(Muon1_idx));
         muonTriplet.push_back(Ntp->Muon_P4(Muon2_idx));
         muonTriplet.push_back(Ntp->Track_P4(Track_idx));

         std::pair<bool, std::vector<float>> triggerCheck;
         triggerCheck.first = false;
         if (trigobjTriplet.size()>=3) triggerCheck = Ntp->triggerMatchTriplet(muonTriplet, trigobjTriplet);
         value.at(TriggerMatch) = triggerCheck.first;

         value.at(DsMassCut) = mumutrkmass;


         pass.at(OSMuons) = value.at(OSMuons) < cut.at(OSMuons);
         pass.at(Mu1PtCut) = (value.at(Mu1PtCut) >= cut.at(Mu1PtCut) );
         pass.at(Mu2PtCut) = (value.at(Mu2PtCut) >= cut.at(Mu2PtCut) );
         pass.at(MuonID) = ( value.at(MuonID) == cut.at(MuonID) );
         pass.at(MuMuMassCut) = ( ( value.at(MuMuMassCut)< PhiMassHigh) && (value.at(MuMuMassCut) >= PhiMassLow) );
         pass.at(TrackPtCut) = ( value.at(TrackPtCut) > cut.at(TrackPtCut) );
         pass.at(NTrackHits) = ( value.at(NTrackHits) > cut.at(NTrackHits) );
         //pass.at(ChiSqCut) = ( value.at(ChiSqCut) < cut.at(ChiSqCut) );
         //pass.at(TriggerMatch) = ( value.at(TriggerMatch) == cut.at(TriggerMatch) );
         pass.at(TriggerMatch) = true;
         pass.at(DsMassCut) = ( value.at(DsMassCut) < dsMassMax && value.at(DsMassCut) >= dsMassMin );

         unsigned int score = 0;
         for (unsigned int k=0; k<NCuts; ++k) if (pass.at(k)) score++;

         if (score==NCuts) selectedIndices.push_back(j);
         candidateRank.push_back(score);
      }

      if (selectedIndices.size()==0) final_idx=std::distance(candidateRank.begin(), std::max_element(candidateRank.begin(), candidateRank.end()));
      else{
         for (size_t j=0; j<selectedIndices.size(); ++j){
            int _ = selectedIndices.at(j);
            double tmpchi = Ntp->TwoMuonsTrack_SV_Chi2(_);
            if (tmpchi<minChiSq) {
               minChiSq = tmpchi;
               final_idx = _;
            }
         }
      }

      unsigned int Muon1_idx = Ntp->TwoMuonsTrackMuonIndices(final_idx).at(0);
      unsigned int Muon2_idx = Ntp->TwoMuonsTrackMuonIndices(final_idx).at(1);
      unsigned int Track_idx = Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0);

      unsigned int Mu1_pt_idx = (Ntp->Muon_P4(Muon1_idx).Pt() > Ntp->Muon_P4(Muon2_idx).Pt() ? Muon1_idx:Muon2_idx);
      unsigned int Mu2_pt_idx = (Ntp->Muon_P4(Muon1_idx).Pt() > Ntp->Muon_P4(Muon2_idx).Pt() ? Muon2_idx:Muon1_idx);

      TLorentzVector mu1_p4 = Ntp->Muon_P4(Mu1_pt_idx);
      TLorentzVector mu2_p4 = Ntp->Muon_P4(Mu2_pt_idx);
      TLorentzVector track_p4 = Ntp->Track_P4(Track_idx);

      double mumutrkmass = (mu1_p4+mu2_p4+track_p4).M();

      value.at(Mu1PtCut) = mu1_p4.Pt();
      value.at(Mu2PtCut) = mu2_p4.Pt();
      value.at(MuonID) = ( Ntp->Muon_isPFMuon(Mu1_pt_idx) && Ntp->Muon_isPFMuon(Mu2_pt_idx) &&
            Ntp->Muon_isGlobalMuon(Mu1_pt_idx) && Ntp->Muon_isGlobalMuon(Mu2_pt_idx) );

      value.at(OSMuons) = abs(Ntp->Muon_charge(Mu1_pt_idx)+Ntp->Muon_charge(Mu2_pt_idx));

      float PhiMass = (mu1_p4+mu2_p4).M();
      value.at(MuMuMassCut) = PhiMass;
      value.at(TrackPtCut) = track_p4.Pt();
      value.at(NTrackHits) = Ntp->Track_numberOfValidHits(Track_idx);
      //value.at(ChiSqCut) = Ntp->TwoMuonsTrack_SV_Chi2(final_idx);
      value.at(DsMassCut) = mumutrkmass;
      // Trigger matching
      vector<TLorentzVector> trigobjTriplet;

      for (int i=0; i<Ntp->NTriggerObjects(); i++){
         TString name = Ntp->TriggerObject_name(i);
         if (!(name.Contains("hltTau3muTkVertexFilter"))) continue;
         TLorentzVector tmp;
         tmp.SetPtEtaPhiM(Ntp->TriggerObject_pt(i), Ntp->TriggerObject_eta(i), Ntp->TriggerObject_phi(i), PDG_Var::Muon_mass());
         trigobjTriplet.push_back(tmp);
      }

      //cout<<trigobjTriplet.size()<<endl;
      std::vector<TLorentzVector> muonTriplet;
      muonTriplet.push_back(Ntp->Muon_P4(Muon1_idx));
      muonTriplet.push_back(Ntp->Muon_P4(Muon2_idx));
      muonTriplet.push_back(Ntp->Track_P4(Track_idx));

      std::pair<bool, std::vector<float>> triggerCheck;
      triggerCheck.first = false;
      if (trigobjTriplet.size()>=3) triggerCheck = Ntp->triggerMatchTriplet(muonTriplet, trigobjTriplet);
      value.at(TriggerMatch) = triggerCheck.first;

      if (triggerCheck.first) {
      Muon1_TriggerMatchdR = triggerCheck.second.at(0);
      Muon2_TriggerMatchdR = triggerCheck.second.at(1);
      Track_TriggerMatchdR = triggerCheck.second.at(2);
      }
      else{
      Muon1_TriggerMatchdR = 999.0;
      Muon2_TriggerMatchdR = 999.0;
      Track_TriggerMatchdR = 999.0;
      }

      pass.at(OSMuons) = value.at(OSMuons) < cut.at(OSMuons);
      pass.at(Mu1PtCut) = (value.at(Mu1PtCut) >= cut.at(Mu1PtCut) );
      pass.at(Mu2PtCut) = (value.at(Mu2PtCut) >= cut.at(Mu2PtCut) );
      pass.at(MuonID) = ( value.at(MuonID) == cut.at(MuonID) );
      pass.at(MuMuMassCut) = ( ( value.at(MuMuMassCut)< PhiMassHigh) && (value.at(MuMuMassCut) >= PhiMassLow) );
      pass.at(TrackPtCut) = ( value.at(TrackPtCut) > cut.at(TrackPtCut) );
      pass.at(NTrackHits) = ( value.at(NTrackHits) > cut.at(NTrackHits) );
      //pass.at(ChiSqCut) = ( value.at(ChiSqCut) < cut.at(ChiSqCut) );
      pass.at(TriggerMatch) = ( value.at(TriggerMatch) == cut.at(TriggerMatch) );
      pass.at(DsMassCut) = ( value.at(DsMassCut) >= dsMassMin && value.at(DsMassCut) < dsMassMax );
   }

   status = AnalysisCuts(t,w,wobs);

   if(status){

      unsigned int mu1_idx = Ntp->TwoMuonsTrackMuonIndices(final_idx).at(0);
      unsigned int mu2_idx = Ntp->TwoMuonsTrackMuonIndices(final_idx).at(1);
      Track_idx = Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0);

      Muon_1_idx = (Ntp->Muon_P4(mu1_idx).Pt()>Ntp->Muon_P4(mu1_idx).Pt() ? mu1_idx:mu2_idx);
      Muon_2_idx = (Ntp->Muon_P4(mu2_idx).Pt()>Ntp->Muon_P4(mu2_idx).Pt() ? mu2_idx:mu1_idx);
      Track_idx = Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0);

      TLorentzVector mu1_p4 = Ntp->Muon_P4(Muon_1_idx);
      TLorentzVector mu2_p4 = Ntp->Muon_P4(Muon_2_idx);
      TLorentzVector track_p4 = Ntp->Track_P4(Track_idx);

      TLorentzVector DsLV = mu1_p4 + mu2_p4 + track_p4;
      TLorentzVector DsRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0,true)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1,true)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2,true);

      double phimass = (mu1_p4+mu2_p4).M();
      double dsmass = (mu1_p4+mu2_p4+track_p4).M();
      TVector3 SVPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(final_idx,true),Ntp->Vertex_MatchedPrimaryVertex(final_idx,true));
      double DecayL = SVPV.Mag()*dsmass/DsLV.E();

      int Nvertices(0);

      for(unsigned int l=0; l < Ntp->NSecondaryVertices(); l++){
         TVector3 SVsignalPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(final_idx,true),Ntp->Vertex_MatchedPrimaryVertex(final_idx,true));
         TVector3 SVfakePV = Ntp->SVPVDirection(Ntp->SecondaryVertexPosition(l),Ntp->Vertex_MatchedPrimaryVertex(final_idx,true));
         if(SVfakePV.DeltaR(SVsignalPV) < 1 && (Ntp->Vertex_Signal_KF_pos(final_idx,true) - Ntp->SecondaryVertexPosition(l)).Mag() > 0.05){
            Nvertices++;
         }
      }

      float MinD0Significance = std::min({Ntp->Vertex_d0sig_reco(final_idx,0,true),
            Ntp->Vertex_d0sig_reco(final_idx,1,true),
            Ntp->Vertex_d0sig_reco(final_idx,2,true)});

      float MaxD0Significance = std::max({Ntp->Vertex_d0sig_reco(final_idx,0,true),
            Ntp->Vertex_d0sig_reco(final_idx,1,true),
            Ntp->Vertex_d0sig_reco(final_idx,2,true)});

      float MinD0SVSignificance = std::min({Ntp->Vertex_d0sigSV_reco(final_idx,0,true),
            Ntp->Vertex_d0sigSV_reco(final_idx,1,true),
            Ntp->Vertex_d0sigSV_reco(final_idx,2,true)});

      float MaxD0SVSignificance = std::max({Ntp->Vertex_d0sigSV_reco(final_idx,0,true),
            Ntp->Vertex_d0sigSV_reco(final_idx,1,true),
            Ntp->Vertex_d0sigSV_reco(final_idx,2,true)});

      double maxMuondR = std::max({mu1_p4.DeltaR(DsLV), mu2_p4.DeltaR(DsLV), track_p4.DeltaR(DsLV)});
      double minMuonPt = std::min({mu1_p4.Pt(), mu2_p4.Pt(), track_p4.Pt()});

      ds_e = DsLV.E();
      ds_pt = DsLV.Pt();
      ds_eta = DsLV.Eta();
      ds_phi = DsLV.Phi();

      if (id!=1){
      ds_motherPdgId = Ntp->DsMotherPdgId(final_idx);
      genMatcheddR = Ntp->DsGenMatch(final_idx);
      }
      else {
      ds_motherPdgId = 0;
      genMatcheddR = 999;
      }

      // Isolation algorithm
      for (int it=0; it<Ntp->NIsolationTrack(final_idx, true); it++){
         double dxy_track = Ntp->IsolationTrack_dxySV(final_idx, it, true);
         double dz_track = Ntp->IsolationTrack_dzSV(final_idx, it, true);
         TLorentzVector TrackLV = Ntp->IsolationTrack_p4(final_idx, it, true);
         double dca_fv = TMath::Sqrt(pow(dxy_track, 2)+ pow(dz_track, 2));

         double dr_ds = DsLV.DeltaR(TrackLV);
         double dr_mu1 = mu1_p4.DeltaR(TrackLV);
         double dr_mu2 = mu2_p4.DeltaR(TrackLV);
         double dr_trk = track_p4.DeltaR(TrackLV);

         // Isolation 1
         if ( dca_fv<0.5 && TrackLV.Pt()<0.33*minMuonPt && dr_ds<3*maxMuondR ){
            sumPtTracks_ds += TrackLV.Pt();
            nTracks_ds++;
            if (dca_fv < mindca_ds) mindca_ds = dca_fv;
         }

         // Isolation 2
         if (TrackLV.Pt()<1.0) {
            continue;
         }

         if (dca_fv < mindca_iso) mindca_iso = dca_fv;

         // Isolation 3 (within dR = 0.5 of tau)
         if (dr_ds<0.5 && dca_fv<0.5){
            sumPtTracks_iso05 += TrackLV.Pt();
            nTracks_iso05++;
            if(dca_fv<mindca_iso05) mindca_iso05 = dca_fv;
         }
      }

      // Relative Pt calculation
      double relPt_iso05 = sumPtTracks_ds/DsLV.Pt();

      dataMCType = id;

      var_relPt_iso05 = relPt_iso05;
      var_nTracks_iso05 = nTracks_iso05;
      var_mindca_iso = mindca_iso;


      Muon1_isGlobal = Ntp->Muon_isGlobalMuon(Muon_1_idx);
      Muon2_isGlobal = Ntp->Muon_isGlobalMuon(Muon_2_idx);
      Muon1_isStandAlone = Ntp->Muon_isStandAloneMuon(Muon_1_idx);
      Muon2_isStandAlone = Ntp->Muon_isStandAloneMuon(Muon_2_idx);
      Muon1_isTracker = Ntp->Muon_isTrackerMuon(Muon_1_idx);
      Muon2_isTracker = Ntp->Muon_isTrackerMuon(Muon_2_idx);
      Muon1_isCalo = Ntp->Muon_isCaloMuon(Muon_1_idx);
      Muon2_isCalo = Ntp->Muon_isCaloMuon(Muon_2_idx);
      Muon1_isIsolationValid = Ntp->Muon_isIsolationValid(Muon_1_idx);
      Muon2_isIsolationValid = Ntp->Muon_isIsolationValid(Muon_2_idx);
      NVtx_woPUWeights = Ntp->NVtx();
      NVtx = Ntp->NVtx();
      PhiMass = phimass;
      TripleMass = dsmass;

      Track_Pt = Ntp->Track_P4(Track_idx).Pt();
      Track_Eta = Ntp->Track_P4(Track_idx).Eta();
      Track_Phi = Ntp->Track_P4(Track_idx).Phi();
      TrackP = track_p4.P();
      Track_vx = Ntp->Track_Poca(Track_idx).X();
      Track_vy = Ntp->Track_Poca(Track_idx).Y();
      Track_vz = Ntp->Track_Poca(Track_idx).Z();
      Track_normalizedChi2 = Ntp->Track_normalizedChi2(Track_idx);
      Track_numberOfValidHits = Ntp->Track_numberOfValidHits(Track_idx);
      Track_charge = Ntp->Track_charge(Track_idx);
      Track_dxy = Ntp->Track_dxy(Track_idx);
      Track_dz = Ntp->Track_dz(Track_idx);

      Muon1_Pt = Ntp->Muon_P4(Muon_1_idx).Pt();
      Muon1_Eta = Ntp->Muon_P4(Muon_1_idx).Eta();
      Muon1_Phi = Ntp->Muon_P4(Muon_1_idx).Phi();
      Muon1_P = mu1_p4.P();
      Muon1_vx = Ntp->Muon_Poca(Muon_1_idx).X();
      Muon1_vy = Ntp->Muon_Poca(Muon_1_idx).Y();
      Muon1_vz = Ntp->Muon_Poca(Muon_1_idx).Z();
      Muon1_cLP = Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_1_idx);
      Muon1_SegmentCompatibility = Ntp->Muon_segmentCompatibility(Muon_1_idx);
      Muon1_NumberOfMatches = Ntp->Muon_numberOfMatches(Muon_1_idx);
      Muon1_NumberOfMatchedStations = Ntp->Muon_numberOfMatchedStations(Muon_1_idx);
      Muon1_kink = Ntp->Muon_combinedQuality_trkKink(Muon_1_idx);

      Muon2_Pt = Ntp->Muon_P4(Muon_2_idx).Pt();
      Muon2_Eta = Ntp->Muon_P4(Muon_2_idx).Eta();
      Muon2_Phi = Ntp->Muon_P4(Muon_2_idx).Phi();
      Muon2_P = mu2_p4.P();
      Muon2_vx = Ntp->Muon_Poca(Muon_2_idx).X();
      Muon2_vy = Ntp->Muon_Poca(Muon_2_idx).Y();
      Muon2_vz = Ntp->Muon_Poca(Muon_2_idx).Z();
      Muon2_cLP = Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_2_idx);
      Muon2_SegmentCompatibility = Ntp->Muon_segmentCompatibility(Muon_2_idx);
      Muon2_NumberOfMatches = Ntp->Muon_numberOfMatches(Muon_2_idx);
      Muon2_NumberOfMatchedStations = Ntp->Muon_numberOfMatchedStations(Muon_2_idx);
      Muon2_kink = Ntp->Muon_combinedQuality_trkKink(Muon_2_idx);

      DimuondR = mu1_p4.DeltaR(mu2_p4);
      Muon1TrkdR = mu1_p4.DeltaR(track_p4);
      Muon2TrkdR = mu2_p4.DeltaR(track_p4);

      VertexKFChi2 = Ntp->Vertex_signal_KF_Chi2(final_idx,true);
      SVPVDsDirAngle = SVPV.Angle(DsLV.Vect());
      NtracksClose = nTracks_ds;
      NSV = Nvertices;
      MinMuon_chi2LocalPosition = std::min({Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_1_idx),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_2_idx)});
      MinDca = std::min({Ntp->Vertex_DCA12(final_idx,true),Ntp->Vertex_DCA23(final_idx,true),Ntp->Vertex_DCA31(final_idx,true)});
      MinD0SigSV = MinD0SVSignificance;
      MinD0SigPV = MinD0Significance;
      MaxVertexPairQuality = std::max({Ntp->Vertex_pair_quality(final_idx,true),
            Ntp->Vertex_pair_quality(final_idx,true),
            Ntp->Vertex_pair_quality(final_idx,true)});
      MaxdeltaMuZ = std::max({fabs(Ntp->Muon_Poca(Muon_1_idx).Z()  - Ntp->Muon_Poca(Muon_2_idx).Z()),
            fabs(Ntp->Muon_Poca(Muon_1_idx).Z()  - Ntp->Track_Poca(Track_idx).Z()),
            fabs(Ntp->Muon_Poca(Muon_2_idx).Z()  - Ntp->Track_Poca(Track_idx).Z())});
      MaxDca = std::max({Ntp->Vertex_DCA12(final_idx,true),Ntp->Vertex_DCA23(final_idx,true),Ntp->Vertex_DCA31(final_idx,true)});
      MaxD0SigSV = MaxD0SVSignificance;
      MaxD0SigPV = MaxD0Significance;
      TransverseFLSignificance =  Ntp->TransverseFlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(final_idx,true),
            Ntp->Vertex_PrimaryVertex_Covariance(final_idx,true),
            Ntp->Vertex_Signal_KF_pos(final_idx,true),
            Ntp->Vertex_Signal_KF_Covariance(final_idx,true));
      FLSignificance =  Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(final_idx,true),
            Ntp->Vertex_PrimaryVertex_Covariance(final_idx,true),
            Ntp->Vertex_Signal_KF_pos(final_idx,true),
            Ntp->Vertex_Signal_KF_Covariance(final_idx,true));
      DecayLength = DecayL;

      if (id!=1 && Ntp->isPromptDs()) var_isPrompt = 1;
      else if (id!=1 && !Ntp->isPromptDs()) var_isPrompt = 0;
      else var_isPrompt = -1;

      TVector3 PVSV = Ntp->Vertex_Signal_KF_pos(final_idx, true)-Ntp->Vertex_MatchedPrimaryVertex(final_idx, true);
      bssv_dxy = Ntp->Vertex_signal_KF_BS_2Ddistance(final_idx, true);
      bssv_significance = Ntp->Vertex_signal_KF_BS_significance(final_idx, true);
      pssv_dxy = sqrt(pow(PVSV.x(),2)+pow(PVSV.y(),2));
      pvsv_dz = abs(PVSV.z());

      int TrackMuonMatchedIdx = -1;

      Track_isTrackerMuon = false;
      Track_isGlobalMuon = false;
      Track_muonMatched = false;

      for (unsigned int iMuon=0; iMuon<Ntp->NMuons(); iMuon++){
         TLorentzVector tmp_lv = Ntp->Muon_P4(iMuon);
         float tmp_dr = tmp_lv.DeltaR(Ntp->Track_P4(Track_idx));
         if (tmp_dr<0.01 && (abs(tmp_lv.Pt()-Track_Pt)/tmp_lv.Pt())<0.1) {
            Track_muonMatched = true;
            if (Ntp->Muon_isGlobalMuon(iMuon)) Track_isGlobalMuon = true;
            if (Ntp->Muon_isTrackerMuon(iMuon)) Track_isTrackerMuon = true;
            TrackMuonMatchedIdx = iMuon;
         }
      }

      Muon1_InnerNC2 = Ntp->Muon_innerTrack_normalizedChi2(Muon_1_idx);
      Muon1_ValidFraction = Ntp->Muon_innerTrack_validFraction(Muon_1_idx);
      Muon1_InnerNValidHits = Ntp->Muon_innerTrack_numberofValidHits(Muon_1_idx);
      Muon1_NValidPixelHits = Ntp->Muon_numberofValidPixelHits(Muon_1_idx);
      Muon1_NValidTrackerHits = Ntp->Muon_innerTrack_numberOfValidTrackerHits(Muon_1_idx);
      Muon1_NLostTrackerHits = Ntp->Muon_innerTrack_numberOfLostTrackerHits(Muon_1_idx);
      Muon1_NLostTrackerHitsInner = Ntp->Muon_innerTrack_numberOfLostTrackerInnerHits(Muon_1_idx);
      Muon1_NLostTrackerHitsOuter = Ntp->Muon_innerTrack_numberOfLostTrackerOuterHits(Muon_1_idx);
      Muon1_PixelLayers = Ntp->Muon_innerTrack_pixelLayersWithMeasurement(Muon_1_idx);
      Muon1_TrackerLayers = Ntp->Muon_trackerLayersWithMeasurement(Muon_1_idx);
      Muon1_PtErrPt = Ntp->Muon_ptErrOverPt(Muon_1_idx);
      Muon1_CaloComp = Ntp->Muon_caloCompatibility(Muon_1_idx);
      Muon1_HadS9 = Ntp->Muon_calEnergy_hadS9(Muon_1_idx);
      Muon1_Had = Ntp->Muon_calEnergy_had(Muon_1_idx);
      Muon1_EM = Ntp->Muon_calEnergy_em(Muon_1_idx);
      Muon1_EMS9 = Ntp->Muon_calEnergy_emS9(Muon_1_idx);
      Muon1_EMS25 = Ntp->Muon_calEnergy_emS25(Muon_1_idx);
      Muon1_Kink = Ntp->Muon_combinedQuality_trkKink(Muon_1_idx);

      Muon2_InnerNC2 = Ntp->Muon_innerTrack_normalizedChi2(Muon_2_idx);
      Muon2_ValidFraction = Ntp->Muon_innerTrack_validFraction(Muon_2_idx);
      Muon2_InnerNValidHits = Ntp->Muon_innerTrack_numberofValidHits(Muon_2_idx);
      Muon2_NValidPixelHits = Ntp->Muon_numberofValidPixelHits(Muon_2_idx);
      Muon2_NValidTrackerHits = Ntp->Muon_innerTrack_numberOfValidTrackerHits(Muon_2_idx);
      Muon2_NLostTrackerHits = Ntp->Muon_innerTrack_numberOfLostTrackerHits(Muon_2_idx);
      Muon2_NLostTrackerHitsInner = Ntp->Muon_innerTrack_numberOfLostTrackerInnerHits(Muon_2_idx);
      Muon2_NLostTrackerHitsOuter = Ntp->Muon_innerTrack_numberOfLostTrackerOuterHits(Muon_2_idx);
      Muon2_PixelLayers = Ntp->Muon_innerTrack_pixelLayersWithMeasurement(Muon_2_idx);
      Muon2_TrackerLayers = Ntp->Muon_trackerLayersWithMeasurement(Muon_2_idx);
      Muon2_PtErrPt = Ntp->Muon_ptErrOverPt(Muon_2_idx);
      Muon2_CaloComp = Ntp->Muon_caloCompatibility(Muon_2_idx);
      Muon2_HadS9 = Ntp->Muon_calEnergy_hadS9(Muon_2_idx);
      Muon2_Had = Ntp->Muon_calEnergy_had(Muon_2_idx);
      Muon2_EM = Ntp->Muon_calEnergy_em(Muon_2_idx);
      Muon2_EMS9 = Ntp->Muon_calEnergy_emS9(Muon_2_idx);
      Muon2_EMS25 = Ntp->Muon_calEnergy_emS25(Muon_2_idx);
      Muon2_Kink = Ntp->Muon_combinedQuality_trkKink(Muon_2_idx);

      if (Track_muonMatched){
         Track_InnerNC2 = Ntp->Muon_innerTrack_normalizedChi2(TrackMuonMatchedIdx);
         Track_ValidFraction = Ntp->Muon_innerTrack_validFraction(TrackMuonMatchedIdx);
         Track_InnerNValidHits = Ntp->Muon_innerTrack_numberofValidHits(TrackMuonMatchedIdx);
         Track_NValidPixelHits = Ntp->Muon_numberofValidPixelHits(TrackMuonMatchedIdx);
         Track_NValidTrackerHits = Ntp->Muon_innerTrack_numberOfValidTrackerHits(TrackMuonMatchedIdx);
         Track_NLostTrackerHits = Ntp->Muon_innerTrack_numberOfLostTrackerHits(TrackMuonMatchedIdx);
         Track_NLostTrackerHitsInner = Ntp->Muon_innerTrack_numberOfLostTrackerInnerHits(TrackMuonMatchedIdx);
         Track_NLostTrackerHitsOuter = Ntp->Muon_innerTrack_numberOfLostTrackerOuterHits(TrackMuonMatchedIdx);
         Track_PixelLayers = Ntp->Muon_innerTrack_pixelLayersWithMeasurement(TrackMuonMatchedIdx);
         Track_TrackerLayers = Ntp->Muon_trackerLayersWithMeasurement(TrackMuonMatchedIdx);
         Track_PtErrPt = Ntp->Muon_ptErrOverPt(TrackMuonMatchedIdx);
         Track_CaloComp = Ntp->Muon_caloCompatibility(TrackMuonMatchedIdx);
         Track_HadS9 = Ntp->Muon_calEnergy_hadS9(TrackMuonMatchedIdx);
         Track_Had = Ntp->Muon_calEnergy_had(TrackMuonMatchedIdx);
         Track_EM = Ntp->Muon_calEnergy_em(TrackMuonMatchedIdx);
         Track_EMS9 = Ntp->Muon_calEnergy_emS9(TrackMuonMatchedIdx);
         Track_EMS25 = Ntp->Muon_calEnergy_emS25(TrackMuonMatchedIdx);
         Track_Kink = Ntp->Muon_combinedQuality_trkKink(TrackMuonMatchedIdx);
      }
      else {
         Track_InnerNC2 = -99;
         Track_ValidFraction = -99;
         Track_InnerNValidHits = -99;
         Track_NValidPixelHits = -99;
         Track_NValidTrackerHits = -99;
         Track_NLostTrackerHits = -99;
         Track_NLostTrackerHitsInner = -99;
         Track_NLostTrackerHitsOuter = -99;
         Track_PixelLayers = -99;
         Track_TrackerLayers = -99;
         Track_PtErrPt = -99;
         Track_CaloComp = -99;
         Track_HadS9 = -99;
         Track_Had = -99;
         Track_EM = -99;
         Track_EMS9 = -99;
         Track_EMS25 = -99;
         Track_Kink = -99;
      }
      
      Muon1_timeAtIpInOutErr = Ntp->Muon_timeAtIpInOutErr(Muon_1_idx);
      Muon1_timeAtIpInOut = Ntp->Muon_timeAtIpInOut(Muon_1_idx);
      Muon2_timeAtIpInOutErr = Ntp->Muon_timeAtIpInOutErr(Muon_2_idx);
      Muon2_timeAtIpInOut = Ntp->Muon_timeAtIpInOut(Muon_2_idx);
      if (Track_muonMatched){
      Muon3_timeAtIpInOutErr = Ntp->Muon_timeAtIpInOutErr(TrackMuonMatchedIdx);
      Muon3_timeAtIpInOut = Ntp->Muon_timeAtIpInOut(TrackMuonMatchedIdx);
      }
      else{
      Muon3_timeAtIpInOut = -999.;
      Muon3_timeAtIpInOutErr = -999.;
      }

      //----------------------------------------------------------------------
      // count number of stations crossed by muon
      // Segment compatibility variables
      for (int istation=0; istation<8; istation++){
         int station_status = 0;
         if (Ntp->Muon_TrackX(Muon_1_idx, istation)<999999.0 && Ntp->Muon_TrackY(Muon_1_idx, istation)<999999.0) station_status++;
         if (Ntp->Muon_dX(Muon_1_idx, istation)>=999999.0 && station_status) station_status = -1;
         Muon1_station_vars[istation][0] = station_status;
         Muon1_station_vars[istation][1] = Ntp->Muon_TrackX(Muon_1_idx, istation);
         Muon1_station_vars[istation][2] = Ntp->Muon_TrackY(Muon_1_idx, istation);
         Muon1_station_vars[istation][3] = Ntp->Muon_pullX(Muon_1_idx, istation);
         Muon1_station_vars[istation][4] = Ntp->Muon_pullY(Muon_1_idx, istation);
         Muon1_station_vars[istation][5] = Ntp->Muon_pullDxDz(Muon_1_idx, istation);
         Muon1_station_vars[istation][6] = Ntp->Muon_pullDyDz(Muon_1_idx, istation);

         station_status = 0;
         if (Ntp->Muon_TrackX(Muon_2_idx, istation)<999999.0 && Ntp->Muon_TrackY(Muon_2_idx, istation)<999999.0) station_status++;
         if (Ntp->Muon_dX(Muon_2_idx, istation)>=999999.0 && station_status) station_status = -1;
         Muon2_station_vars[istation][0] = station_status;
         Muon2_station_vars[istation][1] = Ntp->Muon_TrackX(Muon_2_idx, istation);
         Muon2_station_vars[istation][2] = Ntp->Muon_TrackY(Muon_2_idx, istation);
         Muon2_station_vars[istation][3] = Ntp->Muon_pullX(Muon_2_idx, istation);
         Muon2_station_vars[istation][4] = Ntp->Muon_pullY(Muon_2_idx, istation);
         Muon2_station_vars[istation][5] = Ntp->Muon_pullDxDz(Muon_2_idx, istation);
         Muon2_station_vars[istation][6] = Ntp->Muon_pullDyDz(Muon_2_idx, istation);
   }

      DsPhiPi_Tree->Fill();
   }
}

void  DsPhiPiTree::Finish(){

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // This function is called after the event loop and you can code here any analysis with already filled analysis histograms
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   file= new TFile("FillDsPhiPiTree_l1DoubleMu0.root","recreate");
   DsPhiPi_Tree->SetDirectory(file);

   file->Write();
   file->Close();


   Selection::Finish();
}
