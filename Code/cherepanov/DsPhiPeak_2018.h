#ifndef DsPhiPeak_2018_h
#define DsPhiPeak_2018_h

#include "Selection.h"
#include <vector>
#include "TString.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

class DsPhiPeak_2018 : public Selection {

 public:
  DsPhiPeak_2018(TString Name_, TString id_);
  virtual ~DsPhiPeak_2018();

  virtual void  Configure();
  virtual void  Finish();

  //  enum cuts {L1TOk=0,HLTOk,is2MuTrk,GlobalMu,MuCharge,Mass2Mu,Mu1dR,Mu2dR,TrkdR,Mu1pt,Mu2pt,Trkpt,NTrackHits,DsMassCut,NCuts}; 
  enum cuts {L1TOk=0,HLTOk,is2MuTrk,GlobalMu,MuCharge,Mass2Mu,TriggerMatch,Mu1pt,Mu2pt,Trkpt,DsMassCut,NCuts};


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

  TFile * file;
  TTree * Sync_tree;

 private:

  // PU Weights
  TFile* PUWeightFile;
  TH1D* puWeights;

  TMVA::Reader *readerMuIDBarrel;
  TMVA::Reader *readerMuIDEndcap;
  TMVA::Reader *reader_trackerMuonId;

  // Selection Variables
  bool RunB, RunC, RunD, RunE, RunF = 0;
  float dsMassMin;
  float dsMassMax;

  TRandom rndm;
  double random_num;

  std::vector<std::vector<TH1D>*> peakCollection;
  std::vector<std::vector<TH1D>*> sidebandCollection;
  std::vector<std::vector<TH1D>*> validationCollection;

  // Initializhere your analysis histograms
  std::vector<TH1D> NVtx;

  std::vector<TH1D> DimuondR;
  std::vector<TH1D> Muon1TrkdR;
  std::vector<TH1D> Muon2TrkdR;
  std::vector<TH1D> PhiMass;
  std::vector<TH1D> PhiPlusTrackMass;
  std::vector<TH2D> PhiMassVsDsMass;
  std::vector<TH1D> DsMass;
  std::vector<TH1D> Category;
  std::vector<TH1D> DsGenMatch;

  std::vector<TH1D> DsPt_peak;
  std::vector<TH1D> DsPt_sideband;
  std::vector<TH1D> Ds_Pt;
  std::vector<TH1D> DsP_peak;
  std::vector<TH1D> DsP_sideband;
  std::vector<TH1D> Ds_P;
  std::vector<TH1D> DsM_peak;
  std::vector<TH1D> DsM_sideband;
  std::vector<TH1D> Ds_M;
  std::vector<TH1D> DsL_peak;
  std::vector<TH1D> DsL_sideband;
  std::vector<TH1D> Ds_L;
  std::vector<TH1D> DsEta_peak;
  std::vector<TH1D> DsEta_sideband;
  std::vector<TH1D> Ds_eta;

  //Muon variables
  
  //Track candidate variables
  std::vector<TH1D> Track_P;
  std::vector<TH1D> Track_E;
  std::vector<TH1D> Track_Pt;
  std::vector<TH1D> Track_Eta;
  std::vector<TH1D> Track_Phi;
  std::vector<TH1D> Track_normalizedChi2;
  std::vector<TH1D> Track_numberOfValidHits;
  std::vector<TH1D> Track_charge;

  std::vector<TH1D> Muon1_Pt;
  std::vector<TH1D> Muon1_E;
  std::vector<TH1D> Muon1_P;
  std::vector<TH1D> Muon1_Eta;
  std::vector<TH1D> Muon1_Phi;

  std::vector<TH1D> Muon2_Pt;
  std::vector<TH1D> Muon2_E;
  std::vector<TH1D> Muon2_P;
  std::vector<TH1D> Muon2_Eta;
  std::vector<TH1D> Muon2_Phi;
		
  std::vector<TH1D> Muon1_isGlobal;
  std::vector<TH1D> Muon2_isGlobal;
  std::vector<TH1D> Muon1_isStandAlone;
  std::vector<TH1D> Muon2_isStandAlone;
  std::vector<TH1D> Muon1_isTracker;
  std::vector<TH1D> Muon2_isTracker;
  std::vector<TH1D> Muon1_isCalo;
  std::vector<TH1D> Muon2_isCalo;
  std::vector<TH1D> Muon1_TriggerMatchdR;
  std::vector<TH1D> Muon2_TriggerMatchdR;
  std::vector<TH1D> Track_TriggerMatchdR;

  std::vector<TH1D> DecayLength_peak;
  std::vector<TH1D> DecayLength_sideband;
  std::vector<TH1D> DecayLength_prompt;
  std::vector<TH1D> DecayLength_non_prompt;
  std::vector<TH1D> DecayLength;

  std::vector<TH1D> Muon1_Pt_peak;
  std::vector<TH1D> Muon1_Pt_sideband;
  std::vector<TH1D> Muon1_Eta_peak;
  std::vector<TH1D> Muon1_Eta_sideband;
  std::vector<TH1D> Muon1_Phi_peak;
  std::vector<TH1D> Muon1_Phi_sideband;
  std::vector<TH1D> control_Muon1_Pt;
  std::vector<TH1D> control_Muon1_Eta;
  std::vector<TH1D> control_Muon1_Phi;

  std::vector<TH1D> Muon2_Pt_peak;
  std::vector<TH1D> Muon2_Pt_sideband;
  std::vector<TH1D> Muon2_Eta_peak;
  std::vector<TH1D> Muon2_Eta_sideband;
  std::vector<TH1D> Muon2_Phi_peak;
  std::vector<TH1D> Muon2_Phi_sideband;
  std::vector<TH1D> control_Muon2_Pt;
  std::vector<TH1D> control_Muon2_Eta;
  std::vector<TH1D> control_Muon2_Phi;

  std::vector<TH1D> Track_Pt_peak;
  std::vector<TH1D> Track_Pt_sideband;
  std::vector<TH1D> Track_Eta_peak;
  std::vector<TH1D> Track_Eta_sideband;
  std::vector<TH1D> Track_Phi_peak;
  std::vector<TH1D> Track_Phi_sideband;
  std::vector<TH1D> control_Track_Pt;
  std::vector<TH1D> control_Track_Eta;
  std::vector<TH1D> control_Track_Phi;

  std::vector<TH1D> VertexKFChi2_peak;
  std::vector<TH1D> VertexKFChi2_sideband;
  std::vector<TH1D> VertexKFChi2;

  std::vector<TH1D> SVPVDsDirAngle_peak;
  std::vector<TH1D> SVPVDsDirAngle_sideband;
  std::vector<TH1D> SVPVDsDirAngle;

  std::vector<TH1D> NtracksClose_peak;
  std::vector<TH1D> NtracksClose_sideband;
  std::vector<TH1D> NtracksClose;
  
  std::vector<TH1D> NSV_peak;
  std::vector<TH1D> NSV_sideband;
  std::vector<TH1D> NSV;


  std::vector<TH1D>  MinMuonImpacAngle_peak;
  std::vector<TH1D>  MaxMuonImpacAngle_peak;

  std::vector<TH1D>  MinMuonImpacAngle_sideband;
  std::vector<TH1D>  MaxMuonImpacAngle_sideband;

  std::vector<TH1D>  MinMuonImpacAngle;
  std::vector<TH1D>  MaxMuonImpacAngle;


  std::vector<TH1D> MinMuon_chi2LocalPosition_peak;
  std::vector<TH1D> MinMuon_chi2LocalPosition_sideband;
  std::vector<TH1D> MinMuon_chi2LocalPosition;
  
  std::vector<TH1D> Mu1Muon_chi2LocalPosition;
  std::vector<TH1D> Mu2Muon_chi2LocalPosition;

  std::vector<TH1D> Mu1Muon_chi2LocalPosition_peak;
  std::vector<TH1D> Mu2Muon_chi2LocalPosition_peak;

  std::vector<TH1D> Mu1Muon_chi2LocalPosition_sideband;
  std::vector<TH1D> Mu2Muon_chi2LocalPosition_sideband;


  std::vector<TH1D> MindcaTrackSV_peak;
  std::vector<TH1D> MindcaTrackSV_sideband;
  std::vector<TH1D> MindcaTrackSV;

  std::vector<TH1D> MinDca_peak;
  std::vector<TH1D> MinDca_sideband;
  std::vector<TH1D> MinDca;

  std::vector<TH1D> MinD0SigSV_peak;
  std::vector<TH1D> MinD0SigSV_sideband;
  std::vector<TH1D> MinD0SigSV;

  std::vector<TH1D> MinD0SigPV_peak;
  std::vector<TH1D> MinD0SigPV_sideband;
  std::vector<TH1D> MinD0SigPV;

  std::vector<TH1D> MaxVertexPairQuality_peak;
  std::vector<TH1D> MaxVertexPairQuality_sideband;
  std::vector<TH1D> MaxVertexPairQuality;

  std::vector<TH1D> MaxdeltaMuZ_peak;
  std::vector<TH1D> MaxdeltaMuZ_sideband;
  std::vector<TH1D> MaxdeltaMuZ;

  std::vector<TH1D> MaxDca_peak;
  std::vector<TH1D> MaxDca_sideband;
  std::vector<TH1D> MaxDca;

  std::vector<TH1D> MaxD0SigSV_peak;
  std::vector<TH1D> MaxD0SigSV_sideband;
  std::vector<TH1D> MaxD0SigSV;

  std::vector<TH1D> MaxD0SigPV_peak;
  std::vector<TH1D> MaxD0SigPV_sideband;
  std::vector<TH1D> MaxD0SigPV;


  std::vector<TH1D>   MaxMuDistanceToSV_peak;
  std::vector<TH1D>   MinMuDistanceToSV_peak;

  std::vector<TH1D>   MaxMuDistanceToSV_sideband;
  std::vector<TH1D>   MinMuDistanceToSV_sideband;

  std::vector<TH1D>   MaxMuDistanceToSV;
  std::vector<TH1D>   MinMuDistanceToSV;


  std::vector<TH1D>   mu1segmentCompatibility;
  std::vector<TH1D>   mu2segmentCompatibility;

  std::vector<TH1D>   mu1segmentCompatibility_peak;
  std::vector<TH1D>   mu2segmentCompatibility_peak;

  std::vector<TH1D>   mu1segmentCompatibility_sideband;
  std::vector<TH1D>   mu2segmentCompatibility_sideband;

  std::vector<TH1D> Muon1ImpactAngle_peak;
  std::vector<TH1D> Muon2ImpactAngle_peak;

  std::vector<TH1D> Muon1ImpactAngle_sideband;
  std::vector<TH1D> Muon2ImpactAngle_sideband;

  std::vector<TH1D> Muon1ImpactAngle;
  std::vector<TH1D> Muon2ImpactAngle;


  std::vector<TH1D> MaxD0SigBS_peak;
  std::vector<TH1D> MaxD0SigBS_sideband;
  std::vector<TH1D> MaxD0SigBS;

  std::vector<TH1D> Iso1_peak;
  std::vector<TH1D> Iso1_sideband;
  std::vector<TH1D> Iso1;

  std::vector<TH1D> Iso1Mu1_peak;
  std::vector<TH1D> Iso1Mu1_sideband;
  std::vector<TH1D> Iso1Mu1;

  std::vector<TH1D> Iso8Mu1_peak;
  std::vector<TH1D> Iso8Mu1_sideband;
  std::vector<TH1D> Iso8Mu1;

  std::vector<TH1D> FLSignificance_peak;
  std::vector<TH1D> FLSignificance_sideband;
  std::vector<TH1D> FLSignificance;


  std::vector<TH1D>  Muon1MVAID_peak;
  std::vector<TH1D>  Muon2MVAID_peak;

  std::vector<TH1D>  Muon1MVAID_sideband;
  std::vector<TH1D>  Muon2MVAID_sideband;

  std::vector<TH1D>  Muon1MVAID;
  std::vector<TH1D>  Muon2MVAID;


  std::vector<TH1D>  SVPosX_peak;
  std::vector<TH1D>  SVPosX_sideband;
  std::vector<TH1D>  SVPosX;

  std::vector<TH1D>  SVPosY_peak;
  std::vector<TH1D>  SVPosY_sideband;
  std::vector<TH1D>  SVPosY;

  std::vector<TH1D>  SVPosZ_peak;
  std::vector<TH1D>  SVPosZ_sideband;
  std::vector<TH1D>  SVPosZ;


  std::vector<TH1D>  PVPosX_peak;
  std::vector<TH1D>  PVPosX_sideband;
  std::vector<TH1D>  PVPosX;

  std::vector<TH1D>  PVPosY_peak;
  std::vector<TH1D>  PVPosY_sideband;
  std::vector<TH1D>  PVPosY;

  std::vector<TH1D>  PVPosZ_peak;
  std::vector<TH1D>  PVPosZ_sideband;
  std::vector<TH1D>  PVPosZ;


  std::vector<TH1D>  SVPVEta_peak;
  std::vector<TH1D>  SVPVPhi_peak;

  std::vector<TH1D>  SVPVEta_sideband;
  std::vector<TH1D>  SVPVPhi_sideband;

  std::vector<TH1D>  SVPVEta;
  std::vector<TH1D>  SVPVPhi;

  std::vector<TH1D>  SVPVTransverse_peak;
  std::vector<TH1D>  SVPVTransverse_sideband;
  std::vector<TH1D>  SVPVTransverse;





  std::vector<TH1D>  Mu1combinedQuality_chi2LocalMomentum_peak;
  std::vector<TH1D>  Mu1combinedQuality_chi2LocalPosition_peak;
  std::vector<TH1D>  Mu1combinedQuality_staRelChi2_peak;
  std::vector<TH1D>  Mu1combinedQuality_trkRelChi2_peak;
  std::vector<TH1D>  Mu1combinedQuality_globalDeltaEtaPhi_peak;
  std::vector<TH1D>  Mu1mu_combinedQuality_trkKink_peak;
  std::vector<TH1D>  Mu1mu_combinedQuality_glbKink_peak;
  std::vector<TH1D>  Mu1combinedQuality_glbTrackProbability_peak;
  std::vector<TH1D>  Mu1Numberofvalidtrackerhits_peak;
  std::vector<TH1D>  Mu1Numberofvalidpixelhits_peak;
  std::vector<TH1D>  Mu1trackerLayersWithMeasurement_peak;
  std::vector<TH1D>  Mu1validMuonHitComb_peak;
  std::vector<TH1D>  Mu1numberOfMatchedStations_peak;
  std::vector<TH1D>  Mu1segmentCompatibility_peak;
  std::vector<TH1D>  Mu1timeAtIpInOutErr_peak;
  std::vector<TH1D>  Mu1timeAtIpInOut_peak;
  std::vector<TH1D>  Mu1GLnormChi2_peak;
  std::vector<TH1D>  Mu1innerTrack_normalizedChi2_peak;
  std::vector<TH1D>  Mu1outerTrack_normalizedChi2_peak;
  std::vector<TH1D>  Mu1innerTrack_validFraction_peak;

  std::vector<TH1D>  Mu1combinedQuality_chi2LocalMomentum_sideband;
  std::vector<TH1D>  Mu1combinedQuality_chi2LocalPosition_sideband;
  std::vector<TH1D>  Mu1combinedQuality_staRelChi2_sideband;
  std::vector<TH1D>  Mu1combinedQuality_trkRelChi2_sideband;
  std::vector<TH1D>  Mu1combinedQuality_globalDeltaEtaPhi_sideband;
  std::vector<TH1D>  Mu1mu_combinedQuality_trkKink_sideband;
  std::vector<TH1D>  Mu1mu_combinedQuality_glbKink_sideband;
  std::vector<TH1D>  Mu1combinedQuality_glbTrackProbability_sideband;
  std::vector<TH1D>  Mu1Numberofvalidtrackerhits_sideband;
  std::vector<TH1D>  Mu1Numberofvalidpixelhits_sideband;
  std::vector<TH1D>  Mu1trackerLayersWithMeasurement_sideband;
  std::vector<TH1D>  Mu1validMuonHitComb_sideband;
  std::vector<TH1D>  Mu1numberOfMatchedStations_sideband;
  std::vector<TH1D>  Mu1segmentCompatibility_sideband;
  std::vector<TH1D>  Mu1timeAtIpInOutErr_sideband;
  std::vector<TH1D>  Mu1timeAtIpInOut_sideband;
  std::vector<TH1D>  Mu1GLnormChi2_sideband;
  std::vector<TH1D>  Mu1innerTrack_normalizedChi2_sideband;
  std::vector<TH1D>  Mu1outerTrack_normalizedChi2_sideband;
  std::vector<TH1D>  Mu1innerTrack_validFraction_sideband;

  std::vector<TH1D>  Mu1combinedQuality_chi2LocalMomentum;
  std::vector<TH1D>  Mu1combinedQuality_chi2LocalPosition;
  std::vector<TH1D>  Mu1combinedQuality_staRelChi2;
  std::vector<TH1D>  Mu1combinedQuality_trkRelChi2;
  std::vector<TH1D>  Mu1combinedQuality_globalDeltaEtaPhi;
  std::vector<TH1D>  Mu1mu_combinedQuality_trkKink;
  std::vector<TH1D>  Mu1mu_combinedQuality_glbKink;
  std::vector<TH1D>  Mu1combinedQuality_glbTrackProbability;
  std::vector<TH1D>  Mu1Numberofvalidtrackerhits;
  std::vector<TH1D>  Mu1Numberofvalidpixelhits;
  std::vector<TH1D>  Mu1trackerLayersWithMeasurement;
  std::vector<TH1D>  Mu1validMuonHitComb;
  std::vector<TH1D>  Mu1numberOfMatchedStations;
  std::vector<TH1D>  Mu1segmentCompatibility;
  std::vector<TH1D>  Mu1timeAtIpInOutErr;
  std::vector<TH1D>  Mu1timeAtIpInOut;
  std::vector<TH1D>  Mu1GLnormChi2;
  std::vector<TH1D>  Mu1innerTrack_normalizedChi2;
  std::vector<TH1D>  Mu1outerTrack_normalizedChi2;
  std::vector<TH1D>  Mu1innerTrack_validFraction;






  std::vector<TH1D>  Mu2combinedQuality_chi2LocalMomentum_peak;
  std::vector<TH1D>  Mu2combinedQuality_chi2LocalPosition_peak;
  std::vector<TH1D>  Mu2combinedQuality_staRelChi2_peak;
  std::vector<TH1D>  Mu2combinedQuality_trkRelChi2_peak;
  std::vector<TH1D>  Mu2combinedQuality_globalDeltaEtaPhi_peak;
  std::vector<TH1D>  Mu2mu_combinedQuality_trkKink_peak;
  std::vector<TH1D>  Mu2mu_combinedQuality_glbKink_peak;
  std::vector<TH1D>  Mu2combinedQuality_glbTrackProbability_peak;
  std::vector<TH1D>  Mu2Numberofvalidtrackerhits_peak;
  std::vector<TH1D>  Mu2Numberofvalidpixelhits_peak;
  std::vector<TH1D>  Mu2trackerLayersWithMeasurement_peak;
  std::vector<TH1D>  Mu2validMuonHitComb_peak;
  std::vector<TH1D>  Mu2numberOfMatchedStations_peak;
  std::vector<TH1D>  Mu2segmentCompatibility_peak;
  std::vector<TH1D>  Mu2timeAtIpInOutErr_peak;
  std::vector<TH1D>  Mu2timeAtIpInOut_peak;
  std::vector<TH1D>  Mu2GLnormChi2_peak;
  std::vector<TH1D>  Mu2innerTrack_normalizedChi2_peak;
  std::vector<TH1D>  Mu2outerTrack_normalizedChi2_peak;
  std::vector<TH1D>  Mu2innerTrack_validFraction_peak;

  std::vector<TH1D>  Mu2combinedQuality_chi2LocalMomentum_sideband;
  std::vector<TH1D>  Mu2combinedQuality_chi2LocalPosition_sideband;
  std::vector<TH1D>  Mu2combinedQuality_staRelChi2_sideband;
  std::vector<TH1D>  Mu2combinedQuality_trkRelChi2_sideband;
  std::vector<TH1D>  Mu2combinedQuality_globalDeltaEtaPhi_sideband;
  std::vector<TH1D>  Mu2mu_combinedQuality_trkKink_sideband;
  std::vector<TH1D>  Mu2mu_combinedQuality_glbKink_sideband;
  std::vector<TH1D>  Mu2combinedQuality_glbTrackProbability_sideband;
  std::vector<TH1D>  Mu2Numberofvalidtrackerhits_sideband;
  std::vector<TH1D>  Mu2trackerLayersWithMeasurement_sideband;
  std::vector<TH1D>  Mu2Numberofvalidpixelhits_sideband;
  std::vector<TH1D>  Mu2validMuonHitComb_sideband;
  std::vector<TH1D>  Mu2numberOfMatchedStations_sideband;
  std::vector<TH1D>  Mu2segmentCompatibility_sideband;
  std::vector<TH1D>  Mu2timeAtIpInOutErr_sideband;
  std::vector<TH1D>  Mu2timeAtIpInOut_sideband;
  std::vector<TH1D>  Mu2GLnormChi2_sideband;
  std::vector<TH1D>  Mu2innerTrack_normalizedChi2_sideband;
  std::vector<TH1D>  Mu2outerTrack_normalizedChi2_sideband;
  std::vector<TH1D>  Mu2innerTrack_validFraction_sideband;

  std::vector<TH1D>  Mu2combinedQuality_chi2LocalMomentum;
  std::vector<TH1D>  Mu2combinedQuality_chi2LocalPosition;
  std::vector<TH1D>  Mu2combinedQuality_staRelChi2;
  std::vector<TH1D>  Mu2combinedQuality_trkRelChi2;
  std::vector<TH1D>  Mu2combinedQuality_globalDeltaEtaPhi;
  std::vector<TH1D>  Mu2mu_combinedQuality_trkKink;
  std::vector<TH1D>  Mu2mu_combinedQuality_glbKink;
  std::vector<TH1D>  Mu2combinedQuality_glbTrackProbability;
  std::vector<TH1D>  Mu2trackerLayersWithMeasurement;
  std::vector<TH1D>  Mu2Numberofvalidtrackerhits;
  std::vector<TH1D>  Mu2Numberofvalidpixelhits;
  std::vector<TH1D>  Mu2validMuonHitComb;
  std::vector<TH1D>  Mu2numberOfMatchedStations;
  std::vector<TH1D>  Mu2segmentCompatibility;
  std::vector<TH1D>  Mu2timeAtIpInOutErr;
  std::vector<TH1D>  Mu2timeAtIpInOut;
  std::vector<TH1D>  Mu2GLnormChi2;
  std::vector<TH1D>  Mu2innerTrack_normalizedChi2;
  std::vector<TH1D>  Mu2outerTrack_normalizedChi2;
  std::vector<TH1D>  Mu2innerTrack_validFraction;

  std::vector<TH1D>  MatchedMuonTrackdPT;


  std::vector<TH1D> TrkIDTrackmuonInnerNC2;
  std::vector<TH1D> TrkIDTrackmuonInnerNValidHits;
  std::vector<TH1D> TrkIDTrackmuonValidFraction;
  std::vector<TH1D> TrkIDTrackmuonNLostTrackerHits;
  std::vector<TH1D> TrkIDTrackmuonNLostTrackerHitsInner;
  std::vector<TH1D> TrkIDTrackmuonNLostTrackerHitsOuter;
  std::vector<TH1D> TrkIDTrackmuonPixelLayers;
  std::vector<TH1D> TrkIDTrackmuonNMatchedStations;
  std::vector<TH1D> TrkIDTrackmuonPtErrPt;
  std::vector<TH1D> TrkIDTrackmuonCaloComp;
  std::vector<TH1D> TrkIDTrackmuonSegComp;
  std::vector<TH1D> TrkIDTrackmuonHad;
  std::vector<TH1D> TrkIDTrackmuonEM;


  std::vector<TH1D> TrkIDTrackmuonInnerNC2_sideband;
  std::vector<TH1D> TrkIDTrackmuonInnerNValidHits_sideband;
  std::vector<TH1D> TrkIDTrackmuonValidFraction_sideband;
  std::vector<TH1D> TrkIDTrackmuonNLostTrackerHits_sideband;
  std::vector<TH1D> TrkIDTrackmuonNLostTrackerHitsInner_sideband;
  std::vector<TH1D> TrkIDTrackmuonNLostTrackerHitsOuter_sideband;
  std::vector<TH1D> TrkIDTrackmuonPixelLayers_sideband;
  std::vector<TH1D> TrkIDTrackmuonNMatchedStations_sideband;
  std::vector<TH1D> TrkIDTrackmuonPtErrPt_sideband;
  std::vector<TH1D> TrkIDTrackmuonCaloComp_sideband;
  std::vector<TH1D> TrkIDTrackmuonSegComp_sideband;
  std::vector<TH1D> TrkIDTrackmuonHad_sideband;
  std::vector<TH1D> TrkIDTrackmuonEM_sideband;


  std::vector<TH1D> TrkIDTrackmuonInnerNC2_peak;
  std::vector<TH1D> TrkIDTrackmuonInnerNValidHits_peak;
  std::vector<TH1D> TrkIDTrackmuonValidFraction_peak;
  std::vector<TH1D> TrkIDTrackmuonNLostTrackerHits_peak;
  std::vector<TH1D> TrkIDTrackmuonNLostTrackerHitsInner_peak;
  std::vector<TH1D> TrkIDTrackmuonNLostTrackerHitsOuter_peak;
  std::vector<TH1D> TrkIDTrackmuonPixelLayers_peak;
  std::vector<TH1D> TrkIDTrackmuonNMatchedStations_peak;
  std::vector<TH1D> TrkIDTrackmuonPtErrPt_peak;
  std::vector<TH1D> TrkIDTrackmuonCaloComp_peak;
  std::vector<TH1D> TrkIDTrackmuonSegComp_peak;
  std::vector<TH1D> TrkIDTrackmuonHad_peak;
  std::vector<TH1D> TrkIDTrackmuonEM_peak;

  std::vector<TH1D> TrkIDTrackBDT_peak;
  std::vector<TH1D> TrkIDTrackBDT_sideband;
  std::vector<TH1D> TrkIDTrackBDT;



  //--- muonID variables


  Float_t mu_combinedQuality_chi2LocalMomentum;
  Float_t mu_combinedQuality_chi2LocalPosition;
  Float_t mu_combinedQuality_staRelChi2;
  Float_t mu_combinedQuality_trkRelChi2;
  Float_t mu_combinedQuality_globalDeltaEtaPhi;
  Float_t mu_combinedQuality_trkKink;
  Float_t mu_combinedQuality_glbKink;
  Float_t mu_combinedQuality_glbTrackProbability;
  Float_t mu_Numberofvalidtrackerhits;
  Float_t mu_trackerLayersWithMeasurement;
  Float_t mu_Numberofvalidpixelhits;
  Float_t mu_validMuonHitComb;
  Float_t mu_numberOfMatchedStations;
  Float_t mu_segmentCompatibility;
  Float_t mu_timeAtIpInOutErr;
  Float_t mu_GLnormChi2;
  Float_t mu_innerTrack_normalizedChi2;
  Float_t mu_outerTrack_normalizedChi2;
  Float_t mu_innerTrack_validFraction;

  Float_t mu_eta;
  Float_t mu_pt;
  Float_t mu_phi;
  Float_t mu_SoftMVA;


  Float_t Muon1DetID;
  Float_t Muon2DetID;
  Float_t Muon3DetID;


  //   tracker vs mu id

  float muonPt ;
  float muonEta ;
  float muonPhi ;

  float fake;
  float muonInnerNC2 ;
  float muonValidFraction;
  float muonInnerNValidHits ;
  float muonNLostTrackerHits ;
  float muonNLostTrackerHitsInner ;
  float muonNLostTrackerHitsOuter ;
  float muonPixelLayers ;
  float muonNMatchedStations ;
  float muonPtErrPt ;
  float muonSegComp ;
  float muonCaloComp ;
  float muonHad ;
  float muonEM ; 
  float TrackMuBDT_Track;
  float TrackMuBDT_Muon;

};
#endif
