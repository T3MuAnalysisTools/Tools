#ifndef T3MSelectionTree_h
#define T3MSelectionTree_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "TRandom.h"
#include "TLorentzVector.h"

class T3MSelectionTree : public Selection {

   public:
      T3MSelectionTree(TString Name_, TString id_);
      virtual ~T3MSelectionTree();

      virtual void  Configure();
      virtual void  Finish();

      double GetReferenceFrameAngle(TLorentzVector, TLorentzVector, TLorentzVector);
      double GetReferenceFrameCosine(TLorentzVector, TLorentzVector, TLorentzVector);
      
      enum cuts {SignalCandidate=0, L1Fired, HLTFired, KFChi2, BSSVSig, PFMuons, MuonID, Mu1PtCut, Mu2PtCut, Mu3PtCut, TauMassCut, TriggerMatch, PhiVetoOS1, OmegaVetoOS1, PhiVetoOS2,  OmegaVetoOS2, TrackerLayers, NCuts};

   protected:
      virtual void doEvent();  
      virtual void Store_ExtraDist();

      template <typename T>
         int minQuantityIndex(std::vector<T>& vec);

      template <typename T>
         int maxQuantityIndex(std::vector<T>& vec);



   private:

      // random number generator
      TRandom rndm;
      float random_num;
      int l1FailedRandom;
      ULong64_t eventNumber;
      int run;
      int lumi;
      // PU Weights
      TFile* PUWeightFile;
      TH1D* puWeights;
      
      TFile * file;
      TTree * TMVA_Tree;

      // Selection Variables
      double tauMinMass_, tauMaxMass_;
      double tauMinSideBand_,tauMaxSideBand_;
      double tauMassResCutLow, tauMassResCutHigh;
      double phiVetoSigma, omegaVetoSigma;
      double M_osss1, M_osss2;

      unsigned int Muon_index_1, Muon_index_2, Muon_index_3;

      // categorization variables
      bool MC;
      float category;
      int DataMCType;
      float EventWeight;
      bool threeGlobal;
      int l1seed;

      //commmon variables (2016 + 2017)
      float var_vertexKFChi2; //  KF chi2
      float var_svpvTauAngle; 
      float var_flightLenSig;
      float var_nMatches_mu3;
      float var_max_cLP;
      float var_max_tKink;
      float var_segCompMuMin;
      float var_MinMIPLikelihood;
      float var_sumMuTrkKinkChi2;
      float var_maxdca;
      float var_MuMu_minKFChi2;
      float var_MinD0Significance;
      float var_MaxD0Significance;
      float var_trk_relPt;
      float var_pmin;
      float var_mindca_iso;
      float var_minMatchedStations;
      
      // 2017 variables
      float var_Eta_Tau; 
      float var_tauMassRes;
      float var_tauMass;
      float var_tauMassRefit;

      float var_NtracksClose;
      float var_MindcaTrackSV;
      float var_Mu1TrackMass;
      float var_Mu2TrackMass;
      float var_Mu3TrackMass;
      float var_dcaTrackPV;

      // Muon variables
      float var_Muon1_Pt;
      float var_Muon2_Pt;
      float var_Muon3_Pt;

      float var_Muon1_Eta;
      float var_Muon2_Eta;
      float var_Muon3_Eta;

      float var_Muon1_Phi;
      float var_Muon2_Phi;
      float var_Muon3_Phi;

      float var_Tau_Pt;
      float var_Tau_Eta;
      float var_Tau_Phi;

      float Muon1_station_vars[8][7];
      float Muon2_station_vars[8][7];
      float Muon3_station_vars[8][7];

      float var_Muon1_vx;
      float var_Muon1_vy;
      float var_Muon1_vz;
      bool var_Muon1_IsGlobalMuon;
      bool var_Muon1_IsStandAloneMuon;
      bool var_Muon1_IsTrackerMuon;
      bool var_Muon1_IsCaloMuon;
      bool var_Muon1_IsIsolationValid;
      bool var_Muon1_IsQualityValid;
      bool var_Muon1_IsTimeValid;
      bool var_Muon1_IsPFMuon;
      bool var_Muon1_IsRPCMuon;
      float var_Muon1_emEt03;
      float var_Muon1_emVetoEt03;
      float var_Muon1_hadEt03;
      float var_Muon1_hadVetoEt03;
      float var_Muon1_nJets03;
      float var_Muon1_nTracks03;
      float var_Muon1_StandardSelection;
      float var_Muon1_sumPt03;
      float var_Muon1_trackerVetoPt03;
      float var_Muon1_sumChargedHadronPt03;
      float var_Muon1_sumChargedParticlePt03;
      float var_Muon1_sumNeutralHadronEt03;
      float var_Muon1_sumNeutralHadronEtHighThreshold03;
      float var_Muon1_sumPhotonEt03;
      float var_Muon1_sumPhotonEtHighThreshold03;
      float var_Muon1_sumPUPt03;
      float var_Muon1_numberOfChambers;
      float var_Muon1_Track_idx;
      float var_Muon1_combinedQuality_updatedSta;
      float var_Muon1_combinedQuality_trkKink;
      float var_Muon1_combinedQuality_glbKink;
      float var_Muon1_combinedQuality_trkRelChi2;
      float var_Muon1_combinedQuality_staRelChi2;
      float var_Muon1_combinedQuality_chi2LocalPosition;
      float var_Muon1_combinedQuality_chi2LocalMomentum;
      float var_Muon1_combinedQuality_localDistance;
      float var_Muon1_combinedQuality_globalDeltaEtaPhi;
      float var_Muon1_combinedQuality_tightMatch;
      float var_Muon1_combinedQuality_glbTrackProbability;
      float var_Muon1_prod_inner_outer_charge;
      float var_Muon1_innerTrack_quality;
      float var_Muon1_ptErrOverPt;
      float var_Muon1_calEnergy_hadS9;
      float var_Muon1_calEnergy_had;
      float var_Muon1_calEnergy_emS25;
      float var_Muon1_calEnergy_emS9;
      float var_Muon1_calEnergy_em;
      int var_Muon1_charge;
      float var_Muon1_trackCharge;
      float var_Muon1_hitPattern_pixelLayerwithMeas;
      float var_Muon1_numberOfMatchedStations;
      float var_Muon1_normChi2;
      float var_Muon1_hitPattern_numberOfValidMuonHits;
      float var_Muon1_innerTrack_numberofValidHits;
      float var_Muon1_numberofValidPixelHits;
      float var_Muon1_numberOfMatches;
      float var_Muon1_trackerLayersWithMeasurement;
      float var_Muon1_segmentCompatibility;
      float var_Muon1_caloCompatibility;
      float var_Muon1_innerTrack_validFraction;
      float var_Muon1_innerTrack_pixelLayersWithMeasurement;
      float var_Muon1_innerTrack_numberOfValidTrackerHits;
      float var_Muon1_innerTrack_numberOfLostTrackerHits;
      float var_Muon1_innerTrack_numberOfLostTrackerInnerHits;
      float var_Muon1_innerTrack_numberOfLostTrackerOuterHits;
      float var_Muon1_innerTrack_normalizedChi2;
      float var_Muon1_vmuonhitcomb_reco;
      float var_Muon1_rpchits_reco;
      float var_Muon1_outerTrack_normalizedChi2;
      float var_Muon1_outerTrack_muonStationsWithValidHits;
      bool var_Muon1_isGoodMuon_TM2DCompatibility;
      bool var_Muon1_isGoodMuon_TrackerMuonArbitrated;
      bool var_Muon1_isGoodMuon_TMOneStationTight;
      bool var_Muon1_isGoodMuon_TMOneStationAngTight;
      bool var_Muon1_isGoodMuon_TMLastStationTight;
      bool var_Muon1_isGoodMuon_TMLastStationAngTight;
      bool var_Muon1_isGoodMuon_TMLastStationOptimizedLowPtTight;
      bool var_Muon1_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight;

      float var_Muon2_vx;
      float var_Muon2_vy;
      float var_Muon2_vz;
      bool var_Muon2_IsGlobalMuon;
      bool var_Muon2_IsStandAloneMuon;
      bool var_Muon2_IsTrackerMuon;
      bool var_Muon2_IsCaloMuon;
      bool var_Muon2_IsIsolationValid;
      bool var_Muon2_IsQualityValid;
      bool var_Muon2_IsTimeValid;
      bool var_Muon2_IsPFMuon;
      bool var_Muon2_IsRPCMuon;
      float var_Muon2_emEt03;
      float var_Muon2_emVetoEt03;
      float var_Muon2_hadEt03;
      float var_Muon2_hadVetoEt03;
      float var_Muon2_nJets03;
      float var_Muon2_nTracks03;
      float var_Muon2_StandardSelection;
      float var_Muon2_sumPt03;
      float var_Muon2_trackerVetoPt03;
      float var_Muon2_sumChargedHadronPt03;
      float var_Muon2_sumChargedParticlePt03;
      float var_Muon2_sumNeutralHadronEt03;
      float var_Muon2_sumNeutralHadronEtHighThreshold03;
      float var_Muon2_sumPhotonEt03;
      float var_Muon2_sumPhotonEtHighThreshold03;
      float var_Muon2_sumPUPt03;
      float var_Muon2_numberOfChambers;
      float var_Muon2_Track_idx;
      float var_Muon2_combinedQuality_updatedSta;
      float var_Muon2_combinedQuality_trkKink;
      float var_Muon2_combinedQuality_glbKink;
      float var_Muon2_combinedQuality_trkRelChi2;
      float var_Muon2_combinedQuality_staRelChi2;
      float var_Muon2_combinedQuality_chi2LocalPosition;
      float var_Muon2_combinedQuality_chi2LocalMomentum;
      float var_Muon2_combinedQuality_localDistance;
      float var_Muon2_combinedQuality_globalDeltaEtaPhi;
      float var_Muon2_combinedQuality_tightMatch;
      float var_Muon2_combinedQuality_glbTrackProbability;
      float var_Muon2_prod_inner_outer_charge;
      float var_Muon2_innerTrack_quality;
      float var_Muon2_ptErrOverPt;
      float var_Muon2_calEnergy_hadS9;
      float var_Muon2_calEnergy_had;
      float var_Muon2_calEnergy_emS25;
      float var_Muon2_calEnergy_emS9;
      float var_Muon2_calEnergy_em;
      int var_Muon2_charge;
      float var_Muon2_trackCharge;
      float var_Muon2_hitPattern_pixelLayerwithMeas;
      float var_Muon2_numberOfMatchedStations;
      float var_Muon2_normChi2;
      float var_Muon2_hitPattern_numberOfValidMuonHits;
      float var_Muon2_innerTrack_numberofValidHits;
      float var_Muon2_numberofValidPixelHits;
      float var_Muon2_numberOfMatches;
      float var_Muon2_trackerLayersWithMeasurement;
      float var_Muon2_segmentCompatibility;
      float var_Muon2_caloCompatibility;
      float var_Muon2_innerTrack_validFraction;
      float var_Muon2_innerTrack_pixelLayersWithMeasurement;
      float var_Muon2_innerTrack_numberOfValidTrackerHits;
      float var_Muon2_innerTrack_numberOfLostTrackerHits;
      float var_Muon2_innerTrack_numberOfLostTrackerInnerHits;
      float var_Muon2_innerTrack_numberOfLostTrackerOuterHits;
      float var_Muon2_innerTrack_normalizedChi2;
      float var_Muon2_vmuonhitcomb_reco;
      float var_Muon2_rpchits_reco;
      float var_Muon2_outerTrack_normalizedChi2;
      float var_Muon2_outerTrack_muonStationsWithValidHits;
      bool var_Muon2_isGoodMuon_TM2DCompatibility;
      bool var_Muon2_isGoodMuon_TrackerMuonArbitrated;
      bool var_Muon2_isGoodMuon_TMOneStationTight;
      bool var_Muon2_isGoodMuon_TMOneStationAngTight;
      bool var_Muon2_isGoodMuon_TMLastStationTight;
      bool var_Muon2_isGoodMuon_TMLastStationAngTight;
      bool var_Muon2_isGoodMuon_TMLastStationOptimizedLowPtTight;
      bool var_Muon2_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight;

      float var_Muon3_vx;
      float var_Muon3_vy;
      float var_Muon3_vz;
      bool var_Muon3_IsGlobalMuon;
      bool var_Muon3_IsStandAloneMuon;
      bool var_Muon3_IsTrackerMuon;
      bool var_Muon3_IsCaloMuon;
      bool var_Muon3_IsIsolationValid;
      bool var_Muon3_IsQualityValid;
      bool var_Muon3_IsTimeValid;
      bool var_Muon3_IsPFMuon;
      bool var_Muon3_IsRPCMuon;
      float var_Muon3_emEt03;
      float var_Muon3_emVetoEt03;
      float var_Muon3_hadEt03;
      float var_Muon3_hadVetoEt03;
      float var_Muon3_nJets03;
      float var_Muon3_nTracks03;
      float var_Muon3_StandardSelection;
      float var_Muon3_sumPt03;
      float var_Muon3_trackerVetoPt03;
      float var_Muon3_sumChargedHadronPt03;
      float var_Muon3_sumChargedParticlePt03;
      float var_Muon3_sumNeutralHadronEt03;
      float var_Muon3_sumNeutralHadronEtHighThreshold03;
      float var_Muon3_sumPhotonEt03;
      float var_Muon3_sumPhotonEtHighThreshold03;
      float var_Muon3_sumPUPt03;
      float var_Muon3_numberOfChambers;
      float var_Muon3_Track_idx;
      float var_Muon3_combinedQuality_updatedSta;
      float var_Muon3_combinedQuality_trkKink;
      float var_Muon3_combinedQuality_glbKink;
      float var_Muon3_combinedQuality_trkRelChi2;
      float var_Muon3_combinedQuality_staRelChi2;
      float var_Muon3_combinedQuality_chi2LocalPosition;
      float var_Muon3_combinedQuality_chi2LocalMomentum;
      float var_Muon3_combinedQuality_localDistance;
      float var_Muon3_combinedQuality_globalDeltaEtaPhi;
      float var_Muon3_combinedQuality_tightMatch;
      float var_Muon3_combinedQuality_glbTrackProbability;
      float var_Muon3_prod_inner_outer_charge;
      float var_Muon3_innerTrack_quality;
      float var_Muon3_ptErrOverPt;
      float var_Muon3_calEnergy_hadS9;
      float var_Muon3_calEnergy_had;
      float var_Muon3_calEnergy_emS25;
      float var_Muon3_calEnergy_emS9;
      float var_Muon3_calEnergy_em;
      int var_Muon3_charge;
      float var_Muon3_trackCharge;
      float var_Muon3_hitPattern_pixelLayerwithMeas;
      float var_Muon3_numberOfMatchedStations;
      float var_Muon3_normChi2;
      float var_Muon3_hitPattern_numberOfValidMuonHits;
      float var_Muon3_innerTrack_numberofValidHits;
      float var_Muon3_numberofValidPixelHits;
      float var_Muon3_numberOfMatches;
      float var_Muon3_trackerLayersWithMeasurement;
      float var_Muon3_segmentCompatibility;
      float var_Muon3_caloCompatibility;
      float var_Muon3_innerTrack_validFraction;
      float var_Muon3_innerTrack_pixelLayersWithMeasurement;
      float var_Muon3_innerTrack_numberOfValidTrackerHits;
      float var_Muon3_innerTrack_numberOfLostTrackerHits;
      float var_Muon3_innerTrack_numberOfLostTrackerInnerHits;
      float var_Muon3_innerTrack_numberOfLostTrackerOuterHits;
      float var_Muon3_innerTrack_normalizedChi2;
      float var_Muon3_vmuonhitcomb_reco;
      float var_Muon3_rpchits_reco;
      float var_Muon3_outerTrack_normalizedChi2;
      float var_Muon3_outerTrack_muonStationsWithValidHits;
      bool var_Muon3_isGoodMuon_TM2DCompatibility;
      bool var_Muon3_isGoodMuon_TrackerMuonArbitrated;
      bool var_Muon3_isGoodMuon_TMOneStationTight;
      bool var_Muon3_isGoodMuon_TMOneStationAngTight;
      bool var_Muon3_isGoodMuon_TMLastStationTight;
      bool var_Muon3_isGoodMuon_TMLastStationAngTight;
      bool var_Muon3_isGoodMuon_TMLastStationOptimizedLowPtTight;
      bool var_Muon3_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight;

      // muon time
      float var_Muon1_timeAtIpInOutErr;
      float var_Muon1_timeAtIpInOut;
      float var_Muon2_timeAtIpInOutErr;
      float var_Muon2_timeAtIpInOut;
      float var_Muon3_timeAtIpInOutErr;
      float var_Muon3_timeAtIpInOut;

      // Impact angle
      float var_Muon1ImpactAngle;
      float var_Muon2ImpactAngle;
      float var_Muon3ImpactAngle;

      // Dimuon variables
      float var_NMuTrkOS1;
      float var_NMuTrkOS2;
      float var_NMuTrkOS3;

      float var_MuTrkInvMassOS1Pair1;
      float var_MuKTrkInvMassOS1Pair1;
      float var_MuTrkOS1Pair1dR;
      float var_MuTrkOS1Pair1_CosThetaK;

      float var_MuTrkInvMassOS2Pair1;
      float var_MuKTrkInvMassOS2Pair1;
      float var_MuTrkOS2Pair1dR;
      float var_MuTrkOS2Pair1_CosThetaK;

      float var_MuTrkInvMassOS1Pair2;
      float var_MuKTrkInvMassOS1Pair2;
      float var_MuTrkOS1Pair2dR;
      float var_MuTrkOS1Pair2_CosThetaK;

      float var_MuTrkInvMassOS2Pair2;
      float var_MuKTrkInvMassOS2Pair2;
      float var_MuTrkOS2Pair2dR;
      float var_MuTrkOS2Pair2_CosThetaK;

      float var_MuTrkInvMassOS1Pair3;
      float var_MuKTrkInvMassOS1Pair3;
      float var_MuTrkOS1Pair3dR;
      float var_MuTrkOS1Pair3_CosThetaK;

      float var_MuTrkInvMassOS2Pair3;
      float var_MuKTrkInvMassOS2Pair3;
      float var_MuTrkOS2Pair3dR;
      float var_MuTrkOS2Pair3_CosThetaK;

      // Dimuon variables
      float var_Mosss1;
      float var_Mosss2;
      float var_ntracks;
      float var_nsv;

      float var_VertexMu1D0SigPVReco;
      float var_VertexMu2D0SigPVReco;
      float var_VertexMu3D0SigPVReco;

      float var_VertexMu1D0SigBSReco;
      float var_VertexMu2D0SigBSReco;
      float var_VertexMu3D0SigBSReco;

      float var_VertexMu1D0SigSVReco;
      float var_VertexMu2D0SigSVReco;
      float var_VertexMu3D0SigSVReco;

      // Isolation
      float var_InvMassConeIso0p3;
      float var_InvMassConeIso0p5;
      float var_InvMassConeIso0p8;
      float var_InvMassConeIso1p2;
      float var_InvMassConeIso1p4;
      float var_NTrksIso0p3;
      float var_NTrksIso0p5;
      float var_NTrksIso0p8;
      float var_NTrksIso1p2;
      float var_NTrksIso1p4;

      float var_Iso02;
      float var_Iso04;
      float var_Iso06;
      float var_Iso08;
      float var_Iso1;
      float var_Iso12;
      float var_Iso14;
      float var_Iso16;
      float var_Iso18;
      float var_Iso2;

      float var_Iso02Mu1;
      float var_Iso04Mu1;
      float var_Iso06Mu1;
      float var_Iso08Mu1;
      float var_Iso1Mu1;
      float var_Iso12Mu1;
      float var_Iso14Mu1;
      float var_Iso16Mu1;
      float var_Iso18Mu1;
      float var_Iso2Mu1;

      float var_Iso02Mu2;
      float var_Iso04Mu2;
      float var_Iso06Mu2;
      float var_Iso08Mu2;
      float var_Iso1Mu2;
      float var_Iso12Mu2;
      float var_Iso14Mu2;
      float var_Iso16Mu2;
      float var_Iso18Mu2;
      float var_Iso2Mu2;

      float var_Iso02Mu3;
      float var_Iso04Mu3;
      float var_Iso06Mu3;
      float var_Iso08Mu3;
      float var_Iso1Mu3;
      float var_Iso12Mu3;
      float var_Iso14Mu3;
      float var_Iso16Mu3;
      float var_Iso18Mu3;
      float var_Iso2Mu3;

      float BTagCSVSH;
      float BTagMVASH;
      float BTagCVSBSH;

      float BTagCSVOH;
      float BTagMVAOH;
      float BTagCVSBOH;

      float BTagCSVSHMatchedToTau;
      float BTagMVASHMatchedToTau;
      float BTagCVSBSHMatchedToTau;

      float BTagCSVOHMatchedToTau;
      float BTagMVAOHMatchedToTau;
      float BTagCVSBOHMatchedToTau;

      int Muon1StandardSelectorPass;
      int Muon2StandardSelectorPass;
      int Muon3StandardSelectorPass;

      int Muon1StandardSelectorFail;
      int Muon2StandardSelectorFail;
      int Muon3StandardSelectorFail;

      float EMR_tau_eta;
      float EventMassResolution_PtEtaPhi;

      float var_Muon1InpactAngle;
      float var_Muon2InpactAngle;
      float var_Muon3InpactAngle;

      float deltaMuZ12;
      float deltaMuZ13;
      float deltaMuZ23;

      float sumTracksIso02;
      float sumTracksIso04; 
      float sumTracksIso06;
      float sumTracksIso08;
      float sumTracksIso1;   
      float sumTracksIso12;
      float sumTracksIso14;
      float sumTracksIso16;
      float sumTracksIso18;
      float sumTracksIso2;
      
      float Isolation_NTracks;
      float Isolation_RelPt;
      float Isolation_maxdxy;
      float VertexMu1D0SigPVReco;
      float VertexMu2D0SigPVReco;
      float VertexMu3D0SigPVReco;
      float VertexMu1D0SigBSReco;
      float VertexMu2D0SigBSReco;
      float VertexMu3D0SigBSReco;
      float VertexMu1D0SigSVReco;
      float VertexMu2D0SigSVReco;
      float VertexMu3D0SigSVReco;
      
      float SV_Mass;
      float SVDeltaR;
      float SVDistance;
      int NSV;

      float var_relPt;
      float var_isoMax;

      float mu3_trk_InvMass;
            
};
#endif
