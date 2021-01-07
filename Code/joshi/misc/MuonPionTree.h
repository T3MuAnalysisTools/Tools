#ifndef MuonPionTree_h
#define MuonPionTree_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class MuonPionTree : public Selection {

   public:
      MuonPionTree(TString Name_, TString id_);
      virtual ~MuonPionTree();

      virtual void  Configure();
      virtual void  Finish();

      /////////////////////////////////////////////////////////
      // This is a cut enumerator, for other event cuts please
      // fill the enumerator with put a new enumarator with an 
      // understandanle name, for exmaple  enum cuts {TriggerOk=0,
      // PrimeVts, EventCut1, EventCut2, ..., NCuts};  
      // Do not remove/rename  the last enumerator   NCuts;

      enum cuts {Candidate, NCuts}; 

   protected:
      virtual void doEvent();  
      virtual void Store_ExtraDist();

   private:

      TFile *file;
      TTree* TMVA_Tree;

      // Selection Variables

      int fake;
      int Muon_pdgId;
      int Muon_motherPdgId;
      float Muon_Pt;
      float Muon_Eta;
      float Muon_Phi;
      float Muon_vx;
      float Muon_vy;
      float Muon_vz;
      int Muon_IsGlobalMuon;
      int Muon_IsStandAloneMuon;
      int Muon_IsTrackerMuon;
      int Muon_IsCaloMuon;
      int Muon_IsIsolationValid;
      int Muon_IsQualityValid;
      int Muon_IsTimeValid;
      int Muon_IsPFMuon;
      int Muon_IsRPCMuon;
      float Muon_emEt03;
      float Muon_emVetoEt03;
      float Muon_hadEt03;
      float Muon_hadVetoEt03;
      int Muon_nJets03;
      int Muon_nTracks03;
      float Muon_StandardSelection;
      float Muon_sumPt03;
      float Muon_trackerVetoPt03;
      float Muon_sumChargedHadronPt03;
      float Muon_sumChargedParticlePt03;
      float Muon_sumNeutralHadronEt03;
      float Muon_sumNeutralHadronEtHighThreshold03;
      float Muon_sumPhotonEt03;
      float Muon_sumPhotonEtHighThreshold03;
      float Muon_sumPUPt03;
      int Muon_numberOfChambers;
      int Muon_Track_idx;
      float Muon_combinedQuality_updatedSta;
      float Muon_combinedQuality_trkKink;
      float Muon_combinedQuality_glbKink;
      float Muon_combinedQuality_trkRelChi2;
      float Muon_combinedQuality_staRelChi2;
      float Muon_combinedQuality_chi2LocalPosition;
      float Muon_combinedQuality_chi2LocalMomentum;
      float Muon_combinedQuality_localDistance;
      float Muon_combinedQuality_globalDeltaEtaPhi;
      float Muon_combinedQuality_tightMatch;
      float Muon_combinedQuality_glbTrackProbability;
      float Muon_prod_inner_outer_charge;
      float Muon_innerTrack_quality;
      float Muon_ptErrOverPt;
      float Muon_calEnergy_hadS9;
      float Muon_calEnergy_had;
      float Muon_calEnergy_emS25;
      float Muon_calEnergy_emS9;
      float Muon_calEnergy_em;
      int Muon_charge;
      int Muon_trackCharge;
      int Muon_hitPattern_pixelLayerwithMeas;
      int Muon_numberOfMatchedStations;
      float Muon_normChi2;
      int Muon_hitPattern_numberOfValidMuonHits;
      int Muon_innerTrack_numberofValidHits;
      float Muon_numberofValidPixelHits;
      int Muon_numberOfMatches;
      int Muon_trackerLayersWithMeasurement;
      float Muon_segmentCompatibility;
      float Muon_caloCompatibility;
      float Muon_innerTrack_validFraction;
      int Muon_innerTrack_pixelLayersWithMeasurement;
      int Muon_innerTrack_numberOfValidTrackerHits;
      int Muon_innerTrack_numberOfLostTrackerHits;
      int Muon_innerTrack_numberOfLostTrackerInnerHits;
      int Muon_innerTrack_numberOfLostTrackerOuterHits;
      float Muon_innerTrack_normalizedChi2;
      float Muon_vmuonhitcomb_reco;
      float Muon_rpchits_reco;
      float Muon_outerTrack_normalizedChi2;
      int Muon_outerTrack_muonStationsWithValidHits;
      float Muon_isGoodMuon_TM2DCompatibility;
      float Muon_isGoodMuon_TrackerMuonArbitrated;
      float Muon_isGoodMuon_TMOneStationTight;
      float Muon_isGoodMuon_TMOneStationAngTight;
      float Muon_isGoodMuon_TMLastStationTight;
      float Muon_isGoodMuon_TMLastStationAngTight;
      float Muon_isGoodMuon_TMLastStationOptimizedLowPtTight;
      float Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight;

      float  Muon_station_vars[8][7];

};
#endif
