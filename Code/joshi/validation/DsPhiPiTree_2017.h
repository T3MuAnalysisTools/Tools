#ifndef DsPhiPiTree_h
#define DsPhiPiTree_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class DsPhiPiTree : public Selection {

   public:
      DsPhiPiTree(TString Name_, TString id_);
      virtual ~DsPhiPiTree();

      virtual void  Configure();
      virtual void  Finish();

      /////////////////////////////////////////////////////////
      // This is a cut enumerator, for other event cuts please
      // fill the enumerator with put a new enumarator with an 
      // understandanle name, for exmaple  enum cuts {TriggerOk=0,
      // PrimeVts, EventCut1, EventCut2, ..., NCuts};  
      // Do not remove/rename  the last enumerator   NCuts;

      enum cuts {L1Fired, HLTFired, TwoMuTrkCandidate,OSMuons,Mu1PtCut,Mu2PtCut,MuonID,MuMuMassCut,TrackPtCut,NTrackHits, /*ChiSqCut,*/ DsMassCut,TriggerMatch,NCuts}; 

   protected:
      virtual void doEvent();  
      virtual void Store_ExtraDist();

   private:

      Float_t random_num;
      TRandom rndm;

      // PU Weights
      TFile* PUWeightFile;
      TH1D* puWeights;

      TFile * file;
      TTree * DsPhiPi_Tree;

      // Selection Variables
      float PhiMassHigh;
      float PhiMassLow;
      float sidebandDsMin;
      float sidebandDsMax;
      float peakDsMin;
      float peakDsMax;
      float dsMassMin;
      float dsMassMax;


      double nSidebands;
      double nPeak;

      unsigned int Muon_1_idx, Muon_2_idx, Track_idx;

      float genMatcheddR;

      // Branches
      float eventWeight;
      int era;
      int dataMCType;
      float var_relPt_iso05;
      float var_nTracks_iso05;
      float var_mindca_iso;
      float Muon1_TriggerMatchdR;
      float Muon2_TriggerMatchdR;
      float Track_TriggerMatchdR;
      int Muon1_isGlobal;
      int Muon2_isGlobal;
      int Muon1_isStandAlone;
      int Muon2_isStandAlone;
      int Muon1_isTracker;
      int Muon2_isTracker;
      int Muon1_isCalo;
      int Muon2_isCalo;
      int Muon1_isIsolationValid;
      int Muon2_isIsolationValid;
      int NVtx_woPUWeights;
      int NVtx;
      float PhiMass;
      float TripleMass;
      float Track_Pt;
      float Track_Eta;
      float Track_Phi;
      float TrackP;
      float Track_vx;
      float Track_vy;
      float Track_vz;
      float Track_normalizedChi2;
      float Track_numberOfValidHits;
      float Track_charge;
      float Track_dxy;
      float Track_dz;
      float Muon1_Pt;
      float Muon1_Eta;
      float Muon1_Phi;
      float Muon1_P;
      float Muon1_vx;
      float Muon1_vy;
      float Muon1_vz;
      float Muon1_cLP;
      float Muon1_SegmentCompatibility;
      float Muon1_NumberOfMatches;
      float Muon1_NumberOfMatchedStations;
      float Muon1_kink;
      float Muon2_Pt;
      float Muon2_Eta;
      float Muon2_Phi;
      float Muon2_P;
      float Muon2_vx;
      float Muon2_vy;
      float Muon2_vz;
      float Muon2_cLP;
      float Muon2_SegmentCompatibility;
      float Muon2_NumberOfMatches;
      float Muon2_NumberOfMatchedStations;
      float Muon2_kink;

      float DimuondR;
      float Muon1TrkdR;
      float Muon2TrkdR;
      float VertexKFChi2;
      float SVPVDsDirAngle;
      float NtracksClose;
      float NSV;
      float MinMuon_chi2LocalPosition;
      float MinDca;
      float MinD0SigSV;
      float MinD0SigPV;
      float MaxVertexPairQuality;
      float MaxdeltaMuZ;
      float MaxDca;
      float MaxD0SigSV;
      float MaxD0SigPV;
      float FLSignificance;
      float TransverseFLSignificance;
      float DecayLength;
      int var_isPrompt;

      float bssv_dxy;
      float bssv_significance;
      float pssv_dxy;
      float pvsv_dz;

      // Ds variables
      float ds_pt;
      float ds_eta;
      float ds_phi;
      float ds_e;
      float ds_motherPdgId; 

      // TrackerMuonId quantities
      float Muon1_ValidFraction;
      float Muon1_InnerNValidHits;
      float Muon1_InnerTrackQuality;
      float Muon1_NValidPixelHits;
      float Muon1_NValidTrackerHits;
      float Muon1_NLostTrackerHits;
      float Muon1_NLostTrackerHitsInner;
      float Muon1_NLostTrackerHitsOuter;
      float Muon1_PixelLayers;
      float Muon1_TrackerLayers;
      float Muon1_PtErrPt;
      float Muon1_CaloComp;
      float Muon1_HadS9;
      float Muon1_Had;
      float Muon1_EM;
      float Muon1_EMS9;
      float Muon1_EMS25;
      float Muon1_Kink;
      float Muon1_InnerNC2;

      float Muon2_ValidFraction;
      float Muon2_InnerNValidHits;
      float Muon2_InnerTrackQuality;
      float Muon2_NValidPixelHits;
      float Muon2_NValidTrackerHits;
      float Muon2_NLostTrackerHits;
      float Muon2_NLostTrackerHitsInner;
      float Muon2_NLostTrackerHitsOuter;
      float Muon2_PixelLayers;
      float Muon2_TrackerLayers;
      float Muon2_PtErrPt;
      float Muon2_CaloComp;
      float Muon2_HadS9;
      float Muon2_Had;
      float Muon2_EM;
      float Muon2_EMS9;
      float Muon2_EMS25;
      float Muon2_Kink;
      float Muon2_InnerNC2;

      float Track_ValidFraction;
      float Track_InnerNValidHits;
      float Track_InnerTrackQuality;
      float Track_NValidPixelHits;
      float Track_NValidTrackerHits;
      float Track_NLostTrackerHits;
      float Track_NLostTrackerHitsInner;
      float Track_NLostTrackerHitsOuter;
      float Track_PixelLayers;
      float Track_TrackerLayers;
      float Track_PtErrPt;
      float Track_CaloComp;
      float Track_HadS9;
      float Track_Had;
      float Track_EM;
      float Track_EMS9;
      float Track_EMS25;
      float Track_Kink;
      float Track_InnerNC2;

      // bit to check the status of the track
      bool Track_isTrackerMuon;
      bool Track_isGlobalMuon;
      bool Track_muonMatched;

      // Muon timing variables
      float Muon1_timeAtIpInOutErr;
      float Muon1_timeAtIpInOut;
      float Muon2_timeAtIpInOutErr;
      float Muon2_timeAtIpInOut;
      float Muon3_timeAtIpInOutErr;
      float Muon3_timeAtIpInOut;

      // Muon station variables
      float Muon1_station_vars[8][7];
      float Muon2_station_vars[8][7];

      TMatrixTSym<double> pv_cov;
      TMatrixTSym<double> sv_cov;

};
#endif
