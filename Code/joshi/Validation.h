#ifndef DimuTrk_h
#define DimuTrk_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class DimuTrk : public Selection {

    public:
      DimuTrk(TString Name_, TString id_);
      virtual ~DimuTrk();
		
		virtual double deltaR(double, double, double, double);
      virtual void  Configure();
      virtual void  Finish();

      enum cuts {L1TOk=0,HLTOk,PrimeVtx,NCuts};

    protected:
      virtual void doEvent();
      virtual void Store_ExtraDist();

    private:
      
		// Selection Variables
      std::vector<TH1D> NVtx;
      std::vector<TH1D> NTracks;
		std::vector<TH1D> HLT_Tau3Mu;

		// Track Variables
      std::vector<TH1D> Track_P;
      std::vector<TH1D> Track_E;
      std::vector<TH1D> Track_Pt;
      std::vector<TH1D> Track_Eta;
      std::vector<TH1D> Track_Phi;
      std::vector<TH1D> Track_vx;
      std::vector<TH1D> Track_vy;
      std::vector<TH1D> Track_vz;
      std::vector<TH1D> Track_normalizedChi2;
      std::vector<TH1D> Track_numberOfValidHits;
      std::vector<TH1D> Track_charge;
      std::vector<TH1D> Track_dxy;
      std::vector<TH1D> Track_dz;
      std::vector<TH1D> Track_dxyError;
      std::vector<TH1D> Track_dzError;
		
		// Muon variables
		std::vector<TH1D> Muon_Pt;
      std::vector<TH1D> Muon_Eta;
      std::vector<TH1D> Muon_Phi;
  		std::vector<TH1D> Muon_E;
      std::vector<TH1D> Muon_P;
      std::vector<TH1D> Muon_vx;
      std::vector<TH1D> Muon_vy;
      std::vector<TH1D> Muon_vz;
      std::vector<TH1D> Muon_IsGlobalMuon;
      std::vector<TH1D> Muon_IsStandAloneMuon;
      std::vector<TH1D> Muon_IsTrackerMuon;
      std::vector<TH1D> Muon_IsCaloMuon;
      std::vector<TH1D> Muon_IsIsolationValid;
      std::vector<TH1D> Muon_IsQualityValid;
      std::vector<TH1D> Muon_IsTimeValid;
      std::vector<TH1D> Muon_IsPFMuon;
      std::vector<TH1D> Muon_IsRPCMuon;
      std::vector<TH1D> Muon_emEt03;
      std::vector<TH1D> Muon_emVetoEt03;
      std::vector<TH1D> Muon_hadEt03;
      std::vector<TH1D> Muon_hadVetoEt03;
      std::vector<TH1D> Muon_nJets03;
      std::vector<TH1D> Muon_nTracks03;
      std::vector<TH1D> Muon_StandardSelection;
      std::vector<TH1D> Muon_sumPt03;
      std::vector<TH1D> Muon_trackerVetoPt03;
      std::vector<TH1D> Muon_sumChargedHadronPt03;
      std::vector<TH1D> Muon_sumChargedParticlePt03;
      std::vector<TH1D> Muon_sumNeutralHadronEt03;
      std::vector<TH1D> Muon_sumNeutralHadronEtHighThreshold03;
      std::vector<TH1D> Muon_sumPhotonEt03;
      std::vector<TH1D> Muon_sumPhotonEtHighThreshold03;
      std::vector<TH1D> Muon_sumPUPt03;
      std::vector<TH1D> Muon_numberOfChambers;
      std::vector<TH1D> Muon_Track_idx;
      std::vector<TH1D> Muon_combinedQuality_updatedSta;
      std::vector<TH1D> Muon_combinedQuality_trkKink;
      std::vector<TH1D> Muon_combinedQuality_glbKink;
      std::vector<TH1D> Muon_combinedQuality_trkRelChi2;
      std::vector<TH1D> Muon_combinedQuality_staRelChi2;
      std::vector<TH1D> Muon_combinedQuality_chi2LocalPosition;
      std::vector<TH1D> Muon_combinedQuality_chi2LocalMomentum;
      std::vector<TH1D> Muon_combinedQuality_localDistance;
      std::vector<TH1D> Muon_combinedQuality_globalDeltaEtaPhi;
      std::vector<TH1D> Muon_combinedQuality_tightMatch;
      std::vector<TH1D> Muon_combinedQuality_glbTrackProbability;
      std::vector<TH1D> Muon_prod_inner_outer_charge;
      std::vector<TH1D> Muon_innerTrack_quality;
      std::vector<TH1D> Muon_ptErrOverPt;
      std::vector<TH1D> Muon_calEnergy_hadS9;
      std::vector<TH1D> Muon_calEnergy_had;
      std::vector<TH1D> Muon_calEnergy_emS25;
      std::vector<TH1D> Muon_calEnergy_emS9;
      std::vector<TH1D> Muon_calEnergy_em;
      std::vector<TH1D> Muon_charge;
      std::vector<TH1D> Muon_trackCharge;
      std::vector<TH1D> Muon_pdgid;
      std::vector<TH1D> Muon_B;
      std::vector<TH1D> Muon_M;
      std::vector<TH1D> Muon_par;
      std::vector<TH1D> Muon_cov;
      std::vector<TH1D> Muon_hitPattern_pixelLayerwithMeas;
      std::vector<TH1D> Muon_numberOfMatchedStations;
      std::vector<TH1D> Muon_normChi2;
      std::vector<TH1D> Muon_hitPattern_numberOfValidMuonHits;
      std::vector<TH1D> Muon_innerTrack_numberofValidHits;
      std::vector<TH1D> Muon_numberofValidPixelHits;
      std::vector<TH1D> Muon_numberOfMatches;
      std::vector<TH1D> Muon_trackerLayersWithMeasurement;
      std::vector<TH1D> Muon_segmentCompatibility;
      std::vector<TH1D> Muon_caloCompatibility;
      std::vector<TH1D> Muon_innerTrack_validFraction;
      std::vector<TH1D> Muon_innerTrack_pixelLayersWithMeasurement;
      std::vector<TH1D> Muon_innerTrack_numberOfValidTrackerHits;
      std::vector<TH1D> Muon_innerTrack_numberOfLostTrackerHits;
      std::vector<TH1D> Muon_innerTrack_numberOfLostTrackerInnerHits;
      std::vector<TH1D> Muon_innerTrack_numberOfLostTrackerOuterHits;
      std::vector<TH1D> Muon_innerTrack_normalizedChi2;
      std::vector<TH1D> Muon_vmuonhitcomb_reco;
      std::vector<TH1D> Muon_rpchits_reco;
      std::vector<TH1D> Muon_outerTrack_normalizedChi2;
      std::vector<TH1D> Muon_outerTrack_muonStationsWithValidHits;
      std::vector<TH1D> Muon_isGoodMuon_TM2DCompatibility;
      std::vector<TH1D> Muon_isGoodMuon_TrackerMuonArbitrated;
      std::vector<TH1D> Muon_isGoodMuon_TMOneStationTight;
      std::vector<TH1D> Muon_isGoodMuon_TMOneStationAngTight;
      std::vector<TH1D> Muon_isGoodMuon_TMLastStationTight;
      std::vector<TH1D> Muon_isGoodMuon_TMLastStationAngTight;
      std::vector<TH1D> Muon_isGoodMuon_TMLastStationOptimizedLowPtTight;
      std::vector<TH1D> Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight;

      std::vector<TH1D> PhiMass;
      std::vector<TH1D> TripleMass;
		std::vector<TH1D> MuonsPtRatio;
		std::vector<TH1D> DimuondR;
		std::vector<TH1D> Muon1TrkdR;
		std::vector<TH1D> Muon2TrkdR;
		std::vector<TH1D> Muon1_isGlobal;
      std::vector<TH1D> Muon2_isGlobal;
      std::vector<TH1D> Muon1_isStandAlone;
      std::vector<TH1D> Muon2_isStandAlone;
      std::vector<TH1D> Muon1_isTracker;
      std::vector<TH1D> Muon2_isTracker;
      std::vector<TH1D> Muon1_isCalo;
      std::vector<TH1D> Muon2_isCalo;
      std::vector<TH1D> Muon1_isIsolationValid;
      std::vector<TH1D> Muon2_isIsolationValid;
      std::vector<TH1D> Muon1_isTimeValid;
      std::vector<TH1D> Muon2_isTimeValid;
      std::vector<TH1D> Muon1_emEt03;
      std::vector<TH1D> Muon2_emEt03;
      std::vector<TH1D> Muon1_emVetoEt03;
      std::vector<TH1D> Muon2_emVetoEt03;
      std::vector<TH1D> Muon1_hadEt03;
      std::vector<TH1D> Muon2_hadEt03;
      std::vector<TH1D> Muon1_hadVetoEt03;
      std::vector<TH1D> Muon2_hadVetoEt03;
      std::vector<TH1D> Muon1_nJets03;
      std::vector<TH1D> Muon2_nJets03;
      std::vector<TH1D> Muon1_nTracks03;
      std::vector<TH1D> Muon2_nTracks03;
      std::vector<TH1D> Muon1_sumPt03;
      std::vector<TH1D> Muon2_sumPt03;
      std::vector<TH1D> Muon1_trackerVetoPt03;
      std::vector<TH1D> Muon2_trackerVetoPt03;
      std::vector<TH1D> Muon1_emEt05;
      std::vector<TH1D> Muon2_emEt05;
      std::vector<TH1D> Muon1_emVetoEt05;
      std::vector<TH1D> Muon2_emVetoEt05;
      std::vector<TH1D> Muon1_hadEt05;
      std::vector<TH1D> Muon2_hadEt05;
      std::vector<TH1D> Muon1_hadVetoEt05;
      std::vector<TH1D> Muon2_hadVetoEt05;
      std::vector<TH1D> Muon1_nJets05;
      std::vector<TH1D> Muon2_nJets05;
      std::vector<TH1D> Muon1_nTracks05;
      std::vector<TH1D> Muon2_nTracks05;
      std::vector<TH1D> Muon1_sumPt05;
      std::vector<TH1D> Muon2_sumPt05;
      std::vector<TH1D> Muon1_trackerVetoPt05;
      std::vector<TH1D> Muon2_trackerVetoPt05;
      std::vector<TH1D> Muon1_sumChargedHadronPt03;
      std::vector<TH1D> Muon2_sumChargedHadronPt03;
      std::vector<TH1D> Muon1_sumChargedParticlePt03;
      std::vector<TH1D> Muon2_sumChargedParticlePt03;
      std::vector<TH1D> Muon1_sumNeutralHadronEt03;
      std::vector<TH1D> Muon2_sumNeutralHadronEt03;
      std::vector<TH1D> Muon1_sumNeutralHadronEtHighThreshold03;
      std::vector<TH1D> Muon2_sumNeutralHadronEtHighThreshold03;
      std::vector<TH1D> Muon1_sumPhotonEt03;
      std::vector<TH1D> Muon2_sumPhotonEt03;
      std::vector<TH1D> Muon1_sumPhotonEtHighThreshold03;
      std::vector<TH1D> Muon2_sumPhotonEtHighThreshold03;
      std::vector<TH1D> Muon1_sumPUPt03;
      std::vector<TH1D> Muon2_sumPUPt03;
      std::vector<TH1D> Muon1_sumChargedHadronPt04;
      std::vector<TH1D> Muon2_sumChargedHadronPt04;
      std::vector<TH1D> Muon1_sumChargedParticlePt04;
      std::vector<TH1D> Muon2_sumChargedParticlePt04;
      std::vector<TH1D> Muon1_sumNeutralHadronEt04;
      std::vector<TH1D> Muon2_sumNeutralHadronEt04;
      std::vector<TH1D> Muon1_sumNeutralHadronEtHighThreshold04;
      std::vector<TH1D> Muon2_sumNeutralHadronEtHighThreshold04;
      std::vector<TH1D> Muon1_sumPhotonEt04;
      std::vector<TH1D> Muon2_sumPhotonEt04;
      std::vector<TH1D> Muon1_sumPhotonEtHighThreshold04;
      std::vector<TH1D> Muon2_sumPhotonEtHighThreshold04;
      std::vector<TH1D> Muon1_sumPUPt04;
      std::vector<TH1D> Muon2_sumPUPt04;

      // 2D histograms
      std::vector<TH2D> PhiMassVsDsMass;

};
#endif
