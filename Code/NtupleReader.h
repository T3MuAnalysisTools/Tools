//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 17 13:12:55 2019 by ROOT version 6.10/05
// from TTree t3mtree/
// found on file: DsT3MNtuple_39.root
//////////////////////////////////////////////////////////

#ifndef NtupleReader_h
#define NtupleReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class NtupleReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          Event_EventNumber;
   UInt_t          Event_RunNumber;
   Int_t           Event_bunchCrossing;
   Int_t           Event_orbitNumber;
   UInt_t          Event_luminosityBlock;
   Bool_t          Event_isRealData;
   UInt_t          Event_nsignal_candidates;
   UInt_t          Event_ndsphipi_candidate;
   Int_t           Event_DataMC_Type;
   Double_t        puN;
   std::vector<std::vector<double> > *Track_p4;
   std::vector<double>  *Track_normalizedChi2;
   std::vector<double>  *Track_numberOfValidHits;
   std::vector<double>  *Track_charge;
   std::vector<double>  *Track_dxy;
   std::vector<double>  *Track_dz;
   std::vector<std::vector<double> > *Track_poca;
   std::vector<double>  *Track_dxyError;
   std::vector<double>  *Track_dzError;
   std::vector<std::vector<double> > *Muon_p4;
   std::vector<std::vector<double> > *Muon_Poca;
   std::vector<bool>    *Muon_isGlobalMuon;
   std::vector<bool>    *Muon_isStandAloneMuon;
   std::vector<bool>    *Muon_isTrackerMuon;
   std::vector<bool>    *Muon_isCaloMuon;
   std::vector<bool>    *Muon_isIsolationValid;
   std::vector<bool>    *Muon_isQualityValid;
   std::vector<bool>    *Muon_isTimeValid;
   std::vector<float>   *Muon_emEt03;
   std::vector<float>   *Muon_emVetoEt03;
   std::vector<float>   *Muon_hadEt03;
   std::vector<float>   *Muon_hadVetoEt03;
   std::vector<int>     *Muon_nJets03;
   std::vector<int>     *Muon_nTracks03;
   std::vector<float>   *Muon_sumPt03;
   std::vector<float>   *Muon_trackerVetoPt03;
   std::vector<float>   *Muon_emEt05;
   std::vector<float>   *Muon_emVetoEt05;
   std::vector<float>   *Muon_hadEt05;
   std::vector<float>   *Muon_hadVetoEt05;
   std::vector<int>     *Muon_nJets05;
   std::vector<int>     *Muon_nTracks05;
   std::vector<float>   *Muon_sumPt05;
   std::vector<float>   *Muon_trackerVetoPt05;
   std::vector<float>   *Muon_sumChargedHadronPt03;
   std::vector<float>   *Muon_sumChargedParticlePt03;
   std::vector<float>   *Muon_sumNeutralHadronEt03;
   std::vector<float>   *Muon_sumNeutralHadronEtHighThreshold03;
   std::vector<float>   *Muon_sumPhotonEt03;
   std::vector<float>   *Muon_sumPhotonEtHighThreshold03;
   std::vector<float>   *Muon_sumPUPt03;
   std::vector<float>   *Muon_sumChargedHadronPt04;
   std::vector<float>   *Muon_sumChargedParticlePt04;
   std::vector<float>   *Muon_sumNeutralHadronEt04;
   std::vector<float>   *Muon_sumNeutralHadronEtHighThreshold04;
   std::vector<float>   *Muon_sumPhotonEt04;
   std::vector<float>   *Muon_sumPhotonEtHighThreshold04;
   std::vector<float>   *Muon_sumPUPt04;
   std::vector<unsigned int> *Muon_Track_idx;
   std::vector<int>     *Muon_hitPattern_pixelLayerwithMeas;
   std::vector<int>     *Muon_numberOfMatchedStations;
   std::vector<float>   *Muon_normChi2;
   std::vector<int>     *Muon_hitPattern_numberOfValidMuonHits;
   std::vector<int>     *Muon_innerTrack_numberofValidHits;
   std::vector<int>     *Muon_numberOfMatches;
   std::vector<int>     *Muon_numberOfChambers;
   std::vector<bool>    *Muon_isPFMuon;
   std::vector<bool>    *Muon_isRPCMuon;
   std::vector<int>     *Muon_numberofValidPixelHits;
   std::vector<int>     *Muon_trackerLayersWithMeasurement;
   std::vector<bool>    *Muon_combinedQuality_updatedSta;
   std::vector<double>  *Muon_combinedQuality_trkKink;
   std::vector<double>  *Muon_combinedQuality_glbKink;
   std::vector<double>  *Muon_combinedQuality_trkRelChi2;
   std::vector<double>  *Muon_combinedQuality_staRelChi2;
   std::vector<double>  *Muon_combinedQuality_chi2LocalPosition;
   std::vector<double>  *Muon_combinedQuality_chi2LocalMomentum;
   std::vector<double>  *Muon_combinedQuality_localDistance;
   std::vector<double>  *Muon_combinedQuality_globalDeltaEtaPhi;
   std::vector<bool>    *Muon_combinedQuality_tightMatch;
   std::vector<double>  *Muon_combinedQuality_glbTrackProbability;
   std::vector<double>  *Muon_prod_inner_outer_charge;
   std::vector<double>  *Muon_ptError;
   std::vector<double>  *Muon_phiError;
   std::vector<double>  *Muon_etaError;
   std::vector<std::vector<double> > *Muon_outerTrack_p4;
   std::vector<std::vector<double> > *Muon_innerTrack_p4;
   std::vector<double>  *Muon_innerTrack_quality;
   std::vector<double>  *Muon_ptErrOverPt;
   std::vector<double>  *Muon_calEnergy_hadS9;
   std::vector<double>  *Muon_calEnergy_had;
   std::vector<double>  *Muon_calEnergy_emS25;
   std::vector<double>  *Muon_calEnergy_emS9;
   std::vector<double>  *Muon_calEnergy_em;
   std::vector<double>    *Muon_segmentCompatibility;
   std::vector<double>    *Muon_caloCompatibility;
   std::vector<bool>    *Muon_isGoodMuon_TM2DCompatibility;
   std::vector<double>  *Muon_innerTrack_validFraction;
   std::vector<double>  *Muon_innerTrack_pixelLayersWithMeasurement;
   std::vector<double>  *Muon_innerTrack_numberOfValidTrackerHits;
   std::vector<double>  *Muon_innerTrack_numberOfLostTrackerHits;
   std::vector<double>  *Muon_innerTrack_numberOfLostTrackerInnerHits;
   std::vector<double>  *Muon_innerTrack_numberOfLostTrackerOuterHits;
   std::vector<double>  *Muon_innerTrack_normalizedChi2;
   std::vector<double>  *Muon_outerTrack_normalizedChi2;
   std::vector<double>  *Muon_outerTrack_muonStationsWithValidHits;
   std::vector<bool>    *Muon_isGoodMuon_TrackerMuonArbitrated;
   std::vector<bool>    *Muon_isGoodMuon_TMOneStationTight;
   std::vector<bool>    *Muon_isGoodMuon_TMOneStationAngTight;
   std::vector<bool>    *Muon_isGoodMuon_TMLastStationTight;
   std::vector<bool>    *Muon_isGoodMuon_TMLastStationAngTight;
   std::vector<bool>    *Muon_isGoodMuon_TMLastStationOptimizedLowPtTight;
   std::vector<bool>    *Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight;
   std::vector<double>  *Muon_vmuonhitcomb_reco;
   std::vector<double>  *Muon_rpchits_reco;
   std::vector<int>     *Muon_ID;
   std::vector<int>     *Muon_StandardSelection;
   std::vector<int>     *Muon_charge;
   std::vector<int>     *Muon_trackCharge;
   std::vector<int>     *Muon_pdgid;
   std::vector<double>  *Muon_B;
   std::vector<double>  *Muon_M;
   std::vector<std::vector<double> > *Muon_par;
   std::vector<std::vector<double> > *Muon_cov;
   std::vector<std::vector<float> > *MC_p4;
   std::vector<int>     *MC_pdgid;
   std::vector<int>     *MC_charge;
   std::vector<int>     *MC_midx;
   std::vector<std::vector<int> > *MC_childpdgid;
   std::vector<std::vector<int> > *MC_childidx;
   std::vector<int>     *MC_status;

   std::vector<int> *signalTau_charge;
   std::vector<int> *signalTau_pdgid;
   std::vector<double> *signalTau_B;
   std::vector<double> *signalTau_M;
   std::vector<std::vector<double> > *signalTau_lvp;
   std::vector<std::vector<double> > *signalTau_cov;
   std::vector<int> *signalTau_isLVP;


   Bool_t          MC_isReco;
   std::vector<std::vector<double> > *MCSignalParticle_p4;
   std::vector<int>     *MCSignalParticle_pdgid;
   std::vector<std::vector<int> >   *MCSignalParticle_Sourcepdgid;
   std::vector<std::vector<std::vector<double> > > *MCSignalParticle_Sourcep4;

   std::vector<std::vector<int> > *MCSignalParticle_childpdgid;

   std::vector<std::vector<std::vector<double> > > *MCSignalParticle_childp4;
   std::vector<int>     *MCSignalParticle_charge;
   std::vector<std::vector<unsigned int> > *MCSignalParticle_Tauidx;
   std::vector<std::vector<std::vector<double> > > *MCTauandProd_p4;
   std::vector<std::vector<int> > *MCTauandProd_pdgid;
   std::vector<unsigned int> *MCTauandProd_midx;
   std::vector<std::vector<int> > *MCTauandProd_charge;
   std::vector<std::vector<unsigned int> > *ThreeMuons_index;
   std::vector<double>  *ThreeMuons_SV_Chi2;
   std::vector<double>  *ThreeMuons_SV_NDF;
   std::vector<std::vector<float> > *ThreeMuons_TriggerMatch_dR;
   std::vector<std::vector<unsigned int> > *TwoMuonsTrack_Muonsindex;
   std::vector<std::vector<unsigned int> > *TwoMuonsTrack_Trackindex;
   std::vector<double>  *TwoMuonsTrack_SV_Chi2;
   std::vector<double>  *TwoMuonsTrack_SV_NDF;
   std::vector<std::vector<float> > *TwoMuonsTrack_TriggerMatch_dR;
   std::vector<double>  *Jet_BTagCVSB;
   std::vector<double>  *Jet_BTagMVA;
   std::vector<double>  *Jet_BTagCSV;
   std::vector<std::vector<double> > *Jet_p4;
   Double_t        Vertex_N_primary;
   std::vector<std::vector<double> > *Vertex_signal_dca_reco;
   std::vector<std::vector<double> > *Vertex_signal_KF_pos;
   std::vector<std::vector<double> > *Vertex_signal_KF_cov;
   std::vector<std::vector<std::vector<double> > > *Vertex_signal_KF_refittedTracksP4;
   std::vector<double>  *Vertex_signal_KF_Chi2;
   std::vector<std::vector<double> > *Vertex_signal_AF_pos;
   std::vector<double>  *Vertex_signal_AF_Chi2;
   std::vector<double>  *Vertex_signal_AF_Ndf;
   std::vector<std::vector<double> > *Vertex_pair_quality;
   std::vector<std::vector<double> > *Vertex_pairfit_status;
   std::vector<std::vector<double> > *Vertex_MatchedPrimaryVertex;
   std::vector<bool>    *Vertex_RefitPVisValid;
   std::vector<std::vector<double> > *Vertex_MatchedRefitPrimaryVertex;
   std::vector<std::vector<double> > *Vertex_MatchedRefitPrimaryVertex_covariance;
   std::vector<std::vector<double> > *Vertex_d0_reco;
   std::vector<std::vector<double> > *Vertex_d0sig_reco;
   std::vector<std::vector<double> > *Vertex_2Ddisplacement;
   std::vector<std::vector<double> > *Vertex_3Ddisplacement;
   std::vector<std::vector<float> > *Vertex_Isolation1;
   std::vector<std::vector<float> > *Vertex_Isolation2;
   std::vector<std::vector<float> > *Vertex_Isolation3;
   std::vector<std::vector<float> > *Vertex_Isolation4;
   std::vector<std::string>  *Trigger_l1name;
   std::vector<int>     *Trigger_l1decision;
   std::vector<int>     *Trigger_l1prescale;
   std::vector<std::string>  *Trigger_hltname;
   std::vector<int>     *Trigger_hltdecision;

   // List of branches
   TBranch        *b_Event_EventNumber;   //!
   TBranch        *b_Event_RunNumber;   //!
   TBranch        *b_Event_bunchCrossing;   //!
   TBranch        *b_Event_orbitNumber;   //!
   TBranch        *b_Event_luminosityBlock;   //!
   TBranch        *b_Event_isRealData;   //!
   TBranch        *b_Event_nsignal_candidates;   //!
   TBranch        *b_Event_ndsphipi_candidate;   //!
   TBranch        *b_Event_DataMC_Type;   //!
   TBranch        *b_puN;   //!
   TBranch        *b_Track_p4;   //!
   TBranch        *b_Track_normalizedChi2;   //!
   TBranch        *b_Track_numberOfValidHits;   //!
   TBranch        *b_Track_charge;   //!
   TBranch        *b_Track_dxy;   //!
   TBranch        *b_Track_dz;   //!
   TBranch        *b_Track_poca;   //!
   TBranch        *b_Track_dxyError;   //!
   TBranch        *b_Track_dzError;   //!
   TBranch        *b_Muon_p4;   //!
   TBranch        *b_Muon_Poca;   //!
   TBranch        *b_Muon_isGlobalMuon;   //!
   TBranch        *b_Muon_isStandAloneMuon;   //!
   TBranch        *b_Muon_isTrackerMuon;   //!
   TBranch        *b_Muon_isCaloMuon;   //!
   TBranch        *b_Muon_isIsolationValid;   //!
   TBranch        *b_Muon_isQualityValid;   //!
   TBranch        *b_Muon_isTimeValid;   //!
   TBranch        *b_Muon_emEt03;   //!
   TBranch        *b_Muon_emVetoEt03;   //!
   TBranch        *b_Muon_hadEt03;   //!
   TBranch        *b_Muon_hadVetoEt03;   //!
   TBranch        *b_Muon_nJets03;   //!
   TBranch        *b_Muon_nTracks03;   //!
   TBranch        *b_Muon_sumPt03;   //!
   TBranch        *b_Muon_trackerVetoPt03;   //!
   TBranch        *b_Muon_emEt05;   //!
   TBranch        *b_Muon_emVetoEt05;   //!
   TBranch        *b_Muon_hadEt05;   //!
   TBranch        *b_Muon_hadVetoEt05;   //!
   TBranch        *b_Muon_nJets05;   //!
   TBranch        *b_Muon_nTracks05;   //!
   TBranch        *b_Muon_sumPt05;   //!
   TBranch        *b_Muon_trackerVetoPt05;   //!
   TBranch        *b_Muon_sumChargedHadronPt03;   //!
   TBranch        *b_Muon_sumChargedParticlePt03;   //!
   TBranch        *b_Muon_sumNeutralHadronEt03;   //!
   TBranch        *b_Muon_sumNeutralHadronEtHighThreshold03;   //!
   TBranch        *b_Muon_sumPhotonEt03;   //!
   TBranch        *b_Muon_sumPhotonEtHighThreshold03;   //!
   TBranch        *b_Muon_sumPUPt03;   //!
   TBranch        *b_Muon_sumChargedHadronPt04;   //!
   TBranch        *b_Muon_sumChargedParticlePt04;   //!
   TBranch        *b_Muon_sumNeutralHadronEt04;   //!
   TBranch        *b_Muon_sumNeutralHadronEtHighThreshold04;   //!
   TBranch        *b_Muon_sumPhotonEt04;   //!
   TBranch        *b_Muon_sumPhotonEtHighThreshold04;   //!
   TBranch        *b_Muon_sumPUPt04;   //!
   TBranch        *b_Muon_Track_idx;   //!
   TBranch        *b_Muon_hitPattern_pixelLayerwithMeas;   //!
   TBranch        *b_Muon_numberOfMatchedStations;   //!
   TBranch        *b_Muon_normChi2;   //!
   TBranch        *b_Muon_hitPattern_numberOfValidMuonHits;   //!
   TBranch        *b_Muon_innerTrack_numberofValidHits;   //!
   TBranch        *b_Muon_numberOfMatches;   //!
   TBranch        *b_Muon_numberOfChambers;   //!
   TBranch        *b_Muon_isPFMuon;   //!
   TBranch        *b_Muon_isRPCMuon;   //!
   TBranch        *b_Muon_numberofValidPixelHits;   //!
   TBranch        *b_Muon_trackerLayersWithMeasurement;   //!
   TBranch        *b_Muon_combinedQuality_updatedSta;   //!
   TBranch        *b_Muon_combinedQuality_trkKink;   //!
   TBranch        *b_Muon_combinedQuality_glbKink;   //!
   TBranch        *b_Muon_combinedQuality_trkRelChi2;   //!
   TBranch        *b_Muon_combinedQuality_staRelChi2;   //!
   TBranch        *b_Muon_combinedQuality_chi2LocalPosition;   //!
   TBranch        *b_Muon_combinedQuality_chi2LocalMomentum;   //!
   TBranch        *b_Muon_combinedQuality_localDistance;   //!
   TBranch        *b_Muon_combinedQuality_globalDeltaEtaPhi;   //!
   TBranch        *b_Muon_combinedQuality_tightMatch;   //!
   TBranch        *b_Muon_combinedQuality_glbTrackProbability;   //!
   TBranch        *b_Muon_prod_inner_outer_charge;   //!
   TBranch        *b_Muon_ptError;   //!
   TBranch        *b_Muon_phiError; //!
   TBranch        *b_Muon_etaError;   //!
   TBranch        *b_Muon_outerTrack_p4;   //!
   TBranch        *b_Muon_innerTrack_p4;   //!
   TBranch        *b_Muon_innerTrack_quality;   //!
   TBranch        *b_Muon_ptErrOverPt;   //!
   TBranch        *b_Muon_calEnergy_hadS9;   //!
   TBranch        *b_Muon_calEnergy_had;   //!
   TBranch        *b_Muon_calEnergy_emS25;   //!
   TBranch        *b_Muon_calEnergy_emS9;   //!
   TBranch        *b_Muon_calEnergy_em;   //!
   TBranch        *b_Muon_segmentCompatibility;   //!
   TBranch        *b_Muon_caloCompatibility;   //!
   TBranch        *b_Muon_isGoodMuon_TM2DCompatibility;   //!
   TBranch        *b_Muon_innerTrack_validFraction;   //!
   TBranch        *b_Muon_innerTrack_pixelLayersWithMeasurement;   //!
   TBranch        *b_Muon_innerTrack_numberOfValidTrackerHits;   //!
   TBranch        *b_Muon_innerTrack_numberOfLostTrackerHits;   //!
   TBranch        *b_Muon_innerTrack_numberOfLostTrackerInnerHits;   //!
   TBranch        *b_Muon_innerTrack_numberOfLostTrackerOuterHits;   //!
   TBranch        *b_Muon_innerTrack_normalizedChi2;   //!
   TBranch        *b_Muon_outerTrack_normalizedChi2;   //!
   TBranch        *b_Muon_outerTrack_muonStationsWithValidHits;   //!
   TBranch        *b_Muon_isGoodMuon_TrackerMuonArbitrated;   //!
   TBranch        *b_Muon_isGoodMuon_TMOneStationTight;   //!
   TBranch        *b_Muon_isGoodMuon_TMOneStationAngTight;   //!
   TBranch        *b_Muon_isGoodMuon_TMLastStationTight;   //!
   TBranch        *b_Muon_isGoodMuon_TMLastStationAngTight;   //!
   TBranch        *b_Muon_isGoodMuon_TMLastStationOptimizedLowPtTight;   //!
   TBranch        *b_Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight;   //!
   TBranch        *b_Muon_vmuonhitcomb_reco;   //!
   TBranch        *b_Muon_rpchits_reco;   //!
   TBranch        *b_Muon_ID;   //!
   TBranch        *b_Muon_StandardSelection;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_trackCharge;   //!
   TBranch        *b_Muon_pdgid;   //!
   TBranch        *b_Muon_B;   //!
   TBranch        *b_Muon_M;   //!
   TBranch        *b_Muon_par;   //!
   TBranch        *b_Muon_cov;   //!
   TBranch        *b_MC_p4;   //!
   TBranch        *b_MC_pdgid;   //!
   TBranch        *b_MC_charge;   //!
   TBranch        *b_MC_midx;   //!
   TBranch        *b_MC_childpdgid;   //!
   TBranch        *b_MC_childidx;   //!
   TBranch        *b_MC_status;   //!
   TBranch        *b_MC_isReco;   //!

   TBranch        *b_signalTau_charge;   //!  
   TBranch        *b_signalTau_pdgid;   //!
   TBranch        *b_signalTau_B;   //!
   TBranch        *b_signalTau_M;   //!
   TBranch        *b_signalTau_lvp;   //!
   TBranch        *b_signalTau_cov;   //!
   TBranch        *b_signalTau_isLVP;   //!

   TBranch        *b_MCSignalParticle_p4;   //!
   TBranch        *b_MCSignalParticle_pdgid;   //!

   TBranch        *b_MCSignalParticle_Sourcepdgid;    //!
   TBranch        *b_MCSignalParticle_Sourcep4;    //!

   TBranch        *b_MCSignalParticle_childpdgid;   //!
   TBranch        *b_MCSignalParticle_childp4;   //!
   TBranch        *b_MCSignalParticle_charge;   //!
   TBranch        *b_MCSignalParticle_Tauidx;   //!
   TBranch        *b_MCTauandProd_p4;   //!
   TBranch        *b_MCTauandProd_pdgid;   //!
   TBranch        *b_MCTauandProd_midx;   //!
   TBranch        *b_MCTauandProd_charge;   //!
   TBranch        *b_ThreeMuons_index;   //!
   TBranch        *b_ThreeMuons_SV_Chi2;   //!
   TBranch        *b_ThreeMuons_SV_NDF;   //!
   TBranch        *b_ThreeMuons_TriggerMatch_dR;   //!
   TBranch        *b_TwoMuonsTrack_Muonsindex;   //!
   TBranch        *b_TwoMuonsTrack_Trackindex;   //!
   TBranch        *b_TwoMuonsTrack_SV_Chi2;   //!
   TBranch        *b_TwoMuonsTrack_SV_NDF;   //!
   TBranch        *b_TwoMuonsTrack_TriggerMatch_dR;   //!
   TBranch        *b_Jet_BTagCVSB;   //!
   TBranch        *b_Jet_BTagMVA;   //!
   TBranch        *b_Jet_BTagCSV;   //!
   TBranch        *b_Jet_p4;   //!
   TBranch        *b_Vertex_N_primary;   //!
   TBranch        *b_Vertex_signal_dca_reco;   //!
   TBranch        *b_Vertex_signal_KF_pos;   //!
   TBranch        *b_Vertex_signal_KF_cov;   //!
   TBranch        *b_Vertex_signal_KF_refittedTracksP4;   //!
   TBranch        *b_Vertex_signal_KF_Chi2;   //!
   TBranch        *b_Vertex_signal_AF_pos;   //!
   TBranch        *b_Vertex_signal_AF_Chi2;   //!
   TBranch        *b_Vertex_signal_AF_Ndf;   //!
   TBranch        *b_Vertex_pair_quality;   //!
   TBranch        *b_Vertex_pairfit_status;   //!
   TBranch        *b_Vertex_MatchedPrimaryVertex;   //!
   TBranch        *b_Vertex_RefitPVisValid;   //!
   TBranch        *b_Vertex_MatchedRefitPrimaryVertex;   //!
   TBranch        *b_Vertex_MatchedRefitPrimaryVertex_covariance;   //!
   TBranch        *b_Vertex_d0_reco;   //!
   TBranch        *b_Vertex_d0sig_reco;   //!
   TBranch        *b_Vertex_2Ddisplacement;   //!
   TBranch        *b_Vertex_3Ddisplacement;   //!
   TBranch        *b_Vertex_Isolation1;   //!
   TBranch        *b_Vertex_Isolation2;   //!
   TBranch        *b_Vertex_Isolation3;   //!
   TBranch        *b_Vertex_Isolation4;   //!
   TBranch        *b_Trigger_l1name;   //!
   TBranch        *b_Trigger_l1decision;   //!
   TBranch        *b_Trigger_l1prescale;   //!
   TBranch        *b_Trigger_hltname;   //!
   TBranch        *b_Trigger_hltdecision;   //!

   NtupleReader(TTree *tree=0);
   virtual ~NtupleReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef NtupleReader_cxx
NtupleReader::NtupleReader(TTree *tree) : fChain(0)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {

#ifdef SINGLE_TREE
    // The following code should be used if you want this class to access
    // a single tree instead of a chain
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
    if (!f || !f->IsOpen()) {
      f = new TFile("Memory Directory");
    }
    f->GetObject("T3MTree/t3mtree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
    TChain * chain = new TChain("T3MTree/t3mtree","");
    tree = chain;
#endif // SINGLE_TREE

  }
  Init(tree);
}

NtupleReader::~NtupleReader()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t NtupleReader::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t NtupleReader::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}


void NtupleReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Track_p4 = 0;
   Track_normalizedChi2 = 0;
   Track_numberOfValidHits = 0;
   Track_charge = 0;
   Track_dxy = 0;
   Track_dz = 0;
   Track_poca = 0;
   Track_dxyError = 0;
   Track_dzError = 0;
   Muon_p4 = 0;
   Muon_Poca = 0;
   Muon_isGlobalMuon = 0;
   Muon_isStandAloneMuon = 0;
   Muon_isTrackerMuon = 0;
   Muon_isCaloMuon = 0;
   Muon_isIsolationValid = 0;
   Muon_isQualityValid = 0;
   Muon_isTimeValid = 0;
   Muon_emEt03 = 0;
   Muon_emVetoEt03 = 0;
   Muon_hadEt03 = 0;
   Muon_hadVetoEt03 = 0;
   Muon_nJets03 = 0;
   Muon_nTracks03 = 0;
   Muon_sumPt03 = 0;
   Muon_trackerVetoPt03 = 0;
   Muon_emEt05 = 0;
   Muon_emVetoEt05 = 0;
   Muon_hadEt05 = 0;
   Muon_hadVetoEt05 = 0;
   Muon_nJets05 = 0;
   Muon_nTracks05 = 0;
   Muon_sumPt05 = 0;
   Muon_trackerVetoPt05 = 0;
   Muon_sumChargedHadronPt03 = 0;
   Muon_sumChargedParticlePt03 = 0;
   Muon_sumNeutralHadronEt03 = 0;
   Muon_sumNeutralHadronEtHighThreshold03 = 0;
   Muon_sumPhotonEt03 = 0;
   Muon_sumPhotonEtHighThreshold03 = 0;
   Muon_sumPUPt03 = 0;
   Muon_sumChargedHadronPt04 = 0;
   Muon_sumChargedParticlePt04 = 0;
   Muon_sumNeutralHadronEt04 = 0;
   Muon_sumNeutralHadronEtHighThreshold04 = 0;
   Muon_sumPhotonEt04 = 0;
   Muon_sumPhotonEtHighThreshold04 = 0;
   Muon_sumPUPt04 = 0;
   Muon_Track_idx = 0;
   Muon_hitPattern_pixelLayerwithMeas = 0;
   Muon_numberOfMatchedStations = 0;
   Muon_normChi2 = 0;
   Muon_hitPattern_numberOfValidMuonHits = 0;
   Muon_innerTrack_numberofValidHits = 0;
   Muon_numberOfMatches = 0;
   Muon_numberOfChambers = 0;
   Muon_isPFMuon = 0;
   Muon_isRPCMuon = 0;
   Muon_numberofValidPixelHits = 0;
   Muon_trackerLayersWithMeasurement = 0;
   Muon_combinedQuality_updatedSta = 0;
   Muon_combinedQuality_trkKink = 0;
   Muon_combinedQuality_glbKink = 0;
   Muon_combinedQuality_trkRelChi2 = 0;
   Muon_combinedQuality_staRelChi2 = 0;
   Muon_combinedQuality_chi2LocalPosition = 0;
   Muon_combinedQuality_chi2LocalMomentum = 0;
   Muon_combinedQuality_localDistance = 0;
   Muon_combinedQuality_globalDeltaEtaPhi = 0;
   Muon_combinedQuality_tightMatch = 0;
   Muon_combinedQuality_glbTrackProbability = 0;
   Muon_prod_inner_outer_charge = 0;
   Muon_ptError = 0;
   Muon_phiError = 0;
   Muon_etaError = 0;
   Muon_outerTrack_p4 = 0;
   Muon_innerTrack_p4 = 0;
   Muon_innerTrack_quality = 0;
   Muon_ptErrOverPt = 0;
   Muon_calEnergy_hadS9 = 0;
   Muon_calEnergy_had = 0;
   Muon_calEnergy_emS25 = 0;
   Muon_calEnergy_emS9 = 0;
   Muon_calEnergy_em = 0;
   Muon_segmentCompatibility = 0;
   Muon_caloCompatibility = 0;
   Muon_isGoodMuon_TM2DCompatibility = 0;
   Muon_innerTrack_validFraction = 0;
   Muon_innerTrack_pixelLayersWithMeasurement = 0;
   Muon_innerTrack_numberOfValidTrackerHits = 0;
   Muon_innerTrack_numberOfLostTrackerHits = 0;
   Muon_innerTrack_numberOfLostTrackerInnerHits = 0;
   Muon_innerTrack_numberOfLostTrackerOuterHits = 0;
   Muon_innerTrack_normalizedChi2 = 0;
   Muon_outerTrack_normalizedChi2 = 0;
   Muon_outerTrack_muonStationsWithValidHits = 0;
   Muon_isGoodMuon_TrackerMuonArbitrated = 0;
   Muon_isGoodMuon_TMOneStationTight = 0;
   Muon_isGoodMuon_TMOneStationAngTight = 0;
   Muon_isGoodMuon_TMLastStationTight = 0;
   Muon_isGoodMuon_TMLastStationAngTight = 0;
   Muon_isGoodMuon_TMLastStationOptimizedLowPtTight = 0;
   Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight = 0;
   Muon_vmuonhitcomb_reco = 0;
   Muon_rpchits_reco = 0;
   Muon_ID = 0;
   Muon_StandardSelection = 0;
   Muon_charge = 0;
   Muon_trackCharge = 0;
   Muon_pdgid = 0;
   Muon_B = 0;
   Muon_M = 0;
   Muon_par = 0;
   Muon_cov = 0;
   MC_p4 = 0;
   MC_pdgid = 0;
   MC_charge = 0;
   MC_midx = 0;
   MC_childpdgid = 0;
   MC_childidx = 0;
   MC_status = 0;

   signalTau_charge = 0;
   signalTau_pdgid = 0; 
   signalTau_B = 0; 
   signalTau_M = 0; 
   signalTau_lvp = 0;
   signalTau_cov = 0;
   signalTau_isLVP = 0; 

  
   MCSignalParticle_p4 = 0;
   MCSignalParticle_pdgid = 0;
   MCSignalParticle_Sourcepdgid = 0;
   MCSignalParticle_Sourcep4 = 0;

   MCSignalParticle_childpdgid = 0;
   MCSignalParticle_childp4 = 0;
   MCSignalParticle_charge = 0;
   MCSignalParticle_Tauidx = 0;
   MCTauandProd_p4 = 0;
   MCTauandProd_pdgid = 0;
   MCTauandProd_midx = 0;
   MCTauandProd_charge = 0;
   ThreeMuons_index = 0;
   ThreeMuons_SV_Chi2 = 0;
   ThreeMuons_SV_NDF = 0;
   ThreeMuons_TriggerMatch_dR = 0;
   TwoMuonsTrack_Muonsindex = 0;
   TwoMuonsTrack_Trackindex = 0;
   TwoMuonsTrack_SV_Chi2 = 0;
   TwoMuonsTrack_SV_NDF = 0;
   TwoMuonsTrack_TriggerMatch_dR = 0;
   Jet_BTagCVSB = 0;
   Jet_BTagMVA = 0;
   Jet_BTagCSV = 0;
   Jet_p4 = 0;
   Vertex_signal_dca_reco = 0;
   Vertex_signal_KF_pos = 0;
   Vertex_signal_KF_cov = 0;
   Vertex_signal_KF_refittedTracksP4 = 0;
   Vertex_signal_KF_Chi2 = 0;
   Vertex_signal_AF_pos = 0;
   Vertex_signal_AF_Chi2 = 0;
   Vertex_signal_AF_Ndf = 0;
   Vertex_pair_quality = 0;
   Vertex_pairfit_status = 0;
   Vertex_MatchedPrimaryVertex = 0;
   Vertex_RefitPVisValid = 0;
   Vertex_MatchedRefitPrimaryVertex = 0;
   Vertex_MatchedRefitPrimaryVertex_covariance = 0;
   Vertex_d0_reco = 0;
   Vertex_d0sig_reco = 0;
   Vertex_2Ddisplacement = 0;
   Vertex_3Ddisplacement = 0;
   Vertex_Isolation1 = 0;
   Vertex_Isolation2 = 0;
   Vertex_Isolation3 = 0;
   Vertex_Isolation4 = 0;
   Trigger_l1name = 0;
   Trigger_l1decision = 0;
   Trigger_l1prescale = 0;
   Trigger_hltname = 0;
   Trigger_hltdecision = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event_EventNumber", &Event_EventNumber, &b_Event_EventNumber);
   fChain->SetBranchAddress("Event_RunNumber", &Event_RunNumber, &b_Event_RunNumber);
   fChain->SetBranchAddress("Event_bunchCrossing", &Event_bunchCrossing, &b_Event_bunchCrossing);
   fChain->SetBranchAddress("Event_orbitNumber", &Event_orbitNumber, &b_Event_orbitNumber);
   fChain->SetBranchAddress("Event_luminosityBlock", &Event_luminosityBlock, &b_Event_luminosityBlock);
   fChain->SetBranchAddress("Event_isRealData", &Event_isRealData, &b_Event_isRealData);
   fChain->SetBranchAddress("Event_nsignal_candidates", &Event_nsignal_candidates, &b_Event_nsignal_candidates);
   fChain->SetBranchAddress("Event_ndsphipi_candidate", &Event_ndsphipi_candidate, &b_Event_ndsphipi_candidate);
   fChain->SetBranchAddress("Event_DataMC_Type", &Event_DataMC_Type, &b_Event_DataMC_Type);
   fChain->SetBranchAddress("puN", &puN, &b_puN);
   fChain->SetBranchAddress("Track_p4", &Track_p4, &b_Track_p4);
   fChain->SetBranchAddress("Track_normalizedChi2", &Track_normalizedChi2, &b_Track_normalizedChi2);
   fChain->SetBranchAddress("Track_numberOfValidHits", &Track_numberOfValidHits, &b_Track_numberOfValidHits);
   fChain->SetBranchAddress("Track_charge", &Track_charge, &b_Track_charge);
   fChain->SetBranchAddress("Track_dxy", &Track_dxy, &b_Track_dxy);
   fChain->SetBranchAddress("Track_dz", &Track_dz, &b_Track_dz);
   fChain->SetBranchAddress("Track_poca", &Track_poca, &b_Track_poca);
   fChain->SetBranchAddress("Track_dxyError", &Track_dxyError, &b_Track_dxyError);
   fChain->SetBranchAddress("Track_dzError", &Track_dzError, &b_Track_dzError);
   fChain->SetBranchAddress("Muon_p4", &Muon_p4, &b_Muon_p4);
   fChain->SetBranchAddress("Muon_Poca", &Muon_Poca, &b_Muon_Poca);
   fChain->SetBranchAddress("Muon_isGlobalMuon", &Muon_isGlobalMuon, &b_Muon_isGlobalMuon);
   fChain->SetBranchAddress("Muon_isStandAloneMuon", &Muon_isStandAloneMuon, &b_Muon_isStandAloneMuon);
   fChain->SetBranchAddress("Muon_isTrackerMuon", &Muon_isTrackerMuon, &b_Muon_isTrackerMuon);
   fChain->SetBranchAddress("Muon_isCaloMuon", &Muon_isCaloMuon, &b_Muon_isCaloMuon);
   fChain->SetBranchAddress("Muon_isIsolationValid", &Muon_isIsolationValid, &b_Muon_isIsolationValid);
   fChain->SetBranchAddress("Muon_isQualityValid", &Muon_isQualityValid, &b_Muon_isQualityValid);
   fChain->SetBranchAddress("Muon_isTimeValid", &Muon_isTimeValid, &b_Muon_isTimeValid);
   fChain->SetBranchAddress("Muon_emEt03", &Muon_emEt03, &b_Muon_emEt03);
   fChain->SetBranchAddress("Muon_emVetoEt03", &Muon_emVetoEt03, &b_Muon_emVetoEt03);
   fChain->SetBranchAddress("Muon_hadEt03", &Muon_hadEt03, &b_Muon_hadEt03);
   fChain->SetBranchAddress("Muon_hadVetoEt03", &Muon_hadVetoEt03, &b_Muon_hadVetoEt03);
   fChain->SetBranchAddress("Muon_nJets03", &Muon_nJets03, &b_Muon_nJets03);
   fChain->SetBranchAddress("Muon_nTracks03", &Muon_nTracks03, &b_Muon_nTracks03);
   fChain->SetBranchAddress("Muon_sumPt03", &Muon_sumPt03, &b_Muon_sumPt03);
   fChain->SetBranchAddress("Muon_trackerVetoPt03", &Muon_trackerVetoPt03, &b_Muon_trackerVetoPt03);
   fChain->SetBranchAddress("Muon_emEt05", &Muon_emEt05, &b_Muon_emEt05);
   fChain->SetBranchAddress("Muon_emVetoEt05", &Muon_emVetoEt05, &b_Muon_emVetoEt05);
   fChain->SetBranchAddress("Muon_hadEt05", &Muon_hadEt05, &b_Muon_hadEt05);
   fChain->SetBranchAddress("Muon_hadVetoEt05", &Muon_hadVetoEt05, &b_Muon_hadVetoEt05);
   fChain->SetBranchAddress("Muon_nJets05", &Muon_nJets05, &b_Muon_nJets05);
   fChain->SetBranchAddress("Muon_nTracks05", &Muon_nTracks05, &b_Muon_nTracks05);
   fChain->SetBranchAddress("Muon_sumPt05", &Muon_sumPt05, &b_Muon_sumPt05);
   fChain->SetBranchAddress("Muon_trackerVetoPt05", &Muon_trackerVetoPt05, &b_Muon_trackerVetoPt05);
   fChain->SetBranchAddress("Muon_sumChargedHadronPt03", &Muon_sumChargedHadronPt03, &b_Muon_sumChargedHadronPt03);
   fChain->SetBranchAddress("Muon_sumChargedParticlePt03", &Muon_sumChargedParticlePt03, &b_Muon_sumChargedParticlePt03);
   fChain->SetBranchAddress("Muon_sumNeutralHadronEt03", &Muon_sumNeutralHadronEt03, &b_Muon_sumNeutralHadronEt03);
   fChain->SetBranchAddress("Muon_sumNeutralHadronEtHighThreshold03", &Muon_sumNeutralHadronEtHighThreshold03, &b_Muon_sumNeutralHadronEtHighThreshold03);
   fChain->SetBranchAddress("Muon_sumPhotonEt03", &Muon_sumPhotonEt03, &b_Muon_sumPhotonEt03);
   fChain->SetBranchAddress("Muon_sumPhotonEtHighThreshold03", &Muon_sumPhotonEtHighThreshold03, &b_Muon_sumPhotonEtHighThreshold03);
   fChain->SetBranchAddress("Muon_sumPUPt03", &Muon_sumPUPt03, &b_Muon_sumPUPt03);
   fChain->SetBranchAddress("Muon_sumChargedHadronPt04", &Muon_sumChargedHadronPt04, &b_Muon_sumChargedHadronPt04);
   fChain->SetBranchAddress("Muon_sumChargedParticlePt04", &Muon_sumChargedParticlePt04, &b_Muon_sumChargedParticlePt04);
   fChain->SetBranchAddress("Muon_sumNeutralHadronEt04", &Muon_sumNeutralHadronEt04, &b_Muon_sumNeutralHadronEt04);
   fChain->SetBranchAddress("Muon_sumNeutralHadronEtHighThreshold04", &Muon_sumNeutralHadronEtHighThreshold04, &b_Muon_sumNeutralHadronEtHighThreshold04);
   fChain->SetBranchAddress("Muon_sumPhotonEt04", &Muon_sumPhotonEt04, &b_Muon_sumPhotonEt04);
   fChain->SetBranchAddress("Muon_sumPhotonEtHighThreshold04", &Muon_sumPhotonEtHighThreshold04, &b_Muon_sumPhotonEtHighThreshold04);
   fChain->SetBranchAddress("Muon_sumPUPt04", &Muon_sumPUPt04, &b_Muon_sumPUPt04);
   fChain->SetBranchAddress("Muon_Track_idx", &Muon_Track_idx, &b_Muon_Track_idx);
   fChain->SetBranchAddress("Muon_hitPattern_pixelLayerwithMeas", &Muon_hitPattern_pixelLayerwithMeas, &b_Muon_hitPattern_pixelLayerwithMeas);
   fChain->SetBranchAddress("Muon_numberOfMatchedStations", &Muon_numberOfMatchedStations, &b_Muon_numberOfMatchedStations);
   fChain->SetBranchAddress("Muon_normChi2", &Muon_normChi2, &b_Muon_normChi2);
   fChain->SetBranchAddress("Muon_hitPattern_numberOfValidMuonHits", &Muon_hitPattern_numberOfValidMuonHits, &b_Muon_hitPattern_numberOfValidMuonHits);
   fChain->SetBranchAddress("Muon_innerTrack_numberofValidHits", &Muon_innerTrack_numberofValidHits, &b_Muon_innerTrack_numberofValidHits);
   fChain->SetBranchAddress("Muon_numberOfMatches", &Muon_numberOfMatches, &b_Muon_numberOfMatches);
   fChain->SetBranchAddress("Muon_numberOfChambers", &Muon_numberOfChambers, &b_Muon_numberOfChambers);
   fChain->SetBranchAddress("Muon_isPFMuon", &Muon_isPFMuon, &b_Muon_isPFMuon);
   fChain->SetBranchAddress("Muon_isRPCMuon", &Muon_isRPCMuon, &b_Muon_isRPCMuon);
   fChain->SetBranchAddress("Muon_numberofValidPixelHits", &Muon_numberofValidPixelHits, &b_Muon_numberofValidPixelHits);
   fChain->SetBranchAddress("Muon_trackerLayersWithMeasurement", &Muon_trackerLayersWithMeasurement, &b_Muon_trackerLayersWithMeasurement);
   fChain->SetBranchAddress("Muon_combinedQuality_updatedSta", &Muon_combinedQuality_updatedSta, &b_Muon_combinedQuality_updatedSta);
   fChain->SetBranchAddress("Muon_combinedQuality_trkKink", &Muon_combinedQuality_trkKink, &b_Muon_combinedQuality_trkKink);
   fChain->SetBranchAddress("Muon_combinedQuality_glbKink", &Muon_combinedQuality_glbKink, &b_Muon_combinedQuality_glbKink);
   fChain->SetBranchAddress("Muon_combinedQuality_trkRelChi2", &Muon_combinedQuality_trkRelChi2, &b_Muon_combinedQuality_trkRelChi2);
   fChain->SetBranchAddress("Muon_combinedQuality_staRelChi2", &Muon_combinedQuality_staRelChi2, &b_Muon_combinedQuality_staRelChi2);
   fChain->SetBranchAddress("Muon_combinedQuality_chi2LocalPosition", &Muon_combinedQuality_chi2LocalPosition, &b_Muon_combinedQuality_chi2LocalPosition);
   fChain->SetBranchAddress("Muon_combinedQuality_chi2LocalMomentum", &Muon_combinedQuality_chi2LocalMomentum, &b_Muon_combinedQuality_chi2LocalMomentum);
   fChain->SetBranchAddress("Muon_combinedQuality_localDistance", &Muon_combinedQuality_localDistance, &b_Muon_combinedQuality_localDistance);
   fChain->SetBranchAddress("Muon_combinedQuality_globalDeltaEtaPhi", &Muon_combinedQuality_globalDeltaEtaPhi, &b_Muon_combinedQuality_globalDeltaEtaPhi);
   fChain->SetBranchAddress("Muon_combinedQuality_tightMatch", &Muon_combinedQuality_tightMatch, &b_Muon_combinedQuality_tightMatch);
   fChain->SetBranchAddress("Muon_combinedQuality_glbTrackProbability", &Muon_combinedQuality_glbTrackProbability, &b_Muon_combinedQuality_glbTrackProbability);
   fChain->SetBranchAddress("Muon_prod_inner_outer_charge", &Muon_prod_inner_outer_charge, &b_Muon_prod_inner_outer_charge);
   fChain->SetBranchAddress("Muon_ptError",&Muon_ptError,&b_Muon_ptError);
   fChain->SetBranchAddress("Muon_phiError",&Muon_phiError, &b_Muon_phiError);
   fChain->SetBranchAddress("Muon_etaError",&Muon_etaError, &b_Muon_etaError);
   fChain->SetBranchAddress("Muon_outerTrack_p4", &Muon_outerTrack_p4, &b_Muon_outerTrack_p4);
   fChain->SetBranchAddress("Muon_innerTrack_p4", &Muon_innerTrack_p4, &b_Muon_innerTrack_p4);
   fChain->SetBranchAddress("Muon_innerTrack_quality", &Muon_innerTrack_quality, &b_Muon_innerTrack_quality);
   fChain->SetBranchAddress("Muon_ptErrOverPt", &Muon_ptErrOverPt, &b_Muon_ptErrOverPt);
   fChain->SetBranchAddress("Muon_calEnergy_hadS9", &Muon_calEnergy_hadS9, &b_Muon_calEnergy_hadS9);
   fChain->SetBranchAddress("Muon_calEnergy_had", &Muon_calEnergy_had, &b_Muon_calEnergy_had);
   fChain->SetBranchAddress("Muon_calEnergy_emS25", &Muon_calEnergy_emS25, &b_Muon_calEnergy_emS25);
   fChain->SetBranchAddress("Muon_calEnergy_emS9", &Muon_calEnergy_emS9, &b_Muon_calEnergy_emS9);
   fChain->SetBranchAddress("Muon_calEnergy_em", &Muon_calEnergy_em, &b_Muon_calEnergy_em);
   fChain->SetBranchAddress("Muon_segmentCompatibility", &Muon_segmentCompatibility, &b_Muon_segmentCompatibility);
   fChain->SetBranchAddress("Muon_caloCompatibility", &Muon_caloCompatibility, &b_Muon_caloCompatibility);
   fChain->SetBranchAddress("Muon_isGoodMuon_TM2DCompatibility", &Muon_isGoodMuon_TM2DCompatibility, &b_Muon_isGoodMuon_TM2DCompatibility);
   fChain->SetBranchAddress("Muon_innerTrack_validFraction", &Muon_innerTrack_validFraction, &b_Muon_innerTrack_validFraction);
   fChain->SetBranchAddress("Muon_innerTrack_pixelLayersWithMeasurement", &Muon_innerTrack_pixelLayersWithMeasurement, &b_Muon_innerTrack_pixelLayersWithMeasurement);
   fChain->SetBranchAddress("Muon_innerTrack_numberOfValidTrackerHits", &Muon_innerTrack_numberOfValidTrackerHits, &b_Muon_innerTrack_numberOfValidTrackerHits);
   fChain->SetBranchAddress("Muon_innerTrack_numberOfLostTrackerHits", &Muon_innerTrack_numberOfLostTrackerHits, &b_Muon_innerTrack_numberOfLostTrackerHits);
   fChain->SetBranchAddress("Muon_innerTrack_numberOfLostTrackerInnerHits", &Muon_innerTrack_numberOfLostTrackerInnerHits, &b_Muon_innerTrack_numberOfLostTrackerInnerHits);
   fChain->SetBranchAddress("Muon_innerTrack_numberOfLostTrackerOuterHits", &Muon_innerTrack_numberOfLostTrackerOuterHits, &b_Muon_innerTrack_numberOfLostTrackerOuterHits);
   fChain->SetBranchAddress("Muon_innerTrack_normalizedChi2", &Muon_innerTrack_normalizedChi2, &b_Muon_innerTrack_normalizedChi2);
   fChain->SetBranchAddress("Muon_outerTrack_normalizedChi2", &Muon_outerTrack_normalizedChi2, &b_Muon_outerTrack_normalizedChi2);
   fChain->SetBranchAddress("Muon_outerTrack_muonStationsWithValidHits", &Muon_outerTrack_muonStationsWithValidHits, &b_Muon_outerTrack_muonStationsWithValidHits);
   fChain->SetBranchAddress("Muon_isGoodMuon_TrackerMuonArbitrated", &Muon_isGoodMuon_TrackerMuonArbitrated, &b_Muon_isGoodMuon_TrackerMuonArbitrated);
   fChain->SetBranchAddress("Muon_isGoodMuon_TMOneStationTight", &Muon_isGoodMuon_TMOneStationTight, &b_Muon_isGoodMuon_TMOneStationTight);
   fChain->SetBranchAddress("Muon_isGoodMuon_TMOneStationAngTight", &Muon_isGoodMuon_TMOneStationAngTight, &b_Muon_isGoodMuon_TMOneStationAngTight);
   fChain->SetBranchAddress("Muon_isGoodMuon_TMLastStationTight", &Muon_isGoodMuon_TMLastStationTight, &b_Muon_isGoodMuon_TMLastStationTight);
   fChain->SetBranchAddress("Muon_isGoodMuon_TMLastStationAngTight", &Muon_isGoodMuon_TMLastStationAngTight, &b_Muon_isGoodMuon_TMLastStationAngTight);
   fChain->SetBranchAddress("Muon_isGoodMuon_TMLastStationOptimizedLowPtTight", &Muon_isGoodMuon_TMLastStationOptimizedLowPtTight, &b_Muon_isGoodMuon_TMLastStationOptimizedLowPtTight);
   fChain->SetBranchAddress("Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight", &Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight, &b_Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight);
   fChain->SetBranchAddress("Muon_vmuonhitcomb_reco", &Muon_vmuonhitcomb_reco, &b_Muon_vmuonhitcomb_reco);
   fChain->SetBranchAddress("Muon_rpchits_reco", &Muon_rpchits_reco, &b_Muon_rpchits_reco);
   fChain->SetBranchAddress("Muon_ID", &Muon_ID, &b_Muon_ID);
   fChain->SetBranchAddress("Muon_StandardSelection", &Muon_StandardSelection, &b_Muon_StandardSelection);
   fChain->SetBranchAddress("Muon_charge", &Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_trackCharge", &Muon_trackCharge, &b_Muon_trackCharge);
   fChain->SetBranchAddress("Muon_pdgid", &Muon_pdgid, &b_Muon_pdgid);
   fChain->SetBranchAddress("Muon_B", &Muon_B, &b_Muon_B);
   fChain->SetBranchAddress("Muon_M", &Muon_M, &b_Muon_M);
   fChain->SetBranchAddress("Muon_par", &Muon_par, &b_Muon_par);
   fChain->SetBranchAddress("Muon_cov", &Muon_cov, &b_Muon_cov);
   fChain->SetBranchAddress("MC_p4", &MC_p4, &b_MC_p4);
   fChain->SetBranchAddress("MC_pdgid", &MC_pdgid, &b_MC_pdgid);
   fChain->SetBranchAddress("MC_charge", &MC_charge, &b_MC_charge);
   fChain->SetBranchAddress("MC_midx", &MC_midx, &b_MC_midx);
   fChain->SetBranchAddress("MC_childpdgid", &MC_childpdgid, &b_MC_childpdgid);
   fChain->SetBranchAddress("MC_childidx", &MC_childidx, &b_MC_childidx);
   fChain->SetBranchAddress("MC_status", &MC_status, &b_MC_status);
   fChain->SetBranchAddress("MC_isReco", &MC_isReco, &b_MC_isReco);

   fChain->SetBranchAddress("signalTau_charge",&signalTau_charge,&b_signalTau_charge);
   fChain->SetBranchAddress("signalTau_pdgid",&signalTau_pdgid,&b_signalTau_pdgid);
   fChain->SetBranchAddress("signalTau_B",&signalTau_B, &b_signalTau_B);
   fChain->SetBranchAddress("signalTau_M",&signalTau_M, &b_signalTau_M);
   fChain->SetBranchAddress("signalTau_lvp",&signalTau_lvp, &b_signalTau_lvp);
   fChain->SetBranchAddress("signalTau_cov",&signalTau_cov,&b_signalTau_cov);
   fChain->SetBranchAddress("signalTau_isLVP",&signalTau_isLVP, &b_signalTau_isLVP);



   fChain->SetBranchAddress("MCSignalParticle_p4", &MCSignalParticle_p4, &b_MCSignalParticle_p4);
   fChain->SetBranchAddress("MCSignalParticle_pdgid", &MCSignalParticle_pdgid, &b_MCSignalParticle_pdgid);

   fChain->SetBranchAddress("MCSignalParticle_Sourcepdgid",&MCSignalParticle_Sourcepdgid, &b_MCSignalParticle_Sourcepdgid);
   fChain->SetBranchAddress("MCSignalParticle_Sourcep4",&MCSignalParticle_Sourcep4, &b_MCSignalParticle_Sourcep4);




   fChain->SetBranchAddress("MCSignalParticle_childpdgid", &MCSignalParticle_childpdgid, &b_MCSignalParticle_childpdgid);
   fChain->SetBranchAddress("MCSignalParticle_childp4", &MCSignalParticle_childp4, &b_MCSignalParticle_childp4);
   fChain->SetBranchAddress("MCSignalParticle_charge", &MCSignalParticle_charge, &b_MCSignalParticle_charge);
   fChain->SetBranchAddress("MCSignalParticle_Tauidx", &MCSignalParticle_Tauidx, &b_MCSignalParticle_Tauidx);
   fChain->SetBranchAddress("MCTauandProd_p4", &MCTauandProd_p4, &b_MCTauandProd_p4);
   fChain->SetBranchAddress("MCTauandProd_pdgid", &MCTauandProd_pdgid, &b_MCTauandProd_pdgid);
   fChain->SetBranchAddress("MCTauandProd_midx", &MCTauandProd_midx, &b_MCTauandProd_midx);
   fChain->SetBranchAddress("MCTauandProd_charge", &MCTauandProd_charge, &b_MCTauandProd_charge);
   fChain->SetBranchAddress("ThreeMuons_index", &ThreeMuons_index, &b_ThreeMuons_index);
   fChain->SetBranchAddress("ThreeMuons_SV_Chi2", &ThreeMuons_SV_Chi2, &b_ThreeMuons_SV_Chi2);
   fChain->SetBranchAddress("ThreeMuons_SV_NDF", &ThreeMuons_SV_NDF, &b_ThreeMuons_SV_NDF);
   fChain->SetBranchAddress("ThreeMuons_TriggerMatch_dR", &ThreeMuons_TriggerMatch_dR, &b_ThreeMuons_TriggerMatch_dR);
   fChain->SetBranchAddress("TwoMuonsTrack_Muonsindex", &TwoMuonsTrack_Muonsindex, &b_TwoMuonsTrack_Muonsindex);
   fChain->SetBranchAddress("TwoMuonsTrack_Trackindex", &TwoMuonsTrack_Trackindex, &b_TwoMuonsTrack_Trackindex);
   fChain->SetBranchAddress("TwoMuonsTrack_SV_Chi2", &TwoMuonsTrack_SV_Chi2, &b_TwoMuonsTrack_SV_Chi2);
   fChain->SetBranchAddress("TwoMuonsTrack_SV_NDF", &TwoMuonsTrack_SV_NDF, &b_TwoMuonsTrack_SV_NDF);
   fChain->SetBranchAddress("TwoMuonsTrack_TriggerMatch_dR", &TwoMuonsTrack_TriggerMatch_dR, &b_TwoMuonsTrack_TriggerMatch_dR);
   fChain->SetBranchAddress("Jet_BTagCVSB", &Jet_BTagCVSB, &b_Jet_BTagCVSB);
   fChain->SetBranchAddress("Jet_BTagMVA", &Jet_BTagMVA, &b_Jet_BTagMVA);
   fChain->SetBranchAddress("Jet_BTagCSV", &Jet_BTagCSV, &b_Jet_BTagCSV);
   fChain->SetBranchAddress("Jet_p4", &Jet_p4, &b_Jet_p4);
   fChain->SetBranchAddress("Vertex_N_primary", &Vertex_N_primary, &b_Vertex_N_primary);
   fChain->SetBranchAddress("Vertex_signal_dca_reco", &Vertex_signal_dca_reco, &b_Vertex_signal_dca_reco);
   fChain->SetBranchAddress("Vertex_signal_KF_pos", &Vertex_signal_KF_pos, &b_Vertex_signal_KF_pos);
   fChain->SetBranchAddress("Vertex_signal_KF_cov", &Vertex_signal_KF_cov, &b_Vertex_signal_KF_cov);
   fChain->SetBranchAddress("Vertex_signal_KF_refittedTracksP4", &Vertex_signal_KF_refittedTracksP4, &b_Vertex_signal_KF_refittedTracksP4);
   fChain->SetBranchAddress("Vertex_signal_KF_Chi2", &Vertex_signal_KF_Chi2, &b_Vertex_signal_KF_Chi2);
   fChain->SetBranchAddress("Vertex_signal_AF_pos", &Vertex_signal_AF_pos, &b_Vertex_signal_AF_pos);
   fChain->SetBranchAddress("Vertex_signal_AF_Chi2", &Vertex_signal_AF_Chi2, &b_Vertex_signal_AF_Chi2);
   fChain->SetBranchAddress("Vertex_signal_AF_Ndf", &Vertex_signal_AF_Ndf, &b_Vertex_signal_AF_Ndf);
   fChain->SetBranchAddress("Vertex_pair_quality", &Vertex_pair_quality, &b_Vertex_pair_quality);
   fChain->SetBranchAddress("Vertex_pairfit_status", &Vertex_pairfit_status, &b_Vertex_pairfit_status);
   fChain->SetBranchAddress("Vertex_MatchedPrimaryVertex", &Vertex_MatchedPrimaryVertex, &b_Vertex_MatchedPrimaryVertex);
   fChain->SetBranchAddress("Vertex_RefitPVisValid", &Vertex_RefitPVisValid, &b_Vertex_RefitPVisValid);
   fChain->SetBranchAddress("Vertex_MatchedRefitPrimaryVertex", &Vertex_MatchedRefitPrimaryVertex, &b_Vertex_MatchedRefitPrimaryVertex);
   fChain->SetBranchAddress("Vertex_MatchedRefitPrimaryVertex_covariance", &Vertex_MatchedRefitPrimaryVertex_covariance, &b_Vertex_MatchedRefitPrimaryVertex_covariance);
   fChain->SetBranchAddress("Vertex_d0_reco", &Vertex_d0_reco, &b_Vertex_d0_reco);
   fChain->SetBranchAddress("Vertex_d0sig_reco", &Vertex_d0sig_reco, &b_Vertex_d0sig_reco);
   fChain->SetBranchAddress("Vertex_2Ddisplacement", &Vertex_2Ddisplacement, &b_Vertex_2Ddisplacement);
   fChain->SetBranchAddress("Vertex_3Ddisplacement", &Vertex_3Ddisplacement, &b_Vertex_3Ddisplacement);
   fChain->SetBranchAddress("Vertex_Isolation1", &Vertex_Isolation1, &b_Vertex_Isolation1);
   fChain->SetBranchAddress("Vertex_Isolation2", &Vertex_Isolation2, &b_Vertex_Isolation2);
   fChain->SetBranchAddress("Vertex_Isolation3", &Vertex_Isolation3, &b_Vertex_Isolation3);
   fChain->SetBranchAddress("Vertex_Isolation4", &Vertex_Isolation4, &b_Vertex_Isolation4);
   fChain->SetBranchAddress("Trigger_l1name", &Trigger_l1name, &b_Trigger_l1name);
   fChain->SetBranchAddress("Trigger_l1decision", &Trigger_l1decision, &b_Trigger_l1decision);
   fChain->SetBranchAddress("Trigger_l1prescale", &Trigger_l1prescale, &b_Trigger_l1prescale);
   fChain->SetBranchAddress("Trigger_hltname", &Trigger_hltname, &b_Trigger_hltname);
   fChain->SetBranchAddress("Trigger_hltdecision", &Trigger_hltdecision, &b_Trigger_hltdecision);
   Notify();
}

Bool_t NtupleReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NtupleReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NtupleReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef NtupleReader_cxx
