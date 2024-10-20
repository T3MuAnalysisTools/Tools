#ifndef BackgroundSelector_h
#define BackgroundSelector_h

#include "Selection.h"
#include <vector>
#include "TString.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


class BackgroundSelector : public Selection {

 public:
  BackgroundSelector(TString Name_, TString id_);
  virtual ~BackgroundSelector();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {L1TOk=0,HLTOk,SignalCandidate, Mu1PtCut, Mu2PtCut, Mu3PtCut, MuonID,PhiVeto1, OmegaVeto1, PhiVeto2, OmegaVeto2, TriggerMatch, TauMassCut,NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:
  double tauMinMass_, tauMaxMass_;
  double tauMinSideBand_,tauMaxSideBand_;

  bool RunB, RunC, RunD, RunE, RunF = 0;

  TTree *T3MMiniTree;
  TFile *T3MFMiniTree;

  TRandom rndm;
  double random_num;

  // Selection Variables
  // Initializhere your analysis histograms

  std::vector<TH1D> Muon1Pt;
  std::vector<TH1D> Muon2Pt;
  std::vector<TH1D> Muon3Pt;

  std::vector<TH1D> Muon1Eta;
  std::vector<TH1D> Muon2Eta;
  std::vector<TH1D> Muon3Eta;


  std::vector<TH1D> TauEta;
  std::vector<TH1D> TauPt;
  std::vector<TH1D> TauP;


  std::vector<TH1D> VertexPairDistance1;
  std::vector<TH1D> VertexPairDistance2;
  std::vector<TH1D> VertexPairDistance3;




  std::vector<TH1D> TauMassResolution;
  std::vector<TH1D> EventMassResolution_PtEtaPhi;

  std::vector<TH1D> SVPVTauDirAngle;

  std::vector<TH1D> TauMassRefitA1;
  std::vector<TH1D> TauMassRefitA1MassCut;
  std::vector<TH1D> TauMassRefitA2MassCut;
  std::vector<TH1D> TauMassRefitA1HalfMassCut;
  std::vector<TH1D> TauMassRefitA2HalfMassCut;
  std::vector<TH1D> TauMassRefitA1FullEtaVetoCut;
  std::vector<TH1D> TauMassRefitA2FullEtaVetoCut;
  std::vector<TH1D> TauMassA1;
  std::vector<TH1D> TauMassRefitB1;
  std::vector<TH1D> TauMassRefitB1MassCut;
  std::vector<TH1D> TauMassRefitB2MassCut;
  std::vector<TH1D> TauMassRefitB1HalfMassCut;
  std::vector<TH1D> TauMassRefitB2HalfMassCut;

  std::vector<TH1D> TauMassRefitB1FullEtaVetoCut;
  std::vector<TH1D> TauMassRefitB2FullEtaVetoCut;
  std::vector<TH1D> TauMassB1;
  std::vector<TH1D> TauMassRefitC1;
  std::vector<TH1D> TauMassRefitC1MassCut;
  std::vector<TH1D> TauMassRefitC2MassCut;

  std::vector<TH1D> TauMassRefitC1HalfMassCut;
  std::vector<TH1D> TauMassRefitC2HalfMassCut;

  std::vector<TH1D> TauMassRefitC1FullEtaVetoCut;
  std::vector<TH1D> TauMassRefitC2FullEtaVetoCut;
  std::vector<TH1D> TauMassC1;

  std::vector<TH1D> TauMassRefitA2;
  std::vector<TH1D> TauMassA2;
  std::vector<TH1D> TauMassRefitB2;
  std::vector<TH1D> TauMassB2;
  std::vector<TH1D> TauMassRefitC2;
  std::vector<TH1D> TauMassC2;



  std::vector<TH1D> TauMassRefitABC1;
  std::vector<TH1D> TauMassRefitABC2;

  std::vector<TH1D> dRNearestPair;
  std::vector<TH1D> dRFarestPair;
  std::vector<TH1D> CloseDRMuMuMass;
  std::vector<TH1D> FarDRMuMuMass;


  std::vector<TH1D> SV_Mass_postselection;
  std::vector<TH1D> SV_Mass_preselection;

  std::vector<TH1D> Mu1TrackMass;
  std::vector<TH1D> Mu2TrackMass;
  std::vector<TH1D> Mu3TrackMass;


  std::vector<TH1D> Mu1AllTrackMass;
  std::vector<TH1D> Mu2AllTrackMass;
  std::vector<TH1D> Mu3AllTrackMass;




  std::vector<TH1D>   SV_Mass_postselection_disp_1mu;
  std::vector<TH1D>   SV_Mass_preselection_disp_1mu;



  std::vector<TH1D> TauMassRefitBarrel1;
  std::vector<TH1D> TauMassBarrel1;

  std::vector<TH1D> TauMassRefitEndcap1;
  std::vector<TH1D> TauMassEndcap1;

  std::vector<TH1D> TauMassRefitBarrel2;
  std::vector<TH1D> TauMassBarrel2;

  std::vector<TH1D> TauMassRefitEndcap2;
  std::vector<TH1D> TauMassEndcap2;

  std::vector<TH1D> TauMassResolutionRefit;

  std::vector<TH2D> TauMass_all_nophiVeto;
  std::vector<TH1D> TauMass_all;
  std::vector<TH2D> TauMass_allVsBDTA;
  std::vector<TH2D> TauMass_allVsBDTB;
  std::vector<TH2D> TauMass_allVsBDTC;
  std::vector<TH2D> TauMass_allVsBDTBarrel;
  std::vector<TH2D> TauMass_allVsBDTEndcap;
  std::vector<TH2D> EMR_tau_eta;
  std::vector<TH2D> L1Triggers;


  std::vector<TH2D> L1TriggersB;
  std::vector<TH2D> L1TriggersC;
  std::vector<TH2D> L1TriggersD;
  std::vector<TH2D> L1TriggersE;
  std::vector<TH2D> L1TriggersF;

  std::vector<TH2D> PairMass;
  std::vector<TH2D> PairMassFinalSel;
  std::vector<TH1D> PairMass1;
  std::vector<TH1D> PairMass2;
  std::vector<TH1D> EtaMuMuGammaMass;
  std::vector<TH1D> NonEtaMuMuGammaMass;
  std::vector<TH2D> PairMassWithCut;
  std::vector<TH2D> PairMassEta;
  std::vector<TH2D> PairMassEtaPrime;
  std::vector<TH2D> CategoryOverlap;
  std::vector<TH2D> Mu3IdOverlap;
  std::vector<TH1D> IDOriginOfOSMuon;

  std::vector<TH1D> Muon1isGlob;
  std::vector<TH1D> Muon2isGlob;
  std::vector<TH1D> Muon3isGlob;

  std::vector<TH1D> Muon1isStand;
  std::vector<TH1D> Muon2isStand;
  std::vector<TH1D> Muon3isStand;


  std::vector<TH1D> Muon1isTrack;
  std::vector<TH1D> Muon2isTrack;
  std::vector<TH1D> Muon3isTrack;



  std::vector<TH1D> Muon1PtResolution;
  std::vector<TH1D> Muon2PtResolution;
  std::vector<TH1D> Muon3PtResolution;

  std::vector<TH1D> Muon1EtaResolution;
  std::vector<TH1D> Muon2EtaResolution;
  std::vector<TH1D> Muon3EtaResolution;

  std::vector<TH1D> Muon1DRToTruth;
  std::vector<TH1D> Muon2DRToTruth;
  std::vector<TH1D> Muon3DRToTruth;


  std::vector<TH1D> TriggerMatchdR1;
  std::vector<TH1D> TriggerMatchdR2;
  std::vector<TH1D> TriggerMatchdR3;

  std::vector<TH1D> dR12;
  std::vector<TH1D> dR23;
  std::vector<TH1D> dR31;
  std::vector<TH1D> dR1Tau;
  std::vector<TH1D> dR2Tau;
  std::vector<TH1D> dR3Tau;

  std::vector<TH1D> VertexChi2KF;
  std::vector<TH1D> FLSignificance;
  std::vector<TH1D> BDTOutputA;
  std::vector<TH1D> BDTOutputB;
  std::vector<TH1D> BDTOutputC;
  std::vector<TH1D> BDTOutputBarrel;
  std::vector<TH1D> BDTOutputEndcap;
  std::vector<TH1D> NSignalCandidates;


  std::vector<TH1D>  Muon1MVAID;
  std::vector<TH1D>  Muon2MVAID;
  std::vector<TH1D>  Muon3MVAID;




  TMVA::Reader *readerA;
  TMVA::Reader *readerB;
  TMVA::Reader *readerC;

  TMVA::Reader *readerMuIDBarrel;
  TMVA::Reader *readerMuIDEndcap;



  Float_t var_vertexKFChi2;// (chi sq of the fit of the secondary vertex)
  Float_t var_svpvTauAngle;// (The angle between PV-SV vector and the tau vector)
  Float_t var_flightLenSig;// (Flight length significance of the tau candidate)
  Float_t var_sumMuTrkKinkChi2;// (sum of chi sq of the kink of all three muons)
  Float_t var_segCompMuMin;// (Minimum of the segment compatibility of the three muons)
  Float_t var_MinMIPLikelihood;// (Minimum of the calorimeter compatibility of the three muons)
  Float_t var_tauMass;
  Float_t var_MuMu_mindR;
  Float_t var_RelPt_Mu1Tau;
  Float_t var_Eta_au;
  Float_t var_MuMu_minKFChi2;
  Float_t var_maxdca;
  Float_t var_MuTau_maxdR;
  Float_t var_MaxD0Significance;
  Float_t var_IsolationMinDist;
  Float_t var_Muon1DetID;
  Float_t var_Muon2DetID;
  Float_t var_Muon3DetID;
  Float_t var_mass12_dRsorting;
  Float_t var_mass13_drSorting;


  

  float var_MaxtrkKink;
  float var_MaxD0SigSV;
  float var_MindcaTrackSV;
  float var_maxMuonsDca;
  float var_nsv;

  float var_MaxMuon_chi2LocalPosition;
  float var_MaxVertexPairQuality;
  float var_MuonglbkinkSum;
  float var_MaxMuon_chi2LocalMomentum;



  Double_t m3m;
  Double_t dataMCtype;
  Double_t event_weight;
  Double_t bdt;
  Double_t category;
  Double_t m12;
  Double_t m13;
  Double_t LumiScale;
  Double_t mvaA1;
  Double_t mvaA2;
  Double_t mvaB1;
  Double_t mvaB2;
  Double_t mvaC1;
  Double_t mvaC2;

  Float_t mu_combinedQuality_chi2LocalMomentum;
  Float_t mu_combinedQuality_chi2LocalPosition;
  Float_t mu_combinedQuality_staRelChi2;
  Float_t mu_combinedQuality_trkRelChi2;
  Float_t mu_combinedQuality_globalDeltaEtaPhi;
  Float_t mu_combinedQuality_trkKink;
  Float_t mu_combinedQuality_glbKink;
  Float_t mu_combinedQuality_glbTrackProbability;
  Float_t mu_Numberofvalidtrackerhits;
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



};
#endif
