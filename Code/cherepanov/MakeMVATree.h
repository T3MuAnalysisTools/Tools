#ifndef MakeMVATree_h
#define MakeMVATree_h

#include "Selection.h"
#include <vector>
#include "TString.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"




class MakeMVATree : public Selection {

 public:
  MakeMVATree(TString Name_, TString id_);
  virtual ~MakeMVATree();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {L1=0,HLT,SignalCandidate, Mu1PtCut, Mu2PtCut, Mu3PtCut, MuonID, TRKLWithM, PhiVeto1, OmegaVeto1, PhiVeto2, OmegaVeto2,TriggerMatch, ThreeMuMass, CutCategory,NCuts};

 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();
  TFile * file;
  TTree * TMVA_Tree;

 private:
  // Selection Variables
  double tauMinMass_, tauMaxMass_;
  double tauMinSideBand_,tauMaxSideBand_;
  double tauMassResCutLow, tauMassResCutHigh;

  // random number generator
  TRandom rndm;
  float random_num;
  int l1FailedRandom;
  int eventNumber;

   
  // PU Weights
  TFile* PUWeightFile;
  TH1D* puWeights;

  TMVA::Reader *readerMuIDBarrel;
  TMVA::Reader *readerMuIDEndcap;
  TMVA::Reader *readerBvsD;

  // Initializhere your analysis histograms


  std::vector<TH1D> pTMu1OverMass_TRF;
  std::vector<TH1D> cTheta_MuonOS_TauPol_TRF;
  std::vector<TH1D> OSSS1Angle_TRF;
  std::vector<TH1D> OSSS2Angle_TRF;
  std::vector<TH1D> OSSS1Angle_RRF;
  std::vector<TH1D> OSSS2Angle_RRF;
  std::vector<TH1D> cTheta_TRF_SSSS;
  std::vector<TH1D> cTheta_TRF_OSSS;

  std::vector<TH1D> nPVx;
  std::vector<TH1D> MuonsPt;
  std::vector<TH1D> MuonsEta;
  std::vector<TH1D> MuonsPhi;
  std::vector<TH1D> TripleMass;

  std::vector<TH1D> SVPVTauDirAngle;
  std::vector<TH1D> FLSignificance;
  std::vector<TH1D> FL;
  std::vector<TH1D> VertexChi2KF;
  std::vector<TH1D> Vertex2muTrkKF;
  std::vector<TH1D> Dist2muTrkKF3Mu;
  std::vector<TH1D> VertexQualitySeparator;
  std::vector<TH1D> Muon_segmentCompatibility_mu1;
  std::vector<TH1D> Muon_segmentCompatibility_mu2;
  std::vector<TH1D> Muon_segmentCompatibility_mu3;
  std::vector<TH1D> Muon_segmentCompatibility_min;
  std::vector<TH1D> Muon_segmentCompatibility_max;
  std::vector<TH1D> Muon_ECALCompatibility_mu1;
  std::vector<TH1D> Muon_ECALCompatibility_mu2;
  std::vector<TH1D> Muon_ECALCompatibility_mu3;
  std::vector<TH1D> Muon_ECALCompatibility_min;
  std::vector<TH1D> Muon_ECALCompatibility_max;
  std::vector<TH1D> Isolation_NTracks;
  std::vector<TH1D> Isolation_RelPt;
  std::vector<TH1D> Isolation_maxdxy;

  std::vector<TH1D> VertexDCAMax;
 


  std::vector<TH1D> Vertex2muTrkKFChi2;
  std::vector<TH1D> Vertex2muTrkKFToSignalVertexChi2;
  std::vector<TH1D> Vertex2muTrkKFToSignalVertexDistance;
  




  std::vector<TH2D> EMR_tau_eta;
  std::vector<TH1D> EventMassResolution_PtEtaPhi;

  std::vector<TH1D> deltaMuZ12;
  std::vector<TH1D> deltaMuZ13;
  std::vector<TH1D> deltaMuZ23;
  std::vector<TH1D> MaxdeltaMuZ;
  std::vector<TH1D> MindeltaMuZ;


  std::vector<TH1D> MaxMuonsDca;
  std::vector<TH1D> MinMuonsDca;

  std::vector<TH1D> NSV;
  std::vector<TH1D> SVDeltaR;
  std::vector<TH1D> SVDistance;
  std::vector<TH1D> SV_Mass;
  std::vector<TH1D> NtracksClose;


  std::vector<TH1D>   NTracksCloseToPV;
  std::vector<TH1D>   NTracksCloseToPVTauDR;




  std::vector<TH1D> VertexMu1D0SigPVReco;
  std::vector<TH1D> VertexMu2D0SigPVReco;
  std::vector<TH1D> VertexMu3D0SigPVReco;
  std::vector<TH1D> MaxD0SigPV;
  std::vector<TH1D> MinD0SigPV;


  std::vector<TH1D> VertexMu1D0SigBSReco;
  std::vector<TH1D> VertexMu2D0SigBSReco;
  std::vector<TH1D> VertexMu3D0SigBSReco;
  std::vector<TH1D> MaxD0SigBS;
  std::vector<TH1D> MinD0SigBS;


  std::vector<TH1D> VertexMu1D0SigSVReco;
  std::vector<TH1D> VertexMu2D0SigSVReco;
  std::vector<TH1D> VertexMu3D0SigSVReco;
  std::vector<TH1D> MaxD0SigSV;
  std::vector<TH1D> MinD0SigSV;


  std::vector<TH1D> VertexMu1D0DistSVReco;
  std::vector<TH1D> VertexMu2D0DistSVReco;
  std::vector<TH1D> VertexMu3D0DistSVReco;
  std::vector<TH1D> MaxD0DistSV;
  std::vector<TH1D> MinD0DistSV;
  std::vector<TH1D> VertexMu1DZDistSVReco;
  std::vector<TH1D> VertexMu2DZDistSVReco;
  std::vector<TH1D> VertexMu3DZDistSVReco;


  std::vector<TH1D> MaxDZDistSV;
  std::vector<TH1D> MinDZDistSV;

  std::vector<TH1D> VertexMu1DistanceToSV;
  std::vector<TH1D> VertexMu2DistanceToSV;
  std::vector<TH1D> VertexMu3DistanceToSV;

  std::vector<TH1D> MaxMuDistanceToSV;
  std::vector<TH1D> MinMuDistanceToSV;


  std::vector<TH1D> IsoPhiKKMass_Mu3;
  std::vector<TH1D> IsoPhiKKMass_Mu2;
  std::vector<TH1D> IsoPhiKKMass_Mu1;


  std::vector<TH1D> IsoKStarMass_Mu3;
  std::vector<TH1D> IsoKStarMass_Mu2;
  std::vector<TH1D> IsoKStarMass_Mu1;


  std::vector<TH1D> IsoMuMuMass_Mu3_wideRange;
  std::vector<TH1D> IsoMuMuMass_Mu2_wideRange;
  std::vector<TH1D> IsoMuMuMass_Mu1_wideRange;

  std::vector<TH1D> IsoPhiKKMass_Mu3_wideRange;
  std::vector<TH1D> IsoPhiKKMass_Mu2_wideRange;
  std::vector<TH1D> IsoPhiKKMass_Mu1_wideRange;

  std::vector<TH1D> IsoPhiKKMass_Mu3_midRange;
  std::vector<TH1D> IsoPhiKKMass_Mu2_midRange;
  std::vector<TH1D> IsoPhiKKMass_Mu1_midRange;


  std::vector<TH1D> IsoKStarMass_Mu3_wideRange;
  std::vector<TH1D> IsoKStarMass_Mu2_wideRange;
  std::vector<TH1D> IsoKStarMass_Mu1_wideRange;


  std::vector<TH1D> IsoMuMuMass_Mu3;
  std::vector<TH1D> IsoMuMuMass_Mu2;
  std::vector<TH1D> IsoMuMuMass_Mu1;




  std::vector<TH1D> IsolationCombinatorialMass_pipi;


  std::vector<TH1D> MinMuon_chi2LocalMomentum;
  std::vector<TH1D> MaxMuon_chi2LocalMomentum;

  std::vector<TH1D> Muon1_chi2LocalMomentum;
  std::vector<TH1D> Muon2_chi2LocalMomentum;
  std::vector<TH1D> Muon3_chi2LocalMomentum;


  std::vector<TH1D> MaxMuon_chi2LocalPosition;
  std::vector<TH1D> MinMuon_chi2LocalPosition;

  std::vector<TH1D> Muon1_chi2LocalPosition;
  std::vector<TH1D> Muon2_chi2LocalPosition;
  std::vector<TH1D> Muon3_chi2LocalPosition;

  std::vector<TH1D> MintrkKink;
  std::vector<TH1D> MaxtrkKink;
  std::vector<TH1D> MinglbKink;
  std::vector<TH1D> MaxglbKink;

  std::vector<TH1D> MuonglbkinkSum;

  std::vector<TH1D> MaxVertexPairQuality;
  std::vector<TH1D> MinVertexPairQuality;
	

  std::vector<TH1D> Muon1Pt;
  std::vector<TH1D> Muon2Pt;
  std::vector<TH1D> Muon3Pt;

  std::vector<TH1D> Muon1P;
  std::vector<TH1D> Muon2P;
  std::vector<TH1D> Muon3P;

  std::vector<TH1D> TauEta;
  std::vector<TH1D> TauPt;
  std::vector<TH1D> TauP;


  std::vector<TH1D> Mu1TrackMass;
  std::vector<TH1D> Mu2TrackMass;
  std::vector<TH1D> Mu3TrackMass;

  std::vector<TH1D> MindcaTrackSV;
  std::vector<TH1D> dcaTrackPV;
  std::vector<TH1D> Iso02;
  std::vector<TH1D> Iso04;
  std::vector<TH1D> Iso06;
  std::vector<TH1D> Iso08;
  std::vector<TH1D> Iso1;
  std::vector<TH1D> Iso12;
  std::vector<TH1D> Iso14;
  std::vector<TH1D> Iso16;
  std::vector<TH1D> Iso18;
  std::vector<TH1D> Iso2;


  std::vector<TH1D> Iso02Mu1;
  std::vector<TH1D> Iso04Mu1;
  std::vector<TH1D> Iso06Mu1;
  std::vector<TH1D> Iso08Mu1;
  std::vector<TH1D> Iso1Mu1;
  std::vector<TH1D> Iso12Mu1;
  std::vector<TH1D> Iso14Mu1;
  std::vector<TH1D> Iso16Mu1;
  std::vector<TH1D> Iso18Mu1;
  std::vector<TH1D> Iso2Mu1;

  std::vector<TH1D> Iso02Mu2;
  std::vector<TH1D> Iso04Mu2;
  std::vector<TH1D> Iso06Mu2;
  std::vector<TH1D> Iso08Mu2;
  std::vector<TH1D> Iso1Mu2;
  std::vector<TH1D> Iso12Mu2;
  std::vector<TH1D> Iso14Mu2;
  std::vector<TH1D> Iso16Mu2;
  std::vector<TH1D> Iso18Mu2;
  std::vector<TH1D> Iso2Mu2;

  std::vector<TH1D> Iso02Mu3;
  std::vector<TH1D> Iso04Mu3;
  std::vector<TH1D> Iso06Mu3;
  std::vector<TH1D> Iso08Mu3;
  std::vector<TH1D> Iso1Mu3;
  std::vector<TH1D> Iso12Mu3;
  std::vector<TH1D> Iso14Mu3;
  std::vector<TH1D> Iso16Mu3;
  std::vector<TH1D> Iso18Mu3;
  std::vector<TH1D> Iso2Mu3;

  std::vector<TH1D> Iso08MuMax;
  std::vector<TH1D> Iso08MuMin;

  std::vector<TH1D> MinMatchedStations;
  std::vector<TH1D> MaxMatchedStations;
  std::vector<TH1D> Mu1MatchedStations;
  std::vector<TH1D> Mu2MatchedStations;
  std::vector<TH1D> Mu3MatchedStations;


  std::vector<TH1D> MinMuon_numberOfChambers;
  std::vector<TH1D> MaxMuon_numberOfChambers;
  std::vector<TH1D> Mu1Muon_numberOfChambers;
  std::vector<TH1D> Mu2Muon_numberOfChambers;
  std::vector<TH1D> Mu3Muon_numberOfChambers;
	  


  std::vector<TH1D> MinMuon_numberOfMatches;
  std::vector<TH1D> MaxMuon_numberOfMatches;
  std::vector<TH1D> Mu1Muon_numberOfMatches;
  std::vector<TH1D> Mu2Muon_numberOfMatches;
  std::vector<TH1D> Mu3Muon_numberOfMatches;
  
  std::vector<TH1D> MinMuon_hitPattern_numberOfValidMuonHits;
  std::vector<TH1D> MaxMuon_hitPattern_numberOfValidMuonHits;
  std::vector<TH1D> Mu1Muon_hitPattern_numberOfValidMuonHits;
  std::vector<TH1D> Mu2Muon_hitPattern_numberOfValidMuonHits;
  std::vector<TH1D> Mu3Muon_hitPattern_numberOfValidMuonHits;

  std::vector<TH1D> Separation_BvsD;
  std::vector<TH1D> NBJet4pi;
  std::vector<TH1D> BTagCSV;
  std::vector<TH1D> BTagMVA;
  std::vector<TH1D> BTagCVSB;


  std::vector<TH1D> NBJet4piSH;
  std::vector<TH1D> NBJet4piOH;

  std::vector<TH1D> BTagCSVSH;
  std::vector<TH1D> BTagMVASH;
  std::vector<TH1D> BTagCVSBSH;

  std::vector<TH1D> BTagCSVOH;
  std::vector<TH1D> BTagMVAOH;
  std::vector<TH1D> BTagCVSBOH;


  std::vector<TH2D> BTagCSVSHVsNJets;
  std::vector<TH2D> BTagMVASHVsNJets;
  std::vector<TH2D> BTagCVSBSHVsNJets;

  std::vector<TH2D> BTagCSVOHVsNJets;
  std::vector<TH2D> BTagMVAOHVsNJets;
  std::vector<TH2D> BTagCVSBOHVsNJets;

  std::vector<TH1D> BTagCSVSHMatchedToTau;
  std::vector<TH1D> BTagMVASHMatchedToTau;
  std::vector<TH1D> BTagCVSBSHMatchedToTau;

  std::vector<TH1D> BTagCSVOHMatchedToTau;
  std::vector<TH1D> BTagMVAOHMatchedToTau;
  std::vector<TH1D> BTagCVSBOHMatchedToTau;




  std::vector<TH1D> Muon1ImpactAngle;
  std::vector<TH1D> Muon2ImpactAngle;
  std::vector<TH1D> Muon3ImpactAngle;
  std::vector<TH1D> MinMuonImpacAngle;
  std::vector<TH1D> MaxMuonImpacAngle;


  std::vector<TH1D> Muon1StandardSelectorPass;
  std::vector<TH1D> Muon2StandardSelectorPass;
  std::vector<TH1D> Muon3StandardSelectorPass;

  std::vector<TH1D> Muon1StandardSelectorFail;
  std::vector<TH1D> Muon2StandardSelectorFail;
  std::vector<TH1D> Muon3StandardSelectorFail;


  std::vector<TH1D>  Muon1MVAID;
  std::vector<TH1D>  Muon2MVAID;
  std::vector<TH1D>  Muon3MVAID;



  std::vector<TH1D>  PairMass1AllignedSorting;
  std::vector<TH1D>  PairMass2AllignedSorting;
  std::vector<TH2D>  MuMuMassAllignedSorting;


 std::vector<TH1D>  PairMass1PTSorting;
 std::vector<TH1D>  PairMass2PTSorting;
 std::vector<TH2D>  MuMuMassPTSorting;


 std::vector<TH1D>  PairMass1NoSorting;
 std::vector<TH1D>  PairMass2NoSorting;
 std::vector<TH2D>  MuMuMassNoSorting;







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
  Float_t mu_Numberofvalidpixelhits;
  Float_t mu_trackerLayersWithMeasurement;
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


  //------------ Mini MVA tree variables


  bool MC;
  int category;
  float var_EventNumber;
  float var_id;
  float var_vertexKFChi2 ;
  float var_svpvTauAngle ;
  float var_flightLenSig ;
  float var_flightLenDist ;
  float var_sumMuTrkKinkChi2 ;
  float var_segCompMuMin ;
  float var_segCompMuMax ;
  float var_segCompMu1 ;
  float var_segCompMu2 ;
  float var_segCompMu3 ;
  float var_caloCompMin ;
  float var_caloCompMax ;
  float var_caloCompMu1 ;
  float var_caloCompMu2 ;
  float var_caloCompMu3 ;
  float var_MinMIPLikelihood ;
  float var_tauMass ;
  float var_tauEta ;
  float var_ntracks;
  float var_relPt;
  float var_isoMax;

  float var_MaxdeltaMuZ;
  float var_MindeltaMuZ;

  float var_deltaMuZ12;
  float	var_deltaMuZ13;
  float	var_deltaMuZ23;


  float var_maxMuonsDca;
  float var_minMuonsDca;
  float var_nsv;


  float var_VertexMu1D0SigPVReco;
  float var_VertexMu2D0SigPVReco;
  float var_VertexMu3D0SigPVReco;
	
  float var_MaxD0SigPV;
  float var_MinD0SigPV;


  float var_VertexMu1D0SigBSReco;
  float var_VertexMu2D0SigBSReco;
  float var_VertexMu3D0SigBSReco;


  float var_MaxD0SigBS;
  float var_MinD0SigBS;
	

  float var_VertexMu1D0SigSVReco;
  float var_VertexMu2D0SigSVReco;
  float var_VertexMu3D0SigSVReco;


  float var_MaxD0SigSV;
  float var_MinD0SigSV;
	
  float var_VertexMu1DistanceToSV;
  float var_VertexMu2DistanceToSV;
  float var_VertexMu3DistanceToSV;

  float var_MaxMuDistanceToSV;
  float var_MinMuDistanceToSV;



  float var_MinMuon_chi2LocalPosition;
  float var_MaxMuon_chi2LocalPosition;

  float var_Muon1_chi2LocalPosition;
  float var_Muon2_chi2LocalPosition;
  float var_Muon3_chi2LocalPosition;
	
	
  float var_MinMuon_chi2LocalMomentum;
  float var_MaxMuon_chi2LocalMomentum;

  float var_Muon1_chi2LocalMomentum;
  float var_Muon2_chi2LocalMomentum;
  float var_Muon3_chi2LocalMomentum;
	
  float var_MintrkKink;
  float var_MaxtrkKink;

  float var_Muon1trkKink;
  float var_Muon2trkKink;
  float var_Muon3trkKink;

  float var_MinglbKink;
  float var_MaxglbKink;
	
  float var_MuonglbkinkSum;

  float var_MaxVertexPairQuality;
  float var_MinVertexPairQuality;
	


  float var_Iso02;
  float var_Iso04;
  float var_Iso06;
  float var_Iso08;
  float var_Iso1;
  float var_Iso12;

  float var_Iso02Mu1;
  float var_Iso04Mu1;
  float var_Iso06Mu1;
  float var_Iso08Mu1;
  float var_Iso1Mu1;
  float var_Iso12Mu1;



  float var_Iso02Mu2;
  float var_Iso04Mu2;
  float var_Iso06Mu2;
  float var_Iso08Mu2;
  float var_Iso1Mu2;
  float var_Iso12Mu2;

  float var_Iso02Mu3;
  float var_Iso04Mu3;
  float var_Iso06Mu3;
  float var_Iso08Mu3;
  float var_Iso1Mu3;
  float var_Iso12Mu3;

  float var_Iso08MuMax;
  float var_Iso08MuMin;

  float var_NtracksClose;
  float var_NTracksCloseToPV;
  float var_NTracksCloseToPVTauDR;


  float var_Muon1Pt;
  float var_Muon2Pt;
  float var_Muon3Pt;

  float var_Muon1P;
  float var_Muon2P;
  float var_Muon3P;
  float var_MindcaTrackSV;

  float var_Mu1TrackMass;
  float var_Mu2TrackMass;
  float var_Mu3TrackMass;
 
  float var_MinMatchedStations;
  float var_MaxMatchedStations;
  float var_Mu1MatchedStations;
  float var_Mu2MatchedStations;
  float var_Mu3MatchedStations;


  float var_MinMuon_numberOfChambers;
  float var_MaxMuon_numberOfChambers;
  float var_Mu1Muon_numberOfChambers;
  float var_Mu2Muon_numberOfChambers;
  float var_Mu3Muon_numberOfChambers;

  float  var_Muon1ImpactAngle;
  float  var_Muon2ImpactAngle;
  float  var_Muon3ImpactAngle;
  float  var_MinMuonImpactAngle;
  float  var_MaxMuonImpactAngle;


  float var_MinMuon_numberOfMatches;
  float var_MaxMuon_numberOfMatches;
  float var_Mu1Muon_numberOfMatches;
  float var_Mu2Muon_numberOfMatches;
  float var_Mu3Muon_numberOfMatches;
	
  float var_MinMuon_hitPattern_numberOfValidMuonHits;
  float var_MaxMuon_hitPattern_numberOfValidMuonHits;
  float var_Mu1Muon_hitPattern_numberOfValidMuonHits;
  float var_Mu2Muon_hitPattern_numberOfValidMuonHits;
  float var_Mu3Muon_hitPattern_numberOfValidMuonHits;

  float var_dcaTrackPV;

  float var_Muon1LooseId;
  float var_Muon1MediumId;
  float var_Muon1TightId;

  float var_Muon2LooseId;
  float var_Muon2MediumId;
  float var_Muon2TightId;

  float var_Muon3LooseId;
  float var_Muon3MediumId;
  float var_Muon3TightId;

  float var_Muon1PFIsoLoose;
  float var_Muon1PFIsoMedium;
  float var_Muon1PFIsoTight;
  float var_Muon1PFIsoVTight;

  float var_Muon2PFIsoLoose;
  float var_Muon2PFIsoMedium;
  float var_Muon2PFIsoTight;
  float var_Muon2PFIsoVTight;

  float var_Muon3PFIsoLoose;
  float var_Muon3PFIsoMedium;
  float var_Muon3PFIsoTight;
  float var_Muon3PFIsoVTight;

  float var_mass12;
  float var_mass13;

  float var_mass12_dRsorting;
  float var_mass13_dRsorting;



  float var_Muon1DetID;
  float var_Muon2DetID;
  float var_Muon3DetID;

  float var_VertexMu1D0DistSVReco;
  float var_VertexMu2D0DistSVReco;
  float var_VertexMu3D0DistSVReco;
  float var_MaxD0DistSV;
  float var_MinD0DistSV;

  float var_VertexMu1DZDistSVReco;
  float var_VertexMu2DZDistSVReco;
  float var_VertexMu3DZDistSVReco;
  float var_MaxDZDistSV;
  float var_MinDZDistSV;


  float var_IsoPhiKKMass_Mu1;
  float var_IsoKStarMass_Mu1;
  float var_IsoMuMuMass_Mu1;

  float var_IsoPhiKKMass_Mu2;
  float var_IsoKStarMass_Mu2;
  float var_IsoMuMuMass_Mu2;

  float var_IsoPhiKKMass_Mu3;
  float var_IsoKStarMass_Mu3;
  float var_IsoMuMuMass_Mu3;

  float var_BvsDSeprator;

  float var_Vertex2muTrkKF;
  float var_Dist2muTrkKF3Mu;
  float var_VertexQualitySeparator;

  float var_cTheta_TRF_SSSS;
  float var_cTheta_TRF_OSSS;
  float var_cTheta_MuonOS_TauPol_TRF;
  float var_pTMu1OverMass_TRF;
  float var_costheta_TRF_SSSS;
  float var_costheta_TRF_OSS;
  float var_OSSS1Angle_TRF;
  float var_OSSS2Angle_TRF;

  float var_OSSS1Angle_RRF;
  float var_OSSS2Angle_RRF;





};
#endif
