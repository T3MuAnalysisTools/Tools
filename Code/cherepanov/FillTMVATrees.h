#ifndef FillTMVATrees_h
#define FillTMVATrees_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class FillTMVATrees : public Selection {

 public:
  FillTMVATrees(TString Name_, TString id_);
  virtual ~FillTMVATrees();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,SignalCandidate, Mu1PtCut, Mu2PtCut, Mu3PtCut, MuonID, PhiVeto, OmegaVeto, TriggerMatch, ThreeMuMass,NCuts};

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
  // Initializhere your analysis histograms
  std::vector<TH1D> NVtx;
  std::vector<TH1D> MuonsPt;
  std::vector<TH1D> MuonsEta;
  std::vector<TH1D> MuonsPhi;
  std::vector<TH1D> TripleMass;

  std::vector<TH1D> SVPVTauDirAngle;
  std::vector<TH1D> FLSignificance;
  std::vector<TH1D> FL;
  std::vector<TH1D> VertexChi2KF;
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


  std::vector<TH1D> MinMuon_chi2LocalMomentum;
  std::vector<TH1D> MaxMuon_chi2LocalMomentum;

  std::vector<TH1D> MaxMuon_chi2LocalPosition;
  std::vector<TH1D> MinMuon_chi2LocalPosition;

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

  //------------ Mini MVA tree variables


  bool MC;
  int category;
  float var_vertexKFChi2 ;
  float var_svpvTauAngle ;
  float var_flightLenSig ;
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
  float var_ntracks;
  float var_relPt;
  float var_isoMax;

  float var_MaxdeltaMuZ;
  float var_MindeltaMuZ;
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
	

  float var_MinMuon_chi2LocalPosition;
  float var_MaxMuon_chi2LocalPosition;
	
	
  float var_MinMuon_chi2LocalMomentum;
  float var_MaxMuon_chi2LocalMomentum;
	
	
	
  float var_MintrkKink;
  float var_MaxtrkKink;
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
  float var_Muon1Pt;
  float var_Muon2Pt;
  float var_Muon3Pt;
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
};
#endif
