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
  float var_mu3d0VertexSig;
  float var_maxdca;

  float var_MaxdeltaMuZ;
  float var_MindeltaMuZ;
  float var_maxMuonsDca;
  float var_minMuonsDca;
  float var_nsv;

};
#endif
