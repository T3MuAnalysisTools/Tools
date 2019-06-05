#ifndef softMuIdCuts_h
#define softMuIdCuts_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class softMuIdCuts : public Selection {

 public:
  softMuIdCuts(TString Name_, TString id_);
  virtual ~softMuIdCuts();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,SignalCandidate, Mu1PtCut, Mu2PtCut, Mu3PtCut, MuonID, PhiVeto, OmegaVeto, TriggerMatch, ThreeMuMass,NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:
  double tauMinMass_, tauMaxMass_;
  double tauMinSideBand_,tauMaxSideBand_;


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
  std::vector<TH1D> TauMass;
  std::vector<TH1D> TauMassResolution;

  std::vector<TH1D> TauMassRefit;
  std::vector<TH1D> TauMassResolutionRefit;

  std::vector<TH1D> Muon1StandardSelector;
  std::vector<TH1D> Muon2StandardSelector;
  std::vector<TH1D> Muon3StandardSelector;

  std::vector<TH1D> Muon1isGlob;
  std::vector<TH1D> Muon2isGlob;
  std::vector<TH1D> Muon3isGlob;

  std::vector<TH1D> Muon1isTrack;
  std::vector<TH1D> Muon2isTrack;
  std::vector<TH1D> Muon3isTrack;

  std::vector<TH1D> Muon1kink;
  std::vector<TH1D> Muon2kink;
  std::vector<TH1D> Muon3kink;
  std::vector<TH1D> Muon1InOutTrackMatch;
  std::vector<TH1D> Muon2InOutTrackMatch;
  std::vector<TH1D> Muon3InOutTrackMatch;

  std::vector<TH1D> Muon1PtResolution;
  std::vector<TH1D> Muon2PtResolution;
  std::vector<TH1D> Muon3PtResolution;

  std::vector<TH1D> Muon1EtaResolution;
  std::vector<TH1D> Muon2EtaResolution;
  std::vector<TH1D> Muon3EtaResolution;

  std::vector<TH1D> Muon1DRToTruth;
  std::vector<TH1D> Muon2DRToTruth;
  std::vector<TH1D> Muon3DRToTruth;

  std::vector<TH1D> MuPair1_vertex_chi2;
  std::vector<TH1D> MuPair2_vertex_chi2;
  std::vector<TH1D> MuPair3_vertex_chi2;


  std::vector<TH1D> Pair1Mass;
  std::vector<TH1D> Pair2Mass;
  std::vector<TH1D> Pair3Mass;

  std::vector<TH1D> TriggerMatchdR1;
  std::vector<TH1D> TriggerMatchdR2;
  std::vector<TH1D> TriggerMatchdR3;

  std::vector<TH1D> dR12;
  std::vector<TH1D> dR23;
  std::vector<TH1D> dR31;
  std::vector<TH1D> dR1Tau;
  std::vector<TH1D> dR2Tau;
  std::vector<TH1D> dR3Tau;

  std::vector<TH1D> Isolation_NTracks;
  std::vector<TH1D> Isolation_RelPt;
  std::vector<TH1D> Isolation_MinDist;
  std::vector<TH1D> Isolation05_RelPt;
  std::vector<TH1D> Isolation05_NTracks;
  std::vector<TH1D> Isolation05_MinDist;
  std::vector<TH1D> Isolation_Ntrk1;
  std::vector<TH1D> Isolation_Ntrk2;
  std::vector<TH1D> Isolation_Ntrk3;
  std::vector<TH1D> Isolation_Ntrk0p1;
  std::vector<TH1D> Isolation_Ntrk0p2;
  std::vector<TH1D> Isolation_Ntrk0p5;
  std::vector<TH1D> Isolation_maxdxy;
  
  std::vector<TH1D> VertexChi2KF;
  std::vector<TH1D> VertexChi2AF;
  std::vector<TH1D> VertexDCA12;
  std::vector<TH1D> VertexDCA23;
  std::vector<TH1D> VertexDCA31;
  std::vector<TH1D> VertexDCAMax;
  
  std::vector<TH1D> VertexSignalKFRefittedMu1P;
  std::vector<TH1D> VertexSignalKFRefittedMu1Pt;
  std::vector<TH1D> VertexSignalKFRefittedMu1Eta;
  std::vector<TH1D> VertexSignalKFRefittedMu1Phi;
  std::vector<TH1D> VertexSignalKFRefittedMu2P;
  std::vector<TH1D> VertexSignalKFRefittedMu2Pt;
  std::vector<TH1D> VertexSignalKFRefittedMu2Eta;
  std::vector<TH1D> VertexSignalKFRefittedMu2Phi;
  
  std::vector<TH1D> VertexSignalKFRefittedMu3P;
  std::vector<TH1D> VertexSignalKFRefittedMu3Pt;
  std::vector<TH1D> VertexSignalKFRefittedMu3Eta;
  std::vector<TH1D> VertexSignalKFRefittedMu3Phi;

  std::vector<TH1D> VertexMu1D0Reco;
  std::vector<TH1D> VertexMu1D0SigReco;
  std::vector<TH1D> VertexMu2D0Reco;
  std::vector<TH1D> VertexMu2D0SigReco;
  std::vector<TH1D> VertexMu3D0Reco;
  std::vector<TH1D> VertexMu3D0SigReco;
  std::vector<TH1D> Vertex2DDisplacement;
  std::vector<TH1D> Vertex3DDisplacement;
  std::vector<TH1D> VertexPairQuality;
  std::vector<TH1D> VertexPairRefitStatus;
  std::vector<TH1D> VertexPairfitStatus;
  std::vector<TH1D> VertexSignalKFChi2;
  std::vector<TH1D> VertexSignalAFX;
  std::vector<TH1D> VertexSignalAFY;
  std::vector<TH1D> VertexSignalAFZ;
  std::vector<TH1D> VertexSignalAFChi2;
  std::vector<TH1D> VertexSignalAFNdf;
  std::vector<TH1D> VertexMatchedPrimaryVertexX;
  std::vector<TH1D> VertexMatchedPrimaryVertexY;
  std::vector<TH1D> VertexMatchedPrimaryVertexZ;
  std::vector<TH1D> VertexRefitPVisValid;
  std::vector<TH1D> VertexMatchedRefitPrimaryVertexX;
  std::vector<TH1D> VertexMatchedRefitPrimaryVertexY;
  std::vector<TH1D> VertexMatchedRefitPrimaryVertexZ;
  
  std::vector<TH1D> Isolation_Mu1RelPt;
  std::vector<TH1D> Isolation_Mu2RelPt;
  std::vector<TH1D> Isolation_Mu3RelPt;

  std::vector<TH1D> FlightLengthSig;

};
#endif
