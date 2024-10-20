#ifndef ThreeMu_h
#define ThreeMu_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class ThreeMu : public Selection {

 public:
  ThreeMu(TString Name_, TString id_);
  virtual ~ThreeMu();

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

  std::vector<TH1D> SVPVTauDirAngle;

  std::vector<TH1D> TauMassRefit;
  std::vector<TH1D> TauMassResolutionRefit;

  std::vector<TH1D> Muon1StandardSelector;
  std::vector<TH1D> Muon2StandardSelector;
  std::vector<TH1D> Muon3StandardSelector;

  std::vector<TH1D> Muon1isGlob;
  std::vector<TH1D> Muon2isGlob;
  std::vector<TH1D> Muon3isGlob;

  std::vector<TH1D> Muon1isStand;
  std::vector<TH1D> Muon2isStand;
  std::vector<TH1D> Muon3isStand;


  std::vector<TH1D> Muon1isTrack;
  std::vector<TH1D> Muon2isTrack;
  std::vector<TH1D> Muon3isTrack;

  std::vector<TH1D> Muon1kink;
  std::vector<TH1D> Muon2kink;
  std::vector<TH1D> Muon3kink;
  std::vector<TH1D> MuonkinkMax;
  std::vector<TH1D> MuonkinkMin;
  std::vector<TH1D> MuonkinkSum;

  std::vector<TH1D> Muon1glbkink;
  std::vector<TH1D> Muon2glbkink;
  std::vector<TH1D> Muon3glbkink;
  std::vector<TH1D> MuonglbkinkMax;
  std::vector<TH1D> MuonglbkinkMin;
  std::vector<TH1D> MuonglbkinkSum;

  std::vector<TH1D> Muon1InOutTrackMatch;
  std::vector<TH1D> Muon2InOutTrackMatch;
  std::vector<TH1D> Muon3InOutTrackMatch;

  std::vector<TH1D> MuonInOutTrackMatchMax;

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
  std::vector<TH1D> MuPair_vertex_chi2_min;

  std::vector<TH1D> Pair1Mass;
  std::vector<TH1D> Pair2Mass;
  std::vector<TH1D> Pair3Mass;

  std::vector<TH1D> Pair1Mass_OS1;
  std::vector<TH1D> Pair2Mass_OS2;
  std::vector<TH1D> Pair3Mass_SS;

  std::vector<TH1D> TriggerMatchdR1;
  std::vector<TH1D> TriggerMatchdR2;
  std::vector<TH1D> TriggerMatchdR3;

  std::vector<TH1D> dR12;
  std::vector<TH1D> dR23;
  std::vector<TH1D> dR31;
  std::vector<TH1D> dR1Tau;
  std::vector<TH1D> dR2Tau;
  std::vector<TH1D> dR3Tau;

  std::vector<TH1D> Mu1TauPTRatio;
  std::vector<TH1D> Mu2TauPTRatio;
  std::vector<TH1D> Mu3TauPTRatio;

  std::vector<TH1D> Mu1TauPRatio;
  std::vector<TH1D> Mu2TauPRatio;
  std::vector<TH1D> Mu3TauPRatio;

  std::vector<TH1D> maxMudR;
  std::vector<TH1D> minMudR;

  std::vector<TH1D> dRMaxMuTau;
  std::vector<TH1D> dRMinMuTau;

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

  std::vector<TH1D> Isolation_RelativePtMu1;  
  std::vector<TH1D> Isolation_RelativePtMu2;  
  std::vector<TH1D> Isolation_RelativePtMu3;  
  std::vector<TH1D> Isolation_RelativePtMaxMu;  

  std::vector<TH1D> Isolation_Muon_hadVetoEt03;
  std::vector<TH1D> Isolation_Muon_hadEt03;
  std::vector<TH1D> Isolation_Muon_emVetoEt03;
  std::vector<TH1D> Isolation_Muon_emEt03;
  std::vector<TH1D> Isolation_Muon_nJets03;
  std::vector<TH1D> Isolation_Muon_nTracks03;
  std::vector<TH1D> Isolation_Muon_sumPt03;
  std::vector<TH1D> Isolation_Muon_trackerVetoPt03;
  std::vector<TH1D> Isolation_Muon_hadVetoEt05;
  std::vector<TH1D> Isolation_Muon_hadEt05;
  std::vector<TH1D> Isolation_Muon_emVetoEt05;
  std::vector<TH1D> Isolation_Muon_emEt05;
  std::vector<TH1D> Isolation_Muon_nJets05;
  std::vector<TH1D> Isolation_Muon_nTracks05;
  std::vector<TH1D> Isolation_Muon_sumPt05;
  std::vector<TH1D> Isolation_Muon_trackerVetoPt05;
  std::vector<TH1D> Isolation_Muon_sumChargedHadronPt03;
  std::vector<TH1D> Isolation_Muon_sumChargedParticlePt03;
  std::vector<TH1D> Isolation_Muon_sumNeutralHadronEt03;
  std::vector<TH1D> Isolation_Muon_sumPhotonEt03;
  std::vector<TH1D> Isolation_Muon_sumPUPt03;
  std::vector<TH1D> Isolation_Muon_sumNeutralHadronEtHighThreshold03;
  std::vector<TH1D> Isolation_Muon_sumPhotonEtHighThreshold03;
  std::vector<TH1D> Isolation_Muon_sumChargedHadronPt04;
  std::vector<TH1D> Isolation_Muon_sumChargedParticlePt04;
  std::vector<TH1D> Isolation_Muon_sumNeutralHadronEt04;
  std::vector<TH1D> Isolation_Muon_sumPhotonEt04;
  std::vector<TH1D> Isolation_Muon_sumPUPt04;
  std::vector<TH1D> Isolation_Muon_sumNeutralHadronEtHighThreshold04;
  std::vector<TH1D> Isolation_Muon_sumPhotonEtHighThreshold04;

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
  std::vector<TH1D> FLSignificance;
  std::vector<TH1D> FLSignificance2;

  std::vector<TH1D> VertexMu1D0Reco;
  std::vector<TH1D> VertexMu1D0SigReco;
  std::vector<TH1D> VertexMu2D0Reco;
  std::vector<TH1D> VertexMu2D0SigReco;
  std::vector<TH1D> VertexMu3D0Reco;
  std::vector<TH1D> VertexMu3D0SigReco;
  std::vector<TH1D> VertexMuMinD0SigReco;
  std::vector<TH1D> VertexMuMaxD0SigReco;

  std::vector<TH1D> Vertex2DDisplacement;
  std::vector<TH1D> Vertex3DDisplacement;
  std::vector<TH1D> Vertex2DDisplacementSignificance;
  std::vector<TH1D> Vertex3DDisplacementSignificance;
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

  std::vector<TH1D> Muon1_globalDeltaEtaPhi;
  std::vector<TH1D> Muon2_globalDeltaEtaPhi;
  std::vector<TH1D> Muon3_globalDeltaEtaPhi;



};
#endif
