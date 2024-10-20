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

  /////////////////////////////////////////////////////////
  // This is a cut enumerator, for other event cuts please
  // fill the enumerator with put a new enumarator with an 
  // understandanle name, for exmaple  enum cuts {TriggerOk=0,
  // PrimeVts, EventCut1, EventCut2, ..., NCuts};  
  // Do not remove/rename  the last enumerator   NCuts;

  enum cuts {HLTOk=0, ThreeMuCandidate, TriggerMatch, ThreeMuMass, MuID, PhiVeto, OmegaVeto, NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:
  // Selection Variables
  // Initializhere your analysis histograms

  double tauMinMass_, tauMaxMass_;

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
  std::vector<TH1D> VertexChi2KF;
  std::vector<TH1D> VertexChi2AF;

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


};
#endif
