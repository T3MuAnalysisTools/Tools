#ifndef TripleMuIsolation_h
#define TripleMuIsolation_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class TripleMuIsolation : public Selection {

 public:
  TripleMuIsolation(TString Name_, TString id_);
  virtual ~TripleMuIsolation();

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

  // calorimeter realted quantities
  std::vector<TH1D> Muon1_Isolation_emEt03;
  std::vector<TH1D> Muon1_Isolation_emVetoEt03;
  std::vector<TH1D> Muon1_Isolation_hadEt03;
  std::vector<TH1D> Muon1_Isolation_hadVetoEt03;
  std::vector<TH1D> Muon1_Isolation_nJets03;
  std::vector<TH1D> Muon1_Isolation_nTracks03;
  std::vector<TH1D> Muon1_Isolation_sumPt03;
  std::vector<TH1D> Muon1_Isolation_trackerVetoPt03;
  std::vector<TH1D> Muon1_Isolation_emEt05;
  std::vector<TH1D> Muon1_Isolation_emVetoEt05;
  std::vector<TH1D> Muon1_Isolation_hadEt05;
  std::vector<TH1D> Muon1_Isolation_hadVetoEt05;
  std::vector<TH1D> Muon1_Isolation_nJets05;
  std::vector<TH1D> Muon1_Isolation_nTracks05;
  std::vector<TH1D> Muon1_Isolation_sumPt05;
  std::vector<TH1D> Muon1_Isolation_trackerVetoPt05;
  std::vector<TH1D> Muon1_Isolation_sumChargedHadronPt03;
  std::vector<TH1D> Muon1_Isolation_sumChargedParticlePt03;
  std::vector<TH1D> Muon1_Isolation_sumNeutralHadronEt03;
  std::vector<TH1D> Muon1_Isolation_sumNeutralHadronEtHighThreshold03;
  std::vector<TH1D> Muon1_Isolation_sumPhotonEt03;
  std::vector<TH1D> Muon1_Isolation_sumPhotonEtHighThreshold03;
  std::vector<TH1D> Muon1_Isolation_sumPUPt03;
  std::vector<TH1D> Muon1_Isolation_sumChargedHadronPt04;
  std::vector<TH1D> Muon1_Isolation_sumChargedParticlePt04;
  std::vector<TH1D> Muon1_Isolation_sumNeutralHadronEt04;
  std::vector<TH1D> Muon1_Isolation_sumNeutralHadronEtHighThreshold04;
  std::vector<TH1D> Muon1_Isolation_sumPhotonEt04;
  std::vector<TH1D> Muon1_Isolation_sumPhotonEtHighThreshold04;
  std::vector<TH1D> Muon1_Isolation_sumPUPt04;

  std::vector<TH1D> Muon2_Isolation_emEt03;
  std::vector<TH1D> Muon2_Isolation_emVetoEt03;
  std::vector<TH1D> Muon2_Isolation_hadEt03;
  std::vector<TH1D> Muon2_Isolation_hadVetoEt03;
  std::vector<TH1D> Muon2_Isolation_nJets03;
  std::vector<TH1D> Muon2_Isolation_nTracks03;
  std::vector<TH1D> Muon2_Isolation_sumPt03;
  std::vector<TH1D> Muon2_Isolation_trackerVetoPt03;
  std::vector<TH1D> Muon2_Isolation_emEt05;
  std::vector<TH1D> Muon2_Isolation_emVetoEt05;
  std::vector<TH1D> Muon2_Isolation_hadEt05;
  std::vector<TH1D> Muon2_Isolation_hadVetoEt05;
  std::vector<TH1D> Muon2_Isolation_nJets05;
  std::vector<TH1D> Muon2_Isolation_nTracks05;
  std::vector<TH1D> Muon2_Isolation_sumPt05;
  std::vector<TH1D> Muon2_Isolation_trackerVetoPt05;
  std::vector<TH1D> Muon2_Isolation_sumChargedHadronPt03;
  std::vector<TH1D> Muon2_Isolation_sumChargedParticlePt03;
  std::vector<TH1D> Muon2_Isolation_sumNeutralHadronEt03;
  std::vector<TH1D> Muon2_Isolation_sumNeutralHadronEtHighThreshold03;
  std::vector<TH1D> Muon2_Isolation_sumPhotonEt03;
  std::vector<TH1D> Muon2_Isolation_sumPhotonEtHighThreshold03;
  std::vector<TH1D> Muon2_Isolation_sumPUPt03;
  std::vector<TH1D> Muon2_Isolation_sumChargedHadronPt04;
  std::vector<TH1D> Muon2_Isolation_sumChargedParticlePt04;
  std::vector<TH1D> Muon2_Isolation_sumNeutralHadronEt04;
  std::vector<TH1D> Muon2_Isolation_sumNeutralHadronEtHighThreshold04;
  std::vector<TH1D> Muon2_Isolation_sumPhotonEt04;
  std::vector<TH1D> Muon2_Isolation_sumPhotonEtHighThreshold04;
  std::vector<TH1D> Muon2_Isolation_sumPUPt04;

  std::vector<TH1D> Muon3_Isolation_emEt03;
  std::vector<TH1D> Muon3_Isolation_emVetoEt03;
  std::vector<TH1D> Muon3_Isolation_hadEt03;
  std::vector<TH1D> Muon3_Isolation_hadVetoEt03;
  std::vector<TH1D> Muon3_Isolation_nJets03;
  std::vector<TH1D> Muon3_Isolation_nTracks03;
  std::vector<TH1D> Muon3_Isolation_sumPt03;
  std::vector<TH1D> Muon3_Isolation_trackerVetoPt03;
  std::vector<TH1D> Muon3_Isolation_emEt05;
  std::vector<TH1D> Muon3_Isolation_emVetoEt05;
  std::vector<TH1D> Muon3_Isolation_hadEt05;
  std::vector<TH1D> Muon3_Isolation_hadVetoEt05;
  std::vector<TH1D> Muon3_Isolation_nJets05;
  std::vector<TH1D> Muon3_Isolation_nTracks05;
  std::vector<TH1D> Muon3_Isolation_sumPt05;
  std::vector<TH1D> Muon3_Isolation_trackerVetoPt05;
  std::vector<TH1D> Muon3_Isolation_sumChargedHadronPt03;
  std::vector<TH1D> Muon3_Isolation_sumChargedParticlePt03;
  std::vector<TH1D> Muon3_Isolation_sumNeutralHadronEt03;
  std::vector<TH1D> Muon3_Isolation_sumNeutralHadronEtHighThreshold03;
  std::vector<TH1D> Muon3_Isolation_sumPhotonEt03;
  std::vector<TH1D> Muon3_Isolation_sumPhotonEtHighThreshold03;
  std::vector<TH1D> Muon3_Isolation_sumPUPt03;
  std::vector<TH1D> Muon3_Isolation_sumChargedHadronPt04;
  std::vector<TH1D> Muon3_Isolation_sumChargedParticlePt04;
  std::vector<TH1D> Muon3_Isolation_sumNeutralHadronEt04;
  std::vector<TH1D> Muon3_Isolation_sumNeutralHadronEtHighThreshold04;
  std::vector<TH1D> Muon3_Isolation_sumPhotonEt04;
  std::vector<TH1D> Muon3_Isolation_sumPhotonEtHighThreshold04;
  std::vector<TH1D> Muon3_Isolation_sumPUPt04;

  std::vector<TH1D> Min_Isolation_emEt03;
  std::vector<TH1D> Min_Isolation_emVetoEt03;
  std::vector<TH1D> Min_Isolation_hadEt03;
  std::vector<TH1D> Min_Isolation_hadVetoEt03;
  std::vector<TH1D> Min_Isolation_nJets03;
  std::vector<TH1D> Min_Isolation_nTracks03;
  std::vector<TH1D> Min_Isolation_sumPt03;
  std::vector<TH1D> Min_Isolation_trackerVetoPt03;
  std::vector<TH1D> Min_Isolation_emEt05;
  std::vector<TH1D> Min_Isolation_emVetoEt05;
  std::vector<TH1D> Min_Isolation_hadEt05;
  std::vector<TH1D> Min_Isolation_hadVetoEt05;
  std::vector<TH1D> Min_Isolation_nJets05;
  std::vector<TH1D> Min_Isolation_nTracks05;
  std::vector<TH1D> Min_Isolation_sumPt05;
  std::vector<TH1D> Min_Isolation_trackerVetoPt05;
  std::vector<TH1D> Min_Isolation_sumChargedHadronPt03;
  std::vector<TH1D> Min_Isolation_sumChargedParticlePt03;
  std::vector<TH1D> Min_Isolation_sumNeutralHadronEt03;
  std::vector<TH1D> Min_Isolation_sumNeutralHadronEtHighThreshold03;
  std::vector<TH1D> Min_Isolation_sumPhotonEt03;
  std::vector<TH1D> Min_Isolation_sumPhotonEtHighThreshold03;
  std::vector<TH1D> Min_Isolation_sumPUPt03;
  std::vector<TH1D> Min_Isolation_sumChargedHadronPt04;
  std::vector<TH1D> Min_Isolation_sumChargedParticlePt04;
  std::vector<TH1D> Min_Isolation_sumNeutralHadronEt04;
  std::vector<TH1D> Min_Isolation_sumNeutralHadronEtHighThreshold04;
  std::vector<TH1D> Min_Isolation_sumPhotonEt04;
  std::vector<TH1D> Min_Isolation_sumPhotonEtHighThreshold04;
  std::vector<TH1D> Min_Isolation_sumPUPt04;

  std::vector<TH1D> Max_Isolation_emEt03;
  std::vector<TH1D> Max_Isolation_emVetoEt03;
  std::vector<TH1D> Max_Isolation_hadEt03;
  std::vector<TH1D> Max_Isolation_hadVetoEt03;
  std::vector<TH1D> Max_Isolation_nJets03;
  std::vector<TH1D> Max_Isolation_nTracks03;
  std::vector<TH1D> Max_Isolation_sumPt03;
  std::vector<TH1D> Max_Isolation_trackerVetoPt03;
  std::vector<TH1D> Max_Isolation_emEt05;
  std::vector<TH1D> Max_Isolation_emVetoEt05;
  std::vector<TH1D> Max_Isolation_hadEt05;
  std::vector<TH1D> Max_Isolation_hadVetoEt05;
  std::vector<TH1D> Max_Isolation_nJets05;
  std::vector<TH1D> Max_Isolation_nTracks05;
  std::vector<TH1D> Max_Isolation_sumPt05;
  std::vector<TH1D> Max_Isolation_trackerVetoPt05;
  std::vector<TH1D> Max_Isolation_sumChargedHadronPt03;
  std::vector<TH1D> Max_Isolation_sumChargedParticlePt03;
  std::vector<TH1D> Max_Isolation_sumNeutralHadronEt03;
  std::vector<TH1D> Max_Isolation_sumNeutralHadronEtHighThreshold03;
  std::vector<TH1D> Max_Isolation_sumPhotonEt03;
  std::vector<TH1D> Max_Isolation_sumPhotonEtHighThreshold03;
  std::vector<TH1D> Max_Isolation_sumPUPt03;
  std::vector<TH1D> Max_Isolation_sumChargedHadronPt04;
  std::vector<TH1D> Max_Isolation_sumChargedParticlePt04;
  std::vector<TH1D> Max_Isolation_sumNeutralHadronEt04;
  std::vector<TH1D> Max_Isolation_sumNeutralHadronEtHighThreshold04;
  std::vector<TH1D> Max_Isolation_sumPhotonEt04;
  std::vector<TH1D> Max_Isolation_sumPhotonEtHighThreshold04;
  std::vector<TH1D> Max_Isolation_sumPUPt04;

};
#endif
