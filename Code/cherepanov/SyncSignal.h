#ifndef SyncSignal_h
#define SyncSignal_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class SyncSignal : public Selection {

 public:
  SyncSignal(TString Name_, TString id_);
  virtual ~SyncSignal();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,SignalCandidate, Mu1PtCut, Mu2PtCut, Mu3PtCut, MuonID, PhiVeto, OmegaVeto, TriggerMatch, TauMassCut,NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();
  TFile * file;
  TTree * Sync_tree;


 private:
  double tauMinMass_, tauMaxMass_;
  double tauMinSideBand_,tauMaxSideBand_;


  // Selection Variables
  // Initializhere your analysis histograms

  double evt,run,lumi;

  double sync_pt_1;
  double sync_pt_2;
  double sync_pt_3;

  double sync_eta_1;
  double sync_eta_2;
  double sync_eta_3;

  double muon_1_isGlob;
  double muon_2_isGlob;
  double muon_3_isGlob;

  double muon_1_isTrack;
  double muon_2_isTrack;
  double muon_3_isTrack;

  double nSignalCandidates;


};
#endif
