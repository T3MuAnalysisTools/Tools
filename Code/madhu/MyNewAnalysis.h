#ifndef MyNewAnalysis_h
#define MyNewAnalysis_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class MyNewAnalysis : public Selection {

 public:
  MyNewAnalysis(TString Name_, TString id_);
  virtual ~MyNewAnalysis();

  virtual void  Configure();
  virtual void  Finish();

  /////////////////////////////////////////////////////////
  // This is a cut enumerator, for other event cuts please
  // fill the enumerator with put a new enumarator with an 
  // understandanle name, for exmaple  enum cuts {TriggerOk=0,
  // PrimeVts, EventCut1, EventCut2, ..., NCuts};  
  // Do not remove/rename  the last enumerator   NCuts;

  enum cuts {TriggerOk=0,PrimeVtx,SignalCandidate,LeadingMuonPt,LeadingMuonPt1,LeadingMuonPt2,NCuts};
  
  //InvMuon1G,InvMuon1T,InvMuon1S,InvMuon2G,InvMuon2T,InvMuon2S,InvMuon3G,InvMuon3T,InvMuon3S,InvTauG,InvTauT,InvTauS


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:
  // Selection Variables
  // Initializhere your analysis histograms
  std::vector<TH1D> NVtx;

  std::vector<TH1D> LeadMuonPt;
  std::vector<TH1D> LeadMuonEta;
  std::vector<TH1D> LeadMuonPhi;
  
  std::vector<TH1D> LeadMuonPt1;
  std::vector<TH1D> LeadMuonEta1;
  std::vector<TH1D> LeadMuonPhi1;
  
  std::vector<TH1D> LeadMuonPt2;
  std::vector<TH1D> LeadMuonEta2;
  std::vector<TH1D> LeadMuonPhi2;
  
  std::vector<TH1D> TauPt;
  std::vector<TH1D> TauEta;
  std::vector<TH1D> TauPhi;
  
  std::vector<TH1D> InvMu1G;
  std::vector<TH1D> InvMu1T;
  std::vector<TH1D> InvMu1S;
  
  std::vector<TH1D> InvMu2G;
  std::vector<TH1D> InvMu2T;
  std::vector<TH1D> InvMu2S;
  
  std::vector<TH1D> InvMu3G;
  std::vector<TH1D> InvMu3T;
  std::vector<TH1D> InvMu3S;
  
  std::vector<TH1D> InvTG;
  std::vector<TH1D> InvTT;
  std::vector<TH1D> InvTS;


};
#endif
