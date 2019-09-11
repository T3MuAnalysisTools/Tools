#ifndef MyTest_h
#define MyTest_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class MyTest : public Selection {

 public:
  MyTest(TString Name_, TString id_);
  virtual ~MyTest();

  virtual void  Configure();
  virtual void  Finish();

  /////////////////////////////////////////////////////////
  // This is a cut enumerator, for other event cuts please
  // fill the enumerator with put a new enumarator with an 
  // understandanle name, for exmaple  enum cuts {TriggerOk=0,
  // PrimeVts, EventCut1, EventCut2, ..., NCuts};  
  // Do not remove/rename  the last enumerator   NCuts;

  enum cuts {TriggerOk=0,PrimeVtx,SignalCandidate,LeadingMuonPt,NCuts}; 


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


};
#endif
