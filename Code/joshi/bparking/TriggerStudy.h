#ifndef TriggerStudy_h
#define TriggerStudy_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class TriggerStudy : public Selection {

 public:
  TriggerStudy(TString Name_, TString id_);
  virtual ~TriggerStudy();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0, SignalCandidate, NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:
  double tauMinMass_, tauMaxMass_;
  double tauMinSideBand_,tauMaxSideBand_;
  int n_doubleMu = 0;
  int n_singleMu = 0;
  int n_overlap = 0;

  // Selection Variables
  // Initializhere your analysis histograms

};
#endif
