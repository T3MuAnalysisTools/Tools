#ifndef SimpleTauSelector_h
#define SimpleTauSelector_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class SimpleTauSelector : public Selection {

 public:
  SimpleTauSelector(TString Name_, TString id_);
  virtual ~SimpleTauSelector();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,
	     PrimeVtx, 
	     SignalCandidate,
	     MuonCandidate,
	     OSCharge,
	     nTaus,
	     NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:

  std::vector<TH1D>   NumberOfTaus;
  std::vector<TH1D>   Mu3MuVisibleMass;
  std::vector<TH1D>   Mu3MudPhi;
  std::vector<TH2D>   NumTausvsNumMuons;


};
#endif
