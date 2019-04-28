#ifndef AnalysisWithTMVA_h
#define AnalysisWithTMVA_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class AnalysisWithTMVA : public Selection {

 public:
  AnalysisWithTMVA(TString Name_, TString id_);
  virtual ~AnalysisWithTMVA();

  virtual void  Configure();
  virtual void  Finish();



  enum cuts {TriggerOk=0,PrimeVtx,NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();


 private:
  // Selection Variables
  // Initializhere your analysis histograms
  std::vector<TH1D> NVtx;
  std::vector<TH1D> MuonsPt;
  std::vector<TH1D> MuonsEta;
  std::vector<TH1D> MuonsPhi;
  std::vector<TH1D> PhiMass;
  std::vector<TH1D> TripleMass;
  std::vector<TH2D> PhiMassVsDsMass;
  std::vector<TH1D>  MuonsPtRatio;

};
#endif
