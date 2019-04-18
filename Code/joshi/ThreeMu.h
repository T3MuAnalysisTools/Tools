#ifndef ThreeMu_h
#define ThreeMu_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class ThreeMu : public Selection {

 public:
  ThreeMu(TString Name_, TString id_);
  virtual ~ThreeMu();

  virtual double deltaR(double, double, double, double);
  virtual void  Configure();
  virtual void  Finish();

  /////////////////////////////////////////////////////////
  // This is a cut enumerator, for other event cuts please
  // fill the enumerator with put a new enumarator with an 
  // understandanle name, for exmaple  enum cuts {TriggerOk=0,
  // PrimeVts, EventCut1, EventCut2, ..., NCuts};  
  // Do not remove/rename  the last enumerator   NCuts;

  enum cuts {L1SeedOk=0,HLTOk,PrimeVtx,isThreeMu,fitVtxChiSq,NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:

 // Selection Variables
  // Initializhere your analysis histograms
  std::vector<TH1D> NVtx;
  
  std::vector<TH1D> Muon1Muon2dR;
  std::vector<TH1D> Muon2Muon3dR;
  std::vector<TH1D> Muon1Muon3dR;
  std::vector<TH1D> TripleMass;
  std::vector<TH1D> Category;

  //3 Muon variables
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
