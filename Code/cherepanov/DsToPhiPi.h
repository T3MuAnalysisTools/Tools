#ifndef DsToPhiPi_h
#define DsToPhiPi_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class DsToPhiPi : public Selection {

 public:
  DsToPhiPi(TString Name_, TString id_);
  virtual ~DsToPhiPi();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {HLTOk=0,L1Ok,is2MuTrk,GlobalMu,Chi2Cut,Mass2Mu,Mu1dR,Mu2dR,TrkdR,Mu1pt,Mu2pt,Trkpt,NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:
  // Selection Variables
  // Initializhere your analysis histograms
  double RelLumiB,RelLumiC, RelLumiD, RelLumiE, RelLumiF;
  std::vector<TH1D> NVtx;

  std::vector<TH1D> DimuondR;
  std::vector<TH1D> Muon1TrkdR;
  std::vector<TH1D> Muon2TrkdR;
  std::vector<TH1D> PhiMass;
  std::vector<TH1D> PhiPlusTrackMass;
  std::vector<TH2D> PhiMassVsDsMass;
  std::vector<TH1D> DsMass;
  std::vector<TH1D> DsMassB;
  std::vector<TH1D> DsMassC;
  std::vector<TH1D> DsMassD;
  std::vector<TH1D> DsMassE;
  std::vector<TH1D> DsMassF;

  std::vector<TH1D> Category;
  std::vector<TH1D> DsGenMatch;

  //Muon variables
  
  //Track candidate variables
  std::vector<TH1D> Track_P;
  std::vector<TH1D> Track_E;
  std::vector<TH1D> Track_Pt;
  std::vector<TH1D> Track_Eta;
  std::vector<TH1D> Track_Phi;
  std::vector<TH1D> Track_normalizedChi2;
  std::vector<TH1D> Track_numberOfValidHits;
  std::vector<TH1D> Track_charge;

  std::vector<TH1D> Muon1_Pt;
  std::vector<TH1D> Muon1_E;
  std::vector<TH1D> Muon1_P;
  std::vector<TH1D> Muon1_Eta;
  std::vector<TH1D> Muon1_Phi;

  std::vector<TH1D> Muon2_Pt;
  std::vector<TH1D> Muon2_E;
  std::vector<TH1D> Muon2_P;
  std::vector<TH1D> Muon2_Eta;
  std::vector<TH1D> Muon2_Phi;
		
  std::vector<TH1D> Muon1_isGlobal;
  std::vector<TH1D> Muon2_isGlobal;
  std::vector<TH1D> Muon1_isStandAlone;
  std::vector<TH1D> Muon2_isStandAlone;
  std::vector<TH1D> Muon1_isTracker;
  std::vector<TH1D> Muon2_isTracker;
  std::vector<TH1D> Muon1_isCalo;
  std::vector<TH1D> Muon2_isCalo;
  std::vector<TH1D> Muon1_TriggerMatchdR;
  std::vector<TH1D> Muon2_TriggerMatchdR;
  std::vector<TH1D> Track_TriggerMatchdR;



  std::vector<TH1D> Track_PtB;
  std::vector<TH1D> Track_EtaB;
  std::vector<TH1D> Muon1_PtB;
  std::vector<TH1D> Muon1_EtaB;
  std::vector<TH1D> Muon2_PtB;
  std::vector<TH1D> Muon2_EtaB;
  

  std::vector<TH1D> Track_PtC;
  std::vector<TH1D> Track_EtaC;
  std::vector<TH1D> Muon1_PtC;
  std::vector<TH1D> Muon1_EtaC;
  std::vector<TH1D> Muon2_PtC;
  std::vector<TH1D> Muon2_EtaC;
  

  std::vector<TH1D> Track_PtD;
  std::vector<TH1D> Track_EtaD;
  std::vector<TH1D> Muon1_PtD;
  std::vector<TH1D> Muon1_EtaD;
  std::vector<TH1D> Muon2_PtD;
  std::vector<TH1D> Muon2_EtaD;
 
 
  std::vector<TH1D> Track_PtE;
  std::vector<TH1D> Track_EtaE;
  std::vector<TH1D> Muon1_PtE;
  std::vector<TH1D> Muon1_EtaE;
  std::vector<TH1D> Muon2_PtE;
  std::vector<TH1D> Muon2_EtaE;
 
 
  std::vector<TH1D> Track_PtF;
  std::vector<TH1D> Track_EtaF;
  std::vector<TH1D> Muon1_PtF;
  std::vector<TH1D> Muon1_EtaF;
  std::vector<TH1D> Muon2_PtF;
  std::vector<TH1D> Muon2_EtaF;
 
  std::vector<TH1D> Muon1_PtF_peak;
  std::vector<TH1D> Muon1_PtF_sideband;
  std::vector<TH1D> Muon1_PtF_substracted;



 
};
#endif
