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

  enum cuts {TriggerOk=0, nTaus , SignalCandidate, TauHIsolation, OSCharge, DeepTauVsJ, DeepTauVsMu, DeepTauVsE, pTCut1, pTCut2, pTCut3, pTCutTauX, DM,  NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:

  std::vector<TH1D>   NumberOfTaus;
  std::vector<TH2D>   NumTausvsNumMuons;
  
  std::vector<TH1D>   TauX_TauCand_Inv_Mass_1;
  std::vector<TH1D>   TauX_TauCand_Inv_Mass_2;
  std::vector<TH1D>   TauX_TauCand_Inv_Mass_3;
  std::vector<TH1D>   DecayModes;
  std::vector<TH1D>   MET_Et;
  std::vector<TH1D>   MET_Phi;
  std::vector<TH1D>   Tau3muMass;

  
  std::vector<TH1D>   Kinematics;
  std::vector<TH1D>   Kinematics_1;
  std::vector<TH1D>   Kinematics_MissingTrMass;
  std::vector<TH2D>   Kinematics_TauXPtEta;
  
  TRandom rndm;
  double random_num;

};
#endif
