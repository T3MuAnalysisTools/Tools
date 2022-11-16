#ifndef SimpleTauESelector_h
#define SimpleTauESelector_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class SimpleTauESelector : public Selection {

 public:
  SimpleTauESelector(TString Name_, TString id_);
  virtual ~SimpleTauESelector();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0, nElectrons , SignalCandidate, OSCharge, pTCut1, pTCut2, pTCut3, pTCutTauE,  NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:

  std::vector<TH1D>   NumberOfElectrons;
  std::vector<TH2D>   NumTausvsNumElectrons;
  
  std::vector<TH1D>   El_TauCand_Inv_Mass_1;
  std::vector<TH1D>   El_TauCand_Inv_Mass_2;
  std::vector<TH1D>   El_TauCand_Inv_Mass_3;
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
