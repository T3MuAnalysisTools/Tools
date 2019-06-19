#ifndef TMVA_Tree_1_looseMuon_h
#define TMVA_Tree_1_looseMuon_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class TMVA_Tree_1_looseMuon : public Selection {

 public:
  TMVA_Tree_1_looseMuon(TString Name_, TString id_);
  virtual ~TMVA_Tree_1_looseMuon();

  virtual void  Configure();
  virtual void  Finish();

enum cuts {TriggerOk=0,SignalCandidate, Mu1PtCut, Mu2PtCut, Mu3PtCut, MuonID, PhiVeto, OmegaVeto, TriggerMatch, ThreeMuMass,NCuts};

 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();
  TFile * file;
  TTree * TMVA_Tree;


 private:
  // Selection Variables
  double tauMinMass_, tauMaxMass_;
  double tauMinSideBand_,tauMaxSideBand_;
  double massRes_;
  // Initializhere your analysis histograms
  std::vector<TH1D> NVtx;
  std::vector<TH1D> MuonsPt;
  std::vector<TH1D> MuonsEta;
  std::vector<TH1D> MuonsPhi;
  std::vector<TH1D> TripleMass;

  float var_KFV_chiSq;
  float var_mu3d0sig;
  float var_mu3d0;
  float var_relPt;
  float var_ntrkMu3Iso;
  float var_dcaMu1Mu3;
  float var_mu3InOutMatch;
  float var_mu3Kink;
  float var_vertex2d;
};
#endif