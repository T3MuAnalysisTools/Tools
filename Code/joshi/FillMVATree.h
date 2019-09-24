#ifndef FillMVATree_h
#define FillMVATree_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class FillMVATree : public Selection {

 public:
  FillMVATree(TString Name_, TString id_);
  virtual ~FillMVATree();

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
  double tauMassResCutLow, tauMassResCutHigh;
  // Initializhere your analysis histograms
  std::vector<TH1D> NVtx;
  std::vector<TH1D> MuonsPt;
  std::vector<TH1D> MuonsEta;
  std::vector<TH1D> MuonsPhi;
  std::vector<TH1D> TripleMass;

  std::vector<TH1D> SVPVTauDirAngle;
  std::vector<TH1D> FLSignificance;
  std::vector<TH1D> VertexChi2KF;
  std::vector<TH1D> MuonglbkinkSum;
  std::vector<TH1D> Muon_segmentCompatibility_min;
  std::vector<TH1D> Muon_HCALCompatibility_min;
  //New variables
  std::vector<TH1D> minMudR;
  std::vector<TH1D> dRMaxMuTau;
  std::vector<TH1D> Mu1TauPTRatio;
  std::vector<TH1D> MuPair_vertex_chi2_min;
  std::vector<TH1D> TauEta;
  std::vector<TH1D> VertexDCAMax;


 bool MC; 
 int category;
 float var_vertexKFChi2 ;
 float var_svpvTauAngle ;
 float var_flightLenSig ;
 float var_sumMuTrkKinkChi2 ;
 float var_segCompMuMin ;
 float var_MinMIPLikelihood ;
 float var_tauMass;
 // new variables
 float var_MuMu_mindR;
 float var_RelPt_Mu1Tau;
 float var_Eta_Tau;
 float var_MuMu_minKFChi2;
 float var_maxdca;
 float var_MuTau_maxdR;

};
#endif
