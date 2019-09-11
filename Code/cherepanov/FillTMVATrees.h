#ifndef FillTMVATrees_h
#define FillTMVATrees_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class FillTMVATrees : public Selection {

 public:
  FillTMVATrees(TString Name_, TString id_);
  virtual ~FillTMVATrees();

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
  std::vector<TH1D> Muon_segmentCompatibility_mu1;
  std::vector<TH1D> Muon_segmentCompatibility_mu2;
  std::vector<TH1D> Muon_segmentCompatibility_mu3;
  std::vector<TH1D> Muon_segmentCompatibility_min;
  std::vector<TH1D> Muon_segmentCompatibility_max;
  std::vector<TH1D> Muon_ECALCompatibility_mu1;
  std::vector<TH1D> Muon_ECALCompatibility_mu2;
  std::vector<TH1D> Muon_ECALCompatibility_mu3;
  std::vector<TH1D> Muon_ECALCompatibility_min;
  std::vector<TH1D> Muon_ECALCompatibility_max;
  std::vector<TH1D> Isolation_NTracks;
  std::vector<TH1D> Isolation_RelPt;
  std::vector<TH1D> Isolation_maxdxy;
  std::vector<TH1D> VertexMu3D0SigReco;
  std::vector<TH1D> VertexDCAMax;
 
 bool MC;
 int category;
 float var_vertexKFChi2 ;
 float var_svpvTauAngle ;
 float var_flightLenSig ;
 float var_sumMuTrkKinkChi2 ;
 float var_segCompMuMin ;
 float var_segCompMuMax ;
 float var_segCompMu1 ;
 float var_segCompMu2 ;
 float var_segCompMu3 ;
 float var_caloCompMin ;
 float var_caloCompMax ;
 float var_caloCompMu1 ;
 float var_caloCompMu2 ;
 float var_caloCompMu3 ;
 float var_MinMIPLikelihood ;
 float var_tauMass ;
 float var_ntracks;
 float var_relPt;
 float var_isoMax;
 float var_mu3d0VertexSig;
 float var_maxdca;
};
#endif
