#ifndef TMVA_Signal_glbMuon_h
#define TMVA_Signal_glbMuon_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class TMVA_Signal_glbMuon : public Selection {

 public:
  TMVA_Signal_glbMuon(TString Name_, TString id_);
  virtual ~TMVA_Signal_glbMuon();

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
  // Initializhere your analysis histograms
  std::vector<TH1D> NVtx;
  std::vector<TH1D> MuonsPt;
  std::vector<TH1D> MuonsEta;
  std::vector<TH1D> MuonsPhi;
  std::vector<TH1D> TripleMass;
 	
float var_muMinPt ;
float	var_inOutTrackMatch_Chi2 ;
float	var_mu3kink ;
float	var_vertexKFChi2 ;
float	var_VertexMu3D0Sig ;
float	var_VertexMu3D0 ;
float	var_mindca_iso ;
float	var_iso_relpt ;
float var_fv_cosdphi3d;

};
#endif
