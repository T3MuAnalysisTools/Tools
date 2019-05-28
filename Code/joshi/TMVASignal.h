#ifndef TMVASignal_h
#define TMVASignal_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class TMVASignal : public Selection {

 public:
  TMVASignal(TString Name_, TString id_);
  virtual ~TMVASignal();

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
