#ifndef TMVA_Tree_3_twoGlbOneLooseMuon_h
#define TMVA_Tree_3_twoGlbOneLooseMuon_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class TMVA_Tree_3_twoGlbOneLooseMuon : public Selection {

 public:
  TMVA_Tree_3_twoGlbOneLooseMuon(TString Name_, TString id_);
  virtual ~TMVA_Tree_3_twoGlbOneLooseMuon();

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
 
float var_KFV_chiSq ;
float var_flightLenSig ; // Add flight length significance
float var_mindca_iso ;
float var_muMinPt ;
float var_inOutTrackMatch_Chi2 ;
float var_VertexMu3D0 ;
float var_relPt ;
float var_osMmumu1 ;
float var_osMmumu2 ;

};
#endif
