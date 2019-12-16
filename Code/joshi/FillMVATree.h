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

enum cuts {TriggerOk=0,SignalCandidate, Mu1PtCut, Mu2PtCut, Mu3PtCut, MuonID, PhiVeto, OmegaVeto, TriggerMatch, TauMassCut, GenMatch, NCuts};

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
  std::vector<TH1D> Isolation_MinDist;
  std::vector<TH1D> VertexMuMaxD0SigReco;
  std::vector<TH1D> EventMassResolution_PtEtaPhi;


 // categorization variables
 bool MC;
 float category;

 //commmon variables (2016 + 2017)
 float var_vertexKFChi2; // <= should be changed to normalized KF chi2
 float var_svpvTauAngle; 
 float var_flightLenSig;
 float var_segCompMuMin;
	 
 // 2016 variables
 float var_pmin; // Minimum p of the three muons
 float var_max_cLP; // Maximum chi square of the STA-TRK matching
 float var_max_tKink; // Maximum of the track kink of the 3 muons
 float var_MinD0Significance; // Minimum of the transverse IP significance of the 3 muons
 float var_mindca_iso; // Minimum DCA of tracks to muons with pT > 1 GeV (which muon?)
 float var_trk_relPt; // Ratio of sum of Pt of the tracks in muon isolation to muon (max value) [trk_pt>1 GeV, dR<0.03, dca<1 mm]
 float var_MinMIPLikelihood;


 // 2017 variables
 float var_MuMu_minKFChi2;
 float var_MuTau_maxdR;
 float var_sumMuTrkKinkChi2; // sum of chi square of STA-TRK matching of 3 muons
 float var_MaxD0Significance; // Minimum of the transverse IP significance of the 3 muons
 
 float var_Eta_Tau; 
 float var_tauMassRes;
 float var_tauMass;
 float var_MuMu_mindR;
 float var_RelPt_Mu1Ta;
 float var_maxdca;
 float var_RelPt_Mu1Tau;
};
#endif
