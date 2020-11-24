#ifndef SignalSelectorArun_h
#define SignalSelectorArun_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class SignalSelectorArun : public Selection {

 public:
  SignalSelectorArun(TString Name_, TString id_);
  virtual ~SignalSelectorArun();

  virtual void  Configure();
  virtual void  Finish();

  /////////////////////////////////////////////////////////
  // This is a cut enumerator, for other event cuts please
  // fill the enumerator with put a new enumarator with an 
  // understandanle name, for exmaple  enum cuts {TriggerOk=0,
  // PrimeVts, EventCut1, EventCut2, ..., NCuts};  
  // Do not remove/rename  the last enumerator   NCuts;

  enum cuts {TriggerOk=0,PrimeVtx,SignalCandidate,LeadingMuonPt,LeadingMuonPt1,LeadingMuonPt2,NCuts, Mu1PtCut, Mu2PtCut, Mu3PtCut, MuonID, PVRefit, PhiVeto, OmegaVeto, TriggerMatchMu1, TriggerMatchMu2,TriggerMatchMu3,TauMassCut,GenMatch};
  
  //InvMuon1G,InvMuon1T,InvMuon1S,InvMuon2G,InvMuon2T,InvMuon2S,InvMuon3G,InvMuon3T,InvMuon3S,InvTauG,InvTauT,InvTauS
  



 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();
  
 private:
  double tauMinMass_, tauMaxMass_;
  double tauMinSideBand_,tauMaxSideBand_;
  double tauMassResCutLow, tauMassResCutHigh;
  double phiVetoSigma, omegaVetoSigma;
  
  
  float bdt_cutA1_3glb;  float bdt_cutA2_3glb;  float bdt_cutB1_3glb;  float bdt_cutB2_3glb;  float bdt_cutC1_3glb;  float bdt_cutC2_3glb;  float bdt_cutA1_2glbTrk;  float bdt_cutA2_2glbTrk;  float bdt_cutB1_2glbTrk;  float bdt_cutB2_2glbTrk;  float bdt_cutC1_2glbTrk;  float bdt_cutC2_2glbTrk;

 // categorization variables
 bool MC;
 float category;
 Int_t threeGlobal;

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
 float var_tauMassRefit;

 // MiniTree variables
 Float_t m3m;
  Float_t dataMCtype;
  Float_t event_weight;
  Float_t bdt;
  Float_t rapidity;
  Float_t LumiScale;

  // file and ttree pointer
  TTree *T3MMiniTree;
  TFile *T3MFMiniTree;

 private:
  // Selection Variables
  // Initializhere your analysis histograms
  std::vector<TH1D> NVtx;

  std::vector<TH1D> LeadMuonPt;
  std::vector<TH1D> LeadMuonEta;
  std::vector<TH1D> LeadMuonPhi;
  
  std::vector<TH1D> LeadMuonPt1;
  std::vector<TH1D> LeadMuonEta1;
  std::vector<TH1D> LeadMuonPhi1;
  
  std::vector<TH1D> LeadMuonPt2;
  std::vector<TH1D> LeadMuonEta2;
  std::vector<TH1D> LeadMuonPhi2;
  
  std::vector<TH1D> TauPt;
  std::vector<TH1D> TauEta;
  std::vector<TH1D> TauPhi;
  
  std::vector<TH1D> InvMu1G;
  std::vector<TH1D> InvMu1T;
  std::vector<TH1D> InvMu1S;
  
  std::vector<TH1D> InvMu2G;
  std::vector<TH1D> InvMu2T;
  std::vector<TH1D> InvMu2S;
  
  std::vector<TH1D> InvMu3G;
  std::vector<TH1D> InvMu3T;
  std::vector<TH1D> InvMu3S;
  
  std::vector<TH1D> InvTG;
  std::vector<TH1D> InvTT;
  std::vector<TH1D> InvTS;
  
  std::vector<TH2D> TauMass_all_nophiVeto;
  std::vector<TH1D> NSignalCandidates;


};
#endif
