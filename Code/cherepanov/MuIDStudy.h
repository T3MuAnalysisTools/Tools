#ifndef MuIDStudy_h
#define MuIDStudy_h

#include "Selection.h"
#include <vector>
#include "TString.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


class MuIDStudy : public Selection {

 public:
  MuIDStudy(TString Name_, TString id_);
  virtual ~MuIDStudy();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,SignalCandidate, Mu1PtCut, Mu2PtCut, Mu3PtCut, MuonID, PhiVeto, OmegaVeto, TriggerMatch, TauMassCut,NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:
  double tauMinMass_, tauMaxMass_;
  double tauMinSideBand_,tauMaxSideBand_;


  // Selection Variables
  // Initializhere your analysis histograms

  std::vector<TH1D> Muon1Pt;
  std::vector<TH1D> Muon2Pt;
  std::vector<TH1D> Muon3Pt;

  std::vector<TH1D> Muon1Eta;
  std::vector<TH1D> Muon2Eta;
  std::vector<TH1D> Muon3Eta;

  std::vector<TH1D> Significance3C;
  std::vector<TH1D> Significance2C;
  std::vector<TH1D> Significance1C;

  std::vector<TH1D> Significance3B;
  std::vector<TH1D> Significance2B;
  std::vector<TH1D> Significance1B;

  std::vector<TH1D> Significance3A;
  std::vector<TH1D> Significance2A;
  std::vector<TH1D> Significance1A;

  std::vector<TH1D> Significance3;
  std::vector<TH1D> Significance2;
  std::vector<TH1D> Significance1;

  std::vector<TH1D> TauEta;
  std::vector<TH1D> TauPt;
  std::vector<TH1D> TauP;

  std::vector<TH1D> TauMassResolution;
  std::vector<TH1D> EventMassResolution_PtEtaPhi;

  std::vector<TH1D> SVPVTauDirAngle;


  std::vector<TH2D> EMR_tau_eta;
  std::vector<TH2D> L1Triggers;


  std::vector<TH2D> L1TriggersB;
  std::vector<TH2D> L1TriggersC;
  std::vector<TH2D> L1TriggersD;
  std::vector<TH2D> L1TriggersE;
  std::vector<TH2D> L1TriggersF;


  std::vector<TH1D> Muon1isGlob;
  std::vector<TH1D> Muon2isGlob;
  std::vector<TH1D> Muon3isGlob;

  std::vector<TH1D> Muon1isStand;
  std::vector<TH1D> Muon2isStand;
  std::vector<TH1D> Muon3isStand;


  std::vector<TH1D> Muon1isTrack;
  std::vector<TH1D> Muon2isTrack;
  std::vector<TH1D> Muon3isTrack;

  std::vector<TH1D> Muon1PtResolution;
  std::vector<TH1D> Muon2PtResolution;
  std::vector<TH1D> Muon3PtResolution;

  std::vector<TH1D> Muon1EtaResolution;
  std::vector<TH1D> Muon2EtaResolution;
  std::vector<TH1D> Muon3EtaResolution;

  std::vector<TH1D> Muon1DRToTruth;
  std::vector<TH1D> Muon2DRToTruth;
  std::vector<TH1D> Muon3DRToTruth;

  std::vector<TH1D> CustMuonIdMC1;
  std::vector<TH1D> CustMuonIdData1;

  std::vector<TH1D> CustMuonIdMC2;
  std::vector<TH1D> CustMuonIdData2;

  std::vector<TH1D> CustMuonIdMC3;
  std::vector<TH1D> CustMuonIdData3;


  std::vector<TH1D> CustMuonIdMC1A;
  std::vector<TH1D> CustMuonIdData1A;

  std::vector<TH1D> CustMuonIdMC2A;
  std::vector<TH1D> CustMuonIdData2A;

  std::vector<TH1D> CustMuonIdMC3A;
  std::vector<TH1D> CustMuonIdData3A;




  std::vector<TH1D> CustMuonIdMC1B;
  std::vector<TH1D> CustMuonIdData1B;

  std::vector<TH1D> CustMuonIdMC2B;
  std::vector<TH1D> CustMuonIdData2B;

  std::vector<TH1D> CustMuonIdMC3B;
  std::vector<TH1D> CustMuonIdData3B;




  std::vector<TH1D> CustMuonIdMC1C;
  std::vector<TH1D> CustMuonIdData1C;

  std::vector<TH1D> CustMuonIdMC2C;
  std::vector<TH1D> CustMuonIdData2C;

  std::vector<TH1D> CustMuonIdMC3C;
  std::vector<TH1D> CustMuonIdData3C;








  std::vector<TH1D> TriggerMatchdR1;
  std::vector<TH1D> TriggerMatchdR2;
  std::vector<TH1D> TriggerMatchdR3;

  std::vector<TH1D> dR12;
  std::vector<TH1D> dR23;
  std::vector<TH1D> dR31;
  std::vector<TH1D> dR1Tau;
  std::vector<TH1D> dR2Tau;
  std::vector<TH1D> dR3Tau;

  std::vector<TH1D> VertexChi2KF;
  std::vector<TH1D> FLSignificance;
  std::vector<TH1D> NSignalCandidates;
  std::vector<TH1D> TauMassResolutionRefit;



  Float_t var_vertexKFChi2;// (chi sq of the fit of the secondary vertex)
  Float_t var_svpvTauAngle;// (The angle between PV-SV vector and the tau vector)
  Float_t var_flightLenSig;// (Flight length significance of the tau candidate)
  Float_t var_sumMuTrkKinkChi2;// (sum of chi sq of the kink of all three muons)
  Float_t var_segCompMuMin;// (Minimum of the segment compatibility of the three muons)
  Float_t var_MinMIPLikelihood;// (Minimum of the calorimeter compatibility of the three muons)
  Float_t var_tauMass;
  Float_t var_MuMu_mindR;
  Float_t var_RelPt_Mu1Tau;
  Float_t var_Eta_au;
  Float_t var_MuMu_minKFChi2;
  Float_t var_maxdca;
  Float_t var_MuTau_maxdR;
  Float_t var_MaxD0Significance;
  Float_t var_IsolationMinDist;



};
#endif
