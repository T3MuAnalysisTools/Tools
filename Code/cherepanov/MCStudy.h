#ifndef MCStudy_h
#define MCStudy_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class MCStudy : public Selection {

 public:
  MCStudy(TString Name_, TString id_);
  virtual ~MCStudy();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,SignalCandidate, Mu1PtCut, Mu2PtCut, Mu3PtCut, MuonID, PhiVeto, OmegaVeto, TriggerMatch, ThreeMuMass,NCuts}; 


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

  std::vector<TH1D>   Muon1Eta_EtaSort;
  std::vector<TH1D>   Muon2Eta_EtaSort;
  std::vector<TH1D>   Muon3Eta_EtaSort;
  std::vector<TH1D>   Muon1Pt_EtaSort;
  std::vector<TH1D>   Muon2Pt_EtaSort;
  std::vector<TH1D>   Muon3Pt_EtaSort;


  std::vector<TH1D> TauMassResolution;
  std::vector<TH1D> TauMassResolutionRefit;

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


  std::vector<TH2D> Mu1PtvsEta;
  std::vector<TH2D> Mu2PtvsEta;
  std::vector<TH2D> Mu3PtvsEta;
  std::vector<TH2D> TauPtVsEta;



  std::vector<TH2D> TauMassResolutionVsEta;
  std::vector<TH2D> TauMassResolutionVsPt;

  std::vector<TH2D> TauMassResolutionVsMu1Eta;
  std::vector<TH2D> TauMassResolutionVsMu1Pt;

  std::vector<TH2D> TauMassResolutionVsMu2Eta;
  std::vector<TH2D> TauMassResolutionVsMu2Pt;

  std::vector<TH2D> TauMassResolutionVsMu3Eta;
  std::vector<TH2D> TauMassResolutionVsMu3Pt;


  std::vector<TH1D> PETauMassResolution_Pt;
  std::vector<TH1D> PETauMassResolution_PtEtaPhi;
  std::vector<TH1D> PETauMassResolution_PtEtaPhi_RefitTracks;

  std::vector<TH2D> PETauMassResVsMu1Pt;
  std::vector<TH2D> PETauMassResVsMu2Pt;
  std::vector<TH2D> PETauMassResVsMu3Pt;


  std::vector<TH2D> PETauMassResVsMu1Eta;
  std::vector<TH2D> PETauMassResVsMu2Eta;
  std::vector<TH2D> PETauMassResVsMu3Eta;


  std::vector<TH2D> Mu1PtvsPtError;
  std::vector<TH2D> Mu2PtvsPtError;
  std::vector<TH2D> Mu3PtvsPtError;

  std::vector<TH2D> Mu1EtavsEtaError;
  std::vector<TH2D> Mu2EtavsEtaError;
  std::vector<TH2D> Mu3EtavsPetError;

  std::vector<TH2D> Mu1PhivsPhiError;
  std::vector<TH2D> Mu2PhivsPhiError;
  std::vector<TH2D> Mu3PhivsPhiError;
  std::vector<TH3F> TauMassResolutionVsPtEta;

};
#endif
