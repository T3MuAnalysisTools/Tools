/*
 * GlobalEventFit.h
 *
 *  Created on: Apr 16, 2015
 *      Author: zotz
 */

#ifndef GlobalEventFit_h
#define GlobalEventFit_h

#include "SimpleFits/FitSoftware/interface/TPTRObject.h"
#include "SimpleFits/FitSoftware/interface/GEFObject.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/PTObject.h"
#include "SimpleFits/FitSoftware/interface/LagrangeMultipliersFitter.h"

class GlobalEventFit{
 public:
  GlobalEventFit(TrackParticle Muon, LorentzVectorParticle A1, double Phi_Res, TVector3 PV, TMatrixTSym<double> PVCov);
  GlobalEventFit(TrackParticle Muon, LorentzVectorParticle A1, PTObject METminusNeutrino, TVector3 PV, TMatrixTSym<double> PVCov);
  GlobalEventFit(std::vector<LorentzVectorParticle> A1s, PTObject METminusNeutrino, TVector3 PV, TMatrixTSym<double> PVCov);
  GlobalEventFit(TrackParticle ChargedPion, LorentzVectorParticle NeutralPion, LorentzVectorParticle A1, PTObject MET, TVector3 PV, TMatrixTSym<double> PVCov);
  GlobalEventFit(TrackParticle OneProng,  LorentzVectorParticle MuonsTriplet, PTObject MET, TVector3 PV, TMatrixTSym<double> PVCov, bool Tau3Muons);
  virtual ~GlobalEventFit();

  GEFObject Fit();
  bool AmbiguitySolverByChi2(std::vector< std::vector<bool> > A1Fit, std::vector<bool> EventFit, std::vector<TVectorD> Chi2Vecs, int &IndexToReturn);
  bool AmbiguitySolverByChi2Minuit(std::vector<bool> EventFit, std::vector<double> Chi2s, int &IndexToReturn);
  std::vector<LorentzVectorParticle> FitDaughtersCorr(std::vector<LorentzVectorParticle> FitDaughters);
  const LorentzVectorParticle& getA1() const{ return A1_;}
  const std::vector<LorentzVectorParticle> getA1s() const{ return A1s_;}
  const GEFObject& getGEFObject() const{ return GEFObject_;}

  bool isConfigured() const{return isConfigured_;}
  bool isFit() const{return isFit_;}

  const TrackParticle getMuon() const{return Muon_;}
  const TVector3 getPV() const{return PV_;}
  const TMatrixTSym<double> getPVCov() const{return PVCov_;}
  const TVector3 getSV() const{return SV_;}
  const std::vector<TVector3> getSVs() const{return SVs_;}
  const TMatrixTSym<double> getSVCov() const{return SVCov_;}
  const std::vector< TMatrixTSym<double> > getSVCovs() const{return SVCovs_;}
  const TPTRObject getTPTRObject() const{return TPTRObject_;}
  const std::vector<TPTRObject> getTPTRObjects() const{return TPTRObjects_;}
  std::vector<bool> getFitStatuses() const{return fitstatuses_;}

  double getMassConstraint() const{return MassConstraint_;}
  void setMassConstraint(double MassConstraint){
    MassConstraint_ = MassConstraint;
    useDefaultMassConstraint_ = false;
  }
  void SetCorrectPt(bool correct){correctPt_ = correct;}  //default is set to true. Can be set to false inside your analysis to prevent the corrections of the reconstructed pt of both taus to preserve the hard constraints imposed on the fit resonance
  void setUseCollinearityTauMu(bool useCollinearityTauMu){useCollinearityTauMu_ = useCollinearityTauMu;}
  void setMinimizer(std::string mode);
  void setMinimizer(int Minimizer){minimizer_ = Minimizer;};

 protected:
  void Configure(TrackParticle Muon, LorentzVectorParticle A1, TVector3 PV, TMatrixTSym<double> PVCov);
  void Configure(std::vector<LorentzVectorParticle> A1s, TVector3 PV, TMatrixTSym<double> PVCov);
  void Configure(TrackParticle ChargedPion, LorentzVectorParticle NeutralPion, LorentzVectorParticle A1, TVector3 PV, TMatrixTSym<double> PVCov);
  void Configure(TrackParticle OneProng, LorentzVectorParticle MuonsTriplet, TVector3 PV, TMatrixTSym<double> PVCov, bool Tau3Muons);

  void ThreeProngTauReconstruction();
  bool IsAmbiguous(std::vector<bool> recostatus);
  PTObject SubtractNeutrinoFromMET(unsigned Ambiguity);
  PTObject AddA1(PTObject MET);
  PTObject AddA1s(PTObject MET);
  PTObject AddMuon(PTObject MET);
  PTObject AddNeutralPion(PTObject MET);
  PTObject AddTriplet(PTObject MET);
 private:
  int minimizer_;
  bool isConfigured_;
  bool isFit_;
  bool isValid_;
  bool pionDecay_;
  TPTRObject TPTRObject_;
  std::vector<TPTRObject> TPTRObjects_;
  GEFObject GEFObject_;
  TrackParticle Muon_;
  TrackParticle Track_;
  LorentzVectorParticle MuonsTriplet_;
  LorentzVectorParticle NeutralPion_;
  LorentzVectorParticle A1_;
  std::vector<LorentzVectorParticle> A1s_;
  TVector3 PV_, SV_;
  TMatrixTSym<double> PVCov_, SVCov_;
  std::vector<TVector3> SVs_;
  std::vector< TMatrixTSym<double> > SVCovs_;
  PTObject MET_;
  std::vector<PTObject> METminusNeutrino_;
  std::vector<bool> fitstatuses_;
  TMatrixD FitPar_;
  TMatrixDSym FitCov_;
  bool useMassConstraint_;
  double MassConstraint_;
  double Phi_Res_;
  int MaxIterations_;
  double MaxDelta_;
  double Epsilon_;
  bool useDefaultMassConstraint_, useFullRecoil_, correctPt_, useCollinearityTauMu_;
  bool ZTT3Mu_;
};


#endif /* GlobalEventFit_h */
