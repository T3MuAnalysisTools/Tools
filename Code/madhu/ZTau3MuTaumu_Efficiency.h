#ifndef ZTau3MuTaumu_Efficiency_h
#define ZTau3MuTaumu_Efficiency_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "SimpleFits/FitSoftware/interface/Chi2VertexFitter.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/PTObject.h"

#include "SimpleFits/FitSoftware/interface/TPTRObject.h"
#include "SimpleFits/FitSoftware/interface/GEFObject.h"
#include "SimpleFits/FitSoftware/interface/GlobalEventFit.h"
#include "SimpleFits/FitSoftware/interface/TrackHelixVertexFitter.h"
#include "SimpleFits/FitSoftware/interface/PTObject.h"
#include "SimpleFits/FitSoftware/interface/TPTRObject.h"


class ZTau3MuTaumu_Efficiency : public Selection {

 public:
  ZTau3MuTaumu_Efficiency(TString Name_, TString id_);
  virtual ~ZTau3MuTaumu_Efficiency();

  virtual void  Configure();
  virtual void  Finish();
  
  class IsL1 {
  public:
      // returns true if dimuon L1 trigger satisfied
      bool operator()(float pt1, float pt2, float eta1, float eta2) { 
      if ( fabs(eta1-eta2)<1.81&&((pt1>9.9||pt2>9.9)||(fabs(eta1)<1.61&&fabs(eta2)<1.61)) ) return 1;
      return 0;
      }
  };
  
  class IsdR {
  public:
      // returns true if two muons have dR less than maxdR
      double operator()(double phi1, double phi2, double eta1, double eta2) {
          double dphi = phi1 - phi2;
          if(dphi >= TMath::Pi()) dphi = dphi-2*TMath::Pi();// makes dphi between -pi and pi
          if(dphi <=-TMath::Pi()) dphi = dphi+2*TMath::Pi();
          double dR = sqrt(pow(dphi,2) + pow(eta1 - eta2,2));
          
          return dR;
      }
  };
  
  class IsHLT {
  public:
      // returns true if HLT trigger satisfied
      bool operator()(float pt1, float pt2, float maxhlt1, float maxhlt2) { 
      if ( pt1>maxhlt1&&pt2>maxhlt2) return 1;
      return 0;
      }
  };

  /*
  enum cuts {nMuon=0,
             Conditions_Charge,
             Conditions_Global_ID,
             Conditions_dR,
             Conditions_InvMass,
             Conditions_L1,
             Conditions_HLT,
             Conditions_muPt,
             Conditions_3muPt,
             Conditions_Random,
             L1_TriggerOk,
             HLT_TriggerOk,
             SignalCandidate,
	     NCuts}; 
           
  enum cuts {L1_TriggerOk=0,
             HLT_TriggerOk,
             nMuon,
             Conditions_Charge,
             Conditions_Global_ID,
             Conditions_dR,
             Conditions_InvMass,
             Conditions_L1,
             Conditions_HLT,
             Conditions_muPt,
             Conditions_3muPt,
             Conditions_Random,
             SignalCandidate,
	     NCuts}; 

  */
           
  enum cuts {nMuon=0,
             Conditions_Charge,
             Conditions_Global_ID,
             Conditions_dR,
             Conditions_InvMass,
             Conditions_L1,
             Conditions_HLT,
             Conditions_muPt,
             Conditions_3muPt,
             Conditions_Random,
             L1_TriggerOk,
             HLT_TriggerOk,
             SignalCandidate,
	     NCuts}; 

 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

 private:


  std::vector<TH1D>   Tau3MuRelativeIsolation;
  std::vector<TH1D>   OppositeMuRelativeIsolation;
  std::vector<TH1D>   VisibleDiTauMass;
  std::vector<TH1D>   MTT;
  std::vector<TH1D>   TripletMass;
  std::vector<TH1D>   matched_pdgId;
  std::vector<TH1D>   matched_dR;

  std::vector<TH1D>   Muon1DRToTruth;
  std::vector<TH1D>   Muon2DRToTruth;
  std::vector<TH1D>   Muon3DRToTruth;
  std::vector<TH1D>   dR_betweenTruth_VisibleTaus;
  std::vector<TH1D>   PairMass_OppositeSign_dR12;
  std::vector<TH1D>   PairMass_OppositeSign_dR13;

  std::vector<TH1D>   TripletPt;
  std::vector<TH1D>   OppositeMuonPt;

  std::vector<TH1D>   TripletEta;
  std::vector<TH1D>   OppositeMuonEta;
  
  std::vector<TH1D>   Selection_Cut_Sequential_1;
  std::vector<TH1D>   Selection_Cut_Sequential_2;
  std::vector<TH1D>   Selection_Cut_Sequential_3;
  std::vector<TH1D>   Selection_Cut_Sequential_4;
  std::vector<TH1D>   Selection_Cut_Sequential_5;
  std::vector<TH1D>   Selection_Cut_Sequential_6;
  
  TRandom rndm;
  double random_num;

};
#endif
