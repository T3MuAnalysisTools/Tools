#ifndef Validation_h
#define Validation_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class Validation : public Selection {

    public:
      Validation(TString Name_, TString id_);
      virtual ~Validation();
		
		virtual double deltaR(double, double, double, double);
      virtual void  Configure();
      virtual void  Finish();

      enum cuts {L1TOk=0,HLTOk,PrimeVtx,NCuts};

    protected:
      virtual void doEvent();
      virtual void Store_ExtraDist();

    private:
      
		// Selection Variables
      std::vector<TH1D> NVtx;
      std::vector<TH1D> NTracks;
		std::vector<TH1D> HLT_Tau3Mu;

      std::vector<TH1D> TripleMass;
		std::vector<TH1D> Muon1Muon2dR;
		std::vector<TH1D> Muon2Muon3dR;
		std::vector<TH1D> Muon1Muon3dR;
		
		std::vector<TH1D> Muon1_isGlobal;
      std::vector<TH1D> Muon2_isGlobal;
      std::vector<TH1D> Muon3_isGlobal;
      std::vector<TH1D> Muon1_isStandAlone;
      std::vector<TH1D> Muon2_isStandAlone;
      std::vector<TH1D> Muon3_isStandAlone;
      std::vector<TH1D> Muon1_isTracker;
      std::vector<TH1D> Muon2_isTracker;
      std::vector<TH1D> Muon3_isTracker;
      std::vector<TH1D> Muon1_isCalo;
      std::vector<TH1D> Muon2_isCalo;
      std::vector<TH1D> Muon3_isCalo;
      std::vector<TH1D> Muon1_isIsolationValid;
      std::vector<TH1D> Muon2_isIsolationValid;
      std::vector<TH1D> Muon3_isIsolationValid;
      std::vector<TH1D> Muon1_isTimeValid;
      std::vector<TH1D> Muon2_isTimeValid;
      std::vector<TH1D> Muon3_isTimeValid;
      std::vector<TH1D> Muon1_emEt03;
      std::vector<TH1D> Muon2_emEt03;
      std::vector<TH1D> Muon3_emEt03;
      std::vector<TH1D> Muon1_emVetoEt03;
      std::vector<TH1D> Muon2_emVetoEt03;
      std::vector<TH1D> Muon3_emVetoEt03;
      std::vector<TH1D> Muon1_hadEt03;
      std::vector<TH1D> Muon2_hadEt03;
      std::vector<TH1D> Muon3_hadEt03;
      std::vector<TH1D> Muon1_hadVetoEt03;
      std::vector<TH1D> Muon2_hadVetoEt03;
      std::vector<TH1D> Muon3_hadVetoEt03;
      std::vector<TH1D> Muon1_nJets03;
      std::vector<TH1D> Muon2_nJets03;
      std::vector<TH1D> Muon3_nJets03;
      std::vector<TH1D> Muon1_nTracks03;
      std::vector<TH1D> Muon2_nTracks03;
      std::vector<TH1D> Muon3_nTracks03;
      std::vector<TH1D> Muon1_sumPt03;
      std::vector<TH1D> Muon2_sumPt03;
      std::vector<TH1D> Muon3_sumPt03;
      std::vector<TH1D> Muon1_trackerVetoPt03;
      std::vector<TH1D> Muon2_trackerVetoPt03;
      std::vector<TH1D> Muon3_trackerVetoPt03;
      std::vector<TH1D> Muon1_emEt05;
      std::vector<TH1D> Muon2_emEt05;
      std::vector<TH1D> Muon3_emEt05;
      std::vector<TH1D> Muon1_emVetoEt05;
      std::vector<TH1D> Muon2_emVetoEt05;
      std::vector<TH1D> Muon3_emVetoEt05;
      std::vector<TH1D> Muon1_hadEt05;
      std::vector<TH1D> Muon2_hadEt05;
      std::vector<TH1D> Muon3_hadEt05;
      std::vector<TH1D> Muon1_hadVetoEt05;
      std::vector<TH1D> Muon2_hadVetoEt05;
      std::vector<TH1D> Muon3_hadVetoEt05;
      std::vector<TH1D> Muon1_nJets05;
      std::vector<TH1D> Muon2_nJets05;
      std::vector<TH1D> Muon3_nJets05;
      std::vector<TH1D> Muon1_nTracks05;
      std::vector<TH1D> Muon2_nTracks05;
      std::vector<TH1D> Muon3_nTracks05;
      std::vector<TH1D> Muon1_sumPt05;
      std::vector<TH1D> Muon2_sumPt05;
      std::vector<TH1D> Muon3_sumPt05;
      std::vector<TH1D> Muon1_trackerVetoPt05;
      std::vector<TH1D> Muon2_trackerVetoPt05;
      std::vector<TH1D> Muon3_trackerVetoPt05;
      std::vector<TH1D> Muon1_sumChargedHadronPt03;
      std::vector<TH1D> Muon2_sumChargedHadronPt03;
      std::vector<TH1D> Muon3_sumChargedHadronPt03;
      std::vector<TH1D> Muon1_sumChargedParticlePt03;
      std::vector<TH1D> Muon2_sumChargedParticlePt03;
      std::vector<TH1D> Muon3_sumChargedParticlePt03;
      std::vector<TH1D> Muon1_sumNeutralHadronEt03;
      std::vector<TH1D> Muon2_sumNeutralHadronEt03;
      std::vector<TH1D> Muon3_sumNeutralHadronEt03;
      std::vector<TH1D> Muon1_sumNeutralHadronEtHighThreshold03;
      std::vector<TH1D> Muon2_sumNeutralHadronEtHighThreshold03;
      std::vector<TH1D> Muon3_sumNeutralHadronEtHighThreshold03;
      std::vector<TH1D> Muon1_sumPhotonEt03;
      std::vector<TH1D> Muon2_sumPhotonEt03;
      std::vector<TH1D> Muon3_sumPhotonEt03;
      std::vector<TH1D> Muon1_sumPhotonEtHighThreshold03;
      std::vector<TH1D> Muon2_sumPhotonEtHighThreshold03;
      std::vector<TH1D> Muon3_sumPhotonEtHighThreshold03;
      std::vector<TH1D> Muon1_sumPUPt03;
      std::vector<TH1D> Muon2_sumPUPt03;
      std::vector<TH1D> Muon3_sumPUPt03;
      std::vector<TH1D> Muon1_sumChargedHadronPt04;
      std::vector<TH1D> Muon2_sumChargedHadronPt04;
      std::vector<TH1D> Muon3_sumChargedHadronPt04;
      std::vector<TH1D> Muon1_sumChargedParticlePt04;
      std::vector<TH1D> Muon2_sumChargedParticlePt04;
      std::vector<TH1D> Muon3_sumChargedParticlePt04;
      std::vector<TH1D> Muon1_sumNeutralHadronEt04;
      std::vector<TH1D> Muon2_sumNeutralHadronEt04;
      std::vector<TH1D> Muon3_sumNeutralHadronEt04;
      std::vector<TH1D> Muon1_sumNeutralHadronEtHighThreshold04;
      std::vector<TH1D> Muon2_sumNeutralHadronEtHighThreshold04;
      std::vector<TH1D> Muon3_sumNeutralHadronEtHighThreshold04;
      std::vector<TH1D> Muon1_sumPhotonEt04;
      std::vector<TH1D> Muon2_sumPhotonEt04;
      std::vector<TH1D> Muon3_sumPhotonEt04;
      std::vector<TH1D> Muon1_sumPhotonEtHighThreshold04;
      std::vector<TH1D> Muon2_sumPhotonEtHighThreshold04;
      std::vector<TH1D> Muon3_sumPhotonEtHighThreshold04;
      std::vector<TH1D> Muon1_sumPUPt04;
      std::vector<TH1D> Muon2_sumPUPt04;
      std::vector<TH1D> Muon3_sumPUPt04;
};
#endif
