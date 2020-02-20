#ifndef FillMVATree_TwoGlobalTracker_h
#define FillMVATree_TwoGlobalTracker_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "TRandom.h"

class FillMVATree_TwoGlobalTracker : public Selection {

   public:
      FillMVATree_TwoGlobalTracker(TString Name_, TString id_);
      virtual ~FillMVATree_TwoGlobalTracker();

      virtual void  Configure();
      virtual void  Finish();

      enum cuts {TriggerOk=0,SignalCandidate, Mu1PtCut, Mu2PtCut, Mu3PtCut, TriggerMatchMu1, TriggerMatchMu2, TriggerMatchMu3, MuonID, PVRefit, PhiVetoOS1, OmegaVetoOS1, PhiVetoOS2,  OmegaVetoOS2, TauMassCut, DsGenMatch, GenMatch, NCuts};

   protected:
      virtual void doEvent();  
      virtual void Store_ExtraDist();

      template <typename T>
         int minQuantityIndex(std::vector<T>& vec);

      template <typename T>
         int maxQuantityIndex(std::vector<T>& vec);

      TFile * file;
      TTree * TMVA_Tree;

   private:

      // Selection Variables
      double tauMinMass_, tauMaxMass_;
      double tauMinSideBand_,tauMaxSideBand_;
      double tauMassResCutLow, tauMassResCutHigh;
      double phiVetoSigma, omegaVetoSigma;
      double M_osss1, M_osss2;

      unsigned int Muon_index_1, Muon_index_2, Muon_index_3;

      // Initializhere your analysis histograms
      std::vector<TH1D> L1Seed;
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
      bool threeGlobal;

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
      float var_minMatchedStations; // number of minimum stations matched to a muon


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

      // Additional plots
      std::vector<TH1D> Mu3_isGlobal;
      std::vector<TH1D> Mu2_isGlobal;
      std::vector<TH2D> L1TriggerMatch;
      std::vector<TH1D> MuonVsMinSC;
      std::vector<TH1D> MuonVsMaxCLP;
      std::vector<TH1D> MuonVsMinP;
      std::vector<TH1D> MuMuMass_OS1;
      std::vector<TH1D> MuMuMass_OS2;

      std::vector<TH1D> Muon1Pt;
      std::vector<TH1D> Muon2Pt;
      std::vector<TH1D> Muon3Pt;

      std::vector<TH1D> Muon1Eta;
      std::vector<TH1D> Muon2Eta;
      std::vector<TH1D> Muon3Eta;

      std::vector<TH1D> TauPt;
      std::vector<TH1D> TauP;
      std::vector<TH1D> TauMass;
      std::vector<TH1D> TauMassResolution;

      std::vector<TH1D> TauMassRefit;
      std::vector<TH1D> TauMassResolutionRefit;

      std::vector<TH1D> Muon1StandardSelector;
      std::vector<TH1D> Muon2StandardSelector;
      std::vector<TH1D> Muon3StandardSelector;

      std::vector<TH1D> Muon1isGlob;
      std::vector<TH1D> Muon2isGlob;
      std::vector<TH1D> Muon3isGlob;

      std::vector<TH1D> Muon1isTrack;
      std::vector<TH1D> Muon2isTrack;
      std::vector<TH1D> Muon3isTrack;

      std::vector<TH1D> Muon1kink;
      std::vector<TH1D> Muon2kink;
      std::vector<TH1D> Muon3kink;
      std::vector<TH1D> Muon1InOutTrackMatch;
      std::vector<TH1D> Muon2InOutTrackMatch;
      std::vector<TH1D> Muon3InOutTrackMatch;

      std::vector<TH1D> Muon1PtResolution;
      std::vector<TH1D> Muon2PtResolution;
      std::vector<TH1D> Muon3PtResolution;

      std::vector<TH1D> Muon1EtaResolution;
      std::vector<TH1D> Muon2EtaResolution;
      std::vector<TH1D> Muon3EtaResolution;

      std::vector<TH1D> Muon1DRToTruth;
      std::vector<TH1D> Muon2DRToTruth;
      std::vector<TH1D> Muon3DRToTruth;

      std::vector<TH1D> MuPair1_vertex_chi2;
      std::vector<TH1D> MuPair2_vertex_chi2;
      std::vector<TH1D> MuPair3_vertex_chi2;


      std::vector<TH1D> Pair1Mass;
      std::vector<TH1D> Pair2Mass;
      std::vector<TH1D> Pair3Mass;

      std::vector<TH1D> TriggerMatchdR1;
      std::vector<TH1D> TriggerMatchdR2;
      std::vector<TH1D> TriggerMatchdR3;

      std::vector<TH1D> dR12;
      std::vector<TH1D> dR23;
      std::vector<TH1D> dR31;
      std::vector<TH1D> dR1Tau;
      std::vector<TH1D> dR2Tau;
      std::vector<TH1D> dR3Tau;

      std::vector<TH1D> Isolation_NTracks;
      std::vector<TH1D> Isolation_RelPt;
      std::vector<TH1D> Isolation05_RelPt;
      std::vector<TH1D> Isolation05_NTracks;
      std::vector<TH1D> Isolation05_MinDist;

      std::vector<TH1D> VertexChi2AF;
      std::vector<TH1D> VertexDCA12;
      std::vector<TH1D> VertexDCA23;
      std::vector<TH1D> VertexDCA31;

      std::vector<TH1D> VertexSignalKFRefittedMu1P;
      std::vector<TH1D> VertexSignalKFRefittedMu1Pt;
      std::vector<TH1D> VertexSignalKFRefittedMu1Eta;
      std::vector<TH1D> VertexSignalKFRefittedMu1Phi;
      std::vector<TH1D> VertexSignalKFRefittedMu2P;
      std::vector<TH1D> VertexSignalKFRefittedMu2Pt;
      std::vector<TH1D> VertexSignalKFRefittedMu2Eta;
      std::vector<TH1D> VertexSignalKFRefittedMu2Phi;

      std::vector<TH1D> VertexSignalKFRefittedMu3P;
      std::vector<TH1D> VertexSignalKFRefittedMu3Pt;
      std::vector<TH1D> VertexSignalKFRefittedMu3Eta;
      std::vector<TH1D> VertexSignalKFRefittedMu3Phi;

      std::vector<TH1D> VertexMu1D0Reco;
      std::vector<TH1D> VertexMu1D0SigReco;
      std::vector<TH1D> VertexMu2D0Reco;
      std::vector<TH1D> VertexMu2D0SigReco;
      std::vector<TH1D> VertexMu3D0Reco;
      std::vector<TH1D> VertexMu3D0SigReco;
      std::vector<TH1D> Vertex2DDisplacement;
      std::vector<TH1D> Vertex3DDisplacement;
      std::vector<TH1D> VertexPairQuality;
      std::vector<TH1D> VertexPairRefitStatus;
      std::vector<TH1D> VertexPairfitStatus;
      std::vector<TH1D> VertexSignalKFChi2;
      std::vector<TH1D> VertexSignalAFX;
      std::vector<TH1D> VertexSignalAFY;
      std::vector<TH1D> VertexSignalAFZ;
      std::vector<TH1D> VertexSignalAFChi2;
      std::vector<TH1D> VertexSignalAFNdf;
      std::vector<TH1D> VertexMatchedPrimaryVertexX;
      std::vector<TH1D> VertexMatchedPrimaryVertexY;
      std::vector<TH1D> VertexMatchedPrimaryVertexZ;
      std::vector<TH1D> VertexRefitPVisValid;
      std::vector<TH1D> VertexMatchedRefitPrimaryVertexX;
      std::vector<TH1D> VertexMatchedRefitPrimaryVertexY;
      std::vector<TH1D> VertexMatchedRefitPrimaryVertexZ;

      std::vector<TH1D> Isolation_Mu1RelPt;
      std::vector<TH1D> Isolation_Mu2RelPt;
      std::vector<TH1D> Isolation_Mu3RelPt;

      std::vector<TH1D> FlightLengthSig;

      std::vector<TH2D> Muon1PtEta;
      std::vector<TH2D> Muon2PtEta;
      std::vector<TH2D> Muon3PtEta;

      std::vector<TH1D> MinMatchedStations;
      std::vector<TH1D> MaxMatchedStations;
      std::vector<TH1D> Mu1MatchedStations;
      std::vector<TH1D> Mu2MatchedStations;
      std::vector<TH1D> Mu3MatchedStations;

      std::vector<TH1D> NMuons;
      std::vector<TH1D> NTracks;
      std::vector<TH1D> NTwoMuonsTrack;

      std::vector<TH2D> Test;
      std::vector<TH2D> TrackCorrelation;
      std::vector<TH1D> MuonVstKink;
      std::vector<TH1D> MuonVsD0Sig;
      std::vector<TH1D> MuonVsRelPt;

};
#endif
