#ifndef MuonPionTree_h
#define MuonPionTree_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class MuonPionTree : public Selection {

   public:
      MuonPionTree(TString Name_, TString id_);
      virtual ~MuonPionTree();

      virtual void  Configure();
      virtual void  Finish();

      /////////////////////////////////////////////////////////
      // This is a cut enumerator, for other event cuts please
      // fill the enumerator with put a new enumarator with an 
      // understandanle name, for exmaple  enum cuts {TriggerOk=0,
      // PrimeVts, EventCut1, EventCut2, ..., NCuts};  
      // Do not remove/rename  the last enumerator   NCuts;

      enum cuts {Candidate, NCuts}; 

   protected:
      virtual void doEvent();  
      virtual void Store_ExtraDist();

   private:

      TFile *file;
      TTree* TMVATree;

      // Selection Variables
      float PhiMassHigh;
      float PhiMassLow;
      float sidebandDsMin;
      float sidebandDsMax;
      float peakDsMin;
      float peakDsMax;
      float dsMassMin;
      float dsMassMax;

      int fake;
      int isGlobal;
      int isTracker;
      int isPF;
      float muonPt ;
      float muonEta ;
      float muonPhi ;
      float muonInnerNC2 ;
      float muonValidFraction ;
      int muonInnerNValidHits ;
      float muonInnerTrackQuality ;
      int muonNValidPixelHits ;
      int muonNValidTrackerHits ;
      int muonNLostTrackerHits ;
      int muonNLostTrackerHitsInner ;
      int muonNLostTrackerHitsOuter ;
      int muonPixelLayers ;
      int muonTrackerLayers ;
      int muonNMatchedStations ;
      int muonNMatches ;
      bool muonRPC ;
      float muonPtErrPt ;
      float muonSegComp ;
      float muonCaloComp ;
      float muonHadS9 ;
      float muonHad ;
      float muonEM ;
      float muonEMS9 ;
      float muonEMS25 ;
      float muonKink;
      // Initializhere your analysis histograms

      std::vector<TH2D> NDsPions;

};
#endif
