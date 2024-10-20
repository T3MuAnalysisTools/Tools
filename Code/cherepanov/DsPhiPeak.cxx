#include "DsPhiPeak.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

using namespace std;

DsPhiPeak::DsPhiPeak(TString Name_, TString id_):
   Selection(Name_,id_),
   PhiMassLow(1.00),
   PhiMassHigh(1.04),
   sidebandDsMin(1.70),
   sidebandDsMax(1.80),
   peakDsMin(1.92),
   peakDsMax(2.1),
   dsMassMin(1.68),
   dsMassMax(2.1)
{

   TString basedir = "";
   basedir = (TString)std::getenv("workdir")+"/Code/CommonFiles/PileUp/Collisions2018";
   PUWeightFile = new TFile(basedir+"/PUWeights_Run2018.root");
   puWeights = (TH1D*)PUWeightFile->Get("h1_weights");
}

DsPhiPeak::~DsPhiPeak(){
   for(unsigned int j=0; j<Npassed.size(); j++){
      Logger(Logger::Info) << "Selection Summary before: "
         << Npassed.at(j).GetBinContent(1)  << " +/- " << Npassed.at(j).GetBinError(1)  << " after: "
         << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
   }


   Logger(Logger::Info) << "complete." << std::endl;
}

void  DsPhiPeak::Configure(){

   for(int i=0; i<NCuts;i++){
      cut.push_back(0);
      value.push_back(0);
      pass.push_back(false);
      if(i==TriggerOk)          cut.at(TriggerOk)=1;
      if(i==TwoMuTrkCandidate)  cut.at(TwoMuTrkCandidate)=1;
      if(i==OSMuons)            cut.at(OSMuons)=0.1;
      if(i==Mu1PtCut)           cut.at(Mu1PtCut)=3.0;
      if(i==Mu2PtCut)           cut.at(Mu2PtCut)=3.0;
      if(i==TriggerMatchMu1)    cut.at(TriggerMatchMu1)=0.03;
      if(i==TriggerMatchMu2)    cut.at(TriggerMatchMu2)=0.03;
      if(i==TriggerMatchTrack)  cut.at(TriggerMatchTrack)=0.03;
      if(i==MuonID)             cut.at(MuonID)=1;
      if(i==MuMuMassCut)        cut.at(MuMuMassCut)=1;
      if(i==TrackPtCut)         cut.at(TrackPtCut)=2.0;
      if(i==NTrackHits)         cut.at(NTrackHits)=6;
      //if(i==ChiSqCut)         cut.at(ChiSqCut)=15; // Remove chi sqaure cut
      if(i==DsMassCut)          cut.at(DsMassCut)=1;
      if(i==GenMatch)           cut.at(GenMatch)=0.03;
   }
   TString hlabel;
   TString htitle;

   for(int i=0; i<NCuts; i++){
      title.push_back("");
      distindx.push_back(false);
      dist.push_back(std::vector<float>());
      TString c="_Cut_";c+=i;
      if(i==TriggerOk){
         title.at(i)="Pass HLT";
         hlabel="DoubleTrk_Trk_Tau3mu";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
      else if(i==TwoMuTrkCandidate){
         title.at(i)="is dimu+trk candidate";
         hlabel="is 2mu+trk candidate";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TwoMuTrkCandidate_",htitle,19,1.0,20.0,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TwoMuTrkCandidate_",htitle,19,1.0,20.0,hlabel,"Events"));
      }
      else if(i==Mu1PtCut){
         title.at(i)="Mu1 Pt $>$ ";
         title.at(i)+=cut.at(Mu1PtCut);
         title.at(i)+=" GeV";
         hlabel="Muon1 PT, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1PtCut_",htitle,50,0,25,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1PtCut_",htitle,50,0,25,hlabel,"Events"));
      }
      else if(i==OSMuons){
         title.at(i)="abs(Mu1+Mu2 Charge) < 0.1 ";
         hlabel="Sum of muon charges";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OSMuons_",htitle,5,-2.5,2.5,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OSMuons_",htitle,5,-2.5,2.5,hlabel,"Events"));
      }
      else if(i==Mu2PtCut){
         title.at(i)="Mu2 Pt $>$ ";
         title.at(i)+=cut.at(Mu2PtCut);
         title.at(i)+=" GeV";
         hlabel="Muon2 PT, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2PtCut_",htitle,40,0,20,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2PtCut_",htitle,40,0,20,hlabel,"Events"));
      }
      else if(i==MuonID){
         title.at(i)="All mu pass ID";
         hlabel="pass MuonID";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
      else if(i==MuMuMassCut){
         title.at(i)+="MuMuMass";
         title.at(i)+=" GeV";
         hlabel="InvMass(MuMu), GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuMuMassCut_",htitle,40,0.8,1.2,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuMuMassCut_",htitle,40,0.8,1.2,hlabel,"Events"));
      }
      else if(i==TrackPtCut){
         title.at(i)="Trk Pt $>$ ";
         title.at(i)+=cut.at(TrackPtCut);
         title.at(i)+=" GeV";
         hlabel="Track PT, GeV";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TrackPtCut_",htitle,40,0,20,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TrackPtCut_",htitle,40,0,20,hlabel,"Events"));
      }
      else if(i==NTrackHits){
         title.at(i)="Number of hits in tracker $>$";
         title.at(i)+=cut.at(NTrackHits);
         hlabel="N Tracker Hits";
         Nminus1.push_back(HConfig.GetTH1D(Name+"_Nminus1_NTrackHits_",htitle,50,0,50,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+"_Nminus0_NTrackHits_",htitle,50,0,50,hlabel,"Events"));
      }
      /*
         else if(i==ChiSqCut){
         title.at(i)="Chi Sq of vertex fitting $<$ ";
         title.at(i)+=cut.at(ChiSqCut);
         hlabel="Chi Sq of vertex fitting";
         Nminus1.push_back(HConfig.GetTH1D(Name+"_Nminus1_ChiSq_",htitle,50,0,50,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+"_Nminus0_ChiSq_",htitle,50,0,50,hlabel,"Events"));
         }
         */
      else if(i==TriggerMatchMu1){
         title.at(i)="Trigger Matching (Mu1)";
         hlabel="dR of muon 1 matching";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatchMu1_",htitle,50,0,0.1,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatchMu1_",htitle,50,0,0.1,hlabel,"Events"));
      }
      else if(i==TriggerMatchMu2){
         title.at(i)="Trigger Matching (Mu2)";
         hlabel="dR of muon 2 matching";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatchMu2_",htitle,50,0,0.1,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatchMu2_",htitle,50,0,0.1,hlabel,"Events"));
      }
      else if(i==TriggerMatchTrack){
         title.at(i)="Trigger Matching (Track)";
         hlabel="dR of track matching";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatchTrack_",htitle,50,0,0.1,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatchTrack_",htitle,50,0,0.1,hlabel,"Events"));
      }
      else if(i==GenMatch){
         title.at(i)="GEN matching";
         hlabel="GEN match";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_GENMatch_",htitle,30,0,0.03,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_GENMatch_",htitle,30,0,0.03,hlabel,"Events"));
      }
      else if(i==DsMassCut){
         title.at(i)="Ds Mass Cut [1.68, 2.02]";
         hlabel="Ds Mass Cut";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_DsMassCut_",htitle,85,1.4,2.25,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_DsMassCut_",htitle,85,1.4,2.25,hlabel,"Events"));
      }
   }

   // Setup NPassed Histogams
   Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove

   // Setup Extra Histograms

   // Book here your analysis histrogramms, a good style is to follow selfexplanatory convention
   // Unscaled histograms
   PhiMassVsDsMass=HConfig.GetTH2D(Name+"_PhiMassVsDsMass","#mu#mu Invariant mass vs. #mu#mu + track Invariant mass",50,0.77,1.27,42,1.68,2.1,"M_{#mu#mu} GeV","M_{#mu#mu + track}, GeV");
   TripleMass=HConfig.GetTH1D(Name+"_TripleMass","#mu#mu + track mass",42,1.68,2.1,"Invariant mass of the #mu#mu + track","Events");

   // ChargeMisId and Trigger efficiency measurement
   ChargeMisId=HConfig.GetTH1D(Name+"_ChargeMisId","Charge MisId",42,1.68,2.1,"M_{#mu#mu} GeV","Events");
   TriggerEfficiency=HConfig.GetTH1D(Name+"_TriggerEfficiency","Trigger Efficiency",42,1.68,2.1,"M_{#mu#mu}","Events");

   //*****************
   // Validation plots
   //*****************

   // Muon ID (peak) 
   Muon1_isGlobal_peak=HConfig.GetTH1D(Name+"_Muon1_isGlobal_peak","Global muons status (peak)",2,-.5,1.5,"#mu_{1} isGlb","Events"); peakCollection.push_back(&Muon1_isGlobal_peak);
   Muon2_isGlobal_peak=HConfig.GetTH1D(Name+"_Muon2_isGlobal_peak","(peak)",2,-0.5,0.5,"#mu_{2} isGlb","Events"); peakCollection.push_back(&Muon2_isGlobal_peak);
   Muon1_isStandAlone_peak=HConfig.GetTH1D(Name+"_Muon1_isStandAlone_peak","(peak)",2,-0.5,1.5,"#mu_{1} isStandAlone","Events"); peakCollection.push_back(&Muon1_isStandAlone_peak);
   Muon2_isStandAlone_peak=HConfig.GetTH1D(Name+"_Muon2_isStandAlone_peak","(peak)",2,-0.5,1.5,"#mu_{2} isStandAlone","Events"); peakCollection.push_back(&Muon2_isStandAlone_peak);
   Muon1_isTracker_peak=HConfig.GetTH1D(Name+"_Muon1_isTracker_peak","(peak)",2,-0.5,1.5,"#mu_{1} isTracker","Events"); peakCollection.push_back(&Muon1_isTracker_peak);
   Muon2_isTracker_peak=HConfig.GetTH1D(Name+"_Muon2_isTracker_peak","(peak)",2,-0.5,1.5,"#mu_{2} isTracker","Events"); peakCollection.push_back(&Muon2_isTracker_peak);
   Muon1_isCalo_peak=HConfig.GetTH1D(Name+"_Muon1_isCaloMuon_peak","(peak)",2,-0.5,1.5,"#mu_{1} isCalo","Events"); peakCollection.push_back(&Muon1_isCalo_peak);
   Muon2_isCalo_peak=HConfig.GetTH1D(Name+"_Muon2_isCaloMuon_peak","(peak)",2,-0.5,1.5,"#mu_{2} isCalo","Events"); peakCollection.push_back(&Muon2_isCalo_peak);
   Muon1_isIsolationValid_peak=HConfig.GetTH1D(Name+"_Muon1_isIsolationValid_peak","#mu_{1} isIsoValid(peak)",2,-0.5,1.5,"#mu_{1} isIsolationValid","Events"); peakCollection.push_back(&Muon1_isIsolationValid_peak);
   Muon2_isIsolationValid_peak=HConfig.GetTH1D(Name+"_Muon2_isIsolationValid_peak","#mu_{2} isIsoValid(peak)",2,-0.5,1.5,"#mu_{2} isIsolationValid","Events"); peakCollection.push_back(&Muon2_isIsolationValid_peak);

   // Muon ID (sideband) 
   Muon1_isGlobal_sideband=HConfig.GetTH1D(Name+"_Muon1_isGlobal_sideband","Global muons status (sideband)",2,-.5,1.5,"#mu_{1} isGlb","Events"); sidebandCollection.push_back(&Muon1_isGlobal_sideband);
   Muon2_isGlobal_sideband=HConfig.GetTH1D(Name+"_Muon2_isGlobal_sideband","(sideband)",2,-0.5,0.5,"#mu_{2} isGlb","Events"); sidebandCollection.push_back(&Muon2_isGlobal_sideband);
   Muon1_isStandAlone_sideband=HConfig.GetTH1D(Name+"_Muon1_isStandAlone_sideband","(sideband)",2,-0.5,1.5,"#mu_{1} isStandAlone","Events"); sidebandCollection.push_back(&Muon1_isStandAlone_sideband);
   Muon2_isStandAlone_sideband=HConfig.GetTH1D(Name+"_Muon2_isStandAlone_sideband","(sideband)",2,-0.5,1.5,"#mu_{2} isStandAlone","Events"); sidebandCollection.push_back(&Muon2_isStandAlone_sideband);
   Muon1_isTracker_sideband=HConfig.GetTH1D(Name+"_Muon1_isTracker_sideband","(sideband)",2,-0.5,1.5,"#mu_{1} isTracker","Events"); sidebandCollection.push_back(&Muon1_isTracker_sideband);
   Muon2_isTracker_sideband=HConfig.GetTH1D(Name+"_Muon2_isTracker_sideband","(sideband)",2,-0.5,1.5,"#mu_{2} isTracker","Events"); sidebandCollection.push_back(&Muon2_isTracker_sideband);
   Muon1_isCalo_sideband=HConfig.GetTH1D(Name+"_Muon1_isCaloMuon_sideband","(sideband)",2,-0.5,1.5,"#mu_{1} isCalo","Events"); sidebandCollection.push_back(&Muon1_isCalo_sideband);
   Muon2_isCalo_sideband=HConfig.GetTH1D(Name+"_Muon2_isCaloMuon_sideband","(sideband)",2,-0.5,1.5,"#mu_{2} isCalo","Events"); sidebandCollection.push_back(&Muon2_isCalo_sideband);
   Muon1_isIsolationValid_sideband=HConfig.GetTH1D(Name+"_Muon1_isIsolationValid_sideband","#mu_{1} isIsoValid(sideband)",2,-0.5,1.5,"#mu_{1} isIsolationValid","Events"); sidebandCollection.push_back(&Muon1_isIsolationValid_sideband);
   Muon2_isIsolationValid_sideband=HConfig.GetTH1D(Name+"_Muon2_isIsolationValid_sideband","#mu_{2} isIsoValid(sideband)",2,-0.5,1.5,"#mu_{2} isIsolationValid","Events"); sidebandCollection.push_back(&Muon2_isIsolationValid_sideband);

   // Muon ID (peak) 
   Muon1_isGlobal_validation=HConfig.GetTH1D(Name+"_Muon1_isGlobal_validation","Global muons status (validation)",2,-.5,1.5,"#mu_{1} isGlb","Events"); validationCollection.push_back(&Muon1_isGlobal_validation);
   Muon2_isGlobal_validation=HConfig.GetTH1D(Name+"_Muon2_isGlobal_validation","(validation)",2,-0.5,0.5,"#mu_{2} isGlb","Events"); validationCollection.push_back(&Muon2_isGlobal_validation);
   Muon1_isStandAlone_validation=HConfig.GetTH1D(Name+"_Muon1_isStandAlone_validation","(validation)",2,-0.5,1.5,"#mu_{1} isStandAlone","Events"); validationCollection.push_back(&Muon1_isStandAlone_validation);
   Muon2_isStandAlone_validation=HConfig.GetTH1D(Name+"_Muon2_isStandAlone_validation","(validation)",2,-0.5,1.5,"#mu_{2} isStandAlone","Events"); validationCollection.push_back(&Muon2_isStandAlone_validation);
   Muon1_isTracker_validation=HConfig.GetTH1D(Name+"_Muon1_isTracker_validation","(validation)",2,-0.5,1.5,"#mu_{1} isTracker","Events"); validationCollection.push_back(&Muon1_isTracker_validation);
   Muon2_isTracker_validation=HConfig.GetTH1D(Name+"_Muon2_isTracker_validation","(validation)",2,-0.5,1.5,"#mu_{2} isTracker","Events"); validationCollection.push_back(&Muon2_isTracker_validation);
   Muon1_isCalo_validation=HConfig.GetTH1D(Name+"_Muon1_isCaloMuon_validation","(validation)",2,-0.5,1.5,"#mu_{1} isCalo","Events"); validationCollection.push_back(&Muon1_isCalo_validation);
   Muon2_isCalo_validation=HConfig.GetTH1D(Name+"_Muon2_isCaloMuon_validation","(validation)",2,-0.5,1.5,"#mu_{2} isCalo","Events"); validationCollection.push_back(&Muon2_isCalo_validation);
   Muon1_isIsolationValid_validation=HConfig.GetTH1D(Name+"_Muon1_isIsolationValid_validation","#mu_{1} isIsoValid(validation)",2,-0.5,1.5,"#mu_{1} isIsolationValid","Events"); validationCollection.push_back(&Muon1_isIsolationValid_validation);
   Muon2_isIsolationValid_validation=HConfig.GetTH1D(Name+"_Muon2_isIsolationValid_validation","#mu_{2} isIsoValid(validation)",2,-0.5,1.5,"#mu_{2} isIsolationValid","Events"); validationCollection.push_back(&Muon2_isIsolationValid_validation);

   // Trigger Matching (peak)
   Track_TriggerMatchdR_peak=HConfig.GetTH1D(Name+"_Track_TriggerMatchdR_peak","track dR (trigger match)(peak)",30,0,0.03,"track dR (trigger match)","Events"); peakCollection.push_back(&Track_TriggerMatchdR_peak);
   Muon1_TriggerMatchdR_peak=HConfig.GetTH1D(Name+"_Muon1_TriggerMatchdR_peak","#mu_{1} dR (trigger match)(peak)",30,0,0.03,"#mu_{1} dR (trigger match)","Events"); peakCollection.push_back(&Muon1_TriggerMatchdR_peak);
   Muon2_TriggerMatchdR_peak=HConfig.GetTH1D(Name+"_Muon2_TriggerMatchdR_peak","#mu_{2} dR (trigger match)(peak)",30,0,0.03,"#mu_{2} dR (trigger match)","Events"); peakCollection.push_back(&Muon2_TriggerMatchdR_peak);

   PhiMass_peak=HConfig.GetTH1D(Name+"_PhiMass_peak","#mu#mu mass(peak)",50,0.77,1.27,"Invariant mass of the #mu#mu pair","Events"); peakCollection.push_back(&PhiMass_peak);

   // Trigger Matching (sideband)
   Track_TriggerMatchdR_sideband=HConfig.GetTH1D(Name+"_Track_TriggerMatchdR_sideband","track dR (trigger match)(sideband)",30,0,0.03,"track dR (trigger match)","Events"); sidebandCollection.push_back(&Track_TriggerMatchdR_sideband);
   Muon1_TriggerMatchdR_sideband=HConfig.GetTH1D(Name+"_Muon1_TriggerMatchdR_sideband","#mu_{1} dR (trigger match)(sideband)",30,0,0.03,"#mu_{1} dR (trigger match)","Events"); sidebandCollection.push_back(&Muon1_TriggerMatchdR_sideband);
   Muon2_TriggerMatchdR_sideband=HConfig.GetTH1D(Name+"_Muon2_TriggerMatchdR_sideband","#mu_{2} dR (trigger match)(sideband)",30,0,0.03,"#mu_{2} dR (trigger match)","Events"); sidebandCollection.push_back(&Muon2_TriggerMatchdR_sideband);

   PhiMass_sideband=HConfig.GetTH1D(Name+"_PhiMass_sideband","#mu#mu mass(sideband)",50,0.77,1.27,"Invariant mass of the #mu#mu pair","Events"); sidebandCollection.push_back(&PhiMass_sideband);

   // Trigger Matching (validation)
   Track_TriggerMatchdR_validation=HConfig.GetTH1D(Name+"_Track_TriggerMatchdR_validation","track dR (trigger match)(validation)",30,0,0.03,"track dR (trigger match)","Events"); validationCollection.push_back(&Track_TriggerMatchdR_validation);
   Muon1_TriggerMatchdR_validation=HConfig.GetTH1D(Name+"_Muon1_TriggerMatchdR_validation","#mu_{1} dR (trigger match)(validation)",30,0,0.03,"#mu_{1} dR (trigger match)","Events"); validationCollection.push_back(&Muon1_TriggerMatchdR_validation);
   Muon2_TriggerMatchdR_validation=HConfig.GetTH1D(Name+"_Muon2_TriggerMatchdR_validation","#mu_{2} dR (trigger match)(validation)",30,0,0.03,"#mu_{2} dR (trigger match)","Events"); validationCollection.push_back(&Muon2_TriggerMatchdR_validation);

   PhiMass_validation=HConfig.GetTH1D(Name+"_PhiMass_validation","#mu#mu mass(validation)",50,0.77,1.27,"Invariant mass of the #mu#mu pair","Events"); validationCollection.push_back(&PhiMass_validation);

   // Number of vertices (peak)
   NVtx_peak=HConfig.GetTH1D(Name+"_NVtx_peak","NVtx(peak)",100,0,100,"Number of Vertices (after reweighing)","Events"); peakCollection.push_back(&NVtx_peak);
   NVtx_woPUWeights_peak=HConfig.GetTH1D(Name+"_NVtx_woPUWeights_peak","NVtx(peak)",100,0,100,"Number of Vertices (before reweighing)","Events");  peakCollection.push_back(&NVtx_woPUWeights_peak);

   // Number of vertices (sideband)
   NVtx_sideband=HConfig.GetTH1D(Name+"_NVtx_sideband","NVtx(sideband)",100,0,100,"Number of Vertices (after reweighing)","Events"); sidebandCollection.push_back(&NVtx_sideband);
   NVtx_woPUWeights_sideband=HConfig.GetTH1D(Name+"_NVtx_woPUWeights_sideband","NVtx(sideband)",100,0,100,"Number of Vertices (before reweighing)","Events");  sidebandCollection.push_back(&NVtx_woPUWeights_sideband);

   // Number of vertices (validation)
   NVtx_validation=HConfig.GetTH1D(Name+"_NVtx_validation","NVtx(validation)",100,0,100,"Number of Vertices (after reweighing)","Events"); validationCollection.push_back(&NVtx_validation);
   NVtx_woPUWeights_validation=HConfig.GetTH1D(Name+"_NVtx_woPUWeights_validation","NVtx(validation)",100,0,100,"Number of Vertices (before reweighing)","Events");  validationCollection.push_back(&NVtx_woPUWeights_validation);

   //Dimuon Information (Muons from dimuon + track candidates) (validation)
   DimuondR_validation=HConfig.GetTH1D(Name+"_DimuondR_validation","dR between the muon pair (validation)",17,0,0.34,"dR","Events"); validationCollection.push_back(&DimuondR_validation);
   Muon1TrkdR_validation=HConfig.GetTH1D(Name+"_Muon1TrkdR_validation","dR between the highest p muon and the track (validation)",20,0,1,"dR (#mu_{1},track)","Events"); validationCollection.push_back(&Muon1TrkdR_validation);
   Muon2TrkdR_validation=HConfig.GetTH1D(Name+"_Muon2TrkdR_validation","dR between the lowest p muon and the track (validation)",20,0,1,"dR (#mu_{2},track)","Events"); validationCollection.push_back(&Muon2TrkdR_validation);

   //Dimuon Information (Muons from dimuon + track candidates) (sideband)
   DimuondR_sideband=HConfig.GetTH1D(Name+"_DimuondR_sideband","dR between the muon pair (sideband)",17,0,0.34,"dR","Events"); sidebandCollection.push_back(&DimuondR_sideband);
   Muon1TrkdR_sideband=HConfig.GetTH1D(Name+"_Muon1TrkdR_sideband","dR between the highest p muon and the track (sideband)",20,0,1,"dR (#mu_{1},track)","Events"); sidebandCollection.push_back(&Muon1TrkdR_sideband);
   Muon2TrkdR_sideband=HConfig.GetTH1D(Name+"_Muon2TrkdR_sideband","dR between the lowest p muon and the track (sideband)",20,0,1,"dR (#mu_{2},track)","Events");
   sidebandCollection.push_back(&Muon2TrkdR_sideband);

   //Dimuon Information (Muons from dimuon + track candidates) (peak)
   DimuondR_peak=HConfig.GetTH1D(Name+"_DimuondR_peak","dR between the muon pair (peak)",17,0,0.34,"dR","Events"); peakCollection.push_back(&DimuondR_peak);
   Muon1TrkdR_peak=HConfig.GetTH1D(Name+"_Muon1TrkdR_peak","dR between the highest p muon and the track (peak)",20,0,1,"dR (#mu_{1},track)","Events"); peakCollection.push_back(&Muon1TrkdR_peak);
   Muon2TrkdR_peak=HConfig.GetTH1D(Name+"_Muon2TrkdR_peak","dR between the lowest p muon and the track (peak)",20,0,1,"dR (#mu_{2},track)","Events");   peakCollection.push_back(&Muon2TrkdR_peak);

   // Track Candidate Information (peak)
   Track_Pt_peak=HConfig.GetTH1D(Name+"_Track_Pt_peak","Transverse momentum of track (2mu+trk track candidate) (peak)",50,0,50,"p_{T} (track)","Events"); peakCollection.push_back(&Track_Pt_peak);
   Track_Eta_peak=HConfig.GetTH1D(Name+"_Track_Eta_peak","Psuedorapidity of track (2mu+trk track candidate) (peak)",50,-2.5,2.5,"#eta","Events"); peakCollection.push_back(&Track_Eta_peak);
   Track_Phi_peak=HConfig.GetTH1D(Name+"_Track_Phi_peak","Azimuthal angle of track (2mu+trk track candidate) (peak)",63,-3.15,3.15,"#phi","Events"); peakCollection.push_back(&Track_Phi_peak);
   Track_vx_peak=HConfig.GetTH1D(Name+"_Track_vx_peak","X coordinate of the parent vertex (2mu+trk track candidate) (peak)",50,0,0.2,"Parent vertex x coordinate (cm)","Events"); peakCollection.push_back(&Track_vx_peak);
   Track_vy_peak=HConfig.GetTH1D(Name+"_Track_vy_peak","Y coordinate of the parent vertex (2mu+trk track candidate) (peak)",50,0,0.2,"Parent vertex y coordinate (cm)","Events"); peakCollection.push_back(&Track_vy_peak);
   Track_vz_peak=HConfig.GetTH1D(Name+"_Track_vz_peak","Z coordinate of the parent vertex (2mu+trk track candidate) (peak)",50,0,10,"Parent vertex z coordinate (cm)","Events"); peakCollection.push_back(&Track_vz_peak);
   Track_normalizedChi2_peak=HConfig.GetTH1D(Name+"_Track_normalizedChi2_peak","Normalized chi square (peak)",20,0,10,"#chi^{2}/ndf (track fit)","Events"); peakCollection.push_back(&Track_normalizedChi2_peak);
   Track_numberOfValidHits_peak=HConfig.GetTH1D(Name+"_Track_numberOfValidHits_peak","number of valid hits in the tracker (peak)",50,0,50,"n valid track hits","Events"); peakCollection.push_back(&Track_numberOfValidHits_peak);
   Track_charge_peak=HConfig.GetTH1D(Name+"_Track_charge_peak","Charge of the track (peak)",3,-1.0,1.0,"Number of Vertices","Events"); peakCollection.push_back(&Track_charge_peak);
   Track_dxy_peak=HConfig.GetTH1D(Name+"_Track_dxy_peak","Transverse displacement of the parent vertex from the bs (peak)",60,0,0.3,"dxy (cm)","Events"); peakCollection.push_back(&Track_dxy_peak);
   Track_dz_peak=HConfig.GetTH1D(Name+"_Track_dz_peak","Longitudnal displacement of the parent vertex from the bs (peak)",50,0,10,"dz (cm)","Events"); peakCollection.push_back(&Track_dz_peak);

   // Track Candidate Information (sideband)
   Track_Pt_sideband=HConfig.GetTH1D(Name+"_Track_Pt_sideband","Transverse momentum of track (2mu+trk track candidate) (sideband)",50,0,50,"p_{T} (track)","Events"); sidebandCollection.push_back(&Track_Pt_sideband);
   Track_Eta_sideband=HConfig.GetTH1D(Name+"_Track_Eta_sideband","Psuedorapidity of track (2mu+trk track candidate) (sideband)",50,-2.5,2.5,"#eta","Events"); sidebandCollection.push_back(&Track_Eta_sideband);
   Track_Phi_sideband=HConfig.GetTH1D(Name+"_Track_Phi_sideband","Azimuthal angle of track (2mu+trk track candidate) (sideband)",63,-3.15,3.15,"#phi","Events"); sidebandCollection.push_back(&Track_Phi_sideband);
   Track_vx_sideband=HConfig.GetTH1D(Name+"_Track_vx_sideband","X coordinate of the parent vertex (2mu+trk track candidate) (sideband)",50,0,0.2,"Parent vertex x coordinate (cm)","Events"); sidebandCollection.push_back(&Track_vx_sideband);
   Track_vy_sideband=HConfig.GetTH1D(Name+"_Track_vy_sideband","Y coordinate of the parent vertex (2mu+trk track candidate) (sideband)",50,0,0.2,"Parent vertex y coordinate (cm)","Events"); sidebandCollection.push_back(&Track_vy_sideband);
   Track_vz_sideband=HConfig.GetTH1D(Name+"_Track_vz_sideband","Z coordinate of the parent vertex (2mu+trk track candidate) (sideband)",50,0,10,"Parent vertex z coordinate (cm)","Events"); sidebandCollection.push_back(&Track_vz_sideband);
   Track_normalizedChi2_sideband=HConfig.GetTH1D(Name+"_Track_normalizedChi2_sideband","Normalized chi square (sideband)",20,0,10,"#chi^{2}/ndf (track fit)","Events"); sidebandCollection.push_back(&Track_normalizedChi2_sideband);
   Track_numberOfValidHits_sideband=HConfig.GetTH1D(Name+"_Track_numberOfValidHits_sideband","number of valid hits in the tracker (sideband)",50,0,50,"n valid track hits","Events"); sidebandCollection.push_back(&Track_numberOfValidHits_sideband);
   Track_charge_sideband=HConfig.GetTH1D(Name+"_Track_charge_sideband","Charge of the track (sideband)",3,-1.0,1.0,"Number of Vertices","Events"); sidebandCollection.push_back(&Track_charge_sideband);
   Track_dxy_sideband=HConfig.GetTH1D(Name+"_Track_dxy_sideband","Transverse displacement of the parent vertex from the bs (sideband)",60,0,0.3,"dxy (cm)","Events"); sidebandCollection.push_back(&Track_dxy_sideband);
   Track_dz_sideband=HConfig.GetTH1D(Name+"_Track_dz_sideband","Longitudnal displacement of the parent vertex from the bs (sideband)",50,0,10,"dz (cm)","Events"); sidebandCollection.push_back(&Track_dz_sideband);

   // Track Candidate Information (validation)
   Track_Pt_validation=HConfig.GetTH1D(Name+"_Track_Pt_validation","Transverse momentum of track (2mu+trk track candidate) (validation)",50,0,50,"p_{T} (track)","Events"); validationCollection.push_back(&Track_Pt_validation);
   Track_Eta_validation=HConfig.GetTH1D(Name+"_Track_Eta_validation","Psuedorapidity of track (2mu+trk track candidate) (validation)",50,-2.5,2.5,"#eta","Events"); validationCollection.push_back(&Track_Eta_validation);
   Track_Phi_validation=HConfig.GetTH1D(Name+"_Track_Phi_validation","Azimuthal angle of track (2mu+trk track candidate) (validation)",63,-3.15,3.15,"#phi","Events"); validationCollection.push_back(&Track_Phi_validation);
   Track_vx_validation=HConfig.GetTH1D(Name+"_Track_vx_validation","X coordinate of the parent vertex (2mu+trk track candidate) (validation)",50,0,0.2,"Parent vertex x coordinate (cm)","Events"); validationCollection.push_back(&Track_vx_validation);
   Track_vy_validation=HConfig.GetTH1D(Name+"_Track_vy_validation","Y coordinate of the parent vertex (2mu+trk track candidate) (validation)",50,0,0.2,"Parent vertex y coordinate (cm)","Events"); validationCollection.push_back(&Track_vy_validation);
   Track_vz_validation=HConfig.GetTH1D(Name+"_Track_vz_validation","Z coordinate of the parent vertex (2mu+trk track candidate) (validation)",50,0,10,"Parent vertex z coordinate (cm)","Events"); validationCollection.push_back(&Track_vz_validation);
   Track_normalizedChi2_validation=HConfig.GetTH1D(Name+"_Track_normalizedChi2_validation","Normalized chi square (validation)",20,0,10,"#chi^{2}/ndf (track fit)","Events"); validationCollection.push_back(&Track_normalizedChi2_validation);
   Track_numberOfValidHits_validation=HConfig.GetTH1D(Name+"_Track_numberOfValidHits_validation","number of valid hits in the tracker (validation)",50,0,50,"n valid track hits","Events"); validationCollection.push_back(&Track_numberOfValidHits_validation);
   Track_charge_validation=HConfig.GetTH1D(Name+"_Track_charge_validation","Charge of the track (validation)",3,-1.0,1.0,"Number of Vertices","Events"); validationCollection.push_back(&Track_charge_validation);
   Track_dxy_validation=HConfig.GetTH1D(Name+"_Track_dxy_validation","Transverse displacement of the parent vertex from the bs (validation)",60,0,0.3,"dxy (cm)","Events"); validationCollection.push_back(&Track_dxy_validation);
   Track_dz_validation=HConfig.GetTH1D(Name+"_Track_dz_validation","Longitudnal displacement of the parent vertex from the bs (validation)",50,0,10,"dz (cm)","Events"); validationCollection.push_back(&Track_dz_validation);

   // Muon variables (Muons from dimuon + track candidates) (peak)
   Muon1_Pt_peak=HConfig.GetTH1D(Name+"_Muon1_Pt_peak","Transverse Pt (muon 1) (peak)",25,0,50,"#mu_{1} p_{T} (GeV)","Events"); peakCollection.push_back(&Muon1_Pt_peak);
   Muon1_Eta_peak=HConfig.GetTH1D(Name+"_Muon1_Eta_peak","Psuedorapidity (muon 1) (peak)",25,-2.5,2.5,"#mu_{1} #eta","Events"); peakCollection.push_back(&Muon1_Eta_peak);
   Muon1_Phi_peak=HConfig.GetTH1D(Name+"_Muon1_Phi_peak","Azimuthal angle of (muons 1) (peak)",63,-3.15,3.15,"#mu_{1} #phi","Events");  peakCollection.push_back(&Muon1_Phi_peak);
   Muon1_vx_peak=HConfig.GetTH1D(Name+"_Muon1_Vx_peak","X coordinate of the parent vertex all muons (peak)",50,0,0.2,"#mu_{1} vx","Events");  peakCollection.push_back(&Muon1_vx_peak);
   Muon1_vy_peak=HConfig.GetTH1D(Name+"_Muon1_Vy_peak","Y coordinate of the parent vertex all muons (peak)",50,0,0.2,"#mu_{1} vy","Events");  peakCollection.push_back(&Muon1_vy_peak);
   Muon1_vz_peak=HConfig.GetTH1D(Name+"_Muon1_Vz_peak","Z coordinate of the parent vertex all muons (peak)",50,0,10,"#mu_{1} vz","Events"); peakCollection.push_back(&Muon1_vz_peak);
   Muon2_Pt_peak=HConfig.GetTH1D(Name+"_Muon2_Pt_peak","Transverse Pt (muon 2) (peak)",25,0,50,"#mu_{2} p_{T} (GeV)","Events"); peakCollection.push_back(&Muon2_Pt_peak);
   Muon2_Eta_peak=HConfig.GetTH1D(Name+"_Muon2_Eta_peak","Psuedorapidity (muon 2) (peak)",25,-2.5,2.5,"#mu_{2} #eta","Events"); peakCollection.push_back(&Muon2_Eta_peak);
   Muon2_Phi_peak=HConfig.GetTH1D(Name+"_Muon2_Phi_peak","Azimuthal angle of (muons 1) (peak)",63,-3.15,3.15,"#mu_{2} #phi","Events");  peakCollection.push_back(&Muon2_Phi_peak);
   Muon2_vx_peak=HConfig.GetTH1D(Name+"_Muon2_Vx_peak","X coordinate of the parent vertex all muons (peak)",50,0,0.2,"#mu_{2} vx","Events");  peakCollection.push_back(&Muon2_vx_peak);
   Muon2_vy_peak=HConfig.GetTH1D(Name+"_Muon2_Vy_peak","Y coordinate of the parent vertex all muons (peak)",50,0,0.2,"#mu_{2} vy","Events");  peakCollection.push_back(&Muon2_vy_peak);
   Muon2_vz_peak=HConfig.GetTH1D(Name+"_Muon2_Vz_peak","Z coordinate of the parent vertex all muons (peak)",50,0,10,"#mu_{2} vz","Events"); peakCollection.push_back(&Muon2_vz_peak);

   // Muon variables (Muons from dimuon + track candidates) (sideband)
   Muon1_Pt_sideband=HConfig.GetTH1D(Name+"_Muon1_Pt_sideband","Transverse Pt (muon 1) (sideband)",25,0,50,"#mu_{1} p_{T} (GeV)","Events"); sidebandCollection.push_back(&Muon1_Pt_sideband);
   Muon1_Eta_sideband=HConfig.GetTH1D(Name+"_Muon1_Eta_sideband","Psuedorapidity (muon 1) (sideband)",25,-2.5,2.5,"#mu_{1} #eta","Events"); sidebandCollection.push_back(&Muon1_Eta_sideband);
   Muon1_Phi_sideband=HConfig.GetTH1D(Name+"_Muon1_Phi_sideband","Azimuthal angle of (muons 1) (sideband)",63,-3.15,3.15,"#mu_{1} #phi","Events");  sidebandCollection.push_back(&Muon1_Phi_sideband);
   Muon1_vx_sideband=HConfig.GetTH1D(Name+"_Muon1_Vx_sideband","X coordinate of the parent vertex all muons (sideband)",50,0,0.2,"#mu_{1} vx","Events");  sidebandCollection.push_back(&Muon1_vx_sideband);
   Muon1_vy_sideband=HConfig.GetTH1D(Name+"_Muon1_Vy_sideband","Y coordinate of the parent vertex all muons (sideband)",50,0,0.2,"#mu_{1} vy","Events");  sidebandCollection.push_back(&Muon1_vy_sideband);
   Muon1_vz_sideband=HConfig.GetTH1D(Name+"_Muon1_Vz_sideband","Z coordinate of the parent vertex all muons (sideband)",50,0,10,"#mu_{1} vz","Events"); sidebandCollection.push_back(&Muon1_vz_sideband);
   Muon2_Pt_sideband=HConfig.GetTH1D(Name+"_Muon2_Pt_sideband","Transverse Pt (muon 2) (sideband)",25,0,50,"#mu_{2} p_{T} (GeV)","Events"); sidebandCollection.push_back(&Muon2_Pt_sideband);
   Muon2_Eta_sideband=HConfig.GetTH1D(Name+"_Muon2_Eta_sideband","Psuedorapidity (muon 2) (sideband)",25,-2.5,2.5,"#mu_{2} #eta","Events"); sidebandCollection.push_back(&Muon2_Eta_sideband);
   Muon2_Phi_sideband=HConfig.GetTH1D(Name+"_Muon2_Phi_sideband","Azimuthal angle of (muons 1) (sideband)",63,-3.15,3.15,"#mu_{2} #phi","Events");  sidebandCollection.push_back(&Muon2_Phi_sideband);
   Muon2_vx_sideband=HConfig.GetTH1D(Name+"_Muon2_Vx_sideband","X coordinate of the parent vertex all muons (sideband)",50,0,0.2,"#mu_{2} vx","Events");  sidebandCollection.push_back(&Muon2_vx_sideband);
   Muon2_vy_sideband=HConfig.GetTH1D(Name+"_Muon2_Vy_sideband","Y coordinate of the parent vertex all muons (sideband)",50,0,0.2,"#mu_{2} vy","Events");  sidebandCollection.push_back(&Muon2_vy_sideband);
   Muon2_vz_sideband=HConfig.GetTH1D(Name+"_Muon2_Vz_sideband","Z coordinate of the parent vertex all muons (sideband)",50,0,10,"#mu_{2} vz","Events"); sidebandCollection.push_back(&Muon2_vz_sideband);

   // Muon variables (Muons from dimuon + track candidates) (validation)
   Muon1_Pt_validation=HConfig.GetTH1D(Name+"_Muon1_Pt_validation","Transverse Pt (muon 1) (validation)",25,0,50,"#mu_{1} p_{T} (GeV)","Events"); validationCollection.push_back(&Muon1_Pt_validation);
   Muon1_Eta_validation=HConfig.GetTH1D(Name+"_Muon1_Eta_validation","Psuedorapidity (muon 1) (validation)",25,-2.5,2.5,"#mu_{1} #eta","Events"); validationCollection.push_back(&Muon1_Eta_validation);
   Muon1_Phi_validation=HConfig.GetTH1D(Name+"_Muon1_Phi_validation","Azimuthal angle of (muons 1) (validation)",63,-3.15,3.15,"#mu_{1} #phi","Events");  validationCollection.push_back(&Muon1_Phi_validation);
   Muon1_vx_validation=HConfig.GetTH1D(Name+"_Muon1_Vx_validation","X coordinate of the parent vertex all muons (validation)",50,0,0.2,"#mu_{1} vx","Events");  validationCollection.push_back(&Muon1_vx_validation);
   Muon1_vy_validation=HConfig.GetTH1D(Name+"_Muon1_Vy_validation","Y coordinate of the parent vertex all muons (validation)",50,0,0.2,"#mu_{1} vy","Events");  validationCollection.push_back(&Muon1_vy_validation);
   Muon1_vz_validation=HConfig.GetTH1D(Name+"_Muon1_Vz_validation","Z coordinate of the parent vertex all muons (validation)",100,0,5,"#mu_{1} vz","Events"); validationCollection.push_back(&Muon1_vz_validation);
   Muon2_Pt_validation=HConfig.GetTH1D(Name+"_Muon2_Pt_validation","Transverse Pt (muon 2) (validation)",25,0,50,"#mu_{2} p_{T} (GeV)","Events"); validationCollection.push_back(&Muon2_Pt_validation);
   Muon2_Eta_validation=HConfig.GetTH1D(Name+"_Muon2_Eta_validation","Psuedorapidity (muon 2) (validation)",25,-2.5,2.5,"#mu_{2} #eta","Events"); validationCollection.push_back(&Muon2_Eta_validation);
   Muon2_Phi_validation=HConfig.GetTH1D(Name+"_Muon2_Phi_validation","Azimuthal angle of (muons 1) (validation)",63,-3.15,3.15,"#mu_{2} #phi","Events");  validationCollection.push_back(&Muon2_Phi_validation);
   Muon2_vx_validation=HConfig.GetTH1D(Name+"_Muon2_Vx_validation","X coordinate of the parent vertex all muons (validation)",50,0,0.2,"#mu_{2} vx","Events");  validationCollection.push_back(&Muon2_vx_validation);
   Muon2_vy_validation=HConfig.GetTH1D(Name+"_Muon2_Vy_validation","Y coordinate of the parent vertex all muons (validation)",50,0,0.2,"#mu_{2} vy","Events");  validationCollection.push_back(&Muon2_vy_validation);
   Muon2_vz_validation=HConfig.GetTH1D(Name+"_Muon2_Vz_validation","Z coordinate of the parent vertex all muons (validation)",50,0,10,"#mu_{2} vz","Events"); validationCollection.push_back(&Muon2_vz_validation);

   //*****************
   // plots from Nick
   //*****************

   // Vertex (PV, SV) variables (peak)
   VertexKFChi2_peak=HConfig.GetTH1D(Name+"_VertexKFChi2_peak","KF Vertex Chi Squared (peak)",50,0,20,"KF vertex #chi^{2}","Events"); peakCollection.push_back(&VertexKFChi2_peak);
   SVPVDsDirAngle_peak=HConfig.GetTH1D(Name+"_SVPVDsDirAngle_peak","SVPVDsDirAngle (peak)",100,0,0.15,"Angle btw #vec{SV}-#vec{PV} and #vec{Ds}, rad","Events"); peakCollection.push_back(&SVPVDsDirAngle_peak);
   NtracksClose_peak=HConfig.GetTH1D(Name+"_NtracksClose_peak","NtracksClose (peak)",8,-0.5,7.5,"Number of tracks close to SV",""); peakCollection.push_back(&NtracksClose_peak);
   NSV_peak=HConfig.GetTH1D(Name+"_NSV_peak","NSV (peak)",8,-0.5,7.5,"N vertices in the Ds cone",""); peakCollection.push_back(&NSV_peak);
   MinMuon_chi2LocalPosition_peak=HConfig.GetTH1D(Name+"_MinMuon_chi2LocalPosition_peak","MinMuon_chi2LocalPosition (peak)",50,0,5,"Min Inner/Outer track #chi^{2}",""); peakCollection.push_back(&MinMuon_chi2LocalPosition_peak);
   MindcaTrackSV_peak=HConfig.GetTH1D(Name+"_MindcaTrackSV_peak","MindcaTrackSV (peak)",50,0,0.1,"Min distance of track to SV","");  peakCollection.push_back(&MindcaTrackSV_peak);
   MinDca_peak=HConfig.GetTH1D(Name+"_MinDca_peak","MinDca (peak)",50,0,0.02,"Min distance between muons and track",""); peakCollection.push_back(&MinDca_peak);
   MinD0SigSV_peak=HConfig.GetTH1D(Name+"_MinD0SigSV_peak","MinD0SigSV (peak)",30,0,30,"Min Transverse Impact significance w.r.t SV",""); peakCollection.push_back(&MinD0SigSV_peak);
   MinD0SigPV_peak=HConfig.GetTH1D(Name+"_MinD0SigPV_peak","MinD0SigPV (peak)",30,0,30,"Min Transverse Impact significance w.r.t PV",""); peakCollection.push_back(&MinD0SigPV_peak);
   MaxVertexPairQuality_peak=HConfig.GetTH1D(Name+"_MaxVertexPairQuality_peak","MaxVertexPairQuality (peak)",30,0,10,"max vertex pair quality","Events"); peakCollection.push_back(&MaxVertexPairQuality_peak);
   MaxdeltaMuZ_peak= HConfig.GetTH1D(Name+"_MaxdeltaMuZ_peak","MaxdeltaMuZ (peak)",30,0,0.6,"Max #Delta z (#mu-#mu), cm",""); peakCollection.push_back(&MaxdeltaMuZ_peak);
   MaxDca_peak=HConfig.GetTH1D(Name+"_MaxDca_peak","MaxDca (peak)",50,0,0.10,"Max distance between muons",""); peakCollection.push_back(&MaxDca_peak);
   MaxD0SigSV_peak=HConfig.GetTH1D(Name+"_MaxD0SigSV_peak","MaxD0SigSV (peak)",30,0,30,"Max Transverse Impact significance w.r.t SV",""); peakCollection.push_back(&MaxD0SigSV_peak);
   MaxD0SigPV_peak=HConfig.GetTH1D(Name+"_MaxD0SigPV_peak","MaxD0SigPV (peak)",30,0,30,"Max Transverse Impact significance w.r.t PV",""); peakCollection.push_back(&MaxD0SigPV_peak);
   FLSignificance_peak=HConfig.GetTH1D(Name+"_FLSignificance_peak","FLSignificance",60,0,60,"PV - SV distance  significance","Events"); peakCollection.push_back(&FLSignificance_peak);

   // Vertex (PV, SV) variables (sideband)
   VertexKFChi2_sideband=HConfig.GetTH1D(Name+"_VertexKFChi2_sideband","KF Vertex Chi Squared (sideband)",50,0,20,"KF vertex #chi^{2}","Events"); sidebandCollection.push_back(&VertexKFChi2_sideband);
   SVPVDsDirAngle_sideband=HConfig.GetTH1D(Name+"_SVPVDsDirAngle_sideband","SVPVDsDirAngle (sideband)",100,0,0.15,"Angle btw #vec{SV}-#vec{PV} and #vec{Ds}, rad","Events"); sidebandCollection.push_back(&SVPVDsDirAngle_sideband);
   NtracksClose_sideband=HConfig.GetTH1D(Name+"_NtracksClose_sideband","NtracksClose (sideband)",8,-0.5,7.5,"Number of tracks close to SV",""); sidebandCollection.push_back(&NtracksClose_sideband);
   NSV_sideband=HConfig.GetTH1D(Name+"_NSV_sideband","NSV (sideband)",8,-0.5,7.5,"N vertices in the Ds cone",""); sidebandCollection.push_back(&NSV_sideband);
   MinMuon_chi2LocalPosition_sideband=HConfig.GetTH1D(Name+"_MinMuon_chi2LocalPosition_sideband","MinMuon_chi2LocalPosition (sideband)",50,0,5,"Min Inner/Outer track #chi^{2}",""); sidebandCollection.push_back(&MinMuon_chi2LocalPosition_sideband);
   MindcaTrackSV_sideband=HConfig.GetTH1D(Name+"_MindcaTrackSV_sideband","MindcaTrackSV (sideband)",50,0,0.1,"Min distance of track to SV","");  sidebandCollection.push_back(&MindcaTrackSV_sideband);
   MinDca_sideband=HConfig.GetTH1D(Name+"_MinDca_sideband","MinDca (sideband)",50,0,0.02,"Min distance between muons and track",""); sidebandCollection.push_back(&MinDca_sideband);
   MinD0SigSV_sideband=HConfig.GetTH1D(Name+"_MinD0SigSV_sideband","MinD0SigSV (sideband)",30,0,30,"Min Transverse Impact significance w.r.t SV",""); sidebandCollection.push_back(&MinD0SigSV_sideband);
   MinD0SigPV_sideband=HConfig.GetTH1D(Name+"_MinD0SigPV_sideband","MinD0SigPV (sideband)",30,0,30,"Min Transverse Impact significance w.r.t PV",""); sidebandCollection.push_back(&MinD0SigPV_sideband);
   MaxVertexPairQuality_sideband=HConfig.GetTH1D(Name+"_MaxVertexPairQuality_sideband","MaxVertexPairQuality (sideband)",30,0,10,"max vertex pair quality","Events"); sidebandCollection.push_back(&MaxVertexPairQuality_sideband);
   MaxdeltaMuZ_sideband= HConfig.GetTH1D(Name+"_MaxdeltaMuZ_sideband","MaxdeltaMuZ (sideband)",30,0,0.6,"Max #Delta z (#mu-#mu), cm",""); sidebandCollection.push_back(&MaxdeltaMuZ_sideband);
   MaxDca_sideband=HConfig.GetTH1D(Name+"_MaxDca_sideband","MaxDca (sideband)",50,0,0.10,"Max distance between muons",""); sidebandCollection.push_back(&MaxDca_sideband);
   MaxD0SigSV_sideband=HConfig.GetTH1D(Name+"_MaxD0SigSV_sideband","MaxD0SigSV (sideband)",30,0,30,"Max Transverse Impact significance w.r.t SV",""); sidebandCollection.push_back(&MaxD0SigSV_sideband);
   MaxD0SigPV_sideband=HConfig.GetTH1D(Name+"_MaxD0SigPV_sideband","MaxD0SigPV (sideband)",30,0,30,"Max Transverse Impact significance w.r.t PV",""); sidebandCollection.push_back(&MaxD0SigPV_sideband);
   FLSignificance_sideband=HConfig.GetTH1D(Name+"_FLSignificance_sideband","FLSignificance (sideband)",60,0,60,"PV - SV distance  significance","Events"); sidebandCollection.push_back(&FLSignificance_sideband);

   // Vertex (PV, SV) variables (validation)
   VertexKFChi2_validation=HConfig.GetTH1D(Name+"_VertexKFChi2_validation","KF Vertex Chi Squared (validation)",50,0,20,"KF vertex #chi^{2}","Events"); validationCollection.push_back(&VertexKFChi2_validation);
   SVPVDsDirAngle_validation=HConfig.GetTH1D(Name+"_SVPVDsDirAngle_validation","SVPVDsDirAngle (validation)",100,0,0.15,"Angle btw #vec{SV}-#vec{PV} and #vec{Ds}, rad","Events"); validationCollection.push_back(&SVPVDsDirAngle_validation);
   NtracksClose_validation=HConfig.GetTH1D(Name+"_NtracksClose_validation","NtracksClose (validation)",8,-0.5,7.5,"Number of tracks close to SV",""); validationCollection.push_back(&NtracksClose_validation);
   NSV_validation=HConfig.GetTH1D(Name+"_NSV_validation","NSV (validation)",8,-0.5,7.5,"N vertices in the Ds cone",""); validationCollection.push_back(&NSV_validation);
   MinMuon_chi2LocalPosition_validation=HConfig.GetTH1D(Name+"_MinMuon_chi2LocalPosition_validation","MinMuon_chi2LocalPosition (validation)",50,0,5,"Min Inner/Outer track #chi^{2}",""); validationCollection.push_back(&MinMuon_chi2LocalPosition_validation);
   MindcaTrackSV_validation=HConfig.GetTH1D(Name+"_MindcaTrackSV_validation","MindcaTrackSV (validation)",50,0,0.1,"Min distance of track to SV","");  validationCollection.push_back(&MindcaTrackSV_validation);
   MinDca_validation=HConfig.GetTH1D(Name+"_MinDca_validation","MinDca (validation)",50,0,0.02,"Min distance between muons and track",""); validationCollection.push_back(&MinDca_validation);
   MinD0SigSV_validation=HConfig.GetTH1D(Name+"_MinD0SigSV_validation","MinD0SigSV (validation)",30,0,30,"Min Transverse Impact significance w.r.t SV",""); validationCollection.push_back(&MinD0SigSV_validation);
   MinD0SigPV_validation=HConfig.GetTH1D(Name+"_MinD0SigPV_validation","MinD0SigPV (validation)",30,0,30,"Min Transverse Impact significance w.r.t PV",""); validationCollection.push_back(&MinD0SigPV_validation);
   MaxVertexPairQuality_validation=HConfig.GetTH1D(Name+"_MaxVertexPairQuality_validation","MaxVertexPairQuality (validation)",30,0,10,"max vertex pair quality","Events"); validationCollection.push_back(&MaxVertexPairQuality_validation);
   MaxdeltaMuZ_validation= HConfig.GetTH1D(Name+"_MaxdeltaMuZ_validation","MaxdeltaMuZ (validation)",30,0,0.6,"Max #Delta z (#mu-#mu), cm",""); validationCollection.push_back(&MaxdeltaMuZ_validation);
   MaxDca_validation=HConfig.GetTH1D(Name+"_MaxDca_validation","MaxDca (validation)",50,0,0.10,"Max distance between muons",""); validationCollection.push_back(&MaxDca_validation);
   MaxD0SigSV_validation=HConfig.GetTH1D(Name+"_MaxD0SigSV_validation","MaxD0SigSV (validation)",30,0,30,"Max Transverse Impact significance w.r.t SV",""); validationCollection.push_back(&MaxD0SigSV_validation);
   MaxD0SigPV_validation=HConfig.GetTH1D(Name+"_MaxD0SigPV_validation","MaxD0SigPV (validation)",30,0,30,"Max Transverse Impact significance w.r.t PV",""); validationCollection.push_back(&MaxD0SigPV_validation);
   FLSignificance_validation=HConfig.GetTH1D(Name+"_FLSignificance_validation","FLSignificance (validation)",60,0,60,"PV - SV distance  significance","Events"); validationCollection.push_back(&FLSignificance_validation);

   // Decay length (validation)
   DecayLength_prompt_validation=HConfig.GetTH1D(Name+"_DecayLength_prompt_validation","Proper Decay Length of Prompt Ds (validation)",20,0,.1,"Proper Decay Length (cm)","Events"); validationCollection.push_back(&DecayLength_prompt_validation);
   DecayLength_non_prompt_validation=HConfig.GetTH1D(Name+"_DecayLength_non_prompt_validation","Proper Decay Length of Non-Prompt Ds (validation)",20,0,.1,"Proper Decay Length (cm)","Events"); validationCollection.push_back(&DecayLength_non_prompt_validation);
   DecayLength_validation=HConfig.GetTH1D(Name+"_DecayLength_validation","Proper Decay Length of Ds (validation)",20,0,.1,"Proper Decay Length (cm)","Events"); validationCollection.push_back(&DecayLength_validation);

   // Decay length (sideband)
   DecayLength_prompt_sideband=HConfig.GetTH1D(Name+"_DecayLength_prompt_sideband","Proper Decay Length of Prompt Ds (sideband)",20,0,.1,"Proper Decay Length (cm)","Events"); sidebandCollection.push_back(&DecayLength_prompt_sideband);
   DecayLength_non_prompt_sideband=HConfig.GetTH1D(Name+"_DecayLength_non_prompt_sideband","Proper Decay Length of Non-Prompt Ds (sideband)",20,0,.1,"Proper Decay Length (cm)","Events"); sidebandCollection.push_back(&DecayLength_non_prompt_sideband);
   DecayLength_sideband=HConfig.GetTH1D(Name+"_DecayLength_sideband","Proper Decay Length of Ds (sideband)",20,0,.1,"Proper Decay Length (cm)","Events"); sidebandCollection.push_back(&DecayLength_sideband);

   // Decay length (peak)
   DecayLength_prompt_peak=HConfig.GetTH1D(Name+"_DecayLength_prompt_peak","Proper Decay Length of Prompt Ds (peak)",20,0,.1,"Proper Decay Length (cm)","Events"); peakCollection.push_back(&DecayLength_prompt_peak);
   DecayLength_non_prompt_peak=HConfig.GetTH1D(Name+"_DecayLength_non_prompt_peak","Proper Decay Length of Non-Prompt Ds (peak)",20,0,.1,"Proper Decay Length (cm)","Events"); peakCollection.push_back(&DecayLength_non_prompt_peak);
   DecayLength_peak=HConfig.GetTH1D(Name+"_DecayLength_peak","Proper Decay Length of Ds (peak)",20,0,.1,"Proper Decay Length (cm)","Events"); peakCollection.push_back(&DecayLength_peak);

   // Invariant mass plot for validation

   TripleMass_peak=HConfig.GetTH1D(Name+"_TripleMass_peak","#mu#mu + track mass",42,1.68,2.1,"Invariant mass of the #mu#mu + track (peak)","Events");
   TripleMass_sideband=HConfig.GetTH1D(Name+"_TripleMass_sideband","#mu#mu + track mass",42,1.68,2.1,"Invariant mass of the #mu#mu + track (sideband)","Events");
   TripleMass_validation=HConfig.GetTH1D(Name+"_TripleMass_validation","#mu#mu + track mass",42,1.68,2.1,"Invariant mass of the #mu#mu + track (validation)","Events");

   peakCollection.push_back(&TripleMass_peak);
   sidebandCollection.push_back(&TripleMass_sideband);
   validationCollection.push_back(&TripleMass_validation);

   Selection::ConfigureHistograms(); //do not remove
   HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}

void  DsPhiPeak::Store_ExtraDist(){ 

   // Unscaled 
   Extradist1d.push_back(&TripleMass);
   Extradist2d.push_back(&PhiMassVsDsMass);
   // Charge MisId and trigger efficiency
   Extradist1d.push_back(&ChargeMisId);
   Extradist1d.push_back(&TriggerEfficiency);  

   ////////////////////
   // validation plots  
   ////////////////////

   for (int j=0; j<peakCollection.size(); ++j){
      Extradist1d.push_back(peakCollection.at(j));
      Extradist1d.push_back(sidebandCollection.at(j));
      Extradist1d.push_back(validationCollection.at(j));
   }

   //////////////////////////////////////////////////////////////////////////////////////////////////////
   // Here you must push back all analysis histograms, otherwise they wont be propagated to the output //
   //////////////////////////////////////////////////////////////////////////////////////////////////////

}

void  DsPhiPeak::doEvent(){ 
   unsigned int t;
   int id(Ntp->GetMCID());
   if (id==1){
      if(Ntp->WhichEra(2018).Contains("RunA")){ nPeak = 0.89; nSidebands = 1.0;} 
      if(Ntp->WhichEra(2018).Contains("RunB")){ nPeak = 0.90; nSidebands = 1.0;}
      if(Ntp->WhichEra(2018).Contains("RunC")){ nPeak = 0.90; nSidebands = 1.0;}
      if(Ntp->WhichEra(2018).Contains("RunD")){ nPeak = 0.90; nSidebands = 1.0;}
      
   }

   bool DEBUG = false;
   bool HLTOk(false);
   bool L1Ok(false);
   bool DoubleMuFired(false);
   bool TripleMuFired(false);

   double wobs=1;
   double w;
   double w_normalization=1.0;

   // Isolation variables
   double mindca_iso05 = 99.0;
   double mindca_iso = 99.0;
   double mindca_ds = 99.0;

   double sumPtTracks_mu1 = 0;
   double sumPtTracks_mu3 = 0;
   double sumPtTracks_mu2 = 0;

   double sumPtTracks_ds = 0.;
   double sumPtTracks_iso05 = 0.;

   int nTracks_iso05 = 0;
   int nTracks_ds = 0;

   // Initialize all cut values
   value.at(TriggerOk)=0;
   value.at(TwoMuTrkCandidate)=0;
   value.at(Mu1PtCut)=0;
   value.at(Mu2PtCut)=0;
   value.at(OSMuons)=1.0;
   value.at(TrackPtCut)=0;
   value.at(NTrackHits)=0;
   value.at(MuonID)=0;
   value.at(MuMuMassCut)=0;
   value.at(TriggerMatchMu1)=1.0;
   value.at(TriggerMatchMu2)=1.0;
   value.at(TriggerMatchTrack)=1.0;
   value.at(DsMassCut)=0;
   value.at(GenMatch)=1;   
   
   random_num = rndm.Rndm();

   if(!Ntp->isData()){
      w = w_normalization*(puWeights->GetBinContent(Ntp->TruthNumberOfInteraction())); // Weight MC according to truth number of vertices
   } 
   //  No weights to data
   else{w=1;}

   if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}

   // Apply Selection
   for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
      TString HLT = Ntp->HLTName(iTrigger);
      if(HLT.Contains("HLT_DoubleMu3_Trk_Tau3mu_v") && Ntp->HLTDecision(iTrigger)  ) { HLTOk = true;}
   }

   for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
      TString L1TriggerName = Ntp->L1Name(il1);
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4") && Ntp->L1Decision(il1)) { DoubleMuFired = true; }
      if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9") && Ntp->L1Decision(il1)) { TripleMuFired = true; }
      if( id!=1 && random_num>0.3516 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMuFired = true; }
      if( id==1 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMuFired = true; }
   }

   if (DoubleMuFired || TripleMuFired) L1Ok = true; // Remove all the events seeded by triple mu
   if (L1Ok && HLTOk) value.at(TriggerOk) = true; // HLT requirement
   else value.at(TriggerOk) = 0;

   value.at(TwoMuTrkCandidate) = Ntp->NTwoMuonsTrack(); // Number of two muons + track candidates
   pass.at(TriggerOk) = ( value.at(TriggerOk) == cut.at(TriggerOk) );
   pass.at(TwoMuTrkCandidate) = ( value.at(TwoMuTrkCandidate) >= cut.at(TwoMuTrkCandidate) );
   // Selection of the best candidate

   vector<unsigned int> selectedIndices;
   vector<unsigned int> candidateRank;
   unsigned int final_idx = 0;
   double minChiSq = 999.0;
   bool status; 

   if (Ntp->NTwoMuonsTrack()>0){
      for (size_t j=0; j<Ntp->NTwoMuonsTrack(); ++j){ 

         unsigned int Muon1_idx = Ntp->TwoMuonsTrackMuonIndices(j).at(0);
         unsigned int Muon2_idx = Ntp->TwoMuonsTrackMuonIndices(j).at(1);
         unsigned int Track_idx = Ntp->TwoMuonsTrackTrackIndex(j).at(0);

         unsigned int Mu1_pt_idx = (Ntp->Muon_P4(Muon1_idx).Pt() > Ntp->Muon_P4(Muon2_idx).Pt() ? Muon1_idx:Muon2_idx);
         unsigned int Mu2_pt_idx = (Ntp->Muon_P4(Muon1_idx).Pt() > Ntp->Muon_P4(Muon2_idx).Pt() ? Muon2_idx:Muon1_idx);

         TLorentzVector mu1_p4 = Ntp->Muon_P4(Mu1_pt_idx);
         TLorentzVector mu2_p4 = Ntp->Muon_P4(Mu2_pt_idx);
         TLorentzVector track_p4 = Ntp->Track_P4(Track_idx);

         double mumutrkmass = (mu1_p4+mu2_p4+track_p4).M();

         value.at(Mu1PtCut) = mu1_p4.Pt();
         value.at(Mu2PtCut) = mu2_p4.Pt();
         value.at(MuonID) = (Ntp->Muon_isGlobalMuon(Mu1_pt_idx) && Ntp->Muon_isGlobalMuon(Mu2_pt_idx) );
         value.at(OSMuons) = abs(Ntp->Muon_charge(Mu1_pt_idx)+Ntp->Muon_charge(Mu2_pt_idx));

         float PhiMass = (mu1_p4+mu2_p4).M(); 
         value.at(MuMuMassCut) = PhiMass;
         value.at(TrackPtCut) = track_p4.Pt();
         value.at(NTrackHits) = Ntp->Track_numberOfValidHits(Track_idx);
         //value.at(ChiSqCut) = Ntp->TwoMuonsTrack_SV_Chi2(j);
         value.at(TriggerMatchMu1) = Ntp->TwoMuonsTrack_TriggerMatch_dR(j).at(0);
         value.at(TriggerMatchMu2) = Ntp->TwoMuonsTrack_TriggerMatch_dR(j).at(1);
         value.at(TriggerMatchTrack) = Ntp->TwoMuonsTrack_TriggerMatch_dR(j).at(2);
         value.at(DsMassCut) = mumutrkmass;
         if (id==1) value.at(GenMatch) = 0;
         else value.at(GenMatch) = Ntp->DsGenMatch(j); 


         pass.at(OSMuons) = value.at(OSMuons) < cut.at(OSMuons);
         pass.at(Mu1PtCut) = (value.at(Mu1PtCut) >= cut.at(Mu1PtCut) );
         pass.at(Mu2PtCut) = (value.at(Mu2PtCut) >= cut.at(Mu2PtCut) );
         pass.at(MuonID) = ( value.at(MuonID) == cut.at(MuonID) );
         pass.at(MuMuMassCut) = ( ( value.at(MuMuMassCut)< PhiMassHigh) && (value.at(MuMuMassCut) >= PhiMassLow) );
         pass.at(TrackPtCut) = ( value.at(TrackPtCut) > cut.at(TrackPtCut) );
         pass.at(NTrackHits) = ( value.at(NTrackHits) > cut.at(NTrackHits) );
         //pass.at(ChiSqCut) = ( value.at(ChiSqCut) < cut.at(ChiSqCut) );
         //pass.at(TriggerMatchMu1) = ( value.at(TriggerMatchMu1) < cut.at(TriggerMatchMu1) );
         //pass.at(TriggerMatchMu2) = ( value.at(TriggerMatchMu2) < cut.at(TriggerMatchMu2) );
         //pass.at(TriggerMatchTrack) = ( value.at(TriggerMatchTrack) < cut.at(TriggerMatchTrack) );
         pass.at(TriggerMatchMu1) = true;
         pass.at(TriggerMatchMu2) = true;
         pass.at(TriggerMatchTrack) = true;
         pass.at(GenMatch) = ( value.at(GenMatch) < cut.at(GenMatch) );
         pass.at(DsMassCut) = ( value.at(DsMassCut) < dsMassMax && value.at(DsMassCut) >= dsMassMin );

         unsigned int score = 0;
         for (unsigned int k=0; k<NCuts; ++k) if (pass.at(k)) score++;

         if (score==NCuts) selectedIndices.push_back(j);
         candidateRank.push_back(score);
      }

      if (selectedIndices.size()==0) final_idx=std::distance(candidateRank.begin(), std::max_element(candidateRank.begin(), candidateRank.end()));
      else{
         for (size_t j=0; j<selectedIndices.size(); ++j){
            int _ = selectedIndices.at(j);
            double tmpchi = Ntp->TwoMuonsTrack_SV_Chi2(_);
            if (tmpchi<minChiSq) {
               minChiSq = tmpchi;
               final_idx = _;
            }
         }
      }

      unsigned int Muon1_idx = Ntp->TwoMuonsTrackMuonIndices(final_idx).at(0);
      unsigned int Muon2_idx = Ntp->TwoMuonsTrackMuonIndices(final_idx).at(1);
      unsigned int Track_idx = Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0);

      unsigned int Mu1_pt_idx = (Ntp->Muon_P4(Muon1_idx).Pt() > Ntp->Muon_P4(Muon2_idx).Pt() ? Muon1_idx:Muon2_idx);
      unsigned int Mu2_pt_idx = (Ntp->Muon_P4(Muon1_idx).Pt() > Ntp->Muon_P4(Muon2_idx).Pt() ? Muon2_idx:Muon1_idx);

      TLorentzVector mu1_p4 = Ntp->Muon_P4(Mu1_pt_idx);
      TLorentzVector mu2_p4 = Ntp->Muon_P4(Mu2_pt_idx);
      TLorentzVector track_p4 = Ntp->Track_P4(Track_idx);

      double mumutrkmass = (mu1_p4+mu2_p4+track_p4).M();

      value.at(Mu1PtCut) = mu1_p4.Pt();
      value.at(Mu2PtCut) = mu2_p4.Pt();
      value.at(MuonID) = (Ntp->Muon_isGlobalMuon(Mu1_pt_idx) && Ntp->Muon_isGlobalMuon(Mu2_pt_idx) );
      value.at(OSMuons) = abs(Ntp->Muon_charge(Mu1_pt_idx)+Ntp->Muon_charge(Mu2_pt_idx));

      float PhiMass = (mu1_p4+mu2_p4).M(); 
      value.at(MuMuMassCut) = PhiMass;
      value.at(TrackPtCut) = track_p4.Pt();
      value.at(NTrackHits) = Ntp->Track_numberOfValidHits(Track_idx);
      //value.at(ChiSqCut) = Ntp->TwoMuonsTrack_SV_Chi2(final_idx);
      value.at(TriggerMatchMu1) = Ntp->TwoMuonsTrack_TriggerMatch_dR(final_idx).at(0);
      value.at(TriggerMatchMu2) = Ntp->TwoMuonsTrack_TriggerMatch_dR(final_idx).at(1);
      value.at(TriggerMatchTrack) = Ntp->TwoMuonsTrack_TriggerMatch_dR(final_idx).at(2);
      value.at(DsMassCut) = mumutrkmass;
      if (id==1) value.at(GenMatch) = 0;
      else value.at(GenMatch) = Ntp->DsGenMatch(final_idx); 

      pass.at(OSMuons) = value.at(OSMuons) < cut.at(OSMuons);
      pass.at(Mu1PtCut) = (value.at(Mu1PtCut) >= cut.at(Mu1PtCut) );
      pass.at(Mu2PtCut) = (value.at(Mu2PtCut) >= cut.at(Mu2PtCut) );
      pass.at(MuonID) = ( value.at(MuonID) == cut.at(MuonID) );
      pass.at(MuMuMassCut) = ( ( value.at(MuMuMassCut)< PhiMassHigh) && (value.at(MuMuMassCut) >= PhiMassLow) );
      pass.at(TrackPtCut) = ( value.at(TrackPtCut) > cut.at(TrackPtCut) );
      pass.at(NTrackHits) = ( value.at(NTrackHits) > cut.at(NTrackHits) );
      //pass.at(ChiSqCut) = ( value.at(ChiSqCut) < cut.at(ChiSqCut) );
      //pass.at(TriggerMatchMu1) = ( value.at(TriggerMatchMu1) < cut.at(TriggerMatchMu1) );
      //pass.at(TriggerMatchMu2) = ( value.at(TriggerMatchMu2) < cut.at(TriggerMatchMu2) );
      //pass.at(TriggerMatchTrack) = ( value.at(TriggerMatchTrack) < cut.at(TriggerMatchTrack) );
      pass.at(TriggerMatchMu1) = true;
      pass.at(TriggerMatchMu2) = true;
      pass.at(TriggerMatchTrack) = true;
      pass.at(GenMatch) = ( value.at(GenMatch) < cut.at(GenMatch) );
      pass.at(DsMassCut) = ( value.at(DsMassCut) >= dsMassMin && value.at(DsMassCut) < dsMassMax );

      // charge misid
      if(passAllBut(OSMuons)){
         if (value.at(OSMuons)>0.1) ChargeMisId.at(t).Fill(mumutrkmass, w);
      } 

      // trigger matching 
      vector<unsigned int> triggerMatchCuts;
      triggerMatchCuts.push_back(TriggerMatchMu1);
      triggerMatchCuts.push_back(TriggerMatchMu2);
      triggerMatchCuts.push_back(TriggerMatchTrack);

      if(passAllBut(triggerMatchCuts)){
         if ( !pass.at(TriggerMatchMu1) || !pass.at(TriggerMatchMu2) || !pass.at(TriggerMatchTrack) )
            TriggerEfficiency.at(t).Fill(mumutrkmass, w);
      }     
   }
   status = AnalysisCuts(t,w,wobs);


   /*
    * ----------------------
    * Trigger Matching Test
    *-----------------------
    vector<unsigned int> triggerMatchCuts;
    triggerMatchCuts.push_back(TriggerMatchMu1);
    triggerMatchCuts.push_back(TriggerMatchMu2);
    triggerMatchCuts.push_back(TriggerMatchTrack);

    if (passAllBut(triggerMatchCuts) && (Ntp->TwoMuonsTrack_TriggerMatch_dR(final_idx).at(2)>=0.03)){
    cout<<"============ Event ============"<<endl;
    cout<<"Minimum chi square found: "<<Ntp->TwoMuonsTrack_SV_Chi2(final_idx)<<endl;
    cout<<"Number of dsphipi candidates: "<<Ntp->NTwoMuonsTrack()<<endl;
    cout<<"Number of tracks in the event: "<<Ntp->NTracks()<<endl;
    cout<<"Trying to match other candidates..."<<endl;
    for (size_t j=0; j<Ntp->NTwoMuonsTrack(); ++j){
    if (j==final_idx) continue;
    if ( (Ntp->TwoMuonsTrack_TriggerMatch_dR(j).at(0)<0.03) &&
    (Ntp->TwoMuonsTrack_TriggerMatch_dR(j).at(1)<0.03) &&
    (Ntp->TwoMuonsTrack_TriggerMatch_dR(j).at(2)<0.03)){
    cout<<"Candidate found with chi square "<<Ntp->TwoMuonsTrack_SV_Chi2(j)<<endl;
    cout<<( std::find(selectedIndices.begin(), selectedIndices.end(), j)!=selectedIndices.end() ? "Selected":"Not selected")<<endl;
    }
    }
    }
    */

   if(status){

      unsigned int mu1_idx = Ntp->TwoMuonsTrackMuonIndices(final_idx).at(0);
      unsigned int mu2_idx = Ntp->TwoMuonsTrackMuonIndices(final_idx).at(1);
      Track_idx = Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0);

      Muon_1_idx = (Ntp->Muon_P4(mu1_idx).Pt()>Ntp->Muon_P4(mu1_idx).Pt() ? mu1_idx:mu2_idx); 
      Muon_2_idx = (Ntp->Muon_P4(mu2_idx).Pt()>Ntp->Muon_P4(mu2_idx).Pt() ? mu2_idx:mu1_idx); 
      Track_idx = Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0);

      TLorentzVector mu1_p4 = Ntp->Muon_P4(Muon_1_idx);
      TLorentzVector mu2_p4 = Ntp->Muon_P4(Muon_2_idx);
      TLorentzVector track_p4 = Ntp->Track_P4(Track_idx);

      TLorentzVector DsLV = mu1_p4 + mu2_p4 + track_p4;
      TLorentzVector DsRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0,true)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1,true)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2,true);

      double phimass = (mu1_p4+mu2_p4).M();
      double dsmass = (mu1_p4+mu2_p4+track_p4).M();
      TVector3 SVPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(final_idx,true),Ntp->Vertex_MatchedPrimaryVertex(final_idx,true));
      double DecayL = SVPV.Mag()*dsmass/DsLV.E();

      int Nvertices(0);

      for(unsigned int l=0; l < Ntp->NSecondaryVertices(); l++){
         TVector3 SVsignalPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(final_idx,true),Ntp->Vertex_MatchedPrimaryVertex(final_idx,true));
         TVector3 SVfakePV = Ntp->SVPVDirection(Ntp->SecondaryVertexPosition(l),Ntp->Vertex_MatchedPrimaryVertex(final_idx,true));
         if(SVfakePV.DeltaR(SVsignalPV) < 1 && (Ntp->Vertex_Signal_KF_pos(final_idx,true) - Ntp->SecondaryVertexPosition(l)).Mag() > 0.05){
            Nvertices++;
         }
      }

      float MinD0Significance = std::min({Ntp->Vertex_d0sig_reco(final_idx,0,true),
            Ntp->Vertex_d0sig_reco(final_idx,1,true),
            Ntp->Vertex_d0sig_reco(final_idx,2,true)});

      float MaxD0Significance = std::max({Ntp->Vertex_d0sig_reco(final_idx,0,true),
            Ntp->Vertex_d0sig_reco(final_idx,1,true),
            Ntp->Vertex_d0sig_reco(final_idx,2,true)});

      float MinD0SVSignificance = std::min({Ntp->Vertex_d0sigSV_reco(final_idx,0,true),
            Ntp->Vertex_d0sigSV_reco(final_idx,1,true),
            Ntp->Vertex_d0sigSV_reco(final_idx,2,true)});

      float MaxD0SVSignificance = std::max({Ntp->Vertex_d0sigSV_reco(final_idx,0,true),
            Ntp->Vertex_d0sigSV_reco(final_idx,1,true),
            Ntp->Vertex_d0sigSV_reco(final_idx,2,true)});  

      double maxMuondR = std::max({mu1_p4.DeltaR(DsLV), mu2_p4.DeltaR(DsLV), track_p4.DeltaR(DsLV)});
      double minMuonPt = std::min({mu1_p4.Pt(), mu2_p4.Pt(), track_p4.Pt()});

      // Isolation algorithm
      for (int it=0; it<Ntp->NIsolationTrack(final_idx, true); it++){
         double dxy_track = Ntp->IsolationTrack_dxySV(final_idx, it, true);
         double dz_track = Ntp->IsolationTrack_dzSV(final_idx, it, true);
         TLorentzVector TrackLV = Ntp->IsolationTrack_p4(final_idx, it, true);
         double dca_fv = TMath::Sqrt(pow(dxy_track, 2)+ pow(dz_track, 2));

         double dr_ds = DsLV.DeltaR(TrackLV); 
         double dr_mu1 = mu1_p4.DeltaR(TrackLV);
         double dr_mu2 = mu2_p4.DeltaR(TrackLV);
         double dr_trk = track_p4.DeltaR(TrackLV);

         // Isolation 1
         if ( dca_fv<0.5 && TrackLV.Pt()<0.33*minMuonPt && dr_ds<3*maxMuondR ){
            sumPtTracks_ds += TrackLV.Pt();
            nTracks_ds++;
            if (dca_fv < mindca_ds) mindca_ds = dca_fv;
         }

         // Isolation 2
         if (TrackLV.Pt()<1.0) {
            continue;
         }

         if (dca_fv < mindca_iso) mindca_iso = dca_fv;

         // Isolation 3 (within dR = 0.5 of tau)
         if (dr_ds<0.5 && dca_fv<0.5){
            sumPtTracks_iso05 += TrackLV.Pt();
            nTracks_iso05++;
            if(dca_fv<mindca_iso05) mindca_iso05 = dca_fv;
         }
      }

      // Relative Pt calculation
      double relPt_iso05 = sumPtTracks_ds/DsLV.Pt();

      // Fill histograms
      TripleMass.at(t).Fill(dsmass, w);
      PhiMassVsDsMass.at(t).Fill(phimass, dsmass, w);


      //* ----------------
      // Validation block
      //* ----------------

      // fill sideband mass
      if (dsmass<sidebandDsMax && dsmass>=sidebandDsMin){

         Muon1_TriggerMatchdR_sideband.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(final_idx)).at(0),w);
         Muon2_TriggerMatchdR_sideband.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(final_idx)).at(1),w);
         Track_TriggerMatchdR_sideband.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(final_idx)).at(2),w);
         Muon1_isGlobal_sideband.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_1_idx),w);
         Muon2_isGlobal_sideband.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_2_idx),w);
         Muon1_isStandAlone_sideband.at(t).Fill(Ntp->Muon_isStandAloneMuon(Muon_1_idx),w);
         Muon2_isStandAlone_sideband.at(t).Fill(Ntp->Muon_isStandAloneMuon(Muon_2_idx),w);
         Muon1_isTracker_sideband.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_1_idx),w);
         Muon2_isTracker_sideband.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_2_idx),w);
         Muon1_isCalo_sideband.at(t).Fill(Ntp->Muon_isCaloMuon(Muon_1_idx),w);
         Muon2_isCalo_sideband.at(t).Fill(Ntp->Muon_isCaloMuon(Muon_2_idx),w);
         Muon1_isIsolationValid_sideband.at(t).Fill(Ntp->Muon_isIsolationValid(Muon_1_idx),w);
         Muon2_isIsolationValid_sideband.at(t).Fill(Ntp->Muon_isIsolationValid(Muon_2_idx),w);
         NVtx_woPUWeights_sideband.at(t).Fill(Ntp->NVtx(),w_normalization); 
         NVtx_sideband.at(t).Fill(Ntp->NVtx(),w);
         PhiMass_sideband.at(t).Fill(phimass, w);
         TripleMass_sideband.at(t).Fill(dsmass, w);

         Track_Pt_sideband.at(t).Fill(Ntp->Track_P4(Track_idx).Pt(),w);
         Track_Eta_sideband.at(t).Fill(Ntp->Track_P4(Track_idx).Eta(),w);
         Track_Phi_sideband.at(t).Fill(Ntp->Track_P4(Track_idx).Phi(),w);
         Track_vx_sideband.at(t).Fill(Ntp->Track_Poca(Track_idx).X(),w);
         Track_vy_sideband.at(t).Fill(Ntp->Track_Poca(Track_idx).Y(),w);
         Track_vz_sideband.at(t).Fill(Ntp->Track_Poca(Track_idx).Z(),w);
         Track_normalizedChi2_sideband.at(t).Fill(Ntp->Track_normalizedChi2(Track_idx),w);
         Track_numberOfValidHits_sideband.at(t).Fill(Ntp->Track_numberOfValidHits(Track_idx),w);
         Track_charge_sideband.at(t).Fill(Ntp->Track_charge(Track_idx),w);
         Track_dxy_sideband.at(t).Fill(Ntp->Track_dxy(Track_idx),w);
         Track_dz_sideband.at(t).Fill(Ntp->Track_dz(Track_idx),w);

         Muon1_Pt_sideband.at(t).Fill(Ntp->Muon_P4(Muon_1_idx).Pt(),w);
         Muon1_Eta_sideband.at(t).Fill(Ntp->Muon_P4(Muon_1_idx).Eta(),w);
         Muon1_Phi_sideband.at(t).Fill(Ntp->Muon_P4(Muon_1_idx).Phi(),w);
         Muon1_vx_sideband.at(t).Fill(Ntp->Muon_Poca(Muon_1_idx).X(),w);
         Muon1_vy_sideband.at(t).Fill(Ntp->Muon_Poca(Muon_1_idx).Y(),w);
         Muon1_vz_sideband.at(t).Fill(Ntp->Muon_Poca(Muon_1_idx).Z(),w);

         Muon2_Pt_sideband.at(t).Fill(Ntp->Muon_P4(Muon_2_idx).Pt(),w);
         Muon2_Eta_sideband.at(t).Fill(Ntp->Muon_P4(Muon_2_idx).Eta(),w);
         Muon2_Phi_sideband.at(t).Fill(Ntp->Muon_P4(Muon_2_idx).Phi(),w);
         Muon2_vx_sideband.at(t).Fill(Ntp->Muon_Poca(Muon_2_idx).X(),w);
         Muon2_vy_sideband.at(t).Fill(Ntp->Muon_Poca(Muon_2_idx).Y(),w);
         Muon2_vz_sideband.at(t).Fill(Ntp->Muon_Poca(Muon_2_idx).Z(),w);

         DimuondR_sideband.at(t).Fill(mu1_p4.DeltaR(mu2_p4),w);
         Muon1TrkdR_sideband.at(t).Fill(mu1_p4.DeltaR(track_p4),w);
         Muon2TrkdR_sideband.at(t).Fill(mu2_p4.DeltaR(track_p4),w);

         VertexKFChi2_sideband.at(t).Fill(Ntp->Vertex_signal_KF_Chi2(final_idx,true),w);
         SVPVDsDirAngle_sideband.at(t).Fill(SVPV.Angle(DsLV.Vect()),w);
         NtracksClose_sideband.at(t).Fill(nTracks_ds,w);
         NSV_sideband.at(t).Fill(Nvertices,w);
         MinMuon_chi2LocalPosition_sideband.at(t).Fill(std::min({Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_1_idx),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_2_idx)}),w);
         MinDca_sideband.at(t).Fill(std::min({Ntp->Vertex_DCA12(final_idx,true),Ntp->Vertex_DCA23(final_idx,true),Ntp->Vertex_DCA31(final_idx,true)}),w);
         MinD0SigSV_sideband.at(t).Fill(MinD0SVSignificance,w);
         MinD0SigPV_sideband.at(t).Fill(MinD0Significance,w);
         MaxVertexPairQuality_sideband.at(t).Fill(  std::max({Ntp->Vertex_pair_quality(final_idx,true),
                  Ntp->Vertex_pair_quality(final_idx,true),
                  Ntp->Vertex_pair_quality(final_idx,true)}),w);
         MaxdeltaMuZ_sideband.at(t).Fill( std::max({fabs(Ntp->Muon_Poca(Muon_1_idx).Z()  - Ntp->Muon_Poca(Muon_2_idx).Z()),
                  fabs(Ntp->Muon_Poca(Muon_1_idx).Z()  - Ntp->Track_Poca(Track_idx).Z()),
                  fabs(Ntp->Muon_Poca(Muon_2_idx).Z()  - Ntp->Track_Poca(Track_idx).Z())}),w);
         MaxDca_sideband.at(t).Fill(std::max({Ntp->Vertex_DCA12(final_idx,true),Ntp->Vertex_DCA23(final_idx,true),Ntp->Vertex_DCA31(final_idx,true)}),w);
         MaxD0SigSV_sideband.at(t).Fill(MaxD0SVSignificance,w);
         MaxD0SigPV_sideband.at(t).Fill(MaxD0Significance,w);
         FLSignificance_sideband.at(t).Fill(( Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(final_idx,true),
                     Ntp->Vertex_PrimaryVertex_Covariance(final_idx,true),
                     Ntp->Vertex_Signal_KF_pos(final_idx,true),
                     Ntp->Vertex_Signal_KF_Covariance(final_idx,true))),w);
         DecayLength_sideband.at(t).Fill(DecayL,w);
         if (id!=1){
            if(Ntp->isPromptDs()){
               DecayLength_prompt_sideband.at(t).Fill(DecayL,w);
            }
            else DecayLength_non_prompt_sideband.at(t).Fill(DecayL,w);
         }
      }

      // fill peak mass
      else if (dsmass<peakDsMax && dsmass>=peakDsMin){

         Muon1_TriggerMatchdR_peak.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(final_idx)).at(0),w);
         Muon2_TriggerMatchdR_peak.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(final_idx)).at(1),w);
         Track_TriggerMatchdR_peak.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(final_idx)).at(2),w);
         Muon1_isGlobal_peak.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_1_idx),w);
         Muon2_isGlobal_peak.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_2_idx),w);
         Muon1_isStandAlone_peak.at(t).Fill(Ntp->Muon_isStandAloneMuon(Muon_1_idx),w);
         Muon2_isStandAlone_peak.at(t).Fill(Ntp->Muon_isStandAloneMuon(Muon_2_idx),w);
         Muon1_isTracker_peak.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_1_idx),w);
         Muon2_isTracker_peak.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_2_idx),w);
         Muon1_isCalo_peak.at(t).Fill(Ntp->Muon_isCaloMuon(Muon_1_idx),w);
         Muon2_isCalo_peak.at(t).Fill(Ntp->Muon_isCaloMuon(Muon_2_idx),w);
         Muon1_isIsolationValid_peak.at(t).Fill(Ntp->Muon_isIsolationValid(Muon_1_idx),w);
         Muon2_isIsolationValid_peak.at(t).Fill(Ntp->Muon_isIsolationValid(Muon_2_idx),w);
         NVtx_woPUWeights_peak.at(t).Fill(Ntp->NVtx(),w_normalization); 
         NVtx_peak.at(t).Fill(Ntp->NVtx(),w);
         PhiMass_peak.at(t).Fill(phimass, w);
         TripleMass_peak.at(t).Fill(dsmass, w);

         Track_Pt_peak.at(t).Fill(Ntp->Track_P4(Track_idx).Pt(),w);
         Track_Eta_peak.at(t).Fill(Ntp->Track_P4(Track_idx).Eta(),w);
         Track_Phi_peak.at(t).Fill(Ntp->Track_P4(Track_idx).Phi(),w);
         Track_vx_peak.at(t).Fill(Ntp->Track_Poca(Track_idx).X(),w);
         Track_vy_peak.at(t).Fill(Ntp->Track_Poca(Track_idx).Y(),w);
         Track_vz_peak.at(t).Fill(Ntp->Track_Poca(Track_idx).Z(),w);
         Track_normalizedChi2_peak.at(t).Fill(Ntp->Track_normalizedChi2(Track_idx),w);
         Track_numberOfValidHits_peak.at(t).Fill(Ntp->Track_numberOfValidHits(Track_idx),w);
         Track_charge_peak.at(t).Fill(Ntp->Track_charge(Track_idx),w);
         Track_dxy_peak.at(t).Fill(Ntp->Track_dxy(Track_idx),w);
         Track_dz_peak.at(t).Fill(Ntp->Track_dz(Track_idx),w);

         Muon1_Pt_peak.at(t).Fill(Ntp->Muon_P4(Muon_1_idx).Pt(),w);
         Muon1_Eta_peak.at(t).Fill(Ntp->Muon_P4(Muon_1_idx).Eta(),w);
         Muon1_Phi_peak.at(t).Fill(Ntp->Muon_P4(Muon_1_idx).Phi(),w);
         Muon1_vx_peak.at(t).Fill(Ntp->Muon_Poca(Muon_1_idx).X(),w);
         Muon1_vy_peak.at(t).Fill(Ntp->Muon_Poca(Muon_1_idx).Y(),w);
         Muon1_vz_peak.at(t).Fill(Ntp->Muon_Poca(Muon_1_idx).Z(),w);

         Muon2_Pt_peak.at(t).Fill(Ntp->Muon_P4(Muon_2_idx).Pt(),w);
         Muon2_Eta_peak.at(t).Fill(Ntp->Muon_P4(Muon_2_idx).Eta(),w);
         Muon2_Phi_peak.at(t).Fill(Ntp->Muon_P4(Muon_2_idx).Phi(),w);
         Muon2_vx_peak.at(t).Fill(Ntp->Muon_Poca(Muon_2_idx).X(),w);
         Muon2_vy_peak.at(t).Fill(Ntp->Muon_Poca(Muon_2_idx).Y(),w);
         Muon2_vz_peak.at(t).Fill(Ntp->Muon_Poca(Muon_2_idx).Z(),w);

         DimuondR_peak.at(t).Fill(mu1_p4.DeltaR(mu2_p4),w);
         Muon1TrkdR_peak.at(t).Fill(mu1_p4.DeltaR(track_p4),w);
         Muon2TrkdR_peak.at(t).Fill(mu2_p4.DeltaR(track_p4),w);

         VertexKFChi2_peak.at(t).Fill(Ntp->Vertex_signal_KF_Chi2(final_idx,true),w);
         SVPVDsDirAngle_peak.at(t).Fill(SVPV.Angle(DsLV.Vect()),w);
         NtracksClose_peak.at(t).Fill(nTracks_ds,w);
         NSV_peak.at(t).Fill(Nvertices,w);
         MinMuon_chi2LocalPosition_peak.at(t).Fill(std::min({Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_1_idx),
                  Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_2_idx)}),w);
         MinDca_peak.at(t).Fill(std::min({Ntp->Vertex_DCA12(final_idx,true),Ntp->Vertex_DCA23(final_idx,true),Ntp->Vertex_DCA31(final_idx,true)}),w);
         MinD0SigSV_peak.at(t).Fill(MinD0SVSignificance,w);
         MinD0SigPV_peak.at(t).Fill(MinD0Significance,w);
         MaxVertexPairQuality_peak.at(t).Fill(  std::max({Ntp->Vertex_pair_quality(final_idx,true),
                  Ntp->Vertex_pair_quality(final_idx,true),
                  Ntp->Vertex_pair_quality(final_idx,true)}),w);
         MaxdeltaMuZ_peak.at(t).Fill( std::max({fabs(Ntp->Muon_Poca(Muon_1_idx).Z()  - Ntp->Muon_Poca(Muon_2_idx).Z()),
                  fabs(Ntp->Muon_Poca(Muon_1_idx).Z()  - Ntp->Track_Poca(Track_idx).Z()),
                  fabs(Ntp->Muon_Poca(Muon_2_idx).Z()  - Ntp->Track_Poca(Track_idx).Z())}),w);
         MaxDca_peak.at(t).Fill(std::max({Ntp->Vertex_DCA12(final_idx,true),Ntp->Vertex_DCA23(final_idx,true),Ntp->Vertex_DCA31(final_idx,true)}),w);
         MaxD0SigSV_peak.at(t).Fill(MaxD0SVSignificance,w);
         MaxD0SigPV_peak.at(t).Fill(MaxD0Significance,w);
         FLSignificance_peak.at(t).Fill(( Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(final_idx,true),
                     Ntp->Vertex_PrimaryVertex_Covariance(final_idx,true),
                     Ntp->Vertex_Signal_KF_pos(final_idx,true),
                     Ntp->Vertex_Signal_KF_Covariance(final_idx,true))),w);
         DecayLength_peak.at(t).Fill(DecayL,w);
         if (id!=1){    
            if(Ntp->isPromptDs()){
               DecayLength_prompt_peak.at(t).Fill(DecayL,w);
            }
            else DecayLength_non_prompt_peak.at(t).Fill(DecayL,w);
         }
      }
   }
}

void  DsPhiPeak::Finish(){   

   if (mode==ANALYSIS){
      int id(Ntp->GetMCID());
      if (id==1) {   
         for ( unsigned int j=0; j<validationCollection.size(); ++j){
            (validationCollection.at(j)->at(0)).Add(&(peakCollection.at(j)->at(0)),1.0);
            (validationCollection.at(j)->at(0)).Add(&(sidebandCollection.at(j)->at(0)),-(nPeak/nSidebands));
            cout<<"Ratio: "<<nPeak/nSidebands<<endl;
         }
      }
      else if (id==30){
         for ( unsigned int j=0; j<validationCollection.size(); ++j){
            validationCollection.at(j)->at(1).Add(&(peakCollection.at(j)->at(1)),1.0);
         }
      }
   }
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
   
   if(mode == RECONSTRUCT){
   //   for(unsigned int i=1; i<  Nminus0.at(0).size(); i++){
   //      int id(Ntp->GetMCID());
   //      double scaleDsPhiPi(0.2);
   //      if(Nminus0.at(0).at(1).Integral()!=0) scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(1).Integral();
   //   }  
   double scaleDsPhiPi(0.93*0.676);
   for (unsigned int j=0; j<Extradist1d.size(); ++j) Extradist1d.at(j)->at(1).Scale(scaleDsPhiPi);
   }
   
   
   Selection::Finish();
}
