#include "MuonPionTree.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

using namespace std;

MuonPionTree::MuonPionTree(TString Name_, TString id_):
   Selection(Name_,id_),
   PhiMassLow(1.00),
   PhiMassHigh(1.04),
   sidebandDsMin(1.70),
   sidebandDsMax(1.80),
   peakDsMin(1.92),
   peakDsMax(2.01),
   dsMassMin(1.68),
   dsMassMax(2.02)
{
}

MuonPionTree::~MuonPionTree(){
   for(unsigned int j=0; j<Npassed.size(); j++){
      Logger(Logger::Info) << "Selection Summary before: "
         << Npassed.at(j).GetBinContent(1)  << " +/- " << Npassed.at(j).GetBinError(1)  << " after: "
         << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
   }


   Logger(Logger::Info) << "complete." << std::endl;
}

void  MuonPionTree::Configure(){

   TMVATree = new TTree("tree","tree"); 
   TMVATree->Branch("fake",&fake);
   TMVATree->Branch("isGlobal",&isGlobal);
   TMVATree->Branch("isTracker",&isTracker);
   TMVATree->Branch("isPF",&isPF);
   TMVATree->Branch("muonPt",&muonPt);
   TMVATree->Branch("muonEta",&muonEta);
   TMVATree->Branch("muonPhi",&muonPhi);
   TMVATree->Branch("muonInnerNC2",&muonInnerNC2);
   TMVATree->Branch("muonValidFraction",&muonValidFraction);
   TMVATree->Branch("muonInnerNValidHits",&muonInnerNValidHits);
   TMVATree->Branch("muonInnerTrackQuality",&muonInnerTrackQuality);
   TMVATree->Branch("muonNValidPixelHits",&muonNValidPixelHits);
   TMVATree->Branch("muonNValidTrackerHits",&muonNValidTrackerHits);
   TMVATree->Branch("muonNLostTrackerHits",&muonNLostTrackerHits);
   TMVATree->Branch("muonNLostTrackerHitsInner",&muonNLostTrackerHitsInner);
   TMVATree->Branch("muonNLostTrackerHitsOuter",&muonNLostTrackerHitsOuter);
   TMVATree->Branch("muonPixelLayers",&muonPixelLayers);
   TMVATree->Branch("muonTrackerLayers",&muonTrackerLayers);
   TMVATree->Branch("muonNMatchedStations",&muonNMatchedStations);
   TMVATree->Branch("muonNMatches",&muonNMatches);
   TMVATree->Branch("muonRPC",&muonRPC);
   TMVATree->Branch("muonPtErrPt",&muonPtErrPt);
   TMVATree->Branch("muonSegComp",&muonSegComp);
   TMVATree->Branch("muonCaloComp",&muonCaloComp);
   TMVATree->Branch("muonHadS9",&muonHadS9);
   TMVATree->Branch("muonHad",&muonHad);
   TMVATree->Branch("muonEM",&muonEM);
   TMVATree->Branch("muonEMS9",&muonEMS9);
   TMVATree->Branch("muonEMS25",&muonEMS25);
   TMVATree->Branch("muonKink",&muonKink);
   
   for(int i=0; i<NCuts;i++){
      cut.push_back(0);
      value.push_back(0);
      pass.push_back(false);
      if(i==Candidate)  cut.at(Candidate)=1;
   }
   TString hlabel;
   TString htitle;

   for(int i=0; i<NCuts; i++){
      title.push_back("");
      distindx.push_back(false);
      dist.push_back(std::vector<float>());
      TString c="_Cut_";c+=i;
      if(i==Candidate){
         title.at(i)="candidate";
         hlabel="candidate";
         Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Candidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
         Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Candidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
   }

   // Setup NPassed Histogams
   Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
   NDsPions=HConfig.GetTH2D(Name+"NDsPions","Number of Pions from Ds",5,-0.5,4.5,5,-0.5,4.5,"N D_{s}","N #pi");

   // Setup Extra Histograms
   Selection::ConfigureHistograms(); //do not remove
   HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}

void  MuonPionTree::Store_ExtraDist(){ 
   // store extradist
   Extradist2d.push_back(&NDsPions);
}

void  MuonPionTree::doEvent(){ 
   unsigned int t;
   int id(Ntp->GetMCID());
   bool status;

   double wobs=1;
   double w;

   // Initialize all cut values
   value.at(Candidate)=0;

   if(!Ntp->isData()){
      w = 1; // Weight MC according to truth number of vertices
   } 
   //  No weights to data
   else{w=1;}

   if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}

   pass.at(Candidate) = ( ( Ntp->NTwoMuonsTrack() >= cut.at(Candidate) ) || ( Ntp->NThreeMuons() >= cut.at(Candidate)) );

   status = AnalysisCuts(t,w,wobs);

   if(status){

      int nDsCounter = 0;
      int nPiCounter[5];

      for (int j=0; j<5; j++) nPiCounter[j] = 0;

      float mu_dR = 999.0;
      float tau_dR = 999.0;
      int mu_idx = -1;
      int final_idx = -1;
      bool pi_flag = 0;
      bool mu_flag = 0;

      if (id==30 && Ntp->NTwoMuonsTrack()>0){
         for (unsigned int ngen=0; ngen<Ntp->NMCSignalParticles(); ngen++){
            if (fabs(Ntp->MCSignalParticle_pdgid(ngen))==431){
               nDsCounter++;
               for (int nchild=0; nchild<Ntp->MCSignalParticle_Nchilds(ngen); nchild++){
                  //Find a pi that is a decay product of ds
                  if(fabs(Ntp->MCSignalParticle_childpdgid(ngen,nchild))==211){
                     nPiCounter[nDsCounter]++;
                     // Loop over all the muons and GEN match
                     for (unsigned int iMuon=0; iMuon<Ntp->NMuons(); iMuon++){
                        float muon_pt = Ntp->Muon_P4(iMuon).Pt();
                        float dpt = fabs(muon_pt-Ntp->MCSignalParticle_child_p4(ngen,nchild).Pt())/Ntp->MCSignalParticle_child_p4(ngen,nchild).Pt();
                        //bool muId = (Ntp->Muon_isPFMuon(iMuon) && Ntp->Muon_isTrackerMuon(iMuon) && !Ntp->Muon_isGlobalMuon(iMuon));
                        float dR = Ntp->Muon_P4(iMuon).DeltaR(Ntp->MCSignalParticle_child_p4(ngen,nchild));
                        if (dR<0.01 && dR<mu_dR && dpt<0.1 /*&& muId*/) { 
                           mu_dR = dR;
                           pi_flag = 1;
                           mu_idx = iMuon;
                        }
                     }
                  }
               }
            }
         }
      }

      else if (id==40 && Ntp->NThreeMuons()>0){
         for (unsigned int j=0; j<Ntp->NThreeMuons(); j++){
            float tmp_tau_dR = Ntp->TauGenMatch(j);
            if (tmp_tau_dR<tau_dR && tmp_tau_dR<0.01){
               final_idx = j;
               mu_flag = 1;
               tau_dR = tmp_tau_dR;
            }
         } 
      }

      NDsPions.at(t).Fill(nDsCounter,nPiCounter[nDsCounter]);

      // Fill Tree
      if(pi_flag){
         fake=1;
         isGlobal = Ntp->Muon_isGlobalMuon(mu_idx);
         isTracker = Ntp->Muon_isTrackerMuon(mu_idx);
         isPF = Ntp->Muon_isPFMuon(mu_idx);
         muonPt = Ntp->Muon_P4(mu_idx).Pt();
         muonEta = Ntp->Muon_P4(mu_idx).Eta();
         muonPhi = Ntp->Muon_P4(mu_idx).Phi();
         muonInnerNC2 = Ntp->Muon_innerTrack_normalizedChi2(mu_idx);
         muonValidFraction = Ntp->Muon_innerTrack_validFraction(mu_idx);
         muonInnerNValidHits = Ntp->Muon_innerTrack_numberofValidHits(mu_idx);
         muonInnerTrackQuality = Ntp->Muon_innerTrack_quality(mu_idx);
         muonNValidPixelHits = Ntp->Muon_numberofValidPixelHits(mu_idx);
         muonNValidTrackerHits = Ntp->Muon_innerTrack_numberOfValidTrackerHits(mu_idx);
         muonNLostTrackerHits = Ntp->Muon_innerTrack_numberOfLostTrackerHits(mu_idx);
         muonNLostTrackerHitsInner = Ntp->Muon_innerTrack_numberOfLostTrackerInnerHits(mu_idx);
         muonNLostTrackerHitsOuter = Ntp->Muon_innerTrack_numberOfLostTrackerOuterHits(mu_idx);
         muonPixelLayers = Ntp->Muon_innerTrack_pixelLayersWithMeasurement(mu_idx);
         muonTrackerLayers = Ntp->Muon_trackerLayersWithMeasurement(mu_idx);
         muonNMatchedStations = Ntp->Muon_numberOfMatchedStations(mu_idx);
         muonNMatches = Ntp->Muon_numberOfMatches(mu_idx);
         muonRPC = Ntp->Muon_isRPCMuon(mu_idx);
         muonPtErrPt = Ntp->Muon_ptErrOverPt(mu_idx);
         muonSegComp = Ntp->Muon_segmentCompatibility(mu_idx);
         muonCaloComp = Ntp->Muon_caloCompatibility(mu_idx);
         muonHadS9 = Ntp->Muon_calEnergy_hadS9(mu_idx);
         muonHad = Ntp->Muon_calEnergy_had(mu_idx);
         muonEM = Ntp->Muon_calEnergy_em(mu_idx);
         muonEMS9 = Ntp->Muon_calEnergy_emS9(mu_idx);
         muonEMS25 = Ntp->Muon_calEnergy_emS25(mu_idx);
         muonKink = Ntp->Muon_combinedQuality_trkKink(mu_idx);
         TMVATree->Fill();
      }

      if (mu_flag){
         for (unsigned int j=0; j<3; j++){
            unsigned int iMuon = Ntp->ThreeMuonIndices(final_idx).at(j);
            //if (!Ntp->Muon_isPFMuon(iMuon)) continue;
            //if (!Ntp->Muon_isTrackerMuon(iMuon) || Ntp->Muon_isGlobalMuon(iMuon)) continue;
            fake = 0;
            isGlobal = Ntp->Muon_isGlobalMuon(iMuon);
            isTracker = Ntp->Muon_isTrackerMuon(iMuon);
            isPF = Ntp->Muon_isPFMuon(iMuon);
            muonPt = Ntp->Muon_P4(iMuon).Pt();
            muonEta = Ntp->Muon_P4(iMuon).Eta();
            muonPhi = Ntp->Muon_P4(iMuon).Phi();
            muonInnerNC2 = Ntp->Muon_innerTrack_normalizedChi2(iMuon);
            muonValidFraction = Ntp->Muon_hitPattern_numberOfValidMuonHits(iMuon);
            muonInnerNValidHits = Ntp->Muon_innerTrack_numberofValidHits(iMuon);
            muonInnerTrackQuality = Ntp->Muon_innerTrack_quality(iMuon);
            muonNValidPixelHits = Ntp->Muon_numberofValidPixelHits(iMuon);
            muonNValidTrackerHits = Ntp->Muon_innerTrack_numberOfValidTrackerHits(iMuon);
            muonNLostTrackerHits = Ntp->Muon_innerTrack_numberOfLostTrackerHits(iMuon);
            muonNLostTrackerHitsInner = Ntp->Muon_innerTrack_numberOfLostTrackerInnerHits(iMuon);
            muonNLostTrackerHitsOuter = Ntp->Muon_innerTrack_numberOfLostTrackerOuterHits(iMuon);
            muonPixelLayers = Ntp->Muon_innerTrack_pixelLayersWithMeasurement(iMuon);
            muonTrackerLayers = Ntp->Muon_trackerLayersWithMeasurement(iMuon);
            muonNMatchedStations = Ntp->Muon_numberOfMatchedStations(iMuon);
            muonNMatches = Ntp->Muon_numberOfMatches(iMuon);
            muonRPC = Ntp->Muon_isRPCMuon(iMuon);
            muonPtErrPt = Ntp->Muon_ptErrOverPt(iMuon);
            muonSegComp = Ntp->Muon_segmentCompatibility(iMuon);
            muonCaloComp = Ntp->Muon_caloCompatibility(iMuon);
            muonHadS9 = Ntp->Muon_calEnergy_hadS9(iMuon);
            muonHad = Ntp->Muon_calEnergy_had(iMuon);
            muonEM = Ntp->Muon_calEnergy_em(iMuon);
            muonEMS9 = Ntp->Muon_calEnergy_emS9(iMuon);
            muonEMS25 = Ntp->Muon_calEnergy_emS25(iMuon);
            muonKink = Ntp->Muon_combinedQuality_trkKink(iMuon);
            TMVATree->Fill();
         }
      }
   }
}

void  MuonPionTree::Finish(){
   
   file= new TFile("MuonPionTree.root","recreate");
   TMVATree->SetDirectory(file);

   file->Write();
   file->Close();

  Selection::Finish();
}
