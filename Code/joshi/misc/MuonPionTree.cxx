#include "MuonPionTree.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

using namespace std;

MuonPionTree::MuonPionTree(TString Name_, TString id_):
   Selection(Name_,id_)
{
   file= new TFile("MuonPionTree.root","recreate");
   TMVA_Tree = new TTree("tree","tree"); 
   TMVA_Tree->SetDirectory(file);
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

   TMVA_Tree->Branch("fake",&fake);
   TMVA_Tree->Branch("Muon_pdgId",&Muon_pdgId);
   TMVA_Tree->Branch("Muon_motherPdgId", &Muon_motherPdgId);
   TMVA_Tree->Branch("Muon_Pt",&Muon_Pt);
   TMVA_Tree->Branch("Muon_Eta",&Muon_Eta);
   TMVA_Tree->Branch("Muon_Phi",&Muon_Phi);
   TMVA_Tree->Branch("Muon_vx",&Muon_vx);
   TMVA_Tree->Branch("Muon_vy",&Muon_vy);
   TMVA_Tree->Branch("Muon_vz",&Muon_vz);
   TMVA_Tree->Branch("Muon_IsGlobalMuon",&Muon_IsGlobalMuon);
   TMVA_Tree->Branch("Muon_IsStandAloneMuon",&Muon_IsStandAloneMuon);
   TMVA_Tree->Branch("Muon_IsTrackerMuon",&Muon_IsTrackerMuon);
   TMVA_Tree->Branch("Muon_IsCaloMuon",&Muon_IsCaloMuon);
   TMVA_Tree->Branch("Muon_IsIsolationValid",&Muon_IsIsolationValid);
   TMVA_Tree->Branch("Muon_IsQualityValid",&Muon_IsQualityValid);
   TMVA_Tree->Branch("Muon_IsTimeValid",&Muon_IsTimeValid);
   TMVA_Tree->Branch("Muon_IsPFMuon",&Muon_IsPFMuon);
   TMVA_Tree->Branch("Muon_IsRPCMuon",&Muon_IsRPCMuon);
   TMVA_Tree->Branch("Muon_emEt03",&Muon_emEt03);
   TMVA_Tree->Branch("Muon_emVetoEt03",&Muon_emVetoEt03);
   TMVA_Tree->Branch("Muon_hadEt03",&Muon_hadEt03);
   TMVA_Tree->Branch("Muon_hadVetoEt03",&Muon_hadVetoEt03);
   TMVA_Tree->Branch("Muon_nJets03",&Muon_nJets03);
   TMVA_Tree->Branch("Muon_nTracks03",&Muon_nTracks03);
   TMVA_Tree->Branch("Muon_StandardSelection",&Muon_StandardSelection);
   TMVA_Tree->Branch("Muon_trackerVetoPt03",&Muon_trackerVetoPt03);
   TMVA_Tree->Branch("Muon_sumChargedHadronPt03",&Muon_sumChargedHadronPt03);
   TMVA_Tree->Branch("Muon_sumChargedParticlePt03",&Muon_sumChargedParticlePt03);
   TMVA_Tree->Branch("Muon_sumNeutralHadronEt03",&Muon_sumNeutralHadronEt03);
   TMVA_Tree->Branch("Muon_sumNeutralHadronEtHighThreshold03",&Muon_sumNeutralHadronEtHighThreshold03);
   TMVA_Tree->Branch("Muon_sumPhotonEt03",&Muon_sumPhotonEt03);
   TMVA_Tree->Branch("Muon_sumPhotonEtHighThreshold03",&Muon_sumPhotonEtHighThreshold03);
   TMVA_Tree->Branch("Muon_sumPUPt03",&Muon_sumPUPt03);
   TMVA_Tree->Branch("Muon_numberOfChambers",&Muon_numberOfChambers);
   TMVA_Tree->Branch("Muon_Track_idx",&Muon_Track_idx);
   TMVA_Tree->Branch("Muon_combinedQuality_updatedSta",&Muon_combinedQuality_updatedSta);
   TMVA_Tree->Branch("Muon_combinedQuality_trkKink",&Muon_combinedQuality_trkKink);
   TMVA_Tree->Branch("Muon_combinedQuality_glbKink",&Muon_combinedQuality_glbKink);
   TMVA_Tree->Branch("Muon_combinedQuality_trkRelChi2",&Muon_combinedQuality_trkRelChi2);
   TMVA_Tree->Branch("Muon_combinedQuality_staRelChi2",&Muon_combinedQuality_staRelChi2);
   TMVA_Tree->Branch("Muon_combinedQuality_chi2LocalPosition",&Muon_combinedQuality_chi2LocalPosition);
   TMVA_Tree->Branch("Muon_combinedQuality_chi2LocalMomentum",&Muon_combinedQuality_chi2LocalMomentum);
   TMVA_Tree->Branch("Muon_combinedQuality_localDistance",&Muon_combinedQuality_localDistance);
   TMVA_Tree->Branch("Muon_combinedQuality_globalDeltaEtaPhi",&Muon_combinedQuality_globalDeltaEtaPhi);
   TMVA_Tree->Branch("Muon_combinedQuality_tightMatch",&Muon_combinedQuality_tightMatch);
   TMVA_Tree->Branch("Muon_combinedQuality_glbTrackProbability",&Muon_combinedQuality_glbTrackProbability);
   TMVA_Tree->Branch("Muon_prod_inner_outer_charge",&Muon_prod_inner_outer_charge);
   TMVA_Tree->Branch("Muon_innerTrack_quality",&Muon_innerTrack_quality);
   TMVA_Tree->Branch("Muon_ptErrOverPt",&Muon_ptErrOverPt);
   TMVA_Tree->Branch("Muon_calEnergy_hadS9",&Muon_calEnergy_hadS9);
   TMVA_Tree->Branch("Muon_calEnergy_had",&Muon_calEnergy_had);
   TMVA_Tree->Branch("Muon_calEnergy_emS25",&Muon_calEnergy_emS25);
   TMVA_Tree->Branch("Muon_calEnergy_emS9",&Muon_calEnergy_emS9);
   TMVA_Tree->Branch("Muon_calEnergy_em",&Muon_calEnergy_em);
   TMVA_Tree->Branch("Muon_charge",&Muon_charge);
   TMVA_Tree->Branch("Muon_trackCharge",&Muon_trackCharge);
   TMVA_Tree->Branch("Muon_hitPattern_pixelLayerwithMeas",&Muon_hitPattern_pixelLayerwithMeas);
   TMVA_Tree->Branch("Muon_numberOfMatchedStations",&Muon_numberOfMatchedStations);
   TMVA_Tree->Branch("Muon_normChi2",&Muon_normChi2);
   TMVA_Tree->Branch("Muon_hitPattern_numberOfValidMuonHits",&Muon_hitPattern_numberOfValidMuonHits);
   TMVA_Tree->Branch("Muon_innerTrack_numberofValidHits",&Muon_innerTrack_numberofValidHits);
   TMVA_Tree->Branch("Muon_numberofValidPixelHits",&Muon_numberofValidPixelHits);
   TMVA_Tree->Branch("Muon_numberOfMatches",&Muon_numberOfMatches);
   TMVA_Tree->Branch("Muon_trackerLayersWithMeasurement",&Muon_trackerLayersWithMeasurement);
   TMVA_Tree->Branch("Muon_segmentCompatibility",&Muon_segmentCompatibility);
   TMVA_Tree->Branch("Muon_caloCompatibility",&Muon_caloCompatibility);
   TMVA_Tree->Branch("Muon_innerTrack_validFraction",&Muon_innerTrack_validFraction);
   TMVA_Tree->Branch("Muon_innerTrack_pixelLayersWithMeasurement",&Muon_innerTrack_pixelLayersWithMeasurement);
   TMVA_Tree->Branch("Muon_innerTrack_numberOfValidTrackerHits",&Muon_innerTrack_numberOfValidTrackerHits);
   TMVA_Tree->Branch("Muon_innerTrack_numberOfLostTrackerHits",&Muon_innerTrack_numberOfLostTrackerHits);
   TMVA_Tree->Branch("Muon_innerTrack_numberOfLostTrackerInnerHits",&Muon_innerTrack_numberOfLostTrackerInnerHits);
   TMVA_Tree->Branch("Muon_innerTrack_numberOfLostTrackerOuterHits",&Muon_innerTrack_numberOfLostTrackerOuterHits);
   TMVA_Tree->Branch("Muon_innerTrack_normalizedChi2",&Muon_innerTrack_normalizedChi2);
   TMVA_Tree->Branch("Muon_vmuonhitcomb_reco",&Muon_vmuonhitcomb_reco);
   TMVA_Tree->Branch("Muon_rpchits_reco",&Muon_rpchits_reco);
   TMVA_Tree->Branch("Muon_outerTrack_normalizedChi2",&Muon_outerTrack_normalizedChi2);
   TMVA_Tree->Branch("Muon_outerTrack_muonStationsWithValidHits",&Muon_outerTrack_muonStationsWithValidHits);
   TMVA_Tree->Branch("Muon_isGoodMuon_TM2DCompatibility",&Muon_isGoodMuon_TM2DCompatibility);
   TMVA_Tree->Branch("Muon_isGoodMuon_TrackerMuonArbitrated",&Muon_isGoodMuon_TrackerMuonArbitrated);
   TMVA_Tree->Branch("Muon_isGoodMuon_TMOneStationTight",&Muon_isGoodMuon_TMOneStationTight);
   TMVA_Tree->Branch("Muon_isGoodMuon_TMOneStationAngTight",&Muon_isGoodMuon_TMOneStationAngTight);
   TMVA_Tree->Branch("Muon_isGoodMuon_TMLastStationTight",&Muon_isGoodMuon_TMLastStationTight);
   TMVA_Tree->Branch("Muon_isGoodMuon_TMLastStationAngTight",&Muon_isGoodMuon_TMLastStationAngTight);
   TMVA_Tree->Branch("Muon_isGoodMuon_TMLastStationOptimizedLowPtTight",&Muon_isGoodMuon_TMLastStationOptimizedLowPtTight);
   TMVA_Tree->Branch("Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight",&Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight);

   TMVA_Tree->Branch("Muon_station1_status", &Muon_station_vars[0][0]);
   TMVA_Tree->Branch("Muon_station1_TrackX", &Muon_station_vars[0][1]);
   TMVA_Tree->Branch("Muon_station1_TrackY", &Muon_station_vars[0][2]);
   TMVA_Tree->Branch("Muon_station1_pullX", &Muon_station_vars[0][3]);
   TMVA_Tree->Branch("Muon_station1_pullY", &Muon_station_vars[0][4]);
   TMVA_Tree->Branch("Muon_station1_pullDxDz", &Muon_station_vars[0][5]);
   TMVA_Tree->Branch("Muon_station1_pullDyDz", &Muon_station_vars[0][6]);

   TMVA_Tree->Branch("Muon_station2_status", &Muon_station_vars[1][0]);
   TMVA_Tree->Branch("Muon_station2_TrackX", &Muon_station_vars[1][1]);
   TMVA_Tree->Branch("Muon_station2_TrackY", &Muon_station_vars[1][2]);
   TMVA_Tree->Branch("Muon_station2_pullX", &Muon_station_vars[1][3]);
   TMVA_Tree->Branch("Muon_station2_pullY", &Muon_station_vars[1][4]);
   TMVA_Tree->Branch("Muon_station2_pullDxDz", &Muon_station_vars[1][5]);
   TMVA_Tree->Branch("Muon_station2_pullDyDz", &Muon_station_vars[1][6]);

   TMVA_Tree->Branch("Muon_station3_status", &Muon_station_vars[2][0]);
   TMVA_Tree->Branch("Muon_station3_TrackX", &Muon_station_vars[2][1]);
   TMVA_Tree->Branch("Muon_station3_TrackY", &Muon_station_vars[2][2]);
   TMVA_Tree->Branch("Muon_station3_pullX", &Muon_station_vars[2][3]);
   TMVA_Tree->Branch("Muon_station3_pullY", &Muon_station_vars[2][4]);
   TMVA_Tree->Branch("Muon_station3_pullDxDz", &Muon_station_vars[2][5]);
   TMVA_Tree->Branch("Muon_station3_pullDyDz", &Muon_station_vars[2][6]);

   TMVA_Tree->Branch("Muon_station4_status", &Muon_station_vars[3][0]);
   TMVA_Tree->Branch("Muon_station4_TrackX", &Muon_station_vars[3][1]);
   TMVA_Tree->Branch("Muon_station4_TrackY", &Muon_station_vars[3][2]);
   TMVA_Tree->Branch("Muon_station4_pullX", &Muon_station_vars[3][3]);
   TMVA_Tree->Branch("Muon_station4_pullY", &Muon_station_vars[3][4]);
   TMVA_Tree->Branch("Muon_station4_pullDxDz", &Muon_station_vars[3][5]);
   TMVA_Tree->Branch("Muon_station4_pullDyDz", &Muon_station_vars[3][6]);

   TMVA_Tree->Branch("Muon_station5_status", &Muon_station_vars[4][0]);
   TMVA_Tree->Branch("Muon_station5_TrackX", &Muon_station_vars[4][1]);
   TMVA_Tree->Branch("Muon_station5_TrackY", &Muon_station_vars[4][2]);
   TMVA_Tree->Branch("Muon_station5_pullX", &Muon_station_vars[4][3]);
   TMVA_Tree->Branch("Muon_station5_pullY", &Muon_station_vars[4][4]);
   TMVA_Tree->Branch("Muon_station5_pullDxDz", &Muon_station_vars[4][5]);
   TMVA_Tree->Branch("Muon_station5_pullDyDz", &Muon_station_vars[4][6]);

   TMVA_Tree->Branch("Muon_station6_status", &Muon_station_vars[5][0]);
   TMVA_Tree->Branch("Muon_station6_TrackX", &Muon_station_vars[5][1]);
   TMVA_Tree->Branch("Muon_station6_TrackY", &Muon_station_vars[5][2]);
   TMVA_Tree->Branch("Muon_station6_pullX", &Muon_station_vars[5][3]);
   TMVA_Tree->Branch("Muon_station6_pullY", &Muon_station_vars[5][4]);
   TMVA_Tree->Branch("Muon_station6_pullDxDz", &Muon_station_vars[5][5]);
   TMVA_Tree->Branch("Muon_station6_pullDyDz", &Muon_station_vars[5][6]);

   TMVA_Tree->Branch("Muon_station7_status", &Muon_station_vars[6][0]);
   TMVA_Tree->Branch("Muon_station7_TrackX", &Muon_station_vars[6][1]);
   TMVA_Tree->Branch("Muon_station7_TrackY", &Muon_station_vars[6][2]);
   TMVA_Tree->Branch("Muon_station7_pullX", &Muon_station_vars[6][3]);
   TMVA_Tree->Branch("Muon_station7_pullY", &Muon_station_vars[6][4]);
   TMVA_Tree->Branch("Muon_station7_pullDxDz", &Muon_station_vars[6][5]);
   TMVA_Tree->Branch("Muon_station7_pullDyDz", &Muon_station_vars[6][6]);

   TMVA_Tree->Branch("Muon_station8_status", &Muon_station_vars[7][0]);
   TMVA_Tree->Branch("Muon_station8_TrackX", &Muon_station_vars[7][1]);
   TMVA_Tree->Branch("Muon_station8_TrackY", &Muon_station_vars[7][2]);
   TMVA_Tree->Branch("Muon_station8_pullX", &Muon_station_vars[7][3]);
   TMVA_Tree->Branch("Muon_station8_pullY", &Muon_station_vars[7][4]);
   TMVA_Tree->Branch("Muon_station8_pullDxDz", &Muon_station_vars[7][5]);
   TMVA_Tree->Branch("Muon_station8_pullDyDz", &Muon_station_vars[7][6]);


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

   // Setup Extra Histograms
   Selection::ConfigureHistograms(); //do not remove
   HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}

void  MuonPionTree::Store_ExtraDist(){ 
   // store extradist
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

   //pass.at(Candidate) = ( ( Ntp->NTwoMuonsTrack() >= cut.at(Candidate) ) || ( Ntp->NThreeMuons() >= cut.at(Candidate)) );
   pass.at(Candidate) = true;

   status = AnalysisCuts(t,w,wobs);

   if(status){

      for (unsigned int iMuon=0; iMuon<Ntp->NMuons(); iMuon++){
         //Match RECO muons to GEN MC particles
         std::pair<int, int> genId_index_pair = Ntp->GENMatchedPdgId(Ntp->Muon_P4(iMuon));
         Muon_pdgId = genId_index_pair.first;
         int gen_index = genId_index_pair.second;
         if (abs(Muon_pdgId)==211 || abs(Muon_pdgId)==321 || abs(Muon_pdgId)==13){
            if (abs(Muon_pdgId)==211 || abs(Muon_pdgId)==321) fake = true;
            else fake = false;

            TLorentzVector gen_p4 = Ntp->MCParticle_p4(gen_index);
            if (Ntp->MCParticle_hasMother(gen_index)) 
               Muon_motherPdgId = Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(gen_index));
            else Muon_motherPdgId = 0;

            Muon_Pt = Ntp->Muon_P4(iMuon).Pt();
            Muon_Eta = Ntp->Muon_P4(iMuon).Eta();
            Muon_Phi = Ntp->Muon_P4(iMuon).Phi();
            Muon_vx = Ntp->Muon_Poca(iMuon).X();
            Muon_vy = Ntp->Muon_Poca(iMuon).Y();
            Muon_vz = Ntp->Muon_Poca(iMuon).Z();
            Muon_IsGlobalMuon = Ntp->Muon_isGlobalMuon(iMuon);
            Muon_IsStandAloneMuon = Ntp->Muon_isStandAloneMuon(iMuon);
            Muon_IsTrackerMuon = Ntp->Muon_isTrackerMuon(iMuon);
            Muon_IsCaloMuon = Ntp->Muon_isCaloMuon(iMuon);
            Muon_IsIsolationValid = Ntp->Muon_isIsolationValid(iMuon);
            Muon_IsQualityValid = Ntp->Muon_isQualityValid(iMuon);
            Muon_IsTimeValid = Ntp->Muon_isTimeValid(iMuon);
            Muon_IsPFMuon = Ntp->Muon_isPFMuon(iMuon);
            Muon_IsRPCMuon = Ntp->Muon_isRPCMuon(iMuon);
            Muon_emEt03 = Ntp->Muon_emEt03(iMuon);
            Muon_emVetoEt03 = Ntp->Muon_emVetoEt03(iMuon);
            Muon_hadEt03 = Ntp->Muon_hadEt03(iMuon);
            Muon_hadVetoEt03 = Ntp->Muon_hadVetoEt03(iMuon);
            Muon_nJets03 = Ntp->Muon_nJets03(iMuon);
            Muon_nTracks03 = Ntp->Muon_nTracks03(iMuon);
            Muon_StandardSelection = Ntp->Muon_StandardSelection(iMuon);
            Muon_sumPt03 = Ntp->Muon_sumPt03(iMuon);
            Muon_trackerVetoPt03 = Ntp->Muon_trackerVetoPt03(iMuon);
            Muon_sumChargedHadronPt03 = Ntp->Muon_sumChargedHadronPt03(iMuon);
            Muon_sumChargedParticlePt03 = Ntp->Muon_sumChargedParticlePt03(iMuon);
            Muon_sumNeutralHadronEt03 = Ntp->Muon_sumNeutralHadronEt03(iMuon);
            Muon_sumNeutralHadronEtHighThreshold03 = Ntp->Muon_sumNeutralHadronEtHighThreshold03(iMuon);
            Muon_sumPhotonEt03 = Ntp->Muon_sumPhotonEt03(iMuon);
            Muon_sumPhotonEtHighThreshold03 = Ntp->Muon_sumPhotonEtHighThreshold03(iMuon);
            Muon_sumPUPt03 = Ntp->Muon_sumPUPt03(iMuon);
            Muon_numberOfChambers = Ntp->Muon_numberOfChambers(iMuon);
            Muon_Track_idx = Ntp->Muon_Track_idx(iMuon);
            Muon_combinedQuality_updatedSta = Ntp->Muon_combinedQuality_updatedSta(iMuon);
            Muon_combinedQuality_trkKink = Ntp->Muon_combinedQuality_trkKink(iMuon);
            Muon_combinedQuality_glbKink = Ntp->Muon_combinedQuality_glbKink(iMuon);
            Muon_combinedQuality_trkRelChi2 = Ntp->Muon_combinedQuality_trkRelChi2(iMuon);
            Muon_combinedQuality_staRelChi2 = Ntp->Muon_combinedQuality_staRelChi2(iMuon);
            Muon_combinedQuality_chi2LocalPosition = Ntp->Muon_combinedQuality_chi2LocalPosition(iMuon);
            Muon_combinedQuality_chi2LocalMomentum = Ntp->Muon_combinedQuality_chi2LocalMomentum(iMuon);
            Muon_combinedQuality_localDistance = Ntp->Muon_combinedQuality_localDistance(iMuon);
            Muon_combinedQuality_globalDeltaEtaPhi = Ntp->Muon_combinedQuality_globalDeltaEtaPhi(iMuon);
            Muon_combinedQuality_tightMatch = Ntp->Muon_combinedQuality_tightMatch(iMuon);
            Muon_combinedQuality_glbTrackProbability = Ntp->Muon_combinedQuality_glbTrackProbability(iMuon);
            Muon_prod_inner_outer_charge = Ntp->Muon_prod_inner_outer_charge(iMuon);
            Muon_innerTrack_quality = Ntp->Muon_innerTrack_quality(iMuon);
            Muon_ptErrOverPt = Ntp->Muon_ptErrOverPt(iMuon);
            Muon_calEnergy_hadS9 = Ntp->Muon_calEnergy_hadS9(iMuon);
            Muon_calEnergy_had = Ntp->Muon_calEnergy_had(iMuon);
            Muon_calEnergy_emS25 = Ntp->Muon_calEnergy_emS25(iMuon);
            Muon_calEnergy_emS9 = Ntp->Muon_calEnergy_emS9(iMuon);
            Muon_calEnergy_em = Ntp->Muon_calEnergy_em(iMuon);
            Muon_charge = Ntp->Muon_charge(iMuon);
            Muon_trackCharge = Ntp->Muon_trackCharge(iMuon);
            Muon_hitPattern_pixelLayerwithMeas = Ntp->Muon_hitPattern_pixelLayerwithMeas(iMuon);
            Muon_numberOfMatchedStations = Ntp->Muon_numberOfMatchedStations(iMuon);
            Muon_normChi2 = Ntp->Muon_normChi2(iMuon);
            Muon_hitPattern_numberOfValidMuonHits = Ntp->Muon_hitPattern_numberOfValidMuonHits(iMuon);
            Muon_innerTrack_numberofValidHits = Ntp->Muon_innerTrack_numberofValidHits(iMuon);
            Muon_numberofValidPixelHits = Ntp->Muon_numberofValidPixelHits(iMuon);
            Muon_numberOfMatches = Ntp->Muon_numberOfMatches(iMuon);
            Muon_trackerLayersWithMeasurement = Ntp->Muon_trackerLayersWithMeasurement(iMuon);
            Muon_segmentCompatibility = Ntp->Muon_segmentCompatibility(iMuon);
            Muon_caloCompatibility = Ntp->Muon_caloCompatibility(iMuon);
            Muon_innerTrack_validFraction = Ntp->Muon_innerTrack_validFraction(iMuon);
            Muon_innerTrack_pixelLayersWithMeasurement = Ntp->Muon_innerTrack_pixelLayersWithMeasurement(iMuon);
            Muon_innerTrack_numberOfValidTrackerHits = Ntp->Muon_innerTrack_numberOfValidTrackerHits(iMuon);
            Muon_innerTrack_numberOfLostTrackerHits = Ntp->Muon_innerTrack_numberOfLostTrackerHits(iMuon);
            Muon_innerTrack_numberOfLostTrackerInnerHits = Ntp->Muon_innerTrack_numberOfLostTrackerInnerHits(iMuon);
            Muon_innerTrack_numberOfLostTrackerOuterHits = Ntp->Muon_innerTrack_numberOfLostTrackerOuterHits(iMuon);
            Muon_innerTrack_normalizedChi2 = Ntp->Muon_innerTrack_normalizedChi2(iMuon);
            Muon_vmuonhitcomb_reco = Ntp->Muon_vmuonhitcomb_reco(iMuon);
            Muon_rpchits_reco = Ntp->Muon_rpchits_reco(iMuon);
            Muon_outerTrack_normalizedChi2 = Ntp->Muon_outerTrack_normalizedChi2(iMuon);
            Muon_outerTrack_muonStationsWithValidHits = Ntp->Muon_outerTrack_muonStationsWithValidHits(iMuon);
            Muon_isGoodMuon_TM2DCompatibility = Ntp->Muon_isGoodMuon_TM2DCompatibility(iMuon);
            Muon_isGoodMuon_TrackerMuonArbitrated = Ntp->Muon_isGoodMuon_TrackerMuonArbitrated(iMuon);
            Muon_isGoodMuon_TMOneStationTight = Ntp->Muon_isGoodMuon_TMOneStationTight(iMuon);
            Muon_isGoodMuon_TMOneStationAngTight = Ntp->Muon_isGoodMuon_TMOneStationAngTight(iMuon);
            Muon_isGoodMuon_TMLastStationTight = Ntp->Muon_isGoodMuon_TMLastStationTight(iMuon);
            Muon_isGoodMuon_TMLastStationAngTight = Ntp->Muon_isGoodMuon_TMLastStationAngTight(iMuon);
            Muon_isGoodMuon_TMLastStationOptimizedLowPtTight = Ntp->Muon_isGoodMuon_TMLastStationOptimizedLowPtTight(iMuon);
            Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight = Ntp->Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight(iMuon);

            for (int istation=0; istation<8; istation++){
               int station_status = 0;

               if (Ntp->Muon_TrackX(iMuon, istation)<999999.0 && Ntp->Muon_TrackY(iMuon, istation)<999999.0) station_status++;
               if (Ntp->Muon_dX(iMuon, istation)>=999999.0 && station_status) station_status = -1;
               Muon_station_vars[istation][0] = station_status;
               Muon_station_vars[istation][1] = Ntp->Muon_TrackX(iMuon, istation);
               Muon_station_vars[istation][2] = Ntp->Muon_TrackY(iMuon, istation);
               Muon_station_vars[istation][3] = Ntp->Muon_pullX(iMuon, istation);
               Muon_station_vars[istation][4] = Ntp->Muon_pullY(iMuon, istation);
               Muon_station_vars[istation][5] = Ntp->Muon_pullDxDz(iMuon, istation);
               Muon_station_vars[istation][6] = Ntp->Muon_pullDyDz(iMuon, istation);

            }
            TMVA_Tree->Fill();
         }

      }

   }
}

void  MuonPionTree::Finish(){

   file->Write();
   file->Close();

   Selection::Finish();
}
