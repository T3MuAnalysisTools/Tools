#include "MakeMVATree.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


using namespace std;

MakeMVATree::MakeMVATree(TString Name_, TString id_):
    Selection(Name_,id_),
    tauMinMass_(1.731),
    tauMaxMass_(1.823),
    tauMinSideBand_(1.61),
    tauMaxSideBand_(2.05),
    tauMassResCutLow(0.007),
    tauMassResCutHigh(0.01)
{
    // This is a class constructor;

  //  TString basedir = "";
  //  basedir = (TString)std::getenv("workdir")+"/Code/CommonFiles/PileUp/Collisions2018";
  //  PUWeightFile = new TFile(basedir+"/PUWeights_Run2018.root");
  //  puWeights = (TH1D*)PUWeightFile->Get("h1_weights");


}


MakeMVATree::~MakeMVATree(){
    for(unsigned int j=0; j<Npassed.size(); j++){
      Logger(Logger::Info) << "Selection Summary before: "
        << Npassed.at(j).GetBinContent(1)  << " +/- " << Npassed.at(j).GetBinError(1)  << " after: "
        << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
    }
    Logger(Logger::Info) << "complete." << std::endl;
}

void  MakeMVATree::Configure(){
  // Set tree branches


  readerMuIDBarrel= new TMVA::Reader( "!Color:!Silent" );
  TString basedir = "";



  basedir = (TString)std::getenv("workdir")+"Code/CommonFiles/";
  readerMuIDBarrel->AddVariable("mu_combinedQuality_chi2LocalMomentum>250?250:mu_combinedQuality_chi2LocalMomentum" ,&mu_combinedQuality_chi2LocalMomentum );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_chi2LocalPosition>50?50:mu_combinedQuality_chi2LocalPosition" ,&mu_combinedQuality_chi2LocalPosition );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_staRelChi2>50?50:mu_combinedQuality_staRelChi2" ,&mu_combinedQuality_staRelChi2 );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_trkRelChi2>50?50:mu_combinedQuality_trkRelChi2" ,&mu_combinedQuality_trkRelChi2 );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_globalDeltaEtaPhi" ,&mu_combinedQuality_globalDeltaEtaPhi );
  readerMuIDBarrel->AddVariable("log(mu_combinedQuality_trkKink)" ,&mu_combinedQuality_trkKink );
  readerMuIDBarrel->AddVariable("log(mu_combinedQuality_glbKink)" ,&mu_combinedQuality_glbKink );
  readerMuIDBarrel->AddVariable("mu_combinedQuality_glbTrackProbability>150?150:mu_combinedQuality_glbTrackProbability" ,&mu_combinedQuality_glbTrackProbability );
  //  readerMuIDBarrel->AddVariable("mu_Numberofvalidtrackerhits" ,&mu_Numberofvalidtrackerhits );


  readerMuIDBarrel->AddVariable("mu_Numberofvalidpixelhits" ,&mu_Numberofvalidpixelhits );
  readerMuIDBarrel->AddVariable("mu_trackerLayersWithMeasurement" ,&mu_trackerLayersWithMeasurement );
  readerMuIDBarrel->AddVariable("mu_validMuonHitComb" ,&mu_validMuonHitComb );
  readerMuIDBarrel->AddVariable("mu_numberOfMatchedStations" ,&mu_numberOfMatchedStations );
  readerMuIDBarrel->AddVariable("mu_segmentCompatibility" ,&mu_segmentCompatibility );
  readerMuIDBarrel->AddVariable("mu_timeAtIpInOutErr" ,&mu_timeAtIpInOutErr );
  readerMuIDBarrel->AddVariable("mu_GLnormChi2>50?50:mu_GLnormChi2" ,&mu_GLnormChi2 );
  readerMuIDBarrel->AddVariable("mu_innerTrack_normalizedChi2>50?50:mu_innerTrack_normalizedChi2" ,&mu_innerTrack_normalizedChi2 );
  readerMuIDBarrel->AddVariable("mu_outerTrack_normalizedChi2>80?80:mu_outerTrack_normalizedChi2" ,&mu_outerTrack_normalizedChi2 );
  readerMuIDBarrel->AddVariable("mu_innerTrack_validFraction" ,&mu_innerTrack_validFraction );
  readerMuIDBarrel->AddSpectator("mu_eta" ,&mu_eta);
  readerMuIDBarrel->AddSpectator("mu_pt" ,&mu_pt);
  readerMuIDBarrel->AddSpectator("mu_phi" ,&mu_phi);
  readerMuIDBarrel->AddSpectator("mu_SoftMVA" ,&mu_SoftMVA);
  readerMuIDBarrel->BookMVA( "BDT", basedir+"MuonMVA_2018_mar2021_barrel/weights/TMVA_new_BDT.weights.xml" ); // weights xml file after training, place it to CommonFiles



  readerMuIDEndcap= new TMVA::Reader( "!Color:!Silent" );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_chi2LocalMomentum>250?250:mu_combinedQuality_chi2LocalMomentum" ,&mu_combinedQuality_chi2LocalMomentum );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_chi2LocalPosition>50?50:mu_combinedQuality_chi2LocalPosition" ,&mu_combinedQuality_chi2LocalPosition );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_staRelChi2>50?50:mu_combinedQuality_staRelChi2" ,&mu_combinedQuality_staRelChi2 );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_trkRelChi2>50?50:mu_combinedQuality_trkRelChi2" ,&mu_combinedQuality_trkRelChi2 );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_globalDeltaEtaPhi" ,&mu_combinedQuality_globalDeltaEtaPhi );
  readerMuIDEndcap->AddVariable("log(mu_combinedQuality_trkKink)" ,&mu_combinedQuality_trkKink );
  readerMuIDEndcap->AddVariable("log(mu_combinedQuality_glbKink)" ,&mu_combinedQuality_glbKink );
  readerMuIDEndcap->AddVariable("mu_combinedQuality_glbTrackProbability>150?150:mu_combinedQuality_glbTrackProbability" ,&mu_combinedQuality_glbTrackProbability );
  //  readerMuIDEndcap->AddVariable("mu_Numberofvalidtrackerhits" ,&mu_Numberofvalidtrackerhits );
  readerMuIDEndcap->AddVariable("mu_Numberofvalidpixelhits" ,&mu_Numberofvalidpixelhits );
  readerMuIDEndcap->AddVariable("mu_trackerLayersWithMeasurement" ,&mu_trackerLayersWithMeasurement );
  readerMuIDEndcap->AddVariable("mu_validMuonHitComb" ,&mu_validMuonHitComb );
  readerMuIDEndcap->AddVariable("mu_numberOfMatchedStations" ,&mu_numberOfMatchedStations );
  readerMuIDEndcap->AddVariable("mu_segmentCompatibility" ,&mu_segmentCompatibility );
  readerMuIDEndcap->AddVariable("mu_timeAtIpInOutErr" ,&mu_timeAtIpInOutErr );
  readerMuIDEndcap->AddVariable("mu_GLnormChi2>50?50:mu_GLnormChi2" ,&mu_GLnormChi2 );
  readerMuIDEndcap->AddVariable("mu_innerTrack_normalizedChi2>50?50:mu_innerTrack_normalizedChi2" ,&mu_innerTrack_normalizedChi2 );
  readerMuIDEndcap->AddVariable("mu_outerTrack_normalizedChi2>80?80:mu_outerTrack_normalizedChi2" ,&mu_outerTrack_normalizedChi2 );
  readerMuIDEndcap->AddVariable("mu_innerTrack_validFraction" ,&mu_innerTrack_validFraction );
  readerMuIDEndcap->AddSpectator("mu_eta" ,&mu_eta);
  readerMuIDEndcap->AddSpectator("mu_pt" ,&mu_pt);
  readerMuIDEndcap->AddSpectator("mu_phi" ,&mu_phi);
  readerMuIDEndcap->AddSpectator("mu_SoftMVA" ,&mu_SoftMVA);
  readerMuIDEndcap->BookMVA( "BDT", basedir+"MuonMVA_2018_mar2021_endcap/weights/TMVA_new_BDT.weights.xml" ); // weights xml file after training, place it to CommonFiles




  readerBvsD= new TMVA::Reader( "!Color:!Silent" );
  readerBvsD->AddVariable("var_flightLenSig",&var_flightLenSig);
  readerBvsD->AddVariable("var_svpvTauAngle",&var_svpvTauAngle);
  readerBvsD->AddVariable("var_MindcaTrackSV",&var_MindcaTrackSV);
  readerBvsD->AddVariable("var_NtracksClose",&var_NtracksClose);
  readerBvsD->AddVariable("var_nsv",&var_nsv);
  readerBvsD->AddVariable("var_MaxD0SigBS",&var_MaxD0SigBS);
  readerBvsD->AddVariable("var_MinD0SigBS",&var_MinD0SigBS);
  readerBvsD->AddVariable("var_Iso08",&var_Iso08);
  readerBvsD->AddVariable("var_dcaTrackPV",&var_dcaTrackPV);
  readerBvsD->AddVariable("var_MinMuonImpactAngle",&var_MinMuonImpactAngle);
  readerBvsD->AddVariable("var_flightLenDist",&var_flightLenDist);
  readerBvsD->BookMVA( "BDTG", "/afs/cern.ch/work/c/cherepan/Analysis/workdirHFMVAJul_10_2023/Code/CommonUtils/IterativeTrain/output_0_MCTrainA/weights/TMVAClassification_BDTG.weights.xml" );








  //  Error in <TTree::Fill>: Failed filling branch:tree.var_MinMuon_chi2LocalPosition, nbytes=-1, entry=7975

  TString treeprefix;
  if(Ntp->GetInputNtuplePath().Contains("Z2TauTau")) treeprefix="DoubleMuonLowMass";

  if(Ntp->GetInputNtuplePath().Contains("DsToTau")) treeprefix="Ds";
  if(Ntp->GetInputNtuplePath().Contains("DpToTau")) treeprefix="Dp";
  if(Ntp->GetInputNtuplePath().Contains("BdToTau")) treeprefix="B0";
  if(Ntp->GetInputNtuplePath().Contains("BuToTau")) treeprefix="Bp";


 
  TMVA_Tree= new TTree(treeprefix+"_mva_tree","Mini tree for MVA");
  // this is now working version
  TMVA_Tree->Branch("MC",&MC);
  TMVA_Tree->Branch("category",&category);
  TMVA_Tree->Branch("var_id", &var_id);
  TMVA_Tree->Branch("var_EventNumber", &var_EventNumber);


  TMVA_Tree->Branch("var_vertexKFChi2",&var_vertexKFChi2);
  TMVA_Tree->Branch("var_svpvTauAngle",&var_svpvTauAngle);
  TMVA_Tree->Branch("var_flightLenSig",&var_flightLenSig);
  TMVA_Tree->Branch("var_flightLenDist",&var_flightLenDist);
  TMVA_Tree->Branch("var_sumMuTrkKinkChi2",&var_sumMuTrkKinkChi2);
  TMVA_Tree->Branch("var_segCompMuMin",&var_segCompMuMin);
  TMVA_Tree->Branch("var_segCompMuMax",&var_segCompMuMax);
  TMVA_Tree->Branch("var_segCompMu1",&var_segCompMu1);
  TMVA_Tree->Branch("var_segCompMu2",&var_segCompMu2);
  TMVA_Tree->Branch("var_segCompMu3",&var_segCompMu3);
  TMVA_Tree->Branch("var_caloCompMin",&var_caloCompMin);
  TMVA_Tree->Branch("var_caloCompMax",&var_caloCompMax);
  TMVA_Tree->Branch("var_caloCompMu1",&var_caloCompMu1);
  TMVA_Tree->Branch("var_caloCompMu2",&var_caloCompMu2);
  TMVA_Tree->Branch("var_caloCompMu3",&var_caloCompMu3);
  TMVA_Tree->Branch("var_MinMIPLikelihood",&var_MinMIPLikelihood);
  TMVA_Tree->Branch("var_tauMass",&var_tauMass);
  TMVA_Tree->Branch("var_tauEta",&var_tauEta);
  TMVA_Tree->Branch("var_ntracks",&var_ntracks);
  TMVA_Tree->Branch("var_relPt",&var_relPt);
  TMVA_Tree->Branch("var_isoMax",&var_isoMax);
  
  TMVA_Tree->Branch("var_MaxdeltaMuZ",&var_MaxdeltaMuZ);
  TMVA_Tree->Branch("var_MindeltaMuZ",&var_MindeltaMuZ);
  TMVA_Tree->Branch("var_deltaMuZ12",&var_deltaMuZ12);
  TMVA_Tree->Branch("var_deltaMuZ13",&var_deltaMuZ13);
  TMVA_Tree->Branch("var_deltaMuZ23",&var_deltaMuZ23);

  TMVA_Tree->Branch("var_maxMuonsDca",&var_maxMuonsDca);
  TMVA_Tree->Branch("var_minMuonsDca",&var_minMuonsDca);
  TMVA_Tree->Branch("var_nsv",&var_nsv);

  TMVA_Tree->Branch("var_VertexMu1D0SigPVReco",&var_VertexMu1D0SigPVReco);
  TMVA_Tree->Branch("var_VertexMu2D0SigPVReco",&var_VertexMu2D0SigPVReco);
  TMVA_Tree->Branch("var_VertexMu3D0SigPVReco",&var_VertexMu3D0SigPVReco);

  TMVA_Tree->Branch("var_MaxD0SigPV",&var_MaxD0SigPV);
  TMVA_Tree->Branch("var_MinD0SigPV",&var_MinD0SigPV);

  TMVA_Tree->Branch("var_VertexMu1D0SigBSReco",&var_VertexMu1D0SigBSReco);
  TMVA_Tree->Branch("var_VertexMu2D0SigBSReco",&var_VertexMu2D0SigBSReco);
  TMVA_Tree->Branch("var_VertexMu3D0SigBSReco",&var_VertexMu3D0SigBSReco);

  TMVA_Tree->Branch("var_MaxD0SigBS",&var_MaxD0SigBS);
  TMVA_Tree->Branch("var_MinD0SigBS",&var_MinD0SigBS);

  TMVA_Tree->Branch("var_VertexMu1D0SigSVReco",&var_VertexMu1D0SigSVReco);
  TMVA_Tree->Branch("var_VertexMu2D0SigSVReco",&var_VertexMu2D0SigSVReco);
  TMVA_Tree->Branch("var_VertexMu3D0SigSVReco",&var_VertexMu3D0SigSVReco);

  TMVA_Tree->Branch("var_MaxD0SigSV",&var_MaxD0SigSV);
  TMVA_Tree->Branch("var_MinD0SigSV",&var_MinD0SigSV);

  TMVA_Tree->Branch("var_VertexMu1DistanceToSV",&var_VertexMu1DistanceToSV);
  TMVA_Tree->Branch("var_VertexMu2DistanceToSV",&var_VertexMu2DistanceToSV);
  TMVA_Tree->Branch("var_VertexMu3DistanceToSV",&var_VertexMu3DistanceToSV);

  TMVA_Tree->Branch("var_MaxMuDistanceToSV",&var_MaxMuDistanceToSV);
  TMVA_Tree->Branch("var_MinMuDistanceToSV",&var_MinMuDistanceToSV);


  TMVA_Tree->Branch("var_MinMuon_chi2LocalPosition",&var_MinMuon_chi2LocalPosition);
  TMVA_Tree->Branch("var_MaxMuon_chi2LocalPosition",&var_MaxMuon_chi2LocalPosition);

  TMVA_Tree->Branch("var_Muon1_chi2LocalPosition",&var_Muon1_chi2LocalPosition);
  TMVA_Tree->Branch("var_Muon2_chi2LocalPosition",&var_Muon2_chi2LocalPosition);
  TMVA_Tree->Branch("var_Muon3_chi2LocalPosition",&var_Muon3_chi2LocalPosition);

  TMVA_Tree->Branch("var_MinMuon_chi2LocalMomentum",&var_MinMuon_chi2LocalMomentum);
  TMVA_Tree->Branch("var_MaxMuon_chi2LocalMomentum",&var_MaxMuon_chi2LocalMomentum);

  TMVA_Tree->Branch("var_Muon1_chi2LocalMomentum",&var_Muon1_chi2LocalMomentum);
  TMVA_Tree->Branch("var_Muon2_chi2LocalMomentum",&var_Muon2_chi2LocalMomentum);
  TMVA_Tree->Branch("var_Muon3_chi2LocalMomentum",&var_Muon3_chi2LocalMomentum);

  TMVA_Tree->Branch("var_MintrkKink",&var_MintrkKink);
  TMVA_Tree->Branch("var_MaxtrkKink",&var_MaxtrkKink);
  TMVA_Tree->Branch("var_MinglbKink",&var_MinglbKink);
  TMVA_Tree->Branch("var_MaxglbKink",&var_MaxglbKink);

  TMVA_Tree->Branch("var_MuonglbkinkSum",&var_MuonglbkinkSum);

  TMVA_Tree->Branch("var_MaxVertexPairQuality",&var_MaxVertexPairQuality);
  TMVA_Tree->Branch("var_MinVertexPairQuality",&var_MinVertexPairQuality);

  TMVA_Tree->Branch("var_Iso02",&var_Iso02);
  TMVA_Tree->Branch("var_Iso04",&var_Iso04);
  TMVA_Tree->Branch("var_Iso06",&var_Iso06);
  TMVA_Tree->Branch("var_Iso08",&var_Iso08);
  TMVA_Tree->Branch("var_Iso1",&var_Iso1);
  TMVA_Tree->Branch("var_Iso12",&var_Iso12);

  TMVA_Tree->Branch("var_Iso02Mu1",&var_Iso02Mu1);
  TMVA_Tree->Branch("var_Iso04Mu1",&var_Iso04Mu1);
  TMVA_Tree->Branch("var_Iso06Mu1",&var_Iso06Mu1);
  TMVA_Tree->Branch("var_Iso08Mu1",&var_Iso08Mu1);
  TMVA_Tree->Branch("var_Iso1Mu1",&var_Iso1Mu1);
  TMVA_Tree->Branch("var_Iso12Mu1",&var_Iso12Mu1);


  TMVA_Tree->Branch("var_Iso02Mu2",&var_Iso02Mu2);
  TMVA_Tree->Branch("var_Iso04Mu2",&var_Iso04Mu2);
  TMVA_Tree->Branch("var_Iso06Mu2",&var_Iso06Mu2);
  TMVA_Tree->Branch("var_Iso08Mu2",&var_Iso08Mu2);
  TMVA_Tree->Branch("var_Iso1Mu2",&var_Iso1Mu2);
  TMVA_Tree->Branch("var_Iso12Mu2",&var_Iso12Mu2);

  TMVA_Tree->Branch("var_Iso02Mu3",&var_Iso02Mu3);
  TMVA_Tree->Branch("var_Iso04Mu3",&var_Iso04Mu3);
  TMVA_Tree->Branch("var_Iso06Mu3",&var_Iso06Mu3);
  TMVA_Tree->Branch("var_Iso08Mu3",&var_Iso08Mu3);
  TMVA_Tree->Branch("var_Iso1Mu3",&var_Iso1Mu3);
  TMVA_Tree->Branch("var_Iso12Mu3",&var_Iso12Mu3);

  TMVA_Tree->Branch("var_Iso08MuMax",&var_Iso08MuMax);
  TMVA_Tree->Branch("var_Iso08MuMin",&var_Iso08MuMin);

  TMVA_Tree->Branch("var_NtracksClose",&var_NtracksClose);
  TMVA_Tree->Branch("var_NTracksCloseToPV",&var_NTracksCloseToPV);
  TMVA_Tree->Branch("var_NTracksCloseToPVTauDR",&var_NTracksCloseToPVTauDR);




  TMVA_Tree->Branch("var_Muon1Pt",&var_Muon1Pt);
  TMVA_Tree->Branch("var_Muon2Pt",&var_Muon2Pt);
  TMVA_Tree->Branch("var_Muon3Pt",&var_Muon3Pt);


  TMVA_Tree->Branch("var_Muon1P",&var_Muon1P);
  TMVA_Tree->Branch("var_Muon2P",&var_Muon2P);
  TMVA_Tree->Branch("var_Muon3P",&var_Muon3P);


  TMVA_Tree->Branch("var_MindcaTrackSV",&var_MindcaTrackSV);
  TMVA_Tree->Branch("var_Mu1TrackMass",&var_Mu1TrackMass);
  TMVA_Tree->Branch("var_Mu2TrackMass",&var_Mu2TrackMass);
  TMVA_Tree->Branch("var_Mu3TrackMass",&var_Mu3TrackMass);
  TMVA_Tree->Branch("var_dcaTrackPV",&var_dcaTrackPV);


  TMVA_Tree->Branch("var_MinMatchedStations",&var_MinMatchedStations);
  TMVA_Tree->Branch("var_MaxMatchedStations",&var_MaxMatchedStations);
  TMVA_Tree->Branch("var_Mu1MatchedStations",&var_Mu1MatchedStations);
  TMVA_Tree->Branch("var_Mu2MatchedStations",&var_Mu2MatchedStations);
  TMVA_Tree->Branch("var_Mu3MatchedStations",&var_Mu3MatchedStations);


  TMVA_Tree->Branch("var_MinMuon_numberOfChambers",&var_MinMuon_numberOfChambers);
  TMVA_Tree->Branch("var_MaxMuon_numberOfChambers",&var_MaxMuon_numberOfChambers);
  TMVA_Tree->Branch("var_Mu1Muon_numberOfChambers",&var_Mu1Muon_numberOfChambers);
  TMVA_Tree->Branch("var_Mu2Muon_numberOfChambers",&var_Mu2Muon_numberOfChambers);
  TMVA_Tree->Branch("var_Mu3Muon_numberOfChambers",&var_Mu3Muon_numberOfChambers);


  TMVA_Tree->Branch("var_MinMuon_numberOfMatches",&var_MinMuon_numberOfMatches);
  TMVA_Tree->Branch("var_MaxMuon_numberOfMatches",&var_MaxMuon_numberOfMatches);
  TMVA_Tree->Branch("var_Mu1Muon_numberOfMatches",&var_Mu1Muon_numberOfMatches);
  TMVA_Tree->Branch("var_Mu2Muon_numberOfMatches",&var_Mu2Muon_numberOfMatches);
  TMVA_Tree->Branch("var_Mu3Muon_numberOfMatches",&var_Mu3Muon_numberOfMatches);

  TMVA_Tree->Branch("var_Muon1LooseId",&var_Muon1LooseId);
  TMVA_Tree->Branch("var_Muon1MediumId",&var_Muon1MediumId);
  TMVA_Tree->Branch("var_Muon1TightId",&var_Muon1TightId);

  TMVA_Tree->Branch("var_Muon2LooseId",&var_Muon2LooseId);
  TMVA_Tree->Branch("var_Muon2MediumId",&var_Muon2MediumId);
  TMVA_Tree->Branch("var_Muon2TightId",&var_Muon2TightId);

  TMVA_Tree->Branch("var_Muon3LooseId",&var_Muon3LooseId);
  TMVA_Tree->Branch("var_Muon3MediumId",&var_Muon3MediumId);
  TMVA_Tree->Branch("var_Muon3TightId",&var_Muon3TightId);

  TMVA_Tree->Branch("var_Muon1PFIsoLoose",&var_Muon1PFIsoLoose);
  TMVA_Tree->Branch("var_Muon1PFIsoMedium",&var_Muon1PFIsoMedium);
  TMVA_Tree->Branch("var_Muon1PFIsoTight",&var_Muon1PFIsoTight);
  TMVA_Tree->Branch("var_Muon1PFIsoVTight",&var_Muon1PFIsoVTight);

  TMVA_Tree->Branch("var_Muon2PFIsoLoose",&var_Muon2PFIsoLoose);
  TMVA_Tree->Branch("var_Muon2PFIsoMedium",&var_Muon2PFIsoMedium);
  TMVA_Tree->Branch("var_Muon2PFIsoTight",&var_Muon2PFIsoTight);
  TMVA_Tree->Branch("var_Muon2PFIsoVTight",&var_Muon2PFIsoVTight);

  TMVA_Tree->Branch("var_Muon3PFIsoLoose",&var_Muon3PFIsoLoose);
  TMVA_Tree->Branch("var_Muon3PFIsoMedium",&var_Muon3PFIsoMedium);
  TMVA_Tree->Branch("var_Muon3PFIsoTight",&var_Muon3PFIsoTight);
  TMVA_Tree->Branch("var_Muon3PFIsoVTight",&var_Muon3PFIsoVTight);

  TMVA_Tree->Branch("var_Muon1ImpactAngle",&var_Muon1ImpactAngle);
  TMVA_Tree->Branch("var_Muon2ImpactAngle",&var_Muon2ImpactAngle);
  TMVA_Tree->Branch("var_Muon3ImpactAngle",&var_Muon3ImpactAngle);
  TMVA_Tree->Branch("var_MinMuonImpactAngle",&var_MinMuonImpactAngle);
  TMVA_Tree->Branch("var_MaxMuonImpactAngle",&var_MaxMuonImpactAngle);

  TMVA_Tree->Branch("var_Muon1DetID",&var_Muon1DetID);
  TMVA_Tree->Branch("var_Muon2DetID",&var_Muon2DetID);
  TMVA_Tree->Branch("var_Muon3DetID",&var_Muon3DetID);


  TMVA_Tree->Branch("var_mass12",&var_mass12);
  TMVA_Tree->Branch("var_mass13",&var_mass13);

  TMVA_Tree->Branch("var_mass12_dRsorting",&var_mass12_dRsorting);
  TMVA_Tree->Branch("var_mass13_drSorting",&var_mass13_dRsorting);


  TMVA_Tree->Branch("var_VertexMu1D0DistSVReco",&var_VertexMu1D0DistSVReco);
  TMVA_Tree->Branch("var_VertexMu2D0DistSVReco",&var_VertexMu2D0DistSVReco);
  TMVA_Tree->Branch("var_VertexMu3D0DistSVReco",&var_VertexMu3D0DistSVReco);
  TMVA_Tree->Branch("var_MaxD0DistSV",&var_MaxD0DistSV);
  TMVA_Tree->Branch("var_MinD0DistSV",&var_MinD0DistSV);


  TMVA_Tree->Branch("var_VertexMu1DZDistSVReco",&var_VertexMu1DZDistSVReco);
  TMVA_Tree->Branch("var_VertexMu2DZDistSVReco",&var_VertexMu2DZDistSVReco);
  TMVA_Tree->Branch("var_VertexMu3DZDistSVReco",&var_VertexMu3DZDistSVReco);
  TMVA_Tree->Branch("var_MaxDZDistSV",&var_MaxDZDistSV);
  TMVA_Tree->Branch("var_MinDZDistSV",&var_MinDZDistSV);

  TMVA_Tree->Branch("var_IsoPhiKKMass_Mu1",&var_IsoPhiKKMass_Mu1 );
  TMVA_Tree->Branch("var_IsoKStarMass_Mu1",&var_IsoKStarMass_Mu1 );
  TMVA_Tree->Branch("var_IsoMuMuMass_Mu1",&var_IsoMuMuMass_Mu1 );

  TMVA_Tree->Branch("var_IsoPhiKKMass_Mu2",&var_IsoPhiKKMass_Mu2 );
  TMVA_Tree->Branch("var_IsoKStarMass_Mu2",&var_IsoKStarMass_Mu2 );
  TMVA_Tree->Branch("var_IsoMuMuMass_Mu2",&var_IsoMuMuMass_Mu2 );

  TMVA_Tree->Branch("var_IsoPhiKKMass_Mu3",&var_IsoPhiKKMass_Mu3 );
  TMVA_Tree->Branch("var_IsoKStarMass_Mu3",&var_IsoKStarMass_Mu3 );
  TMVA_Tree->Branch("var_IsoMuMuMass_Mu3",&var_IsoMuMuMass_Mu3 );


  TMVA_Tree->Branch("var_Vertex2muTrkKF",&var_Vertex2muTrkKF );
  TMVA_Tree->Branch("var_Dist2muTrkKF3Mu",&var_Dist2muTrkKF3Mu );
  TMVA_Tree->Branch("var_VertexQualitySeparator",&var_VertexQualitySeparator );

  TMVA_Tree->Branch("var_BvsDSeprator",&var_BvsDSeprator );


  TMVA_Tree->Branch("var_cTheta_TRF_SSSS",&var_cTheta_TRF_SSSS );
  TMVA_Tree->Branch("var_cTheta_TRF_OSSS",&var_cTheta_TRF_OSSS );
  TMVA_Tree->Branch("var_cTheta_MuonOS_TauPol_TRF",&var_cTheta_MuonOS_TauPol_TRF );
  TMVA_Tree->Branch("var_pTMu1OverMass_TRF",&var_pTMu1OverMass_TRF );
  TMVA_Tree->Branch("var_costheta_TRF_SSSS",&var_costheta_TRF_SSSS );
  TMVA_Tree->Branch("var_costheta_TRF_OSS",&var_costheta_TRF_OSS );
  TMVA_Tree->Branch("var_OSSS1Angle_TRF",&var_OSSS1Angle_TRF );
  TMVA_Tree->Branch("var_OSSS2Angle_TRF",&var_OSSS2Angle_TRF );

  TMVA_Tree->Branch("var_OSSS1Angle_RRF",&var_OSSS1Angle_RRF );
  TMVA_Tree->Branch("var_OSSS2Angle_RRF",&var_OSSS2Angle_RRF );




  // -----------------
  
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==HLT)                cut.at(HLT)=1;
    if(i==L1)                 cut.at(L1)=1;
    if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
    if(i==Mu1PtCut)           cut.at(Mu1PtCut)=3.0;
    if(i==Mu2PtCut)           cut.at(Mu2PtCut)=2.0;
    if(i==Mu3PtCut)           cut.at(Mu3PtCut)=2.0;
    if(i==MuonID)             cut.at(MuonID)=1;
    if(i==TRKLWithM)          cut.at(TRKLWithM)=1;
    if(i==PhiVeto1)           cut.at(PhiVeto1)=0; // defined below
    if(i==OmegaVeto1)         cut.at(OmegaVeto1)=0; // defined below
    if(i==PhiVeto2)           cut.at(PhiVeto2)=0; // defined below
    if(i==OmegaVeto2)         cut.at(OmegaVeto2)=0; // defined below
    if(i==TriggerMatch)       cut.at(TriggerMatch)=1;
    if(i==ThreeMuMass)        cut.at(ThreeMuMass)=1;// true for MC and mass side band for data
    if(i==CutCategory)        cut.at(CutCategory)=2;// 
  }
  
  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    if(i==L1){
      title.at(i)="Pass L1";
      hlabel="L1";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_L1_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_L1_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==HLT){
      title.at(i)="Pass HLT";
      hlabel="HLT";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLT_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLT_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==SignalCandidate){
      title.at(i)="is signal candidate";
      hlabel="is 3mu candidate";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Mu1PtCut){
      title.at(i)="Mu1 Pt $>$ ";
      title.at(i)+=cut.at(Mu1PtCut);
      title.at(i)+=" GeV";
      hlabel="Muon1 PT, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1PtCut_",htitle,40,2,25,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1PtCut_",htitle,40,2,25,hlabel,"Events"));
    }
    else if(i==Mu2PtCut){
      title.at(i)="Mu2 Pt $>$ ";
      title.at(i)+=cut.at(Mu2PtCut);
      title.at(i)+=" GeV";
      hlabel="Muon2 PT, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2PtCut_",htitle,40,2,20,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2PtCut_",htitle,40,2,20,hlabel,"Events"));
    }
    else if(i==Mu3PtCut){
      title.at(i)="Mu3 Pt $>$ ";
      title.at(i)+=cut.at(Mu3PtCut);
      title.at(i)+=" GeV";
      hlabel="Muon3 PT, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu3PtCut_",htitle,40,2,15,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu3PtCut_",htitle,40,2,15,hlabel,"Events"));
    }
    else if(i==MuonID){
      title.at(i)="All mu pass ID";
      hlabel="pass MuonID";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonID_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }


    else if(i==TRKLWithM){
      title.at(i)="TrkLayerWithM";
      hlabel="pass TRKLWithM";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TRKLWithM_",htitle,20,-0.5,19.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TRKLWithM_",htitle,20,-0.5,19.5,hlabel,"Events"));
    }

    else if(i==PhiVeto1){
      title.at(i)="phi mass veto";
      hlabel="Phi mass Veto 1 pair, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PhiVeto1_",htitle,40,0.8,1.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PhiVeto1_",htitle,40,0.8,1.2,hlabel,"Events"));
    }
    else if(i==OmegaVeto1){
      title.at(i)="omega mass veto";
      hlabel="Omega mass veto 1 pair, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OmegaVeto1_",htitle,50,0.4,0.9,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OmegaVeto1_",htitle,50,0.4,0.9,hlabel,"Events"));
    }
    else if(i==PhiVeto2){
      title.at(i)="phi mass veto";
      hlabel="Phi mass Veto 2nd pair, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PhiVeto2_",htitle,40,0.8,1.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PhiVeto2_",htitle,40,0.8,1.2,hlabel,"Events"));
    }
    else if(i==OmegaVeto2){
      title.at(i)="omega mass veto";
      hlabel="Omega mass veto 2nd pair, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OmegaVeto2_",htitle,50,0.4,0.9,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OmegaVeto2_",htitle,50,0.4,0.9,hlabel,"Events"));
    }
    else if(i==TriggerMatch){
      title.at(i)="Trigger Matching";
      hlabel="trigger matching dR";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==ThreeMuMass){
      title.at(i)="Tau Mass";
      hlabel="three mu mass, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ThreeMuMass_",htitle,80,1.4,2.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ThreeMuMass_",htitle,80,1.4,2.2,hlabel,"Events"));
    }

   else if(i==CutCategory){
      title.at(i)="CutCategory";
      hlabel="CutCategory";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_CutCategory_",htitle,3,0.5,3.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_CutCategory_",htitle,3,0.5,3.5,hlabel,"Events"));
    }
  } 
  
      Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove


      nPVx = HConfig.GetTH1D(Name+"_nPVx","nPVx",70,-0.5,69.5," # primary vertices","Events");

      Muon1Pt =HConfig.GetTH1D(Name+"_Muon1Pt","Muon1Pt",50,0,25,"  #mu_{1} p_{T}, GeV","Events");
      Muon2Pt =HConfig.GetTH1D(Name+"_Muon2Pt","Muon2Pt",50,0,20,"  #mu_{2} p_{T}, GeV","Events");
      Muon3Pt =HConfig.GetTH1D(Name+"_Muon3Pt","Muon3Pt",50,0,15,"  #mu_{3} p_{T}, GeV","Events");


      Muon1P =HConfig.GetTH1D(Name+"_Muon1P","Muon1P",50,0,25,"  #mu_{1} p, GeV","Events");
      Muon2P =HConfig.GetTH1D(Name+"_Muon2P","Muon2P",50,0,20,"  #mu_{2} p, GeV","Events");
      Muon3P =HConfig.GetTH1D(Name+"_Muon3P","Muon3P",50,0,15,"  #mu_{3} p, GeV","Events");
      
      TauEta=HConfig.GetTH1D(Name+"_TauEta","TauEta",30,-2.6,2.6,"#eta(#tau)","Events");
      TauPt =HConfig.GetTH1D(Name+"_TauPt","TauPt",30,0,50,"p_{T}(#tau), GeV","Events");
      TauP =HConfig.GetTH1D(Name+"_TauP","TauP",40,0,70,"|p|(#tau), GeV","Events");
      
      VertexChi2KF=HConfig.GetTH1D(Name+"_VertexChi2KF","VertexChi2KF",50,0,20,"KF vertex #chi^{2}","Events");
      Vertex2muTrkKF=HConfig.GetTH1D(Name+"_Vertex2muTrkKF","Vertex2muTrkKF",50,-1.5,28.5,"KF vertex #chi^{2} 2#mu + tr","Events");

      Vertex2muTrkKFChi2=HConfig.GetTH1D(Name+"_Vertex2muTrkKFChi2","iso vertex chi2",30,-1.5,28.5,"2#mu+iso track #chi^{2}","Events");
      Vertex2muTrkKFToSignalVertexChi2=HConfig.GetTH1D(Name+"_Vertex2muTrkKFToSignalVertexChi2","Vertex Ratios",7,-1.5,5.5,"2#mu+iso track #chi^{2} / 3#mu #chi^{2}","Events");
      Vertex2muTrkKFToSignalVertexDistance=HConfig.GetTH1D(Name+"_Vertex2muTrkKFToSignalVertexDistance","signal - iso vertex distance",50,0,2,"3#mu Vertex - 3#track Vertex dist, cm","Events");




      Dist2muTrkKF3Mu=HConfig.GetTH1D(Name+"_Dist2muTrkKF3Mu","Dist2muTrkKF3Mu",50,0,2,"Dist( 2#mu + tr - 3#mu), cm","Events");
      MuonglbkinkSum  =HConfig.GetTH1D(Name+"_MuonglbkinkSum","MuonglbkinkSum",50,0.,50," #sum  #mu glb kink #chi^{2}","Events");
      FLSignificance=HConfig.GetTH1D(Name+"_FLSignificance","FLSignificance",60,0,60,"PV - SV distance  significance","Events");
      FL=HConfig.GetTH1D(Name+"_FL","FL",60,0,1.5,"Flight length ,cm","Events");

      SVPVTauDirAngle=HConfig.GetTH1D(Name+"_SVPVTauDirAngle","SVPVTauDirAngle",60,0,0.25,"Angle btw #vec{SV}-#vec{PV} and #vec{#tau}, rad","Events");
      Muon_segmentCompatibility_mu1  = HConfig.GetTH1D(Name+"_Muon_segmentCompatibility_mu1","Muon_segmentCompatibility_mu1",50,0.,1,"Inner Track and muon segment match  #mu_{1} ","Events");
      Muon_segmentCompatibility_mu2  = HConfig.GetTH1D(Name+"_Muon_segmentCompatibility_mu2","Muon_segmentCompatibility_mu2",50,0.,1,"Inner Track and muon segment match  #mu_{2} ","Events");
      Muon_segmentCompatibility_mu3  = HConfig.GetTH1D(Name+"_Muon_segmentCompatibility_mu3","Muon_segmentCompatibility_mu3",50,0.,1,"Inner Track and muon segment match  #mu_{3} ","Events");
      
      Muon_segmentCompatibility_min  = HConfig.GetTH1D(Name+"_Muon_segmentCompatibility_min","Muon_segmentCompatibility_min",50,0.,1,"Inner Track and muon segment match min ","Events");
      Muon_segmentCompatibility_max  = HConfig.GetTH1D(Name+"_Muon_segmentCompatibility_max","Muon_segmentCompatibility_max",50,0.,1,"Inner Track and muon segment match max ","Events");
      
      Muon_ECALCompatibility_mu1  = HConfig.GetTH1D(Name+"_Muon_ECALCompatibility_mu1","Muon_ECALCompatibility_mu1",50,0.,1,"MIP Likelihood  #mu_{1} ","Events");
      Muon_ECALCompatibility_mu2  = HConfig.GetTH1D(Name+"_Muon_ECALCompatibility_mu2","Muon_ECALCompatibility_mu2",50,0.,1,"MIP Likelihood  #mu_{2} ","Events");
      Muon_ECALCompatibility_mu3  = HConfig.GetTH1D(Name+"_Muon_ECALCompatibility_mu3","Muon_ECALCompatibility_mu3",50,0.,1,"MIP Likelihood  #mu_{3} ","Events");
      
      Muon_ECALCompatibility_min  = HConfig.GetTH1D(Name+"_Muon_ECALCompatibility_min","Muon_ECALCompatibility_min",50,0.,1,"MIP Likelihood min ","Events");
      Muon_ECALCompatibility_max  = HConfig.GetTH1D(Name+"_Muon_ECALCompatibility_max","Muon_ECALCompatibility_max",50,0.,1,"MIP Likelihood max ","Events");
      
      Isolation_NTracks=HConfig.GetTH1D(Name+"_Isolation_NTracks","Isolation_NTracks",10,-0.5,9.5,"N tracks","Events");
      Isolation_RelPt=HConfig.GetTH1D(Name+"_Isolation_RelPt","Isolation_RelPt",50,0,1,"relative p_{T}","Events");
      Isolation_maxdxy=HConfig.GetTH1D(Name+"_Isolation_maxdxy","Isolation_maxdxy",40,0,15,"Iso maximum transversely displaced track","Events");

      EventMassResolution_PtEtaPhi = HConfig.GetTH1D(Name+"_EventMassResolution_PtEtaPhi","EventMassResolution_PtEtaPhi",50,0,0.02,"#frac{#Delta m}{m} (ptEtaPhi)","Events");
      EMR_tau_eta  =HConfig.GetTH2D(Name+"_EMR_tau_eta","EMR vs eta",50,0,0.02,50,0,2.5,"3#mu mass, GeV","#eta_{#tau}");

      VertexMu1D0SigPVReco=HConfig.GetTH1D(Name+"_VertexMu1D0SigPVReco","VertexMu1D0SigPVReco",50,0,15,"#mu_{1} - PV transverse distance significance","Events");
      VertexMu2D0SigPVReco=HConfig.GetTH1D(Name+"_VertexMu2D0SigPVReco","VertexMu2D0SigPVReco",50,0,15,"#mu_{2} - PV transverse distance significance","Events");
      VertexMu3D0SigPVReco=HConfig.GetTH1D(Name+"_VertexMu3D0SigPVReco","VertexMu3D0SigPVReco",50,0,15,"#mu_{3} - PV transverse distance significance","Events");


      MaxD0SigPV=HConfig.GetTH1D(Name+"_MaxD0SigPV","MaxD0SigPV",50,0,15,"Max Transverse Impact significance w.r.t PV","");
      MinD0SigPV=HConfig.GetTH1D(Name+"_MinD0SigPV","MinD0SigPV",50,0,10,"Min Transverse Impact significance w.r.t PV","");
   

      VertexMu1D0SigBSReco=HConfig.GetTH1D(Name+"_VertexMu1D0SigBSReco","VertexMu1D0SigBSReco",50,0,15,"#mu_{1} - BS transverse distance significance","Events");
      VertexMu2D0SigBSReco=HConfig.GetTH1D(Name+"_VertexMu2D0SigBSReco","VertexMu2D0SigBSReco",50,0,15,"#mu_{2} - BS transverse distance significance","Events");
      VertexMu3D0SigBSReco=HConfig.GetTH1D(Name+"_VertexMu3D0SigBSReco","VertexMu3D0SigBSReco",50,0,15,"#mu_{3} - BS transverse distance significance","Events");
 
   
      MaxD0SigBS=HConfig.GetTH1D(Name+"_MaxD0SigBS","MaxD0SigBS",50,0,15,"Max Transverse Impact significance w.r.t BS","");
      MinD0SigBS=HConfig.GetTH1D(Name+"_MinD0SigBS","MinD0SigBS",50,0,15,"Min Transverse Impact significance w.r.t BS","");
   

      VertexMu1D0SigSVReco=HConfig.GetTH1D(Name+"_VertexMu1D0SigSVReco","VertexMu1D0SigSVReco",50,0,5,"#mu_{1} - SV transverse distance significance","Events");
      VertexMu2D0SigSVReco=HConfig.GetTH1D(Name+"_VertexMu2D0SigSVReco","VertexMu2D0SigSVReco",50,0,5,"#mu_{2} - SV transverse distance significance","Events");
      VertexMu3D0SigSVReco=HConfig.GetTH1D(Name+"_VertexMu3D0SigSVReco","VertexMu3D0SigSVReco",50,0,5,"#mu_{3} - SV transverse distance significance","Events");
 
      VertexMu1D0DistSVReco=HConfig.GetTH1D(Name+"_VertexMu1D0DistSVReco","VertexMu1D0DistSVReco",50,0,0.02,"#mu_{1} - SV transverse distance, cm","Events");
      VertexMu2D0DistSVReco=HConfig.GetTH1D(Name+"_VertexMu2D0DistSVReco","VertexMu2D0DistSVReco",50,0,0.02,"#mu_{2} - SV transverse distance, cm","Events");
      VertexMu3D0DistSVReco=HConfig.GetTH1D(Name+"_VertexMu3D0DistSVReco","VertexMu3D0DistSVReco",50,0,0.02,"#mu_{3} - SV transverse distance, cm","Events");

      MaxD0DistSV=HConfig.GetTH1D(Name+"_MaxD0DistSV","MaxD0DistSV",50,0,0.02,"#mu - SV Max transverse distance, cm","Events");
      MinD0DistSV=HConfig.GetTH1D(Name+"_MinD0DistSV","MinD0DistSV",50,0,0.01,"#mu - SV Min transverse distance, cm","Events");

      VertexMu1DZDistSVReco=HConfig.GetTH1D(Name+"_VertexMu1DZDistSVReco","VertexMu1DZDistSVReco",50,0,0.02,"#mu_{1} - SV long. distance, cm","Events");
      VertexMu2DZDistSVReco=HConfig.GetTH1D(Name+"_VertexMu2DZDistSVReco","VertexMu2DZDistSVReco",50,0,0.02,"#mu_{2} - SV long. distance, cm","Events");
      VertexMu3DZDistSVReco=HConfig.GetTH1D(Name+"_VertexMu3DZDistSVReco","VertexMu3DZDistSVReco",50,0,0.02,"#mu_{3} - SV long. distance, cm","Events");

      MaxDZDistSV=HConfig.GetTH1D(Name+"_MaxDZDistSV","MaxDZDistSV",50,0,0.02,"#mu - SV max long. distance, cm","Events");
      MinDZDistSV=HConfig.GetTH1D(Name+"_MinDZDistSV","MinDZDistSV",50,0,0.02,"#mu - SV min long. distance, cm","Events");

      VertexMu1DistanceToSV=HConfig.GetTH1D(Name+"_VertexMu1DistanceToSV","VertexMu1DistanceToSV",50,0,0.03,"#mu_{1} - SV  distance, cm","Events");
      VertexMu2DistanceToSV=HConfig.GetTH1D(Name+"_VertexMu2DistanceToSV","VertexMu2DistanceToSV",50,0,0.03,"#mu_{2} - SV  distance, cm","Events");
      VertexMu3DistanceToSV=HConfig.GetTH1D(Name+"_VertexMu3DistanceToSV","VertexMu3DistanceToSV",50,0,0.03,"#mu_{3} - SV  distance, cm","Events");

      MaxMuDistanceToSV=HConfig.GetTH1D(Name+"_MaxMuDistanceToSV","MaxMuDistanceToSV",50,0,0.02,"#mu - SV max  distance, cm","Events");
      MinMuDistanceToSV=HConfig.GetTH1D(Name+"_MinMuDistanceToSV","MinMuDistanceToSV",50,0,0.02,"#mu - SV min  distance, cm","Events");



   
      MaxD0SigSV=HConfig.GetTH1D(Name+"_MaxD0SigSV","MaxD0SigSV",50,0,5,"Max Transverse Impact significance w.r.t SV","");
      MinD0SigSV=HConfig.GetTH1D(Name+"_MinD0SigSV","MinD0SigSV",50,0,2,"Min Transverse Impact significance w.r.t SV","");
      
      MintrkKink= HConfig.GetTH1D(Name+"_MintrkKink","MintrkKink",30,0,30,"Min Tracker Kink","");
      MaxtrkKink= HConfig.GetTH1D(Name+"_MaxtrkKink","MaxtrkKink",30,0,30,"Max Tracker Kink","");
      MinglbKink= HConfig.GetTH1D(Name+"_MinglbKink","MinglbKink",30,0,30,"Min Global  Kink","");
      MaxglbKink= HConfig.GetTH1D(Name+"_MaxglbKink","MaxglbKink",30,0,5,"Max Global  Kink","");

      MinMuon_chi2LocalPosition=HConfig.GetTH1D(Name+"_MinMuon_chi2LocalPosition","MinMuon_chi2LocalPosition",50,0,2,"Min Inner/Outer track #chi^{2}","");
      MaxMuon_chi2LocalPosition=HConfig.GetTH1D(Name+"_MaxMuon_chi2LocalPosition","MaxMuon_chi2LocalPosition",50,0,30,"Max Inner/Outer track #chi^{2}","");

      Muon1_chi2LocalPosition=HConfig.GetTH1D(Name+"_Muon1_chi2LocalPosition","Muon1_chi2LocalPosition",50,0,30,"#mu_{1} Inner/Outer track #chi^{2}","");
      Muon2_chi2LocalPosition=HConfig.GetTH1D(Name+"_Muon2_chi2LocalPosition","Muon2_chi2LocalPosition",50,0,30,"#mu_{2} Inner/Outer track #chi^{2}","");
      Muon3_chi2LocalPosition=HConfig.GetTH1D(Name+"_Muon3_chi2LocalPosition","Muon3_chi2LocalPosition",50,0,15,"#mu_{3} Inner/Outer track #chi^{2}","");

      
      MinMuon_chi2LocalMomentum=HConfig.GetTH1D(Name+"_MinMuon_chi2LocalMomentum","MinMuon_chi2LocalMomentum",50,0,5,"Min Inner/Outer track momentum #chi^{2}","");
      MaxMuon_chi2LocalMomentum=HConfig.GetTH1D(Name+"_MaxMuon_chi2LocalMomentum","MaxMuon_chi2LocalMomentum",50,0,50,"Max Inner/Outer track momentum #chi^{2}","");

      Muon1_chi2LocalMomentum=HConfig.GetTH1D(Name+"_Muon1_chi2LocalMomentum","Muon1_chi2LocalMomentum",50,0,50,"#mu_{1} Inner/Outer track momentum #chi^{2}","");
      Muon2_chi2LocalMomentum=HConfig.GetTH1D(Name+"_Muon2_chi2LocalMomentum","Muon2_chi2LocalMomentum",50,0,50,"#mu_{2} Inner/Outer track momentum #chi^{2}","");
      Muon3_chi2LocalMomentum=HConfig.GetTH1D(Name+"_Muon3_chi2LocalMomentum","Muon3_chi2LocalMomentum",50,0,50,"#mu_{3} Inner/Outer track momentum #chi^{2}","");


      MinMatchedStations=HConfig.GetTH1D(Name+"_MinMatchedStations","MinMatchedStations",10,-0.5,9.5,"min matched stations","");
      MaxMatchedStations=HConfig.GetTH1D(Name+"_MaxMatchedStations","MaxMatchedStations",10,-0.5,9.5,"max matched stations ","");
      Mu1MatchedStations=HConfig.GetTH1D(Name+"_Mu1MatchedStations","Mu1MatchedStations",10,-0.5,9.5,"#mu_{1} matched stations","");
      Mu2MatchedStations=HConfig.GetTH1D(Name+"_Mu2MatchedStations","Mu2MatchedStations",10,-0.5,9.5,"#mu_{2} matched stations","");
      Mu3MatchedStations=HConfig.GetTH1D(Name+"_Mu3MatchedStations","Mu3MatchedStations",10,-0.5,9.5,"#mu_{3} matched stations","");


      MinMuon_numberOfChambers=HConfig.GetTH1D(Name+"_MinMuon_numberOfChambers","MinMuon_numberOfChambers",10,-0.5,9.5,"min numberOfChambers","");
      MaxMuon_numberOfChambers=HConfig.GetTH1D(Name+"_MaxMuon_numberOfChambers","MaxMuon_numberOfChambers",10,-0.5,9.5,"max numberOfChambers","");
      Mu1Muon_numberOfChambers=HConfig.GetTH1D(Name+"_Mu1Muon_numberOfChambers","Mu1Muon_numberOfChambers",10,-0.5,9.5,"#mu_{1} numberOfChambers","");
      Mu2Muon_numberOfChambers=HConfig.GetTH1D(Name+"_Mu2Muon_numberOfChambers","Mu2Muon_numberOfChambers",10,-0.5,9.5,"#mu_{2} numberOfChambers","");
      Mu3Muon_numberOfChambers=HConfig.GetTH1D(Name+"_Mu3Muon_numberOfChambers","Mu3Muon_numberOfChambers",10,-0.5,9.5,"#mu_{3} numberOfChambers","");



      MinMuon_numberOfMatches=HConfig.GetTH1D(Name+"_MinMuon_numberOfMatches","MinMuon_numberOfMatches",10,-0.5,9.5,"min numberOfMatches","");
      MaxMuon_numberOfMatches=HConfig.GetTH1D(Name+"_MaxMuon_numberOfMatches","MaxMuon_numberOfMatches",10,-0.5,9.5,"max numberOfMatches","");
      Mu1Muon_numberOfMatches=HConfig.GetTH1D(Name+"_Mu1Muon_numberOfMatches","Mu1Muon_numberOfMatches",10,-0.5,9.5,"#mu_{1} numberOfMatches","");
      Mu2Muon_numberOfMatches=HConfig.GetTH1D(Name+"_Mu2Muon_numberOfMatches","Mu2Muon_numberOfMatches",10,-0.5,9.5,"#mu_{2} numberOfMatches","");
      Mu3Muon_numberOfMatches=HConfig.GetTH1D(Name+"_Mu3Muon_numberOfMatches","Mu3Muon_numberOfMatches",10,-0.5,9.5,"#mu_{3} numberOfMatches","");
	
      MinMuon_hitPattern_numberOfValidMuonHits=HConfig.GetTH1D(Name+"_MinMuon_hitPattern_numberOfValidMuonHits","MinMuon_hitPattern_numberOfValidMuonHits",10,-0.5,9.5,"min numberOfValidMuonHits","");
      MaxMuon_hitPattern_numberOfValidMuonHits=HConfig.GetTH1D(Name+"_MaxMuon_hitPattern_numberOfValidMuonHits","MaxMuon_hitPattern_numberOfValidMuonHits",10,-0.5,9.5,"max numberOfValidMuonHits","");
      Mu1Muon_hitPattern_numberOfValidMuonHits=HConfig.GetTH1D(Name+"_Mu1Muon_hitPattern_numberOfValidMuonHits","Mu1Muon_hitPattern_numberOfValidMuonHits",10,-0.5,9.5,"#mu_{1} numberOfValidMuonHits","");
      Mu2Muon_hitPattern_numberOfValidMuonHits=HConfig.GetTH1D(Name+"_Mu2Muon_hitPattern_numberOfValidMuonHits","Mu2Muon_hitPattern_numberOfValidMuonHits",10,-0.5,9.5,"#mu_{2} numberOfValidMuonHits","");
      Mu3Muon_hitPattern_numberOfValidMuonHits=HConfig.GetTH1D(Name+"_Mu3Muon_hitPattern_numberOfValidMuonHits","Mu3Muon_hitPattern_numberOfValidMuonHits",10,-0.5,9.5,"#mu_{3} numberOfValidMuonHits","");



      MaxVertexPairQuality=HConfig.GetTH1D(Name+"_MaxVertexPairQuality","MaxVertexPairQuality",30,0,10,"max vertex pair quality","Events");
      MinVertexPairQuality=HConfig.GetTH1D(Name+"_MinVertexPairQuality","MinVertexPairQuality",20,0,2,"minvertex pair quality","Events");


      deltaMuZ12 = HConfig.GetTH1D(Name+"_deltaMuZ12","deltaMuZ12",30,0,0.6,"#Delta z (#mu_{1}-#mu_{2}), cm","");
      deltaMuZ13 = HConfig.GetTH1D(Name+"_deltaMuZ13","deltaMuZ13",30,0,0.6,"#Delta z (#mu_{1}-#mu_{3}), cm","");
      deltaMuZ23 = HConfig.GetTH1D(Name+"_deltaMuZ23","deltaMuZ23",30,0,0.6,"#Delta z (#mu_{2}-#mu_{3}), cm","");
      MaxdeltaMuZ = HConfig.GetTH1D(Name+"_MaxdeltaMuZ","MaxdeltaMuZ",30,0,0.6,"Max #Delta z (#mu-#mu), cm","");
      MindeltaMuZ = HConfig.GetTH1D(Name+"_MindeltaMuZ","MindeltaMuZ",30,0,0.2,"Min #Delta z (#mu-#mu), cm","");
      MaxMuonsDca=HConfig.GetTH1D(Name+"_MaxMuonsDca","MaxMuonsDca",50,0,0.10,"Max distance between muons","");
      MinMuonsDca=HConfig.GetTH1D(Name+"_MinMuonsDca","MinMuonsDca",50,0,0.04,"Min distance between muons","");



      NSV=HConfig.GetTH1D(Name+"_NSV","NSV",8,-0.5,7.5,"N vertices in the tau cone","");

      SVDeltaR=HConfig.GetTH1D(Name+"_SVDeltaR","SVDeltaR",50,0,0.3,"#Delta R (#vec{#tau}  - #vec{Vertex-PV})","");
      SVDistance=HConfig.GetTH1D(Name+"_SVDistance","SVDistance",100,0,0.5,"Distance(SV  - Vertex),cm",""); 
      SV_Mass=HConfig.GetTH1D(Name+"_SV_Mass","SV_Mass",50,0.2,2.5,"VertexMass, GeV","");
      NtracksClose=HConfig.GetTH1D(Name+"_NtracksClose","NtracksClose",8,-0.5,7.5,"Number of tracks close to SV (0.03cm)","");
      NTracksCloseToPV=HConfig.GetTH1D(Name+"_NTracksCloseToPV","NTracksCloseToPV",8,-0.5,7.5,"Number of tracks close to PV (0.05 cm)","");
      NTracksCloseToPVTauDR=HConfig.GetTH1D(Name+"_NTracksCloseToPVTauDR","NTracksCloseToPVTauDR",8,-0.5,7.5,"Number of tracks close to PV (0.05cm #Delta R (#vec{#tau}) < 0.8","");



      Mu1TrackMass=HConfig.GetTH1D(Name+"_Mu1TrackMass","Mu1TrackMass",60,0.2,4.5,"M_{#mu1-track}, GeV","");
      Mu2TrackMass=HConfig.GetTH1D(Name+"_Mu2TrackMass","Mu2TrackMass",60,0.2,4.5,"M_{#mu2-track}, GeV","");
      Mu3TrackMass=HConfig.GetTH1D(Name+"_Mu3TrackMass","Mu3TrackMass",60,0.2,4.5,"M_{#mu3-track}, GeV","");



      Iso02=HConfig.GetTH1D(Name+"_Iso02","Iso02",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 0.2 ");
      Iso04=HConfig.GetTH1D(Name+"_Iso04","Iso04",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 0.4 ");
      Iso06=HConfig.GetTH1D(Name+"_Iso06","Iso06",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 0.6 ");
      Iso08=HConfig.GetTH1D(Name+"_Iso08","Iso08",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 0.8 ");
      Iso1=HConfig.GetTH1D(Name+"_Iso1","Iso1",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 1.0");
      Iso12=HConfig.GetTH1D(Name+"_Iso12","Iso12",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 1.2");
      Iso14=HConfig.GetTH1D(Name+"_Iso14","Iso14",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 1.4");
      Iso16=HConfig.GetTH1D(Name+"_Iso16","Iso16",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 1.6 ");
      Iso18=HConfig.GetTH1D(Name+"_Iso18","Iso18",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 1.8 ");
      Iso2=HConfig.GetTH1D(Name+"_Iso2","Iso2",50,0,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T}(tracks))","#Delta R < 2");




      Iso02Mu1=HConfig.GetTH1D(Name+"_Iso02Mu1","Iso02Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 0.2 ");
      Iso04Mu1=HConfig.GetTH1D(Name+"_Iso04Mu1","Iso04Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 0.4 ");
      Iso06Mu1=HConfig.GetTH1D(Name+"_Iso06Mu1","Iso06Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 0.6 ");
      Iso08Mu1=HConfig.GetTH1D(Name+"_Iso08Mu1","Iso08Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 0.8 ");
      Iso1Mu1=HConfig.GetTH1D(Name+"_Iso1Mu1","Iso1Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 1.0");
      Iso12Mu1=HConfig.GetTH1D(Name+"_Iso12Mu1","Iso12Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 1.2");
      Iso14Mu1=HConfig.GetTH1D(Name+"_Iso14Mu1","Iso14Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 1.4");
      Iso16Mu1=HConfig.GetTH1D(Name+"_Iso16Mu1","Iso16Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 1.6 ");
      Iso18Mu1=HConfig.GetTH1D(Name+"_Iso18Mu1","Iso18Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 1.8 ");
      Iso2Mu1=HConfig.GetTH1D(Name+"_Iso2Mu1","Iso2Mu1",50,0,1.1,"I_{#mu1}= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 2");


      Iso02Mu2=HConfig.GetTH1D(Name+"_Iso02Mu2","Iso02Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 0.2 ");
      Iso04Mu2=HConfig.GetTH1D(Name+"_Iso04Mu2","Iso04Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 0.4 ");
      Iso06Mu2=HConfig.GetTH1D(Name+"_Iso06Mu2","Iso06Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 0.6 ");
      Iso08Mu2=HConfig.GetTH1D(Name+"_Iso08Mu2","Iso08Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 0.8 ");
      Iso1Mu2=HConfig.GetTH1D(Name+"_Iso1Mu2","Iso1Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 1.0");
      Iso12Mu2=HConfig.GetTH1D(Name+"_Iso12Mu2","Iso12Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 1.2");
      Iso14Mu2=HConfig.GetTH1D(Name+"_Iso14Mu2","Iso14Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 1.4");
      Iso16Mu2=HConfig.GetTH1D(Name+"_Iso16Mu2","Iso16Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 1.6 ");
      Iso18Mu2=HConfig.GetTH1D(Name+"_Iso18Mu2","Iso18Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 1.8 ");
      Iso2Mu2=HConfig.GetTH1D(Name+"_Iso2Mu2","Iso2Mu2",50,0,1.1,"I_{#mu2}= p_{T}(#mu2)/(p_{T}(#mu2) + #sum p_{T}(tracks))","#Delta R < 2");

      Iso02Mu3=HConfig.GetTH1D(Name+"_Iso02Mu3","Iso02Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 0.2 ");
      Iso04Mu3=HConfig.GetTH1D(Name+"_Iso04Mu3","Iso04Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 0.4 ");
      Iso06Mu3=HConfig.GetTH1D(Name+"_Iso06Mu3","Iso06Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 0.6 ");
      Iso08Mu3=HConfig.GetTH1D(Name+"_Iso08Mu3","Iso08Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 0.8 ");
      Iso1Mu3=HConfig.GetTH1D(Name+"_Iso1Mu3","Iso1Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 1.0");
      Iso12Mu3=HConfig.GetTH1D(Name+"_Iso12Mu3","Iso12Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 1.2");
      Iso14Mu3=HConfig.GetTH1D(Name+"_Iso14Mu3","Iso14Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 1.4");
      Iso16Mu3=HConfig.GetTH1D(Name+"_Iso16Mu3","Iso16Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 1.6 ");
      Iso18Mu3=HConfig.GetTH1D(Name+"_Iso18Mu3","Iso18Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 1.8 ");
      Iso2Mu3=HConfig.GetTH1D(Name+"_Iso2Mu3","Iso2Mu3",50,0,1.1,"I_{#mu3}= p_{T}(#mu3)/(p_{T}(#mu3) + #sum p_{T}(tracks))","#Delta R < 2");


      Iso08MuMax  =HConfig.GetTH1D(Name+"_Iso08MuMax","Iso08MuMax",50,0,1.1,"Max I_{#mu}= p_{T}(#mu)/(p_{T}(#mu) + #sum p_{T}(tracks))","#Delta R < 0.8 ");
      Iso08MuMin  =HConfig.GetTH1D(Name+"_Iso08MuMin","Iso08MuMin",50,0,1.1,"Min I_{#mu}= p_{T}(#mu)/(p_{T}(#mu) + #sum p_{T}(tracks))","#Delta R < 0.8 ");
      
      MindcaTrackSV=HConfig.GetTH1D(Name+"_MindcaTrackSV","MindcaTrackSV",50,0,0.1,"Min distance of track to SV","");
      dcaTrackPV=HConfig.GetTH1D(Name+"_dcaTrackPV","dcaTrackPV",50,0,0.1,"distance of closest approach to PV","");

      NBJet4pi=HConfig.GetTH1D(Name+"_NBJet4pi","NBJet4pi",10,-0.5,9.5,"N Jets","");


      BTagCSV =HConfig.GetTH1D(Name+"_BTagCSV","BTagCSV",50,0.,1.,"BTagCSV ","");
      BTagMVA =HConfig.GetTH1D(Name+"_BTagMVA","BTagMVA",50,-1.,1.,"BTagMVA ","");
      BTagCVSB=HConfig.GetTH1D(Name+"_BTagCVSB","BTagCVSB",50,-1.,1.,"BTagCVSB ","");

      NBJet4piSH=HConfig.GetTH1D(Name+"_NBJet4piSH","NBJet4piSH",10,-0.5,9.5,"N Jets SH","");
      NBJet4piOH=HConfig.GetTH1D(Name+"_NBJet4piOH","NBJet4piOH",10,-0.5,9.5,"N Jets OH","");

      BTagCSVSH =HConfig.GetTH1D(Name+"_BTagCSVSH","BTagCSVSH",50,0.,1.,"BTagCSV SH","");
      BTagMVASH =HConfig.GetTH1D(Name+"_BTagMVASH","BTagMVASH",50,-1.,1.,"BTagMVA SH","");
      BTagCVSBSH=HConfig.GetTH1D(Name+"_BTagCVSBSH","BTagCVSBSH",50,-1.,1.,"BTagCVSB SH","");

      BTagCSVOH =HConfig.GetTH1D(Name+"_BTagCSVOH","BTagCSVOH",50,0.,1.,"BTagCSV OH ","");
      BTagMVAOH =HConfig.GetTH1D(Name+"_BTagMVAOH","BTagMVAOH",50,-1.,1.,"BTagMVA OH","");
      BTagCVSBOH=HConfig.GetTH1D(Name+"_BTagCVSBOH","BTagCVSBOH",50,-1.,1.,"BTagCVSB OH","");




      BTagCSVSHMatchedToTau =HConfig.GetTH1D(Name+"_BTagCSVSHMatchedToTau","BTagCSVSHMatchedToTau",50,0.,1.,"BTagCSV SH","");
      BTagMVASHMatchedToTau =HConfig.GetTH1D(Name+"_BTagMVASHMatchedToTau","BTagMVASHMatchedToTau",50,-1.,1.,"BTagMVA SH","");
      BTagCVSBSHMatchedToTau=HConfig.GetTH1D(Name+"_BTagCVSBSHMatchedToTau","BTagCVSBSHMatchedToTau",50,-1.,1.,"BTagCVSB SH","");

      BTagCSVOHMatchedToTau =HConfig.GetTH1D(Name+"_BTagCSVOHMatchedToTau","BTagCSVOHMatchedToTau",50,0.,1.,"BTagCSV OH ","");
      BTagMVAOHMatchedToTau =HConfig.GetTH1D(Name+"_BTagMVAOHMatchedToTau","BTagMVAOHMatchedToTau",50,-1.,1.,"BTagMVA OH","");
      BTagCVSBOHMatchedToTau=HConfig.GetTH1D(Name+"_BTagCVSBOHMatchedToTau","BTagCVSBOHMatchedToTau",50,-1.,1.,"BTagCVSB OH","");



      BTagCSVSHVsNJets =HConfig.GetTH2D(Name+"_BTagCSVSHVsNJets","BTagCSVSHVsNJets",50,0.,1.,10,-0.5,9.5,"BTagCSV SH","N");
      BTagMVASHVsNJets =HConfig.GetTH2D(Name+"_BTagMVASHVsNJets","BTagMVASHVsNJets",50,-1.,1.,10,-0.5,9.5,"BTagMVA SH","N");
      BTagCVSBSHVsNJets=HConfig.GetTH2D(Name+"_BTagCVSBSHVsNJets","BTagCVSBSHVsNJets",50,-1.,1.,10,-0.5,9.5,"BTagCVSB SH","N");

      BTagCSVOHVsNJets =HConfig.GetTH2D(Name+"_BTagCSVOHVsNJets","BTagCSVOHVsNJets",50,0.,1.,10,-0.5,9.5,"BTagCSV OH ","N");
      BTagMVAOHVsNJets =HConfig.GetTH2D(Name+"_BTagMVAOHVsNJets","BTagMVAOHVsNJets",50,-1.,1.,10,-0.5,9.5,"BTagMVA OH","N");
      BTagCVSBOHVsNJets=HConfig.GetTH2D(Name+"_BTagCVSBOHVsNJets","BTagCVSBOHVsNJets",50,-1.,1.,10,-0.5,9.5,"BTagCVSB OH","N");

      Muon1ImpactAngle =HConfig.GetTH1D(Name+"_Muon1ImpactAngle","Muon1ImpactAngle",30,-1,1,"#mu_{1} impact angle","");
      Muon2ImpactAngle =HConfig.GetTH1D(Name+"_Muon2ImpactAngle","Muon2ImpactAngle",30,-1,1,"#mu_{2} impact angle","");
      Muon3ImpactAngle =HConfig.GetTH1D(Name+"_Muon3ImpactAngle","Muon3ImpactAngle",30,-1,1,"#mu_{3} impact angle","");
      MinMuonImpacAngle =HConfig.GetTH1D(Name+"_MinMuonImpactAngle","MinMuonImpactAngle",30,-1,1,"Min #mu impact angle","");
      MaxMuonImpacAngle =HConfig.GetTH1D(Name+"_MaxMuonImpactAngle","MaxMuonImpactAngle",30,-1,1,"Max #mu impact angle","");


      Muon1StandardSelectorPass=HConfig.GetTH1D(Name+"_Muon1StandardSelectorPass","Muon1StandardSelectorPass",23,-0.5,22.5,"#mu_{1} standard selector(pass); bin 0 - no ID","Events");
      Muon2StandardSelectorPass=HConfig.GetTH1D(Name+"_Muon2StandardSelectorPass","Muon2StandardSelectorPass",23,-0.5,22.5,"#mu_{2} standard selector(pass); bin 0 - no ID","Events");
      Muon3StandardSelectorPass=HConfig.GetTH1D(Name+"_Muon3StandardSelectorPass","Muon3StandardSelectorPass",23,-0.5,22.5,"#mu_{3} standard selector(pass); bin 0 - no ID","Events");

      Muon1StandardSelectorFail=HConfig.GetTH1D(Name+"_Muon1StandardSelectorFail","Muon1StandardSelectorFail",23,-0.5,22.5,"#mu_{1} standard selector(fail); bin 0 - no ID","Events");
      Muon2StandardSelectorFail=HConfig.GetTH1D(Name+"_Muon2StandardSelectorFail","Muon2StandardSelectorFail",23,-0.5,22.5,"#mu_{2} standard selector(fail); bin 0 - no ID","Events");
      Muon3StandardSelectorFail=HConfig.GetTH1D(Name+"_Muon3StandardSelectorFail","Muon3StandardSelectorFail",23,-0.5,22.5,"#mu_{3} standard selector(fail); bin 0 - no ID","Events");

      Muon1MVAID=HConfig.GetTH1D(Name+"_Muon1MVAID","Muon1MVAID",50,0.0,1.0,"#mu_{1} MVA","Events");
      Muon2MVAID=HConfig.GetTH1D(Name+"_Muon2MVAID","Muon2MVAID",50,0.0,1.0,"#mu_{2} MVA","Events");
      Muon3MVAID=HConfig.GetTH1D(Name+"_Muon3MVAID","Muon3MVAID",50,0.0,1.0,"#mu_{3} MVA","Events");

  
      PairMass1AllignedSorting=HConfig.GetTH1D(Name+"_PairMass1AllignedSorting","PairMass1AllignedSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 1st collimated pair), GeV","Events");
      PairMass2AllignedSorting=HConfig.GetTH1D(Name+"_PairMass2AllignedSorting","PairMass2AllignedSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 2nd collimated pair), GeV","Events");
      MuMuMassAllignedSorting=HConfig.GetTH2D(Name+"_MuMuMassAllignedSorting","MuMuMassAllignedSorting",55,0.1,2.0,50,0.1,2.0,"M_{#mu#mu} (OS-SS, 2nd collimated pair), GeV","M_{#mu#mu} (OS-SS, 1st collimated pair) GeV");



      PairMass1PTSorting=HConfig.GetTH1D(Name+"_PairMass1PTSorting","PairMass1PTSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 1st pT sorting pair), GeV","Events");
      PairMass2PTSorting=HConfig.GetTH1D(Name+"_PairMass2PTSorting","PairMass2PTSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 2nd pT sorting pair), GeV","Events");
      MuMuMassPTSorting=HConfig.GetTH2D(Name+"_MuMuMassPTSorting","MuMuMassPTSorting",55,0.1,2.0,50,0.1,2.0,"M_{#mu#mu} (OS-SS, 1st pT sorting pair), GeV","M_{#mu#mu} (OS-SS, 2nd pT sorting pair), GeV");


      PairMass1NoSorting=HConfig.GetTH1D(Name+"_PairMass1NoSorting","PairMass1NoSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 1 no sorting), GeV","Events");
      PairMass2NoSorting=HConfig.GetTH1D(Name+"_PairMass2NoSorting","PairMass2NoSorting",55,0.1,2.0,"M_{#mu#mu} (OS-SS, 2 no sorting), GeV","Events");
      MuMuMassNoSorting=HConfig.GetTH2D(Name+"_MuMuMassNoSorting","MuMuMassNoSorting",55,0.1,2.0,50,0.1,2.0,"M_{#mu#mu} (OS-SS, 1 no sorting), GeV","M_{#mu#mu} (OS-SS, 2 no sorting sorting), GeV");


      IsoPhiKKMass_Mu3=HConfig.GetTH1D(Name+"_IsoPhiKKMass_Mu3","IsoPhiKKMass_Mu3",55,0.95,1.20,"M_{KK},GeV","Events");
      IsoPhiKKMass_Mu2=HConfig.GetTH1D(Name+"_IsoPhiKKMass_Mu2","IsoPhiKKMass_Mu1",55,0.95,1.20,"M_{KK},GeV","Events");
      IsoPhiKKMass_Mu1=HConfig.GetTH1D(Name+"_IsoPhiKKMass_Mu1","IsoPhiKKMass_Mu1",55,0.95,1.20,"M_{KK},GeV","Events");


      IsoKStarMass_Mu3=HConfig.GetTH1D(Name+"_IsoKStarMass_Mu3","IsoKStarMass_Mu3",55,0.65,2.1,"M_{K#pi},GeV","Events");
      IsoKStarMass_Mu2=HConfig.GetTH1D(Name+"_IsoKStarMass_Mu2","IsoKStarMass_Mu1",55,0.65,2.1,"M_{K#pi},GeV","Events");
      IsoKStarMass_Mu1=HConfig.GetTH1D(Name+"_IsoKStarMass_Mu1","IsoKStarMass_Mu1",55,0.65,2.1,"M_{K#pi},GeV","Events");



      IsoMuMuMass_Mu3=HConfig.GetTH1D(Name+"_IsoMuMuMass_Mu3","IsoMuMuMass_Mu3",55,0.25,1.20,"M_{#mu#mu},GeV","Events");
      IsoMuMuMass_Mu2=HConfig.GetTH1D(Name+"_IsoMuMuMass_Mu2","IsoMuMuMass_Mu1",55,0.25,1.20,"M_{#mu#mu},GeV","Events");
      IsoMuMuMass_Mu1=HConfig.GetTH1D(Name+"_IsoMuMuMass_Mu1","IsoMuMuMass_Mu1",55,0.25,1.20,"M_{#mu#mu},GeV","Events");




      IsoPhiKKMass_Mu3_wideRange=HConfig.GetTH1D(Name+"_IsoPhiKKMass_Mu3_wideRange","IsoPhiKKMass_Mu3_wideRange",75,0.95,3.2,"M_{KK},GeV","Events");
      IsoPhiKKMass_Mu2_wideRange=HConfig.GetTH1D(Name+"_IsoPhiKKMass_Mu2_wideRange","IsoPhiKKMass_Mu1_wideRange",75,0.95,3.2,"M_{KK},GeV","Events");
      IsoPhiKKMass_Mu1_wideRange=HConfig.GetTH1D(Name+"_IsoPhiKKMass_Mu1_wideRange","IsoPhiKKMass_Mu1_wideRange",75,0.95,3.2,"M_{KK},GeV","Events");


      IsoPhiKKMass_Mu3_midRange=HConfig.GetTH1D(Name+"_IsoPhiKKMass_Mu3_midRange","IsoPhiKKMass_Mu3_midRange",75,1.2,1.8,"M_{KK},GeV","Events");
      IsoPhiKKMass_Mu2_midRange=HConfig.GetTH1D(Name+"_IsoPhiKKMass_Mu2_midRange","IsoPhiKKMass_Mu1_midRange",75,1.2,1.8,"M_{KK},GeV","Events");
      IsoPhiKKMass_Mu1_midRange=HConfig.GetTH1D(Name+"_IsoPhiKKMass_Mu1_midRange","IsoPhiKKMass_Mu1_midRange",75,1.2,1.8,"M_{KK},GeV","Events");



      IsoKStarMass_Mu3_wideRange=HConfig.GetTH1D(Name+"_IsoKStarMass_Mu3_wideRange","IsoKStarMass_Mu3_wideRange",75,0.65,3.2,"M_{K#pi},GeV","Events");
      IsoKStarMass_Mu2_wideRange=HConfig.GetTH1D(Name+"_IsoKStarMass_Mu2_wideRange","IsoKStarMass_Mu1_wideRange",75,0.65,3.2,"M_{K#pi},GeV","Events");
      IsoKStarMass_Mu1_wideRange=HConfig.GetTH1D(Name+"_IsoKStarMass_Mu1_wideRange","IsoKStarMass_Mu1_wideRange",75,0.65,3.2,"M_{K#pi},GeV","Events");



      IsoMuMuMass_Mu3_wideRange=HConfig.GetTH1D(Name+"_IsoMuMuMass_Mu3_wideRange","IsoMuMuMass_Mu3_wideRange",75,0.25,3.2,"M_{#mu#mu},GeV","Events");
      IsoMuMuMass_Mu2_wideRange=HConfig.GetTH1D(Name+"_IsoMuMuMass_Mu2_wideRange","IsoMuMuMass_Mu1_wideRange",75,0.25,3.2,"M_{#mu#mu},GeV","Events");
      IsoMuMuMass_Mu1_wideRange=HConfig.GetTH1D(Name+"_IsoMuMuMass_Mu1_wideRange","IsoMuMuMass_Mu1_wideRange",75,0.25,3.2,"M_{#mu#mu},GeV","Events");


      Separation_BvsD=HConfig.GetTH1D(Name+"_Separation_BvsD","Separation_BvsD",70,-1.0,1.0,"B vs D BDTG output","Events");



      IsolationCombinatorialMass_pipi=HConfig.GetTH1D(Name+"_IsolationCombinatorialMass_pipi","IsolationCombinatorialMass_pipi",55,0.27,2.0,"M_{#pi#pi},GeV","Events");
      VertexQualitySeparator=HConfig.GetTH1D(Name+"_VertexQualitySeparator","VertexQualitySeparator",2,-0.5,1.5,"0 - #chi^{2}_{3 #mu} > #chi^{2}_{2 #mu-iso track} ; 1 -","Events");


      OSSS1Angle_TRF=HConfig.GetTH1D(Name+"_OSSS1Angle_TRF","OSSS1Angle_TRF",50,-1.1,1.1,"cos#alpha_{1} (OS-SS1), 3#mu RF","");
      OSSS2Angle_TRF=HConfig.GetTH1D(Name+"_OSSS2Angle_TRF","OSSS2Angle_TRF",50,-1.1,1.1,"cos#alpha_{2} (OS-SS2), 3#mu RF","");


      OSSS1Angle_RRF=HConfig.GetTH1D(Name+"_OSSS1Angle_RRF","OSSS1Angle_RRF",50,-1.1,1.1,"cos#theta (n_{#rho} - n_{#mu})","");
      OSSS2Angle_RRF=HConfig.GetTH1D(Name+"_OSSS2Angle_RRF","OSSS2Angle_RRF",50,-1.1,1.1,"cos#theta (n_{#rho} - n_{#mu})","");

      cTheta_TRF_SSSS=HConfig.GetTH1D(Name+"_cTheta_TRF_SSSS","cTheta_TRF_SSSS",50,-1.1,1.1,"cos#theta (SS1-SS2), 3#mu RF","");
      cTheta_TRF_OSSS=HConfig.GetTH1D(Name+"_cTheta_TRF_OSSS","cTheta_TRF_OSSS",50,-1.1,1.1,"cos#theta (OS -SS1), 3#mu RF","");
      cTheta_MuonOS_TauPol_TRF=HConfig.GetTH1D(Name+"_cTheta_MuonOS_TauPol_TRF","cTheta_MuonOS_TauPol_TRF",50,-1.1,1.1,"cos#Theta (OS - #tau), 3#mu RF","");

      pTMu1OverMass_TRF=HConfig.GetTH1D(Name+"_pTMu1OverMass_TRF","pTMu1OverMass_TRF",50,0,0.6,"P_{#mu_{OS}}/M_{3#mu}  ","");



      Selection::ConfigureHistograms(); //do not remove
      HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
      
}

void  MakeMVATree::Store_ExtraDist(){ 


  Extradist1d.push_back(&pTMu1OverMass_TRF);
  Extradist1d.push_back(&OSSS1Angle_TRF);
  Extradist1d.push_back(&OSSS2Angle_TRF);

  Extradist1d.push_back(&OSSS1Angle_RRF);
  Extradist1d.push_back(&OSSS2Angle_RRF);


  Extradist1d.push_back(&cTheta_MuonOS_TauPol_TRF);


  Extradist1d.push_back(&cTheta_TRF_SSSS);
  Extradist1d.push_back(&cTheta_TRF_OSSS);


  Extradist1d.push_back(&nPVx);

  Extradist1d.push_back(&Muon1Pt);
  Extradist1d.push_back(&Muon2Pt);
  Extradist1d.push_back(&Muon3Pt);

  Extradist1d.push_back(&Muon1P);
  Extradist1d.push_back(&Muon2P);
  Extradist1d.push_back(&Muon3P);

  Extradist1d.push_back(&TauEta);
  Extradist1d.push_back(&TauPt);
  Extradist1d.push_back(&TauP);


  Extradist1d.push_back(&SVPVTauDirAngle);
  Extradist1d.push_back(&FLSignificance);
  Extradist1d.push_back(&FL);
  Extradist1d.push_back(&VertexChi2KF);
  Extradist1d.push_back(&Vertex2muTrkKF);
  Extradist1d.push_back(&VertexQualitySeparator);

  Extradist1d.push_back(&Vertex2muTrkKFChi2);
  Extradist1d.push_back(&Vertex2muTrkKFToSignalVertexChi2);
  Extradist1d.push_back(&Vertex2muTrkKFToSignalVertexDistance);





  Extradist1d.push_back(&Dist2muTrkKF3Mu);
  //  Extradist1d.push_back(&MuonglbkinkSum);
  Extradist1d.push_back(&Muon_segmentCompatibility_mu1);
  Extradist1d.push_back(&Muon_segmentCompatibility_mu2);
  Extradist1d.push_back(&Muon_segmentCompatibility_mu3);
  Extradist1d.push_back(&Muon_segmentCompatibility_min);
  Extradist1d.push_back(&Muon_segmentCompatibility_max);
  Extradist1d.push_back(&Muon_ECALCompatibility_mu1);
  Extradist1d.push_back(&Muon_ECALCompatibility_mu2);
  Extradist1d.push_back(&Muon_ECALCompatibility_mu3);
  Extradist1d.push_back(&Muon_ECALCompatibility_min);
  Extradist1d.push_back(&Muon_ECALCompatibility_max);
  Extradist1d.push_back(&Isolation_RelPt);
  Extradist1d.push_back(&Isolation_NTracks);

  Extradist1d.push_back(&Isolation_maxdxy);

  Extradist1d.push_back(&EventMassResolution_PtEtaPhi);
  Extradist2d.push_back(&EMR_tau_eta);

  Extradist1d.push_back(&MaxMuonsDca);
  Extradist1d.push_back(&MinMuonsDca);


  Extradist1d.push_back(&deltaMuZ12);
  Extradist1d.push_back(&deltaMuZ13);
  Extradist1d.push_back(&deltaMuZ23);
  Extradist1d.push_back(&MaxdeltaMuZ);
  Extradist1d.push_back(&MindeltaMuZ);

  Extradist1d.push_back(&NSV);
  Extradist1d.push_back(&SVDeltaR);
  Extradist1d.push_back(&SVDistance);
  Extradist1d.push_back(&SV_Mass);
  Extradist1d.push_back(&NtracksClose);
  Extradist1d.push_back(&NTracksCloseToPV);
  Extradist1d.push_back(&NTracksCloseToPVTauDR);




  Extradist1d.push_back(&VertexMu1D0SigPVReco);
  Extradist1d.push_back(&VertexMu2D0SigPVReco);
  Extradist1d.push_back(&VertexMu3D0SigPVReco);
  Extradist1d.push_back(&MaxD0SigPV);
  Extradist1d.push_back(&MinD0SigPV);


  Extradist1d.push_back(&VertexMu1D0SigBSReco);
  Extradist1d.push_back(&VertexMu2D0SigBSReco);
  Extradist1d.push_back(&VertexMu3D0SigBSReco);
  Extradist1d.push_back(&MaxD0SigBS);
  Extradist1d.push_back(&MinD0SigBS);


  Extradist1d.push_back(&VertexMu1D0SigSVReco);
  Extradist1d.push_back(&VertexMu2D0SigSVReco);
  Extradist1d.push_back(&VertexMu3D0SigSVReco);
  Extradist1d.push_back(&MaxD0SigSV);
  Extradist1d.push_back(&MinD0SigSV);

  Extradist1d.push_back(&VertexMu1D0DistSVReco);
  Extradist1d.push_back(&VertexMu2D0DistSVReco);
  Extradist1d.push_back(&VertexMu3D0DistSVReco);

  Extradist1d.push_back(&MaxD0DistSV);
  Extradist1d.push_back(&MinD0DistSV);

  Extradist1d.push_back(&VertexMu1DZDistSVReco);
  Extradist1d.push_back(&VertexMu2DZDistSVReco);
  Extradist1d.push_back(&VertexMu3DZDistSVReco);

  Extradist1d.push_back(&MaxDZDistSV);
  Extradist1d.push_back(&MinDZDistSV);


  Extradist1d.push_back(&VertexMu1DistanceToSV);
  Extradist1d.push_back(&VertexMu2DistanceToSV);
  Extradist1d.push_back(&VertexMu3DistanceToSV);

  Extradist1d.push_back(&MaxMuDistanceToSV);
  Extradist1d.push_back(&MinMuDistanceToSV);



  Extradist1d.push_back(&MinMuon_chi2LocalPosition);
  Extradist1d.push_back(&MaxMuon_chi2LocalPosition);

  Extradist1d.push_back(&Muon1_chi2LocalPosition);
  Extradist1d.push_back(&Muon2_chi2LocalPosition);
  Extradist1d.push_back(&Muon3_chi2LocalPosition);
  
	
  Extradist1d.push_back(&MinMuon_chi2LocalMomentum);
  Extradist1d.push_back(&MaxMuon_chi2LocalMomentum);

  Extradist1d.push_back(&Muon1_chi2LocalMomentum);
  Extradist1d.push_back(&Muon2_chi2LocalMomentum);
  Extradist1d.push_back(&Muon3_chi2LocalMomentum);
  
  Extradist1d.push_back(&MintrkKink);
  Extradist1d.push_back(&MaxtrkKink);
  Extradist1d.push_back(&MinglbKink);
  Extradist1d.push_back(&MaxglbKink);
  
  Extradist1d.push_back(&MuonglbkinkSum);
  
  Extradist1d.push_back(&MaxVertexPairQuality);
  Extradist1d.push_back(&MinVertexPairQuality);
	

  Extradist1d.push_back(&Mu1TrackMass);
  Extradist1d.push_back(&Mu2TrackMass);
  Extradist1d.push_back(&Mu3TrackMass);

  Extradist1d.push_back(&Iso02);
  Extradist1d.push_back(&Iso04);
  Extradist1d.push_back(&Iso06);
  Extradist1d.push_back(&Iso08);
  Extradist1d.push_back(&Iso1);
  Extradist1d.push_back(&Iso12);
  Extradist1d.push_back(&Iso14);
  Extradist1d.push_back(&Iso16);
  Extradist1d.push_back(&Iso18);
  Extradist1d.push_back(&Iso2);

  Extradist1d.push_back(&Iso02Mu1);
  Extradist1d.push_back(&Iso04Mu1);
  Extradist1d.push_back(&Iso06Mu1);
  Extradist1d.push_back(&Iso08Mu1);
  Extradist1d.push_back(&Iso1Mu1);
  Extradist1d.push_back(&Iso12Mu1);
  Extradist1d.push_back(&Iso14Mu1);
  Extradist1d.push_back(&Iso16Mu1);
  Extradist1d.push_back(&Iso18Mu1);
  Extradist1d.push_back(&Iso2Mu1);

  Extradist1d.push_back(&Iso02Mu2);
  Extradist1d.push_back(&Iso04Mu2);
  Extradist1d.push_back(&Iso06Mu2);
  Extradist1d.push_back(&Iso08Mu2);
  Extradist1d.push_back(&Iso1Mu2);
  Extradist1d.push_back(&Iso12Mu2);
  Extradist1d.push_back(&Iso14Mu2);
  Extradist1d.push_back(&Iso16Mu2);
  Extradist1d.push_back(&Iso18Mu2);
  Extradist1d.push_back(&Iso2Mu2);
  
  Extradist1d.push_back(&Iso02Mu3);
  Extradist1d.push_back(&Iso04Mu3);
  Extradist1d.push_back(&Iso06Mu3);
  Extradist1d.push_back(&Iso08Mu3);
  Extradist1d.push_back(&Iso1Mu3);
  Extradist1d.push_back(&Iso12Mu3);
  Extradist1d.push_back(&Iso14Mu3);
  Extradist1d.push_back(&Iso16Mu3);
  Extradist1d.push_back(&Iso18Mu3);
  Extradist1d.push_back(&Iso2Mu3);
  


  Extradist1d.push_back(&IsoPhiKKMass_Mu3);
  Extradist1d.push_back(&IsoPhiKKMass_Mu2);
  Extradist1d.push_back(&IsoPhiKKMass_Mu1);


  Extradist1d.push_back(&IsoKStarMass_Mu3);
  Extradist1d.push_back(&IsoKStarMass_Mu2);
  Extradist1d.push_back(&IsoKStarMass_Mu1);


  Extradist1d.push_back(&IsoMuMuMass_Mu3);
  Extradist1d.push_back(&IsoMuMuMass_Mu2);
  Extradist1d.push_back(&IsoMuMuMass_Mu1);



  Extradist1d.push_back(&IsoPhiKKMass_Mu3_wideRange);
  Extradist1d.push_back(&IsoPhiKKMass_Mu2_wideRange);
  Extradist1d.push_back(&IsoPhiKKMass_Mu1_wideRange);

  Extradist1d.push_back(&IsoPhiKKMass_Mu3_midRange);
  Extradist1d.push_back(&IsoPhiKKMass_Mu2_midRange);
  Extradist1d.push_back(&IsoPhiKKMass_Mu1_midRange);


  Extradist1d.push_back(&IsoKStarMass_Mu3_wideRange);
  Extradist1d.push_back(&IsoKStarMass_Mu2_wideRange);
  Extradist1d.push_back(&IsoKStarMass_Mu1_wideRange);


  Extradist1d.push_back(&IsoMuMuMass_Mu3_wideRange);
  Extradist1d.push_back(&IsoMuMuMass_Mu2_wideRange);
  Extradist1d.push_back(&IsoMuMuMass_Mu1_wideRange);







  Extradist1d.push_back(&IsolationCombinatorialMass_pipi);

  Extradist1d.push_back(&Iso08MuMax);
  Extradist1d.push_back(&Iso08MuMin);


  Extradist1d.push_back(&MindcaTrackSV);
  Extradist1d.push_back(&dcaTrackPV);

  Extradist1d.push_back(&MinMatchedStations);
  Extradist1d.push_back(&MaxMatchedStations);
  Extradist1d.push_back(&Mu1MatchedStations);
  Extradist1d.push_back(&Mu2MatchedStations);
  Extradist1d.push_back(&Mu3MatchedStations);

  
  Extradist1d.push_back(&MinMuon_numberOfChambers);
  Extradist1d.push_back(&MaxMuon_numberOfChambers);
  Extradist1d.push_back(&Mu1Muon_numberOfChambers);
  Extradist1d.push_back(&Mu2Muon_numberOfChambers);
  Extradist1d.push_back(&Mu3Muon_numberOfChambers);



  Extradist1d.push_back(&MinMuon_numberOfMatches);
  Extradist1d.push_back(&MaxMuon_numberOfMatches);
  Extradist1d.push_back(&Mu1Muon_numberOfMatches);
  Extradist1d.push_back(&Mu2Muon_numberOfMatches);
  Extradist1d.push_back(&Mu3Muon_numberOfMatches);
  
  Extradist1d.push_back(&MinMuon_hitPattern_numberOfValidMuonHits);
  Extradist1d.push_back(&MaxMuon_hitPattern_numberOfValidMuonHits);
  Extradist1d.push_back(&Mu1Muon_hitPattern_numberOfValidMuonHits);
  Extradist1d.push_back(&Mu2Muon_hitPattern_numberOfValidMuonHits);
  Extradist1d.push_back(&Mu3Muon_hitPattern_numberOfValidMuonHits);


  //  Extradist1d.push_back(&NBJet4pi);
  //  Extradist1d.push_back(&BTagCSV);
  //  Extradist1d.push_back(&BTagMVA);
  //  Extradist1d.push_back(&BTagCVSB);

  Extradist1d.push_back(&NBJet4piSH);
  Extradist1d.push_back(&NBJet4piOH);

  Extradist1d.push_back(&BTagCSVSH);
  Extradist1d.push_back(&BTagMVASH);
  Extradist1d.push_back(&BTagCVSBSH);

  Extradist1d.push_back(&BTagCSVOH);
  Extradist1d.push_back(&BTagMVAOH);
  Extradist1d.push_back(&BTagCVSBOH);

  Extradist1d.push_back(&BTagCSVSHMatchedToTau);
  Extradist1d.push_back(&BTagMVASHMatchedToTau);
  Extradist1d.push_back(&BTagCVSBSHMatchedToTau);

  Extradist1d.push_back(&BTagCSVOHMatchedToTau);
  Extradist1d.push_back(&BTagMVAOHMatchedToTau);
  Extradist1d.push_back(&BTagCVSBOHMatchedToTau);


  Extradist2d.push_back(&BTagCSVSHVsNJets);
  Extradist2d.push_back(&BTagMVASHVsNJets);
  Extradist2d.push_back(&BTagCVSBSHVsNJets);
  
  Extradist2d.push_back(&BTagCSVOHVsNJets);
  Extradist2d.push_back(&BTagMVAOHVsNJets);
  Extradist2d.push_back(&BTagCVSBOHVsNJets);

  
  Extradist1d.push_back(&Muon1ImpactAngle);
  Extradist1d.push_back(&Muon2ImpactAngle);
  Extradist1d.push_back(&Muon3ImpactAngle);
  Extradist1d.push_back(&MinMuonImpacAngle);
  Extradist1d.push_back(&MaxMuonImpacAngle);


  Extradist1d.push_back(&Muon1StandardSelectorPass);
  Extradist1d.push_back(&Muon2StandardSelectorPass);
  Extradist1d.push_back(&Muon3StandardSelectorPass);

  Extradist1d.push_back(&Muon1StandardSelectorFail);
  Extradist1d.push_back(&Muon2StandardSelectorFail);
  Extradist1d.push_back(&Muon3StandardSelectorFail);

  Extradist1d.push_back(&Muon1MVAID);
  Extradist1d.push_back(&Muon2MVAID);
  Extradist1d.push_back(&Muon3MVAID);



  Extradist1d.push_back(&PairMass1AllignedSorting);
  Extradist1d.push_back(&PairMass2AllignedSorting);
  Extradist2d.push_back(&MuMuMassAllignedSorting);


  Extradist1d.push_back(&PairMass1PTSorting);
  Extradist1d.push_back(&PairMass2PTSorting);
  Extradist2d.push_back(&MuMuMassPTSorting);


  Extradist1d.push_back(&PairMass1NoSorting);
  Extradist1d.push_back(&PairMass2NoSorting);
  Extradist2d.push_back(&MuMuMassNoSorting);


  Extradist1d.push_back(&Separation_BvsD);

}



void  MakeMVATree::doEvent(){ 
    unsigned int t;
    int id(Ntp->GetMCID());
    bool hlt_pass = false;
    if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
    //    std::cout<<"  id  "<< id << std::endl;
    bool HLTOk(false);
    bool L1Ok(false);
    bool DoubleMu0Fired(false);
    bool DoubleMu4Fired(false);
    bool DoubleMuFired(false);
    bool TripleMuFired(false);
    bool randomFailed(false);

    random_num = rndm.Uniform();

    for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
      TString HLT = Ntp->HLTName(iTrigger);
      
      //      std::cout<<"   HLT:    "<< Ntp->HLTName(iTrigger) <<"    decision    " <<  Ntp->HLTDecision(iTrigger)  << std::endl;
      if(HLT.Contains("HLT_Tau3Mu_Mu7_Mu1_TkMu1_") || HLT.Contains("HLT_DoubleMu3_Trk_Tau3mu") || HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu") )
	{
	  if( Ntp->HLTDecision(iTrigger)  ) 
	    {
	      HLTOk = true;
	    }
	}
    }

    for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
      TString L1TriggerName = Ntp->L1Name(il1);
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4") && Ntp->L1Decision(il1)) { DoubleMu0Fired = true; }
      if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9") && Ntp->L1Decision(il1)) { TripleMuFired = true; }
      if( id!=1 && random_num>0.30769 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMu4Fired = true;}
      if( id==1 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMu4Fired = true; }
      if( id!=1 && random_num<0.30769 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) {
	randomFailed = true;
      }
    }
    if (DoubleMu0Fired || DoubleMu4Fired) {DoubleMuFired = true;}
    if (!DoubleMu0Fired && !TripleMuFired && randomFailed) l1FailedRandom++;
    if (DoubleMuFired || TripleMuFired) L1Ok = true;

    if (HLTOk) value.at(HLT) = true;
    else value.at(HLT) = false;

    if (L1Ok) value.at(L1) = true;
    else value.at(L1) = false;

    //    if (HLTFired && !L1Fired && !randomFailed) cout<<"wrong hlt"<<endl;

    pass.at(HLT) = (value.at(HLT) == cut.at(HLT));
    pass.at(L1)  = true;//(value.at(L1) == cut.at(L1));






    unsigned int final_idx = 0;

    value.at(TRKLWithM) = 0;
    value.at(SignalCandidate)=0;
    value.at(Mu1PtCut)=0;
    value.at(Mu2PtCut)=0;
    value.at(Mu3PtCut)=0;

    value.at(TriggerMatch)=0;
    value.at(MuonID)=0;
    value.at(ThreeMuMass)=0;

    if(Ntp->NThreeMuons()>0){
      value.at(SignalCandidate) = Ntp->NThreeMuons();
      unsigned int mu1_idx = Ntp->ThreeMuonIndices(final_idx).at(0); 
      unsigned int mu2_idx = Ntp->ThreeMuonIndices(final_idx).at(1); 
      unsigned int mu3_idx = Ntp->ThreeMuonIndices(final_idx).at(2);
      // value.at(MuonID) =  (Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu1_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
      //     Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu2_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium) &&
      //     Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(mu3_idx),Ntp->MuonStandardSelectors::CutBasedIdMedium));
      //----------------  alternatively require two leading muons to be global and trailing muon to be tracker 
      unsigned int mu1_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
      unsigned int mu2_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
      unsigned int mu3_pt_idx=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);
      //


      //      std::cout<<"  "<<Ntp->Muon_isGlobalMuon(mu1_pt_idx)<<"  "<< Ntp->Muon_isTrackerMuon(mu1_pt_idx) << std::endl;

      value.at(MuonID) = (Ntp->Muon_isGlobalMuon(mu1_pt_idx) && 
			  Ntp->Muon_isGlobalMuon(mu2_pt_idx) &&
			  Ntp->Muon_isGlobalMuon(mu3_pt_idx) &&
			  Ntp->Muon_isPFMuon(mu1_pt_idx) &&
			  Ntp->Muon_isPFMuon(mu2_pt_idx) &&
			  Ntp->Muon_isPFMuon(mu3_pt_idx));
      value.at(TRKLWithM) = (Ntp->Muon_trackerLayersWithMeasurement(mu1_pt_idx) >= 7 ? 1:0 &&
			     Ntp->Muon_trackerLayersWithMeasurement(mu2_pt_idx) >= 7 ? 1:0 &&
			     Ntp->Muon_trackerLayersWithMeasurement(mu3_pt_idx) >= 7 ? 1:0 );
      //------------------------------------------------------------------------------------------------------

      value.at(Mu1PtCut) = Ntp->Muon_P4(mu1_pt_idx).Pt();
      value.at(Mu2PtCut) = Ntp->Muon_P4(mu2_pt_idx).Pt();
      value.at(Mu3PtCut) = Ntp->Muon_P4(mu3_pt_idx).Pt();

      vector<unsigned int> idx_vec;

      idx_vec.push_back(mu1_idx);
      idx_vec.push_back(mu2_idx);
      idx_vec.push_back(mu3_idx);

      unsigned int os_idx  = Ntp->SortedChargeMuons(idx_vec).at(0);
      unsigned int ss1_idx = Ntp->SortedChargeMuons(idx_vec).at(1);
      unsigned int ss2_idx = Ntp->SortedChargeMuons(idx_vec).at(2);

      TLorentzVector TauLV = Ntp->Muon_P4(mu1_idx)+Ntp->Muon_P4(mu2_idx)+Ntp->Muon_P4(mu3_idx);



      double M_osss1 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss1_idx)).M();
      double M_osss2 = (Ntp->Muon_P4(os_idx)+Ntp->Muon_P4(ss2_idx)).M();



      value.at(PhiVeto1) =  M_osss1;//M_osss1;
      value.at(PhiVeto2) =  M_osss2;
      value.at(OmegaVeto1) = M_osss1;
      value.at(OmegaVeto2) = M_osss2;


      unsigned int Muon_Eta_index_1=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
      unsigned int Muon_Eta_index_2=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
      unsigned int Muon_Eta_index_3=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);

      std::vector<unsigned int> EtaSortedIndices;


      EtaSortedIndices.push_back(Muon_Eta_index_1);
      EtaSortedIndices.push_back(Muon_Eta_index_2);
      EtaSortedIndices.push_back(Muon_Eta_index_3);
      EMR_tau_eta.at(t).Fill(Ntp->TauMassResolution(EtaSortedIndices,1,false),TauLV.Eta());
      if(Ntp->TauMassResolution(EtaSortedIndices,1,false) < 0.007) value.at(CutCategory)=1;
      if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > 0.007 && Ntp->TauMassResolution(EtaSortedIndices,1,false)< 0.01) value.at(CutCategory)=2;
      if(Ntp->TauMassResolution(EtaSortedIndices,1,false) > 0.01) value.at(CutCategory)=3;


      vector<TLorentzVector> trigobjTriplet;
      for (int i=0; i<Ntp->NTriggerObjects(); i++){
	TString name = Ntp->TriggerObject_name(i);
	if (!(name.Contains("tau3muDisplaced3muFltr"))) continue;
	TLorentzVector tmp;
	tmp.SetPtEtaPhiM(Ntp->TriggerObject_pt(i), Ntp->TriggerObject_eta(i), Ntp->TriggerObject_phi(i), PDG_Var::Muon_mass());
	trigobjTriplet.push_back(tmp);
      }

      std::vector<TLorentzVector> muonTriplet;
      muonTriplet.push_back(Ntp->Muon_P4(mu1_pt_idx));
      muonTriplet.push_back(Ntp->Muon_P4(mu2_pt_idx));
      muonTriplet.push_back(Ntp->Muon_P4(mu3_pt_idx));

      bool triggerCheck = false;
      if (trigobjTriplet.size()>=3) triggerCheck = Ntp->triggerMatchTriplet(muonTriplet, trigobjTriplet).first;
      value.at(TriggerMatch) = triggerCheck;


 
      //      value.at(TriggerMatchMu) = Ntp->ThreeMuons_TriggerMatch_dR(final_idx).at(0);
      //      value.at(TriggerMatchMu2) = Ntp->ThreeMuons_TriggerMatch_dR(final_idx).at(1);
      //      value.at(TriggerMatchMu3) = Ntp->ThreeMuons_TriggerMatch_dR(final_idx).at(2);
     
      value.at(ThreeMuMass) = TauLV.M();
    }
    
    pass.at(SignalCandidate) = (value.at(SignalCandidate) == cut.at(SignalCandidate));
    pass.at(Mu1PtCut) = (value.at(Mu1PtCut) > cut.at(Mu1PtCut));
    pass.at(Mu2PtCut) = (value.at(Mu2PtCut) > cut.at(Mu2PtCut));
    pass.at(Mu3PtCut) = (value.at(Mu3PtCut) > cut.at(Mu3PtCut));
    pass.at(MuonID)  =  true;//(value.at(MuonID)     == cut.at(MuonID));
    pass.at(TRKLWithM) = true;//(value.at(TRKLWithM) == cut.at(TRKLWithM));
    pass.at(TriggerMatch) = (value.at(TriggerMatch) == cut.at(TriggerMatch));
    pass.at(PhiVeto1) = true;//(value.at(PhiVeto1) < 0.98 || value.at(PhiVeto1) > 1.06 );
    pass.at(OmegaVeto1) = true;//(value.at(OmegaVeto1) < 0.742 || value.at(OmegaVeto1) > 0.822 );
    pass.at(PhiVeto2) = true;//(value.at(PhiVeto2) < 0.98 || value.at(PhiVeto2) > 1.06 );
    pass.at(OmegaVeto2) = true;//(value.at(OmegaVeto2) < 0.742 || value.at(OmegaVeto2) > 0.822 );
    pass.at(CutCategory) = true;//( value.at(CutCategory) == cut.at(CutCategory) );

    //    if(id!=1) pass.at(ThreeMuMass) = true;
    //    else  pass.at(ThreeMuMass) = true;//( (value.at(ThreeMuMass) > tauMinSideBand_ && value.at(ThreeMuMass) < tauMinMass_)  || (value.at(ThreeMuMass)> tauMaxMass_ && value.at(ThreeMuMass) < tauMaxSideBand_));
    pass.at(ThreeMuMass) = true;//( (value.at(ThreeMuMass) > tauMinSideBand_) &&  (value.at(ThreeMuMass) < tauMaxSideBand_));

    double wobs=1;
    double w; 
    
    //    if(!Ntp->isData()){w = (puWeights->GetBinContent(Ntp->TruthNumberOfInteraction()));   }//w = 1; /*Ntp->PUReweight(); */} //  No weights to data
    if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
    else{w=1;}

    bool status=AnalysisCuts(t,w,wobs);
    
    nPVx.at(t).Fill(Ntp->NVtx(),w);

    if(status){

      /*
      int NJets(){return Ntp->Jet_p4->size();}
      TLorentzVector Jet_P4(unsigned int i){return TLorentzVector(Ntp->Jet_p4->at(i)->at(1), Ntp->Jet_p4->at(i)->at(2), Ntp->Jet_p4->at(i)->at(3), Ntp->Jet_p4->at(i)->at(0));}
      double JetBTagCVSB(unsigned int i){return Ntp->Jet_BTagCVSB->at(i);}
      double JetBTagMVA(unsigned int i){return Ntp->Jet_BTagMVA->at(i);}
      double JetBTagCSV(unsigned int i){return Ntp->Jet_BTagCSV->at(i);}
      */


      // unsigned int Muon_index_1 =  Ntp->ThreeMuonIndices(final_idx).at(0);
      // unsigned int Muon_index_2 =  Ntp->ThreeMuonIndices(final_idx).at(1);
      // unsigned int Muon_index_3 =  Ntp->ThreeMuonIndices(final_idx).at(2);


      //      //      double     Vertex_2MuonsIsoTrack_KF_Chi2(unsigned int i, bool channel=false){
      //      TVector3   Vertex_2MuonsIsoTrack_KF_pos(unsigned int i, bool channel=false){
      
      

      //      std::cout<<   " 2mt  "<< Ntp->Vertex_2MuonsIsoTrack_KF_Chi2(final_idx) << "   3m   "<<Ntp->Vertex_signal_KF_Chi2(final_idx) << std::endl;
      var_id = id;
      var_EventNumber = Ntp->EventNumber();
      var_Vertex2muTrkKF = Ntp->Vertex_2MuonsIsoTrack_KF_Chi2(final_idx);
      //      std::cout<<"  var_Vertex2muTrkKF   ever -1 ?  " << var_Vertex2muTrkKF << std::endl;
      //      std::cout<<" Event NUmber  "<< var_EventNumber << std::endl;
      Vertex2muTrkKF.at(t).Fill(Ntp->Vertex_2MuonsIsoTrack_KF_Chi2(final_idx) ,w);
      unsigned int Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
      unsigned int Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
      unsigned int Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);

      if(Ntp->Vertex_2MuonsIsoTrack_KF_Chi2(final_idx)!=-1){
	Dist2muTrkKF3Mu.at(t).Fill( (Ntp->Vertex_2MuonsIsoTrack_KF_pos(final_idx) - Ntp->Vertex_Signal_KF_pos(final_idx)).Mag() ,w);
	var_Dist2muTrkKF3Mu = (Ntp->Vertex_2MuonsIsoTrack_KF_pos(final_idx) - Ntp->Vertex_Signal_KF_pos(final_idx)).Mag();
      }else{
	var_Dist2muTrkKF3Mu = -1;
      }


      if(Ntp->Vertex_2MuonsIsoTrack_KF_Chi2(final_idx) < Ntp->Vertex_signal_KF_Chi2(final_idx)){
	VertexQualitySeparator.at(t).Fill(0.,w);
	var_VertexQualitySeparator = 0.;
      }else{
	VertexQualitySeparator.at(t).Fill(1.,w);
	var_VertexQualitySeparator = 1.;
      }




      unsigned int SS1RandomIndex(0);
      unsigned int SS2RandomIndex(0);


      vector<unsigned int> idx_vec;

      idx_vec.push_back( Ntp->ThreeMuonIndices(final_idx).at(0));
      idx_vec.push_back( Ntp->ThreeMuonIndices(final_idx).at(1));
      idx_vec.push_back( Ntp->ThreeMuonIndices(final_idx).at(2));


      unsigned int Muon_index_os  = Ntp->SortedChargeMuons(idx_vec).at(0);
      unsigned int Muon_index_ss1 = Ntp->SortedChargeMuons(idx_vec).at(1);
      unsigned int Muon_index_ss2 = Ntp->SortedChargeMuons(idx_vec).at(2);

      TLorentzVector MuonLV_OS  = Ntp->Muon_P4(Muon_index_os);
      TLorentzVector MuonLV_SS1 = Ntp->Muon_P4(Muon_index_ss1);
      TLorentzVector MuonLV_SS2 = Ntp->Muon_P4(Muon_index_ss2);
      //*** With such sorting pT(ss1) > pT(ss2)
      TLorentzVector MuonOS  = Ntp->Muon_P4(Muon_index_os);  
      TLorentzVector MuonSS1 = Ntp->Muon_P4(Muon_index_ss1);
      TLorentzVector MuonSS2 = Ntp->Muon_P4(Muon_index_ss2);

      TLorentzVector MuonOS_Rotated = MuonOS; MuonOS_Rotated.SetVect(Ntp->Rotate(MuonOS_Rotated.Vect(), (MuonOS+MuonSS1+MuonSS2).Vect()));
      TLorentzVector MuonSS1_Rotated = MuonSS1; MuonSS1_Rotated.SetVect(Ntp->Rotate(MuonSS1_Rotated.Vect(), (MuonOS+MuonSS1+MuonSS2).Vect()));
      TLorentzVector MuonSS2_Rotated = MuonSS2; MuonSS2_Rotated.SetVect(Ntp->Rotate(MuonSS2_Rotated.Vect(), (MuonOS+MuonSS1+MuonSS2).Vect()));
      TLorentzVector ThreeMuon_Rotated = MuonOS + MuonSS1 +MuonSS2; ThreeMuon_Rotated.SetVect(Ntp->Rotate(ThreeMuon_Rotated.Vect(), (MuonOS+MuonSS1+MuonSS2).Vect()));


      TLorentzVector MuonOS_TRF  = Ntp->Boost(MuonOS_Rotated, ThreeMuon_Rotated);    // Muon in Tau REst Frame with Z axis alligned on tau directions
      TLorentzVector MuonSS1_TRF = Ntp->Boost(MuonSS1_Rotated, ThreeMuon_Rotated);
      TLorentzVector MuonSS2_TRF = Ntp->Boost(MuonSS2_Rotated, ThreeMuon_Rotated);



      TLorentzVector MuonOS_R1_RF_rotated = MuonOS_TRF; MuonOS_R1_RF_rotated.SetVect(Ntp->Rotate(MuonOS_R1_RF_rotated.Vect(), (MuonOS_TRF + MuonSS1_TRF).Vect())   );
      TLorentzVector MuonSS1_R1_RF_rotated = MuonSS1_TRF; MuonSS1_R1_RF_rotated.SetVect(Ntp->Rotate(MuonSS1_R1_RF_rotated.Vect(), (MuonOS_TRF + MuonSS1_TRF).Vect())   );



      TLorentzVector MuonOS_R2_RF_rotated = MuonOS_TRF; MuonOS_R2_RF_rotated.SetVect(Ntp->Rotate(MuonOS_R2_RF_rotated.Vect(), (MuonOS_TRF + MuonSS2_TRF).Vect())   );
      TLorentzVector MuonSS2_R2_RF_rotated = MuonSS2_TRF; MuonSS2_R2_RF_rotated.SetVect(Ntp->Rotate(MuonSS2_R2_RF_rotated.Vect(), (MuonOS_TRF + MuonSS2_TRF).Vect()  )   );

      TLorentzVector MuonOS_R1_RF = Ntp->Boost(MuonOS_R1_RF_rotated, MuonOS_R1_RF_rotated+MuonSS1_R1_RF_rotated);
      TLorentzVector MuonSS1_R1_RF=Ntp->Boost(MuonSS1_R1_RF_rotated, MuonOS_R1_RF_rotated+MuonSS1_R1_RF_rotated);
    
      TLorentzVector MuonOS_R2_RF = Ntp->Boost(MuonOS_R2_RF_rotated, MuonOS_R2_RF_rotated+MuonSS2_R2_RF_rotated);
      TLorentzVector MuonSS2_R2_RF=Ntp->Boost(MuonSS2_R2_RF_rotated, MuonOS_R2_RF_rotated+MuonSS2_R2_RF_rotated);

      double R1_RF_Theta =     MuonOS_R1_RF.Vect().Dot( (MuonOS_R1_RF_rotated+MuonSS1_R1_RF_rotated).Vect() ) * (1./MuonOS_R1_RF.Vect().Mag()/(MuonOS_R1_RF_rotated+MuonSS1_R1_RF_rotated).Vect().Mag());
      double R2_RF_Theta =     MuonOS_R2_RF.Vect().Dot( (MuonOS_R2_RF_rotated+MuonSS2_R2_RF_rotated).Vect() ) * (1./MuonOS_R2_RF.Vect().Mag()/(MuonOS_R2_RF_rotated+MuonSS2_R2_RF_rotated).Vect().Mag());
      
      OSSS1Angle_RRF.at(t).Fill(R1_RF_Theta,1.);
      OSSS2Angle_RRF.at(t).Fill(R2_RF_Theta,1.);


      OSSS1Angle_TRF.at(t).Fill(( MuonOS_TRF.Vect() * MuonSS1_TRF.Vect() )*(1/MuonOS_TRF.Vect().Mag())*(1/MuonSS1_TRF.Vect().Mag()),1.);
      OSSS2Angle_TRF.at(t).Fill(( MuonOS_TRF.Vect() * MuonSS2_TRF.Vect() )*(1/MuonOS_TRF.Vect().Mag())*(1/MuonSS2_TRF.Vect().Mag()),1.);

      double costheta_TRF_SSSS;
      double costheta_TRF_OSSS;

      TVector3 n_tau(0,0,1);




    //    (Ntp->Rotate(MuonOS.Vect(),MuonOS.Vect())).Print();
      
      if(MuonSS1_TRF.Vect().Mag() > MuonSS2_TRF.Vect().Mag()){

	costheta_TRF_SSSS = MuonSS1_TRF.Vect().Cross(MuonSS2_TRF.Vect()).Dot(n_tau) *(1/ (MuonSS1_TRF.Vect().Cross(MuonSS2_TRF.Vect())).Mag());
	costheta_TRF_OSSS = MuonOS_TRF.Vect().Cross(MuonSS1_TRF.Vect()).Dot(n_tau) *(1/ ( MuonOS_TRF.Vect().Cross(MuonSS1_TRF.Vect())).Mag());

      } else {
	
	costheta_TRF_SSSS = MuonSS2_TRF.Vect().Cross(MuonSS1_TRF.Vect()).Dot(n_tau) *(1/ (MuonSS2_TRF.Vect().Cross(MuonSS1_TRF.Vect())).Mag());
	costheta_TRF_OSSS = MuonOS_TRF.Vect().Cross(MuonSS2_TRF.Vect()).Dot(n_tau) *(1/  (MuonOS_TRF.Vect().Cross(MuonSS2_TRF.Vect())).Mag());
	
      }
      

      cTheta_TRF_SSSS.at(t).Fill(costheta_TRF_SSSS,1.);
      cTheta_TRF_OSSS.at(t).Fill(costheta_TRF_OSSS,1.);

      cTheta_MuonOS_TauPol_TRF.at(t).Fill(MuonOS_TRF.Vect().Dot(n_tau)*(1./MuonOS_TRF.Vect().Mag()),1.);
      pTMu1OverMass_TRF.at(t).Fill(MuonOS_TRF.P()/(MuonOS_TRF + MuonSS1_TRF + MuonSS2_TRF).M(),1.);

      var_cTheta_TRF_SSSS = costheta_TRF_SSSS;
      var_cTheta_TRF_OSSS = costheta_TRF_OSSS;
      var_cTheta_MuonOS_TauPol_TRF = MuonOS_TRF.Vect().Dot(n_tau)*(1./MuonOS_TRF.Vect().Mag());
      var_pTMu1OverMass_TRF = MuonOS_TRF.P()/(MuonOS_TRF + MuonSS1_TRF + MuonSS2_TRF).M();
      var_costheta_TRF_SSSS = costheta_TRF_SSSS;
      var_costheta_TRF_OSS = costheta_TRF_OSSS;
      var_OSSS1Angle_TRF = ( MuonOS_TRF.Vect() * MuonSS1_TRF.Vect() )*(1/MuonOS_TRF.Vect().Mag())*(1/MuonSS1_TRF.Vect().Mag());
      var_OSSS2Angle_TRF = ( MuonOS_TRF.Vect() * MuonSS2_TRF.Vect() )*(1/MuonOS_TRF.Vect().Mag())*(1/MuonSS2_TRF.Vect().Mag());

      var_OSSS1Angle_RRF = R1_RF_Theta;
      var_OSSS2Angle_RRF = R2_RF_Theta;






      std::vector<unsigned int> Indices;
      Indices.push_back(Muon_index_ss1);
      Indices.push_back(Muon_index_ss2);

      float random_muon_index = rndm.Uniform();
      if(random_muon_index >= 0.5 ){SS1RandomIndex =  Indices.at(0); SS2RandomIndex = Indices.at(1) ; }
      if(random_muon_index <  0.5 ){SS1RandomIndex =  Indices.at(1); SS2RandomIndex = Indices.at(0) ; }

      TLorentzVector MuonLV_RandomSS1 = Ntp->Muon_P4(SS1RandomIndex);
      TLorentzVector MuonLV_RandomSS2 = Ntp->Muon_P4(SS2RandomIndex);
      //      std::cout<<" rnd "<< random_muon_index <<std::endl;


      PairMass1NoSorting.at(t).Fill((MuonLV_OS+MuonLV_RandomSS1).M(),w);
      PairMass2NoSorting.at(t).Fill((MuonLV_OS+MuonLV_RandomSS2).M(),w);
      MuMuMassNoSorting.at(t).Fill((MuonLV_OS+MuonLV_RandomSS2).M(),(MuonLV_OS+MuonLV_RandomSS1).M());


      PairMass1PTSorting.at(t).Fill((MuonLV_OS+MuonLV_SS1).M(),w);
      PairMass2PTSorting.at(t).Fill((MuonLV_OS+MuonLV_SS2).M(),w);
      MuMuMassPTSorting.at(t).Fill((MuonLV_OS+MuonLV_SS1).M(),(MuonLV_OS+MuonLV_SS2).M());


      //      std::cout<<" comapre masses "<< (MuonLV_OS+MuonLV_SS1).M() << "  "<<(MuonLV_OS+MuonLV_RandomSS1).M() << std::endl;
      //      std::cout<<" comapre masses "<< (MuonLV_OS+MuonLV_SS2).M() << "  "<<(MuonLV_OS+MuonLV_RandomSS2).M() << std::endl;
      
      //      std::cout<<" print out masses   "<<(MuonLV_OS+MuonLV_RandomSS1).M()  <<"    " <<(MuonLV_OS+MuonLV_RandomSS2).M() <<std::endl;


      //      std::cout<<"  print out muons LV   "<<std::endl; 
      
      //      std::cout<<"  "<<MuonLV_OS.Pt() <<std::endl;
      //      std::cout<<"  "<<MuonLV_RandomSS1.Pt()<<std::endl;
      //      std::cout<<"  "<<MuonLV_RandomSS2.Pt()<<std::endl;
	    


      //      std::cout<<"  print out muons pT Sortes os,ss1,ss2, sorted   "<<std::endl; 
      //      std::cout<<"  "<<MuonLV_OS.Pt()<<std::endl;
      //      std::cout<<"  "<<MuonLV_SS1.Pt()<<std::endl;
      //      std::cout<<"  "<<MuonLV_SS2.Pt()<<std::endl;


      if(MuonLV_OS.DeltaR(MuonLV_SS1) > MuonLV_OS.DeltaR(MuonLV_SS2)){
	PairMass1AllignedSorting.at(t).Fill((MuonLV_OS+MuonLV_SS2).M(),w);
	PairMass2AllignedSorting.at(t).Fill((MuonLV_OS+MuonLV_SS1).M(),w);
	var_mass12_dRsorting = (MuonLV_OS+MuonLV_SS2).M();
	var_mass13_dRsorting = (MuonLV_OS+MuonLV_SS1).M();

	MuMuMassAllignedSorting.at(t).Fill((MuonLV_OS+MuonLV_SS1).M(),(MuonLV_OS+MuonLV_SS2).M());
      }else{
	PairMass1AllignedSorting.at(t).Fill((MuonLV_OS+MuonLV_SS1).M(),w);
	PairMass2AllignedSorting.at(t).Fill((MuonLV_OS+MuonLV_SS2).M(),w);
	var_mass12_dRsorting = (MuonLV_OS+MuonLV_SS1).M();
	var_mass13_dRsorting = (MuonLV_OS+MuonLV_SS2).M();

	MuMuMassAllignedSorting.at(t).Fill((MuonLV_OS+MuonLV_SS2).M(),(MuonLV_OS+MuonLV_SS1).M());
      }
      

      //      std::cout<<"pT  OS,SS1,SS2 "<< MuonLV_OS.Pt() << "  "<< MuonLV_SS1.Pt() << "  "<<  MuonLV_SS2.Pt() << std::endl;
      //      std::cout<<"ch  OS,SS1,SS2 "<< Ntp->Muon_charge(Muon_index_os) << "  "<< Ntp->Muon_charge(Muon_index_ss1) << "  "<< Ntp->Muon_charge(Muon_index_ss2)  << std::endl;

      var_mass12 = (MuonLV_OS+MuonLV_SS1).M();
      var_mass13 = (MuonLV_OS+MuonLV_SS2).M();


      std::vector<unsigned int> indices;
      indices.push_back(Muon_index_1);
      indices.push_back(Muon_index_2);
      indices.push_back(Muon_index_3);
      

      TLorentzVector Muon1LV = Ntp->Muon_P4(Muon_index_1);
      TLorentzVector Muon2LV = Ntp->Muon_P4(Muon_index_2);
      TLorentzVector Muon3LV = Ntp->Muon_P4(Muon_index_3);
      TLorentzVector TauLV = Muon1LV + Muon2LV + Muon3LV;
      TLorentzVector TauRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,0)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,1)+Ntp->Vertex_signal_KF_refittedTracksP4(final_idx,2);



      int OH_jet_count(0);
      int SH_jet_count(0);
      int jet_index_closest_to_tau(-1);
      int jet_index_closest_to_tau_opposite_in_phi(-1);
      double TJDeltaR(99.);
      double OHDeltaPhi(TMath::Pi());
      for(int ijet= 0; ijet <  Ntp->NJets(); ijet++){
	if(TauLV.DeltaR(Ntp->Jet_P4(ijet)) < TJDeltaR){
	  TJDeltaR = TauLV.DeltaR(Ntp->Jet_P4(ijet));
	  jet_index_closest_to_tau = ijet;
	}


	if(fabs(Ntp->Jet_P4(ijet).Phi()  - TauLV.Phi()) < TMath::Pi() ){
	  SH_jet_count++;
	  BTagCSVSH.at(t).Fill(Ntp->JetBTagCSV(ijet),w);
	  BTagMVASH.at(t).Fill(Ntp->JetBTagMVA(ijet),w);
	  BTagCVSBSH.at(t).Fill(Ntp->JetBTagCVSB(ijet),w);
	}


	if(fabs(Ntp->Jet_P4(ijet).Phi()  - TauLV.Phi()) > TMath::Pi() ){
	  OH_jet_count++;  
	  BTagCSVOH.at(t).Fill(Ntp->JetBTagCSV(ijet),w);
	  BTagMVAOH.at(t).Fill(Ntp->JetBTagMVA(ijet),w);
	  BTagCVSBOH.at(t).Fill(Ntp->JetBTagCVSB(ijet),w);
	}


      }
      for(int ijet= 0; ijet <  Ntp->NJets(); ijet++){
	if(jet_index_closest_to_tau!=-1){
	  //	  if(fabs(Ntp->Jet_P4(ijet).Phi() - fabs(Ntp->Jet_P4(jet_index_closest_to_tau).Phi() - TMath::Pi())) < OHDeltaPhi  )  {
	  if(fabs(TMath::Pi() - fabs(Ntp->Jet_P4(ijet).Phi() - Ntp->Jet_P4(jet_index_closest_to_tau).Phi() )) < OHDeltaPhi  )  {
	    OHDeltaPhi = fabs(TMath::Pi() - fabs(Ntp->Jet_P4(ijet).Phi() - Ntp->Jet_P4(jet_index_closest_to_tau).Phi() ));
	    jet_index_closest_to_tau_opposite_in_phi = ijet;
	  }

	}
      }

      NBJet4piSH.at(t).Fill(SH_jet_count,w);
      NBJet4piOH.at(t).Fill(OH_jet_count,w);
      
      if(jet_index_closest_to_tau!=-1){

	BTagCSVSHMatchedToTau.at(t).Fill(Ntp->JetBTagCSV(jet_index_closest_to_tau),w);
	BTagMVASHMatchedToTau.at(t).Fill(Ntp->JetBTagMVA(jet_index_closest_to_tau),w);
	BTagCVSBSHMatchedToTau.at(t).Fill(Ntp->JetBTagCVSB(jet_index_closest_to_tau),w);

      }
      if(jet_index_closest_to_tau_opposite_in_phi!=-1){

	BTagCSVOHMatchedToTau.at(t).Fill(Ntp->JetBTagCSV(jet_index_closest_to_tau_opposite_in_phi),w);
	BTagMVAOHMatchedToTau.at(t).Fill(Ntp->JetBTagMVA(jet_index_closest_to_tau_opposite_in_phi),w);
	BTagCVSBOHMatchedToTau.at(t).Fill(Ntp->JetBTagCVSB(jet_index_closest_to_tau_opposite_in_phi),w);

      }


      //------------- Muon ID

      for(unsigned int iMuSelector=0; iMuSelector< Ntp->MuonStandardSelectorBitMask(Muon_index_1).size(); iMuSelector++ ){
	if(Ntp->MuonStandardSelectorBitMask(Muon_index_1).at(iMuSelector)==1)  Muon1StandardSelectorPass.at(t).Fill(iMuSelector,1);
	if(Ntp->MuonStandardSelectorBitMask(Muon_index_1).at(iMuSelector)!=1)  Muon1StandardSelectorFail.at(t).Fill(iMuSelector,1);
      }

      for(unsigned int iMuSelector=0; iMuSelector< Ntp->MuonStandardSelectorBitMask(Muon_index_2).size(); iMuSelector++ ){
	if(Ntp->MuonStandardSelectorBitMask(Muon_index_2).at(iMuSelector)==1)  Muon2StandardSelectorPass.at(t).Fill(iMuSelector,1);
	if(Ntp->MuonStandardSelectorBitMask(Muon_index_2).at(iMuSelector)!=1)  Muon2StandardSelectorFail.at(t).Fill(iMuSelector,1);
      }

      for(unsigned int iMuSelector=0; iMuSelector< Ntp->MuonStandardSelectorBitMask(Muon_index_3).size(); iMuSelector++ ){
	if(Ntp->MuonStandardSelectorBitMask(Muon_index_3).at(iMuSelector)==1)  Muon3StandardSelectorPass.at(t).Fill(iMuSelector,1);
	if(Ntp->MuonStandardSelectorBitMask(Muon_index_3).at(iMuSelector)!=1)  Muon3StandardSelectorFail.at(t).Fill(iMuSelector,1);
      }

      for(unsigned int imu=0; imu<3;imu++){

	mu_combinedQuality_chi2LocalMomentum=Ntp->Muon_combinedQuality_chi2LocalMomentum(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu));
	mu_combinedQuality_chi2LocalPosition=Ntp->Muon_combinedQuality_chi2LocalPosition(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu));
	mu_combinedQuality_staRelChi2=Ntp->Muon_combinedQuality_staRelChi2(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu));
	mu_combinedQuality_trkRelChi2=Ntp->Muon_combinedQuality_trkRelChi2(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu));
	mu_combinedQuality_globalDeltaEtaPhi=Ntp->Muon_combinedQuality_globalDeltaEtaPhi(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu));
	mu_combinedQuality_trkKink=log(Ntp->Muon_combinedQuality_glbKink(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu)));
	mu_combinedQuality_glbKink=log(Ntp->Muon_combinedQuality_trkKink(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu)));
	mu_combinedQuality_glbTrackProbability=Ntp->Muon_combinedQuality_glbTrackProbability(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu));
	mu_Numberofvalidtrackerhits=Ntp->Muon_numberofValidPixelHits(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu));
	mu_Numberofvalidpixelhits=Ntp->Muon_innerTrack_numberOfValidTrackerHits(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu));
	mu_trackerLayersWithMeasurement = Ntp->Muon_trackerLayersWithMeasurement(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu));
	mu_validMuonHitComb=Ntp->Muon_hitPattern_numberOfValidMuonHits(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu));
	mu_numberOfMatchedStations=Ntp->Muon_numberOfMatchedStations(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu));
	mu_segmentCompatibility=Ntp->Muon_segmentCompatibility(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu));
	mu_timeAtIpInOutErr=Ntp->Muon_timeAtIpInOutErr(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu));
	mu_GLnormChi2=Ntp->Muon_normChi2(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu));
	
	mu_innerTrack_normalizedChi2=Ntp->Muon_innerTrack_normalizedChi2(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu));
	mu_outerTrack_normalizedChi2= Ntp->Muon_outerTrack_normalizedChi2(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu));
	mu_innerTrack_validFraction=Ntp->Muon_innerTrack_validFraction(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu));


	if(fabs(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(final_idx)).at(imu)).Eta()) < 1.2    )
	  {
	    if(imu==0)
	      {
		Muon1DetID = readerMuIDBarrel->EvaluateMVA("BDT");
	      }
	    if(imu==1)
	      {
		Muon2DetID = readerMuIDBarrel->EvaluateMVA("BDT");
	      }

	    if(imu==2)
	      {
		Muon3DetID = readerMuIDBarrel->EvaluateMVA("BDT");
	      }
	  }
	else
	  {
	    if(imu==0)
	      {
		Muon1DetID = readerMuIDEndcap->EvaluateMVA("BDT");
	      }
	    if(imu==1)
	      {
	       Muon2DetID = readerMuIDEndcap->EvaluateMVA("BDT");
	      }
	    if(imu==2)
	      {
	       Muon3DetID= readerMuIDEndcap->EvaluateMVA("BDT");
	      }
	  }
	
      }

      Muon1MVAID.at(t).Fill(Muon1DetID);
      Muon2MVAID.at(t).Fill(Muon2DetID);
      Muon3MVAID.at(t).Fill(Muon3DetID);

      var_Muon1DetID = Muon1DetID; 
      var_Muon2DetID = Muon2DetID;
      var_Muon3DetID = Muon3DetID;




      //------------- Muon ID




      unsigned int Muon_Eta_index_1=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(0);
      unsigned int Muon_Eta_index_2=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(1);
      unsigned int Muon_Eta_index_3=Ntp->SortedEtaMuons(Ntp->ThreeMuonIndices(final_idx)).at(2);
      
      std::vector<unsigned int> EtaSortedIndices;
      
      
      EtaSortedIndices.push_back(Muon_Eta_index_1);
      EtaSortedIndices.push_back(Muon_Eta_index_2);
      EtaSortedIndices.push_back(Muon_Eta_index_3);
      EMR_tau_eta.at(t).Fill(Ntp->TauMassResolution(EtaSortedIndices,1,false),TauLV.Eta());
      EventMassResolution_PtEtaPhi.at(t).Fill(Ntp->TauMassResolution(EtaSortedIndices,1,false),w);

      float tauMassRes = Ntp->TauMassResolution(EtaSortedIndices,1,false);

      Muon1Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Pt(),w);
      Muon2Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_2).Pt(),w);
      Muon3Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_3).Pt(),w);

      Muon1P.at(t).Fill(Ntp->Muon_P4(Muon_index_1).P(),w);
      Muon2P.at(t).Fill(Ntp->Muon_P4(Muon_index_2).P(),w);
      Muon3P.at(t).Fill(Ntp->Muon_P4(Muon_index_3).P(),w);

  
      TauEta.at(t).Fill(TauLV.Eta(),w);
      TauPt.at(t).Fill(TauLV.Pt(),w);
      TauP.at(t).Fill(TauLV.P(),w);
    


      float MinSegmentCompatibility = std::min({Ntp->Muon_segmentCompatibility(Muon_index_1),Ntp->Muon_segmentCompatibility(Muon_index_2),Ntp->Muon_segmentCompatibility(Muon_index_3)});
      float MaxSegmentCompatibility = std::max({Ntp->Muon_segmentCompatibility(Muon_index_1),Ntp->Muon_segmentCompatibility(Muon_index_2),Ntp->Muon_segmentCompatibility(Muon_index_3)});
      float MinMIPLikelihood = std::min({Ntp->Muon_caloCompatibility(Muon_index_1),Ntp->Muon_caloCompatibility(Muon_index_2),Ntp->Muon_caloCompatibility(Muon_index_3)});
      float MaxMIPLikelihood = std::max({Ntp->Muon_caloCompatibility(Muon_index_1),Ntp->Muon_caloCompatibility(Muon_index_2),Ntp->Muon_caloCompatibility(Muon_index_3)});
      //------------------ calculate var_svpvTauAngle ---------------------
      TVector3 vec_sv = Ntp->Vertex_Signal_KF_pos(final_idx);
      TVector3 vec_pv(0,0,0);
      if (Ntp->Vertex_RefitPVisValid(final_idx)==1){
	vec_pv = Ntp->Vertex_MatchedRefitPrimaryVertex(final_idx);
      }
      TVector3 vec_tau = TauLV.Vect();
      TVector3 d_pvsv = vec_sv - vec_pv;

      TVector3 SVPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(final_idx),Ntp->Vertex_MatchedPrimaryVertex(final_idx));

      TVector3 Mu1ImpactPV = Ntp->SVPVDirection(Ntp->Muon_Poca(Muon_index_1),Ntp->Vertex_MatchedPrimaryVertex(final_idx));
      TVector3 Mu2ImpactPV = Ntp->SVPVDirection(Ntp->Muon_Poca(Muon_index_2),Ntp->Vertex_MatchedPrimaryVertex(final_idx));
      TVector3 Mu3ImpactPV = Ntp->SVPVDirection(Ntp->Muon_Poca(Muon_index_3),Ntp->Vertex_MatchedPrimaryVertex(final_idx));



      //---------------------------------------------------------------
      //------------------ calculate var_flightLenSig ---------------------
      TMatrixTSym<double> fls_PVcov = Ntp->Vertex_PrimaryVertex_Covariance(final_idx);
      TMatrixTSym<double> fls_SVcov = Ntp->Vertex_Signal_KF_Covariance(final_idx);
      
 

      //      if (Ntp->Vertex_RefitPVisValid(final_idx)==1)
	{
	
	  Muon1ImpactAngle.at(t).Fill(SVPV*Mu1ImpactPV*(1/Mu1ImpactPV.Mag()/SVPV.Mag()),w );
	  Muon2ImpactAngle.at(t).Fill(SVPV*Mu2ImpactPV*(1/Mu2ImpactPV.Mag()/SVPV.Mag()),w );
	  Muon3ImpactAngle.at(t).Fill(SVPV*Mu3ImpactPV*(1/Mu3ImpactPV.Mag()/SVPV.Mag()),w );
	  MinMuonImpacAngle.at(t).Fill(std::min({SVPV*Mu1ImpactPV*(1/Mu1ImpactPV.Mag()/SVPV.Mag()),SVPV*Mu2ImpactPV*(1/Mu2ImpactPV.Mag()/SVPV.Mag()),SVPV*Mu3ImpactPV*(1/Mu3ImpactPV.Mag()/SVPV.Mag())}),w);
	  MaxMuonImpacAngle.at(t).Fill(std::max({SVPV*Mu1ImpactPV*(1/Mu1ImpactPV.Mag()/SVPV.Mag()),SVPV*Mu2ImpactPV*(1/Mu2ImpactPV.Mag()/SVPV.Mag()),SVPV*Mu3ImpactPV*(1/Mu3ImpactPV.Mag()/SVPV.Mag())}),w);
	
	  var_Muon1ImpactAngle  = SVPV*Mu1ImpactPV*(1/Mu1ImpactPV.Mag()/SVPV.Mag());
	  var_Muon2ImpactAngle  = SVPV*Mu2ImpactPV*(1/Mu2ImpactPV.Mag()/SVPV.Mag());
	  var_Muon3ImpactAngle  = SVPV*Mu3ImpactPV*(1/Mu3ImpactPV.Mag()/SVPV.Mag());
	  var_MinMuonImpactAngle = std::min({SVPV*Mu1ImpactPV*(1/Mu1ImpactPV.Mag()/SVPV.Mag()),SVPV*Mu2ImpactPV*(1/Mu2ImpactPV.Mag()/SVPV.Mag()),SVPV*Mu3ImpactPV*(1/Mu3ImpactPV.Mag()/SVPV.Mag())});
	  var_MaxMuonImpactAngle = std::max({SVPV*Mu1ImpactPV*(1/Mu1ImpactPV.Mag()/SVPV.Mag()),SVPV*Mu2ImpactPV*(1/Mu2ImpactPV.Mag()/SVPV.Mag()),SVPV*Mu3ImpactPV*(1/Mu3ImpactPV.Mag()/SVPV.Mag())});




	// ----- Fill the histograms -----
	SVPVTauDirAngle.at(t).Fill(var_svpvTauAngle,w);
	FLSignificance.at(t).Fill(var_flightLenSig,w);
	FL.at(t).Fill(SVPV.Mag(),w);
	VertexChi2KF.at(t).Fill(var_vertexKFChi2,w);

	Muon_segmentCompatibility_mu1.at(t).Fill(Ntp->Muon_segmentCompatibility(Muon_index_1),w);
	Muon_segmentCompatibility_mu2.at(t).Fill(Ntp->Muon_segmentCompatibility(Muon_index_2),w);
	Muon_segmentCompatibility_mu3.at(t).Fill(Ntp->Muon_segmentCompatibility(Muon_index_3),w);
	Muon_segmentCompatibility_min.at(t).Fill(MinSegmentCompatibility,w);
	Muon_segmentCompatibility_max.at(t).Fill(MaxSegmentCompatibility,w);
	Muon_ECALCompatibility_mu1.at(t).Fill(Ntp->Muon_caloCompatibility(Muon_index_1),w);
	Muon_ECALCompatibility_mu2.at(t).Fill(Ntp->Muon_caloCompatibility(Muon_index_2),w);
	Muon_ECALCompatibility_mu3.at(t).Fill(Ntp->Muon_caloCompatibility(Muon_index_3),w);
	Muon_ECALCompatibility_min.at(t).Fill(MinMIPLikelihood,w);
	Muon_ECALCompatibility_max.at(t).Fill(MaxMIPLikelihood,w);
	Isolation_NTracks.at(t).Fill(Ntp->Isolation_NTracks(final_idx),w);
	Isolation_RelPt.at(t).Fill(Ntp->Isolation_RelPt(final_idx),w);
	Isolation_maxdxy.at(t).Fill(Ntp->Isolation_maxdy(final_idx),w);

        if(Ntp->Vertex_2MuonsIsoTrack_KF_Chi2(final_idx) > 0){
          Vertex2muTrkKFChi2.at(t).Fill(Ntp->Vertex_2MuonsIsoTrack_KF_Chi2(final_idx),1);
          if(Ntp->Vertex_signal_KF_Chi2(final_idx)!=0)  Vertex2muTrkKFToSignalVertexChi2.at(t).Fill(Ntp->Vertex_2MuonsIsoTrack_KF_Chi2(final_idx)/Ntp->Vertex_signal_KF_Chi2(final_idx),1);
          Vertex2muTrkKFToSignalVertexDistance.at(t).Fill((Ntp->Vertex_Signal_KF_pos(final_idx) - Ntp->Vertex_2MuonsIsoTrack_KF_pos(final_idx)).Mag(),1);
        }



	deltaMuZ12.at(t).Fill(fabs(Ntp->Muon_Poca(Muon_index_1).Z()  - Ntp->Muon_Poca(Muon_index_2).Z()),w);
	deltaMuZ13.at(t).Fill(fabs(Ntp->Muon_Poca(Muon_index_1).Z()  - Ntp->Muon_Poca(Muon_index_3).Z()),w);
	deltaMuZ23.at(t).Fill(fabs(Ntp->Muon_Poca(Muon_index_2).Z()  - Ntp->Muon_Poca(Muon_index_3).Z()),w);
	
	
	MaxdeltaMuZ.at(t).Fill( std::max({fabs(Ntp->Muon_Poca(Muon_index_1).Z()  - Ntp->Muon_Poca(Muon_index_2).Z()),
		fabs(Ntp->Muon_Poca(Muon_index_1).Z()  - Ntp->Muon_Poca(Muon_index_3).Z()),
		fabs(Ntp->Muon_Poca(Muon_index_2).Z()  - Ntp->Muon_Poca(Muon_index_3).Z())}),w);
	
	MindeltaMuZ.at(t).Fill( std::min({fabs(Ntp->Muon_Poca(Muon_index_1).Z()  - Ntp->Muon_Poca(Muon_index_2).Z()),
		fabs(Ntp->Muon_Poca(Muon_index_1).Z()  - Ntp->Muon_Poca(Muon_index_3).Z()),
		fabs(Ntp->Muon_Poca(Muon_index_2).Z()  - Ntp->Muon_Poca(Muon_index_3).Z())}),w);
	
	MaxMuonsDca.at(t).Fill(std::max({Ntp->Vertex_DCA12(final_idx),Ntp->Vertex_DCA23(final_idx),Ntp->Vertex_DCA31(final_idx)}),w);
	MinMuonsDca.at(t).Fill(std::min({Ntp->Vertex_DCA12(final_idx),Ntp->Vertex_DCA23(final_idx),Ntp->Vertex_DCA31(final_idx)}),w);

	var_MaxdeltaMuZ = std::max({Ntp->Vertex_DCA12(final_idx),Ntp->Vertex_DCA23(final_idx),Ntp->Vertex_DCA31(final_idx)});
	var_MindeltaMuZ = std::min({Ntp->Vertex_DCA12(final_idx),Ntp->Vertex_DCA23(final_idx),Ntp->Vertex_DCA31(final_idx)});
	var_deltaMuZ12 = fabs(Ntp->Muon_Poca(Muon_index_1).Z()  - Ntp->Muon_Poca(Muon_index_2).Z());
	var_deltaMuZ13 = fabs(Ntp->Muon_Poca(Muon_index_1).Z()  - Ntp->Muon_Poca(Muon_index_3).Z());
	var_deltaMuZ23 = fabs(Ntp->Muon_Poca(Muon_index_2).Z()  - Ntp->Muon_Poca(Muon_index_3).Z());



	MinMatchedStations.at(t).Fill(std::min({Ntp->Muon_numberOfMatchedStations(Muon_index_1) ,Ntp->Muon_numberOfMatchedStations(Muon_index_2) , Ntp->Muon_numberOfMatchedStations(Muon_index_3)}),w);
	MaxMatchedStations.at(t).Fill(std::max({Ntp->Muon_numberOfMatchedStations(Muon_index_1) ,Ntp->Muon_numberOfMatchedStations(Muon_index_2) , Ntp->Muon_numberOfMatchedStations(Muon_index_3)}),w);
	Mu1MatchedStations.at(t).Fill(Ntp->Muon_numberOfMatchedStations(Muon_index_1),w);
	Mu2MatchedStations.at(t).Fill(Ntp->Muon_numberOfMatchedStations(Muon_index_2),w);
	Mu3MatchedStations.at(t).Fill(Ntp->Muon_numberOfMatchedStations(Muon_index_3),w);


	MinMuon_numberOfChambers.at(t).Fill(std::min({Ntp->Muon_numberOfChambers(Muon_index_1) ,Ntp->Muon_numberOfChambers(Muon_index_2) , Ntp->Muon_numberOfChambers(Muon_index_3)}),w);
	MaxMuon_numberOfChambers.at(t).Fill(std::max({Ntp->Muon_numberOfChambers(Muon_index_1) ,Ntp->Muon_numberOfChambers(Muon_index_2) , Ntp->Muon_numberOfChambers(Muon_index_3)}),w);
	Mu1Muon_numberOfChambers.at(t).Fill(Ntp->Muon_numberOfChambers(Muon_index_1),w);
	Mu2Muon_numberOfChambers.at(t).Fill(Ntp->Muon_numberOfChambers(Muon_index_2),w);
	Mu3Muon_numberOfChambers.at(t).Fill(Ntp->Muon_numberOfChambers(Muon_index_3),w);



	MinMuon_numberOfMatches.at(t).Fill(std::min({Ntp->Muon_numberOfMatches(Muon_index_1) ,Ntp->Muon_numberOfMatches(Muon_index_2) , Ntp->Muon_numberOfMatches(Muon_index_3)}),w);
	MaxMuon_numberOfMatches.at(t).Fill(std::max({Ntp->Muon_numberOfMatches(Muon_index_1) ,Ntp->Muon_numberOfMatches(Muon_index_2) , Ntp->Muon_numberOfMatches(Muon_index_3)}),w);
	Mu1Muon_numberOfMatches.at(t).Fill(Ntp->Muon_numberOfMatches(Muon_index_1),w);
	Mu2Muon_numberOfMatches.at(t).Fill(Ntp->Muon_numberOfMatches(Muon_index_2),w);
	Mu3Muon_numberOfMatches.at(t).Fill(Ntp->Muon_numberOfMatches(Muon_index_3),w);
	
	MinMuon_hitPattern_numberOfValidMuonHits.at(t).Fill(std::min({Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_1) ,Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_2) , Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_3)}),w);
	MaxMuon_hitPattern_numberOfValidMuonHits.at(t).Fill(std::max({Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_1) ,Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_2) , Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_3)}),w);
	Mu1Muon_hitPattern_numberOfValidMuonHits.at(t).Fill(Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_1),w);
	Mu2Muon_hitPattern_numberOfValidMuonHits.at(t).Fill(Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_2),w);
	Mu3Muon_hitPattern_numberOfValidMuonHits.at(t).Fill(Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_3),w);



	var_MinMatchedStations = std::min({Ntp->Muon_numberOfMatchedStations(Muon_index_1) ,Ntp->Muon_numberOfMatchedStations(Muon_index_2) , Ntp->Muon_numberOfMatchedStations(Muon_index_3)});
	var_MaxMatchedStations = std::max({Ntp->Muon_numberOfMatchedStations(Muon_index_1) ,Ntp->Muon_numberOfMatchedStations(Muon_index_2) , Ntp->Muon_numberOfMatchedStations(Muon_index_3)});
	var_Mu1MatchedStations = Ntp->Muon_numberOfMatchedStations(Muon_index_1);
	var_Mu2MatchedStations = Ntp->Muon_numberOfMatchedStations(Muon_index_2);
	var_Mu3MatchedStations = Ntp->Muon_numberOfMatchedStations(Muon_index_3);


	var_MinMuon_numberOfChambers = std::min({Ntp->Muon_numberOfChambers(Muon_index_1) ,Ntp->Muon_numberOfChambers(Muon_index_2) , Ntp->Muon_numberOfChambers(Muon_index_3)});
	var_MaxMuon_numberOfChambers = std::max({Ntp->Muon_numberOfChambers(Muon_index_1) ,Ntp->Muon_numberOfChambers(Muon_index_2) , Ntp->Muon_numberOfChambers(Muon_index_3)});
	var_Mu1Muon_numberOfChambers = Ntp->Muon_numberOfChambers(Muon_index_1);
	var_Mu2Muon_numberOfChambers = Ntp->Muon_numberOfChambers(Muon_index_2);
	var_Mu3Muon_numberOfChambers = Ntp->Muon_numberOfChambers(Muon_index_3);



	var_MinMuon_numberOfMatches = std::min({Ntp->Muon_numberOfMatches(Muon_index_1) ,Ntp->Muon_numberOfMatches(Muon_index_2) , Ntp->Muon_numberOfMatches(Muon_index_3)});
	var_MaxMuon_numberOfMatches = std::max({Ntp->Muon_numberOfMatches(Muon_index_1) ,Ntp->Muon_numberOfMatches(Muon_index_2) , Ntp->Muon_numberOfMatches(Muon_index_3)});
	var_Mu1Muon_numberOfMatches = Ntp->Muon_numberOfMatches(Muon_index_1);
	var_Mu2Muon_numberOfMatches = Ntp->Muon_numberOfMatches(Muon_index_2);
	var_Mu3Muon_numberOfMatches = Ntp->Muon_numberOfMatches(Muon_index_3);
	
	var_MinMuon_hitPattern_numberOfValidMuonHits = std::min({Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_1) ,Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_2) , Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_3)});
	var_MaxMuon_hitPattern_numberOfValidMuonHits = std::max({Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_1) ,Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_2) , Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_3)});
	var_Mu1Muon_hitPattern_numberOfValidMuonHits = Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_1);
	var_Mu2Muon_hitPattern_numberOfValidMuonHits = Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_2);
	var_Mu3Muon_hitPattern_numberOfValidMuonHits = Ntp->Muon_hitPattern_numberOfValidMuonHits(Muon_index_3);



	float MaxD0Significance = std::max({Ntp->Vertex_d0sig_reco(final_idx,0),
	      Ntp->Vertex_d0sig_reco(final_idx,1),
	      Ntp->Vertex_d0sig_reco(final_idx,2)});
	
	
	float MinD0Significance = std::min({Ntp->Vertex_d0sig_reco(final_idx,0),
	      Ntp->Vertex_d0sig_reco(final_idx,1),
	      Ntp->Vertex_d0sig_reco(final_idx,2)});
	

	float MaxD0BSSignificance = std::max({Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,0),
	      Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,1),
	      Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,2)});
	
	
	float MinD0BSSignificance = std::min({Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,0),
	      Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,1),
	      Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,2)});
	

	
	float MaxD0SVSignificance = std::max({Ntp->Vertex_d0sigSV_reco(final_idx,0),
	      Ntp->Vertex_d0sigSV_reco(final_idx,1),
	      Ntp->Vertex_d0sigSV_reco(final_idx,2)});
	
	
	float MinD0SVSignificance = std::min({Ntp->Vertex_d0sigSV_reco(final_idx,0),
	      Ntp->Vertex_d0sigSV_reco(final_idx,1),
	      Ntp->Vertex_d0sigSV_reco(final_idx,2)});
	

	VertexMu1D0SigPVReco.at(t).Fill(Ntp->Vertex_d0sig_reco(final_idx,0),w);
	VertexMu2D0SigPVReco.at(t).Fill(Ntp->Vertex_d0sig_reco(final_idx,1),w);
	VertexMu3D0SigPVReco.at(t).Fill(Ntp->Vertex_d0sig_reco(final_idx,2),w);
	
	MaxD0SigPV.at(t).Fill(MaxD0Significance,w);
	MinD0SigPV.at(t).Fill(MinD0Significance,w);



	VertexMu1D0SigBSReco.at(t).Fill(Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,0),w);
	VertexMu2D0SigBSReco.at(t).Fill(Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,1),w);
	VertexMu3D0SigBSReco.at(t).Fill(Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,2),w);

	MaxD0SigBS.at(t).Fill(MaxD0BSSignificance,w);
	MinD0SigBS.at(t).Fill(MinD0BSSignificance,w);
	

	VertexMu1D0SigSVReco.at(t).Fill(Ntp->Vertex_d0sigSV_reco(final_idx,0),w);
	VertexMu2D0SigSVReco.at(t).Fill(Ntp->Vertex_d0sigSV_reco(final_idx,1),w);
	VertexMu3D0SigSVReco.at(t).Fill(Ntp->Vertex_d0sigSV_reco(final_idx,2),w);



	MaxD0SigSV.at(t).Fill(MaxD0SVSignificance,w);
	MinD0SigSV.at(t).Fill(MinD0SVSignificance,w);




	VertexMu1D0DistSVReco.at(t).Fill(Ntp->Vertex_d0SV_reco(final_idx,0),w);
	VertexMu2D0DistSVReco.at(t).Fill(Ntp->Vertex_d0SV_reco(final_idx,1),w);
	VertexMu3D0DistSVReco.at(t).Fill(Ntp->Vertex_d0SV_reco(final_idx,2),w);
	MaxD0DistSV.at(t).Fill(std::max({Ntp->Vertex_d0SV_reco(final_idx,0),
		Ntp->Vertex_d0SV_reco(final_idx,1),
		Ntp->Vertex_d0SV_reco(final_idx,2)}),w);
	MinD0DistSV.at(t).Fill(std::min({Ntp->Vertex_d0SV_reco(final_idx,0),
                Ntp->Vertex_d0SV_reco(final_idx,1),
                Ntp->Vertex_d0SV_reco(final_idx,2)}),w);

	
	VertexMu1DZDistSVReco.at(t).Fill(Ntp->Vertex_dzSV_reco(final_idx,0),w);
	VertexMu2DZDistSVReco.at(t).Fill(Ntp->Vertex_dzSV_reco(final_idx,1),w);
	VertexMu3DZDistSVReco.at(t).Fill(Ntp->Vertex_dzSV_reco(final_idx,2),w);
	
	
	MaxDZDistSV.at(t).Fill(std::max({Ntp->Vertex_dzSV_reco(final_idx,0),
		Ntp->Vertex_dzSV_reco(final_idx,1),
		Ntp->Vertex_dzSV_reco(final_idx,2)}),w);
	MinDZDistSV.at(t).Fill(std::min({Ntp->Vertex_dzSV_reco(final_idx,0),
		Ntp->Vertex_dzSV_reco(final_idx,1),
		Ntp->Vertex_dzSV_reco(final_idx,2)}),w);



	VertexMu1DistanceToSV.at(t).Fill(sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,0),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,0),2) ) ,w);
	VertexMu2DistanceToSV.at(t).Fill(sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,1),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,1),2) ) ,w);
	VertexMu3DistanceToSV.at(t).Fill(sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,2),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,2),2) ) ,w);

	MaxMuDistanceToSV.at(t).Fill(std::max({sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,0),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,0),2)) ,
		sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,1),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,1),2) ),
		sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,2),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,2),2))}),w);
	MinMuDistanceToSV.at(t).Fill(std::min({sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,0),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,0),2)) ,
                sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,1),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,1),2) ),
                sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,2),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,2),2))}),w);
	  
	var_VertexMu1DistanceToSV=sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,0),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,0),2) );
	var_VertexMu2DistanceToSV=sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,1),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,1),2) );
	var_VertexMu3DistanceToSV=sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,2),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,2),2) );
	
	var_MaxMuDistanceToSV=std::max({sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,0),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,0),2)) ,
	      sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,1),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,1),2) ),
	      sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,2),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,2),2))});
	var_MinMuDistanceToSV=std::min({sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,0),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,0),2)) ,
	      sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,1),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,1),2) ),
	      sqrt(pow(Ntp->Vertex_dzSV_reco(final_idx,2),2)   +  pow(Ntp->Vertex_d0SV_reco(final_idx,2),2))});





	var_VertexMu1D0DistSVReco=Ntp->Vertex_d0SV_reco(final_idx,0);
	var_VertexMu2D0DistSVReco=Ntp->Vertex_d0SV_reco(final_idx,1);
	var_VertexMu3D0DistSVReco=Ntp->Vertex_d0SV_reco(final_idx,2);
	
	var_MaxD0DistSV=std::max({Ntp->Vertex_d0SV_reco(final_idx,0),
	      Ntp->Vertex_d0SV_reco(final_idx,1),
	      Ntp->Vertex_d0SV_reco(final_idx,2)});
	
	var_MinD0DistSV=std::min({Ntp->Vertex_d0SV_reco(final_idx,0),
	      Ntp->Vertex_d0SV_reco(final_idx,1),
	      Ntp->Vertex_d0SV_reco(final_idx,2)});
	
	
	var_VertexMu1DZDistSVReco=Ntp->Vertex_dzSV_reco(final_idx,0);
	var_VertexMu2DZDistSVReco=Ntp->Vertex_dzSV_reco(final_idx,1);
	var_VertexMu3DZDistSVReco=Ntp->Vertex_dzSV_reco(final_idx,2);
	
	var_MaxDZDistSV=std::max({Ntp->Vertex_dzSV_reco(final_idx,0),
	      Ntp->Vertex_dzSV_reco(final_idx,1),
	      Ntp->Vertex_dzSV_reco(final_idx,2)});
	var_MinDZDistSV=std::min({Ntp->Vertex_dzSV_reco(final_idx,0),
	      Ntp->Vertex_dzSV_reco(final_idx,1),
	      Ntp->Vertex_dzSV_reco(final_idx,2)});
	


	

	MinMuon_chi2LocalPosition.at(t).Fill(std::min({Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3)  }));
	MaxMuon_chi2LocalPosition.at(t).Fill(std::max({Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3)  }));


	Muon1_chi2LocalPosition.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1),w);
	Muon2_chi2LocalPosition.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2),w);
	Muon3_chi2LocalPosition.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3),w);
	
	
	MinMuon_chi2LocalMomentum.at(t).Fill(std::min({Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_3)  }));
	MaxMuon_chi2LocalMomentum.at(t).Fill(std::max({Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_3)  }));


	Muon1_chi2LocalMomentum.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_1), 1);
	Muon2_chi2LocalMomentum.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_2), 1);
	Muon3_chi2LocalMomentum.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_3), 1);
	
	
	
	MintrkKink.at(t).Fill(std::min({Ntp->Muon_combinedQuality_trkKink(Muon_index_1),Ntp->Muon_combinedQuality_trkKink(Muon_index_2),Ntp->Muon_combinedQuality_trkKink(Muon_index_3)}),w);
	MaxtrkKink.at(t).Fill(std::max({Ntp->Muon_combinedQuality_trkKink(Muon_index_1),Ntp->Muon_combinedQuality_trkKink(Muon_index_2),Ntp->Muon_combinedQuality_trkKink(Muon_index_3)}),w);
	MinglbKink.at(t).Fill(std::min({Ntp->Muon_combinedQuality_glbKink(Muon_index_1),Ntp->Muon_combinedQuality_glbKink(Muon_index_2),Ntp->Muon_combinedQuality_glbKink(Muon_index_3)}),w);
	MaxglbKink.at(t).Fill(log(std::max({Ntp->Muon_combinedQuality_glbKink(Muon_index_1),Ntp->Muon_combinedQuality_glbKink(Muon_index_2),Ntp->Muon_combinedQuality_glbKink(Muon_index_3)})),w);
	
	MuonglbkinkSum.at(t).Fill(var_sumMuTrkKinkChi2,w);

	MaxVertexPairQuality.at(t).Fill(  std::max({Ntp->Vertex_pair_quality(final_idx,0),Ntp->Vertex_pair_quality(final_idx,1),Ntp->Vertex_pair_quality(final_idx,2)}),w);
	MinVertexPairQuality.at(t).Fill(  std::min({Ntp->Vertex_pair_quality(final_idx,0),Ntp->Vertex_pair_quality(final_idx,1),Ntp->Vertex_pair_quality(final_idx,2)}),w);
	
	
	//  ----------------------------- secondary vertices ----------------
	int NumberOfPrimaryVertices(0);
	for(unsigned int iVertex=0; iVertex < Ntp->NSecondaryVertices(); iVertex++){
	  SV_Mass.at(t).Fill(Ntp->SecondaryVertexMass(iVertex),2);
	  TVector3 SVsignalPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(final_idx),Ntp->Vertex_MatchedPrimaryVertex(final_idx));
	  TVector3 SVfakePV = Ntp->SVPVDirection(Ntp->SecondaryVertexPosition(iVertex),Ntp->Vertex_MatchedPrimaryVertex(final_idx));
	  
	  SVDeltaR.at(t).Fill(SVfakePV.DeltaR(SVsignalPV),w);
	  SVDistance.at(t).Fill((Ntp->Vertex_Signal_KF_pos(final_idx) - Ntp->SecondaryVertexPosition(iVertex)).Mag(),w);
	  
	  
	  if(SVfakePV.DeltaR(SVsignalPV) < 0.3 && (Ntp->Vertex_Signal_KF_pos(final_idx) - Ntp->SecondaryVertexPosition(iVertex)).Mag() > 0.015){ // sv in the tau cone but  displaced
	    NumberOfPrimaryVertices++;
	    
	  }
	}
	NSV.at(t).Fill(NumberOfPrimaryVertices,w);


	// -----------------------------------------------------------------




	// --------------------------- Isolation --------------------------------

	int NcloseTracksCount(0);
	int NPVcloseTracksCount(0);
	int NPVcloseDRToTauTracksCount(0);
	double SumPT02(0),SumPT04(0),SumPT06(0),SumPT08(0),SumPT1(0),SumPT12(0),SumPT14(0),SumPT16(0),SumPT18(0),SumPT2(0);
	double Iso08Muon1,Iso08Muon2,Iso08Muon3;
	double Iso05Muon1,Iso05Muon2,Iso05Muon3;
	int TrackIndex(0);
	int TrackIndex_closestToPV(0);
	double TrackPTtreschold(0.4);
	double dca_temp(999.);
	double dcaPV_temp(999.);

	double Chi2IsoTrackVertexToMuon3(999.);
	TLorentzVector IsoTrack_P4_closestToMu3(0,0,0,0);

	double Chi2IsoTrackVertexToMuon2(999.);
	TLorentzVector IsoTrack_P4_closestToMu2(0,0,0,0);

	double Chi2IsoTrackVertexToMuon1(999.);
	TLorentzVector IsoTrack_P4_closestToMu1(0,0,0,0);

	for(int i =0; i< Ntp->NIsolationTrack(final_idx); i++){
	  
	  for(int ii=i+1 ; ii < Ntp->NIsolationTrack(final_idx); ii++){

	    if(Ntp->IsolationTrack_charge(final_idx,i)*Ntp->IsolationTrack_charge(final_idx,ii)==-1){
	      IsolationCombinatorialMass_pipi.at(t).Fill( (Ntp->IsolationTrack_p4(final_idx,i) +Ntp->IsolationTrack_p4(final_idx,ii)).M(),1  );
	    }

	  }



	   if(Ntp->Muon_charge(Ntp->ThreeMuonIndices(final_idx).at(2)) * Ntp->IsolationTrack_charge(final_idx,i) == -1);
	   {
	     if(Ntp->IsolationTrack_VertexWithSignalMuon3Chi2(final_idx,i) < Chi2IsoTrackVertexToMuon3)
	       {
		 Chi2IsoTrackVertexToMuon3= Ntp->IsolationTrack_VertexWithSignalMuon3Chi2(final_idx,i);
		 IsoTrack_P4_closestToMu3 = Ntp->IsolationTrack_p4(final_idx,i);
	       }
	   }


	   if(Ntp->Muon_charge(Ntp->ThreeMuonIndices(final_idx).at(1)) * Ntp->IsolationTrack_charge(final_idx,i) == -1);
	   {
	     if(Ntp->IsolationTrack_VertexWithSignalMuon2Chi2(final_idx,i) < Chi2IsoTrackVertexToMuon2)
	       {
		 Chi2IsoTrackVertexToMuon2= Ntp->IsolationTrack_VertexWithSignalMuon2Chi2(final_idx,i);
		 IsoTrack_P4_closestToMu2 = Ntp->IsolationTrack_p4(final_idx,i);

	       }

	   }


	   if(Ntp->Muon_charge(Ntp->ThreeMuonIndices(final_idx).at(0)) * Ntp->IsolationTrack_charge(final_idx,i) == -1);
	   {
	     if(Ntp->IsolationTrack_VertexWithSignalMuon1Chi2(final_idx,i) < Chi2IsoTrackVertexToMuon1)
	       {
		 Chi2IsoTrackVertexToMuon1= Ntp->IsolationTrack_VertexWithSignalMuon1Chi2(final_idx,i);
		 IsoTrack_P4_closestToMu1 = Ntp->IsolationTrack_p4(final_idx,i);

	       }

	   }



	  
	  if(Ntp->IsolationTrack_p4(final_idx,i).Pt()> 0.5  && sqrt(  pow(Ntp->IsolationTrack_dzSV(final_idx,i),2)   +   
								      pow(Ntp->IsolationTrack_dxySV(final_idx,i),2)) < 0.03)
	    {
	      NcloseTracksCount++;
	    }


	  if(Ntp->IsolationTrack_p4(final_idx,i).Pt()> 0.5  && sqrt(  pow(Ntp->IsolationTrack_dzPV(final_idx,i),2)   +   
								      pow(Ntp->IsolationTrack_dxyPV(final_idx,i),2)) < 0.05)
	    {
	      NPVcloseTracksCount++;
	    }


	  if(Ntp->IsolationTrack_p4(final_idx,i).Pt()> 0.5  && sqrt(  pow(Ntp->IsolationTrack_dzPV(final_idx,i),2)   +   
								      pow(Ntp->IsolationTrack_dxyPV(final_idx,i),2)) < 0.03  && Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauLV) < 0.8)
	    {
	      NPVcloseDRToTauTracksCount++;
	    }





	  if( sqrt(  pow(Ntp->IsolationTrack_dzSV(final_idx,i),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,i),2) ) <  dca_temp){
	    dca_temp = sqrt(  pow(Ntp->IsolationTrack_dzSV(final_idx,i),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,i),2));
	    TrackIndex = i;
	  }

	  if(Ntp->IsolationTrack_p4(final_idx,i).Pt()> TrackPTtreschold && fabs(Ntp->IsolationTrack_dzPV(final_idx,i)) < 0.05 && 
	     sqrt(  pow(Ntp->IsolationTrack_dzSV(final_idx,i),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,i),2)) < 0.05){


	    if( sqrt(  pow(Ntp->IsolationTrack_dzPV(final_idx,i),2)   +   pow(Ntp->IsolationTrack_dxyPV(final_idx,i),2) ) <  dcaPV_temp){
	      dcaPV_temp = sqrt(  pow(Ntp->IsolationTrack_dzPV(final_idx,i),2)   +   pow(Ntp->IsolationTrack_dxyPV(final_idx,i),2));
	      TrackIndex_closestToPV = i;
	    }


	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 0.2){
	      SumPT02 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 0.4){
	      SumPT04 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 0.6){
	      SumPT06 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 0.8){
	      SumPT08 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 1.0){
	      SumPT1 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 1.2){
	      SumPT12 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 1.4){
	      SumPT14 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }  
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 1.6){
	      SumPT16 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 1.8){
	      SumPT18 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(TauRefitLV) < 2.0){
	      SumPT2 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	  }
	}

	IsoTrack_P4_closestToMu3.SetE(sqrt(IsoTrack_P4_closestToMu3.Px()*IsoTrack_P4_closestToMu3.Px() +
					   IsoTrack_P4_closestToMu3.Py()*IsoTrack_P4_closestToMu3.Py() +
					   IsoTrack_P4_closestToMu3.Pz()*IsoTrack_P4_closestToMu3.Pz() + 0.493677*0.493677));

	TLorentzVector IsoTrack_P4_closestToMu3_PiMass = IsoTrack_P4_closestToMu3;
	IsoTrack_P4_closestToMu3_PiMass.SetE(sqrt(IsoTrack_P4_closestToMu3.Px()*IsoTrack_P4_closestToMu3.Px() +
						  IsoTrack_P4_closestToMu3.Py()*IsoTrack_P4_closestToMu3.Py() +
						  IsoTrack_P4_closestToMu3.Pz()*IsoTrack_P4_closestToMu3.Pz() + 0.135*0.135));

	TLorentzVector IsoTrack_P4_closestToMu3_MuMass = IsoTrack_P4_closestToMu3;
	IsoTrack_P4_closestToMu3_MuMass.SetE(sqrt(IsoTrack_P4_closestToMu3.Px()*IsoTrack_P4_closestToMu3.Px() +
						  IsoTrack_P4_closestToMu3.Py()*IsoTrack_P4_closestToMu3.Py() +
						  IsoTrack_P4_closestToMu3.Pz()*IsoTrack_P4_closestToMu3.Pz() + 0.106*0.106));


	TLorentzVector Mu3_WithKMass = Ntp->Muon_P4(Ntp->ThreeMuonIndices(final_idx).at(2));

	Mu3_WithKMass.SetE(sqrt(Mu3_WithKMass.Px()*Mu3_WithKMass.Px() + 
				Mu3_WithKMass.Py()*Mu3_WithKMass.Py() + 
				Mu3_WithKMass.Pz()*Mu3_WithKMass.Pz() + 0.493677*0.493677));


	IsoPhiKKMass_Mu3.at(t).Fill((IsoTrack_P4_closestToMu3+Mu3_WithKMass).M(),1);
	IsoKStarMass_Mu3.at(t).Fill((IsoTrack_P4_closestToMu3_PiMass+Mu3_WithKMass).M(),1);
	IsoMuMuMass_Mu3.at(t).Fill((IsoTrack_P4_closestToMu3_MuMass+ Ntp->Muon_P4(Ntp->ThreeMuonIndices(final_idx).at(2))).M(),1);
	IsoPhiKKMass_Mu3_wideRange.at(t).Fill((IsoTrack_P4_closestToMu3+Mu3_WithKMass).M(),1);
	IsoPhiKKMass_Mu3_midRange.at(t).Fill((IsoTrack_P4_closestToMu3+Mu3_WithKMass).M(),1);
	IsoKStarMass_Mu3_wideRange.at(t).Fill((IsoTrack_P4_closestToMu3_PiMass+Mu3_WithKMass).M(),1);
	IsoMuMuMass_Mu3_wideRange.at(t).Fill((IsoTrack_P4_closestToMu3_MuMass+ Ntp->Muon_P4(Ntp->ThreeMuonIndices(final_idx).at(2))).M(),1);


	var_IsoPhiKKMass_Mu3 = (IsoTrack_P4_closestToMu3+Mu3_WithKMass).M();
	var_IsoKStarMass_Mu3 = (IsoTrack_P4_closestToMu3_PiMass+Mu3_WithKMass).M();
	var_IsoMuMuMass_Mu3 = (IsoTrack_P4_closestToMu3_MuMass+ Ntp->Muon_P4(Ntp->ThreeMuonIndices(final_idx).at(2))).M();




	IsoTrack_P4_closestToMu2.SetE(sqrt(IsoTrack_P4_closestToMu2.Px()*IsoTrack_P4_closestToMu2.Px() +
					   IsoTrack_P4_closestToMu2.Py()*IsoTrack_P4_closestToMu2.Py() +
					   IsoTrack_P4_closestToMu2.Pz()*IsoTrack_P4_closestToMu2.Pz() + 0.493677*0.493677));
	TLorentzVector Mu2_WithKMass = Ntp->Muon_P4(Ntp->ThreeMuonIndices(final_idx).at(1));

	Mu2_WithKMass.SetE(sqrt(Mu2_WithKMass.Px()*Mu2_WithKMass.Px() + 
				Mu2_WithKMass.Py()*Mu2_WithKMass.Py() + 
				Mu2_WithKMass.Pz()*Mu2_WithKMass.Pz() + 0.493677*0.493677));

	TLorentzVector IsoTrack_P4_closestToMu2_PiMass = IsoTrack_P4_closestToMu2;
	IsoTrack_P4_closestToMu2_PiMass.SetE(sqrt(IsoTrack_P4_closestToMu2.Px()*IsoTrack_P4_closestToMu2.Px() +
						  IsoTrack_P4_closestToMu2.Py()*IsoTrack_P4_closestToMu2.Py() +
						  IsoTrack_P4_closestToMu2.Pz()*IsoTrack_P4_closestToMu2.Pz() + 0.135*0.135));

	TLorentzVector IsoTrack_P4_closestToMu2_MuMass = IsoTrack_P4_closestToMu2;
	IsoTrack_P4_closestToMu2_MuMass.SetE(sqrt(IsoTrack_P4_closestToMu2.Px()*IsoTrack_P4_closestToMu2.Px() +
						  IsoTrack_P4_closestToMu2.Py()*IsoTrack_P4_closestToMu2.Py() +
						  IsoTrack_P4_closestToMu2.Pz()*IsoTrack_P4_closestToMu2.Pz() + 0.106*0.106));



	IsoPhiKKMass_Mu2.at(t).Fill((IsoTrack_P4_closestToMu2+Mu2_WithKMass).M(),1);
	IsoKStarMass_Mu2.at(t).Fill((IsoTrack_P4_closestToMu2_PiMass+Mu2_WithKMass).M(),1);
	IsoMuMuMass_Mu2.at(t).Fill((IsoTrack_P4_closestToMu2_MuMass+ Ntp->Muon_P4(Ntp->ThreeMuonIndices(final_idx).at(1))).M(),1);
	IsoPhiKKMass_Mu2_wideRange.at(t).Fill((IsoTrack_P4_closestToMu2+Mu2_WithKMass).M(),1);
	IsoPhiKKMass_Mu2_midRange.at(t).Fill((IsoTrack_P4_closestToMu2+Mu2_WithKMass).M(),1);
	IsoKStarMass_Mu2_wideRange.at(t).Fill((IsoTrack_P4_closestToMu2_PiMass+Mu2_WithKMass).M(),1);
	IsoMuMuMass_Mu2_wideRange.at(t).Fill((IsoTrack_P4_closestToMu2_MuMass+ Ntp->Muon_P4(Ntp->ThreeMuonIndices(final_idx).at(1))).M(),1);


	var_IsoPhiKKMass_Mu2 = (IsoTrack_P4_closestToMu2+Mu2_WithKMass).M();
	var_IsoKStarMass_Mu2 = (IsoTrack_P4_closestToMu2_PiMass+Mu2_WithKMass).M();
	var_IsoMuMuMass_Mu2 = (IsoTrack_P4_closestToMu2_MuMass+ Ntp->Muon_P4(Ntp->ThreeMuonIndices(final_idx).at(1))).M();







	IsoTrack_P4_closestToMu1.SetE(sqrt(IsoTrack_P4_closestToMu1.Px()*IsoTrack_P4_closestToMu1.Px() +
					   IsoTrack_P4_closestToMu1.Py()*IsoTrack_P4_closestToMu1.Py() +
					   IsoTrack_P4_closestToMu1.Pz()*IsoTrack_P4_closestToMu1.Pz() + 0.493677*0.493677));

	TLorentzVector IsoTrack_P4_closestToMu1_PiMass = IsoTrack_P4_closestToMu1;
	IsoTrack_P4_closestToMu1_PiMass.SetE(sqrt(IsoTrack_P4_closestToMu1.Px()*IsoTrack_P4_closestToMu1.Px() +
						  IsoTrack_P4_closestToMu1.Py()*IsoTrack_P4_closestToMu1.Py() +
						  IsoTrack_P4_closestToMu1.Pz()*IsoTrack_P4_closestToMu1.Pz() + 0.135*0.135));

	TLorentzVector IsoTrack_P4_closestToMu1_MuMass = IsoTrack_P4_closestToMu1;
	IsoTrack_P4_closestToMu1_MuMass.SetE(sqrt(IsoTrack_P4_closestToMu1.Px()*IsoTrack_P4_closestToMu1.Px() +
						  IsoTrack_P4_closestToMu1.Py()*IsoTrack_P4_closestToMu1.Py() +
						  IsoTrack_P4_closestToMu1.Pz()*IsoTrack_P4_closestToMu1.Pz() + 0.106*0.106));



	TLorentzVector Mu1_WithKMass = Ntp->Muon_P4(Ntp->ThreeMuonIndices(final_idx).at(0));

	Mu1_WithKMass.SetE(sqrt(Mu1_WithKMass.Px()*Mu1_WithKMass.Px() + 
				Mu1_WithKMass.Py()*Mu1_WithKMass.Py() + 
				Mu1_WithKMass.Pz()*Mu1_WithKMass.Pz() + 0.493677*0.493677));


	IsoPhiKKMass_Mu1.at(t).Fill((IsoTrack_P4_closestToMu1+Mu1_WithKMass).M(),1);
	IsoKStarMass_Mu1.at(t).Fill((IsoTrack_P4_closestToMu1_PiMass+Mu1_WithKMass).M(),1);
	IsoMuMuMass_Mu1.at(t).Fill((IsoTrack_P4_closestToMu1_MuMass+ Ntp->Muon_P4(Ntp->ThreeMuonIndices(final_idx).at(0))).M(),1);
	IsoPhiKKMass_Mu1_wideRange.at(t).Fill((IsoTrack_P4_closestToMu1+Mu1_WithKMass).M(),1);
	IsoPhiKKMass_Mu1_midRange.at(t).Fill((IsoTrack_P4_closestToMu1+Mu1_WithKMass).M(),1);
	IsoKStarMass_Mu1_wideRange.at(t).Fill((IsoTrack_P4_closestToMu1_PiMass+Mu1_WithKMass).M(),1);
	IsoMuMuMass_Mu1_wideRange.at(t).Fill((IsoTrack_P4_closestToMu1_MuMass+ Ntp->Muon_P4(Ntp->ThreeMuonIndices(final_idx).at(0))).M(),1);


	var_IsoPhiKKMass_Mu1 = (IsoTrack_P4_closestToMu1+Mu1_WithKMass).M();
	var_IsoKStarMass_Mu1 = (IsoTrack_P4_closestToMu1_PiMass+Mu1_WithKMass).M();
	var_IsoMuMuMass_Mu1 = (IsoTrack_P4_closestToMu1_MuMass+ Ntp->Muon_P4(Ntp->ThreeMuonIndices(final_idx).at(0))).M();






	NtracksClose.at(t).Fill(NcloseTracksCount,w);
	NTracksCloseToPV.at(t).Fill(NPVcloseTracksCount,w);
	NTracksCloseToPVTauDR.at(t).Fill(NPVcloseDRToTauTracksCount,w);

	


	Iso02.at(t).Fill(TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT02),w);
	Iso04.at(t).Fill(TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT04),w);
	Iso06.at(t).Fill(TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT06),w);
	Iso08.at(t).Fill(TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT08),w);
	Iso1.at(t).Fill(TauRefitLV.Pt()/  (TauRefitLV.Pt() + SumPT1),w);
	Iso12.at(t).Fill(TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT12),w);
	Iso14.at(t).Fill(TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT14),w);
	Iso16.at(t).Fill(TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT16),w);
	Iso18.at(t).Fill(TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT18),w);
	Iso2.at(t).Fill(TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT2),w);


	var_Iso02 = TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT02);
	var_Iso04 = TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT04);
	var_Iso06 = TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT06);
	var_Iso08 = TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT08);
	var_Iso1 = TauRefitLV.Pt()/  (TauRefitLV.Pt() + SumPT1);
	var_Iso12 = TauRefitLV.Pt()/ (TauRefitLV.Pt() + SumPT12);



	if(Ntp->NIsolationTrack(final_idx)!=0){ 
	  MindcaTrackSV.at(t).Fill(sqrt( pow(Ntp->IsolationTrack_dzSV(final_idx,TrackIndex),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,TrackIndex),2)),w);
	}

	if(Ntp->NIsolationTrack(final_idx)!=0){ 
	  dcaTrackPV.at(t).Fill(sqrt(  pow(Ntp->IsolationTrack_dzPV(final_idx,TrackIndex_closestToPV),2)   +   pow(Ntp->IsolationTrack_dxyPV(final_idx,TrackIndex_closestToPV),2)),1 );
	}


	// ---------------------- I_mu1
	SumPT02=0;SumPT04=0;SumPT06=0;SumPT08=0;SumPT1=0;SumPT12=0;SumPT14=0;SumPT16=0;SumPT18=0;SumPT2=0;
	for(int i =0; i< Ntp->NIsolationTrack(final_idx); i++){

	  if(Ntp->IsolationTrack_p4(final_idx,i).Pt()> TrackPTtreschold && fabs(Ntp->IsolationTrack_dzPV(final_idx,i)) < 0.05
	     && Ntp->IsolationTrack_DocaMu1(final_idx,i) < 0.1){
	    if(Ntp->Muon_charge(Ntp->ThreeMuonIndices(final_idx).at(0))* Ntp->IsolationTrack_charge(final_idx,i)==-1 ){
	      Mu1TrackMass.at(t).Fill(  (Muon1LV +Ntp->IsolationTrack_p4(final_idx,i)).M(),1 );
	      var_Mu1TrackMass = (Muon1LV +Ntp->IsolationTrack_p4(final_idx,i)).M();
	    }else var_Mu1TrackMass = -1;

	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 0.2){
	      SumPT02 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 0.4){
	      SumPT04 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 0.6){
	      SumPT06 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 0.8){
	      SumPT08 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 1.0){
	      SumPT1 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 1.2){
	      SumPT12 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 1.4){
	      SumPT14 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 1.6){
	      SumPT16 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 1.8){
	      SumPT18 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon1LV) < 2.0){
	      SumPT2 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	  }
	}

	Iso02Mu1.at(t).Fill(Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT02),w);
	Iso04Mu1.at(t).Fill(Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT04),w);
	Iso06Mu1.at(t).Fill(Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT06),w);
	Iso08Mu1.at(t).Fill(Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT08),w); Iso08Muon1 = Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT08);
	Iso1Mu1.at(t).Fill(Muon1LV.Pt()/  (Muon1LV.Pt() + SumPT1),w);
	Iso12Mu1.at(t).Fill(Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT12),w);
	Iso14Mu1.at(t).Fill(Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT14),w);
	Iso16Mu1.at(t).Fill(Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT16),w);
	Iso18Mu1.at(t).Fill(Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT18),w);
	Iso2Mu1.at(t).Fill(Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT2),w);


	var_Iso02Mu1 = Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT02);
	var_Iso04Mu1 = Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT04);
	var_Iso06Mu1 = Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT06);
	var_Iso08Mu1 = Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT08);
	var_Iso1Mu1 = Muon1LV.Pt()/  (Muon1LV.Pt() + SumPT1);
	var_Iso12Mu1 = Muon1LV.Pt()/ (Muon1LV.Pt() + SumPT12);



	// ---------------------- I_mu2

	SumPT02=0;SumPT04=0;SumPT06=0;SumPT08=0;SumPT1=0;SumPT12=0;SumPT14=0;SumPT16=0;SumPT18=0;SumPT2=0;
	for(int i =0; i< Ntp->NIsolationTrack(final_idx); i++){

	  if(Ntp->IsolationTrack_p4(final_idx,i).Pt()> TrackPTtreschold && fabs(Ntp->IsolationTrack_dzPV(final_idx,i)) < 0.05
	     && Ntp->IsolationTrack_DocaMu2(final_idx,i) < 0.1){


	    if(Ntp->Muon_charge(Ntp->ThreeMuonIndices(final_idx).at(1))* Ntp->IsolationTrack_charge(final_idx,i)==-1 ){

	      Mu2TrackMass.at(t).Fill(  (Muon2LV +Ntp->IsolationTrack_p4(final_idx,i)).M(),1 );
	      var_Mu2TrackMass = (Muon2LV +Ntp->IsolationTrack_p4(final_idx,i)).M();
	    }else var_Mu2TrackMass = -1;


	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 0.2){
	      SumPT02 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 0.4){
	      SumPT04 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 0.6){
	      SumPT06 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 0.8){
	      SumPT08 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 1.0){
	      SumPT1 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 1.2){
	      SumPT12 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 1.4){
	      SumPT14 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 1.6){
	      SumPT16 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 1.8){
	      SumPT18 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon2LV) < 2.0){
	      SumPT2 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	  }
	}

	Iso02Mu2.at(t).Fill(Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT02),1);
	Iso04Mu2.at(t).Fill(Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT04),1);
	Iso06Mu2.at(t).Fill(Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT06),1);
	Iso08Mu2.at(t).Fill(Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT08),1);Iso08Muon2 = Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT08);
	Iso1Mu2.at(t).Fill(Muon2LV.Pt()/  (Muon2LV.Pt() + SumPT1),1);
	Iso12Mu2.at(t).Fill(Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT12),1);
	Iso14Mu2.at(t).Fill(Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT14),1);
	Iso16Mu2.at(t).Fill(Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT16),1);
	Iso18Mu2.at(t).Fill(Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT18),1);
	Iso2Mu2.at(t).Fill(Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT2),1);
    
	var_Iso02Mu2 = Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT02);
	var_Iso04Mu2 = Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT04);
	var_Iso06Mu2 = Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT06);
	var_Iso08Mu2 = Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT08);
	var_Iso1Mu2 = Muon2LV.Pt()/  (Muon2LV.Pt() + SumPT1);
	var_Iso12Mu2 = Muon2LV.Pt()/ (Muon2LV.Pt() + SumPT12);



	// ---------------------- I_mu3

	SumPT02=0;SumPT04=0;SumPT06=0;SumPT08=0;SumPT1=0;SumPT12=0;SumPT14=0;SumPT16=0;SumPT18=0;SumPT2=0;
	for(int i =0; i< Ntp->NIsolationTrack(final_idx); i++){

	  if(Ntp->IsolationTrack_p4(final_idx,i).Pt()> TrackPTtreschold && fabs(Ntp->IsolationTrack_dzPV(final_idx,i)) < 0.05
	     && Ntp->IsolationTrack_DocaMu3(final_idx,i) < 0.1){


	    if(Ntp->Muon_charge(Ntp->ThreeMuonIndices(final_idx).at(2))* Ntp->IsolationTrack_charge(final_idx,i)==-1 ){
	      Mu3TrackMass.at(t).Fill(  (Muon3LV +Ntp->IsolationTrack_p4(final_idx,i)).M(),1 );
	      var_Mu3TrackMass = (Muon3LV +Ntp->IsolationTrack_p4(final_idx,i)).M();
	    }else var_Mu3TrackMass = -1;


	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 0.2){
	      SumPT02 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 0.4){
	      SumPT04 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 0.6){
	      SumPT06 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 0.8){
	      SumPT08 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 1.0){
	      SumPT1 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 1.2){
	      SumPT12 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 1.4){
	      SumPT14 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 1.6){
	      SumPT16 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 1.8){
	      SumPT18 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	    if(Ntp->IsolationTrack_p4(final_idx,i).DeltaR(Muon3LV) < 2.0){
	      SumPT2 += Ntp->IsolationTrack_p4(final_idx,i).Pt();
	    }
	  }
	}
    
	Iso02Mu3.at(t).Fill(Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT02),1);    
	Iso04Mu3.at(t).Fill(Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT04),1);
	Iso06Mu3.at(t).Fill(Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT06),1);
	Iso08Mu3.at(t).Fill(Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT08),1);Iso08Muon3 = Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT08);
	Iso1Mu3.at(t).Fill(Muon3LV.Pt()/  (Muon3LV.Pt() + SumPT1),1);
	Iso12Mu3.at(t).Fill(Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT12),1);
	Iso14Mu3.at(t).Fill(Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT14),1);
	Iso16Mu3.at(t).Fill(Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT16),1);
	Iso18Mu3.at(t).Fill(Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT18),1);
	Iso2Mu3.at(t).Fill(Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT2),1);
    

	var_Iso02Mu3 = Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT02);
	var_Iso04Mu3 = Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT04);
	var_Iso06Mu3 = Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT06);
	var_Iso08Mu3 = Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT08);
	var_Iso1Mu3 = Muon3LV.Pt()/  (Muon3LV.Pt() + SumPT1);
	var_Iso12Mu3 = Muon3LV.Pt()/ (Muon3LV.Pt() + SumPT12);
	var_Iso08MuMax = std::max({Iso08Muon1,Iso08Muon2,Iso08Muon3});
	var_Iso08MuMin = std::min({Iso08Muon1,Iso08Muon2,Iso08Muon3});


	var_Iso08MuMax = std::max({Iso08Muon1,Iso08Muon2,Iso08Muon3});
	var_Iso08MuMin = std::min({Iso08Muon1,Iso08Muon2,Iso08Muon3});


	//	Iso08MuMax.at(t).Fill(std::max({Iso08Muon1,Iso08Muon2,Iso08Muon3}),1);
	//	Iso08MuMin.at(t).Fill(std::min({Iso08Muon1,Iso08Muon2,Iso08Muon3}),1);


	// -----------------------------------------------------------------------

	// -------------------------- Fill MVA mini tree ------------------------

	var_vertexKFChi2 = Ntp->Vertex_signal_KF_Chi2(final_idx);
	var_svpvTauAngle = SVPV.Angle(TauLV.Vect());
	
	var_flightLenSig =  Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(final_idx),Ntp->Vertex_PrimaryVertex_Covariance(final_idx),
							   Ntp->Vertex_Signal_KF_pos(final_idx),Ntp->Vertex_Signal_KF_Covariance(final_idx));
	

	var_flightLenDist =  (Ntp->Vertex_MatchedPrimaryVertex(final_idx) - Ntp->Vertex_Signal_KF_pos(final_idx)).Mag();
	
	var_sumMuTrkKinkChi2 = (Ntp->Muon_combinedQuality_trkKink(Muon_index_1)+Ntp->Muon_combinedQuality_trkKink(Muon_index_2)+Ntp->Muon_combinedQuality_trkKink(Muon_index_3));
	var_segCompMuMin = MinSegmentCompatibility;
	var_segCompMuMax = MaxSegmentCompatibility;
	var_segCompMu1 = Ntp->Muon_segmentCompatibility(Muon_index_1);
	var_segCompMu2 = Ntp->Muon_segmentCompatibility(Muon_index_2);
	var_segCompMu3 = Ntp->Muon_segmentCompatibility(Muon_index_3);
	var_caloCompMin = MinMIPLikelihood;
	var_caloCompMax = MaxMIPLikelihood;
	var_caloCompMu1 = Ntp->Muon_caloCompatibility(Muon_index_1);
	var_caloCompMu2 = Ntp->Muon_caloCompatibility(Muon_index_2);
	var_caloCompMu3 = Ntp->Muon_caloCompatibility(Muon_index_3);
	var_MinMIPLikelihood = MinMIPLikelihood;
	var_tauMass = TauRefitLV.M();
	var_tauEta = TauRefitLV.Eta();
	var_ntracks = Ntp->Isolation_NTracks(final_idx);
	var_relPt = Ntp->Isolation_RelPt(final_idx);
	var_isoMax = Ntp->Isolation_maxdy(final_idx);
	
	var_maxMuonsDca = std::max({Ntp->Vertex_DCA12(final_idx),Ntp->Vertex_DCA23(final_idx),Ntp->Vertex_DCA31(final_idx)});
	var_minMuonsDca = std::min({Ntp->Vertex_DCA12(final_idx),Ntp->Vertex_DCA23(final_idx),Ntp->Vertex_DCA31(final_idx)});
	
	var_nsv = NumberOfPrimaryVertices;


	var_VertexMu1D0SigPVReco = Ntp->Vertex_d0sig_reco(final_idx,0);
	var_VertexMu2D0SigPVReco = Ntp->Vertex_d0sig_reco(final_idx,1);
	var_VertexMu3D0SigPVReco = Ntp->Vertex_d0sig_reco(final_idx,2);
	
	var_MaxD0SigPV = MaxD0Significance;
	var_MinD0SigPV = MinD0Significance;
	

	var_VertexMu1D0SigBSReco = Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,0);
	var_VertexMu2D0SigBSReco = Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,1);
	var_VertexMu3D0SigBSReco = Ntp->Vertex_d0BeamSpot_reco_sig(final_idx,2);

	var_MaxD0SigBS = MaxD0BSSignificance;
	var_MinD0SigBS = MinD0BSSignificance;
	
	
	var_VertexMu1D0SigSVReco = Ntp->Vertex_d0sigSV_reco(final_idx,0);
	var_VertexMu2D0SigSVReco = Ntp->Vertex_d0sigSV_reco(final_idx,1);
	var_VertexMu3D0SigSVReco = Ntp->Vertex_d0sigSV_reco(final_idx,2);
	
	
	var_MaxD0SigSV = MaxD0SVSignificance;
	var_MinD0SigSV = MinD0SVSignificance;
	
	
	var_MinMuon_chi2LocalPosition =   std::min({Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3)  });
	var_MaxMuon_chi2LocalPosition = std::max({Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3)  });


	var_Muon1_chi2LocalPosition = Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1);
	var_Muon2_chi2LocalPosition = Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2);
	var_Muon3_chi2LocalPosition = Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3);
	
	
	var_MinMuon_chi2LocalMomentum = std::min({Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_3)  });
	var_MaxMuon_chi2LocalMomentum = std::max({Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_1),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_2),Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_3)  });



	var_Muon1_chi2LocalMomentum = Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_1);
	var_Muon2_chi2LocalMomentum = Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_2);
	var_Muon3_chi2LocalMomentum = Ntp->Muon_combinedQuality_chi2LocalMomentum(Muon_index_3);
	
	
	
	var_MintrkKink = std::min({Ntp->Muon_combinedQuality_trkKink(Muon_index_1),Ntp->Muon_combinedQuality_trkKink(Muon_index_2),Ntp->Muon_combinedQuality_trkKink(Muon_index_3)});
	var_MaxtrkKink = std::max({Ntp->Muon_combinedQuality_trkKink(Muon_index_1),Ntp->Muon_combinedQuality_trkKink(Muon_index_2),Ntp->Muon_combinedQuality_trkKink(Muon_index_3)});
	var_MinglbKink = std::min({Ntp->Muon_combinedQuality_glbKink(Muon_index_1),Ntp->Muon_combinedQuality_glbKink(Muon_index_2),Ntp->Muon_combinedQuality_glbKink(Muon_index_3)});
	var_MaxglbKink = std::max({Ntp->Muon_combinedQuality_glbKink(Muon_index_1),Ntp->Muon_combinedQuality_glbKink(Muon_index_2),Ntp->Muon_combinedQuality_glbKink(Muon_index_3)});
	
	var_MuonglbkinkSum = var_sumMuTrkKinkChi2;

	var_MaxVertexPairQuality =   std::max({Ntp->Vertex_pair_quality(final_idx,0),Ntp->Vertex_pair_quality(final_idx,1),Ntp->Vertex_pair_quality(final_idx,2)});
	var_MinVertexPairQuality =   std::min({Ntp->Vertex_pair_quality(final_idx,0),Ntp->Vertex_pair_quality(final_idx,1),Ntp->Vertex_pair_quality(final_idx,2)});
	
	//---------------------
	var_NtracksClose = NcloseTracksCount;
	var_NTracksCloseToPV=NPVcloseTracksCount;
	var_NTracksCloseToPVTauDR=NPVcloseDRToTauTracksCount;



	var_Muon1Pt = Ntp->Muon_P4(Muon_index_1).Pt();
	var_Muon2Pt = Ntp->Muon_P4(Muon_index_2).Pt();
	var_Muon3Pt = Ntp->Muon_P4(Muon_index_3).Pt();

 	var_Muon1P = Ntp->Muon_P4(Muon_index_1).P();
	var_Muon2P = Ntp->Muon_P4(Muon_index_2).P();
	var_Muon3P = Ntp->Muon_P4(Muon_index_3).P();



	if(Ntp->NIsolationTrack(final_idx)!=0){ 
	  var_MindcaTrackSV = sqrt( pow(Ntp->IsolationTrack_dzSV(final_idx,TrackIndex),2)   +   pow(Ntp->IsolationTrack_dxySV(final_idx,TrackIndex),2));
	} else var_MindcaTrackSV = -1;

	if(Ntp->NIsolationTrack(final_idx)!=0){ 
	  var_dcaTrackPV = sqrt(  pow(Ntp->IsolationTrack_dzPV(final_idx,TrackIndex_closestToPV),2)   +   pow(Ntp->IsolationTrack_dxyPV(final_idx,TrackIndex_closestToPV),2));
	} else var_dcaTrackPV = -1;



	if (id==1) MC=0;
	else  MC=1;
	
	if (tauMassRes<tauMassResCutLow) category = 1;
	if (tauMassRes>tauMassResCutLow && tauMassRes<tauMassResCutHigh) category = 2;
	if (tauMassRes>tauMassResCutHigh) category = 3;
	
	
	var_Muon1LooseId = Ntp->CHECK_BIT(Ntp->Muon_ID(Muon_index_1),Ntp->MuonQualityBitMask::Bit_MuonLoose);
	var_Muon1MediumId= Ntp->CHECK_BIT(Ntp->Muon_ID(Muon_index_1),Ntp->MuonQualityBitMask::Bit_MuonMedium);
	var_Muon1TightId= Ntp->CHECK_BIT(Ntp->Muon_ID(Muon_index_1),Ntp->MuonQualityBitMask::Bit_MuonTight);

	var_Muon2LooseId= Ntp->CHECK_BIT(Ntp->Muon_ID(Muon_index_2),Ntp->MuonQualityBitMask::Bit_MuonLoose);
	var_Muon2MediumId= Ntp->CHECK_BIT(Ntp->Muon_ID(Muon_index_2),Ntp->MuonQualityBitMask::Bit_MuonMedium);
	var_Muon2TightId= Ntp->CHECK_BIT(Ntp->Muon_ID(Muon_index_2),Ntp->MuonQualityBitMask::Bit_MuonTight);

	var_Muon3LooseId= Ntp->CHECK_BIT(Ntp->Muon_ID(Muon_index_3),Ntp->MuonQualityBitMask::Bit_MuonLoose);
	var_Muon3MediumId= Ntp->CHECK_BIT(Ntp->Muon_ID(Muon_index_3),Ntp->MuonQualityBitMask::Bit_MuonMedium);
	var_Muon3TightId= Ntp->CHECK_BIT(Ntp->Muon_ID(Muon_index_3),Ntp->MuonQualityBitMask::Bit_MuonTight);

	var_Muon1PFIsoLoose = Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_1),Ntp->MuonStandardSelectors::PFIsoLoose);
	var_Muon1PFIsoMedium = Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_1),Ntp->MuonStandardSelectors::PFIsoMedium);
	var_Muon1PFIsoTight = Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_1),Ntp->MuonStandardSelectors::PFIsoTight);
	var_Muon1PFIsoVTight = Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_1),Ntp->MuonStandardSelectors::PFIsoVeryTight);

	var_Muon2PFIsoLoose = Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_2),Ntp->MuonStandardSelectors::PFIsoLoose);
	var_Muon2PFIsoMedium = Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_2),Ntp->MuonStandardSelectors::PFIsoMedium);
	var_Muon2PFIsoTight = Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_2),Ntp->MuonStandardSelectors::PFIsoTight);
	var_Muon2PFIsoVTight = Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_2),Ntp->MuonStandardSelectors::PFIsoVeryTight);

	var_Muon3PFIsoLoose = Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_3),Ntp->MuonStandardSelectors::PFIsoLoose);
	var_Muon3PFIsoMedium = Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_3),Ntp->MuonStandardSelectors::PFIsoMedium);
	var_Muon3PFIsoTight = Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_3),Ntp->MuonStandardSelectors::PFIsoTight);
	var_Muon3PFIsoVTight = Ntp->CHECK_BIT(Ntp->Muon_StandardSelection(Muon_index_3),Ntp->MuonStandardSelectors::PFIsoVeryTight);


	
	// ---------------- Evaluate B vs D mva

	var_BvsDSeprator = readerBvsD->EvaluateMVA("BDTG");
	Separation_BvsD.at(t).Fill(var_BvsDSeprator,1);

	// ---------------- Evaluate B vs D mva


	TMVA_Tree->Fill();

	}
    }
}



void  MakeMVATree::Finish(){
  
  if(mode == RECONSTRUCT)
    {
      //    for(unsigned int i=1; i<  Nminus0.at(0).size(); i++){
      //    int id(Ntp->GetMCID());
      double scale(1.);
      double scaleDsTau(0.637);
      double scaleBpTau(0.262);
      double scaleB0Tau(0.099);
      double scaleDpTau(0.032);
      
      if(Nminus0.at(0).at(1).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(1).Integral();
      ScaleAllHistOfType(1,scale*scaleDsTau);
    
      if(Nminus0.at(0).at(2).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(2).Integral();
      ScaleAllHistOfType(2,scale*scaleBpTau);
      
      if(Nminus0.at(0).at(3).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(3).Integral();
      ScaleAllHistOfType(3,scale*scaleB0Tau);

      if(Nminus0.at(0).at(4).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(4).Integral();
      ScaleAllHistOfType(4,scale*scaleDpTau);
      
      //    }
    }
  file= new TFile("TMVATreesInput.root","recreate");
  TMVA_Tree->SetDirectory(file);
  
  file->Write();
  file->Close();
    
  Selection::Finish();
}
