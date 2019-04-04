#include "MCEfficiency.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

double MCEfficiency::deltaR(double eta1, double phi1, double eta2, double phi2)
{
    double deta = eta1 - eta2;
    double dphi = phi1 - phi2;
    while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
    while (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();

    return sqrt(deta*deta + dphi*dphi);
}

MCEfficiency::MCEfficiency(TString Name_, TString id_):
    Selection(Name_,id_)
{


    // This is a class constructor;
}

MCEfficiency::~MCEfficiency(){
    for(unsigned int j=0; j<Npassed.size(); j++){
      Logger(Logger::Info) << "Selection Summary before: "
        << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
        << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
    }


    Logger(Logger::Info) << "complete." << std::endl;
}

void  MCEfficiency::Configure(){

    for(int i=0; i<NCuts;i++){
      cut.push_back(0);
      value.push_back(0);
      pass.push_back(false);
      if(i==L1SeedOk)     cut.at(L1SeedOk)=1;
      if(i==HLTOk)        cut.at(HLTOk)=1;
      if(i==is2MuTrk)      cut.at(is2MuTrk)=1;
      if(i==PrimeVtx)      cut.at(PrimeVtx)=5; // Here for example we place cut value on number of PVs
      if(i==trkPt)			cut.at(trkPt)=1.2;
      if(i==nTrkHits)		cut.at(nTrkHits)=12;
		if(i==mumuMass)		cut.at(mumuMass)=PDG_Var::Phi_mass();
		if(i==fitVtxChiSq)	cut.at(fitVtxChiSq)=5.0;
    	}

    TString hlabel;
    TString htitle;

    for(int i=0; i<NCuts; i++){
      title.push_back("");
      distindx.push_back(false);
      dist.push_back(std::vector<float>());
      TString c="_Cut_";c+=i;
      // book here the N-1 and N-0 histrogramms for each cut
      if(i==PrimeVtx){
        title.at(i)="Number of Prime Vertices $(N>$";
        title.at(i)+=cut.at(PrimeVtx);
        title.at(i)+=")";
        htitle=title.at(i);
        htitle.ReplaceAll("$","");
        htitle.ReplaceAll("\\","#");
        hlabel="Number of Prime Vertices";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,61,-0.5,60.5,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,61,-0.5,60.5,hlabel,"Events"));
      }

      else if(i==L1SeedOk){
        title.at(i)="L1 seed ";
        hlabel="L1 triggers";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_L1SeedOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_L1SeedOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }

      else if(i==is2MuTrk){
        title.at(i)="is2MuTrk ";
        hlabel="2muon + track category";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_mumuMass_",htitle,2,-0.5,1.5,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_mumuMass_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
      
		else if(i==mumuMass){
        title.at(i)="mumuMass";
        hlabel="mumuMass";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_mumuMass_",htitle,80,0.8,1.6,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_mumuMass_",htitle,80,0.8,1.6,hlabel,"Events"));
      }

      else if(i==HLTOk){
        title.at(i)="HLT trigger ";
        hlabel="DoubleMu3_Trk_Tau3mu";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }

      else if(i==trkPt){
        title.at(i)="pT of the track";
        hlabel="pT of the track";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_trkPt_",htitle,100,-0.5,9.5,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_trkPt_",htitle,100,-0.5,9.5,hlabel,"Events"));
      }

      else if(i==nTrkHits){
        title.at(i)="Number of hits in the tracker";
        hlabel="Number of hits in the tracker";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nTrkHits_",htitle,11,-1,10,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_NMinus0_nTrkHits_",htitle,11,-1,10,hlabel,"Events"));
      }
		
		else if(i==fitVtxChiSq){
        title.at(i)="Normalized chi square of the vertex fit";
        hlabel="Normalized chi square of the vertex fit";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_fitVtxChiSq_",htitle,100,-0.5,9.5,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_NMinus0_fitVtxChiSq_",htitle,100,-0.5,9.5,hlabel,"Events"));
      }

      // Track Candidate Information
      Track_P=HConfig.GetTH1D(Name+"_Track_P","Momentum magnitude of track (all tracks)",66,-0.5,65.5,"Number of Vertices","Events");
      Track_E=HConfig.GetTH1D(Name+"_Track_E","Energy of track (all tracks)",66,-0.5,65.5,"","Events");
      Track_Pt=HConfig.GetTH1D(Name+"_Track_Pt","Transverse momentum of track (all tracks)",66,-0.5,65.5,"","Events");
      Track_Eta=HConfig.GetTH1D(Name+"_Track_Eta","Psuedorapidity of track (all tracks)",66,-0.5,65.5,"Number of Vertices","Events");
      Track_Phi=HConfig.GetTH1D(Name+"_Track_Phi","Azimuthal angle of track (all tracks)",66,-0.5,65.5,"Number of Vertices","Events");
      Track_vx=HConfig.GetTH1D(Name+"_Track_vx","X coordinate of the parent vertex (all tracks)",66,-0.5,65.5,"Number of Vertices","Events");
      Track_vy=HConfig.GetTH1D(Name+"_Track_vy","Y coordinate of the parent vertex (all tracks)",66,-0.5,65.5,"Number of Vertices","Events");
      Track_vz=HConfig.GetTH1D(Name+"_Track_vz","Z coordinate of the parent vertex (all tracks)",66,-0.5,65.5,"Number of Vertices","Events");
      Track_normalizedChi2=HConfig.GetTH1D(Name+"_Track_normalizedChi2","Normalized chi square",66,-0.5,65.5,"Number of Vertices","Events");
      Track_numberOfValidHits=HConfig.GetTH1D(Name+"_Track_numberOfValidHits","number of valid hits in te tracker",66,-0.5,65.5,"Number of Vertices","Events");
      Track_charge=HConfig.GetTH1D(Name+"_Track_charge","Chargeof the track",66,-0.5,65.5,"Number of Vertices","Events");
      Track_dxy=HConfig.GetTH1D(Name+"_Track_dxy","Transverse displacement of the parent vertex from the bs",66,-0.5,65.5,"Number of Vertices","Events");
      Track_dz=HConfig.GetTH1D(Name+"_Track_dz","Longitudnal displacement of the parent vertex from the bs",66,-0.5,65.5,"Number of Vertices","Events");
      Track_dxyError=HConfig.GetTH1D(Name+"_Track_dxyError","dxy Error",66,-0.5,65.5,"Number of Vertices","Events");
      Track_dzError=HConfig.GetTH1D(Name+"_Track_dzError","dz Error",66,-0.5,65.5,"Number of Vertices","Events");

      //Muon variables (Muons from dimuon + track candidates)
      Muon1_isGlobal=HConfig.GetTH1D(Name+"_Muon1_isGlobal","Global muons status ",2,-.5,1.5,"#mu_{1} isGlb","Events");
      Muon2_isGlobal=HConfig.GetTH1D(Name+"_Muon2_isGlobal","",2,-0.5,0.5,"#mu_{2} isGlb","Events");
      Muon1_isStandAlone=HConfig.GetTH1D(Name+"_Muon1_isStandAlone","",2,-0.5,1.5,"#mu_{1} isStandAlone","Events");
      Muon2_isStandAlone=HConfig.GetTH1D(Name+"_Muon2_isStandAlone","",2,-0.5,1.5,"#mu_{2} isStandAlone","Events");
      Muon1_isTracker=HConfig.GetTH1D(Name+"_Muon1_isTracker","",2,-0.5,1.5,"#mu_{1} isTracker","Events");
      Muon2_isTracker=HConfig.GetTH1D(Name+"_Muon2_isTracker","",2,-0.5,1.5,"#mu_{2} isTracker","Events");
      Muon1_isCalo=HConfig.GetTH1D(Name+"_Muon1_isCaloMuon","",2,-0.5,1.5,"#mu_{1} isCalo","Events");
      Muon2_isCalo=HConfig.GetTH1D(Name+"_Muon2_isCaloMuon","",2,-0.5,1.5,"#mu_{2} isCalo","Events");
      Muon1_isIsolationValid=HConfig.GetTH1D(Name+"_Muon1_isIsolationValid","#mu_{1} isIsoValid",2,-0.5,1.5,"","Events");
      Muon2_isIsolationValid=HConfig.GetTH1D(Name+"_Muon2_isIsolationValid","#mu_{2} isIsoValid",2,-0.5,1.5,"","Events");
      Muon1_isTimeValid=HConfig.GetTH1D(Name+"_Muon1_isTimeValid","#mu_{1} isTimevalid",2,-0.5,1.5,"","Events");
      Muon2_isTimeValid=HConfig.GetTH1D(Name+"_Muon2_isTimeValid","#mu_{2} isTimeValid",2,-0.5,1.5,"","Events");
      Muon1_emEt03=HConfig.GetTH1D(Name+"_Muon1_emEt03","",10,0,10,"","Events");
      Muon2_emEt03=HConfig.GetTH1D(Name+"_Muon2_emEt03","",10,0,10,"","Events");
      Muon1_emVetoEt03=HConfig.GetTH1D(Name+"_Muon1_emVetoEt03","",10,0,10,"","Events");
      Muon2_emVetoEt03=HConfig.GetTH1D(Name+"_Muon2_emVetoEt03","",10,0,10,"","Events");
      Muon1_hadEt03=HConfig.GetTH1D(Name+"_Muon1_hadEt03","",10,0,10,"","Events");
      Muon2_hadEt03=HConfig.GetTH1D(Name+"_Muon2_hadEt03","",10,0,10,"","Events");
      Muon1_hadVetoEt03=HConfig.GetTH1D(Name+"_Muon1_hadVetoEt03","",10,0,10,"","Events");
      Muon2_hadVetoEt03=HConfig.GetTH1D(Name+"_Muon2_hadVetoEt03","",10,0,10,"","Events");
      Muon1_nJets03=HConfig.GetTH1D(Name+"_Muon1_nJets03","",10,0,10,"","Events");
      Muon2_nJets03=HConfig.GetTH1D(Name+"_Muon2_nJets03","",10,0,10,"","Events");
      Muon1_nTracks03=HConfig.GetTH1D(Name+"_Muon1_nTracks03","",10,0,10,"","Events");
      Muon2_nTracks03=HConfig.GetTH1D(Name+"_Muon2_nTracks03","",10,0,10,"","Events");
      Muon1_sumPt03=HConfig.GetTH1D(Name+"_Muon1_sumPt03","",10,0,10,"","Events");
      Muon2_sumPt03=HConfig.GetTH1D(Name+"_Muon2_sumPt03","",10,0,10,"","Events");
      Muon1_trackerVetoPt03=HConfig.GetTH1D(Name+"_Muon1_trackerVetoPt03","",10,0,10,"","Events");
      Muon2_trackerVetoPt03=HConfig.GetTH1D(Name+"_Muon2_trackerVetoPt03","",10,0,10,"","Events");
      Muon1_emEt05=HConfig.GetTH1D(Name+"_Muon1_emEt05","",10,0,10,"","Events");
      Muon2_emEt05=HConfig.GetTH1D(Name+"_Muon2_emEt05","",10,0,10,"","Events");
      Muon1_emVetoEt05=HConfig.GetTH1D(Name+"_Muon1_emVetoEt05","",10,0,10,"","Events");
      Muon2_emVetoEt05=HConfig.GetTH1D(Name+"_Muon2_emVetoEt05","",10,0,10,"","Events");
      Muon1_hadEt05=HConfig.GetTH1D(Name+"_Muon1_hadEt05","",10,0,10,"","Events");
      Muon2_hadEt05=HConfig.GetTH1D(Name+"_Muon2_hadEt05","",10,0,10,"","Events");
      Muon1_hadVetoEt05=HConfig.GetTH1D(Name+"_Muon1_hadVetoEt05","",10,0,10,"","Events");
      Muon2_hadVetoEt05=HConfig.GetTH1D(Name+"_Muon2_hadVetoEt05","",10,0,10,"","Events");
      Muon1_nJets05=HConfig.GetTH1D(Name+"_Muon1_nJets05","",10,0,10,"","Events");
      Muon2_nJets05=HConfig.GetTH1D(Name+"_Muon2_nJets05","",10,0,10,"","Events");
      Muon1_nTracks05=HConfig.GetTH1D(Name+"_Muon1_nTracks05","",10,0,10,"","Events");
      Muon2_nTracks05=HConfig.GetTH1D(Name+"_Muon2_nTracks05","",10,0,10,"","Events");
      Muon1_sumPt05=HConfig.GetTH1D(Name+"_Muon1_sumPt05","",10,0,10,"","Events");
      Muon2_sumPt05=HConfig.GetTH1D(Name+"_Muon2_sumPt05","",10,0,10,"","Events");
      Muon1_trackerVetoPt05=HConfig.GetTH1D(Name+"_Muon1_trackerVetoPt05","",10,0,10,"","Events");
      Muon2_trackerVetoPt05=HConfig.GetTH1D(Name+"_Muon2_trackerVetoPt05","",10,0,10,"","Events");
      Muon1_sumChargedHadronPt03=HConfig.GetTH1D(Name+"_Muon1_sumChargedHadronPt03","",10,0,10,"","Events");
      Muon2_sumChargedHadronPt03=HConfig.GetTH1D(Name+"_Muon2_sumChargedHadronPt03","",10,0,10,"","Events");
      Muon1_sumChargedParticlePt03=HConfig.GetTH1D(Name+"_Muon1_sumChargedParticlePt03","",10,0,10,"","Events");
      Muon2_sumChargedParticlePt03=HConfig.GetTH1D(Name+"_Muon2_sumChargedParticlePt03","",10,0,10,"","Events");
      Muon1_sumNeutralHadronEt03=HConfig.GetTH1D(Name+"_Muon1_sumNeutralHadronEt03","",10,0,10,"","Events");
      Muon2_sumNeutralHadronEt03=HConfig.GetTH1D(Name+"_Muon2_sumNeutralHadronEt03","",10,0,10,"","Events");
      Muon1_sumNeutralHadronEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon1_sumNeutralHadronEtHighThreshold03","",10,0,10,"","Events");
      Muon2_sumNeutralHadronEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon2_sumNeutralHadronEtHighThreshold03","",10,0,10,"","Events");
      Muon1_sumPhotonEt03=HConfig.GetTH1D(Name+"_Muon1_sumPhotonEt03","",10,0,10,"","Events");
      Muon2_sumPhotonEt03=HConfig.GetTH1D(Name+"_Muon2_sumPhotonEt03","",10,0,10,"","Events");
      Muon1_sumPhotonEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon1_sumPhotonEtHighThreshold03","",10,0,10,"","Events");
      Muon2_sumPhotonEtHighThreshold03=HConfig.GetTH1D(Name+"_Muon2_sumPhotonEtHighThreshold03","",10,0,10,"","Events");
      Muon1_sumPUPt03=HConfig.GetTH1D(Name+"_Muon1_sumPUPt03","",10,0,10,"","Events");
      Muon2_sumPUPt03=HConfig.GetTH1D(Name+"_Muon2_sumPUPt03","",10,0,10,"","Events");
      Muon1_sumChargedHadronPt04=HConfig.GetTH1D(Name+"_Muon1_sumChargedHadronPt04","",10,0,10,"","Events");
      Muon2_sumChargedHadronPt04=HConfig.GetTH1D(Name+"_Muon2_sumChargedHadronPt04","",10,0,10,"","Events");
      Muon1_sumChargedParticlePt04=HConfig.GetTH1D(Name+"_Muon1_sumChargedParticlePt04","",10,0,10,"","Events");
      Muon2_sumChargedParticlePt04=HConfig.GetTH1D(Name+"_Muon2_sumChargedParticlePt04","",10,0,10,"","Events");
      Muon1_sumNeutralHadronEt04=HConfig.GetTH1D(Name+"_Muon1_sumNeutralHadronEt04","",10,0,10,"","Events");
      Muon2_sumNeutralHadronEt04=HConfig.GetTH1D(Name+"_Muon2_sumNeutralHadronEt04","",10,0,10,"","Events");
      Muon1_sumNeutralHadronEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon1_sumNeutralHadronEtHighThreshold04","",10,0,10,"","Events");
      Muon2_sumNeutralHadronEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon2_sumNeutralHadronEtHighThreshold04","",10,0,10,"","Events");
      Muon1_sumPhotonEt04=HConfig.GetTH1D(Name+"_Muon1_sumPhotonEt04","",10,0,10,"","Events");
      Muon2_sumPhotonEt04=HConfig.GetTH1D(Name+"_Muon2_sumPhotonEt04","",10,0,10,"","Events");
      Muon1_sumPhotonEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon1_sumPhotonEtHighThreshold04","",10,0,10,"","Events");
      Muon2_sumPhotonEtHighThreshold04=HConfig.GetTH1D(Name+"_Muon2_sumPhotonEtHighThreshold04","",10,0,10,"","Events");
      Muon1_sumPUPt04=HConfig.GetTH1D(Name+"_Muon1_sumPUPt04","",10,0,10,"","Events");
      Muon2_sumPUPt04=HConfig.GetTH1D(Name+"_Muon2_sumPUPt04","",10,0,10,"","Events");

      //Dimuon Information (Muons from dimuon + track candidates)
      MuonsPtRatio=HConfig.GetTH1D(Name+"_MuonsPtRatio","Ratio of Pt of two muons",50,0.1,1.2,"Ratio of first and second muon p_{T}","Events");
      DimuondR=HConfig.GetTH1D(Name+"_DimuondR","dR between the muon pair",20,0,1,"dR","Events");
      Muon1TrkdR=HConfig.GetTH1D(Name+"_Muon1TrkdR","dR between the highest p muon and the track",100,0,5,"dR","Events");
      Muon2TrkdR=HConfig.GetTH1D(Name+"_Muon2TrkdR","dR between the lowest p muon and the track",100,0,5,"dR","Events");
      PhiMass=HConfig.GetTH1D(Name+"_PhiMass","#mu#mu mass",50,0.2,1.5,"Mass of the #mu#mu pair","Events");
      TripleMass=HConfig.GetTH1D(Name+"_TripleMass","#mu#mu + track mass",50,1.7,2.1,"Mass of the #mu#mu + track","Events");
      PhiMassVsDsMass=HConfig.GetTH2D(Name+"_PhiMassVsDsMass","#mu#mu Mass vs. #mu#mu + track mass",50,0.2,1.5,50,1.7,2.1,"M_{#mu#mu}, GeV","M_{#mu#mu + track}, GeV");
    } 
    // Setup NPassed Histogams
    Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
    // Setup Extra Histograms
    // Book here your analysis histrogramms, a good style is to follow selfexplanatory convention
    NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",66,-0.5,65.5,"Number of Vertices","Events");

    Selection::ConfigureHistograms(); //do not remove
    HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}

void  MCEfficiency::Store_ExtraDist(){ 

    //Track candidate variables 
    Extradist1d.push_back(&Track_P);
    Extradist1d.push_back(&Track_E);
    Extradist1d.push_back(&Track_Pt);
    Extradist1d.push_back(&Track_Eta);
    Extradist1d.push_back(&Track_Phi);
    Extradist1d.push_back(&Track_vx);
    Extradist1d.push_back(&Track_vy);
    Extradist1d.push_back(&Track_vz);
    Extradist1d.push_back(&Track_normalizedChi2);
    Extradist1d.push_back(&Track_numberOfValidHits);
    Extradist1d.push_back(&Track_charge);
    Extradist1d.push_back(&Track_dxy);
    Extradist1d.push_back(&Track_dz);
    Extradist1d.push_back(&Track_dxyError);
    Extradist1d.push_back(&Track_dzError);

    //Dimuon variables
    Extradist1d.push_back(&MuonsPtRatio);
    Extradist1d.push_back(&DimuondR);
    Extradist1d.push_back(&Muon1TrkdR);
    Extradist1d.push_back(&Muon2TrkdR);
    Extradist1d.push_back(&PhiMass);
    Extradist1d.push_back(&TripleMass);
    Extradist2d.push_back(&PhiMassVsDsMass);
    Extradist1d.push_back(&Muon1_isGlobal);
    Extradist1d.push_back(&Muon2_isGlobal);
    Extradist1d.push_back(&Muon1_isStandAlone);
    Extradist1d.push_back(&Muon2_isStandAlone);
    Extradist1d.push_back(&Muon1_isTracker);
    Extradist1d.push_back(&Muon2_isTracker);
    Extradist1d.push_back(&Muon1_isCalo);
    Extradist1d.push_back(&Muon2_isCalo);
    Extradist1d.push_back(&Muon1_isIsolationValid);
    Extradist1d.push_back(&Muon2_isIsolationValid);
    Extradist1d.push_back(&Muon1_isTimeValid);
    Extradist1d.push_back(&Muon2_isTimeValid);
    Extradist1d.push_back(&Muon1_emEt03);
    Extradist1d.push_back(&Muon2_emEt03);
    Extradist1d.push_back(&Muon1_emVetoEt03);
    Extradist1d.push_back(&Muon2_emVetoEt03);
    Extradist1d.push_back(&Muon1_hadEt03);
    Extradist1d.push_back(&Muon2_hadEt03);
    Extradist1d.push_back(&Muon1_hadVetoEt03);
    Extradist1d.push_back(&Muon2_hadVetoEt03);
    Extradist1d.push_back(&Muon1_nJets03);
    Extradist1d.push_back(&Muon2_nJets03);
    Extradist1d.push_back(&Muon1_nTracks03);
    Extradist1d.push_back(&Muon2_nTracks03);
    Extradist1d.push_back(&Muon1_sumPt03);
    Extradist1d.push_back(&Muon2_sumPt03);
    Extradist1d.push_back(&Muon1_trackerVetoPt03);
    Extradist1d.push_back(&Muon2_trackerVetoPt03);
    Extradist1d.push_back(&Muon1_emEt05);
    Extradist1d.push_back(&Muon2_emEt05);
    Extradist1d.push_back(&Muon1_emVetoEt05);
    Extradist1d.push_back(&Muon2_emVetoEt05);
    Extradist1d.push_back(&Muon1_hadEt05);
    Extradist1d.push_back(&Muon2_hadEt05);
    Extradist1d.push_back(&Muon1_hadVetoEt05);
    Extradist1d.push_back(&Muon2_hadVetoEt05);
    Extradist1d.push_back(&Muon1_nJets05);
    Extradist1d.push_back(&Muon2_nJets05);
    Extradist1d.push_back(&Muon1_nTracks05);
    Extradist1d.push_back(&Muon2_nTracks05);
    Extradist1d.push_back(&Muon1_sumPt05);
    Extradist1d.push_back(&Muon2_sumPt05);
    Extradist1d.push_back(&Muon1_trackerVetoPt05);
    Extradist1d.push_back(&Muon2_trackerVetoPt05);
    Extradist1d.push_back(&Muon1_sumChargedHadronPt03);
    Extradist1d.push_back(&Muon2_sumChargedHadronPt03);
    Extradist1d.push_back(&Muon1_sumChargedParticlePt03);
    Extradist1d.push_back(&Muon2_sumChargedParticlePt03);
    Extradist1d.push_back(&Muon1_sumNeutralHadronEt03);
    Extradist1d.push_back(&Muon2_sumNeutralHadronEt03);
    Extradist1d.push_back(&Muon1_sumNeutralHadronEtHighThreshold03);
    Extradist1d.push_back(&Muon2_sumNeutralHadronEtHighThreshold03);
    Extradist1d.push_back(&Muon1_sumPhotonEt03);
    Extradist1d.push_back(&Muon2_sumPhotonEt03);
    Extradist1d.push_back(&Muon1_sumPhotonEtHighThreshold03);
    Extradist1d.push_back(&Muon2_sumPhotonEtHighThreshold03);
    Extradist1d.push_back(&Muon1_sumPUPt03);
    Extradist1d.push_back(&Muon2_sumPUPt03);
    Extradist1d.push_back(&Muon1_sumChargedHadronPt04);
    Extradist1d.push_back(&Muon2_sumChargedHadronPt04);
    Extradist1d.push_back(&Muon1_sumChargedParticlePt04);
    Extradist1d.push_back(&Muon2_sumChargedParticlePt04);
    Extradist1d.push_back(&Muon1_sumNeutralHadronEt04);
    Extradist1d.push_back(&Muon2_sumNeutralHadronEt04);
    Extradist1d.push_back(&Muon1_sumNeutralHadronEtHighThreshold04);
    Extradist1d.push_back(&Muon2_sumNeutralHadronEtHighThreshold04);
    Extradist1d.push_back(&Muon1_sumPhotonEt04);
    Extradist1d.push_back(&Muon2_sumPhotonEt04);
    Extradist1d.push_back(&Muon1_sumPhotonEtHighThreshold04);
    Extradist1d.push_back(&Muon2_sumPhotonEtHighThreshold04);
    Extradist1d.push_back(&Muon1_sumPUPt04);
    Extradist1d.push_back(&Muon2_sumPUPt04);

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // Here you must push back all analysis histograms, otherwise they wont be propagated to the output //
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    Extradist1d.push_back(&NVtx);
}

void  MCEfficiency::doEvent(){ 
    unsigned int t;
    int id(Ntp->GetMCID());
	 double phidM = 0.02;
	 bool DEBUG = false;
    if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
    
	 // Apply Selection
	 cout<<trkPt<<" "<<nTrkHits<<" "<<fitVtxChiSq<<endl;
    value.at(L1SeedOk) = 0;
    value.at(HLTOk) = 0;
    value.at(trkPt) = 0;
    value.at(nTrkHits) = 0;
    value.at(fitVtxChiSq) = 999.0;

    for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
      TString HLT = Ntp->HLTName(iTrigger);
      if(HLT.Contains("DoubleMu3_Trk_Tau3mu") && Ntp->HLTDecision(iTrigger) == 1)value.at(HLTOk)=Ntp->HLTDecision(iTrigger);
    }
    for(int l1iTrigger=0; l1iTrigger < Ntp->NL1Seeds(); l1iTrigger++){
      TString L1 = Ntp->L1Name(l1iTrigger);
      if(L1.Contains("L1_DoubleMu0er1p4_dEta_Max1p8_OS") && Ntp->L1Decision(l1iTrigger) == 1)value.at(L1SeedOk)=Ntp->L1Decision(l1iTrigger);
      if(L1.Contains("L1_DoubleMu_10_0_dEta_Max1p8") && Ntp->L1Decision(l1iTrigger) == 1)value.at(L1SeedOk)=Ntp->L1Decision(l1iTrigger);
      if(L1.Contains("L1_TripleMu0") && Ntp->L1Decision(l1iTrigger) == 1)value.at(L1SeedOk)=Ntp->L1Decision(l1iTrigger);
    }

	if (DEBUG) cout<<"HLT ok"<<endl;
    vector<int> cut_count;
    vector<int> tmp_nTrkHits;
    vector<double> tmp_trkPt;
    vector<double> tmp_fitVtxChiSq;
	 vector<double> tmp_mumuMass;

    int tmp_idx = -1;
    double tmp_chiSq = 999.0;

    value.at(PrimeVtx)=Ntp->NVtx(); 
    value.at(is2MuTrk) = 0;
    if(Ntp->NTwoMuonsTrack()!=0 && Ntp->NThreeMuons() == 0) value.at(is2MuTrk) = 1;
	 if (DEBUG) cout<<"select dimutrk"<<endl;
    double deltaMass(999.);
    unsigned int pair_index(0);

	 if (value.at(is2MuTrk)){
    for(unsigned int i2M=0; i2M < Ntp->NTwoMuonsTrack(); i2M++){

      cut_count.push_back(0);
		tmp_mumuMass.push_back(0);
      tmp_trkPt.push_back(0);
      tmp_nTrkHits.push_back(0);
      tmp_fitVtxChiSq.push_back(999.0);

      unsigned int mu1_idx =  Ntp-> TwoMuonsTrackMuonIndices(i2M).at(0);
      unsigned int mu2_idx =  Ntp-> TwoMuonsTrackMuonIndices(i2M).at(1);
      unsigned int trk_idx = Ntp->TwoMuonsTrackTrackIndex(i2M).at(0);

		if (DEBUG) cout<<mu1_idx<<" "<<mu2_idx<<" "<<trk_idx<<endl;
		if (fabs((Ntp->Muon_P4(mu1_idx) + Ntp->Muon_P4(mu2_idx)).M()-PDG_Var::Phi_mass())<phidM){
			cut_count[i2M]++;
			tmp_mumuMass[i2M] = (Ntp->Muon_P4(mu1_idx)+Ntp->Muon_P4(mu2_idx)).M();
			if (DEBUG) cout<<"Dimuon Invariant mass: "<<tmp_mumuMass[i2M]<<endl;

		}
		else continue;

      if (Ntp->Track_P4(trk_idx).Pt()>cut.at(trkPt)) {
        cut_count[i2M]++;
        tmp_trkPt[i2M] = Ntp->Track_P4(trk_idx).Pt();
		  if (DEBUG) cout<<"Track Pt:"<<tmp_trkPt[i2M]<<endl;
      }
		else continue;

      if (Ntp->Track_numberOfValidHits(trk_idx)>=cut.at(nTrkHits)){
        cut_count[i2M]++;
        tmp_nTrkHits[i2M] = Ntp->Track_numberOfValidHits(trk_idx);
		  cout<<"Track nHits:"<<tmp_nTrkHits[i2M]<<endl;
      }
		else continue;

      if (Ntp->TwoMuonsTrack_SV_Chi2(i2M)<=cut.at(fitVtxChiSq)) {
        cut_count[i2M]++;
        tmp_fitVtxChiSq[i2M] = Ntp->TwoMuonsTrack_SV_Chi2(i2M);
		  if (DEBUG) cout<<"Chi square fit vertex: "<<tmp_fitVtxChiSq[i2M]<<endl;
      }
		else continue;

      if (tmp_chiSq>tmp_fitVtxChiSq[i2M]) {
        tmp_idx = i2M;
        tmp_chiSq = tmp_fitVtxChiSq[i2M];
      }
    }
	 if (DEBUG){ 
	 cout<<"select the best candidate"<<endl;
	 cout<<"Index of the best candidate"<<tmp_idx<<endl;
	 cout<<tmp_trkPt.size()<<" "<<tmp_nTrkHits.size()<<" "<<tmp_fitVtxChiSq.size()<<endl;
	 }
    if (tmp_idx==-1) tmp_idx = std::distance(cut_count.begin(), std::max_element(cut_count.begin(), cut_count.end()));
    value.at(trkPt) = tmp_trkPt[tmp_idx];
    value.at(nTrkHits) = tmp_nTrkHits[tmp_idx];
    value.at(fitVtxChiSq) = tmp_fitVtxChiSq[tmp_idx];
	 value.at(mumuMass) = tmp_mumuMass[tmp_idx];
	 if (DEBUG) cout<<"Final index"<<tmp_idx<<endl;
	 }

	cout<<"Values of the variables:"<<value.at(trkPt)<<" "<<value.at(nTrkHits)<<" "<<value.at(fitVtxChiSq)<<endl;

    pass.at(is2MuTrk) = (value.at(is2MuTrk) == cut.at(is2MuTrk));
    pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx)); 
    pass.at(L1SeedOk)= (value.at(L1SeedOk)==cut.at(L1SeedOk)); 
    pass.at(HLTOk)= (value.at(HLTOk)==cut.at(HLTOk)); 
    pass.at(nTrkHits)=(value.at(nTrkHits)>=cut.at(nTrkHits));
    pass.at(trkPt)=(value.at(trkPt)>cut.at(trkPt));
	 pass.at(mumuMass)=(value.at(mumuMass)<cut.at(mumuMass)+phidM || value.at(mumuMass)>cut.at(mumuMass)-phidM);
    pass.at(fitVtxChiSq)=(value.at(fitVtxChiSq)<=cut.at(fitVtxChiSq));

    double wobs=1;
    double w;  

    if(!Ntp->isData()){w = 1;} //  No weights to data
    else{w=1;}

    bool status=AnalysisCuts(t,w,wobs);
    if(status){

    unsigned int mu1 = Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(0);
    unsigned int mu2 = Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1);
    unsigned int track = Ntp-> TwoMuonsTrackTrackIndex(tmp_idx).at(0);
    NVtx.at(t).Fill(Ntp->NVtx(),w);

    Track_Pt.at(t).Fill(Ntp->Track_P4(track).Pt(),w);
    Track_Eta.at(t).Fill(Ntp->Track_P4(track).Eta(),w);
    Track_Phi.at(t).Fill(Ntp->Track_P4(track).Phi(),w);
    Track_E.at(t).Fill(Ntp->Track_P4(track).E(),w);
    Track_P.at(t).Fill(Ntp->Track_P4(track).P(),w);
    Track_vx.at(t).Fill(Ntp->Track_Poca(track).X(),w);
    Track_vy.at(t).Fill(Ntp->Track_Poca(track).Y(),w);
    Track_vz.at(t).Fill(Ntp->Track_Poca(track).Z(),w);
    Track_normalizedChi2.at(t).Fill(Ntp->Track_normalizedChi2(track),w);
    Track_numberOfValidHits.at(t).Fill(Ntp->Track_numberOfValidHits(track),w);
    Track_charge.at(t).Fill(Ntp->Track_charge(track),w);
    Track_dxy.at(t).Fill(Ntp->Track_dxy(track),w);
    Track_dz.at(t).Fill(Ntp->Track_dz(track),w);
    Track_dxyError.at(t).Fill(Ntp->Track_dxyError(track),w);
    Track_dzError.at(t).Fill(Ntp->Track_dzError(track),w);

    Muon1_isGlobal.at(t).Fill(Ntp->Muon_isGlobalMuon(mu1),w);
    Muon2_isGlobal.at(t).Fill(Ntp->Muon_isGlobalMuon(mu2),w);
    Muon1_isStandAlone.at(t).Fill(Ntp->Muon_isStandAloneMuon(mu1),w);
    Muon2_isStandAlone.at(t).Fill(Ntp->Muon_isStandAloneMuon(mu2),w);
    Muon1_isTracker.at(t).Fill(Ntp->Muon_isTrackerMuon(mu1),w);
    Muon2_isTracker.at(t).Fill(Ntp->Muon_isTrackerMuon(mu2),w);
    Muon1_isCalo.at(t).Fill(Ntp->Muon_isCaloMuon(mu1),w);
    Muon2_isCalo.at(t).Fill(Ntp->Muon_isCaloMuon(mu2),w);
    Muon1_isIsolationValid.at(t).Fill(Ntp->Muon_isIsolationValid(mu1),w);
    Muon2_isIsolationValid.at(t).Fill(Ntp->Muon_isIsolationValid(mu2),w);
    Muon1_isTimeValid.at(t).Fill(Ntp->Muon_isTimeValid(mu1),w);
    Muon2_isTimeValid.at(t).Fill(Ntp->Muon_isTimeValid(mu2),w);
    Muon1_emEt03.at(t).Fill(Ntp->Muon_emEt03(mu1),w);
    Muon2_emEt03.at(t).Fill(Ntp->Muon_emEt03(mu2),w);
    Muon1_emVetoEt03.at(t).Fill(Ntp->Muon_emVetoEt03(mu1),w);
    Muon2_emVetoEt03.at(t).Fill(Ntp->Muon_emVetoEt03(mu2),w);
    Muon1_hadEt03.at(t).Fill(Ntp->Muon_hadEt03(mu1),w);
    Muon2_hadEt03.at(t).Fill(Ntp->Muon_hadEt03(mu2),w);
    Muon1_hadVetoEt03.at(t).Fill(Ntp->Muon_hadVetoEt03(mu1),w);
    Muon2_hadVetoEt03.at(t).Fill(Ntp->Muon_hadVetoEt03(mu2),w);
    Muon1_nJets03.at(t).Fill(Ntp->Muon_nJets03(mu1),w);
    Muon2_nJets03.at(t).Fill(Ntp->Muon_nJets03(mu2),w);
    Muon1_nTracks03.at(t).Fill(Ntp->Muon_nTracks03(mu1),w);
    Muon2_nTracks03.at(t).Fill(Ntp->Muon_nTracks03(mu2),w);
    Muon1_sumPt03.at(t).Fill(Ntp->Muon_sumPt03(mu1),w);
    Muon2_sumPt03.at(t).Fill(Ntp->Muon_sumPt03(mu2),w);
    Muon1_trackerVetoPt03.at(t).Fill(Ntp->Muon_trackerVetoPt03(mu1),w);
    Muon2_trackerVetoPt03.at(t).Fill(Ntp->Muon_trackerVetoPt03(mu2),w);
    Muon1_emEt05.at(t).Fill(Ntp->Muon_emEt05(mu1),w);
    Muon2_emEt05.at(t).Fill(Ntp->Muon_emEt05(mu2),w);
    Muon1_emVetoEt05.at(t).Fill(Ntp->Muon_emVetoEt05(mu1),w);
    Muon2_emVetoEt05.at(t).Fill(Ntp->Muon_emVetoEt05(mu2),w);
    Muon1_hadEt05.at(t).Fill(Ntp->Muon_hadEt05(mu1),w);
    Muon2_hadEt05.at(t).Fill(Ntp->Muon_hadEt05(mu2),w);
    Muon1_hadVetoEt05.at(t).Fill(Ntp->Muon_hadVetoEt05(mu1),w);
    Muon2_hadVetoEt05.at(t).Fill(Ntp->Muon_hadVetoEt05(mu2),w);
    Muon1_nJets05.at(t).Fill(Ntp->Muon_nJets05(mu1),w);
    Muon2_nJets05.at(t).Fill(Ntp->Muon_nJets05(mu2),w);
    Muon1_nTracks05.at(t).Fill(Ntp->Muon_nTracks05(mu1),w);
    Muon2_nTracks05.at(t).Fill(Ntp->Muon_nTracks05(mu2),w);
    Muon1_sumPt05.at(t).Fill(Ntp->Muon_sumPt05(mu1),w);
    Muon2_sumPt05.at(t).Fill(Ntp->Muon_sumPt05(mu2),w);
    Muon1_trackerVetoPt05.at(t).Fill(Ntp->Muon_trackerVetoPt05(mu1),w);
    Muon2_trackerVetoPt05.at(t).Fill(Ntp->Muon_trackerVetoPt05(mu2),w);
    Muon1_sumChargedHadronPt03.at(t).Fill(Ntp->Muon_sumChargedHadronPt03(mu1),w);
    Muon2_sumChargedHadronPt03.at(t).Fill(Ntp->Muon_sumChargedHadronPt03(mu2),w);
    Muon1_sumChargedParticlePt03.at(t).Fill(Ntp->Muon_sumChargedParticlePt03(mu1),w);
    Muon2_sumChargedParticlePt03.at(t).Fill(Ntp->Muon_sumChargedParticlePt03(mu2),w);
    Muon1_sumNeutralHadronEt03.at(t).Fill(Ntp->Muon_sumNeutralHadronEt03(mu1),w);
    Muon2_sumNeutralHadronEt03.at(t).Fill(Ntp->Muon_sumNeutralHadronEt03(mu2),w);
    Muon1_sumNeutralHadronEtHighThreshold03.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold03(mu1),w);
    Muon2_sumNeutralHadronEtHighThreshold03.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold03(mu2),w);
    Muon1_sumPhotonEt03.at(t).Fill(Ntp->Muon_sumPhotonEt03(mu1),w);
    Muon2_sumPhotonEt03.at(t).Fill(Ntp->Muon_sumPhotonEt03(mu2),w);
    Muon1_sumPhotonEtHighThreshold03.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold03(mu1),w);
    Muon2_sumPhotonEtHighThreshold03.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold03(mu2),w);
    Muon1_sumPUPt03.at(t).Fill(Ntp->Muon_sumPUPt03(mu1),w);
    Muon2_sumPUPt03.at(t).Fill(Ntp->Muon_sumPUPt03(mu2),w);
    Muon1_sumChargedHadronPt04.at(t).Fill(Ntp->Muon_sumChargedHadronPt04(mu1),w);
    Muon2_sumChargedHadronPt04.at(t).Fill(Ntp->Muon_sumChargedHadronPt04(mu2),w);
    Muon1_sumChargedParticlePt04.at(t).Fill(Ntp->Muon_sumChargedParticlePt04(mu1),w);
    Muon2_sumChargedParticlePt04.at(t).Fill(Ntp->Muon_sumChargedParticlePt04(mu2),w);
    Muon1_sumNeutralHadronEt04.at(t).Fill(Ntp->Muon_sumNeutralHadronEt04(mu1),w);
    Muon2_sumNeutralHadronEt04.at(t).Fill(Ntp->Muon_sumNeutralHadronEt04(mu2),w);
    Muon1_sumNeutralHadronEtHighThreshold04.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold04(mu1),w);
    Muon2_sumNeutralHadronEtHighThreshold04.at(t).Fill(Ntp->Muon_sumNeutralHadronEtHighThreshold04(mu2),w);
    Muon1_sumPhotonEt04.at(t).Fill(Ntp->Muon_sumPhotonEt04(mu1),w);
    Muon2_sumPhotonEt04.at(t).Fill(Ntp->Muon_sumPhotonEt04(mu2),w);
    Muon1_sumPhotonEtHighThreshold04.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold04(mu1),w);
    Muon2_sumPhotonEtHighThreshold04.at(t).Fill(Ntp->Muon_sumPhotonEtHighThreshold04(mu2),w);
    Muon1_sumPUPt04.at(t).Fill(Ntp->Muon_sumPUPt04(mu1),w);
    Muon2_sumPUPt04.at(t).Fill(Ntp->Muon_sumPUPt04(mu2),w);

    DimuondR.at(t).Fill(deltaR(Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(0)).Eta(),Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(0)).Phi(),Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(1)).Eta(),Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(1)).Phi()));
    Muon1TrkdR.at(t).Fill(deltaR(Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(0)).Eta(),Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(0)).Phi(),Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(pair_index).at(0)).Eta(),Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(pair_index).at(0)).Phi()));
    Muon2TrkdR.at(t).Fill(deltaR(Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(1)).Eta(),Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(1)).Phi(),Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(pair_index).at(0)).Eta(),Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(pair_index).at(0)).Phi()));
    MuonsPtRatio.at(t).Fill(Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(0)).Pt()/Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1)).Pt(),w );
    PhiMass.at(t).Fill((Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))).M(), w);
    TripleMass.at(t).Fill((Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))+ 
    Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(pair_index).at(0))).M(), w);
    double phimass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))).M();
    double dsmass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))+
    Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(pair_index).at(0))).M();
      PhiMassVsDsMass.at(t).Fill(phimass, dsmass);

      /*
       if(Ntp->NThreeMuons()!=0){
       unsigned int Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(0);
       unsigned int Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(1);
       unsigned int Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(2);

       double pt1 = Ntp->Muon_P4(Muon_index_1).Pt();
       double pt2 = Ntp->Muon_P4(Muon_index_2).Pt();
       double pt3 = Ntp->Muon_P4(Muon_index_3).Pt();
       FirstMuonsPt.at(t).Fill(pt1,1);
       SecondMuonsPt.at(t).Fill(pt2,1);
       ThirdMuonsPt.at(t).Fill(pt3,1);


       FirstMuonsEta.at(t).Fill( Ntp->Muon_P4(Muon_index_1).Eta(),1);
       SecondMuonsEta.at(t).Fill(  Ntp->Muon_P4(Muon_index_2).Eta(),1);
       ThirdMuonsEta.at(t).Fill(  Ntp->Muon_P4(Muon_index_3).Eta(),1);
       }


       double deltaMass(999.);
       unsigned int pair_index(0);
       if(Ntp->NTwoMuonsTrack()!=0){
       for(unsigned int i2M=0; i2M < Ntp->NTwoMuonsTrack(); i2M++){

       unsigned int muon_1 =  Ntp-> TwoMuonsTrackMuonIndices(i2M).at(0);
       unsigned int muon_2 =  Ntp-> TwoMuonsTrackMuonIndices(i2M).at(1);

       if( fabs((Ntp->Muon_P4(muon_1)  + Ntp->Muon_P4(muon_2)).M()  - PDG_Var::Phi_mass())< deltaMass){
       deltaMass =  fabs((Ntp->Muon_P4(muon_1)  + Ntp->Muon_P4(muon_2)).M()  - PDG_Var::Phi_mass());
       pair_index = i2M; // this is an index of the candidate with best mumu mass
       }
       }

       PhiMass.at(t).Fill((Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))).M(), w);
       TripleMass.at(t).Fill((Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))+ 
       Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(pair_index).at(0))).M(), w);

       double phimass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))).M();
       double dsmass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))+
       Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(pair_index).at(0))).M();

       PhiMassVsDsMass.at(t).Fill(phimass, dsmass);

       }*/
    }
}



void  MCEfficiency::Finish(){
    Selection::Finish();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
}





