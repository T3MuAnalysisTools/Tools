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
  PhiMassHigh(1.04)
{
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
	 if(i==OSMuons)				cut.at(OSMuons)=0.1;
    if(i==Mu1PtCut)           cut.at(Mu1PtCut)=3.0;
    if(i==Mu2PtCut)           cut.at(Mu2PtCut)=3.0;
    if(i==MuonID)             cut.at(MuonID)=1;
    if(i==MuMuMassCut)			cut.at(MuMuMassCut)=1;
    if(i==TrackPtCut)			cut.at(TrackPtCut)=2.0;
    if(i==NTrackHits)			cut.at(NTrackHits)=6;
    if(i==ChiSqCut)				cut.at(ChiSqCut)=15;
    if(i==TriggerMatch)       cut.at(TriggerMatch)=1;
	 if(i==GenMatch)	  		  	cut.at(GenMatch)=0.03;
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
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TwoMuTrkCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TwoMuTrkCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
      else if(i==Mu1PtCut){
        title.at(i)="Mu1 Pt $>$ ";
		  title.at(i)+=cut.at(Mu1PtCut);
		  title.at(i)+=" GeV";
        hlabel="Muon1 PT, GeV";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1PtCut_",htitle,40,2,25,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1PtCut_",htitle,40,2,25,hlabel,"Events"));
      }
		else if(i==OSMuons){
        title.at(i)="abs(Mu1+Mu2 Charge) < 0.1 ";
        hlabel="Sum of muon charges";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OSMuons_",htitle,3,-1.5,1.5,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OSMuons_",htitle,3,-1.5,1.5,hlabel,"Events"));
      }
      else if(i==Mu2PtCut){
        title.at(i)="Mu2 Pt $>$ ";
		  title.at(i)+=cut.at(Mu2PtCut);
		  title.at(i)+=" GeV";
        hlabel="Muon2 PT, GeV";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2PtCut_",htitle,40,2,20,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2PtCut_",htitle,40,2,20,hlabel,"Events"));
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
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuMuMassCut_",htitle,40,2,20,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuMuMassCut_",htitle,40,2,20,hlabel,"Events"));
      }
      else if(i==TrackPtCut){
        title.at(i)="Trk Pt $>$ ";
		  title.at(i)+=cut.at(TrackPtCut);
		  title.at(i)+=" GeV";
        hlabel="Track PT, GeV";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TrackPtCut_",htitle,40,2,15,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TrackPtCut_",htitle,40,2,15,hlabel,"Events"));
      }
		else if(i==NTrackHits){
		  title.at(i)="Number of hits in tracker $>$";
		  title.at(i)+=cut.at(NTrackHits);
		  hlabel="N Tracker Hits";
		  Nminus1.push_back(HConfig.GetTH1D(Name+"_Nminus1_NTrackHits_",htitle,10,0,10,hlabel,"Events"));
		  Nminus0.push_back(HConfig.GetTH1D(Name+"_Nminus0_NTrackHits_",htitle,10,0,10,hlabel,"Events"));
		}
      else if(i==ChiSqCut){
		  title.at(i)="Chi Sq of vertex fitting $<$ ";
		  title.at(i)+=cut.at(ChiSqCut);
		  hlabel="Chi Sq of vertex fitting";
		  Nminus1.push_back(HConfig.GetTH1D(Name+"_Nminus1_ChiSq_",htitle,10,0,10,hlabel,"Events"));
		  Nminus0.push_back(HConfig.GetTH1D(Name+"_Nminus0_ChiSq_",htitle,10,0,10,hlabel,"Events"));
		}
      else if(i==TriggerMatch){
        title.at(i)="Trigger Matching";
        hlabel="Sum of dR_{reco-trigger}";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatch_",htitle,40,0,0.05,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatch_",htitle,40,0,0.05,hlabel,"Events"));
      }
		else if(i==GenMatch){
        title.at(i)="GEN matching";
        hlabel="GEN match";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_GENMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_GENMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
    }

    // Track Candidate Information
    Track_P=HConfig.GetTH1D(Name+"_Track_P","Momentum magnitude of track (2mu+trk track candidate)",66,-0.5,65.5,"p (track)","Events");
    Track_E=HConfig.GetTH1D(Name+"_Track_E","Energy of track (2mu+trk track candidate)",66,-0.5,65.5,"E (track)","Events");
    Track_Pt=HConfig.GetTH1D(Name+"_Track_Pt","Transverse momentum of track (2mu+trk track candidate)",66,-0.5,65.5,"p_{T} (track)","Events");
    Track_Eta=HConfig.GetTH1D(Name+"_Track_Eta","Psuedorapidity of track (2mu+trk track candidate)",66,-0.5,65.5,"#eta","Events");
    Track_Phi=HConfig.GetTH1D(Name+"_Track_Phi","Azimuthal angle of track (2mu+trk track candidate)",66,-0.5,65.5,"#phi","Events");
    Track_vx=HConfig.GetTH1D(Name+"_Track_vx","X coordinate of the parent vertex (2mu+trk track candidate)",66,-0.5,65.5,"Parent vertex x coordinate (cm)","Events");
    Track_vy=HConfig.GetTH1D(Name+"_Track_vy","Y coordinate of the parent vertex (2mu+trk track candidate)",66,-0.5,65.5,"Parent vertex y coordinate (cm)","Events");
    Track_vz=HConfig.GetTH1D(Name+"_Track_vz","Z coordinate of the parent vertex (2mu+trk track candidate)",66,-0.5,65.5,"Parent vertex z coordinate (cm)","Events");
    Track_normalizedChi2=HConfig.GetTH1D(Name+"_Track_normalizedChi2","Normalized chi square",66,-0.5,65.5,"#chi^{2}/ndf (track fit)","Events");
    Track_numberOfValidHits=HConfig.GetTH1D(Name+"_Track_numberOfValidHits","number of valid hits in the tracker",66,-0.5,65.5,"n valid track hits","Events");
    Track_charge=HConfig.GetTH1D(Name+"_Track_charge","Chargeof the track",66,-0.5,65.5,"Number of Vertices","Events");
    Track_dxy=HConfig.GetTH1D(Name+"_Track_dxy","Transverse displacement of the parent vertex from the bs",66,-0.5,65.5,"dxy (cm)","Events");
    Track_dz=HConfig.GetTH1D(Name+"_Track_dz","Longitudnal displacement of the parent vertex from the bs",66,-0.5,65.5,"dz (cm)","Events");
    Track_dxyError=HConfig.GetTH1D(Name+"_Track_dxyError","dxy Error",66,-0.5,65.5,"#Deltadxy (cm)","Events");
    Track_dzError=HConfig.GetTH1D(Name+"_Track_dzError","dz Error",66,-0.5,65.5,"#Deltadz (cm)","Events");

    // Muon variables (Muons from dimuon + track candidates)
    Muon1_Pt=HConfig.GetTH1D(Name+"_Muon1_Pt","Transverse Pt (muon 1)",25,0,50,"#mu_{1} p_{T} (GeV)","Events");
    Muon1_Eta=HConfig.GetTH1D(Name+"_Muon1_Eta","Psuedorapidity (muon 1)",25,-2.5,2.5,"#mu_{1} #eta","Events");
    Muon1_Phi=HConfig.GetTH1D(Name+"_Muon1_Phi","Azimuthal angle of (muons 1)",25,-3.15,3.15,"#mu_{1} #phi","Events"); 
    Muon1_E=HConfig.GetTH1D(Name+"_Muon1_E","Energy of all (muon 1)",20,0,40,"#mu_{1} E (GeV)","Events");
    Muon1_P=HConfig.GetTH1D(Name+"_Muon1_P","Magnitude of momentum of (muon 1)",20,0,40,"#mu_{1} p (GeV)","Events");  

    Muon1_vx=HConfig.GetTH1D(Name+"_Muon1_Vx","X coordinate of the parent vertex all muons",100,0,5,"#mu_{1} vx","Events"); 
    Muon1_vy=HConfig.GetTH1D(Name+"_Muon1_Vy","Y coordinate of the parent vertex all muons",100,0,5,"#mu_{1} vy","Events"); 
    Muon1_vz=HConfig.GetTH1D(Name+"_Muon1_Vz","Z coordinate of the parent vertex all muons",100,0,5,"#mu_{1} vz","Events");

    Muon2_Pt=HConfig.GetTH1D(Name+"_Muon2_Pt","Transverse Pt (muon 2)",25,0,50,"#mu_{2} p_{T} (GeV)","Events");
    Muon2_Eta=HConfig.GetTH1D(Name+"_Muon2_Eta","Psuedorapidity (muon 2)",25,-2.5,2.5,"#mu_{2} #eta","Events");
    Muon2_Phi=HConfig.GetTH1D(Name+"_Muon2_Phi","Azimuthal angle of (muons 1)",25,-3.15,3.15,"#mu_{2} #phi","Events"); 
    Muon2_E=HConfig.GetTH1D(Name+"_Muon2_E","Energy of all (muon 2)",20,0,40,"#mu_{2} E (GeV)","Events");
    Muon2_P=HConfig.GetTH1D(Name+"_Muon2_P","Magnitude of momentum of (muon 2)",20,0,40,"#mu_{2} p (GeV)","Events");  
    Muon2_vx=HConfig.GetTH1D(Name+"_Muon2_Vx","X coordinate of the parent vertex all muons",100,0,5,"#mu_{2} vx","Events"); 
    Muon2_vy=HConfig.GetTH1D(Name+"_Muon2_Vy","Y coordinate of the parent vertex all muons",100,0,5,"#mu_{2} vy","Events"); 
    Muon2_vz=HConfig.GetTH1D(Name+"_Muon2_Vz","Z coordinate of the parent vertex all muons",100,0,5,"#mu_{2} vz","Events");

    Muon1_isGlobal=HConfig.GetTH1D(Name+"_Muon1_isGlobal","Global muons status ",2,-.5,1.5,"#mu_{1} isGlb","Events");
    Muon2_isGlobal=HConfig.GetTH1D(Name+"_Muon2_isGlobal","",2,-0.5,0.5,"#mu_{2} isGlb","Events");
    Muon1_isStandAlone=HConfig.GetTH1D(Name+"_Muon1_isStandAlone","",2,-0.5,1.5,"#mu_{1} isStandAlone","Events");
    Muon2_isStandAlone=HConfig.GetTH1D(Name+"_Muon2_isStandAlone","",2,-0.5,1.5,"#mu_{2} isStandAlone","Events");
    Muon1_isTracker=HConfig.GetTH1D(Name+"_Muon1_isTracker","",2,-0.5,1.5,"#mu_{1} isTracker","Events");
    Muon2_isTracker=HConfig.GetTH1D(Name+"_Muon2_isTracker","",2,-0.5,1.5,"#mu_{2} isTracker","Events");
    Muon1_isCalo=HConfig.GetTH1D(Name+"_Muon1_isCaloMuon","",2,-0.5,1.5,"#mu_{1} isCalo","Events");
    Muon2_isCalo=HConfig.GetTH1D(Name+"_Muon2_isCaloMuon","",2,-0.5,1.5,"#mu_{2} isCalo","Events");
    Muon1_isIsolationValid=HConfig.GetTH1D(Name+"_Muon1_isIsolationValid","#mu_{1} isIsoValid",2,-0.5,1.5,"#mu_{1} isIsolationValid","Events");
    Muon2_isIsolationValid=HConfig.GetTH1D(Name+"_Muon2_isIsolationValid","#mu_{2} isIsoValid",2,-0.5,1.5,"#mu_{2} isIsolationValid","Events");
    Track_TriggerMatchdR=HConfig.GetTH1D(Name+"_Track_TriggerMatchdR","track dR (trigger match)",10,-0.5,9.5,"track dR (trigger match)","Events");
    Muon1_TriggerMatchdR=HConfig.GetTH1D(Name+"_Muon1_TriggerMatchdR","#mu_{1} dR (trigger match)",10,-0.5,9.5,"#mu_{1} dR (trigger match)","Events");
    Muon2_TriggerMatchdR=HConfig.GetTH1D(Name+"_Muon2_TriggerMatchdR","#mu_{2} dR (trigger match)",10,-0.5,9.5,"#mu_{2} dR (trigger match)","Events");

    //Dimuon Information (Muons from dimuon + track candidates)
    DimuondR=HConfig.GetTH1D(Name+"_DimuondR","dR between the muon pair",20,0,1,"dR","Events");
    Muon1TrkdR=HConfig.GetTH1D(Name+"_Muon1TrkdR","dR between the highest p muon and the track",100,0,5,"dR (#mu_{1},track)","Events");
    Muon2TrkdR=HConfig.GetTH1D(Name+"_Muon2TrkdR","dR between the lowest p muon and the track",100,0,5,"dR (#mu_{2},track)","Events");
    PhiMass=HConfig.GetTH1D(Name+"_PhiMass","#mu#mu mass",50,0.77,1.27,"Invariant mass of the #mu#mu pair","Events");
    TripleMass=HConfig.GetTH1D(Name+"_TripleMass","#mu#mu + track mass",42,1.68,2.1,"Invariant mass of the #mu#mu + track","Events");
    PhiMassVsDsMass=HConfig.GetTH2D(Name+"_PhiMassVsDsMass","#mu#mu Invariant mass vs. #mu#mu + track Invariant mass",50,0.77,1.27,42,1.68,2.1,"M_{#mu#mu}, GeV","M_{#mu#mu + track}, GeV");

    // Setup NPassed Histogams
    Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
    // Setup Extra Histograms
    // Book here your analysis histrogramms, a good style is to follow selfexplanatory convention
    NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",66,-0.5,65.5,"Number of Vertices","Events");

    Selection::ConfigureHistograms(); //do not remove
    HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}

void  DsPhiPeak::Store_ExtraDist(){ 
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

    //Dimuon variable
    Extradist1d.push_back(&Muon1_P);
    Extradist1d.push_back(&Muon1_E);
    Extradist1d.push_back(&Muon1_Pt);
    Extradist1d.push_back(&Muon1_Phi);
    Extradist1d.push_back(&Muon1_Eta);
    Extradist1d.push_back(&Muon1_vx);
    Extradist1d.push_back(&Muon1_vy);
    Extradist1d.push_back(&Muon1_vz);
    Extradist1d.push_back(&Muon2_P);
    Extradist1d.push_back(&Muon2_E);
    Extradist1d.push_back(&Muon2_Pt);
    Extradist1d.push_back(&Muon2_Phi);
    Extradist1d.push_back(&Muon2_Eta);
    Extradist1d.push_back(&Muon2_vx);
    Extradist1d.push_back(&Muon2_vy);
    Extradist1d.push_back(&Muon2_vz);
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
    Extradist1d.push_back(&Track_TriggerMatchdR);
    Extradist1d.push_back(&Muon1_TriggerMatchdR);
    Extradist1d.push_back(&Muon2_TriggerMatchdR);

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // Here you must push back all analysis histograms, otherwise they wont be propagated to the output //
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    Extradist1d.push_back(&NVtx);
}

void  DsPhiPeak::doEvent(){ 
  unsigned int t;
  int id(Ntp->GetMCID());
  bool DEBUG = false;
  bool HLTOk(false);
  bool L1Ok(false);
  bool DoubleMuFired(false);
  bool TripleMuFired(false);

  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}

  // Apply Selection
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
      TString HLT = Ntp->HLTName(iTrigger);
      if( HLT.Contains("DoubleMu3_Trk_Tau3mu") && Ntp->HLTDecision(iTrigger)==1) HLTOk = true;
      if( HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu") && Ntp->HLTDecision(iTrigger)==1) HLTOk = true;
    }
  
  for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
    TString L1TriggerName = Ntp->L1Name(il1);
    
    if(id==1 && Ntp->WhichEra(2017).Contains("RunB")){
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
      if(L1TriggerName.Contains("L1_TripleMu_5_3_0_DoubleMu_5_3_OS_Mass_Max17"))      TripleMuFired = Ntp-> L1Decision(il1);
    }

    if(id!=1){
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
      if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9"))TripleMuFired = Ntp-> L1Decision(il1);
    }

    if(id==1 && (Ntp->WhichEra(2017).Contains("RunC") or Ntp->WhichEra(2017).Contains("RunD") or Ntp->WhichEra(2017).Contains("RunF")  or Ntp->WhichEra(2017).Contains("RunE"))){
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
      if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9"))TripleMuFired = Ntp-> L1Decision(il1);
    }

  }

	 if (DoubleMuFired || TripleMuFired) L1Ok = true;
	 if (L1Ok && HLTOk) value.at(TriggerOk) = true;
	 else value.at(TriggerOk) = false;

    unsigned int final_idx = 0;
    double tmp_chiSq = 999.0;

	 for (int it=0; it<Ntp->NTwoMuonsTrack(); it++){
	   if (Ntp->TwoMuonsTrack_SV_Chi2(it)<tmp_chiSq){
		  tmp_chiSq = Ntp->TwoMuonsTrack_SV_Chi2(it);
		  final_idx = it;
		}
	 }

	 value.at(TwoMuTrkCandidate) = Ntp->NTwoMuonsTrack();
	 value.at(OSMuons) = 1;
	 if (Ntp->NTwoMuonsTrack()>0){

		unsigned int Muon1_idx = Ntp->TwoMuonsTrackMuonIndices(final_idx).at(0);
		unsigned int Muon2_idx = Ntp->TwoMuonsTrackMuonIndices(final_idx).at(1);
		unsigned int Track_idx = Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0);

		unsigned int Mu1_pt_idx = (Ntp->Muon_P4(Muon1_idx).Pt() > Ntp->Muon_P4(Muon2_idx).Pt() ? Muon1_idx:Muon2_idx);
		unsigned int Mu2_pt_idx = (Ntp->Muon_P4(Muon1_idx).Pt() > Ntp->Muon_P4(Muon2_idx).Pt() ? Muon2_idx:Muon1_idx);

		TLorentzVector mu1_p4 = Ntp->Muon_P4(Mu1_pt_idx);
		TLorentzVector mu2_p4 = Ntp->Muon_P4(Mu2_pt_idx);
		TLorentzVector track_p4 = Ntp->Track_P4(Track_idx);

		value.at(Mu1PtCut) = mu1_p4.Pt();
		value.at(Mu2PtCut) = mu2_p4.Pt();
      value.at(MuonID) = (Ntp->Muon_isGlobalMuon(Mu1_pt_idx) && Ntp->Muon_isGlobalMuon(Mu2_pt_idx) );
		value.at(OSMuons) = abs(Ntp->Muon_charge(Mu1_pt_idx)+Ntp->Muon_charge(Mu2_pt_idx));
	   
		float PhiMass = (mu1_p4+mu2_p4).M(); 
		
		value.at(MuMuMassCut) = ( (PhiMass < PhiMassHigh) && (PhiMass >= PhiMassLow) );
		value.at(TrackPtCut) = track_p4.Pt();
      value.at(NTrackHits) = Ntp->Track_numberOfValidHits(Track_idx);
		value.at(ChiSqCut) = Ntp->TwoMuonsTrack_SV_Chi2(final_idx);
      bool triggerMatch = true;
		for (auto& i: Ntp->TwoMuonsTrack_TriggerMatch_dR(final_idx)){
		   if (i>0.03) { triggerMatch = false; break; }
		}
		value.at(TriggerMatch) = triggerMatch;
		if (id==1) value.at(GenMatch) = 0;
		else value.at(GenMatch) = Ntp->DsGenMatch(final_idx); 
	 }
   
	 pass.at(TriggerOk) = ( value.at(TriggerOk) == cut.at(TriggerOk) );
	 pass.at(TwoMuTrkCandidate) = ( value.at(TwoMuTrkCandidate) >= cut.at(TwoMuTrkCandidate) );
    pass.at(OSMuons) = value.at(OSMuons) < cut.at(OSMuons);
	 pass.at(Mu1PtCut) = (value.at(Mu1PtCut) >= cut.at(Mu1PtCut) );
	 pass.at(Mu2PtCut) = (value.at(Mu2PtCut) >= cut.at(Mu2PtCut) );
	 pass.at(MuonID) = ( value.at(MuonID) == cut.at(MuonID) );
	 pass.at(MuMuMassCut) = ( value.at(MuMuMassCut) == cut.at(MuMuMassCut) );
	 pass.at(TrackPtCut) = ( value.at(TrackPtCut) > cut.at(TrackPtCut) );
	 pass.at(NTrackHits) = ( value.at(NTrackHits) > cut.at(NTrackHits) );
	 pass.at(ChiSqCut) = ( value.at(ChiSqCut) < cut.at(ChiSqCut) );
	 pass.at(TriggerMatch) = ( value.at(TriggerMatch) == cut.at(TriggerMatch) );
    pass.at(GenMatch) = ( value.at(GenMatch) < cut.at(GenMatch) );
    
    double wobs=1;
    double w;

    if(!Ntp->isData()){w = 1;} //  No weights to data
    else{w=1;}
    bool status=AnalysisCuts(t,w,wobs);
    if(status){
	 	
		unsigned int mu1_idx = Ntp->TwoMuonsTrackMuonIndices(final_idx).at(0);
      unsigned int mu2_idx = Ntp->TwoMuonsTrackMuonIndices(final_idx).at(1);
      unsigned int track = Ntp->TwoMuonsTrackTrackIndex(final_idx).at(0);

		unsigned int mu1 = (Ntp->Muon_P4(mu1_idx).Pt()>Ntp->Muon_P4(mu2_idx).Pt()?mu1_idx:mu2_idx); 
		unsigned int mu2 = (Ntp->Muon_P4(mu1_idx).Pt()>Ntp->Muon_P4(mu2_idx).Pt()?mu2_idx:mu1_idx); 

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

      Muon1_E.at(t).Fill(Ntp->Muon_P4(mu1).E(),w);
      Muon1_P.at(t).Fill(Ntp->Muon_P4(mu1).P(),w);
      Muon1_Pt.at(t).Fill(Ntp->Muon_P4(mu1).Pt(),w);
      Muon1_Eta.at(t).Fill(Ntp->Muon_P4(mu1).Eta(),w);
      Muon1_Phi.at(t).Fill(Ntp->Muon_P4(mu1).Phi(),w);
      Muon1_vx.at(t).Fill(Ntp->Muon_Poca(mu1).X(),w);
      Muon1_vy.at(t).Fill(Ntp->Muon_Poca(mu1).Y(),w);
      Muon1_vz.at(t).Fill(Ntp->Muon_Poca(mu1).Z(),w);

      Muon2_E.at(t).Fill(Ntp->Muon_P4(mu2).E(),w);
      Muon2_P.at(t).Fill(Ntp->Muon_P4(mu2).P(),w);
      Muon2_Pt.at(t).Fill(Ntp->Muon_P4(mu2).Pt(),w);
      Muon2_Eta.at(t).Fill(Ntp->Muon_P4(mu2).Eta(),w);
      Muon2_Phi.at(t).Fill(Ntp->Muon_P4(mu2).Phi(),w);
      Muon2_vx.at(t).Fill(Ntp->Muon_Poca(mu2).X(),w);
      Muon2_vy.at(t).Fill(Ntp->Muon_Poca(mu2).Y(),w);
      Muon2_vz.at(t).Fill(Ntp->Muon_Poca(mu2).Z(),w);

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

      Muon1_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(final_idx)).at(0),w);
      Muon2_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(final_idx)).at(1),w);
      Track_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(final_idx)).at(2),w);

		TLorentzVector mu1_p4 = Ntp->Muon_P4(mu1);
		TLorentzVector mu2_p4 = Ntp->Muon_P4(mu2);
		TLorentzVector track_p4 = Ntp->Track_P4(track);

      DimuondR.at(t).Fill(mu1_p4.DeltaR(mu2_p4));
      Muon1TrkdR.at(t).Fill(mu1_p4.DeltaR(track_p4));
      Muon2TrkdR.at(t).Fill(mu2_p4.DeltaR(track_p4));

      double phimass = (mu1_p4+mu2_p4).M();
      double dsmass = (mu1_p4+mu2_p4+track_p4).M();

		PhiMass.at(t).Fill(phimass, w);
      TripleMass.at(t).Fill(dsmass, w);
      PhiMassVsDsMass.at(t).Fill(phimass, dsmass);
    }
}

void  DsPhiPeak::Finish(){
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /* 
  if(mode == RECONSTRUCT){
    //    for(unsigned int i=1; i<  Nminus0.at(0).size(); i++){
    int id(Ntp->GetMCID());
    double scale(1.);
    double scaleDsTau(0.637);
    double scaleBpTau(0.262);
    double scaleB0Tau(0.099);

    if(Nminus0.at(0).at(2).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(2).Integral();
    ScaleAllHistOfType(2,scale*scaleDsTau);
    
    if(Nminus0.at(0).at(3).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(3).Integral();
    ScaleAllHistOfType(3,scale*scaleB0Tau);

    if(Nminus0.at(0).at(4).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(4).Integral();
    ScaleAllHistOfType(4,scale*scaleBpTau);

    //    }
  }
  */
    Selection::Finish();
}
