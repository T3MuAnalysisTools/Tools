#include "DsToPhiPi.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"



DsToPhiPi::DsToPhiPi(TString Name_, TString id_):
  Selection(Name_,id_)
{


  // This is a class constructor;
}

DsToPhiPi::~DsToPhiPi(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }


  Logger(Logger::Info) << "complete." << std::endl;
}

void  DsToPhiPi::Configure(){

  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==HLTOk)           cut.at(HLTOk)=1;
    if(i==is2MuTrk)        cut.at(is2MuTrk)=1;
    if(i==PhiMassCut)      cut.at(PhiMassCut)=1;
    if(i==Mu1PtCut)        cut.at(Mu1PtCut)=1;
    if(i==Mu2PtCut)        cut.at(Mu2PtCut)=1;
    if(i==TrkPtCut)        cut.at(TrkPtCut)=1;
    if(i==Mu1DeltaR)       cut.at(Mu1DeltaR)=1;
    if(i==Mu2DeltaR)       cut.at(Mu2DeltaR)=1;
    if(i==NumValidHits)    cut.at(NumValidHits)=1;
  }

  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    // book here the N-1 and N-0 histrogramms for each cut
    if(i==is2MuTrk){
      title.at(i)="Category: 2Mu+Trk ";
      hlabel="2muon + track category";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_is2MuTrk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_is2MuTrk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==HLTOk){
      title.at(i)="HLT trigger ";
      hlabel="DoubleMu3_Trk_Tau3mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }	
    else if(i==PhiMassCut){
      title.at(i)= "0.95 $<$ $M_{\\mu\\mu}$ $<$ 1.1 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="OS muon pair mass, GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PhiMassCut_",htitle,40,0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PhiMassCut_",htitle,40,0.5,1.5,hlabel,"Events"));
    }
    else if(i==Mu1PtCut){
      title.at(i)= "Muon 1 Pt $>$ 2.5 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      hlabel="Pt of Muon 1";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1PtCut_",htitle,50,0,30,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1PtCut_",htitle,50,0,30,hlabel,"Events"));
    } 
    else if(i==Mu2PtCut){
      title.at(i)= "Muon 2 Pt $>$ 2.5 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      hlabel="Pt of Muon 2";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2PtCut_",htitle,50,0,30,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2PtCut_",htitle,50,0,30,hlabel,"Events"));
    } 
    else if(i==TrkPtCut){
      title.at(i)= "Track Pt $>$ 2.5 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      hlabel="Pt of Track";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TrkPtCut_",htitle,50,0,30,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TrkPtCut_",htitle,50,0,30,hlabel,"Events"));
    }
    else if(i==Mu1DeltaR){
      title.at(i)= "Muon 1 DeltaR no cut";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      hlabel="#Delta R of #mu_{1}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1DeltaR_",htitle,100,0,.1,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1DeltaR_",htitle,100,0,.1,hlabel,"Events"));
    }
    else if(i==Mu2DeltaR){
      title.at(i)= "Muon 2 DeltaR no cut";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      hlabel="#Delta R of #mu_{2}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2DeltaR_",htitle,100,0,.1,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2DeltaR_",htitle,100,0,.1,hlabel,"Events"));
    }
    else if(i==NumValidHits){
      title.at(i)= "Number of Valid Track Hits $>$ 7";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      hlabel="N Valid Track Hits";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NumValidHits_",htitle,35,0,35,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NumValidHits_",htitle,35,0,35,hlabel,"Events"));
    } 

  }
  // Track Candidate Information
  Track_P=HConfig.GetTH1D(Name+"_Track_P","Momentum magnitude of track (2mu+trk track candidate)",36,-0.5,35.5,"p (track)","Events");
  Track_E=HConfig.GetTH1D(Name+"_Track_E","Energy of track (2mu+trk track candidate)",36,-0.5,35.5,"E (track)","Events");
  Track_Pt=HConfig.GetTH1D(Name+"_Track_Pt","Transverse momentum of track (2mu+trk track candidate)",26,-0.5,25.5,"p_{T} (track)","Events");
  Track_Eta=HConfig.GetTH1D(Name+"_Track_Eta","Psuedorapidity of track (2mu+trk track candidate)",30,-2.5,2.5,"#eta","Events");
  Track_Phi=HConfig.GetTH1D(Name+"_Track_Phi","Azimuthal angle of track (2mu+trk track candidate)",30,-3.5,3.5,"#phi","Events");
  Track_normalizedChi2=HConfig.GetTH1D(Name+"_Track_normalizedChi2","Normalized chi square",20,-0.5,4.5,"#chi^{2} (track fit)","Events");
  Track_numberOfValidHits=HConfig.GetTH1D(Name+"_Track_numberOfValidHits","number of valid hits in te tracker",36,-0.5,35.5,"n valid track hits","Events");
  Track_charge=HConfig.GetTH1D(Name+"_Track_charge","Chargeof the track",3,-1.5,1.5,"Track charge","Events");
  
  Muon1_Pt=HConfig.GetTH1D(Name+"_Muon1_Pt","Transverse Pt (muon 1)",25,0,30,"#mu_{1} p_{T} (GeV)","Events");
  Muon1_Eta=HConfig.GetTH1D(Name+"_Muon1_Eta","Psuedorapidity (muon 1)",25,-2.5,2.5,"#mu_{1} #eta","Events");
  Muon1_Phi=HConfig.GetTH1D(Name+"_Muon1_Phi","Azimuthal angle of (muons 1)",25,-3.4,3.4,"#mu_{1} #phi","Events"); 
  Muon1_E=HConfig.GetTH1D(Name+"_Muon1_E","Energy of all (muon 1)",20,0,40,"#mu_{1} E (GeV)","Events");
  Muon1_P=HConfig.GetTH1D(Name+"_Muon1_P","Magnitude of momentum of (muon 1)",20,0,40,"#mu_{1} p (GeV)","Events");  
  
  Muon2_Pt=HConfig.GetTH1D(Name+"_Muon2_Pt","Transverse Pt (muon 2)",25,0,30,"#mu_{2} p_{T} (GeV)","Events");
  Muon2_Eta=HConfig.GetTH1D(Name+"_Muon2_Eta","Psuedorapidity (muon 2)",25,-2.5,2.5,"#mu_{2} #eta","Events");
  Muon2_Phi=HConfig.GetTH1D(Name+"_Muon2_Phi","Azimuthal angle of (muons 1)",25,-3.4,3.4,"#mu_{2} #phi","Events"); 
  Muon2_E=HConfig.GetTH1D(Name+"_Muon2_E","Energy of all (muon 2)",20,0,40,"#mu_{2} E (GeV)","Events");
  Muon2_P=HConfig.GetTH1D(Name+"_Muon2_P","Magnitude of momentum of (muon 2)",20,0,40,"#mu_{2} p (GeV)","Events");  
  
  Muon1_TriggerMatchdR=HConfig.GetTH1D(Name+"_Muon1_TriggerMatchdR","Trigger Matching mu1",50,0,0.05,"#Delta R Trigger Match #mu_{1}","Events");
  Muon2_TriggerMatchdR=HConfig.GetTH1D(Name+"_Muon2_TriggerMatchdR","Trigger Matching mu2",50,0,0.05,"#Delta R Trigger Match #mu_{2}","Events");
  Track_TriggerMatchdR=HConfig.GetTH1D(Name+"_Track_TriggerMatchdR","Trigger Matching track",50,0,0.05,"#Delta R Trigger Match track","Events");
  
  Muon1_isGlobal=HConfig.GetTH1D(Name+"_Muon1_isGlobal","Global muon status ",2,-.5,1.5,"#mu_{1} isGlb","Events");
  Muon2_isGlobal=HConfig.GetTH1D(Name+"_Muon2_isGlobal","Global muon status",2,-0.5,0.5,"#mu_{2} isGlb","Events");
  Muon1_isStandAlone=HConfig.GetTH1D(Name+"_Muon1_isStandAlone","",2,-0.5,1.5,"#mu_{1} isStandAlone","Events");
  Muon2_isStandAlone=HConfig.GetTH1D(Name+"_Muon2_isStandAlone","",2,-0.5,1.5,"#mu_{2} isStandAlone","Events");
  Muon1_isTracker=HConfig.GetTH1D(Name+"_Muon1_isTracker","",2,-0.5,1.5,"#mu_{1} isTracker","Events");
  Muon2_isTracker=HConfig.GetTH1D(Name+"_Muon2_isTracker","",2,-0.5,1.5,"#mu_{2} isTracker","Events");
  Muon1_isCalo=HConfig.GetTH1D(Name+"_Muon1_isCaloMuon","",2,-0.5,1.5,"#mu_{1} isCalo","Events");
  Muon2_isCalo=HConfig.GetTH1D(Name+"_Muon2_isCaloMuon","",2,-0.5,1.5,"#mu_{2} isCalo","Events");
  
  
  DimuondR=HConfig.GetTH1D(Name+"_DimuondR","dR between the muon pair",20,0,1,"dR","Events");
  Muon1TrkdR=HConfig.GetTH1D(Name+"_Muon1TrkdR","dR between the highest p muon and the track",100,0,5,"dR","Events");
  Muon2TrkdR=HConfig.GetTH1D(Name+"_Muon2TrkdR","dR between the lowest p muon and the track",100,0,5,"dR","Events");
  PhiMass=HConfig.GetTH1D(Name+"_PhiMass","#mu#mu invariant mass",50,0.2,1.5,"Mass of the #mu#mu pair","Events");
  PhiPlusTrackMass=HConfig.GetTH1D(Name+"_PhiPlusTrackMass","#mu#mu + track invariant mass",50,1.7,2.1,"Mass of the #mu#mu + track","Events");
  PhiMassVsDsMass=HConfig.GetTH2D(Name+"_PhiMassVsDsMass","#mu#mu invariant Mass vs. #mu#mu + track invariant mass",50,0.2,1.5,50,1.7,2.1,"M_{#mu#mu}, GeV","M_{#mu#mu + track}, GeV");
  
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms
  // Book here your analysis histrogramms, a good style is to follow selfexplanatory convention
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",66,-0.5,65.5,"Number of Vertices","Events");
 
  DsMass=HConfig.GetTH1D(Name+"_DsMass","Ds invariant mass",50,0,2.1,"M_{Ds} (GeV)", "Events"); 

  DsGenMatch=HConfig.GetTH1D(Name+"_DsGenMatch","dR between Gen Ds to Track",50,0,.1,"dR","Events");

  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}

void  DsToPhiPi::Store_ExtraDist(){ 
  
  Extradist1d.push_back(&Track_P);
  Extradist1d.push_back(&Track_E);
  Extradist1d.push_back(&Track_Pt);
  Extradist1d.push_back(&Track_Eta);
  Extradist1d.push_back(&Track_Phi);
  Extradist1d.push_back(&Track_normalizedChi2);
  Extradist1d.push_back(&Track_numberOfValidHits);
  Extradist1d.push_back(&Track_charge);
  
  Extradist1d.push_back(&Muon1_E);
  Extradist1d.push_back(&Muon1_Pt);
  Extradist1d.push_back(&Muon1_Phi);
  Extradist1d.push_back(&Muon1_Eta);
  Extradist1d.push_back(&Muon2_E);
  Extradist1d.push_back(&Muon2_Pt);
  Extradist1d.push_back(&Muon2_Phi);
  Extradist1d.push_back(&Muon2_Eta);
  Extradist1d.push_back(&DimuondR);
  Extradist1d.push_back(&Muon1TrkdR);
  Extradist1d.push_back(&Muon2TrkdR);
  Extradist1d.push_back(&PhiMass);
  Extradist1d.push_back(&PhiPlusTrackMass);
  Extradist2d.push_back(&PhiMassVsDsMass);
  Extradist1d.push_back(&Muon1_isGlobal);
  Extradist1d.push_back(&Muon2_isGlobal);
  Extradist1d.push_back(&Muon1_isStandAlone);
  Extradist1d.push_back(&Muon2_isStandAlone);
  Extradist1d.push_back(&Muon1_isTracker);
  Extradist1d.push_back(&Muon2_isTracker);
  Extradist1d.push_back(&Muon1_isCalo);
  Extradist1d.push_back(&Muon2_isCalo);
  Extradist1d.push_back(&Track_TriggerMatchdR);
  Extradist1d.push_back(&Muon1_TriggerMatchdR);
  Extradist1d.push_back(&Muon2_TriggerMatchdR);
  Extradist1d.push_back(&NVtx);
  Extradist1d.push_back(&DsMass);
  Extradist1d.push_back(&DsGenMatch);
	 
}


void  DsToPhiPi::doEvent(){ 
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection
  
  value.at(HLTOk) = 0;
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);
    if((HLT.Contains("DoubleMu3_Trk_Tau3mu") || HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu") ) && Ntp->HLTDecision(iTrigger) == 1)value.at(HLTOk)=Ntp->HLTDecision(iTrigger);
  }
    
  int mu1=-1, mu2=-1, track=-1;
  int tmp_idx = -1;
  double tmp_chisq = 999.0;
  
  value.at(is2MuTrk) = 0;
  value.at(PhiMassCut) = 0;
  value.at(Mu1PtCut) = 0;
  value.at(Mu2PtCut) = 0;
  value.at(TrkPtCut) = 0;
  value.at(Mu1DeltaR) = 0;
  value.at(Mu2DeltaR) = 0;
  value.at(NumValidHits) = 0;
  if(Ntp->NTwoMuonsTrack()!=0 && Ntp->NThreeMuons() == 0) value.at(is2MuTrk) = 1;
  
  if (value.at(is2MuTrk)==1){
    for(unsigned int i2M=0; i2M < Ntp->NTwoMuonsTrack(); i2M++){
      tmp_idx = i2M;
      int tmp_mu1 = Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0);
      int tmp_mu2 = Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(1);
      int tmp_track = Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0);
      if (tmp_chisq>Ntp->TwoMuonsTrack_SV_Chi2(i2M)){
	tmp_chisq = Ntp->TwoMuonsTrack_SV_Chi2(i2M);
	mu1 = tmp_mu1;
	mu2 = tmp_mu2;
	track = tmp_track;
      }
    }
    value.at(PhiMassCut) =  (Ntp->Muon_P4(mu1) + Ntp->Muon_P4(mu2)).M();
    value.at(Mu1PtCut) = Ntp->Muon_P4(mu1).Pt();
    value.at(Mu2PtCut) = Ntp->Muon_P4(mu2).Pt();
    value.at(TrkPtCut) = Ntp->Track_P4(track).Pt();
    value.at(Mu1DeltaR) = Ntp->Muon_P4(mu1).DeltaR(Ntp->Track_P4(track));
    value.at(Mu2DeltaR) = Ntp->Muon_P4(mu2).DeltaR(Ntp->Track_P4(track));
    value.at(NumValidHits) = Ntp->Track_numberOfValidHits(track);
  }

  pass.at(PhiMassCut) = ( value.at(PhiMassCut) > 0.95 && value.at(PhiMassCut) < 1.1);
  pass.at(is2MuTrk) = (value.at(is2MuTrk)==cut.at(is2MuTrk));
  pass.at(HLTOk)= (value.at(HLTOk)==cut.at(HLTOk)); 
  pass.at(Mu1PtCut)= (value.at(Mu1PtCut) > 2.5);
  pass.at(Mu2PtCut)= (value.at(Mu2PtCut) > 2.5);
  pass.at(TrkPtCut)= (value.at(TrkPtCut) > 2.5);
  pass.at(Mu1DeltaR)= true;
  pass.at(Mu2DeltaR)= true;
  pass.at(NumValidHits)= (value.at(NumValidHits) > 7);

  double wobs=1;
  double w;  
             
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  bool status=AnalysisCuts(t,w,wobs);
  if(status){
    NVtx.at(t).Fill(Ntp->NVtx(),w);

    Track_Pt.at(t).Fill(Ntp->Track_P4(track).Pt(),w);
    Track_Eta.at(t).Fill(Ntp->Track_P4(track).Eta(),w);
    Track_Phi.at(t).Fill(Ntp->Track_P4(track).Phi(),w);
    Track_P.at(t).Fill(Ntp->Track_P4(track).P(),w);

    Track_normalizedChi2.at(t).Fill(Ntp->Track_normalizedChi2(track),w);
    Track_numberOfValidHits.at(t).Fill(Ntp->Track_numberOfValidHits(track),w);
    Track_charge.at(t).Fill(Ntp->Track_charge(track),w);
    Muon1_isGlobal.at(t).Fill(Ntp->Muon_isGlobalMuon(mu1),w);
    Muon2_isGlobal.at(t).Fill(Ntp->Muon_isGlobalMuon(mu2),w);
    Muon1_isStandAlone.at(t).Fill(Ntp->Muon_isStandAloneMuon(mu1),w);
    Muon2_isStandAlone.at(t).Fill(Ntp->Muon_isStandAloneMuon(mu2),w);
    Muon1_isTracker.at(t).Fill(Ntp->Muon_isTrackerMuon(mu1),w);
    Muon2_isTracker.at(t).Fill(Ntp->Muon_isTrackerMuon(mu2),w);
 

    Muon1_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(0),w);
    Muon2_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(1),w);
    Track_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(2),w);


    Muon1_Pt.at(t).Fill(Ntp->Muon_P4(mu1).Pt(),w);
    Muon1_Eta.at(t).Fill(Ntp->Muon_P4(mu1).Eta(),w);
    Muon1_Phi.at(t).Fill(Ntp->Muon_P4(mu1).Phi(),w);
    Muon1_E.at(t).Fill(Ntp->Muon_P4(mu1).E(),w);


    Muon2_Pt.at(t).Fill(Ntp->Muon_P4(mu2).Pt(),w);
    Muon2_Eta.at(t).Fill(Ntp->Muon_P4(mu2).Eta(),w);
    Muon2_Phi.at(t).Fill(Ntp->Muon_P4(mu2).Phi(),w);
    Muon2_E.at(t).Fill(Ntp->Muon_P4(mu2).E(),w);

    
    DimuondR.at(t).Fill(Ntp->Muon_P4(mu1).DeltaR(Ntp->Muon_P4(mu2)),w);
    Muon1TrkdR.at(t).Fill(Ntp->Muon_P4(mu1).DeltaR(Ntp->Track_P4(track)),w);
    Muon2TrkdR.at(t).Fill(Ntp->Muon_P4(mu2).DeltaR(Ntp->Track_P4(track)),w);
    
    PhiMass.at(t).Fill((Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))).M(), w);
    PhiPlusTrackMass.at(t).Fill((Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))+ 
			   Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0))).M(), w);
    
    double phimass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))).M();
    double dsmass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))+
		     Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0))).M();
    PhiMassVsDsMass.at(t).Fill(phimass, dsmass);

    DsMass.at(t).Fill(dsmass);

    if(id==30){

      DsGenMatch.at(t).Fill(Ntp->DsGenMatch(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0) + Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(1) + Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0)));

    }

  }
}

void  DsToPhiPi::Finish(){

  if(mode == RECONSTRUCT){
    for(unsigned int i=1; i<  Nminus0.at(0).size(); i++){
      double scale(1.);
      if(Nminus0.at(0).at(i).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(i).Integral()/1;
      ScaleAllHistOfType(i,scale);
    }
  }
  Selection::Finish();

}



