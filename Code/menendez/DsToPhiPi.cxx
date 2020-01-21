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
    if(i==L1TOk)           cut.at(L1TOk)=1;
    if(i==HLTOk)           cut.at(HLTOk)=1;
    if(i==is2MuTrk)        cut.at(is2MuTrk)=1;
    if(i==GlobalMu)        cut.at(GlobalMu)=1;
    if(i==MuCharge)        cut.at(MuCharge)=1;
    if(i==Mass2Mu)         cut.at(Mass2Mu)=1;
    if(i==Mu1dR)           cut.at(Mu1dR)=1;
    if(i==Mu2dR)           cut.at(Mu2dR)=1;
    if(i==TrkdR)           cut.at(TrkdR)=1;
    if(i==Mu1pt)           cut.at(Mu1pt)=1;
    if(i==Mu2pt)           cut.at(Mu2pt)=1;
    if(i==Trkpt)           cut.at(Trkpt)=1;

  }

  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    // book here the N-1 and N-0 histrogramms for each cut
    if(i==L1TOk){
      title.at(i)="L1T trigger ";
      hlabel="Level 1 Trigger";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_L1TOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_L1TOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==HLTOk){
      title.at(i)="HLT trigger ";
      hlabel="DoubleMu3_Trk_Tau3mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==is2MuTrk){
      title.at(i)="Category: 2Mu+Trk ";
      hlabel="2muon + track category";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_is2MuTrk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_is2MuTrk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==GlobalMu){
      title.at(i)="Muons are Global";
      hlabel="Muons are Global";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_GlobalMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_GlobalMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
    } 
    else if(i==MuCharge){
      title.at(i)="Muons opposite charge";
      hlabel="Muons have opposite charge";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Mass2Mu){
      title.at(i)="1.00 $<$ $M_{2\\mu}$ $<$ 1.04 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Invariant mass of 2 muons";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mass2Mu_",htitle,200,0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mass2Mu_",htitle,200,0.5,1.5,hlabel,"Events"));
    }
    else if(i==Mu1dR){
      title.at(i)="Mu01 dRtriggerMatch $<$ 0.03";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      hlabel="Trigger Match dR of Muon 1";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1dR_",htitle,50,0,.05,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1dR_",htitle,50,0,.05,hlabel,"Events"));
    }
    else if(i==Mu2dR){
      title.at(i)="Mu02 dRtriggerMatch $<$ 0.03";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      hlabel="Trigger Match dR of Muon 2";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2dR_",htitle,50,0,.05,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2dR_",htitle,50,0,.05,hlabel,"Events"));
    }
    else if(i==TrkdR){
      title.at(i)="Tr dRtriggerMatch $<$ 0.03";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      hlabel="Trigger Match dR of Track";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TrkdR_",htitle,50,0,.05,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TrkdR_",htitle,50,0,.05,hlabel,"Events"));
    }
    else if(i==Mu1pt){
      title.at(i)="$\\mu_{1}$ Pt $>$ 3 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Pt of Muon 1";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1pt_",htitle,80,0,20,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1pt_",htitle,80,0,20,hlabel,"Events"));
    }
    else if(i==Mu2pt){
      title.at(i)="$\\mu_{2}$ Pt $>$ 3 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Pt of Muon 2";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2pt_",htitle,80,0,20,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2pt_",htitle,80,0,20,hlabel,"Events"));
    }
    else if(i==Trkpt){
      title.at(i)="Track Pt $>$ 2 GeV";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Pt of Track";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Trkpt_",htitle,80,0,20,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Trkpt_",htitle,80,0,20,hlabel,"Events"));
    }


  }

  // Track Candidate Information
  //Track_P=HConfig.GetTH1D(Name+"_Track_P","Momentum magnitude of track (2mu+trk track candidate)",36,-0.5,35.5,"p (track)","Events");
  //Track_E=HConfig.GetTH1D(Name+"_Track_E","Energy of track (2mu+trk track candidate)",36,-0.5,35.5,"E (track)","Events");
  //Track_Pt=HConfig.GetTH1D(Name+"_Track_Pt","Transverse momentum of track (2mu+trk track candidate)",26,-0.5,25.5,"p_{T} (track)","Events");
  //Track_Eta=HConfig.GetTH1D(Name+"_Track_Eta","Psuedorapidity of track (2mu+trk track candidate)",30,-2.5,2.5,"#eta","Events");
  //Track_Phi=HConfig.GetTH1D(Name+"_Track_Phi","Azimuthal angle of track (2mu+trk track candidate)",30,-3.5,3.5,"#phi","Events");
  Track_normalizedChi2=HConfig.GetTH1D(Name+"_Track_normalizedChi2","Normalized chi square",20,-0.5,4.5,"#chi^{2} (track fit)","Events");
  //Track_numberOfValidHits=HConfig.GetTH1D(Name+"_Track_numberOfValidHits","number of valid hits in te tracker",36,-0.5,35.5,"n valid track hits","Events");
  //Track_charge=HConfig.GetTH1D(Name+"_Track_charge","Chargeof the track",3,-1.5,1.5,"Track charge","Events");
  
  //Muon1_Pt=HConfig.GetTH1D(Name+"_Muon1_Pt","Transverse Pt (muon 1)",25,0,30,"#mu_{1} p_{T} (GeV)","Events");
  //Muon1_Eta=HConfig.GetTH1D(Name+"_Muon1_Eta","Psuedorapidity (muon 1)",25,-2.5,2.5,"#mu_{1} #eta","Events");
  //Muon1_Phi=HConfig.GetTH1D(Name+"_Muon1_Phi","Azimuthal angle of (muons 1)",25,-3.4,3.4,"#mu_{1} #phi","Events"); 
  //Muon1_E=HConfig.GetTH1D(Name+"_Muon1_E","Energy of all (muon 1)",20,0,40,"#mu_{1} E (GeV)","Events");
  //Muon1_P=HConfig.GetTH1D(Name+"_Muon1_P","Magnitude of momentum of (muon 1)",20,0,40,"#mu_{1} p (GeV)","Events");  
  //
  //Muon2_Pt=HConfig.GetTH1D(Name+"_Muon2_Pt","Transverse Pt (muon 2)",25,0,30,"#mu_{2} p_{T} (GeV)","Events");
  //Muon2_Eta=HConfig.GetTH1D(Name+"_Muon2_Eta","Psuedorapidity (muon 2)",25,-2.5,2.5,"#mu_{2} #eta","Events");
  //Muon2_Phi=HConfig.GetTH1D(Name+"_Muon2_Phi","Azimuthal angle of (muons 1)",25,-3.4,3.4,"#mu_{2} #phi","Events"); 
  //Muon2_E=HConfig.GetTH1D(Name+"_Muon2_E","Energy of all (muon 2)",20,0,40,"#mu_{2} E (GeV)","Events");
  //Muon2_P=HConfig.GetTH1D(Name+"_Muon2_P","Magnitude of momentum of (muon 2)",20,0,40,"#mu_{2} p (GeV)","Events");  
  
  Muon1_TriggerMatchdR=HConfig.GetTH1D(Name+"_Muon1_TriggerMatchdR","Trigger Matching mu1",50,0,0.05,"#Delta R Trigger Match #mu_{1}","Events");
  Muon2_TriggerMatchdR=HConfig.GetTH1D(Name+"_Muon2_TriggerMatchdR","Trigger Matching mu2",50,0,0.05,"#Delta R Trigger Match #mu_{2}","Events");
  Track_TriggerMatchdR=HConfig.GetTH1D(Name+"_Track_TriggerMatchdR","Trigger Matching track",50,0,0.05,"#Delta R Trigger Match track","Events");
  
  //Muon1_isGlobal=HConfig.GetTH1D(Name+"_Muon1_isGlobal","Global muon status ",2,-.5,1.5,"#mu_{1} isGlb","Events");
  //Muon2_isGlobal=HConfig.GetTH1D(Name+"_Muon2_isGlobal","Global muon status",2,-0.5,0.5,"#mu_{2} isGlb","Events");
  //Muon1_isStandAlone=HConfig.GetTH1D(Name+"_Muon1_isStandAlone","",2,-0.5,1.5,"#mu_{1} isStandAlone","Events");
  //Muon2_isStandAlone=HConfig.GetTH1D(Name+"_Muon2_isStandAlone","",2,-0.5,1.5,"#mu_{2} isStandAlone","Events");
  //Muon1_isTracker=HConfig.GetTH1D(Name+"_Muon1_isTracker","",2,-0.5,1.5,"#mu_{1} isTracker","Events");
  //Muon2_isTracker=HConfig.GetTH1D(Name+"_Muon2_isTracker","",2,-0.5,1.5,"#mu_{2} isTracker","Events");
  //Muon1_isCalo=HConfig.GetTH1D(Name+"_Muon1_isCaloMuon","",2,-0.5,1.5,"#mu_{1} isCalo","Events");
  //Muon2_isCalo=HConfig.GetTH1D(Name+"_Muon2_isCaloMuon","",2,-0.5,1.5,"#mu_{2} isCalo","Events");
  
  
  //DimuondR=HConfig.GetTH1D(Name+"_DimuondR","dR between the muon pair",20,0,1,"dR","Events");
  //Muon1TrkdR=HConfig.GetTH1D(Name+"_Muon1TrkdR","dR between the highest p muon and the track",100,0,1,"dR","Events");
  //Muon2TrkdR=HConfig.GetTH1D(Name+"_Muon2TrkdR","dR between the lowest p muon and the track",100,0,1,"dR","Events");
  PhiMass=HConfig.GetTH1D(Name+"_PhiMass","#mu#mu invariant mass",50,0.2,1.5,"Mass of the #mu#mu pair","Events");
  PhiPlusTrackMass=HConfig.GetTH1D(Name+"_PhiPlusTrackMass","#mu#mu + track invariant mass",100,1.7,2.1,"Mass of the #mu#mu + track","Events");
  PhiMassVsDsMass=HConfig.GetTH2D(Name+"_PhiMassVsDsMass","#mu#mu invariant Mass vs. #mu#mu + track invariant mass",50,0.2,1.5,50,1.7,2.1,"M_{#mu#mu}, GeV","M_{#mu#mu + track}, GeV");
  
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms
  // Book here your analysis histrogramms, a good style is to follow selfexplanatory convention
  //NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",66,-0.5,65.5,"Number of Vertices","Events");
 
  DsMass=HConfig.GetTH1D(Name+"_DsMass","Ds invariant mass",100,1.7,2.1,"M_{Ds} (GeV)", "Events"); 
  DsPt_peak=HConfig.GetTH1D(Name+"_DsPt_peak","Transverse Pt (Ds) in Ds Peak",51,-0.5,50.5,"Ds p_{T} (GeV)", "Events");
  DsPt_sideband=HConfig.GetTH1D(Name+"_DsPt_sideband","Transverse Pt (Ds) in Ds Sideband",51,-0.5,50.5,"Ds p_{T} (GeV)", "Events");
  Ds_Pt=HConfig.GetTH1D(Name+"_DsPt","Transverse Pt (Ds)",51,-0.5,50.5,"Ds p_{T} (GeV)", "Events");
  DsP_peak=HConfig.GetTH1D(Name+"_DsP_peak","Transverse Pt (Ds) in Ds Peak",51,-0.5,50.5,"Ds p (GeV)", "Events");
  DsP_sideband=HConfig.GetTH1D(Name+"_DsP_sideband","Transverse Pt (Ds) in Ds Sideband",51,-0.5,50.5,"Ds p (GeV)", "Events");
  Ds_P=HConfig.GetTH1D(Name+"_Ds_P","Transverse Pt (Ds)",51,-0.5,50.5,"Ds p (GeV)", "Events");
  DsM_peak=HConfig.GetTH1D(Name+"_DsM_peak","Ds invariant mass in Ds Peak",100,1.7,2.1,"M_{Ds} (GeV)", "Events");
  DsM_sideband=HConfig.GetTH1D(Name+"_DsM_sideband","Ds invariant mass in Ds Sideband",100,1.7,2.1,"M_{Ds} (GeV)", "Events");
  Ds_M=HConfig.GetTH1D(Name+"_Ds_M","Ds invariant mass",100,1.7,2.1,"M_{Ds} (GeV)", "Events");
  DsL_peak=HConfig.GetTH1D(Name+"_DsL_peak","Length of Ds Decay in Ds Peak",25,0,1,"L_{Ds} (cm)", "Events");
  DsL_sideband=HConfig.GetTH1D(Name+"_DsL_sideband","Length of Ds Decay in Ds Sideband",25,0,1,"L_{Ds} (cm)", "Events");
  Ds_L=HConfig.GetTH1D(Name+"_Ds_L","Length of Ds Decay",25,0,1,"L_{Ds} (cm)", "Events");

  DsGenMatch=HConfig.GetTH1D(Name+"_DsGenMatch","dR between Gen Ds to Track",50,0,.1,"dR","Events");


  DecayLength_peak=HConfig.GetTH1D(Name+"_DecayLength_peak","Proper Decay Length of Ds in Ds Peak",20,0,.1,"Proper Decay Length (cm)","Events");
  DecayLength_sideband=HConfig.GetTH1D(Name+"_DecayLength_sideband","Proper Decay Length of Ds in sideband",20,0,.1,"Proper Decay Length (cm)","Events");
  DecayLength_prompt=HConfig.GetTH1D(Name+"_DecayLength_prompt","Proper Decay Length of Prompt Ds",20,0,.1,"Proper Decay Length (cm)","Events");
  DecayLength_non_prompt=HConfig.GetTH1D(Name+"_DecayLength_non_prompt","Proper Decay Length of Non-Prompt Ds",20,0,.1,"Proper Decay Length (cm)","Events");
  DecayLength=HConfig.GetTH1D(Name+"_DecayLength","Proper Decay Length of Ds",20,0,.1,"Proper Decay Length (cm)","Events");

  Muon1_Pt_peak=HConfig.GetTH1D(Name+"_Muon1_Pt_peak","Transverse Pt in Ds Peak (muon 1)",25,0,30,"#mu_{1} p_{T} (GeV)","Events");
  Muon1_Pt_sideband=HConfig.GetTH1D(Name+"_Muon1_Pt_sideband","Transverse Pt in Ds Sideband (muon 1)",25,0,30,"#mu_{1} p_{T} (GeV)","Events");
  Muon1_Eta_peak=HConfig.GetTH1D(Name+"_Muon1_Eta_peak","Psuedorapidity in Ds Peak (muon 1)",60,-2,2,"#mu_{1} #eta","Events");
  Muon1_Eta_sideband=HConfig.GetTH1D(Name+"_Muon1_Eta_sideband","Psuedorapidity in Ds Sideband (muon 1)",60,-2,2,"#mu_{1} #eta","Events");
  Muon1_Phi_peak=HConfig.GetTH1D(Name+"_Muon1_Phi_peak","Azimuthal angle in Ds Peak of (muons 1)",25,-3.4,3.4,"#mu_{1} #phi","Events");
  Muon1_Phi_sideband=HConfig.GetTH1D(Name+"_Muon1_Phi_sideband","Azimuthal angle in Ds Sideband of (muons 1)",25,-3.4,3.4,"#mu_{1} #phi","Events");
  control_Muon1_Pt=HConfig.GetTH1D(Name+"_control_Muon1_Pt","Transverse Pt (muon 1)",25,0,30,"#mu_{1} p_{T} (GeV)","Events");
  control_Muon1_Eta=HConfig.GetTH1D(Name+"_control_Muon1_Eta","Psuedorapidity (muon 1)",60,-2,2,"#mu_{1} #eta","Events");
  control_Muon1_Phi=HConfig.GetTH1D(Name+"_control_Muon1_Phi","Azimuthal angle of (muons 1)",25,-3.4,3.4,"#mu_{1} #phi","Events");

  Muon2_Pt_peak=HConfig.GetTH1D(Name+"_Muon2_Pt_peak","Transverse Pt in Ds Peak (muon 2)",25,0,30,"#mu_{2} p_{T} (GeV)","Events");
  Muon2_Pt_sideband=HConfig.GetTH1D(Name+"_Muon2_Pt_sideband","Transverse Pt in Ds Sideband (muon 2)",25,0,30,"#mu_{2} p_{T} (GeV)","Events");
  Muon2_Eta_peak=HConfig.GetTH1D(Name+"_Muon2_Eta_peak","Psuedorapidity in Ds Peak (muon 2)",60,-2,2,"#mu_{2} #eta","Events");
  Muon2_Eta_sideband=HConfig.GetTH1D(Name+"_Muon2_Eta_sideband","Psuedorapidity in Ds Sideband (muon 2)",60,-2,2,"#mu_{2} #eta","Events");
  Muon2_Phi_peak=HConfig.GetTH1D(Name+"_Muon2_Phi_peak","Azimuthal angle in Ds Peak of (muons 2)",25,-3.4,3.4,"#mu_{2} #phi","Events");
  Muon2_Phi_sideband=HConfig.GetTH1D(Name+"_Muon2_Phi_sideband","Azimuthal angle in Ds Sideband of (muons 2)",25,-3.4,3.4,"#mu_{2} #phi","Events");
  control_Muon2_Pt=HConfig.GetTH1D(Name+"_control_Muon2_Pt","Transverse Pt (muon 2)",25,0,30,"#mu_{2} p_{T} (GeV)","Events");
  control_Muon2_Eta=HConfig.GetTH1D(Name+"_control_Muon2_Eta","Psuedorapidity (muon 2)",60,-2,2,"#mu_{2} #eta","Events");
  control_Muon2_Phi=HConfig.GetTH1D(Name+"_control_Muon2_Phi","Azimuthal angle of (muons 2)",25,-3.4,3.4,"#mu_{2} #phi","Events");

  Track_Pt_peak=HConfig.GetTH1D(Name+"_Track_Pt_peak","Transverse Pt in Ds Peak (Track)",25,0,30,"Track p_{T} (GeV)","Events");
  Track_Pt_sideband=HConfig.GetTH1D(Name+"_Track_Pt_sideband","Transverse Pt in Ds Sideband (Track)",25,0,30,"Track p_{T} (GeV)","Events");
  Track_Eta_peak=HConfig.GetTH1D(Name+"_Track_Eta_peak","Psuedorapidity in Ds Peak (Track)",60,-2,2,"Track #eta","Events");
  Track_Eta_sideband=HConfig.GetTH1D(Name+"_Track_Eta_sideband","Psuedorapidity in Ds Sideband (Track)",60,-2,2,"Track #eta","Events");
  Track_Phi_peak=HConfig.GetTH1D(Name+"_Track_Phi_peak","Azimuthal angle in Ds Peak of (Track)",25,-3.4,3.4,"Track #phi","Events");
  Track_Phi_sideband=HConfig.GetTH1D(Name+"_Track_Phi_sideband","Azimuthal angle in Ds Sideband of (Track)",25,-3.4,3.4,"Track #phi","Events");
  control_Track_Pt=HConfig.GetTH1D(Name+"_control_Track_Pt","Transverse Pt (Track)",25,0,30,"Track p_{T} (GeV)","Events");
  control_Track_Eta=HConfig.GetTH1D(Name+"_control_Track_Eta","Psuedorapidity (Track)",60,-2,2,"Track #eta","Events");
  control_Track_Phi=HConfig.GetTH1D(Name+"_control_Track_Phi","Azimuthal angle of (Track)",25,-3.4,3.4,"Track #phi","Events");

  

  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}

void  DsToPhiPi::Store_ExtraDist(){ 
  
  //Extradist1d.push_back(&Track_P);
  //Extradist1d.push_back(&Track_E);
  //Extradist1d.push_back(&Track_Pt);
  //Extradist1d.push_back(&Track_Eta);
  //Extradist1d.push_back(&Track_Phi);
  Extradist1d.push_back(&Track_normalizedChi2);
  //Extradist1d.push_back(&Track_numberOfValidHits);
  //Extradist1d.push_back(&Track_charge);
  
  //Extradist1d.push_back(&Muon1_E);
  //Extradist1d.push_back(&Muon1_Pt);
  //Extradist1d.push_back(&Muon1_Phi);
  //Extradist1d.push_back(&Muon1_Eta);
  //Extradist1d.push_back(&Muon2_E);
  //Extradist1d.push_back(&Muon2_Pt);
  //Extradist1d.push_back(&Muon2_Phi);
  //Extradist1d.push_back(&Muon2_Eta);
  //Extradist1d.push_back(&DimuondR);
  //Extradist1d.push_back(&Muon1TrkdR);
  //Extradist1d.push_back(&Muon2TrkdR);
  Extradist1d.push_back(&PhiMass);
  Extradist1d.push_back(&PhiPlusTrackMass);
  Extradist2d.push_back(&PhiMassVsDsMass);
  //Extradist1d.push_back(&Muon1_isGlobal);
  //Extradist1d.push_back(&Muon2_isGlobal);
  //Extradist1d.push_back(&Muon1_isStandAlone);
  //Extradist1d.push_back(&Muon2_isStandAlone);
  //Extradist1d.push_back(&Muon1_isTracker);
  //Extradist1d.push_back(&Muon2_isTracker);
  //Extradist1d.push_back(&Muon1_isCalo);
  //Extradist1d.push_back(&Muon2_isCalo);
  Extradist1d.push_back(&Track_TriggerMatchdR);
  Extradist1d.push_back(&Muon1_TriggerMatchdR);
  Extradist1d.push_back(&Muon2_TriggerMatchdR);
  //Extradist1d.push_back(&NVtx);
  Extradist1d.push_back(&DsMass);
  Extradist1d.push_back(&Ds_Pt);
  Extradist1d.push_back(&Ds_P);
  Extradist1d.push_back(&Ds_M);
  Extradist1d.push_back(&Ds_L);
  Extradist1d.push_back(&DsGenMatch);

  Extradist1d.push_back(&DecayLength_prompt);
  Extradist1d.push_back(&DecayLength_non_prompt);
  Extradist1d.push_back(&DecayLength);
  Extradist1d.push_back(&DecayLength_peak);
  Extradist1d.push_back(&DecayLength_sideband);

  Extradist1d.push_back(&control_Muon1_Pt);
  Extradist1d.push_back(&control_Muon1_Eta);
  Extradist1d.push_back(&control_Muon1_Phi);	 
  Extradist1d.push_back(&control_Muon2_Pt);
  Extradist1d.push_back(&control_Muon2_Eta);
  Extradist1d.push_back(&control_Muon2_Phi);
  Extradist1d.push_back(&control_Track_Pt);
  Extradist1d.push_back(&control_Track_Eta);
  Extradist1d.push_back(&control_Track_Phi);

  Extradist1d.push_back(&Muon1_Pt_peak);
  Extradist1d.push_back(&Muon1_Pt_sideband);
  Extradist1d.push_back(&Muon1_Eta_peak);
  Extradist1d.push_back(&Muon1_Eta_sideband);
  Extradist1d.push_back(&Muon1_Phi_peak);
  Extradist1d.push_back(&Muon1_Phi_sideband);

  Extradist1d.push_back(&Muon2_Pt_peak);
  Extradist1d.push_back(&Muon2_Pt_sideband);
  Extradist1d.push_back(&Muon2_Eta_peak);
  Extradist1d.push_back(&Muon2_Eta_sideband);
  Extradist1d.push_back(&Muon2_Phi_peak);
  Extradist1d.push_back(&Muon2_Phi_sideband);

  Extradist1d.push_back(&Track_Pt_peak);
  Extradist1d.push_back(&Track_Pt_sideband);
  Extradist1d.push_back(&Track_Eta_peak);
  Extradist1d.push_back(&Track_Eta_sideband);
  Extradist1d.push_back(&Track_Phi_peak);
  Extradist1d.push_back(&Track_Phi_sideband);
}


void  DsToPhiPi::doEvent(){ 
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection

  value.at(HLTOk) = 0;
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);
    if((HLT.Contains("DoubleMu3_Trk_Tau3mu") || HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu") ) && Ntp->HLTDecision(iTrigger) == 1){
      value.at(HLTOk)=Ntp->HLTDecision(iTrigger);
    }
  }

  value.at(L1TOk) = 0;
  bool DoubleMuFired(0);
  for(unsigned int il1=0; il1 < Ntp->NL1Seeds(); il1++){
    TString L1TriggerName = Ntp->L1Name(il1);
    
    if(id==1 && Ntp->WhichEra(2017).Contains("RunB")){
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
    }

    if(id!=1){
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
    }

    if(id==1 && (Ntp->WhichEra(2017).Contains("RunC") || Ntp->WhichEra(2017).Contains("RunD") || Ntp->WhichEra(2017).Contains("RunF") || Ntp->WhichEra(2017).Contains("RunE"))){
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
    }
  }
  if(DoubleMuFired) value.at(L1TOk)=1;
   
  int mu1=-1, mu2=-1, track=-1;
  int tmp_idx = -1;
  double tmp_chisq = 999.0;
  double check_PhiMass = 999.0; 

  value.at(is2MuTrk) = 0; 
  value.at(Mass2Mu) = 0;
  value.at(MuCharge) = 0;
  value.at(Mu1dR) = 0;
  value.at(Mu2dR) = 0;
  value.at(TrkdR) = 0;
  if(Ntp->NTwoMuonsTrack()!=0 && Ntp->NThreeMuons() == 0){
    value.at(is2MuTrk) = 1;
  }

  if (value.at(is2MuTrk)==1){
    for(unsigned int i2M=0; i2M < Ntp->NTwoMuonsTrack(); i2M++){
      int tmp_mu1 = Ntp->TwoMuonsTrackMuonIndices(i2M).at(0);
      int tmp_mu2 = Ntp->TwoMuonsTrackMuonIndices(i2M).at(1);
      int tmp_track = Ntp->TwoMuonsTrackTrackIndex(i2M).at(0);
      double tmp_PhiMass = (Ntp->Muon_P4(tmp_mu1)+Ntp->Muon_P4(tmp_mu2)).M();

      if (abs(tmp_PhiMass-1.01)<=check_PhiMass || (tmp_PhiMass > .95 && tmp_PhiMass < 1.1)) {
      if (tmp_chisq>Ntp->TwoMuonsTrack_SV_Chi2(i2M+Ntp->NThreeMuons())){
	tmp_chisq = Ntp->TwoMuonsTrack_SV_Chi2(i2M+Ntp->NThreeMuons());
        check_PhiMass = abs(tmp_PhiMass-1.01);
	mu1 = tmp_mu1;
	mu2 = tmp_mu2;
	track = tmp_track;
        tmp_idx = i2M;
      }
      }
    }

    value.at(GlobalMu) = Ntp->Muon_isGlobalMuon(mu1)==1 && Ntp->Muon_isGlobalMuon(mu2)==1;
    value.at(Mass2Mu) = (Ntp->Muon_P4(mu1) + Ntp->Muon_P4(mu2)).M();
    value.at(MuCharge) = Ntp->Muon_charge(mu1)!=Ntp->Muon_charge(mu2);
    value.at(Mu1dR) = (Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(0);
    value.at(Mu2dR) = (Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(1);
    value.at(TrkdR) = (Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(2);
    value.at(Mu1pt) = Ntp->Muon_P4(mu1).Pt();
    value.at(Mu2pt) = Ntp->Muon_P4(mu2).Pt();
    value.at(Trkpt) = Ntp->Track_P4(track).Pt();

  }

  pass.at(is2MuTrk) = (value.at(is2MuTrk)==cut.at(is2MuTrk));
  pass.at(L1TOk)= (value.at(L1TOk)/*==cut.at(L1TOk)*/);
  pass.at(HLTOk)= (value.at(HLTOk)/*==cut.at(HLTOk)*/);
  pass.at(GlobalMu) = value.at(GlobalMu)==cut.at(GlobalMu);
  pass.at(Mass2Mu) = value.at(Mass2Mu) >= 1 && value.at(Mass2Mu) <= 1.04;
  pass.at(MuCharge) = value.at(MuCharge)==cut.at(MuCharge);
  pass.at(Mu1dR) = value.at(Mu1dR) < .03;
  pass.at(Mu2dR) = value.at(Mu2dR) < .03;
  pass.at(TrkdR) = value.at(TrkdR) < .03;
  pass.at(Mu1pt) = value.at(Mu1pt) > 3;
  pass.at(Mu2pt) = value.at(Mu2pt) > 3;
  pass.at(Trkpt) = value.at(Trkpt) > 2;

  double wobs=1;
  double w;  
  double w_peak;     

  if(!Ntp->isData()){w = 1; w_peak = .82;}//Ntp->PUReweight(); } //  No weights to data
  else{w=1; w_peak=1;}
  bool status=AnalysisCuts(t,w,wobs);
  if(status){
    //NVtx.at(t).Fill(Ntp->NVtx(),w);

    //Track_Pt.at(t).Fill(Ntp->Track_P4(track).Pt(),w);
    //Track_Eta.at(t).Fill(Ntp->Track_P4(track).Eta(),w);
    //Track_Phi.at(t).Fill(Ntp->Track_P4(track).Phi(),w);
    //Track_P.at(t).Fill(Ntp->Track_P4(track).P(),w);

    Track_normalizedChi2.at(t).Fill(Ntp->Track_normalizedChi2(track),w);
    //Track_numberOfValidHits.at(t).Fill(Ntp->Track_numberOfValidHits(track),w);
    //Track_charge.at(t).Fill(Ntp->Track_charge(track),w);
    //Muon1_isGlobal.at(t).Fill(Ntp->Muon_isGlobalMuon(mu1),w);
    //Muon2_isGlobal.at(t).Fill(Ntp->Muon_isGlobalMuon(mu2),w);
    //Muon1_isStandAlone.at(t).Fill(Ntp->Muon_isStandAloneMuon(mu1),w);
    //Muon2_isStandAlone.at(t).Fill(Ntp->Muon_isStandAloneMuon(mu2),w);
    //Muon1_isTracker.at(t).Fill(Ntp->Muon_isTrackerMuon(mu1),w);
    //Muon2_isTracker.at(t).Fill(Ntp->Muon_isTrackerMuon(mu2),w);
 

    Muon1_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(0),w);
    Muon2_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(1),w);
    Track_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(2),w);


    //Muon1_Pt.at(t).Fill(Ntp->Muon_P4(mu1).Pt(),w);
    //Muon1_Eta.at(t).Fill(Ntp->Muon_P4(mu1).Eta(),w);
    //Muon1_Phi.at(t).Fill(Ntp->Muon_P4(mu1).Phi(),w);
    //Muon1_E.at(t).Fill(Ntp->Muon_P4(mu1).E(),w);


    //Muon2_Pt.at(t).Fill(Ntp->Muon_P4(mu2).Pt(),w);
    //Muon2_Eta.at(t).Fill(Ntp->Muon_P4(mu2).Eta(),w);
    //Muon2_Phi.at(t).Fill(Ntp->Muon_P4(mu2).Phi(),w);
    //Muon2_E.at(t).Fill(Ntp->Muon_P4(mu2).E(),w);

    //
    //DimuondR.at(t).Fill(Ntp->Muon_P4(mu1).DeltaR(Ntp->Muon_P4(mu2)),w);
    //Muon1TrkdR.at(t).Fill(Ntp->Muon_P4(mu1).DeltaR(Ntp->Track_P4(track)),w);
    //Muon2TrkdR.at(t).Fill(Ntp->Muon_P4(mu2).DeltaR(Ntp->Track_P4(track)),w);
    
    PhiMass.at(t).Fill((Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))).M(), w);
    PhiPlusTrackMass.at(t).Fill((Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))+ 
			   Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0))).M(), w);
    
    double phimass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))).M();
    double dsmass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))+
		     Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0))).M();
    double dsPt = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))+
                     Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0))).Pt();
    double dsP = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))+
                     Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0))).P();
    PhiMassVsDsMass.at(t).Fill(phimass, dsmass);

    DsMass.at(t).Fill(dsmass,w);

    if(id==30){

      DsGenMatch.at(t).Fill(Ntp->DsGenMatch(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0) + Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(1) + Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0)));

    }

    int vertex_idx = tmp_idx + Ntp->NThreeMuons();
    TVector3 SVPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(vertex_idx),Ntp->Vertex_MatchedPrimaryVertex(vertex_idx));
    double DecayL = SVPV.Mag()*dsmass/dsP;
    if(dsmass > 1.93 && dsmass < 2.01){
      DecayLength_peak.at(t).Fill(DecayL,w);
      Muon1_Pt_peak.at(t).Fill(Ntp->Muon_P4(mu1).Pt(),w);
      Muon1_Eta_peak.at(t).Fill(Ntp->Muon_P4(mu1).Eta(),w);
      Muon1_Phi_peak.at(t).Fill(Ntp->Muon_P4(mu1).Phi(),w);
      Muon2_Pt_peak.at(t).Fill(Ntp->Muon_P4(mu2).Pt(),w);
      Muon2_Eta_peak.at(t).Fill(Ntp->Muon_P4(mu2).Eta(),w);
      Muon2_Phi_peak.at(t).Fill(Ntp->Muon_P4(mu2).Phi(),w);
      Track_Pt_peak.at(t).Fill(Ntp->Track_P4(track).Pt(),w);
      Track_Eta_peak.at(t).Fill(Ntp->Track_P4(track).Eta(),w);
      Track_Phi_peak.at(t).Fill(Ntp->Track_P4(track).Phi(),w);
      
      DsPt_peak.at(t).Fill(dsPt,w);
      DsP_peak.at(t).Fill(dsP,w);
      DsM_peak.at(t).Fill(dsmass,w);
      DsL_peak.at(t).Fill(SVPV.Mag(),w);
    }
    if(dsmass > 1.70 && dsmass < 1.80){
      DecayLength_sideband.at(t).Fill(DecayL,w);
      Muon1_Pt_sideband.at(t).Fill(Ntp->Muon_P4(mu1).Pt(),w);
      Muon1_Eta_sideband.at(t).Fill(Ntp->Muon_P4(mu1).Eta(),w);
      Muon1_Phi_sideband.at(t).Fill(Ntp->Muon_P4(mu1).Phi(),w);
      Muon2_Pt_sideband.at(t).Fill(Ntp->Muon_P4(mu2).Pt(),w);
      Muon2_Eta_sideband.at(t).Fill(Ntp->Muon_P4(mu2).Eta(),w);
      Muon2_Phi_sideband.at(t).Fill(Ntp->Muon_P4(mu2).Phi(),w);
      Track_Pt_sideband.at(t).Fill(Ntp->Track_P4(track).Pt(),w);
      Track_Eta_sideband.at(t).Fill(Ntp->Track_P4(track).Eta(),w);
      Track_Phi_sideband.at(t).Fill(Ntp->Track_P4(track).Phi(),w);

      DsPt_sideband.at(t).Fill(dsPt,w);
      DsP_sideband.at(t).Fill(dsP,w);
      DsM_sideband.at(t).Fill(dsmass,w);
      DsL_sideband.at(t).Fill(SVPV.Mag(),w);
    }

    bool isPrompt(true);
    if(id!=1){

      control_Muon1_Pt.at(t).Fill(Ntp->Muon_P4(mu1).Pt(),w*w_peak);
      control_Muon1_Eta.at(t).Fill(Ntp->Muon_P4(mu1).Eta(),w*w_peak);
      control_Muon1_Phi.at(t).Fill(Ntp->Muon_P4(mu1).Phi(),w*w_peak);
      control_Muon2_Pt.at(t).Fill(Ntp->Muon_P4(mu2).Pt(),w*w_peak);
      control_Muon2_Eta.at(t).Fill(Ntp->Muon_P4(mu2).Eta(),w*w_peak);
      control_Muon2_Phi.at(t).Fill(Ntp->Muon_P4(mu2).Phi(),w*w_peak);
      control_Track_Pt.at(t).Fill(Ntp->Track_P4(track).Pt(),w*w_peak);
      control_Track_Eta.at(t).Fill(Ntp->Track_P4(track).Eta(),w*w_peak);
      control_Track_Phi.at(t).Fill(Ntp->Track_P4(track).Phi(),w*w_peak);

      Ds_Pt.at(t).Fill(dsPt,w*w_peak);
      Ds_P.at(t).Fill(dsP,w*w_peak);
      Ds_M.at(t).Fill(dsmass,w*w_peak);
      Ds_L.at(t).Fill(SVPV.Mag(),w*w_peak);

      for (unsigned int isigp=0; isigp<Ntp->NMCSignalParticles(); isigp++){
        for (int is=0; is<Ntp->NMCSignalParticleSources(isigp); is++){
          if (abs(Ntp->MCSignalParticle_Sourcepdgid(isigp,is))>500){
            isPrompt=false;
          }
        }
      }
      
      if (isPrompt) {DecayLength_prompt.at(t).Fill(DecayL,w*w_peak);}
      else if (!isPrompt) {DecayLength_non_prompt.at(t).Fill(DecayL,w*w_peak);}
      DecayLength.at(t).Fill(DecayL,w*w_peak);
    }


  }
}

void  DsToPhiPi::Finish(){

  int id(Ntp->GetMCID());
  if (id==1) {

    std::vector<double> scaleRun;
    //
    //if(Ntp->WhichEra(2017).Contains("RunB") ){
    //  scaleRun.push_back(2265.07);scaleRun.push_back(4586.64);
    //}
    //if(Ntp->WhichEra(2017).Contains("RunC") ){
    //  scaleRun.push_back(15965.4);scaleRun.push_back(24220.4);
    //}
    //if(Ntp->WhichEra(2017).Contains("RunD") ){
      scaleRun.push_back(5952.73);scaleRun.push_back(11303.6);
    //}
    //if(Ntp->WhichEra(2017).Contains("RunE") ){
    //  scaleRun.push_back(10661.2);scaleRun.push_back(19461.1);
    //}
    //if(Ntp->WhichEra(2017).Contains("RunF") ){
    //  scaleRun.push_back(10093.0);scaleRun.push_back(19046.7);
    //}
    DecayLength_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);//DecayLength_sideband.at(0).Integral());
    Muon1_Pt_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    Muon1_Eta_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    Muon1_Phi_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    Muon2_Pt_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    Muon2_Eta_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    Muon2_Phi_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    Track_Pt_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    Track_Eta_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    Track_Phi_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    DsPt_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    DsP_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    DsM_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    DsL_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);

    DecayLength_prompt.at(0).Add(&DecayLength_peak.at(0));
    DecayLength_prompt.at(0).Add(&DecayLength_sideband.at(0),-1);
    DecayLength_non_prompt.at(0).Add(&DecayLength_peak.at(0));
    DecayLength_non_prompt.at(0).Add(&DecayLength_sideband.at(0),-1);
    DecayLength.at(0).Add(&DecayLength_peak.at(0));
    DecayLength.at(0).Add(&DecayLength_sideband.at(0),-1); 

    control_Muon1_Pt.at(0).Add(&Muon1_Pt_peak.at(0));
    control_Muon1_Pt.at(0).Add(&Muon1_Pt_sideband.at(0),-1);
    control_Muon1_Eta.at(0).Add(&Muon1_Eta_peak.at(0));
    control_Muon1_Eta.at(0).Add(&Muon1_Eta_sideband.at(0),-1);   
    control_Muon1_Phi.at(0).Add(&Muon1_Phi_peak.at(0));
    control_Muon1_Phi.at(0).Add(&Muon1_Phi_sideband.at(0),-1);
    control_Muon2_Pt.at(0).Add(&Muon2_Pt_peak.at(0));
    control_Muon2_Pt.at(0).Add(&Muon2_Pt_sideband.at(0),-1);
    control_Muon2_Eta.at(0).Add(&Muon2_Eta_peak.at(0));
    control_Muon2_Eta.at(0).Add(&Muon2_Eta_sideband.at(0),-1);
    control_Muon2_Phi.at(0).Add(&Muon2_Phi_peak.at(0));
    control_Muon2_Phi.at(0).Add(&Muon2_Phi_sideband.at(0),-1);
    control_Track_Pt.at(0).Add(&Track_Pt_peak.at(0));
    control_Track_Pt.at(0).Add(&Track_Pt_sideband.at(0),-1);
    control_Track_Eta.at(0).Add(&Track_Eta_peak.at(0));
    control_Track_Eta.at(0).Add(&Track_Eta_sideband.at(0),-1);
    control_Track_Phi.at(0).Add(&Track_Phi_peak.at(0));
    control_Track_Phi.at(0).Add(&Track_Phi_sideband.at(0),-1);

    Ds_Pt.at(0).Add(&DsPt_peak.at(0));
    Ds_Pt.at(0).Add(&DsPt_sideband.at(0),-1);
    Ds_P.at(0).Add(&DsP_peak.at(0));
    Ds_P.at(0).Add(&DsP_sideband.at(0),-1);
    Ds_M.at(0).Add(&DsM_peak.at(0));
    Ds_M.at(0).Add(&DsM_sideband.at(0),-1);
    Ds_L.at(0).Add(&DsL_peak.at(0));
    Ds_L.at(0).Add(&DsL_sideband.at(0),-1);

  }

  //if(mode == RECONSTRUCT){
  //  for(unsigned int i=1; i<  Nminus0.at(0).size(); i++){
  //    double scale(1.);
  //    if(Nminus0.at(0).at(i).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(i).Integral()/1;
  //    ScaleAllHistOfType(i,scale);
  //  }
  //}
  Selection::Finish();

}



