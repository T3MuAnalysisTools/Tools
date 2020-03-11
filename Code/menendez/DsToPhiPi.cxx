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
  DsP_peak=HConfig.GetTH1D(Name+"_DsP_peak","Transverse P (Ds) in Ds Peak",51,-0.5,50.5,"Ds p (GeV)", "Events");
  DsP_sideband=HConfig.GetTH1D(Name+"_DsP_sideband","Transverse P (Ds) in Ds Sideband",51,-0.5,50.5,"Ds p (GeV)", "Events");
  Ds_P=HConfig.GetTH1D(Name+"_Ds_P","Transverse P (Ds)",51,-0.5,50.5,"Ds p (GeV)", "Events");
  DsM_peak=HConfig.GetTH1D(Name+"_DsM_peak","Ds invariant mass in Ds Peak",100,1.7,2.1,"M_{Ds} (GeV)", "Events");
  DsM_sideband=HConfig.GetTH1D(Name+"_DsM_sideband","Ds invariant mass in Ds Sideband",100,1.7,2.1,"M_{Ds} (GeV)", "Events");
  Ds_M=HConfig.GetTH1D(Name+"_Ds_M","Ds invariant mass",100,1.7,2.1,"M_{Ds} (GeV)", "Events");
  DsL_peak=HConfig.GetTH1D(Name+"_DsL_peak","Length of Ds Decay in Ds Peak",25,0,1,"L_{Ds} (cm)", "Events");
  DsL_sideband=HConfig.GetTH1D(Name+"_DsL_sideband","Length of Ds Decay in Ds Sideband",25,0,1,"L_{Ds} (cm)", "Events");
  Ds_L=HConfig.GetTH1D(Name+"_Ds_L","Length of Ds Decay",25,0,1,"L_{Ds} (cm)", "Events");
  DsEta_peak=HConfig.GetTH1D(Name+"_DsEta_peak","Psuedorapidity (Ds)",30,-2,2,"Ds #eta","Events");
  DsEta_sideband=HConfig.GetTH1D(Name+"_DsEta_sideband","Psuedorapidity (Ds)",30,-2,2,"Ds #eta","Events");
  Ds_eta=HConfig.GetTH1D(Name+"_Ds_eta","Psuedorapidity (Ds)",30,-2,2,"Ds #eta","Events");

  DsGenMatch=HConfig.GetTH1D(Name+"_DsGenMatch","dR between Gen Ds to Track",50,0,.1,"dR","Events");


  DecayLength_peak=HConfig.GetTH1D(Name+"_DecayLength_peak","Proper Decay Length of Ds in Ds Peak",20,0,.1,"Proper Decay Length (cm)","Events");
  DecayLength_sideband=HConfig.GetTH1D(Name+"_DecayLength_sideband","Proper Decay Length of Ds in sideband",20,0,.1,"Proper Decay Length (cm)","Events");
  DecayLength_prompt=HConfig.GetTH1D(Name+"_DecayLength_prompt","Proper Decay Length of Prompt Ds",20,0,.1,"Proper Decay Length (cm)","Events");
  DecayLength_non_prompt=HConfig.GetTH1D(Name+"_DecayLength_non_prompt","Proper Decay Length of Non-Prompt Ds",20,0,.1,"Proper Decay Length (cm)","Events");
  DecayLength=HConfig.GetTH1D(Name+"_DecayLength","Proper Decay Length of Ds",20,0,.1,"Proper Decay Length (cm)","Events");

  Muon1_Pt_peak=HConfig.GetTH1D(Name+"_Muon1_Pt_peak","Transverse Pt in Ds Peak (muon 1)",25,0,20,"#mu_{1} p_{T} (GeV)","Events");
  Muon1_Pt_sideband=HConfig.GetTH1D(Name+"_Muon1_Pt_sideband","Transverse Pt in Ds Sideband (muon 1)",25,0,20,"#mu_{1} p_{T} (GeV)","Events");
  Muon1_Eta_peak=HConfig.GetTH1D(Name+"_Muon1_Eta_peak","Psuedorapidity in Ds Peak (muon 1)",30,-2,2,"#mu_{1} #eta","Events");
  Muon1_Eta_sideband=HConfig.GetTH1D(Name+"_Muon1_Eta_sideband","Psuedorapidity in Ds Sideband (muon 1)",30,-2,2,"#mu_{1} #eta","Events");
  Muon1_Phi_peak=HConfig.GetTH1D(Name+"_Muon1_Phi_peak","Azimuthal angle in Ds Peak of (muons 1)",25,-3.4,3.4,"#mu_{1} #phi","Events");
  Muon1_Phi_sideband=HConfig.GetTH1D(Name+"_Muon1_Phi_sideband","Azimuthal angle in Ds Sideband of (muons 1)",25,-3.4,3.4,"#mu_{1} #phi","Events");
  control_Muon1_Pt=HConfig.GetTH1D(Name+"_control_Muon1_Pt","Transverse Pt (muon 1)",25,0,20,"#mu_{1} p_{T} (GeV)","Events");
  control_Muon1_Eta=HConfig.GetTH1D(Name+"_control_Muon1_Eta","Psuedorapidity (muon 1)",30,-2,2,"#mu_{1} #eta","Events");
  control_Muon1_Phi=HConfig.GetTH1D(Name+"_control_Muon1_Phi","Azimuthal angle of (muons 1)",25,-3.4,3.4,"#mu_{1} #phi","Events");

  Muon2_Pt_peak=HConfig.GetTH1D(Name+"_Muon2_Pt_peak","Transverse Pt in Ds Peak (muon 2)",25,0,20,"#mu_{2} p_{T} (GeV)","Events");
  Muon2_Pt_sideband=HConfig.GetTH1D(Name+"_Muon2_Pt_sideband","Transverse Pt in Ds Sideband (muon 2)",25,0,20,"#mu_{2} p_{T} (GeV)","Events");
  Muon2_Eta_peak=HConfig.GetTH1D(Name+"_Muon2_Eta_peak","Psuedorapidity in Ds Peak (muon 2)",30,-2,2,"#mu_{2} #eta","Events");
  Muon2_Eta_sideband=HConfig.GetTH1D(Name+"_Muon2_Eta_sideband","Psuedorapidity in Ds Sideband (muon 2)",30,-2,2,"#mu_{2} #eta","Events");
  Muon2_Phi_peak=HConfig.GetTH1D(Name+"_Muon2_Phi_peak","Azimuthal angle in Ds Peak of (muons 2)",25,-3.4,3.4,"#mu_{2} #phi","Events");
  Muon2_Phi_sideband=HConfig.GetTH1D(Name+"_Muon2_Phi_sideband","Azimuthal angle in Ds Sideband of (muons 2)",25,-3.4,3.4,"#mu_{2} #phi","Events");
  control_Muon2_Pt=HConfig.GetTH1D(Name+"_control_Muon2_Pt","Transverse Pt (muon 2)",25,0,20,"#mu_{2} p_{T} (GeV)","Events");
  control_Muon2_Eta=HConfig.GetTH1D(Name+"_control_Muon2_Eta","Psuedorapidity (muon 2)",30,-2,2,"#mu_{2} #eta","Events");
  control_Muon2_Phi=HConfig.GetTH1D(Name+"_control_Muon2_Phi","Azimuthal angle of (muons 2)",25,-3.4,3.4,"#mu_{2} #phi","Events");

  Track_Pt_peak=HConfig.GetTH1D(Name+"_Track_Pt_peak","Transverse Pt in Ds Peak (Track)",25,0,20,"Track p_{T} (GeV)","Events");
  Track_Pt_sideband=HConfig.GetTH1D(Name+"_Track_Pt_sideband","Transverse Pt in Ds Sideband (Track)",25,0,20,"Track p_{T} (GeV)","Events");
  Track_Eta_peak=HConfig.GetTH1D(Name+"_Track_Eta_peak","Psuedorapidity in Ds Peak (Track)",30,-2,2,"Track #eta","Events");
  Track_Eta_sideband=HConfig.GetTH1D(Name+"_Track_Eta_sideband","Psuedorapidity in Ds Sideband (Track)",30,-2,2,"Track #eta","Events");
  Track_Phi_peak=HConfig.GetTH1D(Name+"_Track_Phi_peak","Azimuthal angle in Ds Peak of (Track)",25,-3.4,3.4,"Track #phi","Events");
  Track_Phi_sideband=HConfig.GetTH1D(Name+"_Track_Phi_sideband","Azimuthal angle in Ds Sideband of (Track)",25,-3.4,3.4,"Track #phi","Events");
  control_Track_Pt=HConfig.GetTH1D(Name+"_control_Track_Pt","Transverse Pt (Track)",25,0,20,"Track p_{T} (GeV)","Events");
  control_Track_Eta=HConfig.GetTH1D(Name+"_control_Track_Eta","Psuedorapidity (Track)",30,-2,2,"Track #eta","Events");
  control_Track_Phi=HConfig.GetTH1D(Name+"_control_Track_Phi","Azimuthal angle of (Track)",25,-3.4,3.4,"Track #phi","Events");

  VertexKFChi2_peak=HConfig.GetTH1D(Name+"VertexKFChi2_peak","KF Vertex Chi Squared",50,0,20,"KF vertex #chi^{2}","Events");
  VertexKFChi2_sideband=HConfig.GetTH1D(Name+"VertexKFChi2_sideband","KF Vertex Chi Squared",50,0,20,"KF vertex #chi^{2}","Events");
  VertexKFChi2=HConfig.GetTH1D(Name+"VertexKFChi2","KF Vertex Chi Squared",50,0,20,"KF vertex #chi^{2}","Events");

  SVPVDsDirAngle_peak=HConfig.GetTH1D(Name+"_SVPVDsDirAngle_peak","SVPVDsDirAngle",100,0,0.15,"Angle btw #vec{SV}-#vec{PV} and #vec{Ds}, rad","Events");
  SVPVDsDirAngle_sideband=HConfig.GetTH1D(Name+"_SVPVDsDirAngle_sideband","SVPVDsDirAngle",100,0,0.15,"Angle btw #vec{SV}-#vec{PV} and #vec{Ds}, rad","Events");
  SVPVDsDirAngle=HConfig.GetTH1D(Name+"_SVPVDsDirAngle","SVPVDsDirAngle",100,0,0.15,"Angle btw #vec{SV}-#vec{PV} and #vec{Ds}, rad","Events");

  NtracksClose_peak=HConfig.GetTH1D(Name+"_NtracksClose_peak","NtracksClose",8,-0.5,7.5,"Number of tracks close to SV","");
  NtracksClose_sideband=HConfig.GetTH1D(Name+"_NtracksClose_sideband","NtracksClose",8,-0.5,7.5,"Number of tracks close to SV","");
  NtracksClose=HConfig.GetTH1D(Name+"_NtracksClose","NtracksClose",8,-0.5,7.5,"Number of tracks close to SV","");

  NSV_peak=HConfig.GetTH1D(Name+"_NSV_peak","NSV",8,-0.5,7.5,"N vertices in the Ds cone","");
  NSV_sideband=HConfig.GetTH1D(Name+"_NSV_sideband","NSV",8,-0.5,7.5,"N vertices in the Ds cone","");
  NSV=HConfig.GetTH1D(Name+"_NSV","NSV",8,-0.5,7.5,"N vertices in the Ds cone","");

  MinMuon_chi2LocalPosition_peak=HConfig.GetTH1D(Name+"_MinMuon_chi2LocalPosition_peak","MinMuon_chi2LocalPosition",50,0,5,"Min Inner/Outer track #chi^{2}","");
  MinMuon_chi2LocalPosition_sideband=HConfig.GetTH1D(Name+"_MinMuon_chi2LocalPosition_sideband","MinMuon_chi2LocalPosition",50,0,5,"Min Inner/Outer track #chi^{2}","");
  MinMuon_chi2LocalPosition=HConfig.GetTH1D(Name+"_MinMuon_chi2LocalPosition","MinMuon_chi2LocalPosition",50,0,5,"Min Inner/Outer track #chi^{2}","");
 
  MindcaTrackSV_peak=HConfig.GetTH1D(Name+"_MindcaTrackSV_peak","MindcaTrackSV",50,0,0.1,"Min distance of track to SV","");
  MindcaTrackSV_sideband=HConfig.GetTH1D(Name+"_MindcaTrackSV_sideband","MindcaTrackSV",50,0,0.1,"Min distance of track to SV","");
  MindcaTrackSV=HConfig.GetTH1D(Name+"_MindcaTrackSV","MindcaTrackSV",50,0,0.1,"Min distance of track to SV",""); 

  MinDca_peak=HConfig.GetTH1D(Name+"_MinDca_peak","MinDca",50,0,0.02,"Min distance between muons and track","");
  MinDca_sideband=HConfig.GetTH1D(Name+"_MinDca_sideband","MinDca",50,0,0.02,"Min distance between muons and track","");
  MinDca=HConfig.GetTH1D(Name+"_MinDca","MinDca",50,0,0.02,"Min distance between muons and track","");

  MinD0SigSV_peak=HConfig.GetTH1D(Name+"_MinD0SigSV_peak","MinD0SigSV",30,0,1.5,"Min Transverse Impact significance w.r.t SV","");
  MinD0SigSV_sideband=HConfig.GetTH1D(Name+"_MinD0SigSV_sideband","MinD0SigSV",30,0,1.5,"Min Transverse Impact significance w.r.t SV","");
  MinD0SigSV=HConfig.GetTH1D(Name+"_MinD0SigSV","MinD0SigSV",30,0,1.5,"Min Transverse Impact significance w.r.t SV","");

  MinD0SigPV_peak=HConfig.GetTH1D(Name+"_MinD0SigPV_peak","MinD0SigPV",20,0,20,"Min Transverse Impact significance w.r.t PV","");
  MinD0SigPV_sideband=HConfig.GetTH1D(Name+"_MinD0SigPV_sideband","MinD0SigPV",20,0,20,"Min Transverse Impact significance w.r.t PV","");
  MinD0SigPV=HConfig.GetTH1D(Name+"_MinD0SigPV","MinD0SigPV",20,0,20,"Min Transverse Impact significance w.r.t PV","");

  MaxVertexPairQuality_peak=HConfig.GetTH1D(Name+"_MaxVertexPairQuality_peak","MaxVertexPairQuality",30,0,10,"max vertex pair quality","Events");
  MaxVertexPairQuality_sideband=HConfig.GetTH1D(Name+"_MaxVertexPairQuality_sideband","MaxVertexPairQuality",30,0,10,"max vertex pair quality","Events");
  MaxVertexPairQuality=HConfig.GetTH1D(Name+"_MaxVertexPairQuality","MaxVertexPairQuality",30,0,10,"max vertex pair quality","Events");

  MaxdeltaMuZ_peak = HConfig.GetTH1D(Name+"_MaxdeltaMuZ_peak","MaxdeltaMuZ",30,0,0.6,"Max #Delta z (#mu-#mu), cm","");
  MaxdeltaMuZ_sideband = HConfig.GetTH1D(Name+"_MaxdeltaMuZ_sideband","MaxdeltaMuZ",30,0,0.6,"Max #Delta z (#mu-#mu), cm","");
  MaxdeltaMuZ = HConfig.GetTH1D(Name+"_MaxdeltaMuZ","MaxdeltaMuZ",30,0,0.6,"Max #Delta z (#mu-#mu), cm","");

  MaxDca_peak=HConfig.GetTH1D(Name+"_MaxDca_peak","MaxDca",50,0,0.10,"Max distance between muons","");
  MaxDca_sideband=HConfig.GetTH1D(Name+"_MaxDca_sideband","MaxDca",50,0,0.10,"Max distance between muons","");
  MaxDca=HConfig.GetTH1D(Name+"_MaxDca","MaxDca",50,0,0.10,"Max distance between muons","");

  MaxD0SigSV_peak=HConfig.GetTH1D(Name+"_MaxD0SigSV_peak","MaxD0SigSV",20,0,5,"Max Transverse Impact significance w.r.t SV","");
  MaxD0SigSV_sideband=HConfig.GetTH1D(Name+"_MaxD0SigSV_sideband","MaxD0SigSV",20,0,5,"Max Transverse Impact significance w.r.t SV","");
  MaxD0SigSV=HConfig.GetTH1D(Name+"_MaxD0SigSV","MaxD0SigSV",20,0,5,"Max Transverse Impact significance w.r.t SV","");

  MaxD0SigPV_peak=HConfig.GetTH1D(Name+"_MaxD0SigPV_peak","MaxD0SigPV",20,0,20,"Max Transverse Impact significance w.r.t PV","");
  MaxD0SigPV_sideband=HConfig.GetTH1D(Name+"_MaxD0SigPV_sideband","MaxD0SigPV",20,0,20,"Max Transverse Impact significance w.r.t PV","");
  MaxD0SigPV=HConfig.GetTH1D(Name+"_MaxD0SigPV","MaxD0SigPV",20,0,20,"Max Transverse Impact significance w.r.t PV","");

  Iso1_peak=HConfig.GetTH1D(Name+"_Iso1_peak","Iso1",30,0,1.1,"I= p_{T}(Ds)/(p_{T}(Ds) + #sum p_{T}(tracks))","#Delta R < 1.0");
  Iso1_sideband=HConfig.GetTH1D(Name+"_Iso1_sideband","Iso1",30,0,1.1,"I= p_{T}(Ds)/(p_{T}(Ds) + #sum p_{T}(tracks))","#Delta R < 1.0");
  Iso1=HConfig.GetTH1D(Name+"_Iso1","Iso1",30,0,1.1,"I= p_{T}(Ds)/(p_{T}(Ds) + #sum p_{T}(tracks))","#Delta R < 1.0");

  Iso1Mu1_peak=HConfig.GetTH1D(Name+"_Iso1Mu1_peak","Iso1Mu1",30,0,1.1,"I= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 1.0");
  Iso1Mu1_sideband=HConfig.GetTH1D(Name+"_Iso1Mu1_sideband","Iso1Mu1",30,0,1.1,"I= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 1.0");
  Iso1Mu1=HConfig.GetTH1D(Name+"_Iso1Mu1","Iso1Mu1",30,0,1.1,"I= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 1.0");

  Iso8Mu1_peak=HConfig.GetTH1D(Name+"_Iso8Mu1_peak","Iso8Mu1",30,0,1.1,"I= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 0.8");
  Iso8Mu1_sideband=HConfig.GetTH1D(Name+"_Iso8Mu1_sideband","Iso8Mu1",30,0,1.1,"I= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 0.8");
  Iso8Mu1=HConfig.GetTH1D(Name+"_Iso8Mu1","Iso8Mu1",30,0,1.1,"I= p_{T}(#mu1)/(p_{T}(#mu1) + #sum p_{T}(tracks))","#Delta R < 0.8");

  FLSignificance_peak=HConfig.GetTH1D(Name+"_FLSignificance_peak","FLSignificance",60,0,60,"PV - SV distance  significance","Events");
  FLSignificance_sideband=HConfig.GetTH1D(Name+"_FLSignificance_sideband","FLSignificance",60,0,60,"PV - SV distance  significance","Events");
  FLSignificance=HConfig.GetTH1D(Name+"_FLSignificance","FLSignificance",60,0,60,"PV - SV distance  significance","Events");

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
  //Extradist1d.push_back(&Track_TriggerMatchdR);
  //Extradist1d.push_back(&Muon1_TriggerMatchdR);
  //Extradist1d.push_back(&Muon2_TriggerMatchdR);
  //Extradist1d.push_back(&NVtx);
  Extradist1d.push_back(&DsMass);
  Extradist1d.push_back(&Ds_Pt);
  Extradist1d.push_back(&Ds_P);
  Extradist1d.push_back(&Ds_M);
  Extradist1d.push_back(&Ds_L);
  Extradist1d.push_back(&Ds_eta);
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

  //Extradist1d.push_back(&Muon1_Pt_peak);
  //Extradist1d.push_back(&Muon1_Pt_sideband);
  //Extradist1d.push_back(&Muon1_Eta_peak);
  //Extradist1d.push_back(&Muon1_Eta_sideband);
  //Extradist1d.push_back(&Muon1_Phi_peak);
  //Extradist1d.push_back(&Muon1_Phi_sideband);

  //Extradist1d.push_back(&Muon2_Pt_peak);
  //Extradist1d.push_back(&Muon2_Pt_sideband);
  //Extradist1d.push_back(&Muon2_Eta_peak);
  //Extradist1d.push_back(&Muon2_Eta_sideband);
  //Extradist1d.push_back(&Muon2_Phi_peak);
  //Extradist1d.push_back(&Muon2_Phi_sideband);

  //Extradist1d.push_back(&Track_Pt_peak);
  //Extradist1d.push_back(&Track_Pt_sideband);
  //Extradist1d.push_back(&Track_Eta_peak);
  //Extradist1d.push_back(&Track_Eta_sideband);
  //Extradist1d.push_back(&Track_Phi_peak);
  //Extradist1d.push_back(&Track_Phi_sideband);

  Extradist1d.push_back(&VertexKFChi2);
  Extradist1d.push_back(&SVPVDsDirAngle);
  Extradist1d.push_back(&NtracksClose);
  Extradist1d.push_back(&NSV);
  Extradist1d.push_back(&MinMuon_chi2LocalPosition);
  Extradist1d.push_back(&MindcaTrackSV);
  Extradist1d.push_back(&MinDca);
  Extradist1d.push_back(&MinD0SigSV);
  Extradist1d.push_back(&MinD0SigPV);
  Extradist1d.push_back(&MaxVertexPairQuality);
  Extradist1d.push_back(&MaxdeltaMuZ);
  Extradist1d.push_back(&MaxDca);
  Extradist1d.push_back(&MaxD0SigSV);
  Extradist1d.push_back(&MaxD0SigPV);
  Extradist1d.push_back(&Iso1);
  Extradist1d.push_back(&Iso1Mu1);
  Extradist1d.push_back(&Iso8Mu1);
  Extradist1d.push_back(&FLSignificance);

}


void  DsToPhiPi::doEvent(){ 
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection

  //bool RunB, RunC, RunD, RunE, RunF = 0;
  if(id==1 && Ntp->WhichEra(2017).Contains("RunB") ){ RunB=1;} 
  if(id==1 && Ntp->WhichEra(2017).Contains("RunC") ){ RunC=1;} 
  if(id==1 && Ntp->WhichEra(2017).Contains("RunD") ){ RunD=1;} 
  if(id==1 && Ntp->WhichEra(2017).Contains("RunE") ){ RunE=1;}
  if(id==1 && Ntp->WhichEra(2017).Contains("RunF") ){ RunF=1;}

  value.at(HLTOk) = 0;
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);
    if(HLT.Contains("DoubleMu3_Trk_Tau3mu")/* || HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu") ) && Ntp->HLTDecision(iTrigger) == 1*/){
      value.at(HLTOk)=Ntp->HLTDecision(iTrigger);
    }
  }

  value.at(L1TOk) = 0;
  bool DoubleMuFired(0);
  for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
    TString L1TriggerName = Ntp->L1Name(il1);
    
    if(id==1 && Ntp->WhichEra(2017).Contains("RunB")){
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
    }

    if(id!=1){
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4") && Ntp->L1Decision(il1)) { DoubleMuFired = true; }
      if( random_num>0.3516 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMuFired = true; }
    }

    if(id==1 && (Ntp->WhichEra(2017).Contains("RunC") || Ntp->WhichEra(2017).Contains("RunD") || Ntp->WhichEra(2017).Contains("RunF") || Ntp->WhichEra(2017).Contains("RunE"))){
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
    }
    if (id==1 && Ntp->WhichEra(2018).Contains("Run")){
      if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"))                 DoubleMuFired = Ntp-> L1Decision(il1);
      if(L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2"))                      DoubleMuFired = Ntp-> L1Decision(il1);
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
  if(Ntp->NTwoMuonsTrack()!=0/* && Ntp->NThreeMuons() == 0*/){
    value.at(is2MuTrk) = 1;
  }

  if (value.at(is2MuTrk)==1){
    for(unsigned int i2M=0; i2M < Ntp->NTwoMuonsTrack(); i2M++){
      int tmp_mu1 = Ntp->TwoMuonsTrackMuonIndices(i2M).at(0);
      int tmp_mu2 = Ntp->TwoMuonsTrackMuonIndices(i2M).at(1);
      int tmp_track = Ntp->TwoMuonsTrackTrackIndex(i2M).at(0);
      double tmp_PhiMass = (Ntp->Muon_P4(tmp_mu1)+Ntp->Muon_P4(tmp_mu2)).M();
      
      if (abs(tmp_PhiMass-1.01)<=check_PhiMass || (tmp_PhiMass > .95 && tmp_PhiMass < 1.1)) {
      if (tmp_chisq>Ntp->TwoMuonsTrack_SV_Chi2(i2M)){
	tmp_chisq = Ntp->TwoMuonsTrack_SV_Chi2(i2M);
        check_PhiMass = abs(tmp_PhiMass-1.01);
        
        if (Ntp->Muon_P4(tmp_mu1).Pt() > Ntp->Muon_P4(tmp_mu2).Pt()) {
	  mu1 = tmp_mu1;
	  mu2 = tmp_mu2;
        } else {
          mu1 = tmp_mu2;
          mu2 = tmp_mu1;
        }
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
  pass.at(Mu1dR) = 1; //value.at(Mu1dR) < .03;
  pass.at(Mu2dR) = 1; //value.at(Mu2dR) < .03;
  pass.at(TrkdR) = 1; //value.at(TrkdR) < .03;
  pass.at(Mu1pt) = value.at(Mu1pt) > 3;
  pass.at(Mu2pt) = value.at(Mu2pt) > 3;
  pass.at(Trkpt) = value.at(Trkpt) > 2;

  double wobs=1;
  double w;  
  double w_peak;     

  if(!Ntp->isData()){w = 1; w_peak = /*.083;} //*/.358;}//Ntp->PUReweight(); } //  No weights to data
  else{w=1; w_peak=1;}
  bool status=AnalysisCuts(t,w,wobs);
  if(status){

    int Nvertices(0);
    for(unsigned int l=0; l < Ntp->NSecondaryVertices(); l++){
      TVector3 SVsignalPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(tmp_idx,true),Ntp->Vertex_MatchedPrimaryVertex(tmp_idx,true));
      TVector3 SVfakePV = Ntp->SVPVDirection(Ntp->SecondaryVertexPosition(l),Ntp->Vertex_MatchedPrimaryVertex(tmp_idx,true));
      if(SVfakePV.DeltaR(SVsignalPV) < 1 && (Ntp->Vertex_Signal_KF_pos(tmp_idx,true) - Ntp->SecondaryVertexPosition(l)).Mag() > 0.05){
	    Nvertices++;
      }
    }

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
    double dsEta = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))+
                     Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0))).Eta();
    TLorentzVector DsLV = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1))+
                     Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0)));
    TLorentzVector DsRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(tmp_idx,0,true)+Ntp->Vertex_signal_KF_refittedTracksP4(tmp_idx,1,true)+Ntp->Vertex_signal_KF_refittedTracksP4(tmp_idx,2,true);
    TLorentzVector Muon1LV = Ntp->Muon_P4(mu1);
    TLorentzVector Muon2LV = Ntp->Muon_P4(mu2);


    PhiMassVsDsMass.at(t).Fill(phimass, dsmass);

    DsMass.at(t).Fill(dsmass,w);

    double SumPT1(0), SumPT8(0);
    int NcloseTracksCount(0);
    int TrackIndex(0);
    double dca_temp(999.);
    for(int i =0; i< Ntp->NIsolationTrack(tmp_idx); i++){
      if(Ntp->IsolationTrack_p4(tmp_idx,i).Pt()> 0.5  && sqrt(  pow(Ntp->IsolationTrack_dzSV(tmp_idx,i),2)   +   pow(Ntp->IsolationTrack_dxySV(tmp_idx,i),2)) < 0.03){
	  NcloseTracksCount++;
      }

      if( sqrt(  pow(Ntp->IsolationTrack_dzSV(tmp_idx,i),2)   +   pow(Ntp->IsolationTrack_dxySV(tmp_idx,i),2) ) <  dca_temp){
	dca_temp = sqrt(  pow(Ntp->IsolationTrack_dzSV(tmp_idx,i),2)   +   pow(Ntp->IsolationTrack_dxySV(tmp_idx,i),2));
	TrackIndex = i;
      }
      if(Ntp->IsolationTrack_p4(tmp_idx,i).Pt()> 0.7 && fabs(Ntp->IsolationTrack_dzPV(tmp_idx,i)) < 0.05 && 
	 sqrt(  pow(Ntp->IsolationTrack_dzSV(tmp_idx,i),2)   +   pow(Ntp->IsolationTrack_dxySV(tmp_idx,i),2)) < 0.05){
        if(Ntp->IsolationTrack_p4(tmp_idx,i).DeltaR(DsRefitLV) < 1.0){
	  SumPT1 += Ntp->IsolationTrack_p4(tmp_idx,i).Pt();
	}
        if(Ntp->IsolationTrack_p4(tmp_idx,i).DeltaR(DsRefitLV) < 0.8){
	  SumPT8 += Ntp->IsolationTrack_p4(tmp_idx,i).Pt();
	}
      }
    }

    float MinD0Significance = std::min({Ntp->Vertex_d0sig_reco(tmp_idx,0,true),
	  Ntp->Vertex_d0sig_reco(tmp_idx,1,true),
	  Ntp->Vertex_d0sig_reco(tmp_idx,2,true)});

    float MaxD0Significance = std::max({Ntp->Vertex_d0sig_reco(tmp_idx,0,true),
	  Ntp->Vertex_d0sig_reco(tmp_idx,1,true),
	  Ntp->Vertex_d0sig_reco(tmp_idx,2,true)});

    float MinD0SVSignificance = std::min({Ntp->Vertex_d0sigSV_reco(tmp_idx,0,true),
          Ntp->Vertex_d0sigSV_reco(tmp_idx,1,true),
          Ntp->Vertex_d0sigSV_reco(tmp_idx,2,true)});
  
    float MaxD0SVSignificance = std::max({Ntp->Vertex_d0sigSV_reco(tmp_idx,0,true),
          Ntp->Vertex_d0sigSV_reco(tmp_idx,1,true),
          Ntp->Vertex_d0sigSV_reco(tmp_idx,2,true)});    

    if(id==30){

      DsGenMatch.at(t).Fill(Ntp->DsGenMatch(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0) + Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(1) + Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0)));

    }

    int vertex_idx = tmp_idx;// + Ntp->NThreeMuons();
    TVector3 SVPV = Ntp->SVPVDirection(Ntp->Vertex_Signal_KF_pos(vertex_idx,true),Ntp->Vertex_MatchedPrimaryVertex(vertex_idx,true));
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
      DsEta_peak.at(t).Fill(dsEta,w);

      VertexKFChi2_peak.at(t).Fill(Ntp->Vertex_signal_KF_Chi2(tmp_idx,true));
      SVPVDsDirAngle_peak.at(t).Fill(SVPV.Angle(DsLV.Vect()),w);
      NtracksClose_peak.at(t).Fill(NcloseTracksCount,1);
      NSV_peak.at(t).Fill(Nvertices,1);
      MinMuon_chi2LocalPosition_peak.at(t).Fill(std::min({Ntp->Muon_combinedQuality_chi2LocalPosition(mu1),Ntp->Muon_combinedQuality_chi2LocalPosition(mu2)/*,Ntp->Muon_combinedQuality_chi2LocalPosition(track)*/  }),w);
      if(Ntp->NIsolationTrack(tmp_idx)!=0) MindcaTrackSV_peak.at(t).Fill(sqrt( pow(Ntp->IsolationTrack_dzSV(tmp_idx,TrackIndex),2)   +   pow(Ntp->IsolationTrack_dxySV(tmp_idx,TrackIndex),2)),1);
      MinDca_peak.at(t).Fill(std::min({Ntp->Vertex_DCA12(tmp_idx,true),Ntp->Vertex_DCA23(tmp_idx,true),Ntp->Vertex_DCA31(tmp_idx,true)}),1);
      MinD0SigSV_peak.at(t).Fill(MinD0SVSignificance,1);
      MinD0SigPV_peak.at(t).Fill(MinD0Significance,1);
      MaxVertexPairQuality_peak.at(t).Fill(  std::max({Ntp->Vertex_pair_quality(tmp_idx,0,true),Ntp->Vertex_pair_quality(tmp_idx,1,true),Ntp->Vertex_pair_quality(tmp_idx,2,true)}),1);
      MaxdeltaMuZ_peak.at(t).Fill( std::max({fabs(Ntp->Muon_Poca(mu1).Z()  - Ntp->Muon_Poca(mu2).Z()),fabs(Ntp->Muon_Poca(mu1).Z()  - Ntp->Track_Poca(track).Z()),
	    fabs(Ntp->Muon_Poca(mu1).Z()  - Ntp->Track_Poca(track).Z()),
	    fabs(Ntp->Muon_Poca(mu2).Z()  - Ntp->Track_Poca(track).Z())}),1);
      MaxDca_peak.at(t).Fill(std::max({Ntp->Vertex_DCA12(tmp_idx,true),Ntp->Vertex_DCA23(tmp_idx,true),Ntp->Vertex_DCA31(tmp_idx,true)}),1);
      MaxD0SigSV_peak.at(t).Fill(MaxD0SVSignificance,1);
      MaxD0SigPV_peak.at(t).Fill(MaxD0Significance,1);
      Iso1_peak.at(t).Fill(DsRefitLV.Pt()/  (DsRefitLV.Pt() + SumPT1),1);
      Iso1Mu1_peak.at(t).Fill(Muon1LV.Pt()/  (Muon1LV.Pt() + SumPT1),1);
      Iso8Mu1_peak.at(t).Fill(Muon1LV.Pt()/  (Muon1LV.Pt() + SumPT8),1);
      //FLSignificance_peak.at(t).Fill(( Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(tmp_idx,true),Ntp->Vertex_PrimaryVertex_Covariance(tmp_idx,true),
      //							   Ntp->Vertex_Signal_KF_pos(tmp_idx,true),Ntp->Vertex_Signal_KF_Covariance(tmp_idx,true))),w);

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
      DsEta_sideband.at(t).Fill(dsEta,w);

      VertexKFChi2_sideband.at(t).Fill(Ntp->Vertex_signal_KF_Chi2(tmp_idx,true));
      SVPVDsDirAngle_sideband.at(t).Fill(SVPV.Angle(DsLV.Vect()),w);
      NtracksClose_sideband.at(t).Fill(NcloseTracksCount,1);
      NSV_sideband.at(t).Fill(Nvertices,1);
      MinMuon_chi2LocalPosition_sideband.at(t).Fill(std::min({Ntp->Muon_combinedQuality_chi2LocalPosition(mu1),Ntp->Muon_combinedQuality_chi2LocalPosition(mu2)/*,Ntp->Muon_combinedQuality_chi2LocalPosition(track)*/  }),w);
      if(Ntp->NIsolationTrack(tmp_idx)!=0) MindcaTrackSV_sideband.at(t).Fill(sqrt( pow(Ntp->IsolationTrack_dzSV(tmp_idx,TrackIndex),2)   +   pow(Ntp->IsolationTrack_dxySV(tmp_idx,TrackIndex),2)),1);
      MinDca_sideband.at(t).Fill(std::min({Ntp->Vertex_DCA12(tmp_idx,true),Ntp->Vertex_DCA23(tmp_idx,true),Ntp->Vertex_DCA31(tmp_idx,true)}),1);
      MinD0SigSV_sideband.at(t).Fill(MinD0SVSignificance,1);
      MinD0SigPV_sideband.at(t).Fill(MinD0Significance,1);
      MaxVertexPairQuality_sideband.at(t).Fill(  std::max({Ntp->Vertex_pair_quality(tmp_idx,0,true),Ntp->Vertex_pair_quality(tmp_idx,1,true),Ntp->Vertex_pair_quality(tmp_idx,2,true)}),1);
      MaxdeltaMuZ_sideband.at(t).Fill( std::max({fabs(Ntp->Muon_Poca(mu1).Z()  - Ntp->Muon_Poca(mu2).Z()),fabs(Ntp->Muon_Poca(mu1).Z()  - Ntp->Track_Poca(track).Z()),
            fabs(Ntp->Muon_Poca(mu1).Z()  - Ntp->Track_Poca(track).Z()),
            fabs(Ntp->Muon_Poca(mu2).Z()  - Ntp->Track_Poca(track).Z())}),1);
      MaxDca_sideband.at(t).Fill(std::max({Ntp->Vertex_DCA12(tmp_idx,true),Ntp->Vertex_DCA23(tmp_idx,true),Ntp->Vertex_DCA31(tmp_idx,true)}),1);
      MaxD0SigSV_sideband.at(t).Fill(MaxD0SVSignificance,1);
      MaxD0SigPV_sideband.at(t).Fill(MaxD0Significance,1);
      Iso1_sideband.at(t).Fill(DsRefitLV.Pt()/  (DsRefitLV.Pt() + SumPT1),1);
      Iso1Mu1_sideband.at(t).Fill(Muon1LV.Pt()/  (Muon1LV.Pt() + SumPT1),1);
      Iso8Mu1_sideband.at(t).Fill(Muon1LV.Pt()/  (Muon1LV.Pt() + SumPT8),1);
      //FLSignificance_sideband.at(t).Fill(( Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(tmp_idx,true),Ntp->Vertex_PrimaryVertex_Covariance(tmp_idx,true),
      //                                                             Ntp->Vertex_Signal_KF_pos(tmp_idx,true),Ntp->Vertex_Signal_KF_Covariance(tmp_idx,true))),w);

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
      Ds_eta.at(t).Fill(dsEta,w*w_peak);

      VertexKFChi2.at(t).Fill(Ntp->Vertex_signal_KF_Chi2(tmp_idx,true),w*w_peak);
      SVPVDsDirAngle.at(t).Fill(SVPV.Angle(DsLV.Vect()),w*w_peak);
      NtracksClose.at(t).Fill(NcloseTracksCount,1*w_peak);
      NSV.at(t).Fill(Nvertices,1*w_peak);
      MinMuon_chi2LocalPosition.at(t).Fill(std::min({Ntp->Muon_combinedQuality_chi2LocalPosition(mu1),Ntp->Muon_combinedQuality_chi2LocalPosition(mu2)/*,Ntp->Muon_combinedQuality_chi2LocalPosition(track)*/  }),w*w_peak);
      if(Ntp->NIsolationTrack(tmp_idx)!=0) MindcaTrackSV.at(t).Fill(sqrt( pow(Ntp->IsolationTrack_dzSV(tmp_idx,TrackIndex),2)   +   pow(Ntp->IsolationTrack_dxySV(tmp_idx,TrackIndex),2)),1*w_peak);
      MinDca.at(t).Fill(std::min({Ntp->Vertex_DCA12(tmp_idx,true),Ntp->Vertex_DCA23(tmp_idx,true),Ntp->Vertex_DCA31(tmp_idx,true)}),1*w_peak);
      MinD0SigSV.at(t).Fill(MinD0SVSignificance,1*w_peak);
      MinD0SigPV.at(t).Fill(MinD0Significance,1*w_peak);
      MaxVertexPairQuality.at(t).Fill(  std::max({Ntp->Vertex_pair_quality(tmp_idx,0,true),Ntp->Vertex_pair_quality(tmp_idx,1,true),Ntp->Vertex_pair_quality(tmp_idx,2,true)}),1*w_peak);
      MaxdeltaMuZ.at(t).Fill( std::max({fabs(Ntp->Muon_Poca(mu1).Z()  - Ntp->Muon_Poca(mu2).Z()),fabs(Ntp->Muon_Poca(mu1).Z()  - Ntp->Track_Poca(track).Z()),
            fabs(Ntp->Muon_Poca(mu1).Z()  - Ntp->Track_Poca(track).Z()),
            fabs(Ntp->Muon_Poca(mu2).Z()  - Ntp->Track_Poca(track).Z())}),1*w_peak);
      MaxDca.at(t).Fill(std::max({Ntp->Vertex_DCA12(tmp_idx,true),Ntp->Vertex_DCA23(tmp_idx,true),Ntp->Vertex_DCA31(tmp_idx,true)}),1*w_peak);
      MaxD0SigSV.at(t).Fill(MaxD0SVSignificance,1*w_peak);
      MaxD0SigPV.at(t).Fill(MaxD0Significance,1*w_peak);
      Iso1.at(t).Fill(DsRefitLV.Pt()/  (DsRefitLV.Pt() + SumPT1),1*w_peak);
      Iso1Mu1.at(t).Fill(Muon1LV.Pt()/  (Muon1LV.Pt() + SumPT1),1*w_peak);
      Iso8Mu1.at(t).Fill(Muon1LV.Pt()/  (Muon1LV.Pt() + SumPT8),1*w_peak);
      //FLSignificance.at(t).Fill(( Ntp->FlightLength_significance(Ntp->Vertex_MatchedPrimaryVertex(tmp_idx,true),Ntp->Vertex_PrimaryVertex_Covariance(tmp_idx,true),
      //                                                             Ntp->Vertex_Signal_KF_pos(tmp_idx,true),Ntp->Vertex_Signal_KF_Covariance(tmp_idx,true))),w*w_peak);

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

    //if(RunB){scaleRun.push_back(2265.07);scaleRun.push_back(4586.64);} 
    //if(RunC){scaleRun.push_back(15965.4);scaleRun.push_back(24220.4);}
    //if(RunD){scaleRun.push_back(5952.73);scaleRun.push_back(11303.6);}
    //if(RunE){scaleRun.push_back(10661.2);scaleRun.push_back(19461.1);}
    //if(RunF){scaleRun.push_back(10093.0);scaleRun.push_back(19046.7);}
    //scaleRun.push_back(1122.09);scaleRun.push_back(2246.3);
    scaleRun.push_back(6875.5);scaleRun.push_back(11641.8);

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
    DsEta_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    VertexKFChi2_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    SVPVDsDirAngle_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    NtracksClose_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    NSV_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    MinMuon_chi2LocalPosition_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    MindcaTrackSV_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    MinDca_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    MinD0SigSV_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    MinD0SigPV_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    MaxVertexPairQuality_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    MaxdeltaMuZ_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    MaxDca_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    MaxD0SigSV_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    MaxD0SigPV_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    Iso1_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    Iso1Mu1_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    Iso8Mu1_sideband.at(0).Scale(scaleRun[0]/scaleRun[1]);
    FLSignificance.at(0).Scale(scaleRun[0]/scaleRun[1]);

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
    VertexKFChi2.at(0).Add(&VertexKFChi2_peak.at(0));
    VertexKFChi2.at(0).Add(&VertexKFChi2_sideband.at(0),-1);
    SVPVDsDirAngle.at(0).Add(&SVPVDsDirAngle_peak.at(0));
    SVPVDsDirAngle.at(0).Add(&SVPVDsDirAngle_sideband.at(0),-1);
    NtracksClose.at(0).Add(&NtracksClose_peak.at(0));
    NtracksClose.at(0).Add(&NtracksClose_sideband.at(0),-1);
    NSV.at(0).Add(&NSV_peak.at(0));
    NSV.at(0).Add(&NSV_sideband.at(0),-1);
    MinMuon_chi2LocalPosition.at(0).Add(&MinMuon_chi2LocalPosition_peak.at(0));
    MinMuon_chi2LocalPosition.at(0).Add(&MinMuon_chi2LocalPosition_sideband.at(0),-1);
    MindcaTrackSV.at(0).Add(&MindcaTrackSV_peak.at(0));
    MindcaTrackSV.at(0).Add(&MindcaTrackSV_sideband.at(0),-1);
    MinDca.at(0).Add(&MinDca_peak.at(0));
    MinDca.at(0).Add(&MinDca_sideband.at(0),-1);
    MinD0SigSV.at(0).Add(&MinD0SigSV_peak.at(0));
    MinD0SigSV.at(0).Add(&MinD0SigSV_sideband.at(0),-1);
    MinD0SigPV.at(0).Add(&MinD0SigPV_peak.at(0));
    MinD0SigPV.at(0).Add(&MinD0SigPV_sideband.at(0),-1);
    MaxVertexPairQuality.at(0).Add(&MaxVertexPairQuality_peak.at(0));
    MaxVertexPairQuality.at(0).Add(&MaxVertexPairQuality_sideband.at(0),-1);
    MaxdeltaMuZ.at(0).Add(&MaxdeltaMuZ_peak.at(0));
    MaxdeltaMuZ.at(0).Add(&MaxdeltaMuZ_sideband.at(0),-1);
    MaxDca.at(0).Add(&MaxDca_peak.at(0));
    MaxDca.at(0).Add(&MaxDca_sideband.at(0),-1);
    MaxD0SigSV.at(0).Add(&MaxD0SigSV_peak.at(0));
    MaxD0SigSV.at(0).Add(&MaxD0SigSV_sideband.at(0),-1);
    MaxD0SigPV.at(0).Add(&MaxD0SigPV_peak.at(0));
    MaxD0SigPV.at(0).Add(&MaxD0SigPV_sideband.at(0),-1);
    Iso1.at(0).Add(&Iso1_peak.at(0));
    Iso1.at(0).Add(&Iso1_sideband.at(0),-1);
    Iso1Mu1.at(0).Add(&Iso1Mu1_peak.at(0));
    Iso1Mu1.at(0).Add(&Iso1Mu1_sideband.at(0),-1);
    Iso8Mu1.at(0).Add(&Iso8Mu1_peak.at(0));
    Iso8Mu1.at(0).Add(&Iso8Mu1_sideband.at(0),-1);
    FLSignificance.at(0).Add(&FLSignificance_peak.at(0));
    FLSignificance.at(0).Add(&FLSignificance_sideband.at(0),-1);

    Ds_Pt.at(0).Add(&DsPt_peak.at(0));
    Ds_Pt.at(0).Add(&DsPt_sideband.at(0),-1);
    Ds_P.at(0).Add(&DsP_peak.at(0));
    Ds_P.at(0).Add(&DsP_sideband.at(0),-1);
    Ds_M.at(0).Add(&DsM_peak.at(0));
    Ds_M.at(0).Add(&DsM_sideband.at(0),-1);
    Ds_L.at(0).Add(&DsL_peak.at(0));
    Ds_L.at(0).Add(&DsL_sideband.at(0),-1);
    Ds_eta.at(0).Add(&DsEta_peak.at(0));
    Ds_eta.at(0).Add(&DsEta_sideband.at(0),-1);

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



