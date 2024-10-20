#include "ZTau3MuTauh_PreFC.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
//#include "Logger.h"

ZTau3MuTauh_PreFC::ZTau3MuTauh_PreFC(TString Name_, TString id_):
  Selection(Name_,id_)
{
  // This is a class constructor;
}

ZTau3MuTauh_PreFC::~ZTau3MuTauh_PreFC(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  ZTau3MuTauh_PreFC::Configure(){

  //  This mini tree is for limit extraction

  T3MMiniTree= new TTree("T3MMiniTree","T3MMiniTree");

  T3MMiniTree->Branch("m3m",&m3m);
  T3MMiniTree->Branch("dataMCtype",&dataMCtype);
  T3MMiniTree->Branch("event_weight",&event_weight);
  T3MMiniTree->Branch("m12",&m12);
  T3MMiniTree->Branch("m13",&m13);
  T3MMiniTree->Branch("mDr1",&mDr1);
  T3MMiniTree->Branch("mDr2",&mDr2);
  T3MMiniTree->Branch("LumiScale",&LumiScale);


  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==WhetherDecayFound)        cut.at(WhetherDecayFound)=1;
    if(i==Mu1_Candidate_p)          cut.at(Mu1_Candidate_p)=2.49;
    if(i==Mu1_Candidate_eta)        cut.at(Mu1_Candidate_eta)=2.41;
    if(i==Mu2_Candidate_p)          cut.at(Mu2_Candidate_p)=2.49;
    if(i==Mu2_Candidate_eta)        cut.at(Mu2_Candidate_eta)=2.41;
    if(i==Mu3_Candidate_p)          cut.at(Mu3_Candidate_p)=2.49;
    if(i==Mu3_Candidate_eta)        cut.at(Mu3_Candidate_eta)=2.41;
    if(i==Tau_h_Candidate_p)        cut.at(Tau_h_Candidate_p)=18.0;
    if(i==Tau_h_Candidate_eta)      cut.at(Tau_h_Candidate_eta)=2.41;
    if(i==Mu1_Candidate_recod)      cut.at(Mu1_Candidate_recod)=1;
    if(i==Mu2_Candidate_recod)      cut.at(Mu2_Candidate_recod)=1;
    if(i==Mu3_Candidate_recod)      cut.at(Mu3_Candidate_recod)=1;
    if(i==Tau_h_Candidate_recod)    cut.at(Tau_h_Candidate_recod)=1;
    if(i==L1_TriggerOk)       cut.at(L1_TriggerOk)=1;
    if(i==HLT_TriggerOk)      cut.at(HLT_TriggerOk)=1;
    if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
    if(i==TripletPT)          cut.at(TripletPT)=20;
    if(i==DeepTauJets)        cut.at(DeepTauJets)=1;
    if(i==DeepTauMuons)       cut.at(DeepTauMuons)=1;
    if(i==DeepTauElectrons)   cut.at(DeepTauElectrons)=1;
    if(i==OSCharge)           cut.at(OSCharge)=1;
    if(i==nTaus_pT)           cut.at(nTaus_pT)=17.5;
    if(i==nTaus_eta)          cut.at(nTaus_eta)=2.4;
    if(i==nTaus_dR)           cut.at(nTaus_dR)=1.1;
    if(i==Tau3MuIsolation)    cut.at(Tau3MuIsolation)=0.55;
    if(i==TriggerMatch)       cut.at(TriggerMatch)=1;
    if(i==VisMass)            cut.at(VisMass)=1;


  }

  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    // book here the N-1 and N-0 histrogramms for each cut
    if(i==L1_TriggerOk){
      title.at(i)="L1 Trigger ";
      hlabel="L1 Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_L1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_L1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==HLT_TriggerOk){
      title.at(i)="HLT Trigger ";
      hlabel="HLT Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLT_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLT_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==WhetherDecayFound){
      title.at(i)="$Z \\rightarrow \\tau_{h}, \\tau_{3\\mu}$ decay information found in ntuple";
      hlabel="Decay information found in ntuple ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_WhetherDecayFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_WhetherDecayFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Mu1_Candidate_p){
      title.at(i)=" Whether GEN level $\\mu_{1}$ has $p>2.49 GeV$ ";
      hlabel="$\\mu_{1}$ $p, GeV$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1_Candidate_p_",htitle,160,0.0,80.0,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1_Candidate_p_",htitle,160,0.0,80.0,hlabel,"Events"));
    }
    else if(i==Mu1_Candidate_eta){
      title.at(i)=" Whether GEN level $\\mu_{1}$ has $|\\eta| < 2.41$ ";
      hlabel="$\\mu_{1}$ $|\\eta|$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1_Candidate_eta_",htitle,30,0,3.14,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1_Candidate_eta_",htitle,30,0,3.14,hlabel,"Events"));
    }
    else if(i==Mu2_Candidate_p){
      title.at(i)=" Whether GEN level $\\mu_{2}$ has $p>2.49 GeV$ ";
      hlabel="$\\mu_{2}$ $p, GeV$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2_Candidate_p_",htitle,160,0.0,80.0,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2_Candidate_p_",htitle,160,0.0,80.0,hlabel,"Events"));
    }
    else if(i==Mu2_Candidate_eta){
      title.at(i)=" Whether GEN level $\\mu_{2}$ has $|\\eta| < 2.41$ ";
      hlabel="$\\mu_{1}$ $|\\eta|$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2_Candidate_eta_",htitle,30,0,3.14,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2_Candidate_eta_",htitle,30,0,3.14,hlabel,"Events"));
    }
    else if(i==Mu3_Candidate_p){
      title.at(i)=" Whether GEN level $\\mu_{3}$ has $p>2.49 GeV$ ";
      hlabel="$\\mu_{3}$ $p, GeV$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu3_Candidate_p_",htitle,160,0.0,80.0,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu3_Candidate_p_",htitle,160,0.0,80.0,hlabel,"Events"));
    }
    else if(i==Mu3_Candidate_eta){
      title.at(i)=" Whether GEN level $\\mu_{3}$ has $|\\eta| < 2.41$ ";
      hlabel="$\\mu_{3}$ $|\\eta|$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu3_Candidate_eta_",htitle,30,0,3.14,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu3_Candidate_eta_",htitle,30,0,3.14,hlabel,"Events"));
    }
    else if(i==Tau_h_Candidate_p){
      title.at(i)=" Whether GEN level $\\tau_{h}$ has $pT>18 GeV$ ";
      hlabel="$\\tau_{h}$ $p, GeV$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau_h_Candidate_p_",htitle,40,0.0,80.0,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau_h_Candidate_p_",htitle,40,0.0,80.0,hlabel,"Events"));
    }
    else if(i==Tau_h_Candidate_eta){
      title.at(i)=" Whether GEN level $\\tau_{h}$ has $|\\eta| < 2.41$ ";
      hlabel="$\\tau_{h}$ $|\\eta|$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau_h_Candidate_eta_",htitle,30,0,3.14,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau_h_Candidate_eta_",htitle,30,0,3.14,hlabel,"Events"));
    }
    else if(i==Mu1_Candidate_recod){
      title.at(i)=" Whether $\\mu_{1}$ is reconstructed in MC ";
      hlabel="If $\\mu_{1}$ is reconstructed in MC";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu1_Candidate_recod_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu1_Candidate_recod_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Mu2_Candidate_recod){
      title.at(i)=" Whether $\\mu_{2}$ is reconstructed in MC ";
      hlabel="If $\\mu_{2}$ is reconstructed in MC";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu2_Candidate_recod_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu2_Candidate_recod_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Mu3_Candidate_recod){
      title.at(i)=" Whether $\\mu_{3}$ is reconstructed in MC ";
      hlabel="If $\\mu_{3}$ is reconstructed in MC";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mu3_Candidate_recod_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mu3_Candidate_recod_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Tau_h_Candidate_recod){
      title.at(i)=" Whether $\\tau_{h}$ is reconstructed in MC ";
      hlabel="If $\\tau_{h}$ is reconstructed in MC";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau_h_Candidate_recod_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau_h_Candidate_recod_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==nTaus_pT){
      title.at(i)=" At least one $\\tau_{h}$, $pT>17.5 GeV$";
      hlabel="pT, GeV";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nTaus_pT_",htitle,80,0.0,40,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nTaus_pT_",htitle,80,0.0,40,hlabel,"Events"));
    }
    else if(i==nTaus_eta){
      title.at(i)=" At least one $\\tau_{h}$, $|\\eta| < 2.4$";
      hlabel="$|\\eta|$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nTaus_eta_",htitle,20,2.0,3.0,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nTaus_eta_",htitle,20,2.0,3.0,hlabel,"Events"));
    }
    else if(i==nTaus_dR){
      title.at(i)=" At least one $\\tau_{h}$, $\\Delta R (\\tau_{h}-3\\mu) >$ 1.1";
      hlabel="$\\Delta R (\\tau_{h}-3\\mu) $";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nTaus_dR_",htitle,20,0.0,3.14,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nTaus_dR_",htitle,20,0.0,3.14,hlabel,"Events"));
    }
    else if(i==DeepTauJets){
      title.at(i)=" At least one $\\tau_{h}$ pass DeepTauVsJets (loose WP) ";
      hlabel=" deep tau jets";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_DeepTauJets_",htitle,4,-0.5,3.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_DeepTauJets_",htitle,4,-0.5,3.5,hlabel,"Events"));
    }

    else if(i==DeepTauMuons){
      title.at(i)="$\\tau_{h}$ pass DeepTauVsMuons (loose WP) ";
      hlabel=" deep tau muons";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_DeepTauMuons_",htitle,4,-0.5,3.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_DeepTauMuons_",htitle,4,-0.5,3.5,hlabel,"Events"));
    }

    else if(i==DeepTauElectrons){
      title.at(i)="$\\tau_{h}$ pass DeepTauVsElectrons (loose WP) ";
      hlabel=" deep tau electrons";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_DeepTauElectrons_",htitle,4,-0.5,3.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_DeepTauElectrons_",htitle,4,-0.5,3.5,hlabel,"Events"));
    }

    else if(i==TripletPT){
      title.at(i)="pT(3$\\mu$)  $>$ 18 GeV";
      htitle=title.at(i);
      hlabel="pT(#tau_{3#mu}) , GeV ";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TripletPT_",htitle,50,5,80,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TripletPT_",htitle,50,5,80,hlabel,"Events"));
    }

    else if(i==SignalCandidate){
      title.at(i)="At least one $\\tau_{3\\mu}$ candidate (3,3,2 GeV,  $|\\eta| < 2.4$, dz($\\mu_{i} , \\mu_{j}$)$<$0.5, dR($\\mu_{i} , \\mu_{j}$)$<$0.8, $\\Sigma \\mu_{charge}$ = +-1)";
      htitle=title.at(i);
      hlabel="N $3\\mu$ candidates";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,5,-0.5,4.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }
    else if(i==OSCharge){
      title.at(i)="Charge $\\tau_{h}$ * $\\tau_{3\\mu}$ =  -1; ";
      title.at(i)+=" (at least one)";
      htitle=title.at(i);
      hlabel="Opposite charge? ";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OSCharge_",htitle,4,-0.5,3.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OSCharge_",htitle,4,-0.5,3.5,hlabel,"Events"));
    }
    else if(i==Tau3MuIsolation){
      title.at(i)="$ 3\\mu $ Relative Isolation  $ > $ 0.55";
      //      title.at(i)+= cut.at(MuonIsolation);
      htitle=title.at(i);
      hlabel="I(3#mu)= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T})";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau3MuIsolation_",htitle,50,0,1.1,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau3MuIsolation_",htitle,50,0,1.1,hlabel,"Events"));
    }

    else if(i==VisMass){
      title.at(i)="50 GeV $< M(\\tau(h) + \\tau(3\\mu))  < $ 100 GeV";
      htitle=title.at(i);
      hlabel="M_{#tau(h) - #tau(3#mu)}, GeV (visible mass)";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_VisMass_",htitle,70,0,150,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_VisMass_",htitle,70,0,150,hlabel,"Events"));
    }
    else if(i==TriggerMatch){
      title.at(i)="Selected (by $\\chi^2$) 3$\\mu$ matched to trg";
      hlabel="Trigger Matched ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }



  }
  // Setup NPassed Histogams

  NumberOfTaus=HConfig.GetTH1D(Name+"_NumberOfTaus","NumberOfTaus",5,-0.5,4.5,"Number of #tau ","Events");



  Tau3MuRelativeIsolation=HConfig.GetTH1D(Name+"_Tau3MuRelativeIsolation","Tau3MuRelativeIsolation",50,0.,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T})","#Delta R < 0.4 ");
  TauHDecayMode=HConfig.GetTH1D(Name+"_TauHDecayMode","TauHDecayMode",12,-0.5,11.5,"HPS #tau_{h} decay mode","Events");
  VisibleDiTauMass=HConfig.GetTH1D(Name+"_VisibleDiTauMass","VisibleDiTauMass",70,0.,150,"M_{#tau(h) - #tau(3#mu)}, GeV (visible mass)","Events");
  MTT=HConfig.GetTH1D(Name+"_MTT","MTT",70,0.,140,"M_{#tau(h) - #tau(3#mu)}, GeV (collinear approximation)","Events");
  TripletMass=HConfig.GetTH1D(Name+"_TripletMass","TripletMass",30,1.1,2.2,"M_{3#mu}, GeV","Events");
  PairMass_OppositeSign_dR12=HConfig.GetTH1D(Name+"_PairMass_OppositeSign_dR12","PairMass_OppositeSign_dR12",40,0.2,2.,"M_{1}, GeV (OS - SS dR sorted)","Events");
  PairMass_OppositeSign_dR13=HConfig.GetTH1D(Name+"_PairMass_OppositeSign_dR13","PairMass_OppositeSign_dR13",40,0.2,2.,"M_{2}, GeV (OS - SS dR sorted)","Events");

  TripletPt=HConfig.GetTH1D(Name+"_TripletPt","TripletPt",50,2,80,"pT(3#mu), GeV ","Events");
  OppositeTauPt=HConfig.GetTH1D(Name+"_OppositeTauPt","OppositeTauPt",50,2,40,"pT(#tau), GeV ","Events");

 
  matched_pdgId=HConfig.GetTH1D(Name+"_matched_pdgId","matched_pdgId",25,-0.5,24.5,"pdgID MC matched","Events");
  matched_dR=HConfig.GetTH1D(Name+"_matched_dR","matched_dR",50,-0.1,0.5,"#Delta R(MC-RECO) Object opposite to #tau_{3#mu}","Events");

  Muon1DRToTruth=HConfig.GetTH1D(Name+"_Muon1DRToTruth","Muon1DRToTruth",20,0,0.02,"#Delta R (reco - mc) #mu_{1} ","Events");
  Muon2DRToTruth=HConfig.GetTH1D(Name+"_Muon2DRToTruth","Muon2DRToTruth",20,0,0.02,"#Delta R (reco - mc) #mu_{2} ","Events");
  Muon3DRToTruth=HConfig.GetTH1D(Name+"_Muon3DRToTruth","Muon3DRToTruth",20,0,0.02,"#Delta R (reco - mc) #mu_{3} ","Events");
  dR_betweenTruth_VisibleTaus=HConfig.GetTH1D(Name+"_dR_betweenTruth_VisibleTaus","dR_betweenTruth_VisibleTaus",20,0,0.5,"#Delta R Truth #tau's prods ","Events");
  
  dR_betweenTruth_NeutrinoGuess=HConfig.GetTH1D(Name+"_dR_betweenTruth_NeutrinoGuess","dR_betweenTruth_NeutrinoGuess",20,0,0.5,"#Delta R: Truth to Guessed #nu","Events");
  dR_betweenTruth_Tau=HConfig.GetTH1D(Name+"_dR_betweenTruth_Tau","dR_betweenTruth_Tau",20,0,0.5,"#Delta R: Truth to Guessed #tau","Events");
  Z_Pt=HConfig.GetTH1D(Name+"_Z_Pt","Z_Pt",50,0,70,"Z_{pT}","Events");
  OS_vs_3mu_trigger=HConfig.GetTH2D(Name+"_OS_vs_3mu_trigger","OS_vs_3mu_trigger",2,-0.5,1.5,2,-0.5,1.5,"Whether 3mu Triggered","Whether OS #tau Triggered");
  
  
  Selection_Cut_3mu_Pt=HConfig.GetTH1D(Name+"_Selection_Cut_3mu_Pt","Selection_Cut_3mu_Pt",100,0,50.0,"3#mu p_{T}, GeV","Events");
  Selection_Cut_3mu_Rel_Iso=HConfig.GetTH1D(Name+"_Selection_Cut_3mu_Rel_Iso","Selection_Cut_3mu_Rel_Iso",50,0,1.1,"3 #mu Relative Isolation, p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T})","Events");
  Selection_Cut_tauh_Pt=HConfig.GetTH1D(Name+"_Selection_Cut_tauh_Pt","Selection_Cut_tauh_Pt",100,0,50.0,"#tau_{h} p_{T}, GeV","Events");
  Selection_Cut_tauh_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_tauh_Eta","Selection_Cut_tauh_Eta",30,0,3.0,"#tau_{h} |#eta|","Events");
  Selection_Cut_tauh_DeltaR_3mu=HConfig.GetTH1D(Name+"_Selection_Cut_tauh_DeltaR_3mu","Selection_Cut_tauh_DeltaR_3mu",60,0,1.2,"#Delta R (#tau_{h}-3#mu)","Events");
  Selection_Cut_Vis_InvM=HConfig.GetTH1D(Name+"_Selection_Cut_Vis_InvM","Selection_Cut_Vis_InvM",75,0,150.0,"M_{#tau(h) + #tau(3#mu)}, GeV (visible mass)","Events");
  
  Selection_Cut_Mu1_P=HConfig.GetTH1D(Name+"_Selection_Cut_Mu1_P","Selection_Cut_Mu1_P",200,0.0,100.0,"#mu_{1} p, GeV","Events");
  Selection_Cut_Mu1_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_Mu1_Eta","Selection_Cut_Mu1_Eta",30,0,3.14,"#mu_{1} |#eta|","Events");
  Selection_Cut_Mu1_p_eta_before=HConfig.GetTH2D(Name+"_Selection_Cut_Mu1_p_eta_before","Selection_Cut_Mu1_p_eta_before",200,0.0,100.0,100,0,5.0,"#mu_{1} p, GeV","#mu_{1} |#eta|");
  Selection_Cut_Mu1_p_eta_after=HConfig.GetTH2D(Name+"_Selection_Cut_Mu1_p_eta_after","Selection_Cut_Mu1_p_eta_after",200,0.0,100.0,100,0,5.0,"#mu_{1} p, GeV","#mu_{1} |#eta|");
  Selection_Cut_Mu1_p_eta_after_reco=HConfig.GetTH2D(Name+"_Selection_Cut_Mu1_p_eta_after_reco","Selection_Cut_Mu1_p_eta_after_reco",200,0.0,100.0,100,0,5.0,"#mu_{1} p, GeV","#mu_{1} |#eta|");
  
  Selection_Cut_Mu2_P=HConfig.GetTH1D(Name+"_Selection_Cut_Mu2_P","Selection_Cut_Mu2_P",200,0.0,100.0,"#mu_{2} p, GeV","Events");
  Selection_Cut_Mu2_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_Mu2_Eta","Selection_Cut_Mu2_Eta",30,0,3.14,"#mu_{2} |#eta|","Events");
  Selection_Cut_Mu2_p_eta_before=HConfig.GetTH2D(Name+"_Selection_Cut_Mu2_p_eta_before","Selection_Cut_Mu2_p_eta_before",200,0.0,100.0,100,0,5.0,"#mu_{2} p, GeV","#mu_{2} |#eta|");
  Selection_Cut_Mu2_p_eta_after=HConfig.GetTH2D(Name+"_Selection_Cut_Mu2_p_eta_after","Selection_Cut_Mu2_p_eta_after",200,0.0,100.0,100,0,5.0,"#mu_{2} p, GeV","#mu_{2} |#eta|");
  Selection_Cut_Mu2_p_eta_after_reco=HConfig.GetTH2D(Name+"_Selection_Cut_Mu2_p_eta_after_reco","Selection_Cut_Mu2_p_eta_after_reco",200,0.0,100.0,100,0,5.0,"#mu_{2} p, GeV","#mu_{2} |#eta|");
  
  Selection_Cut_Mu3_P=HConfig.GetTH1D(Name+"_Selection_Cut_Mu3_P","Selection_Cut_Mu3_P",200,0.0,100.0,"#mu_{3} p, GeV","Events");
  Selection_Cut_Mu3_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_Mu3_Eta","Selection_Cut_Mu3_Eta",100,0,5.0,"#mu_{3} |#eta|","Events");
  Selection_Cut_Mu3_p_eta_before=HConfig.GetTH2D(Name+"_Selection_Cut_Mu3_p_eta_before","Selection_Cut_Mu3_p_eta_before",200,0.0,100.0,100,0,5.0,"#mu_{3} p, GeV","#mu_{3} |#eta|");
  Selection_Cut_Mu3_p_eta_after=HConfig.GetTH2D(Name+"_Selection_Cut_Mu3_p_eta_after","Selection_Cut_Mu3_p_eta_after",200,0.0,100.0,100,0,5.0,"#mu_{3} p, GeV","#mu_{3} |#eta|");
  Selection_Cut_Mu3_p_eta_after_reco=HConfig.GetTH2D(Name+"_Selection_Cut_Mu3_p_eta_after_reco","Selection_Cut_Mu3_p_eta_after_reco",200,0.0,100.0,100,0,5.0,"#mu_{3} p, GeV","#mu_{3} |#eta|");
  
  Selection_Cut_h_Pt=HConfig.GetTH1D(Name+"_Selection_Cut_h_Pt","Selection_Cut_h_Pt",40,0.0,80.0,"#tau_{h}, p_{T}, GeV","Events");
  Selection_Cut_h_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_h_Eta","Selection_Cut_h_Eta",30,0,3.14,"#tau_{h}, |#eta|","Events");
  Selection_Cut_h_pt_eta_before=HConfig.GetTH2D(Name+"_Selection_Cut_h_pt_eta_before","Selection_Cut_h_pt_eta_before",200,0.0,100.0,100,0,5.0,"#tau_{h} pT, GeV","#tau_{h} |#eta|");
  Selection_Cut_h_pt_eta_after=HConfig.GetTH2D(Name+"_Selection_Cut_h_pt_eta_after","Selection_Cut_h_pt_eta_after",200,0.0,100.0,100,0,5.0,"#tau_{h} pT, GeV","#tau_{h} |#eta|");
  Selection_Cut_h_pt_eta_after_reco=HConfig.GetTH2D(Name+"_Selection_Cut_h_pt_eta_after_reco","Selection_Cut_h_pt_eta_after_reco",200,0.0,100.0,100,0,5.0,"#tau_{h} pT, GeV","#tau_{h} |#eta|");
  
  Selection_Cut_Mu1_dR=HConfig.GetTH1D(Name+"Selection_Cut_Mu1_dR","Selection_Cut_Mu1_dR",200,0,0.002,"#Delta R","Events");
  Selection_Cut_Mu2_dR=HConfig.GetTH1D(Name+"Selection_Cut_Mu2_dR","Selection_Cut_Mu2_dR",200,0,0.002,"#Delta R","Events");
  Selection_Cut_Mu3_dR=HConfig.GetTH1D(Name+"Selection_Cut_Mu3_dR","Selection_Cut_Mu3_dR",200,0,0.002,"#Delta R","Events");
  Selection_Cut_h_dR=HConfig.GetTH1D(Name+"Selection_Cut_h_dR","Selection_Cut_h_dR",200,0,0.003,"#Delta R","Events");
  Selection_Cut_Mu1_dR_large_scale=HConfig.GetTH1D(Name+"Selection_Cut_Mu1_dR_large_scale","Selection_Cut_Mu1_dR_large_scale",100,0,2.0,"#Delta R","Events");
  Selection_Cut_Mu2_dR_large_scale=HConfig.GetTH1D(Name+"Selection_Cut_Mu2_dR_large_scale","Selection_Cut_Mu2_dR_large_scale",100,0,2.0,"#Delta R","Events");
  Selection_Cut_Mu3_dR_large_scale=HConfig.GetTH1D(Name+"Selection_Cut_Mu3_dR_large_scale","Selection_Cut_Mu3_dR_large_scale",100,0,2.0,"#Delta R","Events");
  Selection_Cut_h_dR_large_scale=HConfig.GetTH1D(Name+"Selection_Cut_h_dR_large_scale","Selection_Cut_h_dR_large_scale",100,0,2.0,"#Delta R","Events");
  
  Selection_Cut_RecoMu_P=HConfig.GetTH1D(Name+"_Selection_Cut_RecoMu_P","Selection_Cut_RecoMu_P",100,0.0,5.0,"#mu p, GeV","Events");
  Selection_Cut_RecoMu_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_RecoMu_Eta","Selection_Cut_RecoMu_Eta",50,2,3.0,"#mu |#eta|","Events");
  Selection_Cut_RecoH_Pt=HConfig.GetTH1D(Name+"_Selection_Cut_RecoH_Pt","Selection_Cut_RecoH_Pt",100,10.0,20.0,"p_{T}, GeV","Events");
  Selection_Cut_RecoH_Eta=HConfig.GetTH1D(Name+"_Selection_Cut_RecoH_Eta","Selection_Cut_RecoH_Eta",50,2,3.0,"|#eta|","Events");

  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms



  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  ZTau3MuTauh_PreFC::Store_ExtraDist(){ 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output

  Extradist1d.push_back(&NumberOfTaus);
  Extradist1d.push_back(&Tau3MuRelativeIsolation);
  Extradist1d.push_back(&TauHDecayMode);
  Extradist1d.push_back(&VisibleDiTauMass);
  Extradist1d.push_back(&MTT);
  Extradist1d.push_back(&TripletMass);

  Extradist1d.push_back(&PairMass_OppositeSign_dR12);
  Extradist1d.push_back(&PairMass_OppositeSign_dR13);


  Extradist1d.push_back(&matched_pdgId);
  Extradist1d.push_back(&matched_dR);


  Extradist1d.push_back(&Muon1DRToTruth);
  Extradist1d.push_back(&Muon2DRToTruth);
  Extradist1d.push_back(&Muon3DRToTruth);

  Extradist1d.push_back(&dR_betweenTruth_VisibleTaus);

  Extradist1d.push_back(&TripletPt);
  Extradist1d.push_back(&OppositeTauPt);
  
  Extradist1d.push_back(&dR_betweenTruth_NeutrinoGuess);
  Extradist1d.push_back(&dR_betweenTruth_Tau);
  Extradist1d.push_back(&Z_Pt);
  Extradist2d.push_back(&OS_vs_3mu_trigger);
  
  Extradist1d.push_back(&Selection_Cut_3mu_Pt);
  Extradist1d.push_back(&Selection_Cut_3mu_Rel_Iso);
  Extradist1d.push_back(&Selection_Cut_tauh_Pt);
  Extradist1d.push_back(&Selection_Cut_tauh_Eta);
  Extradist1d.push_back(&Selection_Cut_tauh_DeltaR_3mu);
  Extradist1d.push_back(&Selection_Cut_Vis_InvM);
  
  Extradist1d.push_back(&Selection_Cut_Mu1_P);
  Extradist1d.push_back(&Selection_Cut_Mu1_Eta);
  Extradist1d.push_back(&Selection_Cut_Mu2_P);
  Extradist1d.push_back(&Selection_Cut_Mu2_Eta);
  Extradist1d.push_back(&Selection_Cut_Mu3_P);
  Extradist1d.push_back(&Selection_Cut_Mu3_Eta);
  Extradist1d.push_back(&Selection_Cut_h_Pt);
  Extradist1d.push_back(&Selection_Cut_h_Eta);
  
  Extradist2d.push_back(&Selection_Cut_Mu1_p_eta_before);
  Extradist2d.push_back(&Selection_Cut_Mu1_p_eta_after);
  Extradist2d.push_back(&Selection_Cut_Mu1_p_eta_after_reco);
  Extradist2d.push_back(&Selection_Cut_Mu2_p_eta_before);
  Extradist2d.push_back(&Selection_Cut_Mu2_p_eta_after);
  Extradist2d.push_back(&Selection_Cut_Mu2_p_eta_after_reco);
  Extradist2d.push_back(&Selection_Cut_Mu3_p_eta_before);
  Extradist2d.push_back(&Selection_Cut_Mu3_p_eta_after);
  Extradist2d.push_back(&Selection_Cut_Mu3_p_eta_after_reco);
  Extradist2d.push_back(&Selection_Cut_h_pt_eta_before);
  Extradist2d.push_back(&Selection_Cut_h_pt_eta_after);
  Extradist2d.push_back(&Selection_Cut_h_pt_eta_after_reco);
  
  Extradist1d.push_back(&Selection_Cut_Mu1_dR);
  Extradist1d.push_back(&Selection_Cut_Mu2_dR);
  Extradist1d.push_back(&Selection_Cut_Mu3_dR);
  Extradist1d.push_back(&Selection_Cut_h_dR);
  
  Extradist1d.push_back(&Selection_Cut_Mu1_dR_large_scale);
  Extradist1d.push_back(&Selection_Cut_Mu2_dR_large_scale);
  Extradist1d.push_back(&Selection_Cut_Mu3_dR_large_scale);
  Extradist1d.push_back(&Selection_Cut_h_dR_large_scale);
  
  Extradist1d.push_back(&Selection_Cut_RecoMu_P);
  Extradist1d.push_back(&Selection_Cut_RecoMu_Eta);
  Extradist1d.push_back(&Selection_Cut_RecoH_Pt);
  Extradist1d.push_back(&Selection_Cut_RecoH_Eta);



}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  ZTau3MuTauh_PreFC::doEvent(){ 

  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection
  
  
  bool HLTOk(false);
  bool L1Ok(false);
  bool DoubleMu0Fired(false);
  bool DoubleMu4Fired(false);
  bool DoubleMuFired(false);
  bool TripleMuFired(false);
  bool randomFailed(false);
  bool HLT_OppositeSide(false);

  
  
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLTName = Ntp->HLTName(iTrigger);
    //std::cout<<"HLT:   "  << Ntp->HLTName(iTrigger)  << "  fires  "<< Ntp->HLTDecision(iTrigger)<< std::endl;
    if(HLTName.Contains("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_v") && Ntp->HLTDecision(iTrigger) ) { HLTOk = true;}
    if(HLTName.Contains("HLT_IsoMu24") && Ntp->HLTDecision(iTrigger) ) { HLT_OppositeSide = true;}
  }
  
  OS_vs_3mu_trigger.at(t).Fill(HLTOk,HLT_OppositeSide,1 );
  
  random_num = rndm.Rndm();
  
  for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
    TString L1TriggerName = Ntp->L1Name(il1);
    //std::cout<<" l1 name  "<< Ntp->L1Name(il1) << std::endl;
    
    if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4") && Ntp->L1Decision(il1)) { DoubleMu0Fired = true; }
    if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9") && Ntp->L1Decision(il1)) { TripleMuFired = true; }
    if( id==1 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMu4Fired = true; }
    if( id!=1 && random_num>0.30769 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMu4Fired = true;}
    if( id!=1 && random_num<0.30769 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) {
      randomFailed = true;
    }
  }

  value.at(L1_TriggerOk)=0;
  value.at(HLT_TriggerOk)=0;  
  if (DoubleMu0Fired || DoubleMu4Fired) {DoubleMuFired = true;}
  if (DoubleMuFired || TripleMuFired) L1Ok = true;



  value.at(L1_TriggerOk)=(L1Ok);
  pass.at(L1_TriggerOk)=(value.at(L1_TriggerOk)==cut.at(L1_TriggerOk));
  
  value.at(HLT_TriggerOk)=(HLTOk);
  pass.at(HLT_TriggerOk)=(value.at(HLT_TriggerOk)==cut.at(HLT_TriggerOk));



  value.at(SignalCandidate) = Ntp->NThreeMuons();

  int  signal_idx=-1;
  double min_chi2(99.);

  for(int i_idx =0; i_idx < Ntp->NThreeMuons(); i_idx++){
    if(Ntp->Vertex_Signal_KF_Chi2(i_idx) < min_chi2){
      min_chi2 = Ntp->Vertex_Signal_KF_Chi2(i_idx);
      signal_idx = i_idx;
    }
  }
  
  //  std::cout << "Test 1. " << std::endl;
  bool WhetherSignalMC = id==210||id==210231||id==210232||id==210233;
  if(WhetherSignalMC){
  
  int Whether_decay_found(0);
  int TausFromZ_Count(0);
  int tau_3mu_idx(-1);
  int tau_mu_idx(-1);
  int tau_e_idx(-1);
  int tau_h_idx(-1);
  for(unsigned int imc =0; imc< Ntp->NMCParticles(); imc++){
    if(abs(Ntp->MCParticle_pdgid(imc)) == 15 && Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(imc) ) == 23){
      TausFromZ_Count++;
          
          //correcting for cases where tau decays multiple times
          int ChildIdx = imc;
          bool Whether_Another_Tau_Decay_In_Chain(false);
          for(unsigned int i =0; i< Ntp->MCParticle_childpdgid(ChildIdx).size(); i++){
            if(abs(Ntp->MCParticle_childpdgid(ChildIdx).at(i))==15) Whether_Another_Tau_Decay_In_Chain=true;
          }
          while (Whether_Another_Tau_Decay_In_Chain){
            Whether_Another_Tau_Decay_In_Chain=false;
            for(unsigned int i =0; i< Ntp->MCParticle_childpdgid(ChildIdx).size(); i++){
              if(abs(Ntp->MCParticle_childpdgid(ChildIdx).at(i))==15){
                ChildIdx = Ntp->MCParticle_childidx(ChildIdx).at(i);
                Whether_Another_Tau_Decay_In_Chain=true;
              }
            }
          }
          
          int nmuons_temp(0);
          int nelectrons_temp(0);
          for(unsigned int i =0; i< Ntp->MCParticle_childpdgid(ChildIdx).size(); i++){
            //std::cout<<"Child "<< i <<" of tau is :  "<< Ntp->MCParticle_childpdgid(ChildIdx).at(i) << std::endl;
            if(abs(Ntp->MCParticle_childpdgid(ChildIdx).at(i))==13){
              nmuons_temp++;
            }
            if(abs(Ntp->MCParticle_childpdgid(ChildIdx).at(i))==11){
              nelectrons_temp++;
            }
          }
          
          if(nmuons_temp==3){
            tau_3mu_idx=ChildIdx;
          }
          else if(nmuons_temp==1){
            tau_mu_idx=ChildIdx;
          }
          else if(nelectrons_temp==1){
            tau_e_idx=ChildIdx;
          }
          else{
            tau_h_idx=ChildIdx;
          }
          
          
      }
  }
    
  if(tau_3mu_idx>-0.5 && tau_h_idx>-0.5) Whether_decay_found=1;
  int Whether_decay_found_temp(0);
  if(tau_3mu_idx>-0.5) Whether_decay_found_temp=1;
  
  //std::cout<<"  No of taus from Z:  "<< TausFromZ_Count << std::endl;
  //if(!Whether_decay_found_temp&&id==210233)      Ntp->printMCDecayChainOfEvent(true,true,true,false);
  
  TLorentzVector Mu1_LV;
  TLorentzVector Mu2_LV;
  TLorentzVector Mu3_LV;
  TLorentzVector Tau_nu_LV;
  TLorentzVector Tau_h_LV;
  if(Whether_decay_found==1){
    std::vector<int> Sorted_MC_Indices = Ntp->SortedPtMuons_MC(Ntp->MCParticle_childidx(tau_3mu_idx));
    
    Mu1_LV=Ntp->MCParticle_p4(Sorted_MC_Indices.at(0));
    Mu2_LV=Ntp->MCParticle_p4(Sorted_MC_Indices.at(1));
    Mu3_LV=Ntp->MCParticle_p4(Sorted_MC_Indices.at(2));
    
    for(int i = 0; i < Ntp->MCParticle_childpdgid(tau_h_idx).size(); i++){
      if(abs(Ntp->MCParticle_childpdgid(tau_h_idx).at(i))==16){
        Tau_nu_LV=Ntp->MCParticle_p4(Ntp->MCParticle_childidx(tau_h_idx).at(i));
      }
    }
    Tau_h_LV = Ntp->MCParticle_p4(tau_h_idx) - Tau_nu_LV;
  }
  
  
  value.at(WhetherDecayFound)=(Whether_decay_found);
  pass.at(WhetherDecayFound)=(value.at(WhetherDecayFound)==cut.at(WhetherDecayFound));
  
  value.at(Mu1_Candidate_p)=(Mu1_LV.Vect().Mag());
  pass.at(Mu1_Candidate_p)=value.at(Mu1_Candidate_p)>cut.at(Mu1_Candidate_p);
  
  value.at(Mu1_Candidate_eta)=(abs(Mu1_LV.Eta()));
  pass.at(Mu1_Candidate_eta)=value.at(Mu1_Candidate_eta)<cut.at(Mu1_Candidate_eta);
  
  value.at(Mu2_Candidate_p)=(Mu2_LV.Vect().Mag());
  pass.at(Mu2_Candidate_p)=value.at(Mu2_Candidate_p)>cut.at(Mu2_Candidate_p);
  
  value.at(Mu2_Candidate_eta)=(abs(Mu2_LV.Eta()));
  pass.at(Mu2_Candidate_eta)=value.at(Mu2_Candidate_eta)<cut.at(Mu2_Candidate_eta);
  
  value.at(Mu3_Candidate_p)=(Mu3_LV.Vect().Mag());
  pass.at(Mu3_Candidate_p)=value.at(Mu3_Candidate_p)>cut.at(Mu3_Candidate_p);
  
  value.at(Mu3_Candidate_eta)=(abs(Mu3_LV.Eta()));
  pass.at(Mu3_Candidate_eta)=value.at(Mu3_Candidate_eta)<cut.at(Mu3_Candidate_eta);
  
  value.at(Tau_h_Candidate_p)=(Tau_h_LV.Pt());
  pass.at(Tau_h_Candidate_p)=value.at(Tau_h_Candidate_p)>cut.at(Tau_h_Candidate_p);
  
  value.at(Tau_h_Candidate_eta)=(abs(Tau_h_LV.Eta()));
  pass.at(Tau_h_Candidate_eta)=value.at(Tau_h_Candidate_eta)<cut.at(Tau_h_Candidate_eta);
  
  double dR1_max(99.0);
  double dR2_max(99.0);
  double dR3_max(99.0);
  double dR4_max(99.0);
  
  
  
  for(unsigned int imu=0; imu < Ntp->NMuons(); imu++)
    {
      TLorentzVector T1=Ntp->Muon_P4(imu);
      if(Ntp->Muon_P4(imu).DeltaR(Mu1_LV)<dR1_max){
        dR1_max=Ntp->Muon_P4(imu).DeltaR(Mu1_LV);
      }
      if(Ntp->Muon_P4(imu).DeltaR(Mu2_LV)<dR2_max){
        dR2_max=Ntp->Muon_P4(imu).DeltaR(Mu2_LV);
      }
      if(Ntp->Muon_P4(imu).DeltaR(Mu3_LV)<dR3_max){
        dR3_max=Ntp->Muon_P4(imu).DeltaR(Mu3_LV);
      }
      
      Selection_Cut_RecoMu_P.at(t).Fill(Ntp->Muon_P4(imu).Vect().Mag(),1 );
      Selection_Cut_RecoMu_Eta.at(t).Fill(abs(Ntp->Muon_P4(imu).Eta()),1 );

    }
  
  for(unsigned int itau=0; itau < Ntp->NTaus(); itau++)
    {
      if(Ntp->Tau_P4(itau).DeltaR(Tau_h_LV)<dR4_max){
        dR4_max=Ntp->Tau_P4(itau).DeltaR(Tau_h_LV);
      }
      
      Selection_Cut_RecoH_Pt.at(t).Fill(Ntp->Tau_P4(itau).Vect().Mag(),1 );
      Selection_Cut_RecoH_Eta.at(t).Fill(abs(Ntp->Tau_P4(itau).Eta()),1 );
    }
  
  value.at(Mu1_Candidate_recod)=(dR1_max<0.01);
  pass.at(Mu1_Candidate_recod)=value.at(Mu1_Candidate_recod);
  value.at(Mu2_Candidate_recod)=(dR2_max<0.01);
  pass.at(Mu2_Candidate_recod)=value.at(Mu2_Candidate_recod);
  value.at(Mu3_Candidate_recod)=(dR3_max<0.01);
  pass.at(Mu3_Candidate_recod)=value.at(Mu3_Candidate_recod);
  
  value.at(Tau_h_Candidate_recod)=(dR4_max<0.05);
  pass.at(Tau_h_Candidate_recod)=value.at(Tau_h_Candidate_recod);
  
  // This is to print out selected event content
  if(id==210233){
          
          //std::cout<<"Whether tau e decay found: "<< Whether_decay_found << ". Event id: "<<id<< std::endl;
          
          if(Whether_decay_found){
                  /*
                  std::cout<<"------------------------------- "<< std::endl;
                  std::cout<<"Event Content "<< std::endl;
                  std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Mu1_LV) << std::endl;
                  std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Mu2_LV) << std::endl;
                  std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Mu3_LV) << std::endl;
                  std::cout<<" idx4:  "<<Ntp->getMatchTruthIndex(Tau_h_LV) << std::endl;
                  Ntp->printMCDecayChainOfEvent(true, true, true, true);
                  std::cout<< "\n\n\n\n\n\n";
                  
                  for(unsigned int i=0; i < Ntp->NMCTaus(); i++){
                    for(unsigned int j=0; j < Ntp->NMCTauDecayProducts(i); j++){
                      std::cout<<"Product coming from index: "<<i<<" with product index: "<<j<<" : "<<Ntp->MCTauandProd_pdgid(i, j)<< std::endl;
                    }
                  }
                  */
                  Selection_Cut_Mu1_P.at(t).Fill(Mu1_LV.Vect().Mag(),1 );
                  Selection_Cut_Mu1_Eta.at(t).Fill(abs(Mu1_LV.Eta()),1 );
                  Selection_Cut_Mu2_P.at(t).Fill(Mu2_LV.Vect().Mag(),1 );
                  Selection_Cut_Mu2_Eta.at(t).Fill(abs(Mu2_LV.Eta()),1 );
                  Selection_Cut_Mu3_P.at(t).Fill(Mu3_LV.Vect().Mag(),1 );
                  Selection_Cut_Mu3_Eta.at(t).Fill(abs(Mu3_LV.Eta()),1 );
                  Selection_Cut_h_Pt.at(t).Fill(Tau_h_LV.Pt(),1 );
                  Selection_Cut_h_Eta.at(t).Fill(abs(Tau_h_LV.Eta()),1 );
                  
                  Selection_Cut_Mu1_p_eta_before.at(t).Fill(Mu1_LV.Vect().Mag(),abs(Mu1_LV.Eta()));
                  if(pass.at(Mu1_Candidate_p)&&pass.at(Mu1_Candidate_eta)){
                    Selection_Cut_Mu1_p_eta_after.at(t).Fill(Mu1_LV.Vect().Mag(),abs(Mu1_LV.Eta()));
                    if(!pass.at(Mu1_Candidate_recod)){
                      Selection_Cut_Mu1_p_eta_after_reco.at(t).Fill(Mu1_LV.Vect().Mag(),abs(Mu1_LV.Eta()));
                    }
                  }
                  
                  if(passAllUntil(Mu1_Candidate_eta)){
                  Selection_Cut_Mu2_p_eta_before.at(t).Fill(Mu2_LV.Vect().Mag(),abs(Mu2_LV.Eta()));
                  if(pass.at(Mu2_Candidate_p)&&pass.at(Mu2_Candidate_eta)){
                    Selection_Cut_Mu2_p_eta_after.at(t).Fill(Mu2_LV.Vect().Mag(),abs(Mu2_LV.Eta()));
                    if(!pass.at(Mu2_Candidate_recod)){
                      Selection_Cut_Mu2_p_eta_after_reco.at(t).Fill(Mu2_LV.Vect().Mag(),abs(Mu2_LV.Eta()));
                    }
                  }
                  }
                  
                  if(passAllUntil(Mu2_Candidate_eta)){
                  Selection_Cut_Mu3_p_eta_before.at(t).Fill(Mu3_LV.Vect().Mag(),abs(Mu3_LV.Eta()));
                  if(pass.at(Mu3_Candidate_p)&&pass.at(Mu3_Candidate_eta)){
                    Selection_Cut_Mu3_p_eta_after.at(t).Fill(Mu3_LV.Vect().Mag(),abs(Mu3_LV.Eta()));
                    if(!pass.at(Mu3_Candidate_recod)){
                      Selection_Cut_Mu3_p_eta_after_reco.at(t).Fill(Mu3_LV.Vect().Mag(),abs(Mu3_LV.Eta()));
                    }
                  }
                  }
                  
                  if(passAllUntil(Mu3_Candidate_eta)){
                  Selection_Cut_h_pt_eta_before.at(t).Fill(Tau_h_LV.Pt(),abs(Tau_h_LV.Eta()));
                  if(pass.at(Tau_h_Candidate_p)&&pass.at(Tau_h_Candidate_eta)){
                    Selection_Cut_h_pt_eta_after.at(t).Fill(Tau_h_LV.Pt(),abs(Tau_h_LV.Eta()));
                    if(!pass.at(Tau_h_Candidate_recod)){
                      Selection_Cut_h_pt_eta_after_reco.at(t).Fill(Tau_h_LV.Pt(),abs(Tau_h_LV.Eta()));
                    }
                  }
                  }
                  
                  Selection_Cut_Mu1_dR.at(t).Fill(dR1_max,1 );
                  Selection_Cut_Mu2_dR.at(t).Fill(dR2_max,1 );
                  Selection_Cut_Mu3_dR.at(t).Fill(dR3_max,1 );
                  Selection_Cut_h_dR.at(t).Fill(dR4_max,1 );
                  
                  Selection_Cut_Mu1_dR_large_scale.at(t).Fill(dR1_max,1 );
                  Selection_Cut_Mu2_dR_large_scale.at(t).Fill(dR2_max,1 );
                  Selection_Cut_Mu3_dR_large_scale.at(t).Fill(dR3_max,1 );
                  Selection_Cut_h_dR_large_scale.at(t).Fill(dR4_max,1 );
                  
          }
          
  }
  
  }//if(id!=1)
  //  std::cout << "Test 2. " << std::endl;
  if(!WhetherSignalMC){
    pass.at(WhetherDecayFound)=1;
    pass.at(Mu1_Candidate_p)=1;
    pass.at(Mu1_Candidate_eta)=1;
    pass.at(Mu2_Candidate_p)=1;
    pass.at(Mu2_Candidate_eta)=1;
    pass.at(Mu3_Candidate_p)=1;
    pass.at(Mu3_Candidate_eta)=1;
    pass.at(Tau_h_Candidate_p)=1;
    pass.at(Tau_h_Candidate_eta)=1;
    pass.at(Mu1_Candidate_recod)=1;
    pass.at(Mu2_Candidate_recod)=1;
    pass.at(Mu3_Candidate_recod)=1;
    pass.at(Tau_h_Candidate_recod)=1;
  }
  //  std::cout << "Test 3. " << std::endl;

  TLorentzVector Tau3MuLV(0,0,0,0);
  pass.at(SignalCandidate) = (value.at(SignalCandidate) >= cut.at(SignalCandidate));

  value.at(TripletPT)=0;
  if(signal_idx!=-1)
    {
      Tau3MuLV = Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0))+
        Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1))+
        Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2));
      
      value.at(TripletPT) = Tau3MuLV.Pt();
      Selection_Cut_3mu_Pt.at(t).Fill(value.at(TripletPT));
    }
  
  pass.at(TripletPT) = (value.at(TripletPT) >= cut.at(TripletPT));



  //  std::cout << "Test 4. " << std::endl;
  std::vector<int> Taus;
  std::vector<int> Taus_OppositeHemisphere_pT;
  std::vector<int> Taus_OppositeHemisphere_eta;
  std::vector<int> Taus_OppositeHemisphere_dR;
  std::vector<int> Taus_OppositeHemisphere;
  std::vector<int> Taus_OppositeHemisphere_OppositeCharge;
  std::vector<int> Taus_DeepTauJetsMedium;
  std::vector<int> Taus_DeepTauJetsTight;
  std::vector<int> Taus_DeepTauJetsLoose;



  value.at(nTaus_pT)  = -1;
  value.at(nTaus_eta)  = -1;
  value.at(nTaus_dR)  = 99.0;
  for(unsigned int itau=0; itau < Ntp->NTaus(); itau++)
    {
      if(signal_idx!=-1)
	{
	  //value.at(nTaus_pT)  = Ntp->Tau_P4(itau).Pt();
          //value.at(nTaus_eta)  = fabs(Ntp->Tau_P4(itau).Eta());
          //value.at(nTaus_dR)  = Ntp->Tau_P4(itau).DeltaR(Tau3MuLV);
          
          if(Ntp->Tau_P4(itau).Pt() > cut.at(nTaus_pT) ) Taus_OppositeHemisphere_pT.push_back(itau);
          if(fabs(Ntp->Tau_P4(itau).Eta()) < cut.at(nTaus_eta) ) Taus_OppositeHemisphere_eta.push_back(itau);
          if(Ntp->Tau_P4(itau).DeltaR(Tau3MuLV) > cut.at(nTaus_dR) ) Taus_OppositeHemisphere_dR.push_back(itau);
          
          if(Ntp->Tau_P4(itau).Pt() > cut.at(nTaus_pT) && fabs(Ntp->Tau_P4(itau).Eta()) < cut.at(nTaus_eta) && 
	     Ntp->Tau_P4(itau).DeltaR(Tau3MuLV) > cut.at(nTaus_dR) ) Taus_OppositeHemisphere.push_back(itau);
	       
    Selection_Cut_tauh_Pt.at(t).Fill(Ntp->Tau_P4(itau).Pt());
    Selection_Cut_tauh_Eta.at(t).Fill(fabs(Ntp->Tau_P4(itau).Eta()));
    Selection_Cut_tauh_DeltaR_3mu.at(t).Fill(Ntp->Tau_P4(itau).DeltaR(Tau3MuLV));
	  
      
	}
    }

  pass.at(nTaus_pT)  = ( Taus_OppositeHemisphere_pT.size() > 0 );
  
  pass.at(nTaus_eta) = ( Taus_OppositeHemisphere_eta.size() > 0 );
  
  pass.at(nTaus_dR)  = ( Taus_OppositeHemisphere_dR.size() > 0 );
  
  //  std::cout << "Test 5. " << std::endl;


  value.at(OSCharge)    =0;
  value.at(Tau3MuIsolation)   = -1;
  value.at(VisMass)           = -1;
  value.at(TriggerMatch) = 0;


  //  std::cout << "Test 6. " << std::endl;
  if(signal_idx!=-1)
    {
      int index_mu_1 = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);
      int index_mu_2 = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);
      int index_mu_3 = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);

      TLorentzVector TripletmuLV = Ntp->Muon_P4(index_mu_1) +  Ntp->Muon_P4(index_mu_2) +  Ntp->Muon_P4(index_mu_3);

      value.at(Tau3MuIsolation) = TripletmuLV.Pt()/  (TripletmuLV.Pt()  + Ntp->Muon_RelIso(index_mu_1) +
						                          Ntp->Muon_RelIso(index_mu_2) +
						                          Ntp->Muon_RelIso(index_mu_3) );
      
      Selection_Cut_3mu_Rel_Iso.at(t).Fill(value.at(Tau3MuIsolation));
    
    
      //Trigger Matching
      bool triggerCheck = 0.1;
      if(pass.at(HLT_TriggerOk))
        {
          vector<TLorentzVector> trigobjTriplet;
          for (int i=0; i<Ntp->NTriggerObjects(); i++)
            {
              TString name = Ntp->TriggerObject_name(i);
              //        if (!(name.Contains("tau3muDisplaced3muFltr"))) continue;
              TLorentzVector tmp;
              tmp.SetPtEtaPhiM(Ntp->TriggerObject_pt(i), Ntp->TriggerObject_eta(i), Ntp->TriggerObject_phi(i), PDG_Var::Muon_mass());
              trigobjTriplet.push_back(tmp);
            }
	  std::vector<TLorentzVector> muonTriplet;
          muonTriplet.push_back(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0)));
          muonTriplet.push_back(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1)));
          muonTriplet.push_back(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2)));

          if (trigobjTriplet.size()>=3) triggerCheck = Ntp->triggerMatchTriplet(muonTriplet, trigobjTriplet).first;
        }
      value.at(TriggerMatch) = triggerCheck;
    
      if(Taus_OppositeHemisphere.size()>0)
	{
	  for(auto i : Taus_OppositeHemisphere)
	    {
	      
	      int Tau3MuCharge = Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)) +
		Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)) +
		Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2));
	      if(Ntp->Tau_charge(i)*Tau3MuCharge == -1) Taus_OppositeHemisphere_OppositeCharge.push_back(i);
	    }
	  value.at(OSCharge) = Taus_OppositeHemisphere_OppositeCharge.size();
	}
    }
    
  //    std::cout << "Test 7. " << std::endl;
    pass.at(TriggerMatch) = (value.at(TriggerMatch)  ==  cut.at(TriggerMatch));    
    pass.at(OSCharge) = (value.at(OSCharge) >= cut.at(OSCharge));  
    


  for(auto i : Taus_OppositeHemisphere_OppositeCharge)
    {
      if(Ntp->Tau_byLooseDeepTau2017v2p1VSjet(i))  Taus_DeepTauJetsLoose.push_back(i);
      if(Ntp->Tau_byMediumDeepTau2017v2p1VSjet(i)) Taus_DeepTauJetsMedium.push_back(i);
      if(Ntp->Tau_byTightDeepTau2017v2p1VSjet(i))  Taus_DeepTauJetsTight.push_back(i);
      
      //      Logger(Logger::Info)<<"Tau:  " << i << " Loose/Med/Tight      DM/NewDM   "<<Ntp->Tau_byLooseDeepTau2017v2p1VSjet(i)  << "   "
      //			  << Ntp->Tau_byMediumDeepTau2017v2p1VSjet(i) <<"  "
      //			  << Ntp->Tau_byTightDeepTau2017v2p1VSjet(i)  << "           "
      //			  << Ntp->Tau_DecayModeFinding(i) << "   " 
	//			  << Ntp-> Tau_NewDecayModeFinding(i)<< "   DM=  "
      //			  << Ntp->Tau_DecayMode(i) <<std::endl;
    }


  value.at(DeepTauJets) = Taus_DeepTauJetsLoose.size(); // Loose
  //value.at(DeepTauJets) = Taus_DeepTauJetsMedium.size(); // Medium
  //value.at(DeepTauJets) = Taus_DeepTauJetsTight.size(); // Tight

  pass.at(DeepTauJets) = (value.at(DeepTauJets) >= cut.at(DeepTauJets));



  value.at(DeepTauMuons) = 0;
  value.at(DeepTauElectrons) = 0;
  std::vector<int> PassedDeepMuonsLoose;
  std::vector<int> PassedDeepMuonsMedium;
  std::vector<int> PassedDeepMuonsTight;

  for(auto i : Taus_DeepTauJetsMedium)
    {
      if(Ntp->Tau_byLooseDeepTau2017v2p1VSmu(i)) PassedDeepMuonsLoose.push_back(i);
      if(Ntp->Tau_byMediumDeepTau2017v2p1VSmu(i)) PassedDeepMuonsMedium.push_back(i);
      if(Ntp->Tau_byTightDeepTau2017v2p1VSmu(i)) PassedDeepMuonsTight.push_back(i);
      
    }
    
    value.at(DeepTauMuons) = PassedDeepMuonsLoose.size();      // Loose
  //  value.at(DeepTauMuons) = PassedDeepMuonsMedium.size(); // Medium
  //  value.at(DeepTauMuons) = PassedDeepMuonsTight.size();  // Tight


  std::vector<int> PassedDeepElectronsLoose;
  std::vector<int> PassedDeepElectronsMedium;
  std::vector<int> PassedDeepElectronsTight;
  for(auto i : PassedDeepMuonsLoose)
    {
      if(Ntp->Tau_byLooseDeepTau2017v2p1VSe(i)) PassedDeepElectronsLoose.push_back(i);
      if(Ntp->Tau_byMediumDeepTau2017v2p1VSe(i)) PassedDeepElectronsMedium.push_back(i);
      if(Ntp->Tau_byTightDeepTau2017v2p1VSe(i)) PassedDeepElectronsTight.push_back(i);

    }


    value.at(DeepTauElectrons) = PassedDeepElectronsLoose.size();       // Loose
  //  value.at(DeepTauElectrons) = PassedDeepElectronsMedium.size();  // Medium
  //  value.at(DeepTauElectrons) = PassedDeepElectronsTight.size();   // Tight


  pass.at(DeepTauMuons) = (value.at(DeepTauMuons) >= cut.at(DeepTauMuons));
  pass.at(DeepTauElectrons) = (value.at(DeepTauElectrons) >= cut.at(DeepTauElectrons));

  if(pass.at(DeepTauElectrons))
    {
      unsigned int tau_index = PassedDeepElectronsLoose.at(0);
      value.at(VisMass) = (Tau3MuLV + Ntp->Tau_P4(tau_index)).M();
      
      Selection_Cut_Vis_InvM.at(t).Fill(value.at(VisMass));
    }

  pass.at(Tau3MuIsolation) = (value.at(Tau3MuIsolation) > cut.at(Tau3MuIsolation));
  pass.at(VisMass)         = (value.at(VisMass) > 50 && value.at(VisMass) < 100);
  





  double wobs=1;
  double w;  
             
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  
  bool presel_cuts_tau;
  
  bool status=AnalysisCuts(t,w,wobs);

  if(status){ 

    //   std::cout << "Test 1. " << std::endl;
    
    unsigned int tau_h_idx = PassedDeepElectronsLoose.at(0);
    NumberOfTaus.at(t).Fill(Ntp->NTaus());

    TLorentzVector TauHLV = Ntp->Tau_P4(tau_h_idx);

    unsigned int muon_1_idx = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);
    unsigned int muon_2_idx = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);
    unsigned int muon_3_idx = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);

    ////////////////////////   sort muons by charge and dR and fill pair masses :
    /////
    vector<unsigned int> idx_vec;
    idx_vec.push_back(muon_1_idx);
    idx_vec.push_back(muon_2_idx);
    idx_vec.push_back(muon_3_idx);

    unsigned int os_mu_idx  = Ntp->SortedChargeMuons(idx_vec).at(0);
    unsigned int ss1_mu_idx = Ntp->SortedChargeMuons(idx_vec).at(1);
    unsigned int ss2_mu_idx = Ntp->SortedChargeMuons(idx_vec).at(2);

    TLorentzVector MuonOS  = Ntp->Muon_P4(os_mu_idx);
    TLorentzVector MuonSS1 = Ntp->Muon_P4(ss1_mu_idx);
    TLorentzVector MuonSS2 = Ntp->Muon_P4(ss2_mu_idx);



    if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){

      PairMass_OppositeSign_dR12.at(t).Fill((MuonOS+MuonSS2).M(),1 );
      PairMass_OppositeSign_dR13.at(t).Fill((MuonOS+MuonSS1).M(),1 );


    }else{

      PairMass_OppositeSign_dR12.at(t).Fill((MuonOS+MuonSS1).M(),1 );
      PairMass_OppositeSign_dR13.at(t).Fill((MuonOS+MuonSS2).M(),1 );

    }
    //////
    ///////////////////////////


    TLorentzVector Muon1LV = Ntp->Muon_P4(muon_1_idx);
    TLorentzVector Muon2LV = Ntp->Muon_P4(muon_2_idx);
    TLorentzVector Muon3LV = Ntp->Muon_P4(muon_3_idx);

    TLorentzVector Tau3muLV = Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)) + 
      Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)) + 
      Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2));

    TLorentzVector TauRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,0) +
      Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,1) +
      Ntp->Vertex_signal_KF_refittedTracksP4(signal_idx,2);



    LorentzVectorParticle Tau3MuLVP = Ntp->Tau3mu_LVP(  signal_idx );
    TVector3 Neutrino_Vect(Ntp->METEt()*TMath::Cos(Ntp->METPhi()),Ntp->METEt()*TMath::Sin(Ntp->METPhi()),Ntp->METEt()/TMath::Tan(TauHLV.Theta()));
    TLorentzVector Neutrino_LV(Neutrino_Vect,Neutrino_Vect.Mag());



    float RelativeIsolationMu1 = Ntp->Muon_RelIso(muon_1_idx);
    float RelativeIsolationMu2 = Ntp->Muon_RelIso(muon_2_idx);
    float RelativeIsolationMu3 = Ntp->Muon_RelIso(muon_3_idx);

    Tau3MuRelativeIsolation.at(t).Fill(    Tau3muLV.Pt()/(RelativeIsolationMu1 + RelativeIsolationMu2 + RelativeIsolationMu3 + Tau3muLV.Pt()),1);
    TauHDecayMode.at(t).Fill(Ntp->Tau_DecayMode(tau_h_idx), 1);
    VisibleDiTauMass.at(t).Fill((TauHLV + Tau3muLV).M(), 1);
    MTT.at(t).Fill( (Tau3muLV + TauHLV  + Neutrino_LV).M(), 1);

    bool PlotMCOnly(false);  // and blind for data
    if(id!=1) PlotMCOnly = true;
    if(id==1 && ( (TauRefitLV.M() > 1.1 && TauRefitLV.M() < 1.7233333) or (TauRefitLV.M() > 1.8333333 && TauRefitLV.M() < 2.2)) ) PlotMCOnly=true;


    if(PlotMCOnly)    TripletMass.at(t).Fill(TauRefitLV.M(),1);

    //////////// kinematics 
    TripletPt.at(t).Fill(Tau3muLV.Pt(),1);
    OppositeTauPt.at(t).Fill(TauHLV.Pt(),1);
    //////////// kinematics 

    //    std::cout << "Test 1.5 " << std::endl;



    ////////////////////////////////////////
    ///
    TLorentzVector OppositeSideLV = TauHLV;
    TLorentzVector Neutrino_LV_Guess_Result;//Guessed Neutrino LV
    //    if(id == WhetherSignalMC)
    if(false)
      {
	matched_pdgId.at(t).Fill(Ntp->matchTruth(OppositeSideLV),1);
	matched_dR.at(t).Fill(Ntp->matchToTruthTauDecay(OppositeSideLV).DeltaR(OppositeSideLV),1);

        TLorentzVector MCMuon1LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)));
        TLorentzVector MCMuon2LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)));
        TLorentzVector MCMuon3LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2)));
        TLorentzVector MCTauLV  = MCMuon1LV+ MCMuon2LV + MCMuon3LV;

        Muon1DRToTruth.at(t).Fill(Muon1LV.DeltaR(MCMuon1LV),1);
        Muon2DRToTruth.at(t).Fill(Muon2LV.DeltaR(MCMuon2LV),1);
        Muon3DRToTruth.at(t).Fill(Muon3LV.DeltaR(MCMuon3LV),1);

	dR_betweenTruth_VisibleTaus.at(t).Fill(MCTauLV.DeltaR(Ntp->matchToTruthTauDecay(OppositeSideLV)),1);
        
        //////// Getting the Neutrino MC
        int HadronicT_MCTau_idx(0);//Getting the index of Spectator Tau in MC
        double Max_dR(99.0);
        for(int i=1; i<Ntp->NMCTaus(); i++){
                if(TauHLV.DeltaR(Ntp->MCTau_p4(i))<Max_dR){// Make HadronicT_MCTau_idx the hadronic tau 
                        HadronicT_MCTau_idx = i;
                }
                HadronicT_MCTau_idx = i;
        }
        TLorentzVector TauMC=Ntp->MCTauandProd_p4(HadronicT_MCTau_idx,0);
        TLorentzVector TauNeutrinoMC(0,0,0,0);
        bool FoundMCNeutrino(false);
        for(int i=1; i<Ntp->NMCTauDecayProducts(HadronicT_MCTau_idx); i++){
                if(abs(Ntp->MCTauandProd_pdgid(HadronicT_MCTau_idx,i))==16){
                        TauNeutrinoMC = Ntp->MCTauandProd_p4(HadronicT_MCTau_idx,i);
                        FoundMCNeutrino = true;
                }
        }
        
        if(FoundMCNeutrino){// Neutrinos satisfying tau invariant mass can be created with any direction of neutrino, PROVING THERE ARE MULTIPLE SOLUTIONS (without taking Z into account). We can construct a neutrino LV consistent with tau and Z mass. The same thing will happen with data, so we can't use this method to distinguish between data and MC? (Or, is it?)
                
                //std::cout << "Found neutrino: " << std::endl;
                double Max_Diff(999.0);
                int Division_No_Thet(1);//200
                int Division_No_Phi(1);//1000
                
                double Phi_Visible_Tau = TauHLV.Phi();
                double Theta_Visible_Tau = TauHLV.Theta();
                TVector3 Visible_Tau_Vec = TauHLV.Vect();
                
                for(int Theta_1_idx=0; Theta_1_idx<Division_No_Thet;Theta_1_idx++){
                        for(int Phi_1_idx=0; Phi_1_idx<Division_No_Phi;Phi_1_idx++){
                                
                                double Theta_1 = 0.0 + Theta_1_idx*(0.5)/Division_No_Thet;
                                double Phi_1 = 0.0 + Phi_1_idx*(2*TMath::Pi())/Division_No_Phi;
                                
                                TVector3 Neutrino_Vect_Guess(TMath::Sin(Theta_1)*TMath::Cos(Phi_1),TMath::Sin(Theta_1)*TMath::Sin(Phi_1),TMath::Cos(Theta_1));
                                if(Phi_Visible_Tau >= TMath::Pi()) Phi_Visible_Tau = Phi_Visible_Tau-2*TMath::Pi();
                                if(Phi_Visible_Tau < -TMath::Pi()) Phi_Visible_Tau = Phi_Visible_Tau+2*TMath::Pi();
                                Neutrino_Vect_Guess.RotateY(Theta_Visible_Tau);
                                Neutrino_Vect_Guess.RotateZ(Phi_Visible_Tau);
                                
                                double scaling_x = (1.77686*1.77686 - TauHLV.M()*TauHLV.M())/(2*(TauHLV.E()-TauHLV.Vect().Mag()*TMath::Cos(Neutrino_Vect_Guess.Angle(Visible_Tau_Vec))));//Imposing tau mass criteria
                                TLorentzVector Neutrino_LV_Guess(scaling_x*Neutrino_Vect_Guess,scaling_x);
                                
                                // Boosting to COM frame
                                //TLorentzVector TauHLV_Boosted_Guess = TauHLV;
                                //TauHLV_Boosted_Guess.Boost(-1*(TauHLV+Neutrino_LV_Guess).BoostVector());
                                //TLorentzVector Neutrino_LV_Boosted_Guess = Neutrino_LV_Guess;
                                //Neutrino_LV_Boosted_Guess.Boost(-1*(TauHLV+Neutrino_LV_Guess).BoostVector());
                                
                                double Diff_Val = abs((TauHLV+Neutrino_LV_Guess+Tau3muLV).M() - 91.1876);// Imposing Z mass
                                if(Diff_Val < Max_Diff){
                                        Max_Diff = Diff_Val;
                                        Neutrino_LV_Guess_Result = Neutrino_LV_Guess;
                                }
        
                        }
                }
                
                //std::cout << "Angle from guess result to actual: " << Neutrino_LV_Guess_Result.Vect().Angle(TauNeutrinoMC.Vect()) << " with Max_Diff: " << Max_Diff << std::endl;
                
                dR_betweenTruth_NeutrinoGuess.at(t).Fill(Neutrino_LV_Guess_Result.DeltaR(TauNeutrinoMC),1);
                dR_betweenTruth_Tau.at(t).Fill((Neutrino_LV_Guess_Result+TauHLV).DeltaR(TauMC),1);
                Z_Pt.at(t).Fill((TauMC+Tau3muLV).Pt(),1);
                
                TVector3 MC_MET = (TauMC-TauNeutrinoMC+Tau3muLV).Vect();
                MC_MET.SetZ(0.);
                MC_MET = -1 * MC_MET;
                
                //std::cout << "MC Delta Phi from actual to MET: " << MC_MET.Phi()-TauNeutrinoMC.Vect().Phi() << std::endl;
                //std::cout << "Z angle: " << (TauMC+Tau3muLV).Vect().Theta() << " Z (Pt): "<< (TauMC+Tau3muLV).Pt() << " Tau vis to Tau3mu angle: " << (TauMC-TauNeutrinoMC).Vect().Phi() - Tau3muLV.Vect().Phi() << " Tau to Tau3mu angle: " << (TauMC).Vect().Phi() - Tau3muLV.Vect().Phi() << std::endl;
                
        }//FoundMCNeutrino
        
        
        


      }//if(id != 1)


        //*** fill up the T3MMiniTree.root for statistical analysis
        
        m3m = TauRefitLV.M();
        
        dataMCtype = id;
        event_weight =1; // 1 for data
        if(dataMCtype == 1){event_weight =1;}
        else if(dataMCtype == 210233){event_weight =5.00e-04;} // event_weight is a value Lumi Scale ; Lumi 45710; Evno 5625
        else if(dataMCtype == 210232){event_weight =4.93e-04;} // Lumi 45710; Evno 5659
        else if(dataMCtype == 210231){event_weight =4.87e-04;} // Lumi 45710; Evno 21314
        
        if(MuonOS.DeltaR(MuonSS1) > MuonOS.DeltaR(MuonSS2)){
          mDr1 = (MuonOS+MuonSS2).M();
          mDr2 = (MuonOS+MuonSS1).M();
        }else{
          mDr1 = (MuonOS+MuonSS1).M();
          mDr2 = (MuonOS+MuonSS2).M();
        }
        
        
        m12 = (MuonOS+MuonSS1).M();
        m13 = (MuonOS+MuonSS2).M();
        LumiScale = 1.;
        T3MMiniTree->Fill();
  
	//  std::cout << "Test 2. " << std::endl;
  }
}


void  ZTau3MuTauh_PreFC::Finish(){

  //*** write down the T3MMiniTree.root for statistical analysis
  T3MFMiniTree = new TFile("T3MMiniTree.root","recreate");
  T3MMiniTree->SetDirectory(T3MFMiniTree);
  T3MFMiniTree->Write();
  T3MFMiniTree->Close();
  
  Selection::Finish();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
}





