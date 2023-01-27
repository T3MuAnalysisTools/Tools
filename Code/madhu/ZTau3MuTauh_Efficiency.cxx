#include "ZTau3MuTauh_Efficiency.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
//#include "Logger.h"

ZTau3MuTauh_Efficiency::ZTau3MuTauh_Efficiency(TString Name_, TString id_):
  Selection(Name_,id_)
{
  // This is a class constructor;
}

ZTau3MuTauh_Efficiency::~ZTau3MuTauh_Efficiency(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  ZTau3MuTauh_Efficiency::Configure(){

  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    
    if(i==L1_TriggerOk)             cut.at(L1_TriggerOk)=1;
    if(i==HLT_TriggerOk)            cut.at(HLT_TriggerOk)=1;
    if(i==Mu1_Candidate)            cut.at(Mu1_Candidate)=1;
    if(i==Mu2_Candidate)            cut.at(Mu2_Candidate)=1;
    if(i==Mu3_Candidate)            cut.at(Mu3_Candidate)=1;
    if(i==Tau_h_Candidate)          cut.at(Tau_h_Candidate)=1;
    if(i==Mu1_Candidate_recod)      cut.at(Mu1_Candidate_recod)=1;
    if(i==Mu2_Candidate_recod)      cut.at(Mu2_Candidate_recod)=1;
    if(i==Mu3_Candidate_recod)      cut.at(Mu3_Candidate_recod)=1;
    if(i==Tau_h_Candidate_recod)    cut.at(Tau_h_Candidate_recod)=1;


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
    else if(i==Mu1_Candidate){
      title.at(i)=" Whether GEN level $\\mu_{1}$ has $p>2.5 GeV, |\\eta| < 2.4$ ";
      hlabel="If GEN level $\\mu_{1}$ has $p>2 GeV, |\\eta| < 2.4$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nTaus_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nTaus_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Mu2_Candidate){
      title.at(i)=" Whether GEN level $\\mu_{2}$ has $p>2.5 GeV, |\\eta| < 2.4$ ";
      hlabel="If GEN level $\\mu_{2}$ has $p>2 GeV, |\\eta| < 2.4$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nTaus_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nTaus_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Mu3_Candidate){
      title.at(i)=" Whether GEN level $\\mu_{3}$ has $p>2.5 GeV, |\\eta| < 2.4$ ";
      hlabel="If GEN level $\\mu_{3}$ has $p>2 GeV, |\\eta| < 2.4$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nTaus_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nTaus_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Tau_h_Candidate){
      title.at(i)=" Whether GEN level $\\tau_{h}$ has $pT>18 GeV, |\\eta| < 2.4$ ";
      hlabel="If GEN level $\\tau_{h}$ has $p>2 GeV, |\\eta| < 2.4$";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nTaus_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nTaus_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Mu1_Candidate_recod){
      title.at(i)=" Whether $\\mu_{1}$ is reconstructed in MC ";
      hlabel="If $\\mu_{1}$ is reconstructed in MC";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nTaus_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nTaus_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Mu2_Candidate_recod){
      title.at(i)=" Whether $\\mu_{2}$ is reconstructed in MC ";
      hlabel="If $\\mu_{2}$ is reconstructed in MC";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nTaus_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nTaus_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Mu3_Candidate_recod){
      title.at(i)=" Whether $\\mu_{3}$ is reconstructed in MC ";
      hlabel="If $\\mu_{3}$ is reconstructed in MC";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nTaus_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nTaus_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Tau_h_Candidate_recod){
      title.at(i)=" Whether $\\tau_{h}$ is reconstructed in MC ";
      hlabel="If $\\tau_{h}$ is reconstructed in MC";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nTaus_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nTaus_",htitle,2,-0.5,1.5,hlabel,"Events"));
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
  dR_betweenTruth_VisibleTaus=HConfig.GetTH1D(Name+"_dR_betweenTruth_VisibleTaus","dR_betweenTruth_VisibleTaus",20,0,0.02,"#Delta R Truth #tau's prods ","Events");
  
  Correlation_MCnu_MET=HConfig.GetTH2D(Name+"_Correlation_MCnu_MET","MET vs MC #nu Energy",100,0.0,100.0,80,0,80.0,"MC #nu Energy, GeV","MET, GeV");
  
  Correlation_Angles=HConfig.GetTH2D(Name+"_Correlation_Angles","MET vs MC #nu Angle",100,0,0.5,100,0,0.5,"MC #nu Angle","Another Angle");
  
  dR_MC_Neutrino_To_Visible_Tau=HConfig.GetTH1D(Name+"_dR_MC_Neutrino_To_Visible_Tau","dR_MC_Neutrino_To_Visible_Tau",50,0,1.0,"#Delta R MC #nu to visible #tau","Events");
  Inv_Mass_VisPrdt_Neutrino=HConfig.GetTH1D(Name+"_Inv_Mass_VisPrdt_Neutrino","Inv_Mass_VisPrdt_Neutrino",60,0,20.0,"Inv mass of #nu + visible #tau","Events");
  
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms



  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  ZTau3MuTauh_Efficiency::Store_ExtraDist(){ 
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
  
  
  
  Extradist1d.push_back(&dR_MC_Neutrino_To_Visible_Tau);
  Extradist1d.push_back(&Inv_Mass_VisPrdt_Neutrino);
  
  Extradist2d.push_back(&Correlation_MCnu_MET);
  Extradist2d.push_back(&Correlation_Angles);



}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  ZTau3MuTauh_Efficiency::doEvent(){ 

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
  
  
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLTName = Ntp->HLTName(iTrigger);
    //std::cout<<"HLT:   "  << Ntp->HLTName(iTrigger)  << "  fires  "<< Ntp->HLTDecision(iTrigger)<< std::endl;
    if(HLTName.Contains("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_v") && Ntp->HLTDecision(iTrigger) ) { HLTOk = true;}
  }
  
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
  
  if (DoubleMu0Fired || DoubleMu4Fired) {DoubleMuFired = true;}
  if (DoubleMuFired || TripleMuFired) L1Ok = true;

  value.at(L1_TriggerOk)=(L1Ok);
  pass.at(L1_TriggerOk)=(value.at(L1_TriggerOk)==cut.at(L1_TriggerOk));
  
  value.at(HLT_TriggerOk)=(HLTOk);
  pass.at(HLT_TriggerOk)=(value.at(HLT_TriggerOk)==cut.at(HLT_TriggerOk));



  unsigned int  signal_idx=0;
  double min_chi2(99.);

  for(unsigned int i_idx =0; i_idx < Ntp->NThreeMuons(); i_idx++){
    if(Ntp->Vertex_Signal_KF_Chi2(i_idx) < min_chi2){
      min_chi2 = Ntp->Vertex_Signal_KF_Chi2(i_idx);
      signal_idx = i_idx;
    }
  }

  TLorentzVector Tau3MuLV(0,0,0,0);
  
  //for(int i = 0; i < Ntp->NMCSignalParticles(); i++){
  //  std::cout << " pdgid of particle" << i+1 << ": " << Ntp->MCSignalParticle_pdgid(i) << std::endl;
  //}
  
  //for(int i = 0; i < Ntp->NMCTauDecayProducts(0); i++){
  //  std::cout << " pdgid of particle" << i+1 << ": " << Ntp->MCTauandProd_pdgid(0, i) << std::endl;
  //}
  //for(int i = 0; i < Ntp->NMCTauDecayProducts(1); i++){
  //  std::cout << " pdgid of particle" << i+1 << ": " << Ntp->MCTauandProd_pdgid(1, i) << std::endl;
  //}
  
  
  int Whether_decay_found(-1);
  int idx_for_3mu(-1);
  
  if(Ntp->NMCTaus()==2&&Ntp->NMCTauDecayProducts(0)>3){
    if(abs(Ntp->MCTauandProd_pdgid(0, 0))==15&&abs(Ntp->MCTauandProd_pdgid(0, 1))==13&&abs(Ntp->MCTauandProd_pdgid(0, 2))==13&&abs(Ntp->MCTauandProd_pdgid(0, 3))==13){
      Whether_decay_found=1;
      idx_for_3mu=0;
    }
  }
  if(Ntp->NMCTaus()==2&&Ntp->NMCTauDecayProducts(1)>3){
    if(abs(Ntp->MCTauandProd_pdgid(1, 0))==15&&abs(Ntp->MCTauandProd_pdgid(1, 1))==13&&abs(Ntp->MCTauandProd_pdgid(1, 2))==13&&abs(Ntp->MCTauandProd_pdgid(1, 3))==13){
      Whether_decay_found=1;
      idx_for_3mu=1;
    }
  }
  
  
  
  TLorentzVector Mu1_LV;
  TLorentzVector Mu2_LV;
  TLorentzVector Mu3_LV;
  TLorentzVector Tau_nu_LV;
  TLorentzVector Tau_h_LV;
  if(Whether_decay_found==1){
    Mu1_LV=Ntp->MCTauandProd_p4(idx_for_3mu, 1);
    Mu2_LV=Ntp->MCTauandProd_p4(idx_for_3mu, 2);
    Mu3_LV=Ntp->MCTauandProd_p4(idx_for_3mu, 3);
    
    for(int i = 0; i < Ntp->NMCTauDecayProducts(1-idx_for_3mu); i++){
      if(abs(Ntp->MCTauandProd_pdgid(1-idx_for_3mu, i))==16){
        Tau_nu_LV=Ntp->MCTauandProd_p4(1-idx_for_3mu, i);
      }
    }
    Tau_h_LV = Ntp->MCTauandProd_p4(1-idx_for_3mu, 0) - Tau_nu_LV;
  }
  
  
  value.at(Mu1_Candidate)=(Mu1_LV.Vect().Mag()>2.5)&&(Mu1_LV.Eta()<2.4);
  pass.at(Mu1_Candidate)=value.at(Mu1_Candidate);
  value.at(Mu2_Candidate)=(Mu2_LV.Vect().Mag()>2.5)&&(Mu2_LV.Eta()<2.4);
  pass.at(Mu2_Candidate)=value.at(Mu2_Candidate);
  value.at(Mu3_Candidate)=(Mu3_LV.Vect().Mag()>2.5)&&(Mu3_LV.Eta()<2.4);
  pass.at(Mu3_Candidate)=value.at(Mu3_Candidate);
  
  value.at(Tau_h_Candidate)=(Tau_h_LV.Pt()>18.0)&&(Tau_h_LV.Eta()<2.4);
  pass.at(Tau_h_Candidate)=value.at(Tau_h_Candidate);
  
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

    }
  
  for(unsigned int itau=0; itau < Ntp->NTaus(); itau++)
    {
      if(Ntp->Tau_P4(itau).DeltaR(Tau_h_LV)<dR4_max){
        dR4_max=Ntp->Tau_P4(itau).DeltaR(Tau_h_LV);
      }
    }
  
  value.at(Mu1_Candidate_recod)=(dR1_max<0.05);
  pass.at(Mu1_Candidate_recod)=value.at(Mu1_Candidate_recod);
  value.at(Mu2_Candidate_recod)=(dR2_max<0.05);
  pass.at(Mu2_Candidate_recod)=value.at(Mu2_Candidate_recod);
  value.at(Mu3_Candidate_recod)=(dR3_max<0.05);
  pass.at(Mu3_Candidate_recod)=value.at(Mu3_Candidate_recod);
  
  value.at(Tau_h_Candidate_recod)=(dR4_max<0.05);
  pass.at(Tau_h_Candidate_recod)=value.at(Tau_h_Candidate_recod);
  
  
  //std::cout << " idx: " << idx_for_3mu << std::endl;
  //std::cout << " pdgid of particle2: " << Ntp->MCTauandProd_pdgid(1, 1) << std::endl;
  //std::cout << " pdgid of particle3: " << Ntp->MCTauandProd_pdgid(1, 2) << std::endl;
  




    double wobs=1;
    double w;  
             
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  

  
  bool status=AnalysisCuts(t,w,wobs);

  if(status){
    
  }
  
  
}


void  ZTau3MuTauh_Efficiency::Finish(){

  Selection::Finish();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
}





