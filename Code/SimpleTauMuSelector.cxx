#include "SimpleTauMuSelector.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

SimpleTauMuSelector::SimpleTauMuSelector(TString Name_, TString id_):
  Selection(Name_,id_)
{
  // This is a class constructor;
}

SimpleTauMuSelector::~SimpleTauMuSelector(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  SimpleTauMuSelector::Configure(){

  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)        cut.at(TriggerOk)=1;
    if(i==nMuons)           cut.at(nMuons)=4;
    if(i==SignalCandidate)  cut.at(SignalCandidate)=1;
    if(i==OSCharge)         cut.at(OSCharge)=0;
    if(i==pTCut1)           cut.at(pTCut1)=15;
    if(i==pTCut2)           cut.at(pTCut2)=8.5;
    if(i==pTCut3)           cut.at(pTCut3)=0;
    if(i==pTCutTauMu)       cut.at(pTCutTauMu)=5;

  }

  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    // book here the N-1 and N-0 histrogramms for each cut
    if(i==TriggerOk){
      title.at(i)="Trigger ";
      hlabel="Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==nMuons){
      title.at(i)="No. of $\\mu$";
      hlabel="No. of $\\mu$ ";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nMuons_",htitle,8,-0.5,7.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nMuons_",htitle,8,-0.5,7.5,hlabel,"Events"));
    }
    else if(i==SignalCandidate){
      title.at(i)="SignalCandidate ";
      hlabel="SignalCandidate ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==OSCharge){
      title.at(i)="OSCharge ";
      hlabel="OSCharge ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OSCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OSCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==pTCut1){
      title.at(i)="$pT(\\mu_{1}) > 15 GeV$";
      hlabel="$pT(\\mu_{1})$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_pTCut1_",htitle,25,5,50,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_pTCut1_",htitle,25,5,50,hlabel,"Events"));
    }
    else if(i==pTCut2){
      title.at(i)="$pT(\\mu_{2}) > 8.5 GeV$";
      hlabel="$pT(\\mu_{2})$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_pTCut2_",htitle,25,5,25,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_pTCut2_",htitle,25,5,25,hlabel,"Events"));
    }
    else if(i==pTCut3){
      title.at(i)="$pT(\\mu_{3}) > 0 GeV$";
      hlabel="$pT(\\mu_{3})$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_pTCut3_",htitle,25,5,25,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_pTCut3_",htitle,25,5,25,hlabel,"Events"));
    }
    else if(i==pTCutTauMu){
      title.at(i)="pT($\\tau_{\\mu}$) $> 5$ GeV";
      htitle=title.at(i);
      hlabel="pT$(\\tau_{\\mu})$, GeV";      
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_pTCutTauMu_",htitle,25,19,60,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_pTCutTauMu_",htitle,25,19,60,hlabel,"Events"));
    }

  }
  // Setup NPassed Histogams

  NumberOfMuons=HConfig.GetTH1D(Name+"_NumberOfMuons","NumberOfMuons",5,-0.5,4.5,"Number of #mu ","Events");
  Mu_TauCand_Inv_Mass_1=HConfig.GetTH1D(Name+"_Mu_TauCand_Inv_Mass_1","Mu_TauCand_Inv_Mass_1",100,0.0,160.0,"#tau_{#mu} + #tau (3#mu) Inv Mass, GeV","Events");
  Mu_TauCand_Inv_Mass_2=HConfig.GetTH1D(Name+"_Mu_TauCand_Inv_Mass_2","Mu_TauCand_Inv_Mass_2",100,0.0,160.0,"#tau_{#mu} + #tau (3#mu) Inv Mass, GeV","Events");
  Mu_TauCand_Inv_Mass_3=HConfig.GetTH1D(Name+"_Mu_TauCand_Inv_Mass_3","Mu_TauCand_Inv_Mass_3",100,0.0,160.0,"#tau_{#mu} + #nu_{#mu} + #tau (3#mu) Inv Mass, GeV","Events");
  MET_Et=HConfig.GetTH1D(Name+"_MET_Et","MET_Et",100,0.0,100.0,"MET Et, GeV","Events");
  MET_Phi=HConfig.GetTH1D(Name+"_MET_Phi","MET_Phi",20,-3.2,3.2,"MET #Phi ","Events");
  Kinematics=HConfig.GetTH1D(Name+"_Kinematics","Kinematics",20,0,3.2,"#tau_{#mu} , #tau (3#mu) opening angle","Events");
  Kinematics_1=HConfig.GetTH1D(Name+"_Kinematics_1","Kinematics_1",40,-3.2,3.2,"#tau_{#mu} , #tau (3#mu) #Delta#Phi","Events");
  Kinematics_MissingTrMass=HConfig.GetTH1D(Name+"_Kinematics_MissingTrMass","Kinematics_MissingTrMass",100,0,500.,"MissingTrMass ","Events");
  NumTausvsNumMuons=HConfig.GetTH2D(Name+"_NumTausvsNumMuons","NumTausvsNumMuons",5,-0.5,4.5,5,-0.5,4.5,"N 3#mu ","N #mu");
  Kinematics_TauXPtEta=HConfig.GetTH2D(Name+"_Kinematics_TauXPtEta","Kinematics_TauXPtEta",20,0.,20.,20,0,2.5,"pT, GeV","eta ");
  Tau3muMass=HConfig.GetTH1D(Name+"_Tau3muMass","Tau3muMass",50,1.4,2.,"M_{3#mu}, GeV","Events");
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms



  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  SimpleTauMuSelector::Store_ExtraDist(){ 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output

  Extradist1d.push_back(&NumberOfMuons);
  Extradist1d.push_back(&Mu_TauCand_Inv_Mass_1);
  Extradist1d.push_back(&Mu_TauCand_Inv_Mass_2);
  Extradist1d.push_back(&Mu_TauCand_Inv_Mass_3);
  Extradist1d.push_back(&MET_Et);
  Extradist1d.push_back(&MET_Phi);
  Extradist1d.push_back(&Kinematics);
  Extradist1d.push_back(&Kinematics_1);
  Extradist1d.push_back(&Kinematics_MissingTrMass);


  Extradist1d.push_back(&Tau3muMass);


  Extradist2d.push_back(&NumTausvsNumMuons);
  Extradist2d.push_back(&Kinematics_TauXPtEta);


}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  SimpleTauMuSelector::doEvent(){ 

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
  
  //random_num = rndm.Rndm();
  random_num = 1.0;
  
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

  value.at(TriggerOk)=(HLTOk && L1Ok);
  pass.at(TriggerOk)=(value.at(TriggerOk)==cut.at(TriggerOk));
  
  value.at(SignalCandidate)=0;
  unsigned int  signal_idx=0;
  
  double min_chi2(99.);
  for(unsigned int i_idx =0; i_idx < Ntp->NThreeMuons(); i_idx++){
    if(Ntp->Vertex_Signal_KF_Chi2(i_idx) < min_chi2){
      min_chi2 = Ntp->Vertex_Signal_KF_Chi2(i_idx);
      signal_idx = i_idx;
    }
  }
  
  TLorentzVector Tau_3mu_LV;
  int  cand_charge=0;
  unsigned int mu1_idx;
  unsigned int mu2_idx;
  unsigned int mu3_idx;

  
  if(Ntp->NThreeMuons()>0){
    value.at(SignalCandidate) = Ntp->NThreeMuons();
    
    mu1_idx = Ntp->ThreeMuonIndices(signal_idx).at(0);
    mu2_idx = Ntp->ThreeMuonIndices(signal_idx).at(1);
    mu3_idx = Ntp->ThreeMuonIndices(signal_idx).at(2);
    
    Tau_3mu_LV = Ntp->Muon_P4(mu1_idx)+Ntp->Muon_P4(mu2_idx)+Ntp->Muon_P4(mu3_idx);
    cand_charge  = Ntp->Muon_charge(mu1_idx)+Ntp->Muon_charge(mu2_idx)+Ntp->Muon_charge(mu3_idx);
  
  }
  
  pass.at(SignalCandidate) = (value.at(SignalCandidate) >= cut.at(SignalCandidate));
  
  value.at(nMuons) = Ntp->NMuons();
  pass.at(nMuons)  = ( value.at(nMuons) >= cut.at(nMuons) );
  
  std::vector<int> muons_idx;
  std::vector<int> muons_os_idx;


  for(unsigned int idx =0; idx < Ntp->NMuons(); idx++){
    //    if(Ntp->Tau_DecayModeFinding(idx))
    muons_idx.push_back(idx);
  }
  
  TLorentzVector MuonFromTauLV_PhiCand;
  int  MuonFromTau_PhiCand_idx=0;
  double max_dPhi(0.);
  
  for(unsigned int i = 0; i < muons_idx.size() && Ntp->NThreeMuons()>0; i++){
  
    //Making sure muon is not from the 3mu candidates
    int muon_matched(0);
    int num = muons_idx.at(i);
    TLorentzVector TMuonLV_i = Ntp->Muon_P4(num);
    for(unsigned int i =0; i < 3; i++){
      
      if(TMuonLV_i.DeltaR(Ntp->Muon_P4(i)) < 0.0001 ){
        muon_matched+=1;
      }
    }
    
    bool Test1 = Ntp->Muon_isGlobalMuon(num);
    bool Test2 = Ntp->Muon_ID(num);
    bool Test3 = abs(Ntp->DeltaPhi(TMuonLV_i.Phi(), Tau_3mu_LV.Phi()))>2.5;
    
    if(Ntp->Muon_charge(muons_idx.at(i)) * cand_charge == -1 && muon_matched==0 && Test1 ){
    
      muons_os_idx.push_back(muons_idx.at(i));
    
      double muTau_tau3mu_dPhi = abs(Ntp->DeltaPhi(TMuonLV_i.Phi(), Tau_3mu_LV.Phi()));
      if(muTau_tau3mu_dPhi > max_dPhi){
        max_dPhi = muTau_tau3mu_dPhi;
        MuonFromTau_PhiCand_idx = muons_idx.at(i);
        MuonFromTauLV_PhiCand=TMuonLV_i;
      }
      
    }

  }
  
  //std::cout<<"   Muon size   "<< muons_idx.size() << "   OS Muon size   " << muons_os_idx.size() << " 3mu Candidates: "<< Ntp->NThreeMuons() <<std::endl;
  
  //  std::cout<<"  dm_os_size   "<<taus_dm_os_idx.size() <<std::endl;
  
  
  value.at(OSCharge) = muons_os_idx.size();
  pass.at(OSCharge) = (value.at(OSCharge) > cut.at(OSCharge));
  
  
  value.at(pTCut1) = 0.0;
  value.at(pTCut2) = 0.0;
  value.at(pTCut3) = 0.0;
  value.at(pTCutTauMu) = 0.0;
  
  if(muons_os_idx.size()>0 && Ntp->NThreeMuons()>0){
    
    value.at(pTCut1) = Ntp->Muon_P4(mu1_idx).Pt();
    value.at(pTCut2) = Ntp->Muon_P4(mu2_idx).Pt();
    value.at(pTCut3) = Ntp->Muon_P4(mu3_idx).Pt();
    value.at(pTCutTauMu) = MuonFromTauLV_PhiCand.Pt();
  }
  pass.at(pTCut1) = (value.at(pTCut1) > cut.at(pTCut1));
  pass.at(pTCut2) = (value.at(pTCut2) > cut.at(pTCut2));
  pass.at(pTCut3) = (value.at(pTCut3) > cut.at(pTCut3));
  pass.at(pTCutTauMu) = (value.at(pTCutTauMu) > cut.at(pTCutTauMu));
  
  double wobs=1;
  double w;  
             
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  
  bool status=AnalysisCuts(t,w,wobs);
  
  NumTausvsNumMuons.at(t).Fill(Ntp->NMuons(), Ntp->NThreeMuons());
  
  
  if(status){
    
    //    std::cout << "Next Event: " <<std::endl;
    
    Tau3muMass.at(t).Fill(Tau_3mu_LV.M(),1);

    NumberOfMuons.at(t).Fill(Ntp->NMuons());
    
    
    //for(unsigned int i =0; i < taus_dm_idx.size(); i++){
      //int idx = taus_dm_idx.at(i);
      int idx = MuonFromTau_PhiCand_idx;
      TLorentzVector MuFromTau_idx = Ntp->Muon_P4(idx);
      
      //if(std::count(taus_dm_os_idx.begin(), taus_dm_os_idx.end(), idx) == 1){
      
        Mu_TauCand_Inv_Mass_1.at(t).Fill( (MuFromTau_idx+Tau_3mu_LV).M() );
        
        if((MuFromTau_idx+Tau_3mu_LV).M()<5.0 && id!=1){
          int pdg_idx = Ntp->getMatchTruthIndex(MuFromTau_idx);
          std::cout << "pdgid: "<< pdg_idx <<std::endl;
        }
        
        Kinematics_TauXPtEta.at(t).Fill(MuFromTau_idx.Pt(),MuFromTau_idx.Eta());
        
        //if((MuFromTau_idx+Tau_3mu_LV).M()>10.0){
        //}
        /*
        for(unsigned int i =0; i < 3; i++){
          if((MuFromTau_idx+Tau_3mu_LV).M()<10.0){
            std::cout << "Tau dR when Inv Mass < 10.0 is: " << MuFromTau_idx.DeltaR(Ntp->Muon_P4(i)) <<std::endl;
          }
          
          if((MuFromTau_idx+Tau_3mu_LV).M()>=10.0){
            std::cout << "Tau dR when Inv Mass > 10.0 is: " << MuFromTau_idx.DeltaR(Ntp->Muon_P4(i)) <<std::endl;
          }
        }
        */
        
        // Kinematics
        
        //Making Bmeson point towards Z axis (and making phi = 0 for the tau: commented): fully reco
        double Phi_init_reco = Tau_3mu_LV.Phi();
        double Theta_init_reco = Tau_3mu_LV.Theta();
        TLorentzVector Tau_3mu_LV_mod = Tau_3mu_LV;
        TLorentzVector MuFromTau_idx_mod = MuFromTau_idx;
        if(Phi_init_reco >= TMath::Pi()) Phi_init_reco = Phi_init_reco-2*TMath::Pi();
        if(Phi_init_reco <=-TMath::Pi()) Phi_init_reco = Phi_init_reco+2*TMath::Pi();
        MuFromTau_idx_mod.RotateZ(-Phi_init_reco);
        MuFromTau_idx_mod.RotateY(-Theta_init_reco);
        Tau_3mu_LV_mod.RotateZ(-Phi_init_reco);
        Tau_3mu_LV_mod.RotateY(-Theta_init_reco);
        
        TVector3 p1 = Tau_3mu_LV.Vect();
        TVector3 p2 = MuFromTau_idx.Vect();
        p1.SetZ(0.);
        p2.SetZ(0.);
        
        Kinematics.at(t).Fill(MuFromTau_idx_mod.Theta());
        Kinematics_1.at(t).Fill(p1.DeltaPhi(p2));
        
        // Missing transverse mass
        Kinematics_MissingTrMass.at(t).Fill(2*Ntp->METEt()*TMath::Sqrt(1.776*1.776+Tau_3mu_LV.Pz()*Tau_3mu_LV.Pz())*(1-TMath::Cos(Ntp->METPhi()-Tau_3mu_LV.Phi()))); //use definition transverse mass for 2 particles
        
        //Approximate neutrino LV
        
        TVector3 Neutrino_Vect(Ntp->METEt()*TMath::Cos(Ntp->METPhi()),Ntp->METEt()*TMath::Sin(Ntp->METPhi()),Ntp->METEt()/TMath::Tan(MuFromTau_idx.Theta()));
        TLorentzVector Neutrino_LV(Neutrino_Vect,Neutrino_Vect.Mag());
        
        //if((MuFromTau_idx+Tau_3mu_LV).M()>10.0){
        Mu_TauCand_Inv_Mass_3.at(t).Fill( (MuFromTau_idx+Tau_3mu_LV+Neutrino_LV).M() );
        //}
        
      //}
    //}// end i
    
    Mu_TauCand_Inv_Mass_2.at(t).Fill( (MuonFromTauLV_PhiCand+Tau_3mu_LV).M() );
    
    MET_Et.at(t).Fill( Ntp->METEt() );
    
    MET_Phi.at(t).Fill( Ntp->METPhi() );
    
  }
}


void  SimpleTauMuSelector::Finish(){

  Selection::Finish();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
}





