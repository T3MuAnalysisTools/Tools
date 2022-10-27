#include "SimpleTauSelector.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

SimpleTauSelector::SimpleTauSelector(TString Name_, TString id_):
  Selection(Name_,id_)
{
  // This is a class constructor;
}

SimpleTauSelector::~SimpleTauSelector(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  SimpleTauSelector::Configure(){

  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)        cut.at(TriggerOk)=1;
    if(i==PrimeVtx)         cut.at(PrimeVtx)=5; // Here for example we place cut value on number of PVs
    if(i==nTaus)            cut.at(nTaus)=1;
    if(i==SignalCandidate)  cut.at(SignalCandidate)=1;
    if(i==OSCharge)        cut.at(OSCharge)=0;

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
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
    }
    else if(i==TriggerOk){
      title.at(i)="Trigger ";
      hlabel="Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==nTaus){
      title.at(i)=" nTaus ";
      hlabel="number of taus";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nTaus_",htitle,3,-0.5,2.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nTaus_",htitle,3,-0.5,2.5,hlabel,"Events"));
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

  }
  // Setup NPassed Histogams

  NumberOfTaus=HConfig.GetTH1D(Name+"_NumberOfTaus","NumberOfTaus",5,-0.5,4.5,"Number of #tau ","Events");
  TauX_TauCand_Inv_Mass_1=HConfig.GetTH1D(Name+"_TauX_TauCand_Inv_Mass_1","TauX_TauCand_Inv_Mass_1",100,0.0,100.0,"#tau (X) + #tau (3#mu) Inv Mass","Events");
  TauX_TauCand_Inv_Mass_2=HConfig.GetTH1D(Name+"_TauX_TauCand_Inv_Mass_2","TauX_TauCand_Inv_Mass_2",100,0.0,100.0,"#tau (X) + #tau (3#mu) Inv Mass","Events");
  DecayModes=HConfig.GetTH1D(Name+"_DecayModes","DecayModes",16,-0.5,15.5,"Decay Modes ","Events");
  MET_Et=HConfig.GetTH1D(Name+"_MET_Et","MET_Et",100,0.0,100.0,"MET Et ","Events");
  MET_Phi=HConfig.GetTH1D(Name+"_MET_Phi","MET_Phi",20,-3.2,3.2,"MET Phi ","Events");
  Kinematics=HConfig.GetTH1D(Name+"_Kinematics","Kinematics",20,0,3.2,"#tau (X) angle to #tau (3#mu) ","Events");
  Kinematics_1=HConfig.GetTH1D(Name+"_Kinematics_1","Kinematics_1",20,-3.2,3.2,"#tau (X) angle to #tau (3#mu) in tranverse plane ","Events");
  NumTausvsNumMuons=HConfig.GetTH2D(Name+"_NumTausvsNumMuons","NumTausvsNumMuons",5,-0.5,4.5,5,-0.5,4.5,"N 3#mu ","N #tau");
  Kinematics_TauXPtEta=HConfig.GetTH2D(Name+"_Kinematics_TauXPtEta","Kinematics_TauXPtEta",100,0.,100.,20,0,4.0,"pT ","eta ");

  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms



  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  SimpleTauSelector::Store_ExtraDist(){ 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output

  Extradist1d.push_back(&NumberOfTaus);
  Extradist1d.push_back(&TauX_TauCand_Inv_Mass_1);
  Extradist1d.push_back(&TauX_TauCand_Inv_Mass_2);
  Extradist1d.push_back(&DecayModes);
  Extradist1d.push_back(&MET_Et);
  Extradist1d.push_back(&MET_Phi);
  Extradist1d.push_back(&Kinematics);
  Extradist1d.push_back(&Kinematics_1);
  Extradist2d.push_back(&NumTausvsNumMuons);
  Extradist2d.push_back(&Kinematics_TauXPtEta);


}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  SimpleTauSelector::doEvent(){ 

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

  value.at(PrimeVtx)=Ntp->NVtx(); 
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx)); 
  
  
  
  value.at(SignalCandidate)=0;
  unsigned int  signal_idx=0;
  
  double min_chi2(99.);
  for(unsigned int i_idx =0; i_idx < Ntp->NThreeMuons(); i_idx++){
    if(Ntp->Vertex_Signal_KF_Chi2(i_idx) < min_chi2){
      min_chi2 = Ntp->Vertex_Signal_KF_Chi2(i_idx);
      signal_idx = i_idx;
    }
  }
  
  TLorentzVector TauLV;
  int  cand_charge=0;
  
  if(Ntp->NThreeMuons()>0){
  
    value.at(SignalCandidate) = Ntp->NThreeMuons();

    unsigned int mu1_idx = Ntp->ThreeMuonIndices(signal_idx).at(0); 
    unsigned int mu2_idx = Ntp->ThreeMuonIndices(signal_idx).at(1); 
    unsigned int mu3_idx = Ntp->ThreeMuonIndices(signal_idx).at(2);
    
    TauLV = Ntp->Muon_P4(mu1_idx)+Ntp->Muon_P4(mu2_idx)+Ntp->Muon_P4(mu3_idx);
    cand_charge  = Ntp->Muon_charge(mu1_idx)+Ntp->Muon_charge(mu2_idx)+Ntp->Muon_charge(mu3_idx);
  
  }
  
  pass.at(SignalCandidate) = (value.at(SignalCandidate) >= cut.at(SignalCandidate));
  
  TLorentzVector TauXLV;
  
  value.at(nTaus) = Ntp->NTaus();
  pass.at(nTaus)  = ( value.at(nTaus) >= cut.at(nTaus) );
  unsigned int  TauX_idx=0;
  value.at(OSCharge) = 10;
  
  double min_dR(99.);
  for(unsigned int idx =0; idx < Ntp->NTaus(); idx++){
    
    TauXLV = Ntp->Tau_P4(idx);
    
    double taux_tau3mu_dR = TauLV.DeltaR(TauXLV);
    if(taux_tau3mu_dR < min_dR){
      min_dR = taux_tau3mu_dR;
      TauX_idx = idx;
    }
    
    if((cand_charge+Ntp->Tau_charge(idx))==0){
      value.at(OSCharge) = 0;
    }
  }
  
  
  
  pass.at(OSCharge) = (value.at(OSCharge) == cut.at(OSCharge));
  
  double wobs=1;
  double w;  
             
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  
  
  
  bool status=AnalysisCuts(t,w,wobs);
  
  NumTausvsNumMuons.at(t).Fill(Ntp->NTaus(), Ntp->NThreeMuons());
  
  
  
  
  if(status){ 
  
    NumberOfTaus.at(t).Fill(Ntp->NTaus());
    
    for(unsigned int idx =0; idx < Ntp->NTaus(); idx++){
      TLorentzVector TauXLV_idx = Ntp->Tau_P4(idx);
      if((Ntp->Tau_charge(idx)+cand_charge) != 0){
      
        TauX_TauCand_Inv_Mass_1.at(t).Fill( (TauXLV_idx+TauLV).M() );
        
        DecayModes.at(t).Fill( Ntp->Tau_DecayMode(idx) );
        
        Kinematics_TauXPtEta.at(t).Fill(TauXLV_idx.Pt(),TauXLV_idx.Eta());
        
        // Kinematics
        
        //Making Bmeson point towards Z axis (and making phi = 0 for the tau: commented): fully reco
        double Phi_init_reco = TauLV.Phi();
        double Theta_init_reco = TauLV.Theta();
        TLorentzVector TauLV_mod = TauLV;
        TLorentzVector TauXLV_idx_mod = TauXLV_idx;
        if(Phi_init_reco >= TMath::Pi()) Phi_init_reco = Phi_init_reco-2*TMath::Pi();
        if(Phi_init_reco <=-TMath::Pi()) Phi_init_reco = Phi_init_reco+2*TMath::Pi();
        TauXLV_idx_mod.RotateZ(-Phi_init_reco);
        TauXLV_idx_mod.RotateY(-Theta_init_reco);
        TauLV_mod.RotateZ(-Phi_init_reco);
        TauLV_mod.RotateY(-Theta_init_reco);
        
        TVector3 p1 = TauLV.Vect();
        TVector3 p2 = TauXLV_idx.Vect();
        p1.SetZ(0.);
        p2.SetZ(0.);
        
        Kinematics.at(t).Fill(TauXLV_idx_mod.Theta());
        Kinematics_1.at(t).Fill(p1.DeltaPhi(p2));
        
        // Missing transverse mass
        
      }
    }
    
    TauX_TauCand_Inv_Mass_2.at(t).Fill( (TauXLV+TauLV).M() );
    
    MET_Et.at(t).Fill( Ntp->METEt() );
    
    MET_Phi.at(t).Fill( Ntp->METPhi() );
    
    

  }
}


void  SimpleTauSelector::Finish(){

  Selection::Finish();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
}





