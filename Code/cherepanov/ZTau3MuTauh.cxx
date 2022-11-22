#include "ZTau3MuTauh.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
//#include "Logger.h"

ZTau3MuTauh::ZTau3MuTauh(TString Name_, TString id_):
  Selection(Name_,id_)
{
  // This is a class constructor;
}

ZTau3MuTauh::~ZTau3MuTauh(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  ZTau3MuTauh::Configure(){

  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
    if(i==TripletKinematics)  cut.at(TripletKinematics)=1;
    if(i==DeepTauJets)        cut.at(DeepTauJets)=1;
    if(i==DeepTauMuons)       cut.at(DeepTauMuons)=1;
    if(i==DeepTauElectrons)   cut.at(DeepTauElectrons)=1;
    if(i==OppositeSide)       cut.at(OppositeSide)=1;
    if(i==OSCharge)           cut.at(OSCharge)=1;
    if(i==nTaus)              cut.at(nTaus)=1;

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
    else if(i==nTaus){
      title.at(i)=" At least one $\\tau_{h}$, $pT>20 GeV, |\\eta| < 2.4$ ";
      hlabel="number of taus";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nTaus_",htitle,4,-0.5,3.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nTaus_",htitle,4,-0.5,3.5,hlabel,"Events"));
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

    else if(i==OppositeSide){
      title.at(i)="At least one  $\\tau_{h}$ is on opposite side $|\\Delta R| > 1/2$";
      hlabel="delta phi ";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OppositeSide_",htitle,4,-0.5,3.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OppositeSide_",htitle,4,-0.5,3.5,hlabel,"Events"));
    }

    else if(i==TripletKinematics){
      title.at(i)="3,3,2 GeV,  $|\\eta| < 2.4$";
      htitle=title.at(i);
      hlabel="3mu kinematics";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TripletKinematics_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TripletKinematics_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }

    else if(i==SignalCandidate){
      title.at(i)="At least one $\\tau_{3\\mu}$ candidate";
      htitle=title.at(i);
      hlabel="N $3\\mu$ candidates";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,5,-0.5,4.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }
    else if(i==OSCharge){
      title.at(i)="Charge $\\tau_{h}$ * $\\tau_{3\\mu}$ =  ";
      title.at(i)+=cut.at(OSCharge);
      title.at(i)+=" (at least one)";
      htitle=title.at(i);
      hlabel="Opposite charge? ";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OSCharge_",htitle,4,-0.5,3.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OSCharge_",htitle,4,-0.5,3.5,hlabel,"Events"));
    }



  }
  // Setup NPassed Histogams

  NumberOfTaus=HConfig.GetTH1D(Name+"_NumberOfTaus","NumberOfTaus",5,-0.5,4.5,"Number of #tau ","Events");



  Tau3MuRelativeIsolation=HConfig.GetTH1D(Name+"_Tau3MuRelativeIsolation","Tau3MuRelativeIsolation",50,0.,1.,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T})","#Delta R < 0.4 ");
  TauHDecayMode=HConfig.GetTH1D(Name+"_TauHDecayMode","TauHDecayMode",12,-0.5,11.5,"HPS #tau_{h} decay mode","Events");
  VisibleDiTauMass=HConfig.GetTH1D(Name+"_VisibleDiTauMass","VisibleDiTauMass",50,20,100,"M_{#tau(h) - #tau(3#mu)}, GeV (visible mass)","Events");
  MTT=HConfig.GetTH1D(Name+"_MTT","MTT",50,60,120,"M_{#tau(h) - #tau(3#mu)}, GeV (collinear approximation)","Events");
  TripletMass=HConfig.GetTH1D(Name+"_TripletMass","TripletMass",30,1.4,2.2,"M_{3#mu}, GeV","Events");



  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms



  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  ZTau3MuTauh::Store_ExtraDist(){ 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output

  Extradist1d.push_back(&NumberOfTaus);
  Extradist1d.push_back(&Tau3MuRelativeIsolation);
  Extradist1d.push_back(&TauHDecayMode);
  Extradist1d.push_back(&VisibleDiTauMass);
  Extradist1d.push_back(&MTT);
  Extradist1d.push_back(&TripletMass);



}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  ZTau3MuTauh::doEvent(){ 

  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection



  bool HLTOk(false);
  bool L1Ok(false);

  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLTName = Ntp->HLTName(iTrigger);
    //std::cout<<"HLT:   "  << Ntp->HLTName(iTrigger)  << "  fires  "<< Ntp->HLTDecision(iTrigger)<< std::endl;
    if(HLTName.Contains("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_v") && Ntp->HLTDecision(iTrigger) ) { HLTOk = true;}
  }
  
  
  value.at(TriggerOk)=HLTOk;
  pass.at(TriggerOk)=(value.at(TriggerOk) == cut.at(TriggerOk)); 



  value.at(SignalCandidate) = Ntp->NThreeMuons();

  unsigned int  signal_idx=0;
  double min_chi2(99.);

  for(unsigned int i_idx =0; i_idx < Ntp->NThreeMuons(); i_idx++){
    if(Ntp->Vertex_Signal_KF_Chi2(i_idx) < min_chi2){
      min_chi2 = Ntp->Vertex_Signal_KF_Chi2(i_idx);
      signal_idx = i_idx;
    }
  }

  TLorentzVector Tau3MuLV(0,0,0,0);
  pass.at(SignalCandidate) = (value.at(SignalCandidate) >= cut.at(SignalCandidate));

  value.at(TripletKinematics) = 0;
  if( pass.at(SignalCandidate) )
    {
      Tau3MuLV = Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0))+
	Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1))+
	Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2));
      
      if(  (Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0)).Pt() >=3 && fabs(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0)).Eta())  < 2.4)  &&
	   (Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1)).Pt() >=3 && fabs(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1)).Eta())  < 2.4)  &&
	   (Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2)).Pt() >=3 && fabs(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2)).Eta())  < 2.4) )
	value.at(TripletKinematics) = 1;

    }

  pass.at(TripletKinematics) = (value.at(TripletKinematics) == cut.at(TripletKinematics));


  std::vector<int> Taus;
  std::vector<int> Taus_OppositeHemisphere;
  std::vector<int> Taus_OppositeHemisphere_OppositeCharge;
  std::vector<int> Taus_DeepTauJetsMedium;
  std::vector<int> Taus_DeepTauJetsTight;
  std::vector<int> Taus_DeepTauJetsLoose;




  for(unsigned int itau=0; itau < Ntp->NTaus(); itau++)
    {
      if(Ntp->Tau_P4(itau).Pt() > 20 && fabs(Ntp->Tau_P4(itau).Eta()) < 2.3) Taus.push_back(itau);
    }

  value.at(nTaus)  = Taus.size();
  pass.at(nTaus)  = ( value.at(nTaus) >= cut.at(nTaus) );
  

  value.at(OppositeSide)=0; 
  value.at(OSCharge)    =0;
  if(pass.at(TripletKinematics))
    {
      for(auto i : Taus)
	{
	  if(  Ntp->Tau_P4(i).DeltaR(Tau3MuLV) > 0.5   ) Taus_OppositeHemisphere.push_back(i);
	  //	  if(fabs(Ntp->DeltaPhi(Ntp->Tau_P4(i).Phi(), Tau3MuLV.Phi() ))  > TMath::Pi() / 2    ) Taus_OppositeHemisphere.push_back(i);

	}
      value.at(OppositeSide)= Taus_OppositeHemisphere.size();

      for(auto i : Taus_OppositeHemisphere)
	{

	  int Tau3MuCharge = Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)) +
	    Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)) +
	    Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2));

	  if(Ntp->Tau_charge(i)*Tau3MuCharge == -1) Taus_OppositeHemisphere_OppositeCharge.push_back(i);
	}
      value.at(OSCharge) = Taus_OppositeHemisphere_OppositeCharge.size();

    }

    pass.at(OppositeSide) = (value.at(OSCharge) >= cut.at(OSCharge));  
    pass.at(OSCharge) = (value.at(OSCharge) >= cut.at(OSCharge));  




  for(auto i : Taus_OppositeHemisphere_OppositeCharge)
    {
      if(Ntp->Tau_byLooseDeepTau2017v2p1VSjet(i))  Taus_DeepTauJetsLoose.push_back(i);
      if(Ntp->Tau_byMediumDeepTau2017v2p1VSjet(i)) Taus_DeepTauJetsMedium.push_back(i);
      if(Ntp->Tau_byTightDeepTau2017v2p1VSjet(i))  Taus_DeepTauJetsTight.push_back(i);
      
      Logger(Logger::Info)<<"Tau:  " << i << " Loose/Med/Tight      DM/NewDM   "<<Ntp->Tau_byLooseDeepTau2017v2p1VSjet(i)  << "   "
			  << Ntp->Tau_byMediumDeepTau2017v2p1VSjet(i) <<"  "
			  << Ntp->Tau_byTightDeepTau2017v2p1VSjet(i)  << "           "
			  << Ntp->Tau_DecayModeFinding(i) << "   " 
	//			  << Ntp-> Tau_NewDecayModeFinding(i)<< "   DM=  "
			  << Ntp->Tau_DecayMode(i) <<std::endl;
    }


  value.at(DeepTauJets) = Taus_DeepTauJetsMedium.size(); // Loose
  //  value.at(DeepTauJets) = Taus_DeepTauJetsMedium.size();// Medium
  //  value.at(DeepTauJets) = Taus_DeepTauJetsMedium.size();// Tight

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

    double wobs=1;
    double w;  
             
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  

  
  bool status=AnalysisCuts(t,w,wobs);

  if(status){ 

    unsigned int tau_h_idx = PassedDeepElectronsLoose.at(0);

    NumberOfTaus.at(t).Fill(Ntp->NTaus());

    TLorentzVector TauHLV = Ntp->Tau_P4(tau_h_idx);

    unsigned int muon_1_idx = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);
    unsigned int muon_2_idx = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);
    unsigned int muon_3_idx = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);


    TLorentzVector Tau3muLV = Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)) + 
      Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)) + 
      Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2));

    int Tau3MuCharge = Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)) +
      Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)) +
      Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2));





    LorentzVectorParticle Tau3MuLVP = Ntp->Tau3mu_LVP(  signal_idx );
    TVector3 Neutrino_Vect(Ntp->METEt()*TMath::Cos(Ntp->METPhi()),Ntp->METEt()*TMath::Sin(Ntp->METPhi()),Ntp->METEt()/TMath::Tan(TauHLV.Theta()));
    TLorentzVector Neutrino_LV(Neutrino_Vect,Neutrino_Vect.Mag());



    float RelativeIsolationMu1 = Ntp->Muon_RelIso(muon_1_idx);
    float RelativeIsolationMu2 = Ntp->Muon_RelIso(muon_2_idx);
    float RelativeIsolationMu3 = Ntp->Muon_RelIso(muon_3_idx);

    std::cout<<"  rel iso  "<<RelativeIsolationMu1<< "   "
	     <<RelativeIsolationMu2 << "    "
	     <<RelativeIsolationMu3 <<std::endl;

    Tau3MuRelativeIsolation.at(t).Fill(    Tau3muLV.Pt()/(RelativeIsolationMu1 + RelativeIsolationMu2 + RelativeIsolationMu3 + Tau3muLV.Pt()),1);
    TauHDecayMode.at(t).Fill(Ntp->Tau_DecayMode(tau_h_idx), 1);
    VisibleDiTauMass.at(t).Fill((TauHLV + Tau3muLV).M(), 1);
    MTT.at(t).Fill( (Tau3muLV + TauHLV  + Neutrino_LV).M(), 1);

    TripletMass.at(t).Fill(Tau3muLV.M(),1);

  }
}


void  ZTau3MuTauh::Finish(){

  Selection::Finish();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
}





