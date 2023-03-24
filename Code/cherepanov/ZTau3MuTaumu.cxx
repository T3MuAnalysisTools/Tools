#include "ZTau3MuTaumu.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
//#include "Logger.h"

ZTau3MuTaumu::ZTau3MuTaumu(TString Name_, TString id_):
  Selection(Name_,id_),
  AnalysisName(Name_)
{
    
  // This is a class constructor;
}

ZTau3MuTaumu::~ZTau3MuTaumu(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  ZTau3MuTaumu::Configure(){



  TString treeprefix;
  if(Ntp->GetInputNtuplePath().Contains("z2tautau")) treeprefix="z2tautau";
  if(Ntp->GetInputNtuplePath().Contains("DoubleMuonLowMass")) treeprefix="DoubleMuonLowMass";




  T3MMiniTree= new TTree(treeprefix + "_" + AnalysisName,"Mini Tree Input for mva");
  T3MMiniTree->Branch("m3m",&mini_m3m);
  T3MMiniTree->Branch("dataMCtype",&mini_dataMCtype);
  T3MMiniTree->Branch("event_weight",&mini_event_weight);
  T3MMiniTree->Branch("tripletpt",&mini_TripletPt);
  T3MMiniTree->Branch("muonpt",&mini_MuonPt);
  T3MMiniTree->Branch("tripletisolation",&mini_Tau3MuRelativeIsolation);
  T3MMiniTree->Branch("oppositemuonisolation",&mini_OppositeMuRelativeIsolation);
  T3MMiniTree->Branch("ditaumass",&mini_ditaumass);


  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==L1_TriggerOk)       cut.at(L1_TriggerOk)=1;
    if(i==HLT_TriggerOk)      cut.at(HLT_TriggerOk)=1;
    if(i==SignalCandidate)    cut.at(SignalCandidate)=1;
    if(i==TripletPT)          cut.at(TripletPT)=20;
    if(i==OSCharge)           cut.at(OSCharge)=1;
    if(i==nMuons)             cut.at(nMuons)=1;
    if(i==Tau3MuIsolation)    cut.at(Tau3MuIsolation)=0.5;
    if(i==TriggerMatch)       cut.at(TriggerMatch)=1;
    if(i==MuonIsolation)      cut.at(MuonIsolation)=0.5;
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
    else if(i==nMuons){
      //      title.at(i)=" At least one extra loose(PF+Gl/Tr) $\\mu$, $pT>15 GeV, |\\eta| < 2.4$ ";
      title.at(i)=" At least one extra (PF+GL) $\\mu$, $pT>10 GeV, |\\eta| < 2.4$ and $\\Delta R (\\mu-3\\mu) >$ 1/2";
      hlabel="At least one extra muon (loose)";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nMuons_",htitle,4,-0.5,3.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nMuons_",htitle,4,-0.5,3.5,hlabel,"Events"));
    }

    else if(i==TripletPT){
      title.at(i)="pT(3$\\mu$)  $>$  30 ";
      htitle=title.at(i);
      hlabel="pT(#tau_{3#mu}) , GeV ";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TripletPT_",htitle,50,5,80,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TripletPT_",htitle,50,5,80,hlabel,"Events"));
    }

    else if(i==SignalCandidate){
      title.at(i)="At least one $\\tau_{3\\mu}$ candidate (3,3,3 GeV,  $|\\eta| < 2.4$) ";
      htitle=title.at(i);
      hlabel="N $3\\mu$ candidates";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,5,-0.5,4.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }
    else if(i==OSCharge){
      title.at(i)="Charge $\\mu$ * $\\tau_{3\\mu}$ = -1; ";
      title.at(i)+=" (at least one)";
      htitle=title.at(i);
      hlabel="Opposite charge ";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OSCharge_",htitle,4,-0.5,3.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OSCharge_",htitle,4,-0.5,3.5,hlabel,"Events"));
    }


    else if(i==Tau3MuIsolation){
      title.at(i)="$3\\mu $ Relative isolation  $ > $ 0.8 ";
      //      title.at(i)+= cut.at(Tau3MuIsolation);
      htitle=title.at(i);
      hlabel="I(#tau_{3#mu})= p_{T}(#tau_{3#mu})/(p_{T}(#tau_{3#mu}) + #sum p_{T})";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau3MuIsolation_",htitle,50,0,1.1,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau3MuIsolation_",htitle,50,0,1.1,hlabel,"Events"));
    }


    else if(i==MuonIsolation){
      title.at(i)="Opposite $\\mu $ Relative isolation  $ > $ 0.85";
      //      title.at(i)+= cut.at(MuonIsolation);
      htitle=title.at(i);
      hlabel="I(#mu)= p_{T}(#mu)/(p_{T}(#mu) + #sum p_{T})";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonIsolation_",htitle,50,0,1.1,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonIsolation_",htitle,50,0,1.1,hlabel,"Events"));
    }

    else if(i==VisMass){
      title.at(i)="40 GeV $< M(\\tau(\\mu) - \\tau(3\\mu))  < $ 85 GeV";
      htitle=title.at(i);
      hlabel="M_{#tau(#mu) - #tau(3#mu)}, GeV (visible mass)";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_VisMass_",htitle,70,0,150,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_VisMass_",htitle,70,0,150,hlabel,"Events"));
    }

    else if(i==TriggerMatch){
      title.at(i)="Selected 3$\\mu$ matched to trg";
      hlabel="Trigger Matched ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }

  }
  // Setup NPassed Histogams




  Tau3MuRelativeIsolation=HConfig.GetTH1D(Name+"_Tau3MuRelativeIsolation","Tau3MuRelativeIsolation",50,0.,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T})","#Delta R < 0.4 ");
  OppositeMuRelativeIsolation=HConfig.GetTH1D(Name+"_OppositeMuRelativeIsolation","OppositeMuRelativeIsolation",50,0.,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T})","#Delta R < 0.4 ");
  
  VisibleDiTauMass=HConfig.GetTH1D(Name+"_VisibleDiTauMass","VisibleDiTauMass",70,0.,150,"M_{#tau(#mu) - #tau(3#mu)}, GeV (visible mass)","Events");
  MTT=HConfig.GetTH1D(Name+"_MTT","MTT",70,0.,140,"M_{#tau(#mu) - #tau(3#mu)}, GeV (collinear approximation)","Events");
  TripletMass=HConfig.GetTH1D(Name+"_TripletMass","TripletMass",30,1.1,2.2,"M_{3#mu}, GeV","Events");
  TripletMassWR=HConfig.GetTH1D(Name+"_TripletMassWR","TripletMassWR",70,1,8,"M_{3#mu}, GeV","Events");


  PairMass_OppositeSign_dR12=HConfig.GetTH1D(Name+"_PairMass_OppositeSign_dR12","PairMass_OppositeSign_dR12",40,0.2,2.,"M_{1}, GeV (OS - SS dR sorted)","Events");
  PairMass_OppositeSign_dR13=HConfig.GetTH1D(Name+"_PairMass_OppositeSign_dR13","PairMass_OppositeSign_dR13",40,0.2,2.,"M_{2}, GeV (OS - SS dR sorted)","Events");




  matched_pdgId=HConfig.GetTH1D(Name+"_matched_pdgId","matched_pdgId",25,-0.5,24.5,"pdgID MC matched","Events");
  matched_dR=HConfig.GetTH1D(Name+"_matched_dR","matched_dR",50,-0.1,0.5,"#Delta R(MC-RECO) Object opposite to #tau_{3#mu}","Events");


  Muon1DRToTruth=HConfig.GetTH1D(Name+"_Muon1DRToTruth","Muon1DRToTruth",20,0,0.02,"#Delta R (reco - mc) #mu_{1} ","Events");
  Muon2DRToTruth=HConfig.GetTH1D(Name+"_Muon2DRToTruth","Muon2DRToTruth",20,0,0.02,"#Delta R (reco - mc) #mu_{2} ","Events");
  Muon3DRToTruth=HConfig.GetTH1D(Name+"_Muon3DRToTruth","Muon3DRToTruth",20,0,0.02,"#Delta R (reco - mc) #mu_{3} ","Events");
  dR_betweenTruth_VisibleTaus=HConfig.GetTH1D(Name+"_dR_betweenTruth_VisibleTaus","dR_betweenTruth_VisibleTaus",20,0,5,"#Delta R Truth #tau's prods ","Events");

  TripletPt=HConfig.GetTH1D(Name+"_TripletPt","TripletPt",50,2,80,"pT(3#mu), GeV ","Events");
  OppositeMuonPt=HConfig.GetTH1D(Name+"_OppositeMuonPt","OppositeMuonPt",50,2,40,"pT(#mu), GeV ","Events");

  TripletEta=HConfig.GetTH1D(Name+"_TripletEta","TripletEta",50,-2.5,2.5,"#eta(3#mu)","Events");
  OppositeMuonEta=HConfig.GetTH1D(Name+"_OppositeMuonEta","OppositeMuonEta",50,-2.5,2.5,"#eta(#mu)","Events");





  SingleVsThreeMuTrigger=HConfig.GetTH2D(Name+"_SingleVsThreeMuTrigger","SingleVsThreeMuTrigger",2,-0.5,1.5,2,-0.5,1.5,"HLT_Tau3Mu ... ","HLT_IsoMu27(24)");


  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms



  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  ZTau3MuTaumu::Store_ExtraDist(){ 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output


  Extradist1d.push_back(&Tau3MuRelativeIsolation);
  Extradist1d.push_back(&OppositeMuRelativeIsolation);
  Extradist1d.push_back(&VisibleDiTauMass);
  Extradist1d.push_back(&MTT);
  Extradist1d.push_back(&TripletMass);
  Extradist1d.push_back(&TripletMassWR);
  Extradist1d.push_back(&matched_pdgId);
  Extradist1d.push_back(&matched_dR);

  Extradist1d.push_back(&Muon1DRToTruth);
  Extradist1d.push_back(&Muon2DRToTruth);
  Extradist1d.push_back(&Muon3DRToTruth);
  Extradist1d.push_back(&dR_betweenTruth_VisibleTaus);

  Extradist1d.push_back(&PairMass_OppositeSign_dR12);
  Extradist1d.push_back(&PairMass_OppositeSign_dR13);

  Extradist1d.push_back(&TripletPt);
  Extradist1d.push_back(&OppositeMuonPt);

  Extradist1d.push_back(&TripletEta);
  Extradist1d.push_back(&OppositeMuonEta);
  Extradist2d.push_back(&SingleVsThreeMuTrigger);

}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  ZTau3MuTaumu::doEvent(){ 

  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection

  //  std::cout<<"   id   "<< id << std::endl;

  bool HLTOk(false);
  bool L1Ok(false);
  bool DoubleMu0Fired(false);
  bool DoubleMu4Fired(false);
  bool DoubleMuFired(false);
  bool TripleMuFired(false);
  bool randomFailed(false);
  
  
  bool Tau3muHLT(false);
  bool TauIsoMu(false);

  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){

    TString HLTName = Ntp->HLTName(iTrigger);
    //    if(Ntp->HLTDecision(iTrigger))    std::cout<<"HLT:   "  << Ntp->HLTName(iTrigger)  << "  fires  "<< Ntp->HLTDecision(iTrigger)<< std::endl;
    //    if(HLTName.Contains("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_v") && Ntp->HLTDecision(iTrigger) ) { HLTOk = true;}


    if(HLTName.Contains("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_v") && Ntp->HLTDecision(iTrigger) ) { Tau3muHLT = true;}
    if(  ( HLTName.Contains("HLT_IsoMu27_v1") || HLTName.Contains("HLT_IsoMu24_v1") ) && Ntp->HLTDecision(iTrigger) ) { HLTOk = true; TauIsoMu = true;}
  }
  



  random_num = rndm.Rndm();
  


  
  for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
    TString L1TriggerName = Ntp->L1Name(il1);
    //    if(Ntp->L1Decision(il1))    std::cout<<" l1 name  "<< Ntp->L1Name(il1) <<  " =   "<<Ntp->L1Decision(il1) <<std::endl;
    //    std::cout<<" random number   "<< random_num << std::endl;
    if(L1TriggerName.Contains("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4") && Ntp->L1Decision(il1)) { DoubleMu0Fired = true; }
    if(L1TriggerName.Contains("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9") && Ntp->L1Decision(il1)) { TripleMuFired = true; }
    if( id==1 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMu4Fired = true; }
    if( id!=1 && random_num>0.30769 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) { DoubleMu4Fired = true;}
    if( id!=1 && random_num<0.30769 && L1TriggerName.Contains("L1_DoubleMu4_SQ_OS_dR_Max1p2") && Ntp->L1Decision(il1)) {
      randomFailed = true;
    }
  }




  bool SingleMu(false);
  for(int il1=0; il1 < Ntp->NL1Seeds(); il1++){
    TString L1TriggerName = Ntp->L1Name(il1);
    //    if(Ntp->L1Decision(il1))    std::cout<<" l1 name  "<< Ntp->L1Name(il1) <<  " =   "<<Ntp->L1Decision(il1) <<std::endl;
    //    std::cout<<" random number   "<< random_num << std::endl;
    if(L1TriggerName.Contains("L1_SingleMu") && Ntp->L1Decision(il1)) { SingleMu = true; }
  }


  value.at(L1_TriggerOk)=0;
  value.at(HLT_TriggerOk)=0;
  if (DoubleMu0Fired || DoubleMu4Fired) {DoubleMuFired = true;}
  if (DoubleMuFired || TripleMuFired) L1Ok = true;

  value.at(L1_TriggerOk)=(L1Ok);
  //  pass.at(L1_TriggerOk)=true;//(value.at(L1_TriggerOk)==cut.at(L1_TriggerOk));



  //  value.at(L1_TriggerOk)=(SingleMu);
  pass.at(L1_TriggerOk)=(value.at(L1_TriggerOk)==cut.at(L1_TriggerOk));
  

  value.at(HLT_TriggerOk)=(Tau3muHLT);
  //  value.at(HLT_TriggerOk)=(HLTOk);
  pass.at(HLT_TriggerOk)=(value.at(HLT_TriggerOk)==cut.at(HLT_TriggerOk));




  SingleVsThreeMuTrigger.at(t).Fill(Tau3muHLT*L1Ok,TauIsoMu*SingleMu);

  //  value.at(SignalCandidate)=0;
  //  if(Ntp->NThreeMuons() !=0)
    value.at(SignalCandidate) = Ntp->NThreeMuons();

  int  signal_idx=-1;
  double min_chi2(99.);

  for(unsigned int i_idx =0; i_idx < Ntp->NThreeMuons(); i_idx++){
    if(Ntp->Vertex_Signal_KF_Chi2(i_idx) < min_chi2){
      min_chi2 = Ntp->Vertex_Signal_KF_Chi2(i_idx);
      signal_idx = i_idx;
    }
  }

  TLorentzVector Tau3MuLV(0,0,0,0);
  pass.at(SignalCandidate) = (value.at(SignalCandidate) >= cut.at(SignalCandidate));


  value.at(TripletPT)=0;
  if(signal_idx!=-1)
    {
      Tau3MuLV = Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0))+
        Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1))+
        Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2));
      value.at(TripletPT) = Tau3MuLV.Pt();

    }

  pass.at(TripletPT) = ( value.at(TripletPT) >= cut.at(TripletPT) );
  //  std::cout<<"   index  "<< signal_idx << "  id  "<< id << std::endl;

  std::vector<int> Muons;
  std::vector<int> Muons_OppositeHemisphere;
  std::vector<int> Muons_OppositeHemisphere_OppositeCharge;



  value.at(nMuons)  = 0;
  for(unsigned int imu=0; imu < Ntp->NMuons(); imu++)
    {
      if(signal_idx!=-1)
	{
	  if(Ntp->Muon_P4(imu).Pt() > 4 && fabs(Ntp->Muon_P4(imu).Eta()) < 2.5 && Ntp->Muon_isPFMuon(imu) 
	     && ( Ntp->Muon_isGlobalMuon(imu) || Ntp->Muon_isTrackerMuon(imu) ) &&
	     Ntp->Muon_P4(imu).DeltaR(Tau3MuLV) > 0.5  )Muons_OppositeHemisphere.push_back(imu);
	}
    }

  value.at(nMuons)  = Muons_OppositeHemisphere.size();
  pass.at(nMuons)   = ( value.at(nMuons) >= cut.at(nMuons) );


  value.at(OSCharge)        =  0;
  value.at(Tau3MuIsolation) = -1;
  value.at(MuonIsolation)   = -1;
  value.at(VisMass)         = -1;

  value.at(TriggerMatch) = 0;
  if(signal_idx!=-1)
    {


      int index_mu_1 = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0);  
      int index_mu_2 = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1);  
      int index_mu_3 = Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2);  

      TLorentzVector TripletmuLV = Ntp->Muon_P4(index_mu_1) +  Ntp->Muon_P4(index_mu_2) +  Ntp->Muon_P4(index_mu_3);
      
      value.at(Tau3MuIsolation) = TripletmuLV.Pt()/  (TripletmuLV.Pt()  + Ntp->Muon_RelIso(index_mu_1) +
 						                          Ntp->Muon_RelIso(index_mu_2) +
						                          Ntp->Muon_RelIso(index_mu_3) );


      if(id ==570){
    
	std::cout<<"------------------------------- "<< std::endl;
	std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Ntp->Muon_P4(index_mu_1)) << "  midx   " << Ntp->MCParticle_midx(Ntp->getMatchTruthIndex(Ntp->Muon_P4(index_mu_1))) <<std::endl;
	std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Ntp->Muon_P4(index_mu_2)) << "  midx   " << Ntp->MCParticle_midx(Ntp->getMatchTruthIndex(Ntp->Muon_P4(index_mu_2))) <<std::endl;
	std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Ntp->Muon_P4(index_mu_3)) << "  midx   " << Ntp->MCParticle_midx(Ntp->getMatchTruthIndex(Ntp->Muon_P4(index_mu_3))) <<std::endl;
      
	Ntp->Muon_P4(index_mu_1).Print(); std::cout<<" idx1:  "<<Ntp->getMatchTruthIndex(Ntp->Muon_P4(index_mu_1)) << std::endl;
	Ntp->Muon_P4(index_mu_2).Print(); std::cout<<" idx2:  "<<Ntp->getMatchTruthIndex(Ntp->Muon_P4(index_mu_2)) << std::endl;
	Ntp->Muon_P4(index_mu_3).Print(); std::cout<<" idx3:  "<<Ntp->getMatchTruthIndex(Ntp->Muon_P4(index_mu_3)) << std::endl;
      


	for(auto i:Ntp->MCParticle_childpdgid(Ntp->MCParticle_midx(Ntp->getMatchTruthIndex(Ntp->Muon_P4(index_mu_1)))))
	  std::cout<<"  daughters of first mother "<< i << std::endl;

	for(auto i:Ntp->MCParticle_childpdgid(Ntp->MCParticle_midx(Ntp->getMatchTruthIndex(Ntp->Muon_P4(index_mu_2)))))
	  std::cout<<"  daughters of second  mother "<< i << std::endl;

	for(auto i:Ntp->MCParticle_childpdgid(Ntp->MCParticle_midx(Ntp->getMatchTruthIndex(Ntp->Muon_P4(index_mu_3)))))
	  std::cout<<"  daughters of third  mother "<< i << std::endl;
	//	for(unsigned int i=0; i < Ntp->MCParticle_childpdgid(); i++)
	//{
	//  std::
	//}

	//	std::cout<<"  Signal  "<<  Ntp->NMCSignalParticles() << std::endl;
	//	for(unsigned int i=0; i < Ntp->NMCParticles(); i++){
	//	  std::cout<<"   Nchilds   "   << Ntp->MCSignalParticle_Nchilds(i) << std::endl;
	//	}

	//	for(unsigned int i=0; i < Ntp->NMCParticles(); i++)
	//	  {

	    //	    std::cout<<"  i   "<<i << " pdg   "<< Ntp->MCParticle_pdgid(i) << std::endl;
	    //	    if(Ntp->MCParticle_pdgid(i) == 23)
	    //	      {
		//		for(unsigned int j =0; j<)

	    //	      }

	//	  }

	Ntp->printMCDecayChainOfEvent(true, true, true, true);
	std::cout<< "\n\n\n\n\n\n";
      }    
    
    //Trigger Matching

      int triggerCheck = 0.1;
      if(pass.at(HLT_TriggerOk))
	{
	  vector<TLorentzVector> trigobjTriplet;
	  for (int i=0; i<Ntp->NTriggerObjects(); i++)
	    {
	      TString name = Ntp->TriggerObject_name(i);
	      //	if (!(name.Contains("tau3muDisplaced3muFltr"))) continue;
	      TLorentzVector tmp;
	      tmp.SetPtEtaPhiM(Ntp->TriggerObject_pt(i), Ntp->TriggerObject_eta(i), Ntp->TriggerObject_phi(i), PDG_Var::Muon_mass());
	      trigobjTriplet.push_back(tmp);
	    }
	  std::vector<TLorentzVector> muonTriplet;
	  muonTriplet.push_back(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0)));
	  muonTriplet.push_back(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1)));
	  muonTriplet.push_back(Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2)));
	  //	  std::cout<<" trigobjTriplet.size()  "<< trigobjTriplet.size() << std::endl;
	  if (trigobjTriplet.size()>=3) triggerCheck = Ntp->triggerMatchTriplet(muonTriplet, trigobjTriplet).first;
	}

      value.at(TriggerMatch) = true;//triggerCheck;      

      
	if( pass.at(nMuons))
	  {
	    for(auto i : Muons_OppositeHemisphere)
	      {
		
		int Tau3MuCharge = Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)) +
		  Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)) +
		  Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2));
		
		if(Ntp->Muon_charge(i)*Tau3MuCharge == -1) Muons_OppositeHemisphere_OppositeCharge.push_back(i);
	      }
	    value.at(OSCharge) = Muons_OppositeHemisphere_OppositeCharge.size();
	  }
	
	
    }
  
    pass.at(TriggerMatch) = (value.at(TriggerMatch)  ==  cut.at(TriggerMatch));
    //    std::cout<<" (value.at(TriggerMatch)   " << value.at(TriggerMatch)  << " cut.at(TriggerMatch)   " << cut.at(TriggerMatch) << "  ??   " 
    //	     <<(value.at(TriggerMatch)  ==  cut.at(TriggerMatch) )<< std::endl;
    pass.at(OSCharge) = (value.at(OSCharge) >= cut.at(OSCharge)); 



    if(pass.at(OSCharge))
      {
	unsigned int muon_index = Muons_OppositeHemisphere_OppositeCharge.at(0);
	value.at(MuonIsolation)   = Ntp->Muon_P4(muon_index).Pt()  /  (Ntp->Muon_P4(muon_index).Pt()  +   Ntp->Muon_RelIso(muon_index) ); 
	value.at(VisMass) = (Tau3MuLV + Ntp->Muon_P4(muon_index)).M();
      }

    pass.at(Tau3MuIsolation) = (value.at(Tau3MuIsolation) > cut.at(Tau3MuIsolation));
    pass.at(MuonIsolation)   = (value.at(MuonIsolation)   > cut.at(MuonIsolation));
    pass.at(VisMass)         = (value.at(VisMass) > 30 && value.at(VisMass) < 100);
    //    std::cout<<"  ----  "<< std::endl;
    /*    for(unsigned int i = 0 ; i< NCuts; i++)
      {

	std::cout<<"  i "<<i<< "    "<<pass.at(i) << std::endl;
	}*/
    
    double wobs=1;
    double w;  
             
    if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
    else{w=1;}
  

  
  bool status=AnalysisCuts(t,w,wobs);

  if(status){ 
    //    std::cout<<"pass nMuons "<< pass.at(nMuons) <<    "  triplet pT    "<< pass.at(TripletPT)  << std::endl;
    //    std::cout<<"   how many muons i have  " << Muons_OppositeHemisphere_OppositeCharge .size() << std::endl;
    unsigned int muon_idx = Muons_OppositeHemisphere_OppositeCharge.at(0);


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


    TLorentzVector MuLV = Ntp->Muon_P4(muon_idx);



    LorentzVectorParticle Tau3MuLVP = Ntp->Tau3mu_LVP(  signal_idx );
    TVector3 Neutrino_Vect(Ntp->METEt()*TMath::Cos(Ntp->METPhi()),Ntp->METEt()*TMath::Sin(Ntp->METPhi()),Ntp->METEt()/TMath::Tan(MuLV.Theta()));
    TLorentzVector Neutrino_LV(Neutrino_Vect,Neutrino_Vect.Mag());



    float RelativeIsolationMu1 = Ntp->Muon_RelIso(muon_1_idx);
    float RelativeIsolationMu2 = Ntp->Muon_RelIso(muon_2_idx);
    float RelativeIsolationMu3 = Ntp->Muon_RelIso(muon_3_idx);


    //////////// kinematics 

    TripletPt.at(t).Fill(Tau3muLV.Pt(),1);
    OppositeMuonPt.at(t).Fill(MuLV.Pt(),1);

    TripletEta.at(t).Fill(Tau3muLV.Eta(),1);
    OppositeMuonEta.at(t).Fill(MuLV.Eta(),1);
    //////////// kinematics 


    Tau3MuRelativeIsolation.at(t).Fill(    Tau3muLV.Pt()/(RelativeIsolationMu1 + RelativeIsolationMu2 + RelativeIsolationMu3 + Tau3muLV.Pt()),1);
    OppositeMuRelativeIsolation.at(t).Fill(    Ntp->Muon_P4(muon_idx).Pt()/(Ntp->Muon_RelIso(muon_idx) +  Ntp->Muon_P4(muon_idx).Pt()),1);

    VisibleDiTauMass.at(t).Fill((MuLV + Tau3muLV).M(), 1);
    MTT.at(t).Fill( (Tau3muLV + MuLV  + Neutrino_LV).M(), 1);


    mini_m3m = TauRefitLV.M();
    mini_dataMCtype = id;
    mini_event_weight = 1;
    mini_TripletPt = TauRefitLV.Pt();
    mini_MuonPt = MuLV.Pt();
    mini_Tau3MuRelativeIsolation = Tau3muLV.Pt()/(RelativeIsolationMu1 + RelativeIsolationMu2 + RelativeIsolationMu3 + Tau3muLV.Pt());
    mini_OppositeMuRelativeIsolation = Ntp->Muon_P4(muon_idx).Pt()/(Ntp->Muon_RelIso(muon_idx) +  Ntp->Muon_P4(muon_idx).Pt());
    mini_ditaumass = (MuLV + Tau3muLV).M();
  
    T3MMiniTree->Fill();


    bool PlotMCOnly(false);  // and blind for data
    if(id!=1) PlotMCOnly = true;
    if(id==1 && ( (TauRefitLV.M() > 1.1 && TauRefitLV.M() < 1.74) or (TauRefitLV.M() > 1.812 && TauRefitLV.M() < 2.2)) ) PlotMCOnly=true;


    if(PlotMCOnly)  TripletMass.at(t).Fill(TauRefitLV.M(),1);
    TripletMassWR.at(t).Fill(TauRefitLV.M(),1);

    TLorentzVector OppositeSideLV = MuLV;
    if(id != 1)
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
      }


  }
}


void  ZTau3MuTaumu::Finish(){

  //*** write down the T3MMiniTree.root for statistical analysis
  TString out_file_name  = "MVA_Mini_Tree_"+ AnalysisName+".root";
  T3MFMiniTree = new TFile(out_file_name,"recreate");
  T3MMiniTree->SetDirectory(T3MFMiniTree);
  T3MFMiniTree->Write();
  T3MFMiniTree->Close();

  Selection::Finish();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
}





