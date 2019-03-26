#include "MCStudy.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

MCStudy::MCStudy(TString Name_, TString id_):
  Selection(Name_,id_)
{


  // This is a class constructor;
}

MCStudy::~MCStudy(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }


  Logger(Logger::Info) << "complete." << std::endl;
}

void  MCStudy::Configure(){

  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==L1SeedOk)    cut.at(L1SeedOk)=1;
    if(i==HLTOk)    cut.at(HLTOk)=1;
    if(i==PrimeVtx)     cut.at(PrimeVtx)=5; // Here for example we place cut value on number of PVs
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
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,61,-0.5,60.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,61,-0.5,60.5,hlabel,"Events"));
    }
    else if(i==L1SeedOk){
      title.at(i)="L1 seed ";
      hlabel="L1 triggers";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_L1SeedOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_L1SeedOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }

    else if(i==HLTOk){
      title.at(i)="HLT trigger ";
      hlabel="DoubleMu3_Trk_Tau3mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms
  // Book here your analysis histrogramms, a good style is to follow selfexplanatory convention
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",66,-0.5,65.5,"Number of Vertices","Events");
  MuonsPt=HConfig.GetTH1D(Name+"_MuonsPt","All Muons Pt",25,0,50,"p_{T} of all stored Muons, GeV","Events");
  Category=HConfig.GetTH1D(Name+"_Category","Category",2,0.5,2.5,"1 - #mu#mu#mu; 2 - #mu#mu+track ","Events");
  MuonsPtRatio=HConfig.GetTH1D(Name+"_MuonsPtRatio","Ratio of 2 muons pt",50,0.1,1.2,"Ratio of first and second muon p_{T}","Events");
  MuonsEta=HConfig.GetTH1D(Name+"_MuonsEta","All Muons #eta",25,-2.5,2.5,"Rapidity of all stored Muons","Events");
  PhiMass=HConfig.GetTH1D(Name+"_PhiMass","#mu#mu mass",75,0.2,1.4,"Mass of the #mu#mu pair, GeV","Events");
  TripleMass=HConfig.GetTH1D(Name+"_TripleMass","#mu#mu + track mass",50,1.7,2.1,"Mass of the #mu#mu + track, GeV","Events");
  PhiMassVsDsMass=HConfig.GetTH2D(Name+"_PhiMassVsDsMass","#mu#mu Mass vs. #mu#mu + track mass",75,0.2,1.4,50,1.7,2.1,"M_{#mu#mu}, GeV","M_{#mu#mu + track}, GeV");

  FirstMuonsPt=HConfig.GetTH1D(Name+"_FirstMuonsPt","FirstMuonPt",25,0,50,"First  #mu p_{T}, GeV","Events");
  SecondMuonsPt=HConfig.GetTH1D(Name+"_SecondMuonsPt","SecondMuonPt",25,0,50,"Second #mu p_{T}, GeV","Events");
  ThirdMuonsPt=HConfig.GetTH1D(Name+"_ThirdMuonsPt","ThirdMuonPt",25,0,50,"Third  #mu p_{T}, GeV","Events");


  FirstMuonsEta=HConfig.GetTH1D(Name+"_FirstMuonsEta","FirstMuonEta",25,-2.5,2.5,"First  #mu  rapidity","Events");
  SecondMuonsEta=HConfig.GetTH1D(Name+"_SecondMuonsEta","SecondMuonEta",25,-2.5,2.5,"Second #mu  rapidity","Events");
  ThirdMuonsEta=HConfig.GetTH1D(Name+"_ThirdMuonsEta","ThirdMuonEta",25,-2.5,2.5,"Third  #mu rapidity","Events");





  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  MCStudy::Store_ExtraDist(){ 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output
  Extradist1d.push_back(&NVtx);
  Extradist1d.push_back(&MuonsPt);
  Extradist1d.push_back(&MuonsEta);
  Extradist1d.push_back(&MuonsPtRatio);
  Extradist1d.push_back(&PhiMass);
  Extradist1d.push_back(&TripleMass);
  Extradist1d.push_back(&Category);
  Extradist2d.push_back(&PhiMassVsDsMass);

  Extradist1d.push_back(&FirstMuonsPt);
  Extradist1d.push_back(&SecondMuonsPt);
  Extradist1d.push_back(&ThirdMuonsPt);

  Extradist1d.push_back(&FirstMuonsEta);
  Extradist1d.push_back(&SecondMuonsEta);
  Extradist1d.push_back(&ThirdMuonsEta);



}


void  MCStudy::doEvent(){ 
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection
  std::cout<<"----------------id "<<id << std::endl;
  value.at(L1SeedOk) = 0;
  value.at(HLTOk) = 0;
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);
    if(HLT.Contains("DoubleMu3_Trk_Tau3mu") && Ntp->HLTDecision(iTrigger) == 1)value.at(HLTOk)=Ntp->HLTDecision(iTrigger);
  }
  for(int l1iTrigger=0; l1iTrigger < Ntp->NL1Seeds(); l1iTrigger++){
    TString L1 = Ntp->L1Name(l1iTrigger);
    if(L1.Contains("L1_DoubleMu0er1p4_dEta_Max1p8_OS") && Ntp->L1Decision(l1iTrigger) == 1)value.at(L1SeedOk)=Ntp->L1Decision(l1iTrigger);
    if(L1.Contains("L1_DoubleMu_10_0_dEta_Max1p8") && Ntp->L1Decision(l1iTrigger) == 1)value.at(L1SeedOk)=Ntp->L1Decision(l1iTrigger);
    if(L1.Contains("L1_TripleMu0") && Ntp->L1Decision(l1iTrigger) == 1)value.at(L1SeedOk)=Ntp->L1Decision(l1iTrigger);
  }

  value.at(PrimeVtx)=Ntp->NVtx(); 
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx)); 
  
  pass.at(L1SeedOk)= (value.at(L1SeedOk)==cut.at(L1SeedOk)); 
  pass.at(HLTOk)= (value.at(HLTOk)==cut.at(HLTOk)); 
  
  double wobs=1;
  double w;  
             
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}

  /*

  unsigned int               NMCSignalParticles(){return Ntp->MCSignalParticle_p4->size();}
  TLorentzVector             MCSignalParticle_p4(unsigned int i){return TLorentzVector(Ntp->MCSignalParticle_p4->at(i).at(1),Ntp->MCSignalParticle_p4->at(i).at(2),Ntp->MCSignalParticle_p4->at(i).at(3),Ntp->MCSignalParticle_p4->at(i).at(0));}
  int                        MCSignalParticle_pdgid(unsigned int i){return Ntp->MCSignalParticle_pdgid->at(i);}
  int                        MCSignalParticle_charge(unsigned int i){return Ntp->MCSignalParticle_charge->at(i);}
  std::vector<unsigned int>  MCSignalParticle_Tauidx(unsigned int i){return Ntp->MCSignalParticle_Tauidx->at(i);}
  // Tau decays (Tau is first element of vector)
  int NMCTaus(){return Ntp->MCTauandProd_p4->size();}
  TLorentzVector MCTau_p4(unsigned int i){return MCTauandProd_p4(i,0);}
  int MCTau_pdgid(unsigned int i){return MCTauandProd_pdgid(i,0);}
  int MCTau_charge(unsigned int i){return MCTauandProd_charge(i,0);}


  //Tau and decay products
  int NMCTauDecayProducts(unsigned int i){if(0<=i && i<(unsigned int)NMCTaus()) return Ntp->MCTauandProd_p4->at(i).size(); return 0;}
  TLorentzVector MCTauandProd_p4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->MCTauandProd_p4->at(i).at(j).at(1),Ntp->MCTauandProd_p4->at(i).at(j).at(2),Ntp->MCTauandProd_p4->at(i).at(j).at(3),Ntp->MCTauandProd_p4->at(i).at(j).at(0));}
  int MCTauandProd_pdgid(unsigned int i, unsigned int j){return Ntp->MCTauandProd_pdgid->at(i).at(j);}
  int MCTauandProd_charge(unsigned int i, unsigned int j){return Ntp->MCTauandProd_charge->at(i).at(j);}
  */


  std::cout<<"N signal particles  "<< Ntp->NMCSignalParticles() << std::endl;
    for(int isig=0; isig< Ntp->NMCSignalParticles() ; isig++){
      std::cout<<"---------------------------- Signal Particle ID:  "<< Ntp->MCSignalParticle_pdgid(isig) << " Mass  "<< Ntp->MCSignalParticle_p4(isig).M() << std::endl;
      std::cout<<" TauIDx size  "<< Ntp->MCSignalParticle_Tauidx(isig).size() <<"   NTaus  "<< Ntp-> NMCTaus()<< std::endl;
      std::cout<<"Signal particle P4 ";Ntp->MCSignalParticle_p4(isig).Print();



      for(int ichi=0; ichi < Ntp->MCSignalParticle_Nchilds(isig); ichi++){
	std::cout<<"child n:  "<< ichi <<"  pdg  "<< Ntp->MCSignalParticle_childpdgid(isig,ichi)<< "  E  "<< Ntp->MCSignalParticle_child_p4(isig,ichi).E() << std::endl;

      }
      for(int itau=0; itau< Ntp-> NMCTaus(); itau++){
	std::cout<<"Tau Number:   "<< itau<< std::endl;
	Ntp->MCTau_p4(itau).Print();

      }
    }


    for(int itau=0; itau< Ntp-> NMCTaus(); itau++){
      std::cout<<"Tau Number:   "<< itau<< " Mother ID    "<< Ntp-> MCTau_midx(itau)<<std::endl;

      for(int ida =0 ; ida < Ntp->NMCTauDecayProducts(itau); ida++){
	std::cout<<"  tau product id   "<< Ntp->MCTauandProd_pdgid(itau,ida)<< std::endl;
	std::cout<<"Momenta ";
	Ntp->MCTauandProd_p4(itau,ida).Print();
      }
    }




    if(Ntp->NThreeMuons()==2){
      unsigned int Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(0);
      unsigned int Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(1);
      unsigned int Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(2);


      unsigned int Muon_index_11=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(1)).at(0);
      unsigned int Muon_index_21=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(1)).at(1);
      unsigned int Muon_index_31=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(1)).at(2);


      std::cout<<"RECO MUONS      11111111111111111111 a"<<std::endl;
      Ntp->Muon_P4(Muon_index_1).Print();
      Ntp->Muon_P4(Muon_index_2).Print();
      Ntp->Muon_P4(Muon_index_3).Print();


      std::cout<<"RECO MUONS      22222222222222222222 a"<<std::endl;
      Ntp->Muon_P4(Muon_index_11).Print();
      Ntp->Muon_P4(Muon_index_21).Print();
      Ntp->Muon_P4(Muon_index_31).Print();

    }





  bool status=AnalysisCuts(t,w,wobs);
  if(status){
    NVtx.at(t).Fill(Ntp->NVtx(),w);

    //    std::cout<<"N 2 mu + track  "<< Ntp->NTwoMuonsTrack() << " N three mu  "<< Ntp->NThreeMuons() << "  Number of Vertices  " << Ntp->NumberOfSVertices() <<std::endl;
    for(unsigned int iMuon=0; iMuon < Ntp->NMuons(); iMuon++){ // loop over muons
      MuonsPt.at(t).Fill(Ntp->Muon_P4(iMuon).Pt(),w);
      MuonsEta.at(t).Fill(Ntp->Muon_P4(iMuon).Eta(),w);
    }
    if(Ntp->NThreeMuons()!=0) Category.at(t).Fill(1);
    if(Ntp->NTwoMuonsTrack()!=0) Category.at(t).Fill(2);

    //    std::cout<<"  "<< Ntp->NThreeMuons() << "   "<< Ntp->NTwoMuonsTrack()<<std::endl;
    //    Ntp->Vertex_signal_KF_refittedTracksP4(0,0).Print();
    //    Ntp->Vertex_signal_KF_refittedTracksP4(0,1).Print();
    //    Ntp->Vertex_signal_KF_refittedTracksP4(0,2).Print();
    if(Ntp->NThreeMuons()!=0){
      unsigned int Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(0);
      unsigned int Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(1);
      unsigned int Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(2);

      double pt1 = Ntp->Muon_P4(Muon_index_1).Pt();
      double pt2 = Ntp->Muon_P4(Muon_index_2).Pt();
      double pt3 = Ntp->Muon_P4(Muon_index_3).Pt();
      FirstMuonsPt.at(t).Fill(pt1,1);
      SecondMuonsPt.at(t).Fill(pt2,1);
      ThirdMuonsPt.at(t).Fill(pt3,1);


      FirstMuonsEta.at(t).Fill( Ntp->Muon_P4(Muon_index_1).Eta(),1);
      SecondMuonsEta.at(t).Fill(  Ntp->Muon_P4(Muon_index_2).Eta(),1);
      ThirdMuonsEta.at(t).Fill(  Ntp->Muon_P4(Muon_index_3).Eta(),1);
    }


    double deltaMass(999.);
    unsigned int pair_index(0);
    if(Ntp->NTwoMuonsTrack()!=0){
    for(unsigned int i2M=0; i2M < Ntp->NTwoMuonsTrack(); i2M++){

      unsigned int muon_1 =  Ntp-> TwoMuonsTrackMuonIndices(i2M).at(0);
      unsigned int muon_2 =  Ntp-> TwoMuonsTrackMuonIndices(i2M).at(1);

      if( fabs((Ntp->Muon_P4(muon_1)  + Ntp->Muon_P4(muon_2)).M()  - PDG_Var::Phi_mass())< deltaMass){
	deltaMass =  fabs((Ntp->Muon_P4(muon_1)  + Ntp->Muon_P4(muon_2)).M()  - PDG_Var::Phi_mass());
	pair_index = i2M; // this is an index of the candidate with best mumu mass
      }
    }
    

    MuonsPtRatio.at(t).Fill(Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(0)).Pt()/Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1)).Pt(),w );
    PhiMass.at(t).Fill((Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))).M(), w);
    TripleMass.at(t).Fill((Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))+ 
			   Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(pair_index).at(0))).M(), w);

    double phimass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))).M();
    double dsmass = (Ntp->Muon_P4( Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(0))  + Ntp->Muon_P4(Ntp-> TwoMuonsTrackMuonIndices(pair_index).at(1))+
		     Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(pair_index).at(0))).M();

    PhiMassVsDsMass.at(t).Fill(phimass, dsmass);

    }
  }
}



void  MCStudy::Finish(){
  Selection::Finish();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
}





