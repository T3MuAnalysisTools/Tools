#include "ThreeMu.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

ThreeMu::ThreeMu(TString Name_, TString id_):
  Selection(Name_,id_)
{


  // This is a class constructor;
}

ThreeMu::~ThreeMu(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }


  Logger(Logger::Info) << "complete." << std::endl;
}

void  ThreeMu::Configure(){

  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==HLTOk)               cut.at(HLTOk)=1;
    if(i==ThreeMuCandidate)    cut.at(ThreeMuCandidate)=1;
  }
  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    // book here the N-1 and N-0 histrogramms for each cut

    if(i==HLTOk){
      title.at(i)="HLT trigger ";
      hlabel="DoubleMu3_Trk_Tau3mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==ThreeMuCandidate){
      title.at(i)="ThreeMuCandidate is found ";
      hlabel="is Three Mu Candidate";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ThreeMuCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ThreeMuCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms  // Book here your analysis histrogramms, a good style is to follow selfexplanatory convention

  Muon1Pt =HConfig.GetTH1D(Name+"_Muon1Pt","Muon1Pt",75,0,25,"  #mu_{1} p_{T}, GeV","Events");
  Muon2Pt =HConfig.GetTH1D(Name+"_Muon2Pt","Muon2Pt",75,0,20,"  #mu_{2} p_{T}, GeV","Events");
  Muon3Pt =HConfig.GetTH1D(Name+"_Muon3Pt","Muon3Pt",75,0,15,"  #mu_{3} p_{T}, GeV","Events");

  Muon1isGlob =HConfig.GetTH1D(Name+"_Muon1isGlob","Muon1isGlob",2,-0.5,1.5,"  #mu_{1} is global muon","Events");
  Muon2isGlob =HConfig.GetTH1D(Name+"_Muon2isGlob","Muon2isGlob",2,-0.5,1.5,"  #mu_{2} is global muon","Events");
  Muon3isGlob =HConfig.GetTH1D(Name+"_Muon3isGlob","Muon3isGlob",2,-0.5,1.5,"  #mu_{3} is global muon","Events");

  Muon1isTrack =HConfig.GetTH1D(Name+"_Muon1isTrack","Muon1isTrack",2,-0.5,1.5,"  #mu_{1} is tracker muon","Events");
  Muon2isTrack =HConfig.GetTH1D(Name+"_Muon2isTrack","Muon2isTrack",2,-0.5,1.5,"  #mu_{2} is tracker muon","Events");
  Muon3isTrack =HConfig.GetTH1D(Name+"_Muon3isTrack","Muon3isTrack",2,-0.5,1.5,"  #mu_{3} is tracker muon","Events");

  Muon1kink =HConfig.GetTH1D(Name+"_Muon1kink","Muon1kink",50,0.,50,"  #mu_{1} kink","Events");
  Muon2kink =HConfig.GetTH1D(Name+"_Muon2kink","Muon2kink",50,0.,50,"  #mu_{2} kink","Events");
  Muon3kink =HConfig.GetTH1D(Name+"_Muon3kink","Muon3kink",50,0.,50,"  #mu_{3} kink","Events");

  Muon1InOutTrackMatch =HConfig.GetTH1D(Name+"_Muon1InOutTrackMatch","Muon1InOutTrackMatch",50,0.,10,"  #mu_{1} inner and outer tracker match","Events");
  Muon2InOutTrackMatch =HConfig.GetTH1D(Name+"_Muon2InOutTrackMatch","Muon2InOutTrackMatch",50,0.,10,"  #mu_{2} inner and outer tracker match","Events");
  Muon3InOutTrackMatch =HConfig.GetTH1D(Name+"_Muon3InOutTrackMatch","Muon3InOutTrackMatch",50,0.,10,"  #mu_{3} inner and outer tracker match","Events");

  Muon1Eta=HConfig.GetTH1D(Name+"_Muon1Eta","Muon1Eta",25,-2.6,2.6,"#mu_{1}  rapidity","Events");
  Muon2Eta=HConfig.GetTH1D(Name+"_Muon2Eta","Muon2Eta",25,-2.6,2.6,"#mu_{2}  rapidity","Events");
  Muon3Eta=HConfig.GetTH1D(Name+"_Muon3Eta","Muon3Eta",25,-2.6,2.6,"#mu_{3}  rapidity","Events");

  Muon1StandardSelector=HConfig.GetTH1D(Name+"_Muon1StandardSelector","Muon1StandardSelector",23,-0.5,22.5,"#mu_{1} standard selector; bin 0 - no ID","Events");
  Muon2StandardSelector=HConfig.GetTH1D(Name+"_Muon2StandardSelector","Muon2StandardSelector",23,-0.5,22.5,"#mu_{2} standard selector; bin 0 - no ID","Events");
  Muon3StandardSelector=HConfig.GetTH1D(Name+"_Muon3StandardSelector","Muon3StandardSelector",23,-0.5,22.5,"#mu_{3} standard selector; bin 0 - no ID","Events");


  TauEta=HConfig.GetTH1D(Name+"_TauEta","TauEta",30,-2.6,2.6,"#tau  rapidity","Events");
  TauPt =HConfig.GetTH1D(Name+"_TauPt","TauPt",30,0,50,"  #tau p_{T}, GeV","Events");
  TauP =HConfig.GetTH1D(Name+"_TauP","TauP",30,0,50,"  #tau |p|, GeV","Events");
  TauMass =HConfig.GetTH1D(Name+"_TauMass","#tau lepton mass",50,1.5,1.9,"  M_{#tau} , GeV","Events");
  TauMassResolution=HConfig.GetTH1D(Name+"_TauMassResolution","TauMassResolution",50,-0.2,0.2,"#Delta M_{#tau}  (reco - mc)/mc ","Events");
  TauMassRefit =HConfig.GetTH1D(Name+"_TauMassRefit","Refit #tau lepton mass",50,1.5,1.9,"KF refit  M_{#tau} , GeV","Events");
  TauMassResolutionRefit=HConfig.GetTH1D(Name+"_TauMassResolutionRefit","TauMassResolutionRefit",50,-0.2,0.2,"KF refit #Delta M_{#tau}  (reco - mc)/mc ","Events");

  VertexChi2KF=HConfig.GetTH1D(Name+"_VertexChi2KF","VertexChi2KF",55,0,10,"#chi^{2} of the KF SV","Events");
  VertexChi2AF=HConfig.GetTH1D(Name+"_VertexChi2AF","VertexChi2AF",55,0,10,"#chi^{2} of the AF SV","Events");

  Muon1PtResolution=HConfig.GetTH1D(Name+"_Muon1PtResolution","Muon1PtResolution",50,-0.2,0.2," #Delta p_{T}(#mu_{1})  (reco - mc)/mc ","Events");
  Muon2PtResolution=HConfig.GetTH1D(Name+"_Muon2PtResolution","Muon2PtResolution",50,-0.2,0.2," #Delta p_{T}(#mu_{2})  (reco - mc)/mc ","Events");
  Muon3PtResolution=HConfig.GetTH1D(Name+"_Muon3PtResolution","Muon3PtResolution",50,-0.2,0.2," #Delta p_{T}(#mu_{3})  (reco - mc)/mc  ","Events");

  Muon1EtaResolution=HConfig.GetTH1D(Name+"_Muon1EtaResolution","Muon1EtaResolution",50,-0.05,0.05," #Delta #eta(#mu_{1})  (reco - mc)/mc  ","Events");
  Muon2EtaResolution=HConfig.GetTH1D(Name+"_Muon2EtaResolution","Muon2EtaResolution",50,-0.05,0.05," #Delta #eta(#mu_{1})  (reco - mc)/mc   ","Events");
  Muon3EtaResolution=HConfig.GetTH1D(Name+"_Muon3EtaResolution","Muon3EtaResolution",50,-0.05,0.05," #Delta #eta(#mu_{1})  (reco - mc)/mc   ","Events");


  Muon1DRToTruth=HConfig.GetTH1D(Name+"_Muon1DRToTruth","Muon1DRToTruth",20,0,0.02,"reco - mc #mu_{1} #Delta R","Events");
  Muon2DRToTruth=HConfig.GetTH1D(Name+"_Muon2DRToTruth","Muon2DRToTruth",20,0,0.02,"reco - mc #mu_{1} #Delta R","Events");
  Muon3DRToTruth=HConfig.GetTH1D(Name+"_Muon3DRToTruth","Muon3DRToTruth",20,0,0.02,"reco - mc #mu_{1} #Delta R","Events");

  MuPair1_vertex_chi2=HConfig.GetTH1D(Name+"_MuPair1_vertex_chi2","MuPair1_vertex_chi2",50,0,5,"KF  #chi^{2} of first #mu pair","Events");
  MuPair2_vertex_chi2=HConfig.GetTH1D(Name+"_MuPair2_vertex_chi2","MuPair2_vertex_chi2",50,0,5,"KF  #chi^{2} of second #mu pair","Events");
  MuPair3_vertex_chi2=HConfig.GetTH1D(Name+"_MuPair3_vertex_chi2","MuPair3_vertex_chi2",50,0,5,"KF  #chi^{2} of third #mu pair","Events");

  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  ThreeMu::Store_ExtraDist(){ 

  Extradist1d.push_back(&Muon1Pt);
  Extradist1d.push_back(&Muon2Pt);
  Extradist1d.push_back(&Muon3Pt);

  Extradist1d.push_back(&Muon1Eta);
  Extradist1d.push_back(&Muon2Eta);
  Extradist1d.push_back(&Muon3Eta);

  Extradist1d.push_back(&Muon1StandardSelector);
  Extradist1d.push_back(&Muon2StandardSelector);
  Extradist1d.push_back(&Muon3StandardSelector);

  Extradist1d.push_back(&Muon1isGlob);
  Extradist1d.push_back(&Muon2isGlob);
  Extradist1d.push_back(&Muon3isGlob);

  Extradist1d.push_back(&Muon1isTrack);
  Extradist1d.push_back(&Muon2isTrack);
  Extradist1d.push_back(&Muon3isTrack);

  Extradist1d.push_back(&TauEta);
  Extradist1d.push_back(&TauPt);
  Extradist1d.push_back(&TauP);
  Extradist1d.push_back(&TauMass);
  Extradist1d.push_back(&TauMassRefit);
  Extradist1d.push_back(&TauMassResolution);
  Extradist1d.push_back(&TauMassResolutionRefit);

  Extradist1d.push_back(&Muon1kink);
  Extradist1d.push_back(&Muon2kink);
  Extradist1d.push_back(&Muon3kink);
  Extradist1d.push_back(&VertexChi2KF);
  Extradist1d.push_back(&VertexChi2AF);

  Extradist1d.push_back(&Muon1InOutTrackMatch);
  Extradist1d.push_back(&Muon2InOutTrackMatch);
  Extradist1d.push_back(&Muon3InOutTrackMatch);

  Extradist1d.push_back(&Muon1PtResolution);
  Extradist1d.push_back(&Muon2PtResolution);
  Extradist1d.push_back(&Muon3PtResolution);

  Extradist1d.push_back(&Muon1EtaResolution);
  Extradist1d.push_back(&Muon2EtaResolution);
  Extradist1d.push_back(&Muon3EtaResolution);

  Extradist1d.push_back(&Muon1DRToTruth);
  Extradist1d.push_back(&Muon2DRToTruth);
  Extradist1d.push_back(&Muon3DRToTruth);

  Extradist1d.push_back(&MuPair1_vertex_chi2);
  Extradist1d.push_back(&MuPair2_vertex_chi2);
  Extradist1d.push_back(&MuPair3_vertex_chi2);

}


void  ThreeMu::doEvent(){ 
  unsigned int t;
  int id(Ntp->GetMCID());
  //  std::cout<<"id  "<<id<< std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection

  value.at(HLTOk) = 0;
  value.at(ThreeMuCandidate) = 0;
  if(Ntp->NThreeMuons()>0)  value.at(ThreeMuCandidate)=1;
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);
    if(HLT.Contains("DoubleMu3_Trk_Tau3mu") && Ntp->HLTDecision(iTrigger) == 1)value.at(HLTOk)=Ntp->HLTDecision(iTrigger);
  }
  pass.at(HLTOk)= (value.at(HLTOk)==cut.at(HLTOk)); 
  pass.at(ThreeMuCandidate)= (value.at(ThreeMuCandidate)==cut.at(ThreeMuCandidate)); 



  double wobs=1;
  double w;  
             
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}

  bool status=AnalysisCuts(t,w,wobs);
  if(status){

    unsigned int Muon_index_1=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(0);
    unsigned int Muon_index_2=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(1);
    unsigned int Muon_index_3=Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(2);

    TLorentzVector Muon1LV = Ntp->Muon_P4(Muon_index_1);
    TLorentzVector Muon2LV = Ntp->Muon_P4(Muon_index_2);
    TLorentzVector Muon3LV = Ntp->Muon_P4(Muon_index_3);
    
    TLorentzVector TauRefitLV = Ntp->Vertex_signal_KF_refittedTracksP4(0,0)+Ntp->Vertex_signal_KF_refittedTracksP4(0,1)+Ntp->Vertex_signal_KF_refittedTracksP4(0,2);

    Muon1Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_1).Pt(),1);
    Muon2Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_2).Pt(),1);
    Muon3Pt.at(t).Fill(Ntp->Muon_P4(Muon_index_3).Pt(),1);
    

    Muon1Eta.at(t).Fill( Ntp->Muon_P4(Muon_index_1).Eta(),1);
    Muon2Eta.at(t).Fill(  Ntp->Muon_P4(Muon_index_2).Eta(),1);
    Muon3Eta.at(t).Fill(  Ntp->Muon_P4(Muon_index_3).Eta(),1);

    TLorentzVector TauLV = Ntp->Muon_P4(Muon_index_1)  + Ntp->Muon_P4(Muon_index_2) + Ntp->Muon_P4(Muon_index_3);
    
    TauEta.at(t).Fill(TauLV.Eta(),1);
    TauPt.at(t).Fill(TauLV.Pt(),1);
    TauP.at(t).Fill(TauLV.P(),1);
    TauMass.at(t).Fill(TauLV.M(),1);
    TauMassRefit.at(t).Fill(TauRefitLV.M(),1);    
    for(unsigned int iMuSelector=0; iMuSelector< Ntp->MuonStandardSelectorBitMask(Muon_index_1).size(); iMuSelector++ ){
      if(Ntp->MuonStandardSelectorBitMask(Muon_index_1).at(iMuSelector)==1)  Muon1StandardSelector.at(t).Fill(iMuSelector,1);
    }

    for(unsigned int iMuSelector=0; iMuSelector< Ntp->MuonStandardSelectorBitMask(Muon_index_2).size(); iMuSelector++ ){
      if(Ntp->MuonStandardSelectorBitMask(Muon_index_2).at(iMuSelector)==1)  Muon2StandardSelector.at(t).Fill(iMuSelector,1);
    }

    for(unsigned int iMuSelector=0; iMuSelector< Ntp->MuonStandardSelectorBitMask(Muon_index_3).size(); iMuSelector++ ){
      if(Ntp->MuonStandardSelectorBitMask(Muon_index_3).at(iMuSelector)==1)  Muon3StandardSelector.at(t).Fill(iMuSelector,1);
    }

    Muon1isGlob.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_index_1),1);
    Muon2isGlob.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_index_2),1);
    Muon3isGlob.at(t).Fill(Ntp->Muon_isGlobalMuon(Muon_index_3),1);

    Muon1isTrack.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_index_1),1);
    Muon2isTrack.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_index_2),1);
    Muon3isTrack.at(t).Fill(Ntp->Muon_isTrackerMuon(Muon_index_3),1);

    Muon1kink.at(t).Fill(Ntp->Muon_combinedQuality_trkKink(Muon_index_1),1);
    Muon2kink.at(t).Fill(Ntp->Muon_combinedQuality_trkKink(Muon_index_2),1);
    Muon3kink.at(t).Fill(Ntp->Muon_combinedQuality_trkKink(Muon_index_3),1);
    VertexChi2KF.at(t).Fill(Ntp->ThreeMuons_SV_Chi2(0),1);
    VertexChi2AF.at(t).Fill(Ntp->Vertex_signal_AF_Chi2(0),1);


    Muon1InOutTrackMatch.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_1),1);
    Muon2InOutTrackMatch.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_2),1);
    Muon3InOutTrackMatch.at(t).Fill(Ntp->Muon_combinedQuality_chi2LocalPosition(Muon_index_3),1);

    MuPair1_vertex_chi2.at(t).Fill(Ntp->Vertex_pair_quality(0,0),1);
    MuPair2_vertex_chi2.at(t).Fill(Ntp->Vertex_pair_quality(0,1),1);
    MuPair3_vertex_chi2.at(t).Fill(Ntp->Vertex_pair_quality(0,2),1);

    //---------------  Fill MC plots 
    if(id==40 || id == 60 || id ==90){
      if(Ntp->MCEventIsReconstructed()){

      TLorentzVector MCMuon1LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(0)));
      TLorentzVector MCMuon2LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(1)));
      TLorentzVector MCMuon3LV= Ntp->matchToTruthTauDecay(Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(0)).at(2)));
      TLorentzVector MCTauLV  = MCMuon1LV+ MCMuon2LV + MCMuon3LV;
      Muon1PtResolution.at(t).Fill((Muon1LV.Pt() - MCMuon1LV.Pt())/MCMuon1LV.Pt(), 1);
      Muon2PtResolution.at(t).Fill((Muon2LV.Pt() - MCMuon2LV.Pt())/MCMuon2LV.Pt(), 1);
      Muon3PtResolution.at(t).Fill((Muon3LV.Pt() - MCMuon3LV.Pt())/MCMuon3LV.Pt(), 1);

      Muon1EtaResolution.at(t).Fill((Muon1LV.Eta() - MCMuon1LV.Eta())/MCMuon1LV.Eta(), 1);
      Muon2EtaResolution.at(t).Fill((Muon2LV.Eta() - MCMuon2LV.Eta())/MCMuon2LV.Eta(), 1);
      Muon3EtaResolution.at(t).Fill((Muon3LV.Eta() - MCMuon3LV.Eta())/MCMuon3LV.Eta(), 1);

      TauMassResolution.at(t).Fill((TauLV.M() - MCTauLV.M())/MCTauLV.M(),1);
      TauMassResolutionRefit.at(t).Fill((TauRefitLV.M() - MCTauLV.M())/MCTauLV.M(),1);

      Muon1DRToTruth.at(t).Fill(Muon1LV.DeltaR(MCMuon1LV),1);
      Muon2DRToTruth.at(t).Fill(Muon2LV.DeltaR(MCMuon2LV),1);
      Muon3DRToTruth.at(t).Fill(Muon3LV.DeltaR(MCMuon3LV),1);


      }
    }

  }
}



void  ThreeMu::Finish(){
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 

  if(mode == RECONSTRUCT){
    for(unsigned int i=1; i<  Nminus0.at(0).size(); i++){
      double scale(1.);
      if(Nminus0.at(0).at(i).Integral()!=0)scale = Nminus0.at(0).at(0).Integral()/Nminus0.at(0).at(i).Integral()/3;
      ScaleAllHistOfType(i,scale);
    }
  }


  Selection::Finish();

}





