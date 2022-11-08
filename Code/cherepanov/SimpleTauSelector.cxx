#include "SimpleTauSelector.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
//#include "Logger.h"

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
    if(i==SignalCandidate)  cut.at(SignalCandidate)=1;
    if(i==MuonCandidate)    cut.at(MuonCandidate)=1;
    if(i==OSCharge)         cut.at(OSCharge)=-1;
    if(i==nTaus)            cut.at(nTaus)=1;

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
      title.at(i)="N 3$\\mu$ candidates";
      htitle=title.at(i);
      hlabel="N $3\\mu$ candidates";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,5,-0.5,4.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }
    else if(i==MuonCandidate){
      title.at(i)=" $\\mu$ candidate";
      htitle=title.at(i);
      hlabel=" $3\\mu$ candidates";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==OSCharge){
      title.at(i)=" $\\mu$ and 3$\\mu$ are of opposite charge";
      htitle=title.at(i);
      hlabel="Opposite charge? ";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OSCharge_",htitle,3,-1.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OSCharge_",htitle,3,-1.5,1.5,hlabel,"Events"));
    }



  }
  // Setup NPassed Histogams

  NumberOfTaus=HConfig.GetTH1D(Name+"_NumberOfTaus","NumberOfTaus",5,-0.5,4.5,"Number of #tau ","Events");
  NumTausvsNumMuons=HConfig.GetTH2D(Name+"_NumTausvsNumMuons","NumTausvsNumMuons",5,-0.5,4.5,5,-0.5,4.5,"N 3#mu ","N #tau");

  Mu3MuVisibleMass=HConfig.GetTH1D(Name+"_Mu3MuVisibleMass","Mu3MuVisibleMass",60,0,120,"M_{#mu-3#mu}, GeV (Visible Mass)","Events");
  Mu3MudPhi=HConfig.GetTH1D(Name+"_Mu3MudPhi","Mu3MudPhi",30,-3.14,3.14,"#Delta #phi{#mu-3#mu}","Events");
  VertexChi2KF_vs_HelixFit=HConfig.GetTH2D(Name+"_VertexChi2KF_vs_HelixFit","VertexChi2KF_vs_HelixFit",50,0,10,50,0,10,"Kalman Vertex #chi^{2}","Helix Vertex  Fitter #chi^{2}");


  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms



  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  SimpleTauSelector::Store_ExtraDist(){ 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output

  Extradist1d.push_back(&NumberOfTaus);
  Extradist2d.push_back(&NumTausvsNumMuons);
  Extradist1d.push_back(&Mu3MuVisibleMass);
  Extradist1d.push_back(&Mu3MudPhi);
  Extradist2d.push_back(&VertexChi2KF_vs_HelixFit);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  SimpleTauSelector::doEvent(){ 

  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  // Apply Selection


  value.at(PrimeVtx)=Ntp->NVtx(); 
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx)); 
  
  value.at(TriggerOk)=(Ntp->EventNumber()%1000)==1;
  pass.at(TriggerOk)=true; 





  value.at(nTaus) = Ntp->NTaus();
  value.at(SignalCandidate) = Ntp->NThreeMuons();

  unsigned int  signal_idx=0;
  double min_chi2(99.);

  for(unsigned int i_idx =0; i_idx < Ntp->NThreeMuons(); i_idx++){
    if(Ntp->Vertex_Signal_KF_Chi2(i_idx) < min_chi2){
      min_chi2 = Ntp->Vertex_Signal_KF_Chi2(i_idx);
      signal_idx = i_idx;
    }
  }

  int muon_idx = -1;
  double dPhi(TMath::Pi()/2);


  pass.at(SignalCandidate) = (value.at(SignalCandidate) >= cut.at(SignalCandidate));
  pass.at(nTaus)  = true;//( value.at(nTaus) >= cut.at(nTaus) );
  value.at(OSCharge) = 1;
  if( pass.at(SignalCandidate) )
    {
      TLorentzVector  Tau3MuLV = Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(0))+
	Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(1))+
	Ntp->Muon_P4(Ntp->ThreeMuonIndices(signal_idx).at(2));


      for(unsigned int i =0; i < Ntp->NMuons(); i++)
	{

	  if(i == Ntp->ThreeMuonIndices(signal_idx).at(0) or 
	     i == Ntp->ThreeMuonIndices(signal_idx).at(1) or
	     i == Ntp->ThreeMuonIndices(signal_idx).at(2)) continue;




	  if( fabs(Ntp->DeltaPhi(Ntp->Muon_P4(i).Phi(), Tau3MuLV.Phi())) > dPhi)
	    {
	      dPhi = fabs(Ntp->DeltaPhi(Ntp->Muon_P4(i).Phi(),Tau3MuLV.Phi()));
	      muon_idx=i;
	    }
	}


    }

  if(muon_idx == -1) value.at(MuonCandidate) =0;
  if(muon_idx >=  0) value.at(MuonCandidate) =1;
  

  pass.at(MuonCandidate) = (value.at(MuonCandidate) == cut.at(MuonCandidate));

  if( pass.at(SignalCandidate)  && pass.at(MuonCandidate))
    {
      int Tau3MuCharge = Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)) +
	Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)) +
	Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2));
      
      value.at(OSCharge) = ( Tau3MuCharge*Ntp->Muon_charge(muon_idx)   );
    }
  pass.at(OSCharge) = (value.at(OSCharge) == cut.at(OSCharge));  

  double wobs=1;
  double w;  
             
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  

  
  bool status=AnalysisCuts(t,w,wobs);
  NumTausvsNumMuons.at(t).Fill(Ntp->NTaus(), Ntp->NThreeMuons());
  if(status){ 


    //      std::cout<<" Nsignal particles  "<< Ntp->NMCSignalParticles() << std::endl;

    //      for(unsigned int iMCsignalParticle = 0 ; iMCsignalParticle<  Ntp->NMCSignalParticles(); iMCsignalParticle++)
    //	{
    //	  std::cout<<" signal particle pdgid   "<< Ntp->MCSignalParticle_pdgid(iMCsignalParticle)<< " and daughtersr   "<< Ntp->MCSignalParticle_Nchilds(iMCsignalParticle)<< std::endl;
    //	}



    NumberOfTaus.at(t).Fill(Ntp->NTaus());


    TLorentzVector Tau3muLV = Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)) + 
      Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)) + 
      Ntp->Muon_P4(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2));

    int Tau3MuCharge = Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(0)) +
      Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(1)) +
      Ntp->Muon_charge(Ntp->SortedPtMuons(Ntp->ThreeMuonIndices(signal_idx)).at(2));

    //    std::cout<<"   n truth taus   : "<< Ntp->NMCTaus() << std::endl;

    //    for(unsigned int itau = 0 ; itau<  Ntp->NTaus(); itau++){
      //      std::cout<<" ------------------- "<< std::endl;
    //      for(unsigned int idautau =0; idautau< Ntp->NMCTauDecayProducts(itau); idautau++){
	//	std::cout<<"   daus pdgid   "<< Ntp->MCTauandProd_pdgid(itau,idautau) << std::endl;
    //	Ntp->MCTauandProd_p4(itau,idautau).Print();
    //      }
    //    }
    //    std::cout<<" recoed muon "<< std::endl;




    TLorentzVector MuLV = Ntp->Muon_P4(muon_idx);
    //    MuLV.Print();
    //    TrackParticle MuTrP = Ntp->Muon_TrackParticle(muon_idx);
    //    MuTrP.getCovMatrix().Print();
    //  bool Muon_TrackParticleHasMomentum(unsigned int i){if(Ntp->Muon_par->at(i).size()!=0)return true; return false;}
    // TrackParticle Muon_TrackParticle(unsigned int i){

    TrackParticle MuonTP = Ntp->Muon_TrackParticle(muon_idx);


    //    std::cout<<"   dxy  "<<     MuonTP.Parameter(TrackParticle::dxy)<< std::endl;
    //    std::cout<<"   phi  "<<     MuonTP.Parameter(TrackParticle::phi)<< std::endl;
    //    std::cout<<"   lambda "<<     MuonTP.Parameter(TrackParticle::lambda)<< std::endl;
    ///    std::cout<<"   dz "<<     MuonTP.Parameter(TrackParticle::dz)<< std::endl;
    //    std::cout<<"   kappa "<<     MuonTP.Parameter(TrackParticle::kappa)<< std::endl;
    //    MuonTP.getParMatrix();

    LorentzVectorParticle Tau3MuLVP = Ntp->Tau3mu_LVP(  signal_idx );
    std::cout<<"  LVP from the ntuple    " << std::endl;
    Tau3MuLVP.LV().Print();


    unsigned int nMetComponents = 2;
    TMatrixT<double> metVector(nMetComponents, 1);
    TMatrixTSym<double> metCovariance; metCovariance.ResizeTo(nMetComponents, nMetComponents);

    metVector[0][0] = Ntp->METEt() * cos( Ntp->METPhi());
    metVector[1][0] = Ntp->METEt() * sin( Ntp->METPhi());

    metCovariance[0][0]= Ntp->METXX();
    metCovariance[1][1]= Ntp->METYY();
    metCovariance[0][1]= Ntp->METXY();
    metCovariance[1][0]= Ntp->METXY();


    PTObject metInput(metVector, metCovariance);

    TVector3 pvInput = Ntp->Vertex_HighestPt_PrimaryVertex();
    TMatrixTSym<double> pvCovarianceInput = Ntp->Vertex_HighestPt_PrimaryVertex_Covariance();

    //    void SetLevel(level _l){l=_l;}
    //    Logger::Instance()->SetLevel(Logger::Warning);

    //    std::cout<<" LVP vertex "<<std::endl;    Tau3MuLVP.Vertex().Print();
    //    std::cout<<"   ntuple SV "<< std::endl;  Ntp->Vertex_Signal_KF_pos(signal_idx).Print();
    
    GlobalEventFit* globalEventFit = nullptr;
    globalEventFit = new GlobalEventFit(MuonTP,  Tau3MuLVP , metInput, pvInput, pvCovarianceInput, true);


    globalEventFit->setMinimizer("Default");
    TPTRObject tauReco = globalEventFit->getTPTRObject();
    GEFObject  fitResult = globalEventFit->Fit();


  std:cout<<"   Tau3MuLVP  "<< Tau3MuLVP.LV().Px() << "  "<<Tau3MuLVP.LV().Py() << "  "<<Tau3MuLVP.LV().Pz() <<std::endl;

    std::cout<<"========================================  Refit Particles  "<< std::endl;
    std::cout<<"==";    fitResult.getTauH().LV().Print();
    std::cout<<"==";    fitResult.getTauMu().LV().Print();
    std::cout<<"========================================  Refit Particles  "<< std::endl;
    std::cout<<"========================================  Initial Particles  "<< std::endl;
    std::cout<<"==";    fitResult.getInitTauH().LV().Print();
    std::cout<<"==";    fitResult.getInitTauMu().LV().Print();
    std::cout<<"========================================  Initial Particles  "<< std::endl;
    //    std::cout<<"  ZTT3MuSimpleFitProducer     fit result  "<< fitResult.isValid()<< std::endl;
    //    std::cout<<"  csum, chi2, nIterations   "<< fitResult.getCsum() << "  "<< fitResult.getChi2() << "   "<< fitResult.getNiterations() << std::endl;
    //    std::cout<<"  Taus Pt  "<< fitResult.getTaus().at(0).LV().Pt() << "  resonance mass  " <<fitResult.getResonance().LV().M() << std::endl;
    





    //    metVector.Print();
    //    metCovariance.Print();
    /*
    TVector3 vguess(0.1,0.1,0.1);
    std::vector<TrackParticle> MuonsTrackParticles;
    MuonsTrackParticles.push_back(    Ntp->Muon_TrackParticle(Ntp->ThreeMuonIndices(signal_idx).at(0)) );
    MuonsTrackParticles.push_back(    Ntp->Muon_TrackParticle(Ntp->ThreeMuonIndices(signal_idx).at(1)) );
    MuonsTrackParticles.push_back(    Ntp->Muon_TrackParticle(Ntp->ThreeMuonIndices(signal_idx).at(2)) );


    Chi2VertexFitter  Fitter(MuonsTrackParticles,vguess);
   
    Fitter.Fit();
    VertexChi2KF_vs_HelixFit.at(t).Fill(Ntp->Vertex_signal_KF_Chi2(signal_idx),Fitter.ChiSquare());
    std::vector<LorentzVectorParticle> ReffitedLVParticles = Fitter.GetReFitLorentzVectorParticles();

    TLorentzVector ReffitedTau = Fitter.GetMother(444).LV();
      
    */

    //  if()


    Mu3MuVisibleMass.at(t).Fill( (Tau3muLV+MuLV).M() , 1.);
    Mu3MudPhi.at(t).Fill( Ntp->DeltaPhi(Tau3muLV.Phi(), MuLV.Phi())  , 1.);
    //    std::cout<<"  masses  "<< Ntp->Tau3mu_LVP(signal_idx).LV().M() << "  "<< Tau3muLV.M() << std::endl;
    //    Ntp->Tau3mu_LVP(signal_idx).LV().Print();

    //    Ntp->Tau3mu_LVP(signal_idx).VertexCov().Print();
    //    Tau3muLV.Print();

    //    Ntp->Tau3mu_LVP(signal_idx).LVCov().Print();
    //    Ntp->Tau3mu_LVP(signal_idx).Vertex().Print();

    //    std::cout<<"  signal KF   Cov  "<< std::endl;
    //    Ntp->Vertex_Signal_KF_Covariance(signal_idx).Print();

    //    Ntp->Vertex_Signal_KF_pos(signal_idx).Print();

  }
}


void  SimpleTauSelector::Finish(){

  Selection::Finish();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
}





