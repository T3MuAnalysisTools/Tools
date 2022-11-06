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
  
  std::cout<<"  N eles "<< Ntp->NElectrons() <<std::endl;
  for(unsigned int ie=0; ie < Ntp->NElectrons(); ie++)
    {
      std::cout<<      Ntp->Electron_P4(ie).Pt()  << std::endl;
	//      std::cout<<"  charge  "<< Ntp->Electron_charge(ie) << " trackiso  "<< Ntp->Electron_puppiPhotonIso(ie) << "  medium  " << Ntp->Electron_cutBasedElectronID_Fall17_94X_V2_medium(ie) << std::endl;
    }

  
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

    for(unsigned int itau = 0 ; itau<  Ntp->NTaus(); itau++){
      //      std::cout<<" ------------------- "<< std::endl;
      for(unsigned int idautau =0; idautau< Ntp->NMCTauDecayProducts(itau); idautau++){
	//	std::cout<<"   daus pdgid   "<< Ntp->MCTauandProd_pdgid(itau,idautau) << std::endl;
	Ntp->MCTauandProd_p4(itau,idautau).Print();
      }
    }
    //    std::cout<<" recoed muon "<< std::endl;




    TLorentzVector MuLV = Ntp->Muon_P4(muon_idx);
    //    MuLV.Print();
    //    TrackParticle MuTrP = Ntp->Muon_TrackParticle(muon_idx);
    //    MuTrP.getCovMatrix().Print();
    //  bool Muon_TrackParticleHasMomentum(unsigned int i){if(Ntp->Muon_par->at(i).size()!=0)return true; return false;}
    // TrackParticle Muon_TrackParticle(unsigned int i){



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





