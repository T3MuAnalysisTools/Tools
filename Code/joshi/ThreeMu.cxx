#include "ThreeMu.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

double ThreeMu::deltaR(double eta1, double phi1, double eta2, double phi2)
{
    double deta = eta1 - eta2;
    double dphi = phi1 - phi2;
    while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
    while (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();

    return sqrt(deta*deta + dphi*dphi);
}
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
    if(i==L1SeedOk)     cut.at(L1SeedOk)=1;
    if(i==HLTOk)        cut.at(HLTOk)=1;
    if(i==isThreeMu)        cut.at(isThreeMu)=1;
    if(i==PrimeVtx)     cut.at(PrimeVtx)=5; // Here for example we place cut value on number of PVs
	 if(i==fitVtxChiSq)	cut.at(fitVtxChiSq)=5.0;
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

    else if(i==isThreeMu){
      title.at(i)="isThreeMu ";
      hlabel="3mu category";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_isThreeMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_isThreeMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }

    else if(i==HLTOk){
      title.at(i)="HLT trigger ";
      hlabel="DoubleMu3_Trk_Tau3mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
	 
	 else if(i==fitVtxChiSq){
      title.at(i)="Normalized chi sq fit vertex $(<5)$ ";
      hlabel="Normalized chi sq fit vertex";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_fitVtxChiSq_",htitle,100,-0.5,19.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_fitVtxChiSq_",htitle,100,-0.5,19.5,hlabel,"Events"));
    }
  } 

// Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms
  // Book here your analysis histrogramms, a good style is to follow selfexplanatory convention
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",66,-0.5,65.5,"Number of Vertices","Events");
  //Muon variables (Muons from dimuon + track candidates)
  Isolation_NTracks=HConfig.GetTH1D(Name+"_Isolation_NTracks","",20,0,10,"Iso ntrks","Events");
  Isolation_RelPt=HConfig.GetTH1D(Name+"_Isolation_RelPt","",20,0,10,"Iso rel p_{T}","Events");
  Isolation_MinDist=HConfig.GetTH1D(Name+"_Isolation_MinDist","",10,0,1,"Iso MinDist","Events");
  Isolation05_RelPt=HConfig.GetTH1D(Name+"_Isolation05_RelPt","",10,0,3,"Iso05 rel p_{T}","Events");
  Isolation05_NTracks=HConfig.GetTH1D(Name+"_Isolation05_NTracks","",20,0,10,"Iso05 ntrks","Events");
  Isolation05_MinDist=HConfig.GetTH1D(Name+"_Isolation05_MinDist","",10,0,1,"Iso05 MinDist","Events");
  Isolation_Ntrk1=HConfig.GetTH1D(Name+"_Isolation_Ntrk1","",10,0,10,"Iso ntrk 1","Events");
  Isolation_Ntrk2=HConfig.GetTH1D(Name+"_Isolation_Ntrk2","",10,0,10,"Iso ntrk 2","Events");
  Isolation_Ntrk3=HConfig.GetTH1D(Name+"_Isolation_Ntrk3","",10,0,10,"Iso ntrk 3","Events");
  Isolation_Ntrk0p1=HConfig.GetTH1D(Name+"_Isolation_Ntrk0p1","",10,0,10,"Iso ntrk0p1","Events");
  Isolation_Ntrk0p2=HConfig.GetTH1D(Name+"_Isolation_Ntrk0p2","",10,0,10,"Iso ntrk0p2","Events");
  Isolation_Ntrk0p5=HConfig.GetTH1D(Name+"_Isolation_Ntrk0p5","",10,0,10,"Iso ntrk0p5","Events");
  Isolation_maxdxy=HConfig.GetTH1D(Name+"_Isolation_maxdxy","",40,0,20,"Iso max(dxy)","Events");
  
  //Dimuon Information (Muons from dimuon + track candidates)
  Muon1Muon3dR=HConfig.GetTH1D(Name+"_Muon1Muon3dR","dR between the highest p muon and the lowest pt muon",100,0,5,"dR","Events");
  Muon2Muon3dR=HConfig.GetTH1D(Name+"_Muon2Muon3dR","dR between the second highest p muon and the lowest pt muon",100,0,5,"dR","Events");
  Muon1Muon2dR=HConfig.GetTH1D(Name+"_Muon1Muon1dR","dR between the highest p muon and the second highest p muon",100,0,5,"dR","Events");
  TripleMass=HConfig.GetTH1D(Name+"_TripleMass","#mu#mu#mu mass",50,1.5,2.0,"Mass of the #mu#mu#mu","Events");

Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}

void  ThreeMu::Store_ExtraDist(){ 
  Extradist1d.push_back(&Muon1Muon3dR);
  Extradist1d.push_back(&Muon2Muon3dR);
  Extradist1d.push_back(&Muon1Muon2dR);
  Extradist1d.push_back(&TripleMass);

  Extradist1d.push_back(&Isolation_NTracks);
  Extradist1d.push_back(&Isolation_RelPt);
  Extradist1d.push_back(&Isolation_MinDist);
  Extradist1d.push_back(&Isolation05_RelPt);
  Extradist1d.push_back(&Isolation05_NTracks);
  Extradist1d.push_back(&Isolation05_MinDist);
  Extradist1d.push_back(&Isolation_Ntrk1);
  Extradist1d.push_back(&Isolation_Ntrk2);
  Extradist1d.push_back(&Isolation_Ntrk3);
  Extradist1d.push_back(&Isolation_Ntrk0p1);
  Extradist1d.push_back(&Isolation_Ntrk0p2);
  Extradist1d.push_back(&Isolation_Ntrk0p5);
  Extradist1d.push_back(&Isolation_maxdxy);
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output //
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  Extradist1d.push_back(&NVtx);
}

void  ThreeMu::doEvent(){ 
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}
  
  // Apply Selection
  value.at(L1SeedOk) = 0;
  value.at(HLTOk) = 0;
  for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
    TString HLT = Ntp->HLTName(iTrigger);
    if(( HLT.Contains("DoubleMu3_Trk_Tau3mu") || HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu_v2")) && Ntp->HLTDecision(iTrigger) == 1)value.at(HLTOk)=Ntp->HLTDecision(iTrigger);
  }
  for(int l1iTrigger=0; l1iTrigger < Ntp->NL1Seeds(); l1iTrigger++){
    TString L1 = Ntp->L1Name(l1iTrigger);
    if(L1.Contains("L1_DoubleMu0er1p4_dEta_Max1p8_OS") && Ntp->L1Decision(l1iTrigger) == 1)value.at(L1SeedOk)=Ntp->L1Decision(l1iTrigger);
    if(L1.Contains("L1_DoubleMu_10_0_dEta_Max1p8") && Ntp->L1Decision(l1iTrigger) == 1)value.at(L1SeedOk)=Ntp->L1Decision(l1iTrigger);
    if(L1.Contains("L1_TripleMu0") && Ntp->L1Decision(l1iTrigger) == 1)value.at(L1SeedOk)=Ntp->L1Decision(l1iTrigger);
  }

  value.at(PrimeVtx)=Ntp->NVtx();
  value.at(fitVtxChiSq)=0;

  value.at(isThreeMu) = 0;
  if(Ntp->NTwoMuonsTrack()==0 && Ntp->NThreeMuons() != 0) value.at(isThreeMu) = 1;

  pass.at(isThreeMu) = (value.at(isThreeMu) == cut.at(isThreeMu));
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx)); 
  pass.at(L1SeedOk)= (value.at(L1SeedOk)==cut.at(L1SeedOk)); 
  pass.at(HLTOk)= (value.at(HLTOk)==cut.at(HLTOk));

  int mu1=-1, mu2=-1, mu3=-1;
  int tmp_idx = -1;
  double tmp_chisq = 999;

  if (value.at(isThreeMu)==1){
  	for (unsigned int i3M=0; i3M<Ntp->NThreeMuons(); i3M++){
		int tmp_mu1 =  Ntp->ThreeMuonIndices(i3M).at(0);
      int tmp_mu2 =  Ntp->ThreeMuonIndices(i3M).at(1);
      int tmp_mu3 =  Ntp->ThreeMuonIndices(i3M).at(2);
		if (tmp_chisq>Ntp->ThreeMuons_SV_Chi2(i3M)){
			tmp_idx = i3M;
			tmp_chisq = Ntp->ThreeMuons_SV_Chi2(i3M);
			mu1 = tmp_mu1;
			mu2 = tmp_mu2;
			mu3 = tmp_mu3;
		}
	}
  }
  if (tmp_idx==-1) value.at(fitVtxChiSq)=999.0;
  else value.at(fitVtxChiSq) = tmp_chisq;
  pass.at(fitVtxChiSq) = (value.at(fitVtxChiSq)<cut.at(fitVtxChiSq));
  
  double wobs=1;
  double w;  
             
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}

  bool status=AnalysisCuts(t,w,wobs);
  if(status){
  if (DEBUG) cout<<mu1<<" "<<mu2<<" "<<mu3<<" "<<Ntp->NMuons()<<endl;
	 NVtx.at(t).Fill(Ntp->NVtx(),w);

	 
	 Isolation_NTracks.at(t).Fill(Ntp->Isolation_NTracks(tmp_idx),w);
    Isolation_RelPt.at(t).Fill(Ntp->Isolation_RelPt(tmp_idx),w);
    Isolation_MinDist.at(t).Fill(Ntp->Isolation_MinDist(tmp_idx),w);
    Isolation05_RelPt.at(t).Fill(Ntp->Isolation05_RelPt(tmp_idx),w);
    Isolation05_NTracks.at(t).Fill(Ntp->Isolation05_NTracks(tmp_idx),w);
    Isolation05_MinDist.at(t).Fill(Ntp->Isolation05_MinDist(tmp_idx),w);
    Isolation_Ntrk1.at(t).Fill(Ntp->Isolation_Ntrk1(tmp_idx),w);
    Isolation_Ntrk2.at(t).Fill(Ntp->Isolation_Ntrk2(tmp_idx),w);
    Isolation_Ntrk3.at(t).Fill(Ntp->Isolation_Ntrk3(tmp_idx),w);
    Isolation_Ntrk0p1.at(t).Fill(Ntp->Isolation_Ntrk0p1(tmp_idx),w);
    Isolation_Ntrk0p2.at(t).Fill(Ntp->Isolation_Ntrk0p2(tmp_idx),w);
    Isolation_Ntrk0p5.at(t).Fill(Ntp->Isolation_Ntrk0p5(tmp_idx),w);
    Isolation_maxdxy.at(t).Fill(Ntp->Isolation_maxdy(tmp_idx),w);
	 
	 TripleMass.at(t).Fill((Ntp->Muon_P4(mu1)+Ntp->Muon_P4(mu2)+Ntp->Muon_P4(mu3)).M(), w);
  }
}


void  ThreeMu::Finish(){
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if(mode == RECONSTRUCT){
    for(unsigned int i=0; i<  Nminus0.at(0).size(); i++){
      double scale(1.);
      if(Nminus0.at(0).at(i).Integral()!=0)scale = 1/Nminus0.at(0).at(i).Integral();
      ScaleAllHistOfType(i,scale);
    }
  }


  Selection::Finish();

}
