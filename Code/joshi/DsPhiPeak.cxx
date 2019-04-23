#include "MCEfficiency.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
#include "Logger.h"

template<typename T>
void printVec(int size, T& vec){
    for (int i = 0; i<size; i++) cout<<vec[i]<<" ";
    cout<<endl;
}

double MCEfficiency::deltaR(double eta1, double phi1, double eta2, double phi2)
{
    double deta = eta1 - eta2;
    double dphi = phi1 - phi2;
    while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
    while (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();

    return sqrt(deta*deta + dphi*dphi);
}

MCEfficiency::MCEfficiency(TString Name_, TString id_):
    Selection(Name_,id_)
{


    // This is a class constructor;
}

MCEfficiency::~MCEfficiency(){
    for(unsigned int j=0; j<Npassed.size(); j++){
      Logger(Logger::Info) << "Selection Summary before: "
        << Npassed.at(j).GetBinContent(1)  << " +/- " << Npassed.at(j).GetBinError(1)  << " after: "
        << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
    }


    Logger(Logger::Info) << "complete." << std::endl;
}

void  MCEfficiency::Configure(){

    for(int i=0; i<NCuts;i++){
      cut.push_back(0);
      value.push_back(0);
      pass.push_back(false);
      if(i==L1SeedOk)  cut.at(L1SeedOk)=1;
      if(i==HLTOk)  cut.at(HLTOk)=1;
      if(i==is2MuTrk)   cut.at(is2MuTrk)=1;
      if(i==PrimeVtx)   cut.at(PrimeVtx)=5; // Here for example we place cut value on number of PVs
      if(i==trkPt)   cut.at(trkPt)=1.2;
      if(i==nTrkHits)  cut.at(nTrkHits)=12;
      if(i==mumuMass)  cut.at(mumuMass)=PDG_Var::Phi_mass();
      if(i==fitVtxChiSq) cut.at(fitVtxChiSq)=5.0;
		if(i==trigObjMatch) cut.at(trigObjMatch)=1;
		if(i==genMatch) cut.at(genMatch)=1;
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

      else if(i==is2MuTrk){
        title.at(i)="Catrgory: 2Mu+Trk ";
        hlabel="2muon + track category";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_mumuMass_",htitle,2,-0.5,1.5,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_mumuMass_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }

      else if(i==mumuMass){
        title.at(i)="MuMu Invariant Mass within (1.0,1.04) GeV";
        hlabel="mumuMass";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_mumuMass_",htitle,80,0.8,1.6,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_mumuMass_",htitle,80,0.8,1.6,hlabel,"Events"));
      }

      else if(i==HLTOk){
        title.at(i)="HLT trigger ";
        hlabel="DoubleMu3_Trk_Tau3mu";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HLTOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }

      else if(i==trkPt){
        title.at(i)="Track pt $(pt>$";
        title.at(i)+=cut.at(trkPt);
        title.at(i)+=")";
        hlabel="pT of the track";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_trkPt_",htitle,100,-0.5,9.5,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_trkPt_",htitle,100,-0.5,9.5,hlabel,"Events"));
      }

      else if(i==nTrkHits){
        title.at(i)="Number of hits in the tracker $(N>$";
        title.at(i)+=cut.at(nTrkHits);
        title.at(i)+=")";
        hlabel="Number of hits in the tracker";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nTrkHits_",htitle,11,-1,10,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_NMinus0_nTrkHits_",htitle,11,-1,10,hlabel,"Events"));
      }

      else if(i==fitVtxChiSq){
        title.at(i)="Normalized chi sq of the vertex fit $(<$";
        title.at(i)+=cut.at(fitVtxChiSq);
        title.at(i)+=")";
        hlabel="Normalized chi square of the vertex fit";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_fitVtxChiSq_",htitle,100,-0.5,9.5,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_NMinus0_fitVtxChiSq_",htitle,100,-0.5,9.5,hlabel,"Events"));
      }
		
		else if(i==trigObjMatch){
        title.at(i)="2mu+trk matched to trigger object";
        title.at(i)+=cut.at(trigObjMatch);
        hlabel="Trigger Object Match";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_trigObjMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_NMinus0_trigObjMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
		
		else if(i==genMatch){
        title.at(i)="2mu+trk matched to GEN particles";
        title.at(i)+=cut.at(genMatch);
        hlabel="GEN particle match";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_genMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_NMinus0_genMatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
    }

    // Track Candidate Information
    Track_P=HConfig.GetTH1D(Name+"_Track_P","Momentum magnitude of track (2mu+trk track candidate)",66,-0.5,65.5,"p (track)","Events");
    Track_E=HConfig.GetTH1D(Name+"_Track_E","Energy of track (2mu+trk track candidate)",66,-0.5,65.5,"E (track)","Events");
    Track_Pt=HConfig.GetTH1D(Name+"_Track_Pt","Transverse momentum of track (2mu+trk track candidate)",66,-0.5,65.5,"p_{T} (track)","Events");
    Track_Eta=HConfig.GetTH1D(Name+"_Track_Eta","Psuedorapidity of track (2mu+trk track candidate)",66,-0.5,65.5,"#eta","Events");
    Track_Phi=HConfig.GetTH1D(Name+"_Track_Phi","Azimuthal angle of track (2mu+trk track candidate)",66,-0.5,65.5,"#phi","Events");
    Track_vx=HConfig.GetTH1D(Name+"_Track_vx","X coordinate of the parent vertex (2mu+trk track candidate)",66,-0.5,65.5,"Parent vertex x coordinate (cm)","Events");
    Track_vy=HConfig.GetTH1D(Name+"_Track_vy","Y coordinate of the parent vertex (2mu+trk track candidate)",66,-0.5,65.5,"Parent vertex y coordinate (cm)","Events");
    Track_vz=HConfig.GetTH1D(Name+"_Track_vz","Z coordinate of the parent vertex (2mu+trk track candidate)",66,-0.5,65.5,"Parent vertex z coordinate (cm)","Events");
    Track_normalizedChi2=HConfig.GetTH1D(Name+"_Track_normalizedChi2","Normalized chi square",66,-0.5,65.5,"#chi^{2}/ndf (track fit)","Events");
    Track_numberOfValidHits=HConfig.GetTH1D(Name+"_Track_numberOfValidHits","number of valid hits in te tracker",66,-0.5,65.5,"n valid track hits","Events");
    Track_charge=HConfig.GetTH1D(Name+"_Track_charge","Chargeof the track",66,-0.5,65.5,"Number of Vertices","Events");
    Track_dxy=HConfig.GetTH1D(Name+"_Track_dxy","Transverse displacement of the parent vertex from the bs",66,-0.5,65.5,"dxy (cm)","Events");
    Track_dz=HConfig.GetTH1D(Name+"_Track_dz","Longitudnal displacement of the parent vertex from the bs",66,-0.5,65.5,"dz (cm)","Events");
    Track_dxyError=HConfig.GetTH1D(Name+"_Track_dxyError","dxy Error",66,-0.5,65.5,"#Deltadxy (cm)","Events");
    Track_dzError=HConfig.GetTH1D(Name+"_Track_dzError","dz Error",66,-0.5,65.5,"#Deltadz (cm)","Events");

    // Muon variables (Muons from dimuon + track candidates)
    Muon1_Pt=HConfig.GetTH1D(Name+"_Muon1_Pt","Transverse Pt (muon 1)",25,0,50,"#mu_{1} p_{T} (GeV)","Events");
    Muon1_Eta=HConfig.GetTH1D(Name+"_Muon1_Eta","Psuedorapidity (muon 1)",25,-2.5,2.5,"#mu_{1} #eta","Events");
    Muon1_Phi=HConfig.GetTH1D(Name+"_Muon1_Phi","Azimuthal angle of (muons 1)",25,-3.15,3.15,"#mu_{1} #phi","Events"); 
    Muon1_E=HConfig.GetTH1D(Name+"_Muon1_E","Energy of all (muon 1)",20,0,40,"#mu_{1} E (GeV)","Events");
    Muon1_P=HConfig.GetTH1D(Name+"_Muon1_P","Magnitude of momentum of (muon 1)",20,0,40,"#mu_{1} p (GeV)","Events");  

    Muon1_vx=HConfig.GetTH1D(Name+"_Muon1_Vx","X coordinate of the parent vertex all muons",100,0,5,"#mu_{1} vx","Events"); 
    Muon1_vy=HConfig.GetTH1D(Name+"_Muon1_Vy","Y coordinate of the parent vertex all muons",100,0,5,"#mu_{1} vy","Events"); 
    Muon1_vz=HConfig.GetTH1D(Name+"_Muon1_Vz","Z coordinate of the parent vertex all muons",100,0,5,"#mu_{1} vz","Events");

    Muon2_Pt=HConfig.GetTH1D(Name+"_Muon2_Pt","Transverse Pt (muon 2)",25,0,50,"#mu_{2} p_{T} (GeV)","Events");
    Muon2_Eta=HConfig.GetTH1D(Name+"_Muon2_Eta","Psuedorapidity (muon 2)",25,-2.5,2.5,"#mu_{2} #eta","Events");
    Muon2_Phi=HConfig.GetTH1D(Name+"_Muon2_Phi","Azimuthal angle of (muons 1)",25,-3.15,3.15,"#mu_{2} #phi","Events"); 
    Muon2_E=HConfig.GetTH1D(Name+"_Muon2_E","Energy of all (muon 2)",20,0,40,"#mu_{2} E (GeV)","Events");
    Muon2_P=HConfig.GetTH1D(Name+"_Muon2_P","Magnitude of momentum of (muon 2)",20,0,40,"#mu_{2} p (GeV)","Events");  
    Muon2_vx=HConfig.GetTH1D(Name+"_Muon2_Vx","X coordinate of the parent vertex all muons",100,0,5,"#mu_{2} vx","Events"); 
    Muon2_vy=HConfig.GetTH1D(Name+"_Muon2_Vy","Y coordinate of the parent vertex all muons",100,0,5,"#mu_{2} vy","Events"); 
    Muon2_vz=HConfig.GetTH1D(Name+"_Muon2_Vz","Z coordinate of the parent vertex all muons",100,0,5,"#mu_{2} vz","Events");

    Muon1_isGlobal=HConfig.GetTH1D(Name+"_Muon1_isGlobal","Global muons status ",2,-.5,1.5,"#mu_{1} isGlb","Events");
    Muon2_isGlobal=HConfig.GetTH1D(Name+"_Muon2_isGlobal","",2,-0.5,0.5,"#mu_{2} isGlb","Events");
    Muon1_isStandAlone=HConfig.GetTH1D(Name+"_Muon1_isStandAlone","",2,-0.5,1.5,"#mu_{1} isStandAlone","Events");
    Muon2_isStandAlone=HConfig.GetTH1D(Name+"_Muon2_isStandAlone","",2,-0.5,1.5,"#mu_{2} isStandAlone","Events");
    Muon1_isTracker=HConfig.GetTH1D(Name+"_Muon1_isTracker","",2,-0.5,1.5,"#mu_{1} isTracker","Events");
    Muon2_isTracker=HConfig.GetTH1D(Name+"_Muon2_isTracker","",2,-0.5,1.5,"#mu_{2} isTracker","Events");
    Muon1_isCalo=HConfig.GetTH1D(Name+"_Muon1_isCaloMuon","",2,-0.5,1.5,"#mu_{1} isCalo","Events");
    Muon2_isCalo=HConfig.GetTH1D(Name+"_Muon2_isCaloMuon","",2,-0.5,1.5,"#mu_{2} isCalo","Events");
    Muon1_isIsolationValid=HConfig.GetTH1D(Name+"_Muon1_isIsolationValid","#mu_{1} isIsoValid",2,-0.5,1.5,"#mu_{1} isIsolationValid","Events");
    Muon2_isIsolationValid=HConfig.GetTH1D(Name+"_Muon2_isIsolationValid","#mu_{2} isIsoValid",2,-0.5,1.5,"#mu_{2} isIsolationValid","Events");
    Track_TriggerMatchdR=HConfig.GetTH1D(Name+"_Track_TriggerMatchdR","track dR (trigger match)",10,-0.5,9.5,"track dR (trigger match)","Events");
    Muon1_TriggerMatchdR=HConfig.GetTH1D(Name+"_Muon1_TriggerMatchdR","#mu_{1} dR (trigger match)",10,-0.5,9.5,"#mu_{1} dR (trigger match)","Events");
    Muon2_TriggerMatchdR=HConfig.GetTH1D(Name+"_Muon2_TriggerMatchdR","#mu_{2} dR (trigger match)",10,-0.5,9.5,"#mu_{2} dR (trigger match)","Events");

    //Dimuon Information (Muons from dimuon + track candidates)
    DimuondR=HConfig.GetTH1D(Name+"_DimuondR","dR between the muon pair",20,0,1,"dR","Events");
    Muon1TrkdR=HConfig.GetTH1D(Name+"_Muon1TrkdR","dR between the highest p muon and the track",100,0,5,"dR (#mu_{1},track)","Events");
    Muon2TrkdR=HConfig.GetTH1D(Name+"_Muon2TrkdR","dR between the lowest p muon and the track",100,0,5,"dR (#mu_{2},track)","Events");
    PhiMass=HConfig.GetTH1D(Name+"_PhiMass","#mu#mu mass",50,0.2,1.5,"Invariant mass of the #mu#mu pair","Events");
    TripleMass=HConfig.GetTH1D(Name+"_TripleMass","#mu#mu + track mass",50,1.7,2.1,"Invariant mass of the #mu#mu + track","Events");
    PhiMassVsDsMass=HConfig.GetTH2D(Name+"_PhiMassVsDsMass","#mu#mu Invariant mass vs. #mu#mu + track Invariant mass",50,0.2,1.5,50,1.7,2.1,"M_{#mu#mu}, GeV","M_{#mu#mu + track}, GeV");

    // Setup NPassed Histogams
    Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
    // Setup Extra Histograms
    // Book here your analysis histrogramms, a good style is to follow selfexplanatory convention
    NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",66,-0.5,65.5,"Number of Vertices","Events");

    Selection::ConfigureHistograms(); //do not remove
    HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}

void  MCEfficiency::Store_ExtraDist(){ 
    //Track candidate variables 
    Extradist1d.push_back(&Track_P);
    Extradist1d.push_back(&Track_E);
    Extradist1d.push_back(&Track_Pt);
    Extradist1d.push_back(&Track_Eta);
    Extradist1d.push_back(&Track_Phi);
    Extradist1d.push_back(&Track_vx);
    Extradist1d.push_back(&Track_vy);
    Extradist1d.push_back(&Track_vz);
    Extradist1d.push_back(&Track_normalizedChi2);
    Extradist1d.push_back(&Track_numberOfValidHits);
    Extradist1d.push_back(&Track_charge);
    Extradist1d.push_back(&Track_dxy);
    Extradist1d.push_back(&Track_dz);
    Extradist1d.push_back(&Track_dxyError);
    Extradist1d.push_back(&Track_dzError);

    //Dimuon variable
    Extradist1d.push_back(&Muon1_P);
    Extradist1d.push_back(&Muon1_E);
    Extradist1d.push_back(&Muon1_Pt);
    Extradist1d.push_back(&Muon1_Phi);
    Extradist1d.push_back(&Muon1_Eta);
    Extradist1d.push_back(&Muon1_vx);
    Extradist1d.push_back(&Muon1_vy);
    Extradist1d.push_back(&Muon1_vz);
    Extradist1d.push_back(&Muon2_P);
    Extradist1d.push_back(&Muon2_E);
    Extradist1d.push_back(&Muon2_Pt);
    Extradist1d.push_back(&Muon2_Phi);
    Extradist1d.push_back(&Muon2_Eta);
    Extradist1d.push_back(&Muon2_vx);
    Extradist1d.push_back(&Muon2_vy);
    Extradist1d.push_back(&Muon2_vz);
    Extradist1d.push_back(&DimuondR);
    Extradist1d.push_back(&Muon1TrkdR);
    Extradist1d.push_back(&Muon2TrkdR);
    Extradist1d.push_back(&PhiMass);
    Extradist1d.push_back(&TripleMass);
    Extradist2d.push_back(&PhiMassVsDsMass);

    Extradist1d.push_back(&Muon1_isGlobal);
    Extradist1d.push_back(&Muon2_isGlobal);
    Extradist1d.push_back(&Muon1_isStandAlone);
    Extradist1d.push_back(&Muon2_isStandAlone);
    Extradist1d.push_back(&Muon1_isTracker);
    Extradist1d.push_back(&Muon2_isTracker);
    Extradist1d.push_back(&Muon1_isCalo);
    Extradist1d.push_back(&Muon2_isCalo);
    Extradist1d.push_back(&Muon1_isIsolationValid);
    Extradist1d.push_back(&Muon2_isIsolationValid);
    Extradist1d.push_back(&Track_TriggerMatchdR);
    Extradist1d.push_back(&Muon1_TriggerMatchdR);
    Extradist1d.push_back(&Muon2_TriggerMatchdR);

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // Here you must push back all analysis histograms, otherwise they wont be propagated to the output //
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    Extradist1d.push_back(&NVtx);
}

void  MCEfficiency::doEvent(){ 
    unsigned int t;
    int id(Ntp->GetMCID());
    double phidM = 0.02;
    bool DEBUG = false;
    if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}

    // Apply Selection
    value.at(L1SeedOk) = 0;
    value.at(HLTOk) = 0;
    value.at(trkPt) = 0;
    value.at(nTrkHits) = 0;
    value.at(fitVtxChiSq) = 999.0;
	 value.at(PrimeVtx)=Ntp->NVtx(); 
    value.at(is2MuTrk) = 0;
	 value.at(trigObjMatch) = 0;
	 value.at(genMatch) = 0;

    for(int iTrigger=0; iTrigger < Ntp->NHLT(); iTrigger++){
      TString HLT = Ntp->HLTName(iTrigger);
      if((HLT.Contains("DoubleMu3_Trk_Tau3mu") || HLT.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu")) && Ntp->HLTDecision(iTrigger) == 1)value.at(HLTOk)=Ntp->HLTDecision(iTrigger);
    }
    for(int l1iTrigger=0; l1iTrigger < Ntp->NL1Seeds(); l1iTrigger++){
      TString L1 = Ntp->L1Name(l1iTrigger);
      if(L1.Contains("L1_DoubleMu0er1p4_dEta_Max1p8_OS") && Ntp->L1Decision(l1iTrigger) == 1)value.at(L1SeedOk)=Ntp->L1Decision(l1iTrigger);
      if(L1.Contains("L1_DoubleMu_10_0_dEta_Max1p8") && Ntp->L1Decision(l1iTrigger) == 1)value.at(L1SeedOk)=Ntp->L1Decision(l1iTrigger);
      if(L1.Contains("L1_TripleMu0") && Ntp->L1Decision(l1iTrigger) == 1)value.at(L1SeedOk)=Ntp->L1Decision(l1iTrigger);
    }

    if (DEBUG) cout<<"HLT ok"<<endl;


    vector<int> cut_count;
    vector<int> tmp_nTrkHits;
    vector<double> tmp_trkPt;
    vector<double> tmp_fitVtxChiSq;
    vector<double> tmp_mumuMass;

    int tmp_idx = -1;
    double tmp_chiSq = 999.0;

    bool dsphipi_flag = false;
    bool trigObj_match = false;


    if(Ntp->NTwoMuonsTrack()!=0 && Ntp->NThreeMuons() == 0) value.at(is2MuTrk) = 1;
    if (DEBUG) cout<<"select dimutrk"<<endl;

    if (value.at(is2MuTrk)){
      for(unsigned int i2M=0; i2M < Ntp->NTwoMuonsTrack(); i2M++){

        cut_count.push_back(0);
        tmp_mumuMass.push_back(0);
        tmp_trkPt.push_back(0);
        tmp_nTrkHits.push_back(0);
        tmp_fitVtxChiSq.push_back(999.0);

        unsigned int mu1_idx =  Ntp-> TwoMuonsTrackMuonIndices(i2M).at(0);
        unsigned int mu2_idx =  Ntp-> TwoMuonsTrackMuonIndices(i2M).at(1);
        unsigned int trk_idx = Ntp->TwoMuonsTrackTrackIndex(i2M).at(0);

        if (DEBUG) cout<<mu1_idx<<" "<<mu2_idx<<" "<<trk_idx<<endl;

        if (fabs((Ntp->Muon_P4(mu1_idx) + Ntp->Muon_P4(mu2_idx)).M()-PDG_Var::Phi_mass())<phidM){
            cut_count[i2M]++;
            tmp_mumuMass[i2M] = (Ntp->Muon_P4(mu1_idx)+Ntp->Muon_P4(mu2_idx)).M();
            if (DEBUG) cout<<"Dimuon Invariant mass: "<<tmp_mumuMass[i2M]<<endl;
        }
        else continue;

        if (Ntp->Track_P4(trk_idx).Pt()>cut.at(trkPt)) {
            cut_count[i2M]++;
            tmp_trkPt[i2M] = Ntp->Track_P4(trk_idx).Pt();
            if (DEBUG) cout<<"Track Pt: "<<tmp_trkPt[i2M]<<endl;
        }
        else continue;

        if (Ntp->Track_numberOfValidHits(trk_idx)>=cut.at(nTrkHits)){
            cut_count[i2M]++;
            tmp_nTrkHits[i2M] = Ntp->Track_numberOfValidHits(trk_idx);
            if (DEBUG) cout<<"Track nHits: "<<tmp_nTrkHits[i2M]<<endl;
        }
        else continue;

        if (Ntp->TwoMuonsTrack_SV_Chi2(i2M)<=cut.at(fitVtxChiSq)) {
            cut_count[i2M]++;
            tmp_fitVtxChiSq[i2M] = Ntp->TwoMuonsTrack_SV_Chi2(i2M);
            if (DEBUG) cout<<"Chi square fit vertex: "<<tmp_fitVtxChiSq[i2M]<<endl;
        }
        else continue;

        if (tmp_chiSq>tmp_fitVtxChiSq[i2M]) {
            tmp_idx = i2M;
            tmp_chiSq = tmp_fitVtxChiSq[i2M];
        }


      }

      if (DEBUG){ 
        cout<<"select the best candidate"<<endl;
        cout<<"Index of the best candidate"<<tmp_idx<<endl;
        cout<<tmp_trkPt.size()<<" "<<tmp_nTrkHits.size()<<" "<<tmp_fitVtxChiSq.size()<<endl;
      }
      if (tmp_idx==-1) tmp_idx = std::distance(cut_count.begin(), std::max_element(cut_count.begin(), cut_count.end()));
      value.at(trkPt) = tmp_trkPt[tmp_idx];
      value.at(nTrkHits) = tmp_nTrkHits[tmp_idx];
      value.at(fitVtxChiSq) = tmp_fitVtxChiSq[tmp_idx];
      value.at(mumuMass) = tmp_mumuMass[tmp_idx];
      if (DEBUG) cout<<"Final index"<<tmp_idx<<endl;
		

      
		//Introduce GEN matching
      unsigned int ds_count = 0, dspi_count=0, dsphi_count = 0;

      float ds_genMatch[2];
      vector<unsigned int> dsphi_idxList;
		
		ds_genMatch[0] = 999.;
		ds_genMatch[1] = 999.;

      for (unsigned int ngen=0; ngen<Ntp->NMCSignalParticles(); ngen++){
        bool dspi_flag = 0, dsphi_flag = 0;

        if (fabs(Ntp->MCSignalParticle_pdgid(ngen))==431){
            if (DEBUG) cout<<"Found a Ds"<<endl;
            ds_count++;
            dspi_flag = 0;
            dsphi_flag = 0;

            for (int nchild=0; nchild<Ntp->MCSignalParticle_Nchilds(ngen); nchild++){
              if (DEBUG) cout<<"Looping over all the child particles"<<endl;


              //Find a pi that is a decay product of ds
              if(fabs(Ntp->MCSignalParticle_childpdgid(ngen,nchild))==211){
                dspi_flag = 1;
                if (DEBUG) cout<<"Found a Pion as a child particle"<<endl;
                dspi_flag = 1;
                dspi_count++;
                float tmp_dR = 999.;
                int tmp_pi_idx = Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0);
                float dR = Ntp->Track_P4(tmp_pi_idx).DeltaR(Ntp->MCSignalParticle_child_p4(ngen,nchild));
                if (dR<0.01 && dR<tmp_dR) { 
                   tmp_dR = dR;
                   ds_genMatch[1] = dR;
                   }
              }

              //Find a phi that is a decay product of Ds
              if (fabs(Ntp->MCSignalParticle_childpdgid(ngen,nchild))==333){
                dsphi_flag = 1;
                if (DEBUG) cout<<"Found a Phi as a child particle"<<endl;
                dsphi_idxList.push_back(ngen);
                dsphi_count++;
                dsphi_flag = 1;
                float tmp_dR = 999.;
                unsigned int tmp_mu1_idx = Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0);
                unsigned int tmp_mu2_idx = Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(1);
                float dR = Ntp->MCSignalParticle_child_p4(ngen,nchild).DeltaR(Ntp->Muon_P4(tmp_mu1_idx)+Ntp->Muon_P4(tmp_mu2_idx));
                if (dR<0.01 && dR<tmp_dR) { 
					 	tmp_dR = dR;
                  ds_genMatch[0] = dR;
                  }
              }

            }
            // Check the decay products if Ds doesn't decay to PhiPi
            if (DEBUG && (!dspi_flag || !dsphi_flag)){
              cout<<"Ds->PhiPi decay not found"<<endl;
              for (int nchild = 0; nchild<Ntp->MCSignalParticle_Nchilds(ngen); nchild++){
                cout<<"Daughter PDG Id: "<<Ntp->MCSignalParticle_childpdgid(ngen,nchild)<<endl;
              }
            } 
        } 
      }

      cout<<"GEN Match Indices"<<endl;
		cout<<tmp_idx<<endl;
      printVec(dsphi_idxList.size(), dsphi_idxList);
      cout<<"------------------"<<endl;
      
		cout<<"Number of Ds produced: "<<ds_count<<endl;
      cout<<"Number of Pis form Ds: "<<dspi_count<<endl;
      cout<<"Number of Phis from Ds: "<<dsphi_count<<endl;

      if (ds_genMatch[0]<0.01 && ds_genMatch[1]<0.01) dsphipi_flag = true;
      cout<<"Ds GEN matched: "<<(dsphipi_flag ? "YES":"NO")<<endl;

      // END of GEN matching
      // Begining of trigger object matching


      float mu1_trObj_dR = Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx).at(0);
      float mu2_trObj_dR = Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx).at(1);
      float trk_trObj_dR = Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx).at(2);

      if (DEBUG){
        if (mu1_trObj_dR<0.01) cout<<"Muon 1 matched to a trigger object"<<endl;
        else cout<<"Muon 1 wasn't matched to a trigger object"<<endl;
        if (mu2_trObj_dR<0.01) cout<<"Muon 2 matched to a trigger object"<<endl;
        else cout<<"Muon 2 wasn't matched to a trigger object"<<endl;
        if (trk_trObj_dR<0.01) cout<<"Track matched to a trigger object"<<endl;
        else cout<<"Track wasn't matched to a trigger object"<<endl;
      }

      if(mu1_trObj_dR<0.01 && mu2_trObj_dR<0.01 && trk_trObj_dR<0.01) trigObj_match = true;

      // End of trigger object matching

      value.at(trigObjMatch)=trigObj_match;
      value.at(genMatch)=dsphipi_flag;
    }

    if (DEBUG) cout<<"Values of the variables:"<<value.at(trkPt)<<" "<<value.at(nTrkHits)<<" "<<value.at(fitVtxChiSq)<<endl;


    pass.at(is2MuTrk) = (value.at(is2MuTrk) == cut.at(is2MuTrk));
    pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx)); 
    pass.at(L1SeedOk)= (value.at(L1SeedOk)==cut.at(L1SeedOk)); 
    pass.at(HLTOk)= (value.at(HLTOk)==cut.at(HLTOk)); 
    pass.at(nTrkHits)=(value.at(nTrkHits)>=cut.at(nTrkHits));
    pass.at(trkPt)=(value.at(trkPt)>cut.at(trkPt));
    pass.at(mumuMass)=(value.at(mumuMass)<cut.at(mumuMass)+phidM && value.at(mumuMass)>cut.at(mumuMass)-phidM);
    pass.at(fitVtxChiSq)=(value.at(fitVtxChiSq)<=cut.at(fitVtxChiSq));
	 pass.at(trigObjMatch)=(value.at(trigObjMatch)==cut.at(trigObjMatch));
	 pass.at(genMatch) = (value.at(genMatch)==cut.at(genMatch));

    double wobs=1;
    double w;

    if(!Ntp->isData()){w = 1;} //  No weights to data
    else{w=1;}

    bool status=AnalysisCuts(t,w,wobs);
    if(status){
	 	
		unsigned int mu1 = Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(0);
      unsigned int mu2 = Ntp-> TwoMuonsTrackMuonIndices(tmp_idx).at(1);
      unsigned int track = Ntp-> TwoMuonsTrackTrackIndex(tmp_idx).at(0);

      NVtx.at(t).Fill(Ntp->NVtx(),w);

      Track_Pt.at(t).Fill(Ntp->Track_P4(track).Pt(),w);
      Track_Eta.at(t).Fill(Ntp->Track_P4(track).Eta(),w);
      Track_Phi.at(t).Fill(Ntp->Track_P4(track).Phi(),w);
      Track_E.at(t).Fill(Ntp->Track_P4(track).E(),w);
      Track_P.at(t).Fill(Ntp->Track_P4(track).P(),w);
      Track_vx.at(t).Fill(Ntp->Track_Poca(track).X(),w);
      Track_vy.at(t).Fill(Ntp->Track_Poca(track).Y(),w);
      Track_vz.at(t).Fill(Ntp->Track_Poca(track).Z(),w);
      Track_normalizedChi2.at(t).Fill(Ntp->Track_normalizedChi2(track),w);
      Track_numberOfValidHits.at(t).Fill(Ntp->Track_numberOfValidHits(track),w);
      Track_charge.at(t).Fill(Ntp->Track_charge(track),w);
      Track_dxy.at(t).Fill(Ntp->Track_dxy(track),w);
      Track_dz.at(t).Fill(Ntp->Track_dz(track),w);
      Track_dxyError.at(t).Fill(Ntp->Track_dxyError(track),w);
      Track_dzError.at(t).Fill(Ntp->Track_dzError(track),w);

      Muon1_E.at(t).Fill(Ntp->Muon_P4(mu1).E(),w);
      Muon1_P.at(t).Fill(Ntp->Muon_P4(mu1).P(),w);
      Muon1_Pt.at(t).Fill(Ntp->Muon_P4(mu1).Pt(),w);
      Muon1_Eta.at(t).Fill(Ntp->Muon_P4(mu1).Eta(),w);
      Muon1_Phi.at(t).Fill(Ntp->Muon_P4(mu1).Phi(),w);
      Muon1_vx.at(t).Fill(Ntp->Muon_Poca(mu1).X(),w);
      Muon1_vy.at(t).Fill(Ntp->Muon_Poca(mu1).Y(),w);
      Muon1_vz.at(t).Fill(Ntp->Muon_Poca(mu1).Z(),w);

      Muon2_E.at(t).Fill(Ntp->Muon_P4(mu2).E(),w);
      Muon2_P.at(t).Fill(Ntp->Muon_P4(mu2).P(),w);
      Muon2_Pt.at(t).Fill(Ntp->Muon_P4(mu2).Pt(),w);
      Muon2_Eta.at(t).Fill(Ntp->Muon_P4(mu2).Eta(),w);
      Muon2_Phi.at(t).Fill(Ntp->Muon_P4(mu2).Phi(),w);
      Muon2_vx.at(t).Fill(Ntp->Muon_Poca(mu2).X(),w);
      Muon2_vy.at(t).Fill(Ntp->Muon_Poca(mu2).Y(),w);
      Muon2_vz.at(t).Fill(Ntp->Muon_Poca(mu2).Z(),w);

      Muon1_isGlobal.at(t).Fill(Ntp->Muon_isGlobalMuon(mu1),w);
      Muon2_isGlobal.at(t).Fill(Ntp->Muon_isGlobalMuon(mu2),w);
      Muon1_isStandAlone.at(t).Fill(Ntp->Muon_isStandAloneMuon(mu1),w);
      Muon2_isStandAlone.at(t).Fill(Ntp->Muon_isStandAloneMuon(mu2),w);
      Muon1_isTracker.at(t).Fill(Ntp->Muon_isTrackerMuon(mu1),w);
      Muon2_isTracker.at(t).Fill(Ntp->Muon_isTrackerMuon(mu2),w);
      Muon1_isCalo.at(t).Fill(Ntp->Muon_isCaloMuon(mu1),w);
      Muon2_isCalo.at(t).Fill(Ntp->Muon_isCaloMuon(mu2),w);
      Muon1_isIsolationValid.at(t).Fill(Ntp->Muon_isIsolationValid(mu1),w);
      Muon2_isIsolationValid.at(t).Fill(Ntp->Muon_isIsolationValid(mu2),w);

      // cout<<(Ntp->TwoMuonsTrack_TriggerMatch_dR).size()<<endl;
      //Muon1_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(0),w);
      //Muon2_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(1),w);
      //Track_TriggerMatchdR.at(t).Fill((Ntp->TwoMuonsTrack_TriggerMatch_dR(tmp_idx)).at(2),w);

      DimuondR.at(t).Fill(deltaR(Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0)).Eta(),Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0)).Phi(),Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(1)).Eta(),Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(1)).Phi()));
      Muon1TrkdR.at(t).Fill(deltaR(Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0)).Eta(),Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(0)).Phi(),Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0)).Eta(),Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0)).Phi()));
      Muon2TrkdR.at(t).Fill(deltaR(Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(1)).Eta(),Ntp->Muon_P4(Ntp->TwoMuonsTrackMuonIndices(tmp_idx).at(1)).Phi(),Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0)).Eta(),Ntp->Track_P4(Ntp->TwoMuonsTrackTrackIndex(tmp_idx).at(0)).Phi()));

      double phimass = (Ntp->Muon_P4(mu1)  + Ntp->Muon_P4(mu2)).M();
      double dsmass = (Ntp->Muon_P4(mu1)  + Ntp->Muon_P4(mu2)+ Ntp->Track_P4(track)).M();

		if(fabs(dsmass-PDG_Var::Ds_mass())>phidM) {
			cout<<"Ds mass: "<<dsmass<<endl;
		   cout<<"δ(π mass) = "<<fabs(Ntp->Track_P4(track).M()-PDG_Var::Pi_mass())<<endl;
		   cout<<"δ(φ mass) = "<<fabs((Ntp->Muon_P4(mu1)+Ntp->Muon_P4(mu2)).M()-PDG_Var::Phi_mass())<<endl;
			}

		PhiMass.at(t).Fill(phimass, w);
      TripleMass.at(t).Fill(dsmass, w);
      PhiMassVsDsMass.at(t).Fill(phimass, dsmass);
    }
}

void  MCEfficiency::Finish(){

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


