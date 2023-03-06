#include "ZTau3MuTaumu_Efficiency.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include "PDG_Var.h"
#include <iostream>
//#include "Logger.h"

ZTau3MuTaumu_Efficiency::ZTau3MuTaumu_Efficiency(TString Name_, TString id_):
  Selection(Name_,id_)
{
  // This is a class constructor;
}

ZTau3MuTaumu_Efficiency::~ZTau3MuTaumu_Efficiency(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  ZTau3MuTaumu_Efficiency::Configure(){

  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    
    if(i==nMuon)                    cut.at(nMuon)=3;
    if(i==Conditions_Charge)        cut.at(Conditions_Charge)=1;
    if(i==Conditions_Global_ID)     cut.at(Conditions_Global_ID)=1;
    if(i==Conditions_dR)            cut.at(Conditions_dR)=1;
    if(i==Conditions_InvMass)       cut.at(Conditions_InvMass)=1;
    if(i==Conditions_L1)            cut.at(Conditions_L1)=1;
    if(i==Conditions_HLT)           cut.at(Conditions_HLT)=1;
    if(i==Conditions_muPt)          cut.at(Conditions_muPt)=1;
    if(i==Conditions_3muPt)         cut.at(Conditions_3muPt)=1;
    if(i==Conditions_Random)        cut.at(Conditions_Random)=1;
    if(i==L1_TriggerOk)             cut.at(L1_TriggerOk)=1;
    if(i==HLT_TriggerOk)            cut.at(HLT_TriggerOk)=1;
    if(i==SignalCandidate)          cut.at(SignalCandidate)=1;
    
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
    else if(i==nMuon){
      title.at(i)="Atleast 3 reco muons ";
      hlabel="Atleast 3 reco muons ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nMuon_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nMuon_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==Conditions_Charge){
      title.at(i)="Charge conditions";
      htitle=title.at(i);
      hlabel="Charge conditions";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Conditions_Charge_",htitle,5,-0.5,4.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Conditions_Charge_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }
    else if(i==Conditions_Global_ID){
      title.at(i)="Global ID conditions";
      htitle=title.at(i);
      hlabel="Global_ID conditions";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Conditions_Global_ID_",htitle,5,-0.5,4.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Conditions_Global_ID_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }
    else if(i==Conditions_dR){
      title.at(i)="dR conditions";
      htitle=title.at(i);
      hlabel="dR conditions";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Conditions_dR_",htitle,5,-0.5,4.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Conditions_dR_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }
    else if(i==Conditions_InvMass){
      title.at(i)="InvMass conditions";
      htitle=title.at(i);
      hlabel="InvMass conditions";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Conditions_InvMass_",htitle,5,-0.5,4.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Conditions_InvMass_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }
    else if(i==Conditions_L1){
      title.at(i)="L1 conditions";
      htitle=title.at(i);
      hlabel="L1 conditions";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Conditions_L1_",htitle,5,-0.5,4.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Conditions_L1_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }
    else if(i==Conditions_HLT){
      title.at(i)="HLT conditions";
      htitle=title.at(i);
      hlabel="HLT conditions";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Conditions_HLT_",htitle,5,-0.5,4.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Conditions_HLT_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }
    else if(i==Conditions_muPt){
      title.at(i)="Single mu pt conditions";
      htitle=title.at(i);
      hlabel="Single mu pt conditions";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Conditions_muPt_",htitle,5,-0.5,4.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Conditions_muPt_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }
    else if(i==Conditions_3muPt){
      title.at(i)="3mu pt conditions";
      htitle=title.at(i);
      hlabel="3mu pt conditions";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Conditions_3muPt_",htitle,5,-0.5,4.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Conditions_3muPt_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }
    else if(i==Conditions_Random){
      title.at(i)="Random conditions";
      htitle=title.at(i);
      hlabel="Random conditions";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Conditions_Random_",htitle,5,-0.5,4.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Conditions_Random_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }
    else if(i==SignalCandidate){
      title.at(i)="At least one $\\tau_{3\\mu}$ candidate (3,3,2 GeV,  $|\\eta| < 2.4$, dz($\\mu_{i} , \\mu_{j}$)$<$0.5, dR($\\mu_{i} , \\mu_{j}$)$<$0.8, $\\Sigma \\mu_{charge}$ = +-1)";
      htitle=title.at(i);
      hlabel="N $3\\mu$ candidates";
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SignalCandidate_",htitle,5,-0.5,4.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SignalCandidate_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }

  }
  // Setup NPassed Histogams




  Tau3MuRelativeIsolation=HConfig.GetTH1D(Name+"_Tau3MuRelativeIsolation","Tau3MuRelativeIsolation",50,0.,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T})","#Delta R < 0.4 ");
  OppositeMuRelativeIsolation=HConfig.GetTH1D(Name+"_OppositeMuRelativeIsolation","OppositeMuRelativeIsolation",50,0.,1.1,"I= p_{T}(#tau)/(p_{T}(#tau) + #sum p_{T})","#Delta R < 0.4 ");
  
  VisibleDiTauMass=HConfig.GetTH1D(Name+"_VisibleDiTauMass","VisibleDiTauMass",70,0.,150,"M_{#tau(#mu) - #tau(3#mu)}, GeV (visible mass)","Events");
  MTT=HConfig.GetTH1D(Name+"_MTT","MTT",70,0.,140,"M_{#tau(#mu) - #tau(3#mu)}, GeV (collinear approximation)","Events");
  TripletMass=HConfig.GetTH1D(Name+"_TripletMass","TripletMass",30,1.1,2.2,"M_{3#mu}, GeV","Events");


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
  
  Selection_Cut_Sequential_1=HConfig.GetTH1D(Name+"_Selection_Cut_Sequential_1","Selection_Cut_Sequential_1",8,-0.5,7.5,"#mu no","Events");
  Selection_Cut_Sequential_2=HConfig.GetTH1D(Name+"_Selection_Cut_Sequential_2","Selection_Cut_Sequential_2",60,0.0,1.5,"dR","Events");
  Selection_Cut_Sequential_3=HConfig.GetTH1D(Name+"_Selection_Cut_Sequential_3","Selection_Cut_Sequential_3",60,1.0,2.5,"Inv Mass","Events");
  Selection_Cut_Sequential_4=HConfig.GetTH1D(Name+"_Selection_Cut_Sequential_4","Selection_Cut_Sequential_4",80,0.0,80,"sum 3mu pT","Events");
  Selection_Cut_Sequential_5=HConfig.GetTH1D(Name+"_Selection_Cut_Sequential_5","Selection_Cut_Sequential_5",100,0.0,5.0,"Tau3MuIsolation","Events");
  Selection_Cut_Sequential_6=HConfig.GetTH1D(Name+"_Selection_Cut_Sequential_6","Selection_Cut_Sequential_6",40,0.0,20,"largest Mu pT","Events");

  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events"); // Do not remove
  // Setup Extra Histograms



  Selection::ConfigureHistograms(); //do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour); // do not remove
}




void  ZTau3MuTaumu_Efficiency::Store_ExtraDist(){ 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Here you must push back all analysis histograms, otherwise they wont be propagated to the output


  Extradist1d.push_back(&Tau3MuRelativeIsolation);
  Extradist1d.push_back(&OppositeMuRelativeIsolation);
  Extradist1d.push_back(&VisibleDiTauMass);
  Extradist1d.push_back(&MTT);
  Extradist1d.push_back(&TripletMass);
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
  
  Extradist1d.push_back(&Selection_Cut_Sequential_1);
  Extradist1d.push_back(&Selection_Cut_Sequential_2);
  Extradist1d.push_back(&Selection_Cut_Sequential_3);
  Extradist1d.push_back(&Selection_Cut_Sequential_4);
  Extradist1d.push_back(&Selection_Cut_Sequential_5);
  Extradist1d.push_back(&Selection_Cut_Sequential_6);


}

////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is called on each event

void  ZTau3MuTaumu_Efficiency::doEvent(){ 
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
  
  random_num = rndm.Rndm();
  
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

  value.at(L1_TriggerOk)=(L1Ok);
  pass.at(L1_TriggerOk)=(value.at(L1_TriggerOk)==cut.at(L1_TriggerOk));
  
  value.at(HLT_TriggerOk)=(HLTOk);
  pass.at(HLT_TriggerOk)=(value.at(HLT_TriggerOk)==cut.at(HLT_TriggerOk));
  
  //std::cout << "Test 1." << std::endl;
  
  unsigned int  signal_idx=0;
  double min_chi2(99.);

  for(unsigned int i_idx =0; i_idx < Ntp->NThreeMuons(); i_idx++){
    if(Ntp->Vertex_Signal_KF_Chi2(i_idx) < min_chi2){
      min_chi2 = Ntp->Vertex_Signal_KF_Chi2(i_idx);
      signal_idx = i_idx;
    }
  }
  
  int RecoMu_Count = Ntp->NMuons();
  value.at(nMuon) = RecoMu_Count;
  pass.at(nMuon) = (value.at(nMuon) >= cut.at(nMuon));
  
  value.at(SignalCandidate) = Ntp->NThreeMuons();
  pass.at(SignalCandidate) = (value.at(SignalCandidate) >= cut.at(SignalCandidate));
  
  std::vector<float> pts; // Transverse momenta
  std::vector<float> ps; // Momenta
  std::vector<float> phis;// signed values of phi
  std::vector<float> etas;// signed values of eta
  std::vector<int> charges;// charges
  
  std::vector<float> es;// muon energy
  std::vector<float> p1;// muon px
  std::vector<float> p2;
  std::vector<float> p3;
  
  for(unsigned int imu=0; imu < Ntp->NMuons(); imu++){
        TLorentzVector Muon_Reco_LV=Ntp->Muon_P4(imu);
        
        pts.push_back(Muon_Reco_LV.Pt());
        ps.push_back(Muon_Reco_LV.Vect().Mag());
        etas.push_back(Muon_Reco_LV.Eta());
        phis.push_back(Muon_Reco_LV.Phi());
        charges.push_back(Ntp->Muon_charge(imu));
        
        es.push_back(Muon_Reco_LV.E());
        p1.push_back(Muon_Reco_LV.Px());
        p2.push_back(Muon_Reco_LV.Py());
        p3.push_back(Muon_Reco_LV.Pz());
        
        
  }
  
  using Row = vector<double>;
  using Matrix = vector<Row>;
  
  Matrix HLT_pass;//stores whether a particular combination of 3 muons would pass the HLT.
                             //Condition for track not applied. Inv Mass considered separately.
  Matrix HLT_restrictive_pass;// restrictive HLT. Checks for common vertex and invariant mass
  Matrix L1_pass;//if the 3 muons would pass L1 Double mu condition. Doesn't take into account events with only 2 muons, so triple mu criteria always satisfied
  Matrix Charge_pass;//if total charge of muons is +/- 1
  Matrix InvMass_pass;//Pass if invariant mass is between 1.55-2.1 GeV
  Matrix dR_pass;// Check if all three pairs of muons are within maxDr_ of each other
  Matrix Global_ID_pass;// Check for pT and eta conditions of global ID
  
  Matrix Single_mu_pt_pass;// other checks
  Matrix Three_mu_pt_pass;// other checks
  Matrix Random_pass;// other checks
  
  int reserve_size=((RecoMu_Count*(RecoMu_Count-1)*(RecoMu_Count-2))/6);
  if(RecoMu_Count>=3){// reserve space for pass vectors so that the code runs faster
      
      HLT_pass.reserve(reserve_size);
      HLT_restrictive_pass.reserve(reserve_size);
      L1_pass.reserve(reserve_size);
      Charge_pass.reserve(reserve_size);
      InvMass_pass.reserve(reserve_size);
      dR_pass.reserve(reserve_size);
      Global_ID_pass.reserve(reserve_size);
      
      Single_mu_pt_pass.reserve(reserve_size);
      Three_mu_pt_pass.reserve(reserve_size);
      Random_pass.reserve(reserve_size);
  }
  
  IsL1 isl1;
  IsHLT ishlt;
  IsdR isdr;
  
  double ptMin_(3.0);
  double pMin_(2.0);
  double etaMax_(2.4);
  double invMassMin_(1.4);
  double invMassMax_(2.1);
  double maxDr_(0.3);
  
  int x(0);// stores number of unique muon combinations
  
  if(RecoMu_Count>=3){// need atleast 3 muons.
    for (unsigned int A = 0; A < RecoMu_Count; A++) {//all possible combinations of muons (indexed by A, B and C)
        for (unsigned int B = 0; B < RecoMu_Count; B++) {
            for (unsigned int C = 0; C < RecoMu_Count; C++) {
                if(A<B&&B<C){//All unique combinations of muons
                    x+=1;
                    std::vector<unsigned int> mu_i = {A,B,C};
                    
                    // Look at their invariant mass
                    TLorentzVector sum_3mu(p1.at(A)+p1.at(B)+p1.at(C),p2.at(A)+p2.at(B)+p2.at(C),p3.at(A)+p3.at(B)+p3.at(C),es.at(A)+es.at(B)+es.at(C));
                    bool three_mu_inv = (sum_3mu.M()>invMassMin_)&&(sum_3mu.M()<invMassMax_);
                    InvMass_pass.push_back({three_mu_inv,sum_3mu.M()});
                    
                    
                    // Check the HLT pTs
                    bool HLT(false);
                    for (int i = 0; i < 3; i++){//sends 3 pairs of muons to the IsHLT class
                        // (b + (a%b)) % b is used to get a positive value for negative modulo
                        if(ishlt(pts.at(mu_i.at((3 + ((i-1)%3)) % 3)),pts.at(mu_i.at((3 + ((i+1)%3)) % 3)),ptMin_,ptMin_)){
                            HLT=true;
                            break;
                        }
                    }
                    HLT_pass.push_back({HLT});
                    
                    // Check the L1 double mu pT and etas
                    bool L1T(false);
                    for (int i = 0; i < 3; i++){//sends 3 pairs of muons to the IsL1 class
                        // (b + (a%b)) % b is used to get a positive value for negative modulo
                        if(isl1(pts.at(mu_i.at((3 + ((i-1)%3)) % 3)),pts.at(mu_i.at((3 + ((i+1)%3)) % 3)),etas.at(mu_i.at((3 + ((i-1)%3)) % 3)),etas.at(mu_i.at((3 + ((i+1)%3)) % 3)))){
                            L1T=true;
                            break;
                        }
                    }
                    L1_pass.push_back({L1T});
                    
                    // Check dR
                    bool dRp(true);
                    double maxdRinTriplet(-1);
                    for (int i = 0; i < 3; i++){//sends 3 pairs of muons to the IsdR class
                        // (b + (a%b)) % b is used to get a positive value for negative modulo
                        double dR_pair = isdr(phis.at(mu_i.at((3 + ((i-1)%3)) % 3)),phis.at(mu_i.at((3 + ((i+1)%3)) % 3)),etas.at(mu_i.at((3 + ((i-1)%3)) % 3)),etas.at(mu_i.at((3 + ((i+1)%3)) % 3)));
                        if(dR_pair>maxdRinTriplet){
                          maxdRinTriplet=dR_pair;
                        }
                        if(dR_pair>maxDr_){
                            dRp=false;// if one pair has greater dR than maxdR, it doesn't pass
                            break;
                        }
                    }
                    dR_pass.push_back({dRp,maxdRinTriplet});
                    
                    
                    // Check the charges
                    Charge_pass.push_back({abs(charges.at(A)+charges.at(B)+charges.at(C))==1});
                    
                    
                    // Check for Global Muon ID pT and eTa restrictions
                    Global_ID_pass.push_back({(pts.at(A)>ptMin_)&&(pts.at(B)>ptMin_)&&(pts.at(C)>ptMin_)&&(fabs(etas.at(A))<etaMax_)&&(fabs(etas.at(B))<etaMax_)&&(fabs(etas.at(C))<etaMax_)});
                    
                    // Other checks tested. Incl. Luca analysis notes.
                    bool Single_Mu_pT(false);
                    double largest_Mu_pT(-1);
                    for (int i = 0; i < 3; i++){
                        if(pts.at(mu_i.at(i))>7.0){
                            Single_Mu_pT=true;// 
                        }
                        if(pts.at(mu_i.at(i))>largest_Mu_pT){
                          largest_Mu_pT=pts.at(mu_i.at(i));
                        }
                    }
                    bool ThreeMu_pT=(sum_3mu.Pt()>18.0);
                    
                    double three_mu_sum_pt = (pts.at(A)+pts.at(B)+pts.at(C));
                    
                    double mu1_pt_sum(0.0);
                    double mu2_pt_sum(0.0);
                    double mu3_pt_sum(0.0);
                    for (int i = 0; i < Ntp->NTracks(); i++){
                      TLorentzVector Track_LV = Ntp->Track_P4(i);
                      if(Ntp->Muon_P4(A).DeltaR(Track_LV)<0.4&&Ntp->Track_dz(i)<0.2){
                        mu1_pt_sum+=Track_LV.Pt();
                      }
                      if(Ntp->Muon_P4(B).DeltaR(Track_LV)<0.4&&Ntp->Track_dz(i)<0.2){
                        mu2_pt_sum+=Track_LV.Pt();
                      }
                      if(Ntp->Muon_P4(C).DeltaR(Track_LV)<0.4&&Ntp->Track_dz(i)<0.2){
                        mu3_pt_sum+=Track_LV.Pt();
                      }
                    }
                    
                    //double mu1_iso1 = Ntp->Muon_sumChargedParticlePt04(A) + std::max(0., Ntp->Muon_sumPhotonEt04(A) -0.2 * Ntp->Muon_sumPUPt04(A));
                    //double mu2_iso1 = Ntp->Muon_sumChargedParticlePt04(B) + std::max(0., Ntp->Muon_sumPhotonEt04(B) -0.2 * Ntp->Muon_sumPUPt04(B));
                    //double mu3_iso1 = Ntp->Muon_sumChargedParticlePt04(C) + std::max(0., Ntp->Muon_sumPhotonEt04(C) -0.2 * Ntp->Muon_sumPUPt04(C));
                    
                    //double mu1_iso_test = mu1_pt_sum + std::max(0., Ntp->Muon_sumPhotonEt04(A) -0.2 * Ntp->Muon_sumPUPt04(A));
                    //double mu2_iso_test = mu2_pt_sum + std::max(0., Ntp->Muon_sumPhotonEt04(B) -0.2 * Ntp->Muon_sumPUPt04(B));
                    //double mu3_iso_test = mu3_pt_sum + std::max(0., Ntp->Muon_sumPhotonEt04(C) -0.2 * Ntp->Muon_sumPUPt04(C));
                    
                    //double mu1_iso_test1 = Ntp->Muon_sumChargedParticlePt04(A) + Ntp->Muon_sumPhotonEt04(A) -0.2 * Ntp->Muon_sumPUPt04(A);
                    //double mu2_iso_test1 = Ntp->Muon_sumChargedParticlePt04(B) + Ntp->Muon_sumPhotonEt04(B) -0.2 * Ntp->Muon_sumPUPt04(B);
                    //double mu3_iso_test1 = Ntp->Muon_sumChargedParticlePt04(C) + Ntp->Muon_sumPhotonEt04(C) -0.2 * Ntp->Muon_sumPUPt04(C);
                    
                    //double mu1_iso_test2 = mu1_pt_sum + Ntp->Muon_sumPhotonEt04(A) -0.2 * Ntp->Muon_sumPUPt04(A);
                    //double mu2_iso_test2 = mu2_pt_sum + Ntp->Muon_sumPhotonEt04(B) -0.2 * Ntp->Muon_sumPUPt04(B);
                    //double mu3_iso_test2 = mu3_pt_sum + Ntp->Muon_sumPhotonEt04(C) -0.2 * Ntp->Muon_sumPUPt04(C);
                    
                    //double mu1_iso2 = Ntp->Muon_sumChargedParticlePt03(A) + std::max(0., Ntp->Muon_sumPhotonEt03(A) -0.2 * Ntp->Muon_sumPUPt03(A));
                    //double mu2_iso2 = Ntp->Muon_sumChargedParticlePt03(B) + std::max(0., Ntp->Muon_sumPhotonEt03(B) -0.2 * Ntp->Muon_sumPUPt03(B));
                    //double mu3_iso2 = Ntp->Muon_sumChargedParticlePt03(C) + std::max(0., Ntp->Muon_sumPhotonEt03(C) -0.2 * Ntp->Muon_sumPUPt03(C));
                    
                    //double Tau3MuIsolation = (Ntp->Muon_RelIso(A) + Ntp->Muon_RelIso(B) + Ntp->Muon_RelIso(C))/  (3* sum_3mu.Pt() );
                    //double Tau3MuIsolation = sum_3mu.Pt()/  (sum_3mu.Pt()  + mu1_iso + mu2_iso + mu3_iso );
                    //double Tau3MuIsolation = three_mu_sum_pt/  (three_mu_sum_pt  + mu1_iso + mu2_iso + mu3_iso - 2*three_mu_sum_pt);
                    //double Tau3MuIsolation = three_mu_sum_pt/  (three_mu_sum_pt  + mu1_iso_test + mu2_iso_test + mu3_iso_test - 2*three_mu_sum_pt);//seems to distinguish events that pass HLT and those that don't
                    //double Tau3MuIsolation = three_mu_sum_pt/  (three_mu_sum_pt  + mu1_iso_test + mu2_iso_test + mu3_iso_test);// why does this give 2 peaks? because some events include the tracks that can be matched to muons? because of std::max?
                    //double Tau3MuIsolation = three_mu_sum_pt/  (three_mu_sum_pt  + mu1_iso_test + mu2_iso_test + mu3_iso_test - three_mu_sum_pt);//seems to distinguish events that pass HLT and those that don't
                    //double Tau3MuIsolation = three_mu_sum_pt/  (three_mu_sum_pt  + std::max(0., mu1_iso_test + mu2_iso_test + mu3_iso_test - three_mu_sum_pt));
                    //double Tau3MuIsolation = three_mu_sum_pt/  (three_mu_sum_pt  + mu1_iso_test + mu2_iso_test + mu3_iso_test);
                    //double Tau3MuIsolation = (mu1_iso_test + mu2_iso_test + mu3_iso_test - 3*three_mu_sum_pt)/three_mu_sum_pt;
                    
                    //double Tau3MuIsolation = three_mu_sum_pt/  (three_mu_sum_pt  + mu1_iso_test1 + mu2_iso_test1 + mu3_iso_test1 - 2*three_mu_sum_pt);
                    //double Tau3MuIsolation = three_mu_sum_pt/  (three_mu_sum_pt  + mu1_iso_test2 + mu2_iso_test2 + mu3_iso_test2 - 3*three_mu_sum_pt);
                    
                    double mu1_iso3 = Ntp->Muon_sumChargedParticlePt03(A);
                    double mu2_iso3 = Ntp->Muon_sumChargedParticlePt03(B);
                    double mu3_iso3 = Ntp->Muon_sumChargedParticlePt03(C);
                    double Tau3MuIsolation = (mu1_iso3 + mu2_iso3 + mu3_iso3)/three_mu_sum_pt;
                    
                    
                    Single_mu_pt_pass.push_back({Single_Mu_pT});
                    Three_mu_pt_pass.push_back({ThreeMu_pT});
                    Random_pass.push_back({(Tau3MuIsolation>1.99)&&(Tau3MuIsolation<2.5),sum_3mu.Pt(),Tau3MuIsolation,largest_Mu_pT,pts.at(A),pts.at(B),pts.at(C)});
                  
                }// end unique muon if statement
            }// C loop
        }// B loop
    }// A loop
  }// end pts.size()>=3
  
  
  bool Final_pass(false);
  bool Final_pass_Charge(false);
  bool Final_pass_Global_ID(false);
  bool Final_pass_dR(false);
  bool Final_pass_InvMass(false);
  bool Final_pass_HLT(false);
  bool Final_pass_L1(false);
  
  bool Final_pass_Single_mu_pt(false);
  bool Final_pass_Three_mu_pt(false);
  bool Final_pass_Random(false);
  
  double dR_min_of_pairMax(99.0);
  double central_3mu_Mass(99.0);
  double largest_mu_pt_final(-1.0);
  double largest_3mu_pt_final(-1.0);
  double smallest_3mu_isolation(99.0);
  
  for (unsigned int comb_idx = 0; comb_idx < Charge_pass.size(); comb_idx++) {//comb_idx indexes all unique combinations of muons
      
      if(Charge_pass[comb_idx][0]){
        Final_pass_Charge=true;
        
        if(Global_ID_pass[comb_idx][0]){
          Final_pass_Global_ID=true;
          
          if(dR_pass[comb_idx][1]<dR_min_of_pairMax){
            dR_min_of_pairMax=dR_pass[comb_idx][1];
          }
          if(dR_pass[comb_idx][0]){
            Final_pass_dR=true;
            
            if(fabs(InvMass_pass[comb_idx][1]-1.776)   < fabs(central_3mu_Mass-1.776)  ){
              central_3mu_Mass=InvMass_pass[comb_idx][1];
            }
            if(InvMass_pass[comb_idx][0]){
              Final_pass_InvMass=true;
              
              if(L1_pass[comb_idx][0]){
                Final_pass_L1=true;
                
                if(HLT_pass[comb_idx][0]){
                  Final_pass_HLT=true;
                  
                  
                  if(Random_pass[comb_idx][3]>largest_mu_pt_final){
                    largest_mu_pt_final=Random_pass[comb_idx][3];
                  }
                  
                  if(Single_mu_pt_pass[comb_idx][0]){
                    Final_pass_Single_mu_pt=true;
                    
                    if(Random_pass[comb_idx][1]>largest_3mu_pt_final){
                      largest_3mu_pt_final=Random_pass[comb_idx][1];
                    }
                    
                    if(Three_mu_pt_pass[comb_idx][0]){
                      Final_pass_Three_mu_pt=true;
                      
                      if(Random_pass[comb_idx][2]<smallest_3mu_isolation){
                        smallest_3mu_isolation=Random_pass[comb_idx][2];
                      }
                      
                      if(Random_pass[comb_idx][0]){
                        Final_pass_Random=true;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
  }
  
  pass.at(Conditions_Charge) = Final_pass_Charge;
  pass.at(Conditions_Global_ID) = Final_pass_Global_ID;
  pass.at(Conditions_dR) = Final_pass_dR;
  pass.at(Conditions_InvMass) = Final_pass_InvMass;
  pass.at(Conditions_L1) = Final_pass_L1;
  pass.at(Conditions_HLT) = Final_pass_HLT;
  pass.at(Conditions_muPt) = Final_pass_Single_mu_pt;
  pass.at(Conditions_3muPt) = Final_pass_Three_mu_pt;
  pass.at(Conditions_Random) = Final_pass_Random;
  
    double wobs=1;
    double w;  
             
  if(!Ntp->isData()){w = 1; /*Ntp->PUReweight(); */} //  No weights to data
  else{w=1;}
  
  if(passAllUntil(HLT_TriggerOk)) Selection_Cut_Sequential_1.at(t).Fill(RecoMu_Count,1);
  
  if(RecoMu_Count>=3){
          
          if(passAllUntil(Conditions_Global_ID)) Selection_Cut_Sequential_2.at(t).Fill(dR_min_of_pairMax,1);
          if(passAllUntil(Conditions_dR)) Selection_Cut_Sequential_3.at(t).Fill(central_3mu_Mass,1);
          if(passAllUntil(Conditions_HLT)) Selection_Cut_Sequential_4.at(t).Fill(largest_3mu_pt_final,1);
          if(passAllUntil(Conditions_HLT)) Selection_Cut_Sequential_5.at(t).Fill(smallest_3mu_isolation,1);
          if(passAllUntil(Conditions_HLT)) Selection_Cut_Sequential_6.at(t).Fill(largest_mu_pt_final,1);
          
  }
  
  HLT_pass.clear();
  L1_pass.clear();
  Charge_pass.clear();
  InvMass_pass.clear();
  dR_pass.clear();
  Global_ID_pass.clear();
  HLT_restrictive_pass.clear();
  Random_pass.clear();
  
  bool status=AnalysisCuts(t,w,wobs);

  if(status){ 
    
  }
}


void  ZTau3MuTaumu_Efficiency::Finish(){

  Selection::Finish();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This function is called after the event loop and you can code here any analysis with already filled analysis histograms 
}





