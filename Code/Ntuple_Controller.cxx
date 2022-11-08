//Ntuple_Controller.cxx IMPLEMENTATION FILE
#include "Ntuple_Controller.h"
#include "Parameters.h"
#include <tuple>
#include "Logger.h"
#include "PDGInfo.h"
using namespace std;


// External code

///////////////////////////////////////////////////////////////////////
//
// Constructor
//
///////////////////////////////////////////////////////////////////////
Ntuple_Controller::Ntuple_Controller(std::vector<TString> RootFiles):
   copyTree(false)
   ,ObjEvent(-1)
   ,isInit(false)
{
   // TChains the ROOTuple file
   TChain *chain = new TChain("T3MTree/t3mtree");
   //  Logger(Logger::Verbose) << "Loading " << RootFiles.size() << " files" << std::endl;
   for(unsigned int i=0; i<RootFiles.size(); i++){
      chain->Add(RootFiles[i]);
   }
   TTree *tree = (TTree*)chain;
   if(chain==0){
      //	Logger(Logger::Error) << "chain points to NULL" << std::endl;
   }
   std::cout << "Number of Events in Ntuple: " << chain->GetEntries() << std::endl;
   Ntp=new NtupleReader(tree);
   nbytes=0; 
   nb=0;
   std::cout << "Ntuple Configured" << std::endl;

   gRandom->SetSeed(1234);
}

///////////////////////////////////////////////////////////////////////
//
// Function: void InitEvent()
//
// Purpose: Initialize variables etc on event base
//
///////////////////////////////////////////////////////////////////////

void Ntuple_Controller::InitEvent(){
   Muon_corrected_p4.clear();
   //	Muon_corrected_p4.resize(NMuons());
   Muon_isCorrected = false;

   // after everything is initialized
   isInit = true;
}

///////////////////////////////////////////////////////////////////////
//
// Function: Int_t Get_Entries()
//
// Purpose: To get the number of events in the Ntuple
//
///////////////////////////////////////////////////////////////////////
Int_t Ntuple_Controller::Get_Entries(){
   return Int_t(Ntp->fChain->GetEntries());
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Get_Event(int _jentry)
//
// Purpose: To get the event _jentry
//
///////////////////////////////////////////////////////////////////////
void Ntuple_Controller::Get_Event(int _jentry){
   jentry=_jentry;
   Ntp->LoadTree(jentry);
   nb = Ntp->fChain->GetEntry(jentry);   nbytes += nb;
   isInit = false;
   InitEvent();
}


///////////////////////////////////////////////////////////////////////
//
// Function: void Get_EventIndex()
//
// Purpose: To get the event index (jentry)
//
///////////////////////////////////////////////////////////////////////
Int_t Ntuple_Controller::Get_EventIndex(){
   return jentry;
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Get_Event(int _jentry)
//
// Purpose: To get the file name of the root file currently being 
//          accesses
//
///////////////////////////////////////////////////////////////////////
TString Ntuple_Controller::Get_File_Name(){
   return Ntp->fChain->GetCurrentFile()->GetName();
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Branch_Setup(TString B_Name, int type)
//
// Purpose: To setup a branch
//
///////////////////////////////////////////////////////////////////////
void Ntuple_Controller::Branch_Setup(TString B_Name, int type){   
   Ntp->fChain->SetBranchStatus(B_Name,type);
}

///////////////////////////////////////////////////////////////////////
//
// destructor
//
///////////////////////////////////////////////////////////////////////
Ntuple_Controller::~Ntuple_Controller() {
   std::cout << "Cleaning up" << std::endl;
   delete Ntp;
   //  delete rmcor;
   std::cout << "Complete." << std::endl;
}


void Ntuple_Controller::CloneTree(TString n){
   if(!copyTree){
      std::cout << "Starting D3PD cloning" << std::endl;
      newfile = new TFile(n+".root","recreate");
      SkimmedTree=Ntp->fChain->CloneTree(0);
      copyTree=true;
   }
}

void Ntuple_Controller::SaveCloneTree(){
   if(copyTree){
      SkimmedTree->AutoSave();
      newfile->Close();
   }
   std::cout << "Done"<< std::endl;
}

void Ntuple_Controller::ThinTree(){
   std::cout << "ThinTree not implemented..." << std::endl;
}

int Ntuple_Controller::SetupSystematics(TString sys){
   return Default;
}


void Ntuple_Controller::ConfigureObjects(){
   if(ObjEvent!=EventNumber()){
      ObjEvent=EventNumber();
   }
}



//Physics get Functions
Long64_t  Ntuple_Controller::GetMCID(){

   Long64_t  DataMCTypeFromTupel =  Ntp->Event_DataMC_Type;

   if(DataMCTypeFromTupel==40 || DataMCTypeFromTupel==60 ||  DataMCTypeFromTupel==90 ){
      //    return 555;
   }
   return DataMCTypeFromTupel;

}

// return DataMCType without mass information
int Ntuple_Controller::GetStrippedMCID(){
   return GetMCID() % 100;
}

// return path of input dataset (as given in the line InputNtuples in Input.txt)
TString Ntuple_Controller::GetInputNtuplePath(){
   Parameters Par; // assumes configured in Analysis.cxx
   TString dsPath;
   Par.GetString("InputNtuples:",dsPath);
   return dsPath;
}

// return name of input dataset
TString Ntuple_Controller::GetInputDatasetName(){
   TString dsPath = GetInputNtuplePath();
   return gSystem->BaseName( gSystem->DirName( gSystem->DirName(dsPath) ) );
}

//
TString Ntuple_Controller::GetInputPublishDataName(){
   TString dsPath = GetInputNtuplePath();
   return gSystem->BaseName( gSystem->DirName(dsPath) );
}

// determine Higgs mass (from Dataset name or fallback options)

double Ntuple_Controller::DeltaPhi(double angle1,double angle2)
{
   double diff=angle1-angle2;
   while(diff>=TMath::Pi())diff=diff-TMath::TwoPi();
   while(diff<=-TMath::Pi())diff=diff+TMath::TwoPi();
   return diff;
}

//-------------------------- print vector -------------------
template<typename T>
void printVec(int size, T& vec){
   for (int i = 0; i<size; i++) cout<<vec[i]<<" ";
   cout<<endl;
}
//-------------------------- --------------------------------- 

//-------------------------- return deltaR -------------------

double Ntuple_Controller::deltaR(double eta1, double phi1, double eta2, double phi2)
{
   double deta = abs(eta1 - eta2);
   double dphi = fabs(phi1 - phi2);

   if (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();

   return sqrt(deta*deta + dphi*dphi);
}
//-------------------------- --------------------------------- 


//-------------------------- return GEN pdgId -------------------
// Returns the pdgId of the GEN particle matched to a RECO particle

std::pair<int, int> Ntuple_Controller::GENMatchedPdgId(TLorentzVector vec){
   std::pair<int, int> tmp;
   float dR_min = 999.0;
   int pdgId = 9999999;
   tmp.first = pdgId;
   tmp.second = -1;
   for (unsigned int it=0; it<NMCParticles(); ++it){
      float tmp_dR = MCParticle_p4(it).DeltaR(vec);
      float dpt = abs(MCParticle_p4(it).Pt()-vec.Pt())/vec.Pt();
      if (tmp_dR<dR_min && tmp_dR<0.03 && dpt<0.1){
         dR_min = tmp_dR;
         pdgId = MCParticle_pdgid(it);
         tmp.second = it;
         tmp.first = pdgId;
      }
   } 
   return tmp;
}

//-------------------------- --------------------------------- -------------------
// Trigger matching
//-------------------------- --------------------------------- -------------------
std::pair<bool, std::vector<float>> Ntuple_Controller::triggerMatchTriplet(std::vector<TLorentzVector> v1, std::vector<TLorentzVector> v2){

   std::pair<bool, std::vector<float>> matched;
   int SIZE = v2.size();
   matched.first = false;

   for (int i=0; i<3; i++) matched.second.push_back(999.);

   double dpt[3][SIZE];
   double dR[3][SIZE];

   for (int i=0; i<3; i++){
      for (int j=0; j<SIZE; j++){
         double p = fabs(v1[i].Pt()-v2[j].Pt())/v1[i].Pt();
         double r = v1[i].DeltaR(v2[j]);
         dpt[i][j] = p;
         dR[i][j] = r;
      }
   }

   // Find if any combination matches the trigger matching criteria
   for (int i=0; i<SIZE; i++){
      for (int j=0; j<SIZE; j++){
         if (i==j) continue;
         for (int k=0; k<SIZE; k++){
            if (k==i || k==j) continue;
            if ( dpt[0][i]<0.1 && dpt[1][j]<0.1 && dpt[2][k]<0.1 &&
                  dR[0][i]<0.03 && dR[1][j]<0.03 && dR[2][k]<0.03 ) {
               matched.first = true;
               matched.second.at(0) = dR[0][i]; 
               matched.second.at(1) = dR[1][j]; 
               matched.second.at(2) = dR[2][k]; 
            }
         }
      }
   }

   return matched;
}

//-------------------------- return GEN LV -------------------
// Returns the lorentz vector of the GEN particle matched to a RECO particle

TLorentzVector Ntuple_Controller::GENMatchedLV(TLorentzVector vec){
   float dR_min = 999.0;
   TLorentzVector tmp_vec(0,0,0,0);
   for (unsigned int it=0; it<NMCSignalParticles(); ++it){
      float tmp_dR = MCSignalParticle_p4(it).DeltaR(vec);
      if (tmp_dR<dR_min && tmp_dR<0.05){
         dR_min = tmp_dR;
         tmp_vec = MCSignalParticle_p4(it);
      }
   } 
   return tmp_vec;
}

//-------------------------- --------------------------------- -------------------
//-------------------------- --------------------------------- -------------------
float Ntuple_Controller::ThreeMuonDsGenMatch(unsigned int tmp_idx){
   float ds_dR = 999.0;

   if (tmp_idx>=NThreeMuons()) {
      cout<<"Index out of range!"<<endl;
      return ds_dR;
   }

   unsigned int Muon_index_1 = ThreeMuonIndices(tmp_idx).at(0);
   unsigned int Muon_index_2 = ThreeMuonIndices(tmp_idx).at(1);
   unsigned int Muon_index_3 = ThreeMuonIndices(tmp_idx).at(2);

   for (unsigned int j=0; j<NTwoMuonsTrack(); ++j){
      unsigned int mu1_idx = TwoMuonsTrackMuonIndices(j).at(0);
      unsigned int mu2_idx = TwoMuonsTrackMuonIndices(j).at(1);
      unsigned int trk_idx = TwoMuonsTrackTrackIndex(j).at(0);

      float dRMu1Trk = Muon_P4(Muon_index_1).DeltaR(Track_P4(trk_idx));
      float dRMu2Trk = Muon_P4(Muon_index_2).DeltaR(Track_P4(trk_idx));
      float dRMu3Trk = Muon_P4(Muon_index_3).DeltaR(Track_P4(trk_idx));

      if ( (((mu1_idx==Muon_index_1 && mu2_idx==Muon_index_2) || (mu1_idx==Muon_index_2 && mu2_idx==Muon_index_1)) && (dRMu3Trk<0.01)) ||
           (((mu1_idx==Muon_index_1 && mu2_idx==Muon_index_2) || (mu1_idx==Muon_index_2 && mu2_idx==Muon_index_1)) && (dRMu2Trk<0.01)) ||
           (((mu1_idx==Muon_index_1 && mu2_idx==Muon_index_2) || (mu1_idx==Muon_index_2 && mu2_idx==Muon_index_1)) && (dRMu1Trk<0.01))
            ) {
        float tmp_ds_dR = DsGenMatch(j);
        if (tmp_ds_dR<ds_dR){
           ds_dR = tmp_ds_dR;
        }
      }
   }
   return ds_dR;
}

//-------------------------- Match GEN Ds to 2mu+trk candidate -------------------
// Returns the dR of best GEN machted Ds
// The function takes the index of the two muons + track candidate as the input
// and iterates over all the GEN Ds candidates. It matches the track to a daughter pion and
// the dimuon to a daughter Phi separately. If both matches are found, then it saves the dR.
//
float Ntuple_Controller::DsGenMatch(unsigned int tmp_idx){
   float ds_dR = 999.0;

   if (tmp_idx>=NTwoMuonsTrack()) {
      cout<<"Index out of range!"<<endl;
      return ds_dR;
   }
   for (unsigned int ngen=0; ngen<NMCSignalParticles(); ngen++){
      float tmp_phi_dR = 999.0, tmp_pi_dR = 999.0;

      int track_idx = -1, mu1_idx = -1, mu2_idx = -1;

      TLorentzVector tmp_pi_p4(0,0,0,0);
      TLorentzVector tmp_phi_p4(0,0,0,0);

      if (fabs(MCSignalParticle_pdgid(ngen))==431){

         tmp_phi_dR = 999.0;
         tmp_pi_dR = 999.0;

         for (int nchild=0; nchild<MCSignalParticle_Nchilds(ngen); nchild++){

            //Find a pi that is a decay product of ds
            if(fabs(MCSignalParticle_childpdgid(ngen,nchild))==211){
               unsigned int tmp_track_idx = TwoMuonsTrackTrackIndex(tmp_idx).at(0);
               float dR = Track_P4(tmp_track_idx).DeltaR(MCSignalParticle_child_p4(ngen,nchild));
               if (dR<0.03 && dR<tmp_pi_dR) {
                  tmp_pi_dR = dR;
                  tmp_pi_p4 = MCSignalParticle_child_p4(ngen,nchild);
                  track_idx = tmp_track_idx;
               }
            }
            //Find a phi that is a decay product of Ds
            if (fabs(MCSignalParticle_childpdgid(ngen,nchild))==333){
               unsigned int tmp_mu1_idx = TwoMuonsTrackMuonIndices(tmp_idx).at(0);
               unsigned int tmp_mu2_idx = TwoMuonsTrackMuonIndices(tmp_idx).at(1);
               float dR = MCSignalParticle_child_p4(ngen,nchild).DeltaR(Muon_P4(tmp_mu1_idx)+Muon_P4(tmp_mu2_idx));
               if (dR<0.03 && dR<tmp_phi_dR) {
                  tmp_phi_dR = dR;
                  tmp_phi_p4 = MCSignalParticle_child_p4(ngen,nchild);
                  mu1_idx = tmp_mu1_idx;
                  mu2_idx = tmp_mu2_idx;
               }
            }
         }

         if (tmp_phi_dR<0.03 && tmp_pi_dR<0.03){
            float tmp_ds_dR = (tmp_phi_p4+tmp_pi_p4).DeltaR(Muon_P4(mu1_idx)+Muon_P4(mu2_idx)+Track_P4(track_idx));
            if (tmp_ds_dR<ds_dR) ds_dR = tmp_ds_dR;
         }
      }
   }
   return ds_dR;
}

//-------------------------- --------------------------------- -------------------

//-------------------------- Match GEN Ds to 2mu+trk candidate -------------------
// Return true if the Ds is a prompt one (doesn't come from B+, B0, Bs decays)
bool Ntuple_Controller::isPromptDs(){
   // Check if the candidate is a prompt Ds
   bool isPrompt = true;

   for (unsigned int ngen=0; ngen<NMCSignalParticles(); ++ngen){
      if (!(abs(MCSignalParticle_pdgid(ngen))==431)) continue;
      for (int nmom=0; nmom<NMCSignalParticleSources(ngen); ++nmom){
         if ( (abs(MCSignalParticle_Sourcepdgid(ngen,nmom))==511) || (abs(MCSignalParticle_Sourcepdgid(ngen,nmom))==521) || (abs(MCSignalParticle_Sourcepdgid(ngen,nmom))==531)) isPrompt = false;
      }
   }

   return isPrompt;
}
//-------------------------- --------------------------------- -------------------

// Get Mother PDGId of the Ds
int Ntuple_Controller::DsMotherPdgId(unsigned int index){
   int pdgid = 0;
   for (unsigned int ngen=0; ngen<NMCSignalParticles(); ++ngen){
      if (!(abs(MCSignalParticle_pdgid(ngen))==431)) continue;
      if (DsGenMatch(index)>0.03) continue;
      for (int nmom=0; nmom<NMCSignalParticleSources(ngen); ++nmom){
         pdgid = abs(MCSignalParticle_Sourcepdgid(ngen,nmom));
      }
   }
   return pdgid; 
}

//-------------------------- Match Tau3Mu candidates to Taus -------------------

float Ntuple_Controller::TauGenMatch(unsigned int tmp_idx){

   float tau_dR = 999.0;

   if (tmp_idx>=NThreeMuons()){
      return tau_dR;
   }

   int mu1_idx = ThreeMuonIndices(tmp_idx).at(0);
   int mu2_idx = ThreeMuonIndices(tmp_idx).at(1);
   int mu3_idx = ThreeMuonIndices(tmp_idx).at(2);

   TLorentzVector tau_p4 = Muon_P4(mu1_idx)+Muon_P4(mu2_idx)+Muon_P4(mu3_idx);

   for (int ngen=0; ngen<NMCTaus(); ++ngen){

      int gen_mu1_idx = -1, gen_mu2_idx = -1, gen_mu3_idx = -1;

      bool mu1_match = false;
      bool mu2_match = false;
      bool mu3_match = false;

      float tmp_dR = 999.0;

      for (int nchild=1; nchild<NMCTauDecayProducts(ngen); nchild++){

         if (Muon_P4(mu1_idx).DeltaR(MCTauandProd_p4(ngen,nchild))<0.03 && !mu1_match) {
            mu1_match = true;
            gen_mu1_idx = nchild;
         }
         if (Muon_P4(mu2_idx).DeltaR(MCTauandProd_p4(ngen,nchild))<0.03 && !mu2_match) {
            mu2_match = true;
            gen_mu2_idx = nchild;
         }
         if (Muon_P4(mu3_idx).DeltaR(MCTauandProd_p4(ngen,nchild))<0.03 && !mu3_match) {
            mu3_match = true;
            gen_mu3_idx = nchild;
         }
      }

      if (mu1_match && mu2_match && mu3_match) {
         TLorentzVector gen_tau_p4 = (MCTauandProd_p4(ngen,gen_mu1_idx)+MCTauandProd_p4(ngen,gen_mu2_idx)+MCTauandProd_p4(ngen,gen_mu3_idx));
         tmp_dR = gen_tau_p4.DeltaR(tau_p4);
         if (tmp_dR<tau_dR) tau_dR = tmp_dR;
      }
   }

   return tau_dR;
}

//-------------------------- --------------------------------- -------------------
std::vector<unsigned int> Ntuple_Controller::SortedPtMuons(std::vector<unsigned int> indices){

   std::vector<unsigned int> out;
   unsigned int i1,i2,i3;
   if(indices.size()!=3) {std::cout<< "Sorting muons by Pt failed, wrong vector size! "<< std::endl; return out;}
   double pt1 = Muon_P4(indices.at(0)).Pt();
   double pt2 = Muon_P4(indices.at(1)).Pt();
   double pt3 = Muon_P4(indices.at(2)).Pt();

   if(pt1>pt2)
   {
      if(pt2>pt3)
      {
         i1=indices.at(0); i2 = indices.at(1); i3 = indices.at(2);
      }
      else if(pt1>pt3)
      {
         i1=indices.at(0); i2 = indices.at(2); i3 = indices.at(1);
      }
      else
      {
         i1=indices.at(2); i2 = indices.at(0); i3 = indices.at(1);
      }
   }
   else
   {
      if(pt1>pt3)
      {
         i1=indices.at(1); i2 = indices.at(0); i3 = indices.at(2);
      }
      else if(pt2>pt3)
      {
         i1=indices.at(1); i2 = indices.at(2); i3 = indices.at(0);
      }
      else
      {
         i1=indices.at(2); i2 = indices.at(1); i3 = indices.at(0);
      }
   }
   out.push_back(i1);
   out.push_back(i2);
   out.push_back(i3);

   return out;
}



std::vector<unsigned int> Ntuple_Controller::SortedEtaMuons(std::vector<unsigned int> indices){

   std::vector<unsigned int> out;
   unsigned int i1,i2,i3;

   double eta1 = fabs(Muon_P4(indices.at(0)).Eta());
   double eta2 = fabs(Muon_P4(indices.at(1)).Eta());
   double eta3 = fabs(Muon_P4(indices.at(2)).Eta());

   if(eta1>eta2)
   {
      if(eta2>eta3)
      {
         i1=indices.at(0); i2 = indices.at(1); i3 = indices.at(2);
      }
      else if(eta1>eta3)
      {
         i1=indices.at(0); i2 = indices.at(2); i3 = indices.at(1);
      }
      else
      {
         i1=indices.at(2); i2 = indices.at(0); i3 = indices.at(1);
      }
   }
   else
   {
      if(eta1>eta3)
      {
         i1=indices.at(1); i2 = indices.at(0); i3 = indices.at(2);
      }
      else if(eta2>eta3)
      {
         i1=indices.at(1); i2 = indices.at(2); i3 = indices.at(0);
      }
      else
      {
         i1=indices.at(2); i2 = indices.at(1); i3 = indices.at(0);
      }
   }
   out.push_back(i1);
   out.push_back(i2);
   out.push_back(i3);

   return out;
}


TLorentzVector Ntuple_Controller::MatchedLV(std::vector<TLorentzVector> list, unsigned int index){

   TLorentzVector ToMatchLV = Muon_P4(index);
   TLorentzVector out(0,0,0,0);
   double dR(0.99);

   for(auto &i : list){
      if(ToMatchLV.DeltaR(i)< dR){
         dR = ToMatchLV.DeltaR(i);
         out = i;
      }
   }

   return out;
}



double Ntuple_Controller::TauMassResolution(std::vector<unsigned int>  indices, int type, bool UseRefited=true){ // type = 0 - only pt propagation; type = 1 - indluding direction


   std::vector<TLorentzVector > listRefited;
   unsigned int final_idx = 0;
   listRefited.push_back(Vertex_signal_KF_refittedTracksP4(final_idx,0));
   listRefited.push_back(Vertex_signal_KF_refittedTracksP4(final_idx,1));
   listRefited.push_back(Vertex_signal_KF_refittedTracksP4(final_idx,2));

   TLorentzVector Mu1(0,0,0,0);
   TLorentzVector Mu2(0,0,0,0);
   TLorentzVector Mu3(0,0,0,0);


   if(UseRefited){
      Mu1=MatchedLV(listRefited,indices.at(0));
      Mu2=MatchedLV(listRefited,indices.at(1));
      Mu3=MatchedLV(listRefited,indices.at(2));
   } else {

      Mu1 = Muon_P4(indices.at(0));
      Mu2 = Muon_P4(indices.at(1));
      Mu3 = Muon_P4(indices.at(2));

   }

   TLorentzVector Tau = Mu1 + Mu2 + Mu3;
   TLorentzVector Mu1_deltaPt(0,0,0,0);
   TLorentzVector Mu2_deltaPt(0,0,0,0);
   TLorentzVector Mu3_deltaPt(0,0,0,0);

   TLorentzVector Mu1_deltaEta(0,0,0,0);
   TLorentzVector Mu2_deltaEta(0,0,0,0);
   TLorentzVector Mu3_deltaEta(0,0,0,0);

   TLorentzVector Mu1_deltaPhi(0,0,0,0);
   TLorentzVector Mu2_deltaPhi(0,0,0,0);
   TLorentzVector Mu3_deltaPhi(0,0,0,0);


   Mu1_deltaPt.SetPtEtaPhiM(Mu1.Pt() + Muon_ptError(indices.at(0)),Mu1.Eta(), Mu1.Phi(), Mu1.M());
   Mu2_deltaPt.SetPtEtaPhiM(Mu2.Pt() + Muon_ptError(indices.at(1)),Mu2.Eta(), Mu2.Phi(), Mu2.M());
   Mu3_deltaPt.SetPtEtaPhiM(Mu3.Pt() + Muon_ptError(indices.at(2)),Mu3.Eta(), Mu3.Phi(), Mu3.M());


   Mu1_deltaEta.SetPtEtaPhiM(Mu1.Pt(),Mu1.Eta()  + Muon_etaError(indices.at(0)), Mu1.Phi(), Mu1.M());
   Mu2_deltaEta.SetPtEtaPhiM(Mu2.Pt(),Mu2.Eta()  + Muon_etaError(indices.at(1)), Mu2.Phi(), Mu2.M());
   Mu3_deltaEta.SetPtEtaPhiM(Mu3.Pt(),Mu3.Eta()  + Muon_etaError(indices.at(2)), Mu3.Phi(), Mu3.M());


   Mu1_deltaPhi.SetPtEtaPhiM(Mu1.Pt(), Mu1.Eta(), Mu1.Phi() + Muon_phiError(indices.at(0)), Mu1.M());
   Mu2_deltaPhi.SetPtEtaPhiM(Mu2.Pt(), Mu2.Eta(), Mu2.Phi() + Muon_phiError(indices.at(1)), Mu2.M());
   Mu3_deltaPhi.SetPtEtaPhiM(Mu3.Pt(), Mu3.Eta(), Mu3.Phi() + Muon_phiError(indices.at(2)), Mu3.M());


   TLorentzVector Tau_Mu1PtUp = Mu1_deltaPt + Mu2 + Mu3;
   TLorentzVector Tau_Mu2PtUp = Mu1 + Mu2_deltaPt + Mu3;
   TLorentzVector Tau_Mu3PtUp = Mu1 + Mu2 + Mu3_deltaPt;

   TLorentzVector Tau_Mu1EtaUp = Mu1_deltaEta + Mu2 + Mu3;
   TLorentzVector Tau_Mu2EtaUp = Mu1 + Mu2_deltaEta + Mu3;
   TLorentzVector Tau_Mu3EtaUp = Mu1 + Mu2 + Mu3_deltaEta;

   TLorentzVector Tau_Mu1PhiUp = Mu1_deltaPhi + Mu2 + Mu3;
   TLorentzVector Tau_Mu2PhiUp = Mu1 + Mu2_deltaPhi + Mu3;
   TLorentzVector Tau_Mu3PhiUp = Mu1 + Mu2 + Mu3_deltaPhi;


   double  deltaM_pt_1 = Tau_Mu1PtUp.M() - Tau.M();
   double  deltaM_pt_2 = Tau_Mu2PtUp.M() - Tau.M();
   double  deltaM_pt_3 = Tau_Mu3PtUp.M() - Tau.M();

   double  deltaM_Eta_1 = Tau_Mu1EtaUp.M() - Tau.M();
   double  deltaM_Eta_2 = Tau_Mu2EtaUp.M() - Tau.M();
   double  deltaM_Eta_3 = Tau_Mu3EtaUp.M() - Tau.M();

   double  deltaM_Phi_1 = Tau_Mu1PhiUp.M() - Tau.M();
   double  deltaM_Phi_2 = Tau_Mu2PhiUp.M() - Tau.M();
   double  deltaM_Phi_3 = Tau_Mu3PhiUp.M() - Tau.M();


   double deltaTauM(0);
   if(type == 0) 
      deltaTauM = sqrt(deltaM_pt_1*deltaM_pt_1 + deltaM_pt_2*deltaM_pt_2 + deltaM_pt_3*deltaM_pt_3);

   if(type == 1)
      deltaTauM = sqrt(deltaM_pt_1*deltaM_pt_1 + deltaM_pt_2*deltaM_pt_2 + deltaM_pt_3*deltaM_pt_3
            +deltaM_Eta_1*deltaM_Eta_1 + deltaM_Eta_2*deltaM_Eta_2 + deltaM_Eta_3*deltaM_Eta_3
            +deltaM_Phi_1*deltaM_Phi_1 + deltaM_Phi_2*deltaM_Phi_2 + deltaM_Phi_3*deltaM_Phi_3);


   return deltaTauM/Tau.M();

}




LorentzVectorParticle Ntuple_Controller::Tau3mu_LVP(unsigned int i){
  TMatrixT<double>    t3m_par(LorentzVectorParticle::NLorentzandVertexPar,1);
  TMatrixTSym<double> t3m_cov(LorentzVectorParticle::NLorentzandVertexPar);
  int l=0;
  if(Ntp->signalTau_lvp->at(i).size()==LorentzVectorParticle::NLorentzandVertexPar){
    for(int k=0; k<LorentzVectorParticle::NLorentzandVertexPar; k++){
      t3m_par(k,0)=Ntp->signalTau_lvp->at(i).at(k);
      for(int j=k; j<LorentzVectorParticle::NLorentzandVertexPar; j++){
        t3m_cov(k,j)=Ntp->signalTau_cov->at(i).at(l);
        t3m_cov(j,k)=Ntp->signalTau_cov->at(i).at(l);
        l++;
      }
    }
  }
  return LorentzVectorParticle(t3m_par,t3m_cov,Ntp->signalTau_pdgid->at(i),Ntp->signalTau_charge->at(i),Ntp->signalTau_B->at(i));
}




TLorentzVector Ntuple_Controller::Boost(TLorentzVector pB, TLorentzVector frame){// make the boost of vector pB to the "frame"
   TMatrixT<double> transform(4,4);
   TMatrixT<double> result(4,1);
   TVectorT<double> vec(4); 
   TVector3 b;
   if(frame.Vect().Mag()==0){  Logger(Logger::Warning)<< "Boost is not set, perfrom calculation in the Lab Frame  "<<std::endl; return pB;}
   if(frame.E()==0){ Logger(Logger::Warning)<<" Caution: Please check that you perform boost correctly!  " <<std::endl; return pB;} 
   else   b=frame.Vect()*(1/frame.E());
   vec(0)  = pB.E();    vec(1)  = pB.Px(); 
   vec(2)  = pB.Py();  vec(3)  = pB.Pz();
   double gamma  = 1/sqrt( 1 - b.Mag2());
   transform(0,0)=gamma; transform(0,1) =- gamma*b.X() ;  transform(0,2) =  - gamma*b.Y();  transform(0,3) = - gamma*b.Z(); 
   transform(1,0)=-gamma*b.X(); transform(1,1) =(1+ (gamma-1)*b.X()*b.X()/b.Mag2()) ;  transform(1,2) = ((gamma-1)*b.X()*b.Y()/b.Mag2());  transform(1,3) = ((gamma-1)*b.X()*b.Z()/b.Mag2());
   transform(2,0)=-gamma*b.Y(); transform(2,1) = ((gamma-1)*b.Y()*b.X()/b.Mag2());  transform(2,2) = (1 + (gamma-1)*b.Y()*b.Y()/b.Mag2());  transform(2,3) =  ((gamma-1)*b.Y()*b.Z()/b.Mag2()); 
   transform(3,0)=-gamma*b.Z(); transform(3,1) =((gamma-1)*b.Z()*b.X()/b.Mag2()) ;  transform(3,2) = ((gamma-1)*b.Z()*b.Y()/b.Mag2());  transform(3,3) = (1 + (gamma-1)*b.Z()*b.Z()/b.Mag2()); 
   result=transform*convertToMatrix(vec);
   return TLorentzVector(result(1,0), result(2,0) ,result(3,0), result(0,0));
}




TVector3 Ntuple_Controller::Rotate(TVector3 LVec, TVector3 Rot){
  TVector3 vec = LVec;
  vec.RotateZ(0.5*TMath::Pi() - Rot.Phi());  
  vec.RotateX(Rot.Theta());
  return vec;
}



TMatrixT<double>  Ntuple_Controller::convertToMatrix(TVectorT<double> V){
   TMatrixT<double> M(V.GetNrows(),1);
   for(int i=0; i < M.GetNrows(); i++){
      M(i,0)=V(i);
   } return M;
}


TMatrixTSym<float> Ntuple_Controller::SecondaryVertexCovariance(unsigned int i){
   TMatrixTSym<float> V_cov(3);
   int l=0;

   for(int j=0;j<3;j++){
      for(int k=j;k<3;k++){
         V_cov(j,k)=Ntp->SV_PosCovariance->at(i).at(l);
         V_cov(k,j)=Ntp->SV_PosCovariance->at(i).at(l);
         l++;
      }
   }
   return  V_cov;
}




TMatrixTSym<double> Ntuple_Controller::Vertex_Signal_KF_Covariance(unsigned int i, bool channel){

   unsigned int index = i + channel*NThreeMuons();
   TMatrixTSym<double> V_cov(3);
   int l=0;

   for(int j=0;j<3;j++){
      for(int k=j;k<3;k++){
         V_cov(j,k)=Ntp->Vertex_signal_KF_cov->at(index).at(l);
         V_cov(k,j)=Ntp->Vertex_signal_KF_cov->at(index).at(l);
         l++;
      }
   }
   return  V_cov;
}


TMatrixTSym<double> Ntuple_Controller::Vertex_2MuonsIsoTrack_KF_Covariance(unsigned int i, bool channel){

  unsigned int index = i + channel*NThreeMuons();
  TMatrixTSym<double> V_cov(3);
  int l=0;

  for(int j=0;j<3;j++){
    for(int k=j;k<3;k++){
      V_cov(j,k)=Ntp->Vertex_2MuonsIsoTrack_KF_cov->at(index).at(l);
      V_cov(k,j)=Ntp->Vertex_2MuonsIsoTrack_KF_cov->at(index).at(l);
      l++;
    }
  }
  return  V_cov;
}



TMatrixTSym<double> Ntuple_Controller::Vertex_PrimaryVertex_Covariance(unsigned int i, bool channel){
   unsigned int index = i + channel*NThreeMuons();
   TMatrixTSym<float> V_cov(3);
   int l=0;
   for(int j=0;j<3;j++){
      for(int k=j;k<3;k++){
         //if(j==k) V_cov(i,j)=pow(0.0001,2.0);
         V_cov(j,k)=Ntp->Vertex_MatchedRefitPrimaryVertex_covariance->at(index).at(l);
         V_cov(k,j)=Ntp->Vertex_MatchedRefitPrimaryVertex_covariance->at(index).at(l);
         l++;
      }
   }
   return V_cov;
}




TMatrixTSym<double> Ntuple_Controller::Vertex_HighestPt_PrimaryVertex_Covariance(){
   TMatrixTSym<float> V_cov(3);
   int l=0;
   for(int j=0;j<3;j++){
      for(int k=j;k<3;k++){
         //if(j==k) V_cov(i,j)=pow(0.0001,2.0);
         V_cov(j,k)=Ntp->Vertex_HighestPt_PrimaryVertex_covariance->at(0).at(l);
         V_cov(k,j)=Ntp->Vertex_HighestPt_PrimaryVertex_covariance->at(0).at(l);
         l++;
      }
   }
   return V_cov;
}




TVector3  Ntuple_Controller::SVPVDirection(TVector3 SV, TVector3 PV){
   return SV-PV;
}


double Ntuple_Controller::FlightLength_significance(TVector3 pv,TMatrixTSym<double> PVcov, TVector3 sv, TMatrixTSym<double> SVcov ){
   TVector3 SVPV = sv - pv;
   TVectorF FD;
   FD.ResizeTo(3);
   FD(0) = SVPV.X();
   FD(1) = SVPV.Y();
   FD(2) = SVPV.Z();

   TMatrixT<double> PVcv;
   PVcv.ResizeTo(3,3);
   for(int nr =0; nr<PVcov.GetNrows(); nr++){
      for(int nc =0; nc<PVcov.GetNcols(); nc++){
         PVcv(nr,nc) = PVcov(nr,nc);
      }
   }
   TMatrixT<double> SVcv;
   SVcv.ResizeTo(3,3);
   for(int nr =0; nr<SVcov.GetNrows(); nr++){
      for(int nc =0; nc<SVcov.GetNcols(); nc++){
         SVcv(nr,nc) = SVcov(nr,nc);
      }
   }

   TMatrixT<double> SVPVMatrix(3,1);
   for(int i=0; i<SVPVMatrix.GetNrows();i++){
      SVPVMatrix(i,0)=FD(i);
   }

   TMatrixT<double> SVPVMatrixT=SVPVMatrix;
   SVPVMatrixT.T();

   TMatrixT<double> lambda2 = SVPVMatrixT*(SVcv + PVcv)*SVPVMatrix;
   double sigmaabs = sqrt(lambda2(0,0))/SVPV.Mag();
   double sign = SVPV.Mag()/sigmaabs;

   return sign;
}

double Ntuple_Controller::TransverseFlightLength_significance(TVector3 pv, TMatrixTSym<double> PVcov, TVector3 sv, TMatrixTSym<double> SVcov ){
   TVector3 SVPV = sv - pv;
   TVectorF FD;
   FD.ResizeTo(3);
   FD(0) = SVPV.X();
   FD(1) = SVPV.Y();
   FD(2) = SVPV.Z();

   TMatrixT<double> PVcv;
   PVcv.ResizeTo(2,2);
   for(int nr =0; nr<PVcov.GetNrows()-1; nr++){
      for(int nc =0; nc<PVcov.GetNcols()-1; nc++){
         PVcv(nr,nc) = PVcov(nr,nc);
      }
   }
   TMatrixT<double> SVcv;
   SVcv.ResizeTo(2,2);
   for(int nr =0; nr<SVcov.GetNrows()-1; nr++){
      for(int nc =0; nc<SVcov.GetNcols()-1; nc++){
         SVcv(nr,nc) = SVcov(nr,nc);
      }
   }

   TMatrixT<double> SVPVMatrix(2,1);
   for(int i=0; i<SVPVMatrix.GetNrows()-1;i++){
      SVPVMatrix(i,0)=FD(i);
   }

   TMatrixT<double> SVPVMatrixT=SVPVMatrix;
   SVPVMatrixT.T();

   TMatrixT<double> lambda2 = SVPVMatrixT*(SVcv + PVcv)*SVPVMatrix;
   double sigmaabs = sqrt(lambda2(0,0))/SVPV.Mag();
   double sign = SVPV.Mag()/sigmaabs;

   return sign;

}



TMatrixD Ntuple_Controller::EigenVectors(TMatrixTSym<double> M){

  TMatrixDSymEigen eigen_matrix(M);
  return eigen_matrix.GetEigenVectors();

}

TVectorD Ntuple_Controller::EigenValues(TMatrixTSym<double> M){
  TMatrixDSymEigen eigen_matrix(M);
  return eigen_matrix.GetEigenValues();
}



//////////////////////////////////////////////////////////////////////
//   If the covariance matrix contains the negative eigen values then 
//   the dioganal elements are recursively increased  by 1%

TMatrixTSym<double>  Ntuple_Controller::RegulariseCovariance(TMatrixTSym<double>  M, double coef){

  TMatrixTSym<double>  M_infl = M;
  for(unsigned int i=0; i<TrackParticle::NHelixPar;i++){
    M_infl(i,i)*=coef;
  }

  TVectorD eigen_val = EigenValues(M_infl);
  for(unsigned int i=0; i<eigen_val.GetNrows(); i++){
    if(eigen_val(i) < 0){
      coef*=1.01;
      return RegulariseCovariance(M_infl,coef);
    }
  }
  return M_infl;
}



TMatrixTSym<double>  
Ntuple_Controller::InvertMatrix(TMatrixTSym<double>   M){
  TDecompBK Inverter(M);
  return  Inverter.Invert();
}





std::vector<unsigned int> Ntuple_Controller::SortedChargeMuons(std::vector<unsigned int> indices){

   std::vector<unsigned int> out;

   double q1 = Muon_charge(indices[0]);
   double q2 = Muon_charge(indices[1]);
   double q3 = Muon_charge(indices[2]);

   double pt1 = Muon_P4(indices[0]).Pt();
   double pt2 = Muon_P4(indices[1]).Pt();
   double pt3 = Muon_P4(indices[2]).Pt();

   unsigned int i1,i2,i3;

   if (q1==q2){
      i1 = indices[2];
      i2 = pt1>pt2?indices[0]:indices[1];
      i3 = pt1>pt2?indices[1]:indices[0];
   }
   else if (q2==q3){
      i1 = indices[0];
      i2 = pt2>pt3?indices[1]:indices[2];
      i3 = pt2>pt3?indices[2]:indices[1];
   }
   else if (q1==q3){
      i1 = indices[1];
      i2 = pt1>pt3?indices[0]:indices[2];
      i3 = pt1>pt3?indices[2]:indices[0];
   }

   if (abs(q1+q2+q3)>1.1) Logger(Logger::Warning)<< "Sum of charges is greater than 1"<<endl;

   out.push_back(i1);
   out.push_back(i2);
   out.push_back(i3);

   return out;
}



TLorentzVector Ntuple_Controller::matchToTruthTauDecay(TLorentzVector vector){
   TLorentzVector out(0,0,0,0);
   double dr(0.3);
   if(NMCTaus()==0)
   {
      Logger(Logger::Warning) << "No truth tau leptons in this event found; return TLorentzVector(0,0,0,0) " << std::endl; return out;
   }
   for(int i=0; i < NMCTaus(); i++ ){
      for(int j =0; j < NMCTauDecayProducts(i); j++){
         if(MCTauandProd_p4(i,j).DeltaR(vector) < dr){
            out=MCTauandProd_p4(i,j);
            dr = MCTauandProd_p4(i,j).DeltaR(vector);
         }
      }
   }
   return out;
}



int Ntuple_Controller::matchTruth(TLorentzVector tvector){
   double testdr=0.3;
   int pdgid = 0;
   for(unsigned int i=0;i<NMCParticles();i++){
      if(MCParticle_p4(i).Pt()>0.){
         if(tvector.DeltaR(MCParticle_p4(i))<testdr){
            testdr = tvector.DeltaR(MCParticle_p4(i));
            pdgid = MCParticle_pdgid(i);
         }
      }
   }
   return pdgid;
}

int Ntuple_Controller::getMatchTruthIndex(TLorentzVector tvector){
  int index = -9;
  double dr=999.;
  for(unsigned int i=0;i<NMCParticles();i++){
    if(tvector.DeltaR(MCParticle_p4(i))<dr){
      index = i;
      dr = tvector.DeltaR(MCParticle_p4(i));
    }
  }
  return index;
}

TLorentzVector Ntuple_Controller::KaonTrack_P4(unsigned int index){
   TLorentzVector p(Track_P4(index).Px(),
                    Track_P4(index).Py(),
                    Track_P4(index).Pz(),
                    sqrt(pow(Track_P4(index).P(),2)+pow(PDG_Var::Kp_mass(),2)));
   return p;
}

TLorentzVector Ntuple_Controller::MuonToKaon(unsigned int index){
  TLorentzVector p(Muon_P4(index).Px(),
                   Muon_P4(index).Py(),
                   Muon_P4(index).Pz(),
                   sqrt(pow(Muon_P4(index).P(),2)+pow(PDG_Var::Kp_mass(),2)));
  return p;
}

TLorentzVector Ntuple_Controller::MuonToPion(unsigned int index){
  TLorentzVector p(Muon_P4(index).Px(),
                   Muon_P4(index).Py(),
                   Muon_P4(index).Pz(),
                   sqrt(pow(Muon_P4(index).P(),2)+pow(PDG_Var::Pi_mass(),2)));
  return p;
}

std::vector<int> Ntuple_Controller::MuonCustomID(unsigned int MuonIndex){

   std::vector<int> out;


   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::CutBasedIdLoose))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::SoftCutBasedId))out.push_back(1);
   else out.push_back(0);


   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::CutBasedIdMedium))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::CutBasedIdTight))out.push_back(1);
   else out.push_back(0);


   if(Muon_isTrackerMuon(MuonIndex))out.push_back(1);
   else out.push_back(0);

   if(Muon_isGlobalMuon(MuonIndex))out.push_back(1);
   else out.push_back(0);


   return out;
}


std::vector<int> Ntuple_Controller::MuonStandardSelectorBitMask(unsigned int MuonIndex){

   std::vector<int> out;
   out.push_back(1);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::CutBasedIdLoose))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::CutBasedIdMedium))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::CutBasedIdMediumPrompt))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::CutBasedIdTight))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::CutBasedIdGlobalHighPt))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::CutBasedIdTrkHighPt))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::PFIsoVeryLoose))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::PFIsoLoose))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::PFIsoMedium))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::PFIsoTight))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::PFIsoVeryTight))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::TkIsoLoose))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::TkIsoTight))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::SoftCutBasedId))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::SoftMvaId))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::MvaLoose))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::MvaMedium))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::MvaTight))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::MiniIsoLoose))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::MiniIsoMedium))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::MiniIsoTight))out.push_back(1);
   else out.push_back(0);

   if(CHECK_BIT(Muon_StandardSelection(MuonIndex),MuonStandardSelectors::MiniIsoVeryTight))out.push_back(1);
   else out.push_back(0);

   return out;
}


TString Ntuple_Controller::WhichEra(int year){
   TString out("");
   if(!isData()) {
      Logger(Logger::Info)<< "You are currently running over MC, no era returned;  "<< std::endl; 
      return out;
   }
   if(year == 2017){
      if(RunNumber() >=297047 and RunNumber() <=299329) return out = "RunB";

      if(RunNumber() >=299368 and RunNumber() <=302029) return out = "RunC";

      if(RunNumber() >=302031 and RunNumber() <=302663) return out = "RunD";

      if(RunNumber() >=303824 and RunNumber() <=304797) return out = "RunE";

      if(RunNumber() >=305040 and RunNumber() <=306460) return out = "RunF";
   }

   if(year == 2018){
      if (RunNumber() >=316239 and RunNumber() <=316944) return out = "RunA";

      if (RunNumber() >=318070 and RunNumber() <=319310) return out = "RunB";

      if (RunNumber() >=319449 and RunNumber() <=319756) return out = "RunC";

      if (RunNumber() >=320500 and RunNumber() <=325175) return out = "RunD";
   }
   return out;
}

//// draw decay chain
void Ntuple_Controller::printMCDecayChainOfEvent(bool printStatus, bool printPt, bool printEtaPhi, bool printQCD){
   Logger(Logger::Info) << "=== Draw MC decay chain of event ===" << std::endl;
   for(unsigned int par = 0; par < NMCParticles(); par++){
      if ( !MCParticle_hasMother(par) )
         printMCDecayChainOfMother(par, printStatus, printPt, printEtaPhi, printQCD);
   }
}

void Ntuple_Controller::printMCDecayChainOfParticle(unsigned int index, bool printStatus, bool printPt, bool printEtaPhi, bool printQCD){
   Logger(Logger::Info) << "=== Draw MC decay chain of particle ===" << std::endl;
   printMCDecayChainOfMother(index, printStatus, printPt, printEtaPhi, printQCD);
}
void Ntuple_Controller::printMCDecayChainOfMother(unsigned int par, bool printStatus, bool printPt, bool printEtaPhi, bool printQCD){
   Logger(Logger::Info) << "Draw decay chain for mother particle at index " << par << " :" << std::endl;
   printMCDecayChain(par, 0, printStatus, printPt, printEtaPhi, printQCD);
}
void Ntuple_Controller::printMCDecayChain(unsigned int par, unsigned int level, bool printStatus, bool printPt, bool printEtaPhi, bool printQCD){
   std::ostringstream out;
   for(unsigned int l = 0; l < level; l++) out << "    ";
   out << MCParticleToString(par, printStatus, printPt, printEtaPhi);
   if (!printQCD && MCParticle_pdgid(par) == 92) {out << "    QCD stuff truncated ...\n"; std::cout << out.str(); return;}
   out << "\n";
   std::cout << out.str();
   for (unsigned int i_dau = 0; i_dau < MCParticle_childidx(par).size(); i_dau++){
      if (MCParticle_childidx(par).at(i_dau) != (int)par) // skip cases where particles are daughters of themselves
         printMCDecayChain(MCParticle_childidx(par).at(i_dau), level+1, printStatus, printPt, printEtaPhi, printQCD);
   }
}

std::string Ntuple_Controller::MCParticleToString(unsigned int par, bool printStatus, bool printPt, bool printEtaPhi){
   std::ostringstream out;
   out << "+-> ";
   out << PDGInfo::pdgIdToName( MCParticle_pdgid(par) );
   if (printStatus) out << " (status = " << MCParticle_status(par) << ", idx = " << par << ")";
   if (printPt || printEtaPhi) out <<  " [";
   if (printPt) out << "pT = " << MCParticle_p4(par).Pt() << "GeV";
   if (printEtaPhi) out << " eta = " << MCParticle_p4(par).Eta() << " phi = " << MCParticle_p4(par).Phi();
   if (printPt || printEtaPhi) out << "]";
   return out.str();
}
