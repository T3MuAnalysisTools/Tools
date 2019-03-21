//Ntuple_Controller.cxx IMPLEMENTATION FILE
#include "Ntuple_Controller.h"
#include "Parameters.h"
#include <tuple>
#include "Logger.h"
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


std::vector<unsigned int> Ntuple_Controller::SortedPtMuons(std::vector<unsigned int> indices){
  
  std::vector<unsigned int> out;
  unsigned int i1,i2,i3;

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

