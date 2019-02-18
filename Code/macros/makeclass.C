
void makeclass(TString file){

  TFile *_file0 = TFile::Open(file);
  _file0->Cd("T3MTree");
  TTree* treePtr = (TTree*)_file0 ->Get("T3MTree/t3mtree");
  //   TChain c("HTauTauTree");
  //   c.Add(file);
   treePtr->MakeClass("NtupleReader");
}
