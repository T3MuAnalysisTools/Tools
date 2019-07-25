#ifndef SyncDsPhiPi_h
#define SyncDsPhiPi_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class SyncDsPhiPi : public Selection {

 public:
  SyncDsPhiPi(TString Name_, TString id_);
  virtual ~SyncDsPhiPi();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {HLTOk=0,is2MuTrk,PhiMassCut,NCuts}; 


 protected:
  virtual void doEvent();  
  virtual void Store_ExtraDist();

  TFile * file;
  TTree * Sync_tree;



 private:
 // Selection Variables
  // Initializhere your analysis histograms

  double sync_pt_1;
  double sync_pt_2;
  double sync_pt_3;

  double sync_eta_1;
  double sync_eta_2;
  double sync_eta_3;

  double muon_1_isGlob;
  double muon_2_isGlob;

  double muon_1_isTrack;
  double muon_2_isTrack;

  double phi_mass;
  double ds_mass;

  double evt,run,lumi;

  double sync_DsPhiPiVtx_x;
  double sync_DsPhiPiVtx_y;
  double sync_DsPhiPiVtx_z;
  double sync_DsPhiPiVtx_Chi2;



};
#endif
