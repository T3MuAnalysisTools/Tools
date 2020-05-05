#include "Selection_Factory.h"
#include "Logger.h"

#include "Example.h"
#ifdef USE_cherepanov
#include "cherepanov/MyTest.h"
#include "cherepanov/MCStudy.h"
#include "cherepanov/ThreeMu.h"
#include "cherepanov/DsToPhiPi.h"

#include "cherepanov/SyncSignal.h"
#include "cherepanov/SyncDsPhiPi.h"
#include "cherepanov/FillTMVATrees.h"
#include "cherepanov/SignalSelector.h"
#include "cherepanov/BackgroundSelector.h"
#include "cherepanov/MuIDStudy.h"
#include "cherepanov/MinBiasSelector.h"
#include "cherepanov/Isolation.h"

#endif


#ifdef USE_joshi
#include "joshi/MyTest.h"
//#include "joshi/Validation.h"
#include "joshi/DimuTrk.h"
#include "joshi/ThreeMu.h"
#include "joshi/TMVASignal.h"
#include "joshi/AnalysisWithTMVA.h"
#include "joshi/MCEfficiency.h"
#include "joshi/FillMVATree.h"
#include "joshi/FillMVATree_ThreeGlobal_TrackerOnly.h"
#include "joshi/FillMVATree_ThreeGlobal.h"
#include "joshi/FillMVATree_TwoGlobalTracker_TrackerOnly.h"
#include "joshi/FillMVATree_TwoGlobalTracker.h"
#include "joshi/DsPhiPeak.h"
#include "joshi/SignalSelector.h"
#include "joshi/MuonPionTree.h"
#endif


#ifdef USE_wang
#include "wang/MyTest.h"

#endif


#ifdef USE_menendez
#include "menendez/MyTest.h"
#include "menendez/DsToPhiPi.h"
#include "menendez/SyncDsPhiPi.h"
#endif


#ifdef USE_madhu
#include "madhu/MyTest.h"

#endif



Selection_Factory::Selection_Factory(){
}

Selection_Factory::~Selection_Factory(){
}

Selection_Base* Selection_Factory::Factory(TString Analysis, TString UncertType,int mode, int runtype, double lumi){
  Selection_Base* s;
  Analysis.ToLower();

  // ensuring code will compile independently of user code
  // WARNING: be aware of the consequences of "Contains". Make sure that Class "foo" is put after "foobar".
  if(Analysis.Contains("example"))s=new Example(Analysis,UncertType);
#ifdef USE_cherepanov
  else if(Analysis.Contains("mytest"))s=new MyTest(Analysis,UncertType);
  else if(Analysis.Contains("mcstudy"))s=new MCStudy(Analysis,UncertType);
  else if(Analysis.Contains("threemu"))s=new ThreeMu(Analysis,UncertType);
  else if(Analysis.Contains("dstophipi"))s=new DsToPhiPi(Analysis,UncertType);
  else if(Analysis.Contains("syncsignal"))s=new SyncSignal(Analysis,UncertType);
  else if(Analysis.Contains("syncdsphipi"))s=new SyncDsPhiPi(Analysis,UncertType);
  else if(Analysis.Contains("filltmvatrees"))s=new FillTMVATrees(Analysis,UncertType);
  else if(Analysis.Contains("signalselector"))s=new SignalSelector(Analysis,UncertType);
  else if(Analysis.Contains("backgroundselector"))s=new BackgroundSelector(Analysis,UncertType);
  else if(Analysis.Contains("muidstudy"))s=new MuIDStudy(Analysis,UncertType);
  else if(Analysis.Contains("isolation"))s=new Isolation(Analysis,UncertType);



#endif

#ifdef USE_joshi
  else if(Analysis.Contains("mytest"))s=new MyTest(Analysis,UncertType);
  //  else if(Analysis.Contains("validation"))s=new Validation(Analysis,UncertType);
  else if(Analysis.Contains("dimutrk"))s=new DimuTrk(Analysis,UncertType);
  else if(Analysis.Contains("threemu"))s=new ThreeMu(Analysis,UncertType); //  three mu id cuts
  else if(Analysis.Contains("mcefficiency"))s=new MCEfficiency(Analysis,UncertType);
  else if(Analysis.Contains("analysiswithtmva"))s=new AnalysisWithTMVA(Analysis,UncertType);
  else if(Analysis.Contains("fillmvatree_threeglobal_trackeronly"))s=new FillMVATree_ThreeGlobal_TrackerOnly(Analysis,UncertType);
  else if(Analysis.Contains("fillmvatree_threeglobal"))s=new FillMVATree_ThreeGlobal(Analysis,UncertType);
  else if(Analysis.Contains("fillmvatree_twoglobaltracker_trackeronly"))s=new FillMVATree_TwoGlobalTracker_TrackerOnly(Analysis,UncertType);
  else if(Analysis.Contains("fillmvatree_twoglobaltracker"))s=new FillMVATree_TwoGlobalTracker(Analysis,UncertType);
  else if(Analysis.Contains("fillmvatree"))s=new FillMVATree(Analysis,UncertType);
  else if(Analysis.Contains("signalselector"))s = new SignalSelector(Analysis,UncertType); 
  else if(Analysis.Contains("dsphipeak"))s = new DsPhiPeak(Analysis,UncertType); 
  else if(Analysis.Contains("SignalSelector"))s=new SignalSelector(Analysis, UncertType);
  else if(Analysis.Contains("muonpiontree"))s=new MuonPionTree(Analysis, UncertType);
  
#endif


#ifdef USE_wang
  else if(Analysis.Contains("mytest"))s=new MyTest(Analysis,UncertType);
#endif


#ifdef USE_madhu
  else if(Analysis.Contains("mytest"))s=new MyTest(Analysis,UncertType);
#endif


#ifdef USE_menendez
  else if(Analysis.Contains("mytest"))s=new MyTest(Analysis,UncertType);
  else if(Analysis.Contains("dstophipi"))s=new DsToPhiPi(Analysis,UncertType);
  else if(Analysis.Contains("syncdsphipi"))s=new SyncDsPhiPi(Analysis,UncertType);
#endif

  else{
	Logger(Logger::Error)<< "Invalid Analysis type \"" << Analysis << "\". Using default <Example.h> " << std::endl;
    s=new Example(Analysis,UncertType);
  }
  s->SetMode(mode);
  s->SetRunType(runtype);
  s->SetLumi(lumi);
  s->Configure();
  return s;
}
