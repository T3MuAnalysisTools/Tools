#include "Selection_Factory.h"
//#include "Logger.h"

#include "Example.h"
#ifdef USE_cherepanov
#include "cherepanov/MyTest.h"
#include "cherepanov/MCStudy.h"
#include "cherepanov/ThreeMu.h"
#include "cherepanov/DsToPhiPi.h"
#include "cherepanov/CommonSelector.h"
#include "cherepanov/SyncSignal.h"
#include "cherepanov/SyncDsPhiPi.h"
#include "cherepanov/FillTMVATrees.h"
#include "cherepanov/MakeMVATree.h"
#include "cherepanov/MakeMVACategoryII.h"
#include "cherepanov/SignalSelector.h"
#include "cherepanov/CategoryIISelector.h"
#include "cherepanov/BackgroundSelector.h"
#include "cherepanov/MuIDStudy.h"
#include "cherepanov/MinBiasSelector.h"
#include "cherepanov/Isolation.h"
#include "cherepanov/MCBackgroundStudy.h"
#include "cherepanov/DsPhiPeak.h"
#include "cherepanov/VertexFits.h"
#include "cherepanov/DebugFit.h"
#include "cherepanov/SignalVertexSelector.h"
#include "cherepanov/SimpleTauSelector.h"
#include "cherepanov/ZTau3MuTauh.h"
#include "cherepanov/ZTau3MuTaue.h"
#include "cherepanov/ZTau3MuTaumu.h"
#endif


#ifdef USE_joshi
#include "joshi/DsPhiPiTree.h"
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
#include "madhu/MCBackgroundStudy.h"
#include "madhu/ThreeMuonDecay.h"
#include "madhu/CommonSelector.h"
#include "madhu/BDTSelector.h"
#include "madhu/MakeMVATree.h"
#include "madhu/SignalVertexSelector.h"
#include "madhu/DebugFit.h"
#include "madhu/SimpleTauSelector.h"

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
  //  else if(Analysis.Contains("mcstudy"))s=new MCStudy(Analysis,UncertType);
  else if(Analysis.Contains("threemu"))s=new ThreeMu(Analysis,UncertType);
  //else if(Analysis.Contains("dstophipi"))s=new DsToPhiPi(Analysis,UncertType);
  //else if(Analysis.Contains("syncsignal"))s=new SyncSignal(Analysis,UncertType);
  //else if(Analysis.Contains("syncdsphipi"))s=new SyncDsPhiPi(Analysis,UncertType);
  else if(Analysis.Contains("filltmvatrees"))s=new FillTMVATrees(Analysis,UncertType);
  else if(Analysis.Contains("makemvatree"))s=new MakeMVATree(Analysis,UncertType);
  else if(Analysis.Contains("makemvacategoryII"))s=new MakeMVACategoryII(Analysis,UncertType);
  //else if(Analysis.Contains("signalselector"))s=new SignalSelector(Analysis,UncertType);
  //else if(Analysis.Contains("categoryiiselector"))s=new CategoryIISelector(Analysis,UncertType);
  //else if(Analysis.Contains("backgroundselector"))s=new BackgroundSelector(Analysis,UncertType);
  //else if(Analysis.Contains("muidstudy"))s=new MuIDStudy(Analysis,UncertType);
  //else if(Analysis.Contains("isolation"))s=new Isolation(Analysis,UncertType);
  else if(Analysis.Contains("mcbackgroundstudy"))s=new MCBackgroundStudy(Analysis,UncertType);
  //  else if(Analysis.Contains("dsphipeak"))s=new DsPhiPeak(Analysis,UncertType);
  else if(Analysis.Contains("commonselector"))s=new CommonSelector(Analysis,UncertType);
  else if(Analysis.Contains("vertexfits"))s=new VertexFits(Analysis,UncertType);
  else if(Analysis.Contains("debugfit"))s=new DebugFit(Analysis,UncertType);
  else if(Analysis.Contains("signalvertexselector"))s=new SignalVertexSelector(Analysis,UncertType);
  else if(Analysis.Contains("simpletauselector"))s=new SimpleTauSelector(Analysis,UncertType);
  else if(Analysis.Contains("ztau3mutauh"))s=new ZTau3MuTauh(Analysis,UncertType);
  else if(Analysis.Contains("ztau3mutaue"))s=new ZTau3MuTaue(Analysis,UncertType);
  else if(Analysis.Contains("ztau3mutaumu"))s=new ZTau3MuTaumu(Analysis,UncertType);




#endif

#ifdef USE_joshi
  else if(Analysis.Contains("dsphipitree"))s=new MuonPionTree(Analysis, UncertType);
#endif


#ifdef USE_wang
  else if(Analysis.Contains("mytest"))s=new MyTest(Analysis,UncertType);
#endif


#ifdef USE_madhu
  else if(Analysis.Contains("mytest"))s=new MyTest(Analysis,UncertType);
  else if(Analysis.Contains("mcbackgroundstudy"))s=new MCBackgroundStudy(Analysis,UncertType);
  else if(Analysis.Contains("threemuondecay"))s=new ThreeMuonDecay(Analysis,UncertType);
  else if(Analysis.Contains("commonselector"))s=new CommonSelector(Analysis,UncertType);
  else if(Analysis.Contains("bdtselector"))s=new BDTSelector(Analysis,UncertType);
  else if(Analysis.Contains("makemvatree"))s=new MakeMVATree(Analysis,UncertType);
  else if(Analysis.Contains("signalvertexselector"))s=new SignalVertexSelector(Analysis,UncertType);
  else if(Analysis.Contains("debugfit"))s=new DebugFit(Analysis,UncertType);
  else if(Analysis.Contains("simpletauselector"))s=new SimpleTauSelector(Analysis,UncertType);
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
