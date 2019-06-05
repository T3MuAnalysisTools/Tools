#include "Selection_Factory.h"
#include "Logger.h"

#include "Example.h"
#ifdef USE_cherepanov
#include "cherepanov/MyTest.h"
#include "cherepanov/Validation.h"
#include "cherepanov/MCStudy.h"
#include "cherepanov/ThreeMu.h"
#include "cherepanov/TwoMuTrack.h"
#include "cherepanov/TMVATree.h"
#include "cherepanov/AnalysisWithTMVA.h"
#endif


#ifdef USE_joshi
#include "joshi/MyTest.h"
//#include "joshi/Validation.h"
#include "joshi/DimuTrk.h"
#include "joshi/ThreeMu.h"
#include "joshi/TMVASignal.h"
#include "joshi/AnalysisWithTMVA.h"
#include "joshi/MCEfficiency.h"
#include "joshi/twoGlbOneLooseMuIdCuts.h"
#include "joshi/looseMuIdCuts.h"
#include "joshi/trkMuIdCuts.h"
#include "joshi/glbMuIdCuts.h"
#include "joshi/softMuIdCuts.h"
#endif


#ifdef USE_wang
#include "wang/MyTest.h"

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
  else if(Analysis.Contains("validation"))s=new Validation(Analysis,UncertType);
  else if(Analysis.Contains("mcstudy"))s=new MCStudy(Analysis,UncertType);
  else if(Analysis.Contains("threemu"))s=new ThreeMu(Analysis,UncertType);
  else if(Analysis.Contains("twomutrack"))s=new TwoMuTrack(Analysis,UncertType);
  else if(Analysis.Contains("tmvatree"))s=new TMVATree(Analysis,UncertType);
  else if(Analysis.Contains("analysiswithtmva"))s=new AnalysisWithTMVA(Analysis,UncertType);



#endif

#ifdef USE_joshi
  else if(Analysis.Contains("mytest"))s=new MyTest(Analysis,UncertType);
  //  else if(Analysis.Contains("validation"))s=new Validation(Analysis,UncertType);
  else if(Analysis.Contains("dimutrk"))s=new DimuTrk(Analysis,UncertType);
  else if(Analysis.Contains("threemu"))s=new ThreeMu(Analysis,UncertType); // two global and one tracker id cuts
  else if(Analysis.Contains("mcefficiency"))s=new MCEfficiency(Analysis,UncertType);
  else if(Analysis.Contains("tmvasignal"))s=new TMVASignal(Analysis,UncertType);
  else if(Analysis.Contains("analysiswithtmva"))s=new AnalysisWithTMVA(Analysis,UncertType);
  else if(Analysis.Contains("glbMuIdCuts"))s=new AnalysisWithTMVA(Analysis,UncertType);
  else if(Analysis.Contains("looseMuIdCuts"))s=new AnalysisWithTMVA(Analysis,UncertType);
  else if(Analysis.Contains("softMuIdCuts"))s=new AnalysisWithTMVA(Analysis,UncertType);
  else if(Analysis.Contains("trkMuIdCuts"))s=new AnalysisWithTMVA(Analysis,UncertType);
  else if(Analysis.Contains("twoGlbOnLooseMuIdCuts"))s=new AnalysisWithTMVA(Analysis,UncertType);

#endif


#ifdef USE_wang
  else if(Analysis.Contains("mytest"))s=new MyTest(Analysis,UncertType);
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
