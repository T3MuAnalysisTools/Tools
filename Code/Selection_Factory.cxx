#include "Selection_Factory.h"
#include "Logger.h"

#include "Example.h"
#ifdef USE_cherepanov
#include "cherepanov/MyTest.h"
#include "cherepanov/Validation.h"
#include "cherepanov/MCStudy.h"
#include "cherepanov/ThreeMu.h"
#include "cherepanov/DsToPhiPi.h"

#include "cherepanov/SignalCategories.h"

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
#include "joshi/TripleMuIsolation.h"

#include "joshi/TMVATree_1.h"
#include "joshi/TMVATree_2.h"
#include "joshi/TMVATree_3.h"

#include "joshi/TMVA_Signal_glbMuon.h"
#include "joshi/TMVA_Signal_looseMuon.h"
#include "joshi/TMVA_Signal_softMuon.h"
#include "joshi/TMVA_Signal_trkMuon.h"
#include "joshi/TMVA_Signal_twoGlbOneLooseMuon.h"

#include "joshi/TMVA_Tree_1_glbMuon.h"
#include "joshi/TMVA_Tree_1_looseMuon.h"
#include "joshi/TMVA_Tree_1_softMuon.h"
#include "joshi/TMVA_Tree_1_trkMuon.h"
#include "joshi/TMVA_Tree_1_twoGlbOneLooseMuon.h"

#include "joshi/TMVA_Tree_2_glbMuon.h"
#include "joshi/TMVA_Tree_2_looseMuon.h"
#include "joshi/TMVA_Tree_2_softMuon.h"
#include "joshi/TMVA_Tree_2_trkMuon.h"
#include "joshi/TMVA_Tree_2_twoGlbOneLooseMuon.h"

#include "joshi/TMVA_Tree_3_glbMuon.h"
#include "joshi/TMVA_Tree_3_looseMuon.h"
#include "joshi/TMVA_Tree_3_softMuon.h"
#include "joshi/TMVA_Tree_3_trkMuon.h"
#include "joshi/TMVA_Tree_3_twoGlbOneLooseMuon.h"
#endif


#ifdef USE_wang
#include "wang/MyTest.h"

#endif



#ifdef USE_menendez
#include "menendez/MyTest.h"
#include "menendez/DsToPhiPi.h"
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
  else if(Analysis.Contains("dstophipi"))s=new DsToPhiPi(Analysis,UncertType);
  else if(Analysis.Contains("signalcategories"))s=new SignalCategories(Analysis,UncertType);



#endif

#ifdef USE_joshi
  else if(Analysis.Contains("mytest"))s=new MyTest(Analysis,UncertType);
  //  else if(Analysis.Contains("validation"))s=new Validation(Analysis,UncertType);
  else if(Analysis.Contains("dimutrk"))s=new DimuTrk(Analysis,UncertType);
  else if(Analysis.Contains("threemu"))s=new ThreeMu(Analysis,UncertType); //  three mu id cuts
  else if(Analysis.Contains("triplemuisolation"))s=new TripleMuIsolation(Analysis,UncertType); // calorimeter isolation variables 
  else if(Analysis.Contains("mcefficiency"))s=new MCEfficiency(Analysis,UncertType);
  else if(Analysis.Contains("analysiswithtmva"))s=new AnalysisWithTMVA(Analysis,UncertType);
  else if(Analysis.Contains("glbmuidcuts"))s=new glbMuIdCuts(Analysis,UncertType);
  else if(Analysis.Contains("loosemuidcuts"))s=new looseMuIdCuts(Analysis,UncertType);
  else if(Analysis.Contains("softmuidcuts"))s=new softMuIdCuts(Analysis,UncertType);
  else if(Analysis.Contains("trkmuidcuts"))s=new trkMuIdCuts(Analysis,UncertType);
  else if(Analysis.Contains("twoglboneloosemuidcuts"))s=new twoGlbOneLooseMuIdCuts(Analysis,UncertType);
  
  else if(Analysis.Contains("tmvasignal"))s=new TMVASignal(Analysis,UncertType);
  else if(Analysis.Contains("tmvatree_1"))s=new TMVATree_1(Analysis,UncertType);
  else if(Analysis.Contains("tmvatree_2"))s=new TMVATree_2(Analysis,UncertType);
  else if(Analysis.Contains("tmvatree_3"))s=new TMVATree_3(Analysis,UncertType);

  else if(Analysis.Contains("tmva_signal_glbmuon"))s=new TMVA_Signal_glbMuon(Analysis,UncertType);
  else if(Analysis.Contains("tmva_signal_loosemuon"))s=new TMVA_Signal_looseMuon(Analysis,UncertType);
  else if(Analysis.Contains("tmva_signal_softmuon"))s=new TMVA_Signal_softMuon(Analysis,UncertType);
  else if(Analysis.Contains("tmva_signal_trkmuon"))s=new TMVA_Signal_trkMuon(Analysis,UncertType);
  else if(Analysis.Contains("tmva_signal_twoglboneloosemuon"))s=new TMVA_Signal_twoGlbOneLooseMuon(Analysis,UncertType);

  else if(Analysis.Contains("tmva_tree_1_glbmuon"))s=new TMVA_Tree_1_glbMuon(Analysis,UncertType);
  else if(Analysis.Contains("tmva_tree_1_loosemuon"))s=new TMVA_Tree_1_looseMuon(Analysis,UncertType);
  else if(Analysis.Contains("tmva_tree_1_softmuon"))s=new TMVA_Tree_1_softMuon(Analysis,UncertType);
  else if(Analysis.Contains("tmva_tree_1_trkmuon"))s=new TMVA_Tree_1_trkMuon(Analysis,UncertType);
  else if(Analysis.Contains("tmva_tree_1_twoglboneloosemuon"))s=new TMVA_Tree_1_twoGlbOneLooseMuon(Analysis,UncertType);

  else if(Analysis.Contains("tmva_tree_2_glbmuon"))s=new TMVA_Tree_2_glbMuon(Analysis,UncertType);
  else if(Analysis.Contains("tmva_tree_2_loosemuon"))s=new TMVA_Tree_2_looseMuon(Analysis,UncertType);
  else if(Analysis.Contains("tmva_tree_2_softmuon"))s=new TMVA_Tree_2_softMuon(Analysis,UncertType);
  else if(Analysis.Contains("tmva_tree_2_trkmuon"))s=new TMVA_Tree_2_trkMuon(Analysis,UncertType);
  else if(Analysis.Contains("tmva_tree_2_twoglboneloosemuon"))s=new TMVA_Tree_2_twoGlbOneLooseMuon(Analysis,UncertType);

  else if(Analysis.Contains("tmva_tree_3_glbmuon"))s=new TMVA_Tree_3_glbMuon(Analysis,UncertType);
  else if(Analysis.Contains("tmva_tree_3_loosemuon"))s=new TMVA_Tree_3_looseMuon(Analysis,UncertType);
  else if(Analysis.Contains("tmva_tree_3_softmuon"))s=new TMVA_Tree_3_softMuon(Analysis,UncertType);
  else if(Analysis.Contains("tmva_tree_3_trkmuon"))s=new TMVA_Tree_3_trkMuon(Analysis,UncertType);
  else if(Analysis.Contains("tmva_tree_3_twoglboneloosemuon"))s=new TMVA_Tree_3_twoGlbOneLooseMuon(Analysis,UncertType);
#endif


#ifdef USE_wang
  else if(Analysis.Contains("mytest"))s=new MyTest(Analysis,UncertType);
#endif


#ifdef USE_menendez
  else if(Analysis.Contains("mytest"))s=new MyTest(Analysis,UncertType);
  else if(Analysis.Contains("dstophipi"))s=new DsToPhiPi(Analysis,UncertType);
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
