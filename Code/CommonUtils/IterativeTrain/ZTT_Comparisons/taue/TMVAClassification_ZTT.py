#!/usr/bin/env python


import sys    # exit
import time   # time accounting
import getopt # command line parser
from TrainConfigs_taue_comparisons import configuration,selection

DEFAULT_INFNAME  = "../../MVA_Mini_Tree_ZTT.root"


#DEFAULT_METHODS  = "Cuts,CutsD,CutsPCA,CutsGA,CutsSA,Likelihood,LikelihoodD,LikelihoodPCA,LikelihoodKDE,LikelihoodMIX,PDERS,PDERSD,PDERSPCA,PDEFoam,PDEFoamBoost,KNN,LD,Fisher,FisherG,BoostedFisher,HMatrix,FDA_GA,FDA_SA,FDA_MC,FDA_MT,FDA_GAMT,FDA_MCMT,MLP,MLPBFGS,MLPBNN,CFMlpANN,TMlpANN,SVM,BDT,BDTD,BDTG,BDTB,RuleFit"


#DEFAULT_METHODS  = "BDT,BDTG,Likelihood"
DEFAULT_METHODS  = "BDT"



def usage():
    print " "
    print "Usage: python %s [options]" % sys.argv[0]
    print "  -m | --methods    : gives methods to be run (default: all methods)" % DEFAULT_METHODS
    print "  -i | --inputfile  : name of input ROOT file (default: '%s')" % DEFAULT_INFNAME
    print "  -v | --verbose"
    print "  -? | --usage      : print this help message"
    print " "





def doTrain(configs,training_cuts,mlist,infname):

    count = 0
    for train in configs:
        for category_wagon in train.keys():

            print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  ",category_wagon
            prefix = ""

            outfname = "category_"+str(count)+"_"+category_wagon+prefix+".root"
            outputFile = TFile(outfname , 'RECREATE' )


            factory = TMVA.Factory( "TMVAClassification", outputFile,
                                    "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" )

            dataloader = TMVA.DataLoader("output"+"_"+str(count)+"_"+category_wagon + prefix)





#            outfname = "Training_N_"+str(count)+"_"+category_wagon+prefix+".root"
#            outputFile = TFile(outfname , 'RECREATE' )
#            factory = TMVA.Factory( "TMVAClassification", outputFile,
#                                    "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" )
#            dataloader = TMVA.DataLoader("Trainin_N_"+str(count)+"_"+category_wagon + prefix)
 

            factory.SetVerbose( True  )




            # Setup dataloader and define cuts
            cuts=''
            cutb='dataMCtype      == 1      && '

            for v in train.get(category_wagon):
                dataloader.AddVariable(v,v,"F")

                if v in selection:

                  cuts = cuts+"("+v+" > "+str(selection[v][0])+ \
                          " && "+v+" < "+str(selection[v][1])+")  && "
                  cutb = cutb+"("+v+" > "+str(selection[v][0])+ \
                          " && "+v+" < "+str(selection[v][1])+")  && "




            input = TFile.Open(  infname  )
            signalWeight  = 1.0;
            backgroundWeight = 1.0;
            categoryCut = ''
            if category_wagon=="ZTT_mu3mu":
                categoryCut=''
                cuts  = cuts + 'dataMCtype      == 210232   '
                signal_ztt_mu3mu= input.Get("z2tautau_ztau3mutaumu")
                background      = input.Get("DoubleMuonLowMass_ztau3mutaumu")
                dataloader.AddSignalTree    ( signal_ztt_mu3mu,      signalWeight          );
                dataloader.AddBackgroundTree( background,            backgroundWeight      );

            if category_wagon=="ZTT_tau3mu":
                categoryCut=''
                cuts  = cuts + 'dataMCtype      == 210233   '
                signal_ztt_tau3mu= input.Get("z2tautau_ztau3mutauh")
                background       = input.Get("DoubleMuonLowMass_ztau3mutauh")
                dataloader.AddSignalTree    ( signal_ztt_tau3mu,     signalWeight          );
                dataloader.AddBackgroundTree( background,            backgroundWeight      );
                
                
            if category_wagon=="ZTT_tau_NoCV_3mu":
                categoryCut=''
                cuts  = cuts + 'dataMCtype      == 210233  && ifCommonCV == 0   '
                signal_ztt_tau3mu= input.Get("z2tautau_ztau3mutauh_A")
                background       = input.Get("DoubleMuonLowMass_ztau3mutauh_A")
                dataloader.AddSignalTree    ( signal_ztt_tau3mu,     signalWeight          );
                dataloader.AddBackgroundTree( background,            backgroundWeight      );
                
                
            if category_wagon=="ZTT_tau_CV_3mu":
                categoryCut=''
                cuts  = cuts + 'dataMCtype      == 210233  && ifCommonCV == 1   '
                signal_ztt_tau3mu= input.Get("z2tautau_ztau3mutauh_B")
                background       = input.Get("DoubleMuonLowMass_ztau3mutauh_B")
                dataloader.AddSignalTree    ( signal_ztt_tau3mu,     signalWeight          );
                dataloader.AddBackgroundTree( background,            backgroundWeight      );


            if category_wagon=="ZTT_e3mu":
                categoryCut=''
                cuts  = cuts + 'dataMCtype      == 210231   '
                signal_ztt_e3mu= input.Get("z2tautau_ztau3mutaue")
                background     = input.Get("DoubleMuonLowMass_ztau3mutaue")
                dataloader.AddSignalTree    ( signal_ztt_e3mu,      signalWeight          );
                dataloader.AddBackgroundTree( background,           backgroundWeight      );


            mycutSig = TCut(  cuts + categoryCut) 
            mycutBkg = TCut(  cutb + categoryCut + " (m3m < 1.75 || m3m > 1.81) ") 


            print "-------------------------------------------------------------------------------------------- Signal Cut", mycutSig
            print "-------------------------------------------------------------------------------------------- Bkg    Cut", mycutBkg


            dataloader.PrepareTrainingAndTestTree( mycutSig, mycutBkg,
                                                   "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" )


            if "BDTG" in mlist:
                factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDTG",
                                    "!H:!V:NTrees=1000:MinNodeSize=1.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:MaxDepth=3" )       
                                    
                                    
                                                     

            if category_wagon=="ZTT_mu3mu":
                    if "BDT" in mlist:
                        factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDT",
                                            "!H:!V:NTrees=800:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.1:SeparationType=GiniIndex:nCuts=100" )
                                            
            if category_wagon=="ZTT_tau3mu":
                    if "BDT" in mlist:
                        factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDT",
                                            "!H:!V:NTrees=800:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.1:SeparationType=GiniIndex:nCuts=100" )
                                            
            if category_wagon=="ZTT_tau_NoCV_3mu":
                    if "BDT" in mlist:
                        factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDT",
                                            "!H:!V:NTrees=800:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.1:SeparationType=GiniIndex:nCuts=100" )
                                            
            if category_wagon=="ZTT_tau_CV_3mu":
                    if "BDT" in mlist:
                        factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDT",
                                            "!H:!V:NTrees=800:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.1:SeparationType=GiniIndex:nCuts=100" )
                                            
            if category_wagon=="ZTT_e3mu":
                    if "BDT" in mlist:
                        factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDT",
                                            "!H:!V:NTrees=800:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.1:SeparationType=GiniIndex:nCuts=100" )




            if "Likelihood" in mlist:
                factory.BookMethod(dataloader, TMVA.Types.kLikelihood, "Likelihood",
                                    "H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" )

            if "MLP" in mlist:
                factory.BookMethod(dataloader, TMVA.Types.kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" )


            factory.TrainAllMethods()
            factory.TestAllMethods()
            factory.EvaluateAllMethods()    

            # Save the output.
            outputFile.Close()
    
            print "=== wrote root file %s\n" % outfname
            print "=== TMVAClassification is done!\n"
    
            # open the GUI for the result macros    
            # TMVA.TMVAGui(outfname)
    
            # keep the ROOT thread running
            #gApplication.Run() 

        count+=1  # counter over trainings defined in TrainConfig




if __name__ == "__main__":

    from ROOT import gSystem, gROOT, gApplication, TFile, TTree, TCut
    if gROOT.GetVersionCode() >= 332288 and gROOT.GetVersionCode() < 332544:
        print "*** You are running ROOT version 5.18, which has problems in PyROOT such that TMVA"
        print "*** does not run properly (function calls with enums in the argument are ignored)."
        print "*** Solution: either use CINT or a C++ compiled version (see TMVA/macros or TMVA/examples),"
        print "*** or use another ROOT version (e.g., ROOT 5.19)."
        sys.exit(1)



    try:
        # retrive command line options
        shortopts  = "m:i:t:o:vh?"
        longopts   = ["methods=", "inputfile=", "verbose", "help", "usage"]
        opts, args = getopt.getopt( sys.argv[1:], shortopts, longopts )

    except getopt.GetoptError:
        # print help information and exit:
        print "ERROR: unknown options in argument %s" % sys.argv[1:]
        usage()
        sys.exit(1)



    infname     = DEFAULT_INFNAME
    methods     = DEFAULT_METHODS
    verbose     = False

    for o, a in opts:
        if o in ("--usage"):
            usage()
            sys.exit(0)
        elif o in ("-m", "--methods"):
            methods = a
        elif o in ("-i", "--inputfile"):
            infname = a
        elif o in ("-v", "--verbose"):
            verbose = True
     


    # Print methods
    mlist = methods.replace(' ',',').split(',')
    print "=== TMVAClassification: use method(s)..."
    for m in mlist:
        if m.strip() != '':
            print "=== - <%s>" % m.strip()

    # Run TMVA with given configs
    from ROOT import TMVA
    doTrain(configuration, selection, mlist,infname)

    

