#!/usr/bin/env python


import sys    # exit
import time   # time accounting
import getopt # command line parser
from TrainConfigs_D import configuration,selection

DEFAULT_INFNAME  = "FillMVATreeInput_2018combined.root"


#DEFAULT_METHODS  = "Cuts,CutsD,CutsPCA,CutsGA,CutsSA,Likelihood,LikelihoodD,LikelihoodPCA,LikelihoodKDE,LikelihoodMIX,PDERS,PDERSD,PDERSPCA,PDEFoam,PDEFoamBoost,KNN,LD,Fisher,FisherG,BoostedFisher,HMatrix,FDA_GA,FDA_SA,FDA_MC,FDA_MT,FDA_GAMT,FDA_MCMT,MLP,MLPBFGS,MLPBNN,CFMlpANN,TMlpANN,SVM,BDT,BDTD,BDTG,BDTB,RuleFit"


#DEFAULT_METHODS  = "BDT,BDTG,Likelihood"
DEFAULT_METHODS  = "BDT,BDTG"



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


            prefix = "_DS"
            outfname = str(count)+"_"+category_wagon+prefix+".root"
            outputFile = TFile(outfname , 'RECREATE' )


            factory = TMVA.Factory( "TMVAClassification", outputFile,
                                    "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" )

            dataloader = TMVA.DataLoader("output"+"_"+str(count)+"_"+category_wagon + prefix)
            factory.SetVerbose( True  )




            # Setup dataloader and define cuts
            cuts=''
            cutb=''

            for v in train.get(category_wagon):
                dataloader.AddVariable(v,v,"F")

                if v in selection:

                  cuts = cuts+"("+v+" > "+str(selection[v][0])+ \
                          " && "+v+" < "+str(selection[v][1])+") && "
                  cutb = cutb+"("+v+" > "+str(selection[v][0])+ \
                          " && "+v+" < "+str(selection[v][1])+") && "




            input = TFile.Open(  infname  )
            signal_ds= input.Get("TreeS_Ds")
            signal_bd= input.Get("TreeS_Bd")
            signal_bu= input.Get("TreeS_Bu")
            background= input.Get("TreeB")



            signalWeight_ds  = 1.0;
            backgroundWeight = 1.0;
  
 
            dataloader.AddSignalTree    ( signal_ds,     signalWeight_ds       );
            dataloader.AddBackgroundTree( background,    backgroundWeight      );



            categoryCut = ''
            if category_wagon=="A":categoryCut='category==1&& (var_mass12_dRsorting<0.994 || var_mass12_dRsorting> 1.044) && (var_mass13_drSorting<0.994 || var_mass13_drSorting> 1.044)'
            if category_wagon=="B":categoryCut='category==2&& (var_mass12_dRsorting<0.985 || var_mass12_dRsorting> 1.053) && (var_mass13_drSorting<0.985 || var_mass13_drSorting> 1.053) '
            if category_wagon=="C":categoryCut='category==3&& (var_mass12_dRsorting<0.974 || var_mass12_dRsorting> 1.064) && (var_mass13_drSorting<0.974 || var_mass13_drSorting> 1.064)'

            mycutSig = TCut( "MC == 1 &&" + cuts + categoryCut) 
            mycutBkg = TCut( "MC == 0 &&" + cutb + categoryCut + "&& (var_tauMass < 1.75 || var_tauMass > 1.81) ") 


            dataloader.PrepareTrainingAndTestTree( mycutSig, mycutBkg,
                                                   "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" )


            if "BDTG" in mlist:
                factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDTG",
                                    "!H:!V:NTrees=1000:MinNodeSize=1.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:MaxDepth=3" )                        

            if "BDT" in mlist:
                factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDT",
                                    "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" )


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

    

