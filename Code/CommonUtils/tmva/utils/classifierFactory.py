#!/usr/bin/env python
import os
import ROOT
import argparse
from varlist import var_limits, varsets

BASE_DIR = ""

# ---------- Argument parser -------------
parser = argparse.ArgumentParser()

#parser.add_argument("--verbose", help="print event by event output", action="store_true")
parser.add_argument("-c", "--classifier", help="name of the ouput classifier; Default: %(default)s", default="TMVAClassification")
parser.add_argument("-i", "--input-file", help="input to the classifier", default="FillMVATreeInput_combined.root")
parser.add_argument("-g", "--glb", help="1: three global category; 0: two global and tracker category", default=1)
args = parser.parse_args()
# ----------------------------------------

input_file = "../../macros/"+args.input_file
# tags
_input_tag = "// Input file"
_var_tag = "// Add variables"
_output_tag = "// Output file"
_cut_tag = "// Add cuts"
_training_tag = "// Prepare training"
_factory_tag = "// Factory name"
_template_filename = "../TMVAClassification_template.cxx"

categoryCut_bkg = {"A": "category==1 && (var_tauMass > 1.825 || var_tauMass < 1.729)", \
        "B": "category==2 && (var_tauMass > 1.853 || var_tauMass < 1.701)", \
        "C":"category==3 && (var_tauMass > 1.877 || var_tauMass < 1.677)"}
categoryCut_signal = {"A": "category==1","B": "category==2","C":"category==3"}

# specify cuts=1 to implement cuts, cuts=0 to not
def makeClassifiers(classifierName, input_file, outputFile, varset, cat, glb, cut):
   
   print "Producing classifier "+classifierName+" ..."
   #-------------------------------
   # open files
   #-------------------------------

   _tmp_filename = "./TMVAClassification.cxx.tmp"
   _dest_filename = classifierName+".cxx"
   templateFile = open(_template_filename, "r")
   oFile = open(_tmp_filename,"w")
   #-------------------------------

   cuts = ""
   cutb = ""
   lines = templateFile.readlines()
   templateFile.close()

   rootfile = ROOT.TFile(input_file, "READ")

   _varlist = varsets[varset]["training"]
   spectatorList = varsets[varset]["spectator"]

   for line in lines:
      if ("TMVAClassification_template" in line):
         line = line.replace("TMVAClassification_template",classifierName)
      
      if _input_tag in line:
         oFile.write("TString fname = \""+input_file+"\";")
         oFile.write(line)

      elif _output_tag in line:
         oFile.write("TString outfileName = \""+outputFile+"\";")
         oFile.write(line)
      
      elif _factory_tag in line:
         oFile.write("TString factoryName = \""+classifierName+"\";")

      elif _var_tag in line:
         for item in _varlist:
            Type = "F"
            if ("category" in item or "MC" in item): continue
            if ("var_trackerMuonId" in item and glb==1): continue
            if ("var_minMatchedStations" in item): Type = "I"
            low_expression = item+" > "+str(var_limits[item][0])+" ? "+item+":"+str(var_limits[item][0])
            expression = "("+low_expression+") < "+str(var_limits[item][1])+\
                    " ? ("+low_expression+"):"+str(var_limits[item][1])
            if item not in spectatorList: oFile.write("   dataloader->AddVariable(\""+ \
                    expression+"\",\"Variable "+item.replace("var_","")+"\",\"units\", \'"+Type+"\' );\n")
         for item in spectatorList: oFile.write("   dataloader->AddSpectator(\""+item+ \
                 "\",\"Variable "+item.replace("var_","")+"\",\"units\", \'F\' );\n")

      elif _cut_tag in line:
         if (cut==1):
            for var in _varlist:
               if (var in var_limits):
                  cuts = cuts+"("+var+" > "+str(var_limits[var][0])+ \
                          " && "+var+" < "+str(var_limits[var][1])+") && "
                  cutb = cutb+"("+var+" > "+str(var_limits[var][0])+ \
                          " && "+var+" < "+str(var_limits[var][1])+") && "
         cuts+=categoryCut_signal[cat]+" && threeGlobal=="+str(glb)
         cutb+=categoryCut_bkg[cat]+" && threeGlobal=="+str(glb)
         oFile.write("TCut mycutb = \""+cutb+"\";\n")
         oFile.write("TCut mycuts = \""+cuts+"\";\n")

      elif _training_tag in line:
         
         TreeB = rootfile.Get("TreeB")
         TreeS_Ds = rootfile.Get("TreeS_Ds")
         TreeS_Bu = rootfile.Get("TreeS_Bu")
         TreeS_Bd = rootfile.Get("TreeS_Bd")
         
         nevents_bkg = TreeB.Draw("category", cutb)
         nevents_signal = TreeS_Ds.Draw("category", cuts) + \
                 TreeS_Bu.Draw("category", cuts) + \
                 TreeS_Bd.Draw("category", cuts)
      
         split_factor = 0.7 # ratio of data used for training and testing; to be added as an argument...
         
         nTrain_bkg = int(split_factor*nevents_bkg)
         nTrain_signal = int(split_factor*nevents_signal)
         nTest_bkg = nevents_bkg - nTrain_bkg
         nTest_signal = nevents_signal - nTrain_signal
         oFile.write("dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,\n\"nTrain_Signal="+str(nTrain_signal)+ \
                 ":nTrain_Background="+str(nTrain_bkg)+":nTest_signal="+str(nTest_signal)+":nTest_Background="+str(nTest_bkg)+ \
                 ":SplitMode=Random:NormMode=NumEvents:!V\");\n")

      else:
         oFile.write(line)
        
   #----------------------------------------
   # close files
   #----------------------------------------
   templateFile.close()
   oFile.close()
   os.rename(_tmp_filename, _dest_filename)
   os.system("chmod +x "+_dest_filename)
   #----------------------------------------

#categories = ["A","B","C","Barrel","Endcap"]
categories = ["A","B","C"]

if __name__== "__main__":
   glb = args.glb
   input_file = args.input_file
   for cat in categories:  
      for varset in varsets:
         if (glb==1): classifierName = args.classifier+"_"+varset+"_"+cat+"_threeGlobal"
         else: classifierName = args.classifier+"_"+varset+"_"+cat+"_twoGlobalTracker"
         print classifierName
         outputFile = classifierName+"_Output.root"
         makeClassifiers(classifierName, input_file, outputFile, varset, cat, glb, 0)
         #makeClassifiers(classifierName, input_file, outputFile, varset, cat, glb, 1)    
         #os.system("cp ../classifiers/"+classifierName+".cxx .")
         os.system("root -q -b -l "+classifierName+".cxx\(\\\"BDT\\\"\) &> ../logs/"+classifierName+".log")
         if os.path.isdir("../weights/"+classifierName):
             os.system("rm -rf ../weights/"+classifierName)
         os.system("mv datasets ../weights/"+classifierName)
         os.system("mv "+classifierName+"_Output.root ../tmva_output")
         os.system("mv "+classifierName+".cxx ../classifiers")
