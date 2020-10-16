#!/usr/bin/env python

import os
from ROOT import *
import argparse
import json

BASE_DIR = ''

def makeClassifiers(classifierName, inputFile, outputFile, variableSet):
   
   #-------------------------------
   # tags
   #-------------------------------
   _input_tag = "// InputFileTag"
   _var_tag = "// VariableTag"
   _output_tag = "// OutputFileTag"
   _cut_tag = "// CutTag"
   _training_tag = "// TrainingTag"
   _factory_tag = "// FactoryNameTag"
   _template_filename = "./TrackerMuonId_TMVATemplate.cxx"
   
   #-------------------------------
   # open files
   #-------------------------------
   _tmp_filename = "./TrackerMuonId_TMVATemplate.cxx.tmp"
   _dest_filename = classifierName+'.cxx'
   templateFile = open(_template_filename, "r")
   classifierFile = open(_tmp_filename,'w')
   
   jsonFile = open("TrackerMuonIdVariables.json")
   variable_sets = json.load(jsonFile)       
   #-------------------------------
   
   lines = templateFile.readlines()
   templateFile.close()

   for line in lines:
      if ('TrackerMuonId_TMVATemplate' in line):
         line = line.replace('TrackerMuonId_TMVATemplate',classifierName)
      
      if _input_tag in line:
         classifierFile.write("TString fname = \""+inputFile+"\";")
         classifierFile.write(line)

      elif _output_tag in line:
         classifierFile.write("TString outfileName = \""+outputFile+"\";")
         classifierFile.write(line)
      
      elif _factory_tag in line:
         classifierFile.write("TMVA::Factory *factory = new TMVA::Factory(\""+classifierName+"\",\"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification\");\n");
         classifierFile.write("TMVA::DataLoader *dataloader=new TMVA::DataLoader(\""+classifierName+"\");\n");

      elif _var_tag in line:
         for item in variable_sets[variableSet]["training_variables"]:
             classifierFile.write("dataloader->AddVariable(\""+item+"\",\"Variable \",\"units\", \'"+variable_sets[variableSet]['training_variables'][item]['type']+"\' );\n")
         for item in variable_sets[variableSet]["spectator_variables"]:
             classifierFile.write("dataloader->AddSpectator(\""+item+"\",\"Variable \",\"units\", \'"+variable_sets[variableSet]['spectator_variables'][item]['type']+"\' );\n")

      #elif _training_tag in line:
      #   classifierFile.write("dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,\n\"nTrain_Signal=0:nTrain_Background=0:nTest_signal=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V\");\n")

      else:
         classifierFile.write(line)
   #----------------------------------------
   # close files
   #----------------------------------------
   jsonFile.close()
   templateFile.close()
   classifierFile.close()
   os.system("vim -c \"execute \'normal!=G\' | :wq\" "+_tmp_filename)
   os.rename(_tmp_filename, _dest_filename)
   os.system("chmod +x "+_dest_filename)
   #----------------------------------------



def main():

    # ---------- Argument parser -------------
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--classifier", help="Name of the ouput classifier;[Default: %(default)s]", default="TMVAClassification")
    parser.add_argument("-i","--input-file", help="Name of the input file for TMVA classification; [Default: %(default)s]", default="../Ntuples/MuonPionTree_combined_2017.root")
    parser.add_argument("-o","--output-file", help="Name of the output file; [Default: %(default)s]", default="TMVAOutput/TrackerMuonId_output_example.root")
    parser.add_argument("-v","--variable-set", help="Set of variables used for training; [Default: %(default)s]", default="2017_standard_variables")
    args = parser.parse_args()
    # ----------------------------------------
    
    classifierName = args.classifier
    inputFile = args.input_file 
    outputFile = args.output_file 
    variableSet = args.variable_set
    makeClassifiers(classifierName, inputFile, outputFile, variableSet)
    os.system("root -q -b -l "+classifierName+".cxx\(\\\"BDT\\\"\) &>> TMVATraining.log")

if __name__== "__main__":
    main()
