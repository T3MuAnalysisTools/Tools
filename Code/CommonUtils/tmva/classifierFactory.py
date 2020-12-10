#!/usr/bin/env python


import os
from ROOT import *
import argparse

BASE_DIR = ''

# ---------- Argument parser -------------
parser = argparse.ArgumentParser()

#parser.add_argument("--verbose", help="print event by event output", action="store_true")
parser.add_argument("-c","--classifier", help="name of the ouput classifier;[Default: %(default)s]", default="MyTMVAClassification")
parser.add_argument("-f","--inputFile", help="name of the of the input root with MVA trees; [Default: %(default)s]",  type = str, action="store", default = 'FillMVATreeInput_combined_August.root')
parser.add_argument("-m","--methods", help="name of methods (coma separated); [Default: %(default)s]",  type = str, action="store", default = "BDT")
args = parser.parse_args()
# ----------------------------------------

# tags
_input_tag = "// Input file"
_var_tag = "// Add variables"
_output_tag = "// Output file"
_cut_tag = "// Add cuts"
_training_tag = "// Prepare training"
_factory_tag = "// Factory name"
_template_filename = "./TMVAClassification_template.cxx"

spectatorList = ['var_tauMass']

# specify cuts=1 to implement cuts, cuts=0 to not
def makeClassifiers(classifierName, inputFile, outputFile,cat, varset, cut):
        
   #-------------------------------
   # open files
   #-------------------------------

   _tmp_filename = "./TMVAClassification.cxx.tmp"
   _dest_filename = classifierName+'.cxx'
   templateFile = open(_template_filename, "r")
   oFile = open(_tmp_filename,'w')
   #-------------------------------

   cutb = ""
   cuts = ""
   lines = templateFile.readlines()
   templateFile.close()

   _varlist = varsets[varset]

   for line in lines:
      if ('TMVAClassification_template' in line):
         line = line.replace('TMVAClassification_template',classifierName)
      
      if _input_tag in line:
         oFile.write("TString fname = \""+inputFile+"\";")
         oFile.write(line)

      elif _output_tag in line:
         oFile.write("TString outfileName = \""+outputFile+"\";")
         oFile.write(line)
      
      elif _factory_tag in line:
         oFile.write("TString factoryName = \""+classifierName+"\";")

      elif _var_tag in line:
         for item in _varlist:
            if ('category' in item or 'MC' in item): continue
            if item not in spectatorList: oFile.write("dataloader->AddVariable(\""+item+"\",\"Variable "+item.replace('var_','')+"\",\"units\", \'F\' );\n")
         for item in spectatorList: oFile.write("dataloader->AddSpectator(\""+item+"\",\"Variable "+item.replace('var_','')+"\",\"units\", \'F\' );\n")

      elif _cut_tag in line:
         if (cut==1):
            for var in _varlist:
               if (var in var_limits):
                  cuts = cuts+"("+var+" > "+str(var_limits[var][0])+" && "+var+" < "+str(var_limits[var][1])+") && "
                  cutb = cutb+"("+var+" > "+str(var_limits[var][0])+" && "+var+" < "+str(var_limits[var][1])+") && "
         if cat=="A": category_cut = " && category == 1  "
         if cat=="B": category_cut = " && category == 2  "
         if cat=="C": category_cut = " && category == 3  "
         if cat=="BC": category_cut = " && (category == 2 || category ==3)  "
         oFile.write("TCut mycutb = \""+cutb+" MC==0 " + category_cut + "  && var_flightLenSig < 100.  && var_maxMuonsDca < 0.3   &&   var_MaxVertexPairQuality < 30. &&  var_svpvTauAngle < 0.2   && var_MaxtrkKink  < 100 && var_MindcaTrackSV > 0 && var_MindcaTrackSV < 0.5 &&  var_MuonglbkinkSum < 100   && (var_tauMass < 1.75 || var_tauMass > 1.81)  \";\n")
         oFile.write("TCut mycuts = \""+cuts+" MC==1 " + category_cut + "  && var_flightLenSig < 100.  && var_maxMuonsDca < 0.3   && var_MaxVertexPairQuality < 30. &&  var_svpvTauAngle < 0.2  && var_MaxtrkKink  < 100 && var_MindcaTrackSV > 0 && var_MindcaTrackSV < 0.5 &&  var_MuonglbkinkSum < 100        \";\n")

      elif _training_tag in line:
         oFile.write("dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,\n\"nTrain_Signal=0:nTrain_Background=0:nTest_signal=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V\");\n")

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

var_limits = {  'var_min_p':[0,10],
                'var_cLP': [0,100],
                'var_max_tKink': [0,100],
                'var_mindca_iso': [0,3],
                'var_trkRel': [0,10],
                'var_vertexKFChi2': [0,30],
                'var_svpvTauAngle': [0,0.2],
                'var_flightLenSig': [0,20],
                'var_segCompMuMin': [0.2,1],
                'var_MaxD0Significance': [0,6],
                'var_MinD0Significance': [0,6],
                'var_MuMu_minKFChi2': [0,3],
                'var_sumMuTrkKinkChi2': [0,99],
                'var_MuMu_mindR': [0,0.15],
                'var_maxdca': [0,0.15],
                'var_MuTau_maxdR': [0,1],
                'var_RelPt_Mu1Tau': [0,100],
                'var_MinMIPLikelihood': [0,1],
                'var_pmin': [0,10],
                'var_max_cLP': [0,60]
}




varsets = {'A':['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig',
                ,'var_maxMuonsDca',
                'var_Muon1DetID','var_Muon2DetID','var_Muon3DetID' ,'var_MaxVertexPairQuality'
                ],


           'B':['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig',
                ,'var_maxMuonsDca',
                'var_Muon1DetID','var_Muon2DetID','var_Muon3DetID' ,'var_MaxVertexPairQuality'
                ],


           'C':['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig',
                'var_maxMuonsDca',
                'var_Muon1DetID','var_Muon2DetID','var_Muon3DetID' ,'var_MaxVertexPairQuality'
                ],



           'BC':['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig',
                 'var_MaxD0SigSV','var_MindcaTrackSV','var_maxMuonsDca',
                 'var_segCompMuMin','var_MaxMuon_chi2LocalPosition','var_MaxVertexPairQuality'
                 ,'var_MaxMuon_chi2LocalMomentum']}





print varsets

categories = ['A','B','C','BC']
datasets = {'Vars':'FillMVATree'}

if __name__== "__main__":
      for varset in varsets:
            print varset
            cat = varset
            classifierName = args.classifier+"_"+cat
            inputFile = args.inputFile
            outputFile = classifierName+"_Output.root"
            makeClassifiers(classifierName, inputFile, outputFile, cat,  varset, 0)
            os.system("rm TMVATraining.log")
            print classifierName
            print args.methods
#            os.system("root -q -b -l "+classifierName+".cxx\(\\\ "+args.methods+"\\\"\) &>> TMVATraining.log")
#            os.system("root -q -b -l "+classifierName+".cxx\(\\\""+args.methods+"\\\"\) &>> TMVATraining.log")

