import os
from ROOT import *
import argparse

BASE_DIR = ''

# ---------- Argument parser -------------
parser = argparse.ArgumentParser()

#parser.add_argument("--verbose", help="print event by event output", action="store_true")
parser.add_argument("--classifier", help="name of the ouput classifier", default="MyTMVAClassification.cxx")

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

spectatorList = ['var_tauMass','var_tauMassRes', 'var_Eta_Tau']

# specify cuts=1 to implement cuts, cuts=0 to not
def makeClassifiers(classifierName, inputFile, outputFile, category, varset, cut):
        
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

   _varlist = varset

   for line in lines:
      if ('TMVAClassification_template' in line):
         print "found it"
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
            if item not in spectatorList: oFile.write("   dataloader->AddVariable(\""+item+"\",\"Variable "+item.replace('var_','')+"\",\"units\", \'F\' );\n")
         for item in spectatorList: oFile.write("   dataloader->AddSpectator(\""+item+"\",\"Variable "+item.replace('var_','')+"\",\"units\", \'F\' );\n")

      elif _cut_tag in line:
         if (cut==1):
            for var in _varlist:
               if (var in var_limits):
                  cuts = cuts+"("+var+" > "+str(var_limits[var][0])+" && "+var+" < "+str(var_limits[var][1])+") && "
                  cutb = cutb+"("+var+" > "+str(var_limits[var][0])+" && "+var+" < "+str(var_limits[var][1])+") && "
         if (category=='A'):
            cuts += "(category == 1) && "
            cutb += "(category == 1 || category==4) && "
         if (category=='B'):
            cuts += "(category == 2) && "
            cutb += "(category == 2 || category==5) && "
         if (category=='C'):
            cuts += "(category == 3) && "
            cutb += "(category == 3 || category==6) && "
         if (category=='Barrel'):
            cuts += "var_Eta_Tau <= 0.8 && "
            cutb += "var_Eta_Tau <= 0.8 && "
         if (category=='Endcap'):
            cuts += "var_Eta_Tau > 0.8 && "
            cutb += "var_Eta_Tau > 0.8 && "
         oFile.write("TCut mycutb = \""+cutb+" MC==0 && abs(var_MaxD0Significance)<100\";\n")
         oFile.write("TCut mycuts = \""+cuts+" MC==1 && abs(var_MaxD0Significance)<100\";\n")

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

var_limits = { 'var_cLP': [0,100],
                                'var_max_tKink': [0,100],
                                'var_mindca_iso': [0,0.05],
                                'var_trkRel': [0,10],
                                'var_vertexKFChi2': [0,30],
                                'var_svpvTauAngle': [0,0.2],
                                'var_flightLenSig': [0,20],
                                'var_segCompMuMin': [0,1],
                                'var_MaxD0Significance': [0,6],
                                'var_MinD0Significance': [0,6],
                                'var_MuMu_minKFChi2': [0,3],
                                'var_sumMuTrkKinkChi2': [0,99],
                                'var_MuMu_mindR': [0,0.25],
                                'var_maxdca': [0,0.15],
                                'var_MuTau_maxdR': [0,1],
                                'var_RelPt_Mu1Tau': [0,100],
                                'var_MinMIPLikelihood': [0,1],
                                'var_pmin': [0,25],
                                'var_max_cLP': [0,60]
                                }

varsets = {'2016vars':['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig', 'var_segCompMuMin', 'var_pmin', 'var_max_cLP', 'var_max_tKink', 'var_MinD0Significance', 'var_mindca_iso', 'var_trk_relPt'],
        '2017biased':['var_vertexKFChi2','var_svpvTauAngle','var_flightLenSig','var_segCompMuMin','var_MuMu_minKFChi2','var_MuTau_maxdR','var_sumMuTrkKinkChi2','var_MaxD0Significance', 'var_MinMIPLikelihood', 'var_MuMu_mindR', 'var_maxdca', 'var_mindca_iso', 'var_pmin'],
        '2017unbiased':['var_vertexKFChi2','var_svpvTauAngle','var_flightLenSig','var_segCompMuMin','var_MuMu_minKFChi2','var_sumMuTrkKinkChi2','var_MaxD0Significance', 'var_MinMIPLikelihood', 'var_maxdca', 'var_mindca_iso','var_pmin']}

varsets_noVeto = {'2016vars':['var_vertexKFChi2', 'var_svpvTauAngle', 'var_flightLenSig', 'var_segCompMuMin', 'var_pmin', 'var_max_cLP', 'var_max_tKink', 'var_MinD0Significance', 'var_mindca_iso', 'var_trk_relPt'],
        '2017biased':['var_vertexKFChi2','var_svpvTauAngle','var_flightLenSig','var_segCompMuMin','var_MuMu_minKFChi2','var_MuTau_maxdR','var_sumMuTrkKinkChi2','var_MaxD0Significance', 'var_MinMIPLikelihood', 'var_MuMu_mindR', 'var_maxdca', 'var_RelPt_Mu1Tau', 'var_mindca_iso'],
        '2017unbiased':['var_vertexKFChi2','var_svpvTauAngle','var_flightLenSig','var_segCompMuMin','var_sumMuTrkKinkChi2','var_MaxD0Significance', 'var_MinMIPLikelihood', 'var_maxdca', 'var_mindca_iso'],
        '2017unbiased_MuMuMass':['var_vertexKFChi2','var_svpvTauAngle','var_flightLenSig','var_segCompMuMin','var_sumMuTrkKinkChi2','var_MaxD0Significance', 'var_MinMIPLikelihood', 'var_maxdca', 'var_mindca_iso','var_phiMass','var_omegaMass']}



#categories = ['A','B','C','Barrel','Endcap']
categories = ['A','B','C']
datasets = {'Veto':'FillMVATree','NoVeto':'MVATree_NoVeto'}

if __name__== "__main__":
   os.system("rm TMVATraining_2glb1trk.log")
   for data in datasets:
        if (data=='Veto') : tmpSets = varsets
        else: tmpSets = varsets_noVeto
        for cat in categories:  
                for varset in tmpSets:
                   if '2017' in varset: continue
                   classifierName = args.classifier+"_TwoGlbMuOneTrk_"+data+"_"+varset+"_"+cat
                   #inputFile =  'Input.root'
                   #outputFile = 'Output.root'
                   if (data=="Veto"): inputFile = '../../macros/TMVATrees/TMVATrees_TwoGlobalMuonOneTracker_BDT_Run2017_MC_combined.root'
                   else: continue
                   #else: inputFile = '../../macros/TMVATrees/TMVATrees_MVATree_NoVeto_Run2017_MC_combined.root'
                   outputFile = classifierName+"_Output.root"
                   makeClassifiers(classifierName, inputFile, outputFile, cat, tmpSets[varset], 1)
                   #makeClassifiers(classifierName, inputFile, outputFile, cat, tmpSets[varset], 1)
                   os.system("rm -rf "+classifierName)
                   os.system("root -q -b -l "+classifierName+".cxx\(\\\"BDT\\\"\) &>> TMVATraining_2glb1trk.log")
                   os.system("mv datasets "+classifierName)

