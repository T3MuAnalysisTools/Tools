import os
from ROOT import *
import argparse

#----- Argument parser --------
parser = argparse.ArgumentParser()

#parser.add_argument("--verbose", help="print event by event output", action="store_true")
parser.add_argument("--inputFile", help="input root file")
parser.add_argument("--outputFile", help="output root file")
parser.add_argument("--classifier", help="name of the ouput classifier", default="MyTMVAClassification.cxx")
parser.add_argument("--category", help="category number [1-9]", default=0)

args = parser.parse_args()
#------------------------------
# tags

_var_tag = "// Add all the variables"
_input_tag = "// input file"
_output_tag = "// output file"
_cut_tag = "// Add cuts"
_training_tag = "// Prepare Training"
_factory_tag = "// Factory name"
_template_filename = "./MyTMVAClassification_template.cxx"
_tmp_filename = "./MyTMVAClassification.cxx.tmp"
_dest_filename = args.classifier+'.cxx'

#-------------------------------
# open files

rootFile = TFile(args.inputFile,"READ")
templateFile = open(_template_filename,"r")
oFile = open(_tmp_filename, "w")

#------------------------------

_treeB = rootFile.Get("TreeB")
_treeS = rootFile.Get("TreeS")
_varlist = []

Nbkg = _treeB.GetEntriesFast()
Nsig = _treeS.GetEntriesFast()

cut_bcategory = TCut("")
cut_scategory = TCut("")

if (args.category>0 and args.category<7): cut_bcategory = TCut("MC==0 && category=="+str(args.category))
if (args.category>0 and args.category>6): cut_bcategory = TCut("MC==0 && ( category=="+str(args.category)+"|| category=="+str((int(args.category)-1)%3+1)+")")
if (args.category>0): cut_scategory = TCut("MC==1 && category=="+str((int(args.category)-1)%3+1))

_treeB.Draw(">>elistb",cut_bcategory)
_treeS.Draw(">>elists",cut_scategory)

elist = gDirectory.Get("elistb");
elist = gDirectory.Get("elists");

NTrainBkg = (int)(elistb.GetN()/2); 
NTrainSig = (int)(elists.GetN()/2);

print "Number of events for training/testing in background sample = "+str(NTrainBkg)
print "Number of events for training/testing in signal sample = "+str(NTrainSig)

for item in _treeB.GetListOfLeaves():
   _varlist.append(item.GetName())

lines = templateFile.readlines()

for line in lines:
   if ('MyTMVAClassification_template' in line):
      line = line.replace('MyTMVAClassification_template',args.classifier)
      line = line.replace('_tmvatree','')
      line = line.replace('.root','')
      oFile.write(line)

   elif _factory_tag in line:
      oFile.write("TString factoryName = \""+args.classifier+"\";")
      
   elif _input_tag in line:
      oFile.write("TString fname = \""+args.inputFile+"\";")
      oFile.write(line)
   
   elif _output_tag in line:
      oFile.write("TString outfileName = \""+args.outputFile+"\";")
      oFile.write(line)
   
   elif _var_tag in line:
      for item in _varlist:
         if ( ('category' in item) or ('MC' in item) or ('Eta_au' in item)): continue
         if 'tauMass' not in item: oFile.write("   dataloader->AddVariable(\""+item+"\",\"Variable "+item.replace('var_','')+"\",\"units\", \'F\' );\n") 
         if 'tauMass' in item: oFile.write("   dataloader->AddSpectator(\""+item+"\",\"Variable "+item.replace('var_','')+"\",\"units\", \'F\' );\n") 
   
   elif _cut_tag in line:
      if (args.category>0 and args.category<7): oFile.write("TCut mycutb = \"category=="+str(args.category)+" && MC==0\";\n")
      else: oFile.write("TCut mycutb = \"( (category=="+str((int(args.category)-1)%3+4)+" || category=="+str((int(args.category)-1)%3+1)+") && MC==0)\";\n")
      #oFile.write("TCut mycuts = \"category==0\";\n")
      oFile.write("TCut mycuts = \"category=="+str((int(args.category)-1)%3+1)+" && MC==1\";\n")
   
   
   elif _training_tag in line:
      oFile.write("dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,\n\"nTrain_Signal="+str(NTrainSig)+":nTrain_Background="+str(NTrainBkg)+":SplitMode=Random:NormMode=NumEvents:!V\");\n")

   else: 
      oFile.write(line) 

#----------------------------------------
# close files
#----------------------------------------
rootFile.Close()
templateFile.close()
oFile.close()
os.rename(_tmp_filename, _dest_filename)
os.system("chmod +x "+_dest_filename)
#----------------------------------------
