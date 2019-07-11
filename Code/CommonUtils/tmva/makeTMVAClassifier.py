import os
from ROOT import *
import argparse

#----- Argument parser --------
parser = argparse.ArgumentParser()

#parser.add_argument("--verbose", help="print event by event output", action="store_true")
parser.add_argument("--inputFile", help="input root file")
parser.add_argument("--outputFile", help="output root file")
parser.add_argument("--classifier", help="name of the ouput classifier", default="MyTMVAClassification.cxx")

args = parser.parse_args()
#------------------------------

_var_tag = "// Add all the variables"
_input_tag = "// input file"
_output_tag = "// output file"
_template_filename = "./MyTMVAClassification_template.cxx"
_tmp_filename = "./MyTMVAClassification.cxx.tmp"
_dest_filename = args.classifier

#-------------------------------
# open files

rootFile = TFile(args.inputFile,"READ")
templateFile = open(_template_filename,"r")
oFile = open(_tmp_filename, "w")

#------------------------------

_tree = rootFile.Get("TreeB")
_varlist = []

for item in _tree.GetListOfLeaves():
   _varlist.append(item.GetName())

lines = templateFile.readlines()

for line in lines:
   if ('MyTMVAClassification_template' in line):
      line = line.replace('MyTMVAClassification_template',args.inputFile)
      line = line.replace('_tmvatree','')
      line = line.replace('.root','')
      oFile.write(line)

   elif _input_tag in line:
      oFile.write("TString fname = \""+args.inputFile+"\";")
      oFile.write(line)
   
   elif _output_tag in line:
      oFile.write("TString outfileName = \""+args.outputFile+"\";")
      oFile.write(line)
   
   elif _var_tag in line:
      for item in _varlist:
         oFile.write("   dataloader->AddVariable(\""+item+"\",\"Variable "+item.replace('var_','')+"\",\"units\", \'F\' );\n") 
   else: 
      oFile.write(line) 

#------------------------------
# close files
rootFile.Close()
templateFile.close()
oFile.close()
os.rename(_tmp_filename, _dest_filename)
os.system("chmod +x "+_dest_filename)
#------------------------------
