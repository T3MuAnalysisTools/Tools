import os
import argparse

#----- Argument parser --------
parser = argparse.ArgumentParser()

parser.add_argument("--classifier", help="name of the ouput classifier", default="MyTMVAClassification.cxx")
parser.add_argument("--methods", help="classification methods", nargs='+')

args = parser.parse_args()
#------------------------------

#CLASSIFIER = args.classifier.replace('.C','')
DIR = 'dataset/weights/'

TMVATreeList = [ 'Run2017', 'Run2017B', 'Run']

#------------------------------
#preparing classes
for i in range(1,10):
"python makeTMVAClassifier.py --inputFile ../../macros/TMVATrees_Run2017_MC.root --outputFile TMVATrees_Run2017_MC_Output.root --category 1 --classifier Run2017Classification_+str(i)
#------------------------------
methodlist = args.methods
print methodlist
for i in range(1,10):
   CLASSIFIER = args.classifier.replace('.cxx','')
   CLASSIFIER+='_'+str(i)
   cmd = 'root -q -b -l ./'+CLASSIFIER+'.cxx\(\\"';
   for num, method in enumerate(methodlist):
      cmd+=method
      if not (num==len(methodlist)-1): cmd+=','
   cmd+='\\"\)'
   print cmd
   os.system(cmd)
