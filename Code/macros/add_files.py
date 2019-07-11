import os, sys
import argparse
import ROOT as r

#----- Argument parser --------
parser = argparse.ArgumentParser()
parser.add_argument('--instanceName', help="Name of the class to be combined.")
parser.add_argument('--inputFile', help="input file.")
args = parser.parse_args()
#------------------------------

SRCDIR = '../../'

mclist = []
datalist = []

setlist = [f for f in os.listdir(SRCDIR) if 'Set_' in f]
for st in setlist:
   ifile = open(SRCDIR+st+"/Input.txt")
   lines = ifile.readlines()
   datasetline = [l for l in lines if 'InputNtuples:' in l]
   datasetline = datasetline[0].strip('\r\n')
   if ( ('DoubleMuonLowMass' in datasetline) and ('2017F' in datasetline) ):
      datatype = 'Data'
   elif ( ('DsTau3Mu' in datasetline) or ('B0Tau3Mu' in datasetline) or ('BpTau3Mu' in datasetline) ): datatype = 'MC'
   #elif ( 'DsTau3Mu' in datasetline ): datatype = 'MC'
   else: datatype = None
   if (datatype=='MC'): mclist.append(st)
   elif (datatype=='Data'): datalist.append(st)
   ifile.close()

if (args.instanceName==None):
   print "Please provide the name of the class."

cmd = 'hadd -f '+args.instanceName+'_MC_merged.root'

#merge mc files
for mcset in mclist:
    path=(SRCDIR+mcset+"/"+args.inputFile+".root")
    if not (os.path.isfile(path)):
        print "File "+path+" does not exist"
        continue
    cmd+=' '+path

os.system(cmd) 

cmd = 'hadd -f '+args.instanceName+'_DATA_merged.root'

#merge data files
for dset in datalist:
   path=(SRCDIR+dset+"/"+args.inputFile+".root")
   if not(os.path.isfile(path)):
      print "File "+path+" does not exist"
      continue
   cmd+=' '+path
os.system(cmd) 
