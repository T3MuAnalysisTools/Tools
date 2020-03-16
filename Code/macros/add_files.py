#!/usr/bin/env python

import os, sys
import argparse
import ROOT as r

#----- Argument parser --------
parser = argparse.ArgumentParser()
parser.add_argument('--instanceName', help="Name of the class to be combined.",default='FillMVATree')
parser.add_argument('--inputFile', help="input file.", default='FillMVATreeInput.root')
args = parser.parse_args()
#------------------------------

SRCDIR = '../../'

setdict = {
'DoubleMuonLowMass':[],
'DsToTau':[],
'BuToTau':[],
'BdToTau':[]
}

setlist = [f for f in os.listdir(SRCDIR) if 'Set_' in f]
for st in setlist:
   ifile = open(SRCDIR+st+"/Input.txt")
   lines = ifile.readlines()
   ifile.close()
   datasetline = [l for l in lines if 'InputNtuples:' in l]
   datasetline = datasetline[0].strip('\r\n')
   for dset in setdict:
      if (dset in datasetline): setdict[dset].append(st)

if (args.instanceName==None):
   print "Please provide the name of the class."

def add_sets(datatype, setlist):
   cmd = 'hadd -f '+args.instanceName+'_'+datatype+'_merged.root'
   for mcset in setlist:
       path=(SRCDIR+mcset+"/"+args.inputFile)
       if not (os.path.isfile(path)):
           print "File "+path+" does not exist"
           continue
       cmd+=' '+path
   os.system(cmd)

print setdict
for dset in setdict:
   add_sets(dset, setdict[dset])
