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

mclistBu = []
mclistBd = []
mclistDs = []
datalistF = []
datalistE = []
datalistD = []
datalistC = []
datalistB = []

setlist = [f for f in os.listdir(SRCDIR) if 'Set_' in f]
for st in setlist:
   ifile = open(SRCDIR+st+"/Input.txt")
   lines = ifile.readlines()
   datasetline = [l for l in lines if 'InputNtuples:' in l]
   datasetline = datasetline[0].strip('\r\n')
   if ('DoubleMuonLowMass' in datasetline):
      if ('2017B' in datasetline): datatype = 'DataB'
      elif ('2017C' in datasetline): datatype = 'DataC'
      elif ('2017D' in datasetline): datatype = 'DataD'
      elif ('2017E' in datasetline): datatype = 'DataE'
      elif ('2017F' in datasetline): datatype = 'DataF'
   elif ('DsToTau' in datasetline): datatype = 'DsMC'
   elif ('BuToTau' in datasetline): datatype = 'BuMC' 
   elif ('BdToTau' in datasetline): datatype = 'BdMC'
   #elif ( 'DsToTau' in datasetline ): datatype = 'MC'
   else: datatype= None
   if (datatype=='DsMC'): mclistDs.append(st)
   elif (datatype=='BuMC'): mclistBu.append(st)
   elif (datatype=='BdMC'): mclistBd.append(st)
   elif (datatype=='DataB'): datalistB.append(st)
   elif (datatype=='DataC'): datalistC.append(st)
   elif (datatype=='DataD'): datalistD.append(st)
   elif (datatype=='DataE'): datalistE.append(st)
   elif (datatype=='DataF'): datalistF.append(st)
   ifile.close()

if (args.instanceName==None):
   print "Please provide the name of the class."

def add_sets(datatype, setlist):
   cmd = 'hadd -f '+args.instanceName+'_'+datatype+'_merged.root'
   for mcset in setlist:
       path=(SRCDIR+mcset+"/"+args.inputFile+".root")
       if not (os.path.isfile(path)):
           print "File "+path+" does not exist"
           continue
       cmd+=' '+path
   os.system(cmd)

setdict = {'DsMC': mclistDs, 'BuMC': mclistBu, 'BdMC': mclistBd,
            '2017B': datalistB, '2017C': datalistC, '2017D': datalistD, '2017E': datalistE, '2017F': datalistF }

for dtype in setdict: add_sets(dtype, setdict[dtype])
