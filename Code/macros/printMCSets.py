import os, sys

mclist = []

setlist = [f for f in os.listdir('.') if 'Set_' in f]
for st in setlist:
   ifile = open(st+"/Input.txt")
   lines = ifile.readlines()
   datasetline = [l for l in lines if 'InputNtuples:' in l]
   datasetline = datasetline[0].strip('\r\n')
   print "Set: "+st+" | Dataset: "+datasetline
   if 'DoubleMuonLowMass' in datasetline:
      datatype = 'Data'
      print datasetline
   elif ( ('DsTau3Mu' in datasetline) or ('B0Tau3Mu' in datasetline) or ('BpTau3Mu' in datasetline) ): datatype = 'MC'
   else: datatype = None
   if (datatype=='MC'): mclist.append(st)

print mclist
