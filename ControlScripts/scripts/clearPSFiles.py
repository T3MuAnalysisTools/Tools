import os, sys
import argparse

#----- Argument parser --------
parser = argparse.ArgumentParser(description="Clear post-script files from a directory and EPS sub-directory.")

parser.add_argument("--inputDir", help="path to the directory")
parser.add_argument("--verbose", help="verbose option", action="store_true")
args = parser.parse_args()
#------------------------------

VERBOSE = args.verbose
DIR_PATH = args.inputDir
EPS_PATH = DIR_PATH+'/EPS'

if not (os.path.exists(args.inputDir)):
    print "Path "+args.inputDir+" doesn't exist!"
    parser.print_help()
    sys.exit()

if (os.path.exists(DIR_PATH)):
   if (VERBOSE): print "Removing post-script files ("+DIR_PATH+')...'
   out = os.popen("ls "+DIR_PATH).read()
   files = [item for item in out.split('\n') if (item!='\n' and '.ps' in item)]
   if (len(files)==0):
      if (VERBOSE): print "No post-script files found!"
   else:
      for f in files: 
         if (VERBOSE): print("* "+f)
         #os.system('rm '+DIR_PATH+'/'+f)

if (os.path.exists(EPS_PATH)):
   if (VERBOSE): print "Removing files from EPS directory ("+EPS_PATH+')...'
   out = os.popen("ls "+EPS_PATH).read()
   files = [item for item in out.split('\n') if (item!='\n' and '.eps' in item)]
   if (len(files)==0):
      if (VERBOSE): print "No files found in EPS directory!"
   else: 
      os.system("rm "+EPS_PATH)

setlst = [ f for f in os.listdir(DIR_PATH) if "Set_" in f]

for s in setlst:    
    # set paths
    DIR_PATH = args.inputDir+'/'+s
    EPS_PATH = DIR_PATH+'/EPS'

    if (os.path.exists(DIR_PATH)):
       if (VERBOSE): print "Removing post-script files from "+DIR_PATH+' ...'
       out = os.popen("ls "+DIR_PATH).read()
       files = [item for item in out.split('\n') if (item!='\n' and '.ps' in item)]
       if (len(files)==0):
          if (VERBOSE): print "No .ps files found!"  
       else:
          for f in files: print("* "+f)
    
    if (os.path.exists(EPS_PATH)):
       if (VERBOSE): print "Removing eps files from "+EPS_PATH+" ..."
       out = os.popen('ls '+EPS_PATH).read()
       files = [item for item in out.split('\n') if (item!='\n' and '.eps' in item)]
       if (len(files)==0):
          if(VERBOSE): print "No .eps files found!"
       else:
          if (VERBOSE): print "Clearing "+EPS_PATH+" directory ..."
          os.system('rm '+EPS_PATH+'/*')
