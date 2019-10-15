#!/usr/bin/env python

import os
import sys
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--input-file",help="input file; [Default: %(default)s] ", action="store", default = 'datasets.dat')

    args = parser.parse_args()

    flist = args.input_file
    if not os.path.isfile(flist):
        print "File %s not found!!!" % flist
        sys.exit()


    with open(flist) as filelist:
        for l in filelist:
            l = l.strip()
            for content in l.split("/"):
                if "Set_" in content:
                    command = "cd "+ l + ";"+ "condor_submit Condor_"+content+"; cd ../;"
                    os.system(command)

