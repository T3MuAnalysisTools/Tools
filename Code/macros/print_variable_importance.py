import os
import argparse

#----- Argument parser --------
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", help="input file", default="logs/TMVAClassification_2016vars_A_twoGlobalTracker.log")
args = parser.parse_args()
#------------------------------

with open(args.input_file) as file_:
    lines = file_.readlines()

start_ = 0
end_ = 0

for iline, line in enumerate(lines):
    if 'Rank : Variable' in line: start_ = iline+2
    elif 'Destroy' in line: end_ = iline-1

print 'Varible\tImportance'
for i in range(start_,end_):
    line = lines[i].strip('\n')
    _ = line.split(':')
    print (_[2].split('>')[0]).replace(' (',''),'\t', _[6]
