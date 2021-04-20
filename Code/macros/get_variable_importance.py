import os
import argparse
import numpy as np
from matplotlib import pyplot as plt


#----- Argument parser --------
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", help="input file", default="TMVAClassification_2016vars_A_twoGlobalTracker.log")
parser.add_argument("-o", "--output-file",help="output file", default="TMVAClassification_2016vars_twoGlobalTracker_importance.png")
args = parser.parse_args()
#------------------------------

plt.rcdefaults()
#plt.tight_layout()
fig, ax = plt.subplots()

def get_variable_importance(inputfile):

    with open(inputfile) as file_:
        lines = file_.readlines()

    start_ = 0
    end_ = 0

    var_list = []
    var_importance = []

    for iline, line in enumerate(lines):
        if 'Rank : Variable' in line: start_ = iline+2
        elif 'Destroy' in line: end_ = iline-1

    for i in range(start_,end_):
        line = lines[i].strip('\n')
        _ = line.split(':')
        varname = (_[2].split('>')[0]).replace(' (','')
        if 'min(Muon1_segmentCompatibility,Muon2_segmentCompatibility)' in varname:
            varname = 'min_muon_seg_comp'
        var_list.append(varname)
        var_importance.append(_[6])
    return (var_list, var_importance)


category_list = ['_A_', '_B_', '_C_']
input_cat = ''
for cat in category_list:
    if cat in args.input_file: input_cat = cat

files = {}

for cat in category_list:
    files[cat.replace('_','')] = args.input_file.replace(input_cat, cat)

importance_list = []
label_list = []
var_list = []

for cat_ in files:
    _, __ = get_variable_importance(files[cat_])
    var_list = _
    importance_list.append(__)
    label_list.append(cat_)

bar_size = 0.25
padding = 0.25

y_pos = np.arange(len(var_list)) * (bar_size * 3 + padding)

ax.barh(y_pos, importance_list[0], align ='center', label=label_list[0], color='blue', height=bar_size)
ax.barh(y_pos+bar_size, importance_list[1], align ='center', label=label_list[1], color='orange', height=bar_size)
ax.barh(y_pos+2*bar_size, importance_list[2], align ='center', label=label_list[2], color='gray', height=bar_size)
ax.set_yticks(y_pos)
ax.set_yticklabels(tuple(var_list))
ax.invert_yaxis()
ax.set(ylim=[0 - 2*padding, len(var_list)+padding])
ax.set_xlabel('Importance')
path_list = args.input_file.split('/')
title = path_list[len(path_list)-1].split('.')[0]
ax.set_title('BDT variable ranking (2017)')
plt.legend(fontsize=14, loc='upper right')
plt.tight_layout()
plt.savefig(args.output_file)
#plt.show()
