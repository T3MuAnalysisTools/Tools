import os.path as path
import os

pathbase='/afs/cern.ch/work/m/mmadhu/Analysis/workdirZTT_PostFCFeb_20_2023/'

samplepath=''

for i in range(283):
  #print(i+1)
  setno=str(i+1)
  skimmed_filepath=pathbase+'Set_'+setno+'/SKIMMED_NTUP.root'
  input_file=pathbase+'Set_'+setno+'/Input.txt'
  ifSkimmedFileExists=path.isfile(skimmed_filepath)
  #print(ifSkimmedFileExists)
  
  fp= open(input_file, 'r')
  current_samplepath=''
  for i, line in enumerate(fp):
        if 'InputNtuples' in line:
            current_samplepath=line.split()[1]
            break
  fp.close()
  
  if samplepath=='' or samplepath!=current_samplepath:
    samplepath=current_samplepath
    
  trimmed_sample_path=samplepath.split('/')[5]
  
  #print(trimmed_sample_path)
  
  gfalcopy_cmd='(eval `scram unsetenv -sh`; gfal-copy '+skimmed_filepath+' davs://cmsio7.rc.ufl.edu:1094/store/user/mmadhu/test/'+trimmed_sample_path+'/Skimmed_'+setno+'.root )'
  #print(gfalcopy_cmd)
  #os.system(gfalcopy_cmd)
  
  if ifSkimmedFileExists:
    os.system(gfalcopy_cmd)
    print('File Copied')
  else:
    print('Skimmed File Doesn\'t Exist')