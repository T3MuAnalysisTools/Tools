
samplepath=''

set_no_int=1
setno=str(set_no_int)
set_dir=pathbase+'Set_'+setno+'/'
ifDirExists=os.path.isdir(set_dir)   #checks whether directory exists

while ifDirExists:
  skimmed_filepath=pathbase+'Set_'+setno+'/SKIMMED_NTUP.root'
  input_file=pathbase+'Set_'+setno+'/Input.txt'
  ifSkimmedFileExists=os.path.isfile(skimmed_filepath)
  
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
  
  gfalcopy_cmd='(eval `scram unsetenv -sh`; gfal-copy '+skimmed_filepath+' davs://cmsio7.rc.ufl.edu:1094/store/user/'+username+'/'+ tag +'/'+trimmed_sample_path+'/Skimmed_'+setno+'.root )'
  
  if ifSkimmedFileExists:
    os.system(gfalcopy_cmd)
    print('File Copied to: '+username+'/'+ tag +'/'+trimmed_sample_path)
  else:
    print('Skimmed file for set ' + setno + ' doesn\'t exist')
  set_no_int=set_no_int+1
  setno=str(set_no_int)
  set_dir=pathbase+'Set_'+setno+'/'
  ifDirExists=os.path.isdir(set_dir)
