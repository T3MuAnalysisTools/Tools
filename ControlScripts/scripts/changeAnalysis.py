import os
for i in range(1,31):
   filename = 'Set_'+str(i)+'/Input.txt'
   ifile = open(filename,'r')
   tmp_file = open('Set_'+str(i)+'/tmp_Input.txt','w+')
   dimu_count = 0
   threemu_count = 0
   for line in ifile:
      line = line.strip('\n\r')
      if 'DimuTrk' in line:
         dimu_count += 1
         if (dimu_count<2): 
            tmp_file.write(line)
            tmp_file.write('\n')
      elif 'ThreeMu' in line:
         threemu_count += 1
         if (threemu_count<2):
            tmp_file.write(line)
            tmp_file.write('\n')
      else: 
         tmp_file.write(line)
         tmp_file.write('\n')
   tmp_file.close()
   ifile.close()
   cmd = 'mv Set_'+str(i)+'/tmp_Input.txt '+filename
   os.system(cmd)
