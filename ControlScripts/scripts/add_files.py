import os

cmd = 'hadd TMVAInput_merged_mc.root'
for i in range(58,66):
    cmd+=(' Set_'+str(i)+'/TMVAInput.root')
os.system(cmd)
