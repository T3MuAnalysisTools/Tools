import fileinput
corruptlist = open("unavailable_due_to_cmsufoss11_issue_cherepan_mmadhu_bjoshi.txt", 'r')

for rono, roline in enumerate(corruptlist):
    for setno in range(1, 396+1):#396 is number of sets
        setinputtext = fileinput.input("Set_"+str(setno)+"/Input.txt", inplace=True)
        setgetsh = fileinput.input("Set_"+str(setno)+"/Set_"+str(setno)+"-get.sh", inplace=True)
        
        inputfile = open("Set_"+str(setno)+"/Input.txt", 'r') #get the sample the set is derived from
        for i, lineinput in enumerate(inputfile):
            j=i+1
            if j == 13:
                sampleline=lineinput
        inputfile.close()
        
        sample=sampleline[14:(len(sampleline)-1)]  #get the sample the set is derived from
        
        #checking if sample of the set is the same as the particular sample of the current roline root file:
        
        if sample in roline:
            rootfilestr=roline[(len(roline)-17):(len(roline)-2)]  #get the root file from roline
            
            for line in setinputtext:
                if rootfilestr in line:
                    print('{} {}'.format("#", line), end='')
                else:
                    print('{}'.format(line), end='')
                    
            for line1 in setgetsh:
                if rootfilestr in line1:
                    print('{} {}'.format("#", line1), end='')
                else:
                    print('{}'.format(line1), end='')
            
    print(rono)
corruptlist.close()
