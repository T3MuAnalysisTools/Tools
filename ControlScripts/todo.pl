#! /usr/bin/perl
use Cwd;
use POSIX;
use POSIX qw(strftime);

#############################################
$numArgs = $#ARGV +1;
$ARGV[$argnum];

$UserID= POSIX::cuserid();
$UserDir="";
$UserName="";
if($UserID eq "cherepan"){
    $UserDir="--cherepanov";
    $UserName="Vladimir";
} 

if($UserID eq "bjoshi"){
    $UserDir="--joshi";
    $UserName="Bhargav";
}


if($UserID eq "wangjian"){
    $UserDir="--wang";
    $UserName="Jian";
}



if($UserID eq "nimenend"){
    $UserDir="--menendez";
    $UserName="Nicholas";
}

if($UserID eq "mmadhu"){
    $UserDir="--madhu";
    $UserName="Arun";
}




$letter = substr($UserID, 0, 1);

#Default values
$InputDir="/afs/cern.ch/work/$letter/$UserID/InputTest";
$OutputDir="/afs/cern.ch/work/$letter/$UserID/Analysis";
$CodeDir="../Code";
$set="ControlSample_";
$CMSSWRel="10_2_18";
$maxdata=20;
$maxmc=5;
$maxemb=20;
$ARCH="slc7_amd64_gcc700";
$PWD=getcwd;
$DsdevBranch = "master";

if($ARGV[0] eq "--help" || $ARGV[0] eq ""){
    printf("\n ===========>>>>  Bonjour $UserName ! <<<<=============\n\n");
    printf("\nThis code requires one input option. The syntax is:./todo.pl [OPTION]");
    printf("\nPlease choose from the following options:\n");
    printf("\n./todo.pl --help                                   Prints this message\n");
    printf("\n./todo.pl --DsTauTo3MNtuple <dir>                  Clone and compile DsToTau ntuple. Example: ./todo.pl --DsTauTo3MNtuple workdir  ");
    printf("\n                                                                        --Branch <branch> developing branch; Default: master ");
    printf("\n./todo.pl                                           --MuonPogNtuple <dir> MuonPogNtuple  ");
    printf("\n                                                    --ARCH  <SCRAM_ARCH>   Setup SCRAM_ARCH; Default: slc7_amd64_gcc700 ");
    printf("\n                                                    --CMSSWRel <X_Y_Z>  Configure CMSSW_X_Y_Z; Default: CMSSW_10_2_18  \n\n");
    printf("\n./todo.pl --Local <Input.txt>                      INTENTED FOR SMALL SCALE TESTS ONLY");  
    printf("\n                                                   Configure a directory to run locally. <InputPar.txt> name of file that");
    printf("\n                                                   contains input command template.");
    printf("\n                                                   Optional commands:  ");
    printf("\n                                                     --InputDir <InputDir>   Default value: $InputDir"); 
    printf("\n                                                       Note: the root files must be inside a $InputDir/<dir>/ ");
    printf("\n                                                       where dir contains user, data or mc as a substring");
    printf("\n                                                     --OutputDir <OutputDir> Default value: $OutputDir");
    printf("\n                                                     --SetName <SetName>     Default value: $set ");
    printf("\n                                                     --NMaxData <Max Number of data files per job >     Default value: $maxdata ");
    printf("\n                                                     --NMaxMC <Max Number of MC files per job >     Default value: $maxmc ");
    printf("\n                                                     --Proxy <path>;   a path to proxy recieved by running grid-init");
    printf("\n./todo.pl --DCache <Input.txt> <ListofDS.txt>      INTENTED FOR REGULAR USE (DEFAULT)");
    printf("\n                                                   Configure a directory to run from. <InputPar.txt> name of file that");
    printf("\n                                                   contains input command template.");
    printf("\n                                                   Please make sure that you have run voms-proxy-init -voms cms -valid 192:00");
    printf("\n                                                   beforehand and copy your received proxy to proxy dir");
    printf("\n                                                   <ListoDS.txt> list of DCache Dataset directories you want to run on.");
    printf("\n                                                   Optional commands:  ");
    printf("\n                                                     --OutputDir <OutputDir> Default value: $OutputDir ");
    printf("\n                                                     --SetName <SetName>     Default value: $set  "); 
    printf("\n                                                     --NMaxData <Max Number of data files per job >     Default value: $maxdata ");
    printf("\n                                                     --NMaxMC <Max Number of MC files per job >     Default value: $maxmc ");
    printf("\n                                                     --Proxy <path>;   a path to proxy recieved by running grid-init");
    printf("\n  ");


    exit(0);  
} 

######################################
$InputFile=$ARGV[1];
$buildRoot=0;
$hasroot=0;
for($l=2;$l<$numArgs; $l++){
    if($ARGV[$l] eq "--InputDir"){
	$l++;
	$InputDir=$ARGV[$l];
    }
    if($ARGV[$l] eq "--OutputDir"){
	$l++;
	$OutputDir=$ARGV[$l];
    }
    if($ARGV[$l] eq "--SetName"){
	$l++;
	$set=$ARGV[$l];
    }
    if($ARGV[$l] eq "--CMSSWRel"){
        $l++;
        $CMSSWRel=$ARGV[$l];
    }
    if($ARGV[$l] eq  "--NMaxData" ){
	$l++;
	$maxdata=$ARGV[$l];
    }
    if($ARGV[$l] eq  "--NMaxMC" ){
	$l++;
	$maxmc=$ARGV[$l];
    }

    if($ARGV[$l] eq  "--ARCH" ){
        $l++;
        $ARCH=$ARGV[$l];
    }
    if($ARGV[$l] eq  "--Proxy" ){
        $l++;
        $Proxy=$ARGV[$l];
    }
}
my $dir = getcwd;

$time= strftime("%h_%d_%Y",localtime);
$temp= $set . $time;
$set=$temp;



if( $ARGV[0] eq "--MuonPogNtuple"){
    $currentdir=getcwd;
    if($ARGV[1] ne ""){
	$basedir=$ARGV[1];
    }
    else{
	printf("\nWorkingDir for CMSSW is required. Please follow the syntax:./todo.pl --MuonPogNtuple <dir> ");
	printf("\nFor more details use: ./todo --help\n"); 
	exit(0);
    }

    printf("\nWorkingDir for CMSSW: $basedir");
    printf("\nCurrentDir is: $currentdir \n");
    $CMSPATH="/CMSSW_$CMSSWRel";
    $CMSSW_BASE="$basedir$CMSPATH";

    system(sprintf("rm Install_MuPoGNtuple_$time"));
    system(sprintf("echo \"export SCRAM_ARCH=\\\"$ARCH\\\"\" >> Install_MuPoGNtuple_$time"));
    system(sprintf("echo \"source /cvmfs/cms.cern.ch/cmsset_default.sh\" >> Install_MuPoGNtuple_$time"));

    system(sprintf("echo \"mkdir $basedir\" >>  Install_MuPoGNtuple_$time")); 
    system(sprintf("echo \"cd $basedir\" >>  Install_MuPoGNtuple_$time")); 
    system(sprintf("echo \"cmsrel CMSSW_$CMSSWRel\" >>  Install_MuPoGNtuple_$time")); 
    system(sprintf("echo \"cd CMSSW_$CMSSWRel/src\" >> Install_MuPoGNtuple_$time")); 
    system(sprintf("echo \"cmsenv\" >> Install_MuPoGNtuple_$time")); 
    system(sprintf("echo \"git clone git\@github.com:T3MuAnalysisTools/MuonPOGtreeProducer.git\" >> Install_MuPoGNtuple_$time"));
    system(sprintf("echo \"cd MuonPOGtreeProducer; git checkout 9X_tau3mu;\" >> Install_MuPoGNtuple_$time"));
    system(sprintf("echo \"cd $currentdir/$CMSSW_BASE/src\" >> Install_MuPoGNtuple_$time")); 
    system(sprintf("echo \"scram b -j 4\" >> Install_MuPoGNtuple_$time")); 

    printf("\n\nInstructions:");
    printf("\nsource  Install_MuPoGNtuple_$time to complete installation, compilation might take some time...  ");
    printf("\n\n\nTo run test job do  'cmsRun muonPogNtuples_cfg.py'  in  $CMSSW_BASE/src/MuonPOGtreeProducer/Tools/test/ \n\n");
}

if( $ARGV[0] eq "--DsTauTo3MNtuple"){

    # User check
    if($UserID ne "cherepan" and $UserID ne "wangjian" and $UserID ne "bjoshi" and $UserID ne "nimenend"  and $UserID ne "mmadhu"){
	printf("\nUnrecognized user. Exit.\n"); 
	exit(0);
    }

    $currentdir=getcwd;
    if($ARGV[1] ne ""){
        $basedir=$ARGV[1];
    }
    else{
        printf("\nWorkingDir for CMSSW is required. Please follow the syntax:./todo.pl --DsTauTo3MNtuple <dir> ");
        printf("\nFor more details use: ./todo --help\n");
        exit(0);
    }
    for($l=2;$l<$numArgs; $l++){
	if($ARGV[$l] eq "--Branch"){
	    $l++;
	    $DsdevBranch=$ARGV[$l];
	}
    }

    printf("\nWorkingDir for CMSSW: $basedir");
    printf("\nCurrentDir is: $currentdir \n");
    $CMSPATH="/CMSSW_$CMSSWRel";
    $CMSSW_BASE="$basedir$CMSPATH";

    system(sprintf("rm Install_DsTNtuple_$time"));
    system(sprintf("echo \"export SCRAM_ARCH=\\\"$ARCH\\\"\" >> Install_DsTNtuple_$time"));
    system(sprintf("echo \"source /cvmfs/cms.cern.ch/cmsset_default.sh\" >> Install_DsTNtuple_$time"));

    system(sprintf("echo \"mkdir $basedir\" >>  Install_DsTNtuple_$time"));
    system(sprintf("echo \"cd $basedir\" >>  Install_DsTNtuple_$time"));
    system(sprintf("echo \"cmsrel CMSSW_$CMSSWRel\" >>  Install_DsTNtuple_$time"));
    system(sprintf("echo \"cd CMSSW_$CMSSWRel/src\" >> Install_DsTNtuple_$time"));
    system(sprintf("echo \"cmsenv\" >> Install_DsTNtuple_$time"));
    system(sprintf("echo \"git clone git\@github.com:T3MuAnalysisTools/DsTau23Mu.git\" >> Install_DsTNtuple_$time"));
    system(sprintf("echo \"cd DsTau23Mu; git checkout $DsdevBranch; \" >> Install_DsTNtuple_$time"));


    system(sprintf("echo \"cd ../\" >> Install_DsTNtuple_$time"));
    system(sprintf("echo \"git clone git\@github.com:T3MuAnalysisTools/SkimProduction.git\" >> Install_DsTNtuple_$time"));
    system(sprintf("echo \"scram b -j 4\" >> Install_DsTNtuple_$time"));

    printf("\n\nInstructions:");
    printf("\nTo complete the installation do the following command (compilation might take some time ...):");
    printf("\nsource  Install_DsTNtuple_$time");
    printf("\n\n\nTo run test job do  'cmsRun analyze.py'  in  $CMSSW_BASE/src/DsTau23Mu/T3MNtuple/test/ \n\n");
}



if( $ARGV[0] eq "--Local" ){

    # User check
    if($UserID ne "cherepan" and $UserID ne "wangjian" and $UserID ne "bjoshi" and $UserID ne "nimenend" and $UserID ne "mmadhu"){
	printf("\nUnrecognized user. Exit.\n"); 
	exit(0);
    }

    # Print out input parameters
    printf("Active directory will be: $OutputDir/workdir$set \n");
    printf("Input directory will be:  $InputDir \n");
    printf("Code Repository is:       $CodeDir \n");
    # Clean Directory in case workdir$set exists
    printf("Cleaning Directories \n");
    system(sprintf("cd $OutputDir"));
    system(sprintf("rm -rf $OutputDir/workdir$set \n"));
    system(sprintf("mkdir $OutputDir/workdir$set "));
    printf("Cleaning complete \n");
    
    # create directory stucture
    system(sprintf("mkdir  $OutputDir/workdir$set/Code "));
    system(sprintf("mkdir  $OutputDir/workdir$set/Code/libs "));
    system(sprintf("cp -r $CodeDir/* $OutputDir/workdir$set/Code/ "));
    system(sprintf("mkdir $OutputDir/workdir$set/EPS "));
    system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/InputData "));
    
    # generate init script 
    system(sprintf("echo \"cd  $OutputDir/workdir$set/Code/CommonUtils/CMSSW_9_3_8/src\" >> $OutputDir/workdir$set/init.sh "));
    system(sprintf("echo \"cmsenv \" >> $OutputDir/workdir$set/init.sh "));




    # generate compile script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd  $OutputDir/workdir$set/Code/\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"./config \\\$\@ \"   >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"source $OutputDir/workdir$set/init.sh \"   >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd  $OutputDir/workdir$set/Code/\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"gmake all \" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/compile")) ;
 
    # Generate Combine script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Combine")) ;
    system(sprintf("echo \"export workdir=\\\"$OutputDir/workdir$set/\\\"\" >> $OutputDir/workdir$set/Combine"));
    system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; ./config \" >> $OutputDir/workdir$set/Combine")); 
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Combine")) ; 
    system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe \" >> $OutputDir/workdir$set/Combine")) ;

    # Generate Combine Input
    system(sprintf("cp $InputFile $OutputDir/workdir$set/Input.txt "));
    system(sprintf("cd $OutputDir/workdir$set; $dir/subs '{SET}' COMBINE Input.txt; cd $dir "));
    system(sprintf("cd $OutputDir/workdir$set; $dir/subs '{FileDir}' NotAvailable Input.txt; cd $dir "));
    system(sprintf("echo \"Mode: RECONSTRUCT\" >> $OutputDir/workdir$set/Input.txt"));
    system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Input.txt"));

    # Setup Condor Combine scripts
    system(sprintf("echo \"universe     = vanilla      \"  >> $OutputDir/workdir$set/Condor_Combine"));
    system(sprintf("echo \"rank         = memory       \"  >> $OutputDir/workdir$set/Condor_Combine"));
    system(sprintf("echo \"executable   = Combine      \"  >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"output       = Combine-Condor_\\\$(cluster)_\\\$(proccess).o  \" >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"error        = Combine-Condor_\\\$(cluster)_\\\$(proccess).e  \" >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"log          = Combine-Condor_\\\$(cluster)_\\\$(proccess).log  \" >> $OutputDir/workdir$set/Condor_Combine")); 			
    system(sprintf("echo \"getenv       = True\" >> $OutputDir/workdir$set/Condor_Combine"));
    system(sprintf("echo \"+JobFlavour = 'microcentury'\" >> $OutputDir/workdir$set/Condor_Combine"));
    system(sprintf("echo \"queue  1 \" >> $OutputDir/workdir$set/Condor_Combine"));

    # Start Submit script
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Submit")) ; 
    system(sprintf("chmod u+x $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"rm Set*/*.o; rm Set*/*.e; rm Set*/*.log; \" >> $OutputDir/workdir$set/Submit")) ;

    opendir(DIR,"$InputDir/");
    @dirs = grep {( /user/ || /data/ || /mc/)} readdir(DIR);
    closedir DIR;
    
    $B=0;
    $max=1;
    for($l=0;$l<2; $l++){
	printf("\n\nStarting Loop $l  maxdata = $maxdata maxmc = $maxmc maxemb = $maxemb \n");
	$A=$maxdata+$maxmc+$maxemb+10;
	foreach $subdir (@dirs){
	    printf("subdir  $subdir\n");
	    if(($l==0 && ($subdir =~ m/data/)) || ($l==1 && !($subdir =~ m/data/))){
		if($l==0){
		    $max=$maxdata;
		}
		else{
		    $max=$maxmc;
			if($subdir =~ m/embed/){
				$max=$maxemb
			}
		}
		printf(" \nAccessing Directory  $subdir \n");
		opendir(SUBDIR,"$InputDir/$subdir/");
		@files = grep { /root/ } readdir(SUBDIR);

		$nfiles = @files;
		$idx=0;
		foreach $file (@files){
		    $idx++;
		    printf("$file Set= $B  Index=  $A   Max.= $max N Files = $nfiles Current File = $idx \n");
		    if($A > $max){
			$A=1;
			$B++;
			# Add Set information to Combining scripts and Input.txt
			system(sprintf("echo \"File: $OutputDir/workdir$set/Set_$B/  \" >>  $OutputDir/workdir$set/Input.txt ")) ;
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B \" >> $OutputDir/workdir$set/Submit")) ;
			system(sprintf("echo \"condor_submit  Condor_Set_$B  \" >> $OutputDir/workdir$set/Submit")) ;

			# Create and configure Set_$B dir
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B ")) ;
			system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/Set_$B/InputData "));
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B/EPS ")) ;

                        # Setup Set_$B.sh
			system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;
			system(sprintf("echo \"echo 'Starting Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"export workdir=\\\"$OutputDir/workdir$set/\\\"\" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; ./config \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;
			system(sprintf("chmod +x $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe 2>&1 | tee >(sed -r \\\"s/\\\\x1B\\\\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g\\\" > Set_$B.output) \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"echo 'Science Completed' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));

			# Setup Input.txt
			system(sprintf("cp   $InputFile $OutputDir/workdir$set/Set_$B/Input.txt ")); 
			system(sprintf("cd $OutputDir/workdir$set/Set_$B; $dir/subs '{SET}' Set_$B Input.txt; cd $dir "));
			system(sprintf("cd $OutputDir/workdir$set/Set_$B; $dir/subs '{FileDir}' $InputDir/$subdir Input.txt; cd $dir ")); 
			system(sprintf("echo \"Mode: ANALYSIS\" >> $OutputDir/workdir$set/Set_$B/Input.txt")); 
			system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Set_$B/Input.txt"));
			# Setup Condor scripts
			system(sprintf("echo \"universe     = vanilla      \"  >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
			system(sprintf("echo \"rank         = memory       \"  >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
			system(sprintf("echo \"executable   = Set_$B.sh      \"  >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"output       = Set_$B-Condor_\\\$(cluster)_\\\$(proccess).o  \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"error        = Set_$B-Condor_\\\$(cluster)_\\\$(proccess).e  \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"log          = Set_$B-Condor_\\\$(cluster)_\\\$(proccess).log  \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"getenv       = True\" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
			system(sprintf("echo \"notification = Error        \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));		
			system(sprintf("echo \"queue  1 \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));

	
		    }
		    system(sprintf("echo \"File: $InputDir/$subdir/$file  \" >> $OutputDir/workdir$set/Set_$B/Input.txt")) ;
		    $A++;
		}
	    }
	}
    }
    
    printf("Setting up root....\n");
    if($buildRoot==1){
        printf("Building custom root $buildrootversion... Enjoy your coffee!!! ");
	system(sprintf("wget ftp://root.cern.ch/root/root_v$buildrootversion.source.tar.gz"));
	system(sprintf("gzip -dc root_v$buildrootversion.source.tar.gz | tar -xf -"));
	system(sprintf("mkdir $OutputDir/workdir$set/root"));
	system(sprintf("cd root_v$buildrootversion; ./configure --enable-python --enable-roofit --enable-minuit2 --disable-xrootd --disable-sqlite --disable-python --disable-mysql --prefix=$OutputDir/workdir$set/root; make & make install "));
    }

    # Finish Submit script
    system(sprintf("echo \"cd  $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit"));

    # print Instructions  
    printf("\n\nInstructions");
    printf("\nPlease make sure you have run:");
    printf("\nvoms-proxy-init -voms cms -valid 192:00");
    printf("\nNow you can run the analysis.");
    printf("\nTo go to the Test workdir: cd  $OutputDir/workdir$set ");
    printf("\nTo compile the code in the workdir: source compile   $UserDir");
    printf("\nTo submit jobs to the batch queue: source Submit ");
    printf("\nTo combine jobs submitted to the batch queue: source Combine \n");
    printf("\nTo test a single job: cd  $OutputDir/workdir$set; source compile   $UserDir; cd $OutputDir/workdir$set/Set_1; ./Set_1 | tee log; cd ..\n");

}





if( $ARGV[0] eq "--DCache" ){

    if($UserID ne "cherepan" and $UserID ne "wangjian" and $UserID ne "bjoshi" and $UserID ne "nimenend"  and $UserID ne "mmadhu" ){
	printf("\nUnrecognized user. Exit.\n"); 
	exit(0);
    }
    $RemoteScrathDir="/afs/cern.ch/work/$letter/$UserID";
    $RemoteDir='\$TMPDIR'; 
    $TempDataSetFile=$ARGV[2];
    # Print out input parameters
    printf("Active directory will be: $OutputDir/workdir$set \n");
    printf("Code Repository is:       $CodeDir \n");
    printf("List of dcache dir:       $TempDataSetFile \n");

    # Open ListofFile.txt
    @DataSets;
    open(DAT, $TempDataSetFile) || die("Could not open file $TempDataSetFile!");
    while ($item = <DAT>) {
	chomp($item);
	push(@DataSets,$item);
    }
    close(DAT);
    # Clean Directory in case workdir$set exists
    printf("Cleaning Directories \n");
    system(sprintf("cd $OutputDir"));
    system(sprintf("rm -rf $OutputDir/workdir$set \n"));
    system(sprintf("mkdir $OutputDir/workdir$set "));
    printf("Cleaning complete \n");
    
    # creat directory stucture
    system(sprintf("mkdir  $OutputDir/workdir$set/Code "));
    system(sprintf("mkdir  $OutputDir/workdir$set/Code/libs "));
    system(sprintf("cp -r $CodeDir/* $OutputDir/workdir$set/Code/ "));
    system(sprintf("mkdir $OutputDir/workdir$set/EPS "));
    system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/InputData "));

    # generate init script
    system(sprintf("echo \"cd  $OutputDir/workdir$set/Code/CommonUtils/CMSSW_9_3_8/src\" >> $OutputDir/workdir$set/init.sh "));
    system(sprintf("echo \"cmsenv \" >> $OutputDir/workdir$set/init.sh "));
    system(sprintf("echo \"cd -\" >> $OutputDir/workdir$set/init.sh "));

    # generate compile script
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd  $OutputDir/workdir$set/Code/\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"./config \\\$\@ \"   >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"source $OutputDir/workdir$set/init.sh \"   >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd  $OutputDir/workdir$set/Code/\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"gmake all \" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/compile")) ;
 
    # copy helpers
    system(sprintf("cp scripts/run_rest.py  $OutputDir/workdir$set/"));


    # Generate Combine script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Combine")) ;
    system(sprintf("echo \"export workdir=\\\"$OutputDir/workdir$set/\\\"\" >> $OutputDir/workdir$set/Combine"));
    system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; ./config \" >> $OutputDir/workdir$set/Combine"));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Combine")) ; 
    system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe \" >> $OutputDir/workdir$set/Combine")) ;

    # Generate Combine Input
    system(sprintf("cp $InputFile $OutputDir/workdir$set/Input.txt "));
    system(sprintf("cd $OutputDir/workdir$set/; $dir/subs '{SET}' COMBINE Input.txt; cd $dir"));
    system(sprintf("cd $OutputDir/workdir$set/; $dir/subs '{FileDir}' COMBINE Input.txt; cd $dir"));
    system(sprintf("echo \"Mode: RECONSTRUCT\" >> $OutputDir/workdir$set/Input.txt"));
    system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Input.txt"));

    # Generate SubmitManual and set_env 
    system(sprintf("cp  subs  $OutputDir/workdir$set/;"));
    system(sprintf("cp  SubmitManual  $OutputDir/workdir$set/;"));

    # Setup Condor Combine scripts
    system(sprintf("echo \"universe     = vanilla      \"  >> $OutputDir/workdir$set/Condor_Combine"));
    system(sprintf("echo \"rank         = memory       \"  >> $OutputDir/workdir$set/Condor_Combine"));
    system(sprintf("echo \"executable   = Combine      \"  >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"output       = Combine-Condor_\\\$(cluster)_\\\$(proccess).o  \" >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"error        = Combine-Condor_\\\$(cluster)_\\\$(proccess).e  \" >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"log          = Combine-Condor_\\\$(cluster)_\\\$(proccess).log  \" >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"getenv       = True\" >> $OutputDir/workdir$set/Condor_Combine"));
    system(sprintf("echo \"queue   1 \" >> $OutputDir/workdir$set/Condor_Combine"));

    # Start Submit script
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Submit")) ; 
    system(sprintf("chmod u+x $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"verbosity=\\\$(grep SetLevel Code/Analysis.cxx | grep -c -e Debug -e Verbose)\" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"  if [[ \\\${verbosity} -ne 0 ]]; then \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"    echo 'ERROR: Please make sure to set the verbosity level to Info in Analysis.cxx, otherwise your log-files will break QSUB! Abort...' \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"    exit \\\${verbosity}\" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"  fi\" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"rm Set*/*.o; rm Set*/*.e; rm Set*/*.log; \" >> $OutputDir/workdir$set/Submit")) ;
 


    $B=0;
    for($l=0;$l<2; $l++){
	system(sprintf("echo \"notification = Error        \" >> $OutputDir/workdir$set/Condor_Combine"));
	$max=1;
	foreach $DS (@DataSets){
	    print $DS; 	       
	    if(($l==0 && ($DS =~ m/DoubleMuonLowMass/)) || ($l==1 && !($DS =~ m/DoubleMuonLowMass/))){
		if($l==0){
		    $max=$maxdata;
		}
		else{
		    $max=$maxmc;
		    if($DS =~ m/embed/){
			$max=$maxemb
		    }
		}
		print "  max  = $max;    maxmc = $maxmc \n";
		printf("\n\nStarting Loop $l \n");
		$A=$maxdata+$maxmc+$maxemb+10;

		# find the root files for the current DataSet (DS)
		printf("Accessing Directory  $DS \n");
		system(sprintf("touch junk0"));
		system(sprintf("touch junk1"));
		system(sprintf("touch junk2"));
		system(sprintf("touch junk"));

		system(sprintf("uberftp cmsio.rc.ufl.edu \"ls /cms/data$DS/ \" | grep .root >&junk0"));
		system(sprintf("cat junk0 | awk '{print \$9}' >& junk1")); 


		system(sprintf("junk1")); 


		# Get list of files in dcache dir
		@files=();
		open(DAT, "junk1");
		while ($item = <DAT>) {
		    $item =~ s/\r\n$/\n/;
		    chomp($item);
		    push(@files,$item);
		}
		close(DAT);


		system(sprintf("rm junk0"));
		system(sprintf("rm junk1"));
		system(sprintf("rm junk2"));
		system(sprintf("rm junk"));

		$nfiles = @files;
		$idx=0;
		printf("A = $A  and B=$B \n ");
		foreach $file (@files){
		    $idx++;
		    printf("$file Set = $B  Index =  $A   Max. = $max    N Files = $nfiles   Current File = $idx \n");
		    if($A > $max ){
			$A=1;
			$B++;
			$Filedir="$RemoteDir/Set_$B";
			# Add Set information to Combining scripts and Input.txt
			system(sprintf("echo \"File: $OutputDir/workdir$set/Set_$B/ \" >>  $OutputDir/workdir$set/Input.txt ")) ;
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B \" >> $OutputDir/workdir$set/Submit")) ;
			system(sprintf("echo \"condor_submit  Condor_Set_$B  \" >> $OutputDir/workdir$set/Submit")) ;


			# Create and configure Set_$B dir
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B ")) ;
			system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/Set_$B/InputData "));
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B/EPS ")) ;
			

			# Setup Set_$B.sh
			system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;
			system(sprintf("echo \"echo 'Starting Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"export workdir=\\\"$OutputDir/workdir$set/\\\"\" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"export X509_USER_PROXY=\\\"$Proxy\\\"\"  >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; ./config \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"source $OutputDir/workdir$set/Set_$B/Set_$B-get.sh \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ; 
			system(sprintf("chmod +x $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"mkdir $RemoteDir/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cp -r * $RemoteDir/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd  $RemoteDir/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe 2>&1 | tee >(sed -r \\\"s/\\\\x1B\\\\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g\\\" > Set_$B.output) \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cp -r *  $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"source $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"rm -r   $RemoteDir/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"export HOME=\\\"/afs/cern.ch/user/c/$UserID\\\"         \"   >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));

			system(sprintf("echo \"echo 'Completed Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));





                        # Setup Set_$B_get.sh and Set_$B_clean.sh
			system(sprintf("echo \"#! /bin/bash\"         >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
			system(sprintf("echo \"mkdir $Filedir \" >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
			system(sprintf("echo \"cd $Filedir \"    >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
			system(sprintf("echo \"#! /bin/bash\"         >> $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh"));
			system(sprintf("echo \"cd $RemoteDir \"    >> $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh")); 
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;
			system(sprintf("echo \"echo 'Cleaning EPS' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"rm EPS/*eps  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));

			# Setup Input.txt
			system(sprintf("cp $InputFile $OutputDir/workdir$set/Set_$B/Input.txt ")); 
			system(sprintf("cd $OutputDir/workdir$set/Set_$B; $dir/subs '{SET}' Set_$B Input.txt; cd $dir "));
			system(sprintf("cd $OutputDir/workdir$set/Set_$B; $dir/subs '{FileDir}' $DS Input.txt; cd $dir "));
			system(sprintf("echo \"Mode: ANALYSIS\" >> $OutputDir/workdir$set/Set_$B/Input.txt")); 
			system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Set_$B/Input.txt"));


			system(sprintf("echo \"universe     = vanilla      \"  >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
			system(sprintf("echo \"rank         = memory       \"  >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
			system(sprintf("echo \"executable   = Set_$B.sh      \"  >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"output       = Set_$B-Condor_\\\$(cluster)_\\\$(proccess).o  \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"error        = Set_$B-Condor_\\\$(cluster)_\\\$(proccess).e  \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"log          = Set_$B-Condor_\\\$(cluster)_\\\$(proccess).log  \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"getenv       = True\" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
			system(sprintf("echo \"notification = Error        \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
			system(sprintf("echo \"queue  1 \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
		    
			
		    }
		    ($a,$b,$c)=split('/',$file);
		    $myfile=$file;
		    if($a =~ m/root/){
			$myfile=$a;
		    }
		    if($b =~ m/root/){
                        $myfile=$b;
			system(sprintf("echo \"notification = Error        \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
                    }
		    if($c =~ m/root/){
                        $myfile=$c;
                    }
		    $myfiletrunc = $myfile;

		    my @wholepath = split /\//, $file;
		    foreach $TreeName (@wholepath){
			if($TreeName =~ m/root/){
			        $myfiletrunc=$TreeName
			}
			$myTestpp +=$TreeName
		    }

		    $ind = rindex($file,"/");
		    $mypath = substr($file,0,$ind);

#uberftp  cmsio.rc.ufl.edu 'cd /cms/data//store/user/cherepan/DoubleMuonLowMass/Prod_07_10_2019_DoubleMuonLowMass__Run2017F-17Nov2017-v1/191007_093213/0000; get DsT3MNtuple_1.root'
		    system(sprintf("echo \"    uberftp  cmsio.rc.ufl.edu 'cd /cms/data/$DS/; get $myfiletrunc ' \"  >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
	#	    system(sprintf("echo \"gfal-copy gsiftp://cmsio.rc.ufl.edu/cms/data/$file . \"  >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
		    system(sprintf("echo \" File:  $Filedir/$myfiletrunc \"     >> $OutputDir/workdir$set/Set_$B/Input.txt")) ;
		    system(sprintf("echo \"rm -rf $Filedir/$myfiletrunc \"    >> $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh"));




		    $A++;
		}
	    }
	}
    }
    system(sprintf("echo \"cd  $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit"));

    # print Instructions
    printf("\n\nInstructions");
    printf("\nPlease make sure you have run:");
    printf("\nvoms-proxy-init -voms cms -valid 192:00"); 
    printf("\nNow you can run the analysis using dcache.");
    printf("\nTo go to the Test workdir: cd  $OutputDir/workdir$set ");
    printf("\nTo compile the code in the workdir: source compile  $UserDir ");
    printf("\nTo submit jobs to the batch queue: source Submit ");
    printf("\nTo combine jobs submitted to the batch queue: source Combine \n");
    printf("\nTo test a single job: cd  $OutputDir/workdir$set; source compile   $UserDir; cd $OutputDir/workdir$set/Set_1; source Set_1 | tee log; cd ..\n");

}
