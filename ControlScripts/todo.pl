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




$letter = substr($UserID, 0, 1);

#Default values
$InputDir="/afs/cern.ch/work/$letter/$UserID/InputTest";
$OutputDir="/afs/cern.ch/work/$letter/$UserID/Test";
$CodeDir="../Code";
$set="ControlSample_";
$CMSSWRel="9_4_4";
$Cleaning ="NO";
$maxdata=20;
$maxmc=5;
$maxemb=20;
$ARCH="slc6_amd64_gcc630";
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
    printf("\n                                                    --ARCH  <SCRAM_ARCH>   Setup SCRAM_ARCH; Default: slc6_amd64_gcc630 ");
    printf("\n                                                    --CMSSWRel <CMSSW_X_Y_Z>  Configure CMSSW_X_Y_Z; Default: CMSSW_9_4_4  \n\n");
    printf("\n./todo.pl --Local <Input.txt>                      INTENTED FOR SMALL SCALE TESTS ONLY");  
    printf("\n                                                   Configure a directory to run locally. <InputPar.txt> name of file that");
    printf("\n                                                   contains input command template.");
    printf("\n                                                   Optional commands:  ");
    printf("\n                                                     --InputDir <InputDir>   Default value: $InputDir"); 
    printf("\n                                                       Note: the root files must be inside a $InputDir/<dir>/ ");
    printf("\n                                                       where dir contains user, data or mc as a substring");
    printf("\n                                                     --OutputDir <OutputDir> Default value: $OutputDir");
    printf("\n                                                     --CodeDir  <CodeDir>    Default value: $CodeDir");
    printf("\n                                                     --SetName <SetName>     Default value: $set ");
    printf("\n                                                     --NMaxData <Max Number of data files per job >     Default value: $maxdata ");
    printf("\n                                                     --NMaxMC <Max Number of MC files per job >     Default value: $maxmc ");
    printf("\n                                                     --NMaxEmbed <Max Number of Embedding files per job > Default value: $maxemb ");
    printf("\n                                                     --ROOTSYS <ROOTSYS> the current ROOTSYS variable if --BuildRoot is not defined");
    printf("\n                                                     --Proxy <path>;   a path to proxy recieved by running grid-init");

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
    if($ARGV[$l] eq "--CodeDir"){
	$l++;
	$CodeDir=$ARGV[$l];
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
    if($ARGV[$l] eq  "--NMaxEmbed" ){
	$l++;
	$maxemb=$ARGV[$l];
    }
    if($ARGV[$l] eq  "--BuildRoot"){
	$l++;
	$buildrootversion=$ARGV[$l];
	$buildRoot=1;
    } 
    if($ARGV[$l] eq  "--ROOTSYS"){
        $l++;
        $MYROOTSYS=$ARGV[$l];
	$hasroot=1;
    }
    if($ARGV[$l] eq  "--LongQueue"){
		$l++;
		$Queue="cream-pbs-cms";
    }
    if($ARGV[$l] eq  "--ARCH" ){
        $l++;
        $ARCH=$ARGV[$l];
    }
    if($ARGV[$l] eq  "--Proxy" ){
        $l++;
        $Proxy=$ARGV[$l];
    }
    if($ARGV[$l] eq  "--QsubQueue" ){
	$l++;
	if($ARGV[$l] eq  "short"){
	    $QsubQue="cms_local_short";
	}
	if($ARGV[$l] eq  "medium"){
	    $QsubQue="sbg_local_mdm";
	}
	if($ARGV[$l] eq  "long"){
	    $QsubQue="cms_local";
	}
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
    if($UserID neq "cherepan" or $UserID neq "wangjian" or $UserID neq "bjoshi" ){
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
    system(sprintf("echo \"cd $currentdir/$CMSSW_BASE/src/DsTau23Mu; git checkout $DsdevBranch; \" >> Install_DsTNtuple_$time"));


    system(sprintf("echo \"cd $currentdir/$CMSSW_BASE/src\" >> Install_DsTNtuple_$time"));
    system(sprintf("echo \"git clone git\@github.com:T3MuAnalysisTools/SkimProduction.git\" >> Install_DsTNtuple_$time"));
    system(sprintf("echo \"cd $currentdir/$CMSSW_BASE/src\" >> Install_DsTNtuple_$time"));
    system(sprintf("echo \"scram b -j 4\" >> Install_DsTNtuple_$time"));

    printf("\n\nInstructions:");
    printf("\nsource  Install_DsTNtuple_$time to complete installation, compilation might take some time...  ");
    printf("\n\n\nTo run test job do  'cmsRun analyze.py'  in  $CMSSW_BASE/src/DsTau23Mu/T3MNtuple/test/ \n\n");
}



if( $ARGV[0] eq "--Local" ){

    # User check
    if($UserID neq "cherepan" or $UserID neq "wangjian" or $UserID neq "bjoshi" ){
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
    system(sprintf("mkdir  $OutputDir/workdir$set/Code/i386_linux "));
    system(sprintf("cp -r $CodeDir/* $OutputDir/workdir$set/Code/ "));
    system(sprintf("mkdir $OutputDir/workdir$set/EPS "));
    system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/InputData "));
    
    # generate init script 
    system(sprintf("echo \"cd  $OutputDir/workdir$set/Code/CommonUtils/CMSSW_9_3_8/src\" >> $OutputDir/workdir$set/init "));
    system(sprintf("echo \"cmsenv \" >> $OutputDir/workdir$set/init "));

    # generate compile script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd  $OutputDir/workdir$set/Code/\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"./config \\\$\@ \"   >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"source $OutputDir/workdir$set/init \"   >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd  $OutputDir/workdir$set/Code/\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"gmake all \" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/compile")) ;
 
    # Generate Combine script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Combine")) ;
    system(sprintf("echo \"export workdir=\\\"$OutputDir/workdir$set/\\\"\" >> $OutputDir/workdir$set/Combine"));
    system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; source config \" >> $OutputDir/workdir$set/Combine")); 
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
    system(sprintf("echo \"queue  1 \" >> $OutputDir/workdir$set/Condor_Combine"));

    # Start Submit script
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Submit")) ; 
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
			system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; source config \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
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

