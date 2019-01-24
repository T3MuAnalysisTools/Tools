#! /usr/bin/perl
use Cwd;
use POSIX;
use POSIX qw(strftime);

#############################################
$numArgs = $#ARGV +1;
$ARGV[$argnum];

$UserID= POSIX::cuserid();
$UserIDCern=$UserID;
$UserDir="";

$PWD=getcwd;



$CMSSWRel="9_4_4";
$ARCH="slc6_amd64_gcc530";
$DsdevBranch = "master";
#export SCRAM_ARCH= slc6_amd64_gcc630
#source /cvmfs/cms.cern.ch/cmsset_default.sh

#printf("\n ---> Your user ID is:   $UserID \n");
if($ARGV[0] eq "--help" || $ARGV[0] eq ""){
    printf("\n========================================================================================");
    printf("\nThis code requires input options. The syntax is:./todo.pl [--OPTION1]  [--OPTION2] ... ");
    printf("\nPlease choose from the following options:\n");
    printf("\n./todo.pl --help                                   Prints this message");
    printf("\n./todo.pl --DsTauTo3MNtuple <dir>                  Clone and compile DsToTau ntuple. Example: ./todo.pl --DsTauTo3MNtuple workdir  ");
    printf("\n                                                                        --Branch <branch> developing branch; Default: master ");
    printf("\n./todo.pl                                           --MuonPogNtuple <dir> MuonPogNtuple  ");
    printf("\n                                                    --ARCH  <SCRAM_ARCH>   Setup SCRAM_ARCH; Default: slc6_amd64_gcc630 ");
    printf("\n                                                    --CMSSWRel <CMSSW_X_Y_Z>  Configure CMSSW_X_Y_Z; Default: CMSSW_9_4_4  \n\n");

    exit(0);  
}


for($l=2;$l<$numArgs; $l++){
    if($ARGV[$l] eq "--CMSSWRel"){
        $l++;
        $CMSSWRel=$ARGV[$l];
    }

    if($ARGV[$l] eq  "--ARCH" ){
        $l++;
        $ARCH=$ARGV[$l];
    }
}
$time= strftime("%h_%d_%Y",localtime);

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

    system(sprintf("rm Install_MuPoGNtuple_$time"));
    system(sprintf("echo \"export SCRAM_ARCH=\\\"$ARCH\\\"\" >> Install_MuPoGNtuple_$time"));
    system(sprintf("echo \"source /cvmfs/cms.cern.ch/cmsset_default.sh\" >> Install_MuPoGNtuple_$time"));

    system(sprintf("echo \"mkdir $basedir\" >>  Install_MuPoGNtuple_$time")); 
    system(sprintf("echo \"cd $basedir\" >>  Install_MuPoGNtuple_$time")); 
    system(sprintf("echo \"cmsrel CMSSW_$CMSSWRel\" >>  Install_MuPoGNtuple_$time")); 
    system(sprintf("echo \"cd CMSSW_$CMSSWRel/src\" >> Install_MuPoGNtuple_$time")); 
    system(sprintf("echo \"cmsenv\" >> Install_MuPoGNtuple_$time")); 
    $CMSPATH="/CMSSW_$CMSSWRel";
    $CMSSW_BASE="$basedir$CMSPATH";
    system(sprintf("echo \"git clone git\@github.com:T3MuAnalysisTools/MuonPOGtreeProducer.git\" >> Install_MuPoGNtuple_$time"));
    system(sprintf("echo \"cd MuonPOGtreeProducer; git checkout 9X_tau3mu;\" >> Install_MuPoGNtuple_$time"));
    system(sprintf("echo \"cd $currentdir/$CMSSW_BASE/src\" >> Install_MuPoGNtuple_$time")); 
    system(sprintf("echo \"scram b -j 4\" >> Install_MuPoGNtuple_$time")); 

    printf("\n\nInstructions:");
    printf("\nsource  Install_MuPoGNtuple_$time to complete installation, compilation might take some time...  ");
    printf("\n\n\nTo run test job do  'cmsRun muonPogNtuples_cfg.py'  in  $CMSSW_BASE/src/MuonPOGtreeProducer/Tools/test/ \n\n");
}

if( $ARGV[0] eq "--DsTauTo3MNtuple"){
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

    system(sprintf("rm Install_DsTNtuple_$time"));
    system(sprintf("echo \"export SCRAM_ARCH=\\\"$ARCH\\\"\" >> Install_DsTNtuple_$time"));
    system(sprintf("echo \"source /cvmfs/cms.cern.ch/cmsset_default.sh\" >> Install_DsTNtuple_$time"));

    system(sprintf("echo \"mkdir $basedir\" >>  Install_DsTNtuple_$time"));
    system(sprintf("echo \"cd $basedir\" >>  Install_DsTNtuple_$time"));
    system(sprintf("echo \"cmsrel CMSSW_$CMSSWRel\" >>  Install_DsTNtuple_$time"));
    system(sprintf("echo \"cd CMSSW_$CMSSWRel/src\" >> Install_DsTNtuple_$time"));
    system(sprintf("echo \"cmsenv\" >> Install_DsTNtuple_$time"));
    $CMSPATH="/CMSSW_$CMSSWRel";
    $CMSSW_BASE="$basedir$CMSPATH";
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
