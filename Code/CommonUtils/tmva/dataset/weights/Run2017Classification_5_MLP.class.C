// Class: ReadMLP
// Automatically generated by MethodBase::MakeClass
//

/* configuration options =====================================================

#GEN -*-*-*-*-*-*-*-*-*-*-*- general info -*-*-*-*-*-*-*-*-*-*-*-

Method         : MLP::MLP
TMVA Release   : 4.2.1         [262657]
ROOT Release   : 6.10/05       [395781]
Creator        : bjoshi
Date           : Wed Sep 25 21:21:23 2019
Host           : Linux cmsbuild09.cern.ch 2.6.32-696.1.1.el6.x86_64 #1 SMP Wed Apr 12 08:44:28 CEST 2017 x86_64 x86_64 x86_64 GNU/Linux
Dir            : /afs/cern.ch/work/b/bjoshi/Analysis/workdirTMVAAnalysisSep_24_2019/Code/CommonUtils/tmva
Training events: 16339
Analysis type  : [Classification]


#OPT -*-*-*-*-*-*-*-*-*-*-*-*- options -*-*-*-*-*-*-*-*-*-*-*-*-

# Set by User:
NCycles: "600" [Number of training cycles]
HiddenLayers: "N+5" [Specification of hidden layer architecture]
NeuronType: "tanh" [Neuron activation function type]
V: "False" [Verbose output (short form of "VerbosityLevel" below - overrides the latter one)]
VarTransform: "N" [List of variable transformations performed before training, e.g., "D_Background,P_Signal,G,N_AllClasses" for: "Decorrelation, PCA-transformation, Gaussianisation, Normalisation, each for the given class of events ('AllClasses' denotes all events of all classes, if no class indication is given, 'All' is assumed)"]
H: "True" [Print method-specific help message]
TestRate: "5" [Test for overtraining performed at each #th epochs]
UseRegulator: "False" [Use regulator to avoid over-training]
# Default:
RandomSeed: "1" [Random seed for initial synapse weights (0 means unique seed for each run; default value '1')]
EstimatorType: "CE" [MSE (Mean Square Estimator) for Gaussian Likelihood or CE(Cross-Entropy) for Bernoulli Likelihood]
NeuronInputType: "sum" [Neuron input function type]
VerbosityLevel: "Default" [Verbosity level]
CreateMVAPdfs: "False" [Create PDFs for classifier outputs (signal and background)]
IgnoreNegWeightsInTraining: "False" [Events with negative weights are ignored in the training (but are included for testing and performance evaluation)]
TrainingMethod: "BP" [Train with Back-Propagation (BP), BFGS Algorithm (BFGS), or Genetic Algorithm (GA - slower and worse)]
LearningRate: "2.000000e-02" [ANN learning rate parameter]
DecayRate: "1.000000e-02" [Decay rate for learning parameter]
EpochMonitoring: "False" [Provide epoch-wise monitoring plots according to TestRate (caution: causes big ROOT output file!)]
Sampling: "1.000000e+00" [Only 'Sampling' (randomly selected) events are trained each epoch]
SamplingEpoch: "1.000000e+00" [Sampling is used for the first 'SamplingEpoch' epochs, afterwards, all events are taken for training]
SamplingImportance: "1.000000e+00" [ The sampling weights of events in epochs which successful (worse estimator than before) are multiplied with SamplingImportance, else they are divided.]
SamplingTraining: "True" [The training sample is sampled]
SamplingTesting: "False" [The testing sample is sampled]
ResetStep: "50" [How often BFGS should reset history]
Tau: "3.000000e+00" [LineSearch "size step"]
BPMode: "sequential" [Back-propagation learning mode: sequential or batch]
BatchSize: "-1" [Batch size: number of events/batch, only set if in Batch Mode, -1 for BatchSize=number_of_events]
ConvergenceImprove: "1.000000e-30" [Minimum improvement which counts as improvement (<0 means automatic convergence check is turned off)]
ConvergenceTests: "-1" [Number of steps (without improvement) required for convergence (<0 means automatic convergence check is turned off)]
UpdateLimit: "10000" [Maximum times of regulator update]
CalculateErrors: "False" [Calculates inverse Hessian matrix at the end of the training to be able to calculate the uncertainties of an MVA value]
WeightRange: "1.000000e+00" [Take the events for the estimator calculations from small deviations from the desired value to large deviations only over the weight range]
##


#VAR -*-*-*-*-*-*-*-*-*-*-*-* variables *-*-*-*-*-*-*-*-*-*-*-*-

NVar 12
var_vertexKFChi2              var_vertexKFChi2              var_vertexKFChi2              Variable vertexKFChi2         units                             'F'    [0.0049333521165,99.0441207886]
var_svpvTauAngle              var_svpvTauAngle              var_svpvTauAngle              Variable svpvTauAngle         units                             'F'    [0.000118228665087,2.82946825027]
var_flightLenSig              var_flightLenSig              var_flightLenSig              Variable flightLenSig         units                             'F'    [0.349412500858,21.5535106659]
var_sumMuTrkKinkChi2          var_sumMuTrkKinkChi2          var_sumMuTrkKinkChi2          Variable sumMuTrkKinkChi2     units                             'F'    [5.5704293251,3133.14746094]
var_segCompMuMin              var_segCompMuMin              var_segCompMuMin              Variable segCompMuMin         units                             'F'    [0,1]
var_MinMIPLikelihood          var_MinMIPLikelihood          var_MinMIPLikelihood          Variable MinMIPLikelihood     units                             'F'    [1.1731815448e-06,0.995340228081]
var_maxdca                    var_maxdca                    var_maxdca                    Variable maxdca               units                             'F'    [0.000125885009766,5.06948852539]
var_MuMu_mindR                var_MuMu_mindR                var_MuMu_mindR                Variable MuMu_mindR           units                             'F'    [0.00179446791299,0.370722234249]
var_RelPt_Mu1Tau              var_RelPt_Mu1Tau              var_RelPt_Mu1Tau              Variable RelPt_Mu1Tau         units                             'F'    [0.0573893263936,0.855991661549]
var_Eta_au                    var_Eta_au                    var_Eta_au                    Variable Eta_au               units                             'F'    [-1.97041010857,1.96270954609]
var_MuMu_minKFChi2            var_MuMu_minKFChi2            var_MuMu_minKFChi2            Variable MuMu_minKFChi2       units                             'F'    [0,26.7062072754]
var_MuTau_maxdR               var_MuTau_maxdR               var_MuTau_maxdR               Variable MuTau_maxdR          units                             'F'    [0.00222029048018,0.379105210304]
NSpec 1
var_tauMass                   var_tauMass                   var_tauMass                   Variable tauMass              units                             'F'    [0.589990794659,3.07679438591]


============================================================================ */

#include <vector>
#include <cmath>
#include <string>
#include <iostream>

#ifndef IClassifierReader__def
#define IClassifierReader__def

class IClassifierReader {

 public:

   // constructor
   IClassifierReader() : fStatusIsClean( true ) {}
   virtual ~IClassifierReader() {}

   // return classifier response
   virtual double GetMvaValue( const std::vector<double>& inputValues ) const = 0;

   // returns classifier status
   bool IsStatusClean() const { return fStatusIsClean; }

 protected:

   bool fStatusIsClean;
};

#endif

class ReadMLP : public IClassifierReader {

 public:

   // constructor
   ReadMLP( std::vector<std::string>& theInputVars ) 
      : IClassifierReader(),
        fClassName( "ReadMLP" ),
        fNvars( 12 ),
        fIsNormalised( false )
   {      
      // the training input variables
      const char* inputVars[] = { "var_vertexKFChi2", "var_svpvTauAngle", "var_flightLenSig", "var_sumMuTrkKinkChi2", "var_segCompMuMin", "var_MinMIPLikelihood", "var_maxdca", "var_MuMu_mindR", "var_RelPt_Mu1Tau", "var_Eta_au", "var_MuMu_minKFChi2", "var_MuTau_maxdR" };

      // sanity checks
      if (theInputVars.size() <= 0) {
         std::cout << "Problem in class \"" << fClassName << "\": empty input vector" << std::endl;
         fStatusIsClean = false;
      }

      if (theInputVars.size() != fNvars) {
         std::cout << "Problem in class \"" << fClassName << "\": mismatch in number of input values: "
                   << theInputVars.size() << " != " << fNvars << std::endl;
         fStatusIsClean = false;
      }

      // validate input variables
      for (size_t ivar = 0; ivar < theInputVars.size(); ivar++) {
         if (theInputVars[ivar] != inputVars[ivar]) {
            std::cout << "Problem in class \"" << fClassName << "\": mismatch in input variable names" << std::endl
                      << " for variable [" << ivar << "]: " << theInputVars[ivar].c_str() << " != " << inputVars[ivar] << std::endl;
            fStatusIsClean = false;
         }
      }

      // initialize min and max vectors (for normalisation)
      fVmin[0] = -1;
      fVmax[0] = 1;
      fVmin[1] = -1;
      fVmax[1] = 1;
      fVmin[2] = -1;
      fVmax[2] = 1;
      fVmin[3] = -1;
      fVmax[3] = 1;
      fVmin[4] = -1;
      fVmax[4] = 1;
      fVmin[5] = -1;
      fVmax[5] = 1;
      fVmin[6] = -1;
      fVmax[6] = 1;
      fVmin[7] = -1;
      fVmax[7] = 1;
      fVmin[8] = -1;
      fVmax[8] = 0.99999988079071;
      fVmin[9] = -1;
      fVmax[9] = 1;
      fVmin[10] = -1;
      fVmax[10] = 1;
      fVmin[11] = -1;
      fVmax[11] = 1;

      // initialize input variable types
      fType[0] = 'F';
      fType[1] = 'F';
      fType[2] = 'F';
      fType[3] = 'F';
      fType[4] = 'F';
      fType[5] = 'F';
      fType[6] = 'F';
      fType[7] = 'F';
      fType[8] = 'F';
      fType[9] = 'F';
      fType[10] = 'F';
      fType[11] = 'F';

      // initialize constants
      Initialize();

      // initialize transformation
      InitTransform();
   }

   // destructor
   virtual ~ReadMLP() {
      Clear(); // method-specific
   }

   // the classifier response
   // "inputValues" is a vector of input values in the same order as the 
   // variables given to the constructor
   double GetMvaValue( const std::vector<double>& inputValues ) const;

 private:

   // method-specific destructor
   void Clear();

   // input variable transformation

   double fMin_1[3][12];
   double fMax_1[3][12];
   void InitTransform_1();
   void Transform_1( std::vector<double> & iv, int sigOrBgd ) const;
   void InitTransform();
   void Transform( std::vector<double> & iv, int sigOrBgd ) const;

   // common member variables
   const char* fClassName;

   const size_t fNvars;
   size_t GetNvar()           const { return fNvars; }
   char   GetType( int ivar ) const { return fType[ivar]; }

   // normalisation of input variables
   const bool fIsNormalised;
   bool IsNormalised() const { return fIsNormalised; }
   double fVmin[12];
   double fVmax[12];
   double NormVariable( double x, double xmin, double xmax ) const {
      // normalise to output range: [-1, 1]
      return 2*(x - xmin)/(xmax - xmin) - 1.0;
   }

   // type of input variable: 'F' or 'I'
   char   fType[12];

   // initialize internal variables
   void Initialize();
   double GetMvaValue__( const std::vector<double>& inputValues ) const;

   // private members (method specific)

   double ActivationFnc(double x) const;
   double OutputActivationFnc(double x) const;

   int fLayers;
   int fLayerSize[3];
   double fWeightMatrix0to1[18][13];   // weight matrix from layer 0 to 1
   double fWeightMatrix1to2[1][18];   // weight matrix from layer 1 to 2

   double * fWeights[3];
};

inline void ReadMLP::Initialize()
{
   // build network structure
   fLayers = 3;
   fLayerSize[0] = 13; fWeights[0] = new double[13]; 
   fLayerSize[1] = 18; fWeights[1] = new double[18]; 
   fLayerSize[2] = 1; fWeights[2] = new double[1]; 
   // weight matrix from layer 0 to 1
   fWeightMatrix0to1[0][0] = 0.159355433149732;
   fWeightMatrix0to1[1][0] = 2.91371941096358;
   fWeightMatrix0to1[2][0] = 1.42607169902996;
   fWeightMatrix0to1[3][0] = 1.70649912012474;
   fWeightMatrix0to1[4][0] = -1.23263378286829;
   fWeightMatrix0to1[5][0] = 1.6996876625311;
   fWeightMatrix0to1[6][0] = -0.212342870149715;
   fWeightMatrix0to1[7][0] = 1.59390772594193;
   fWeightMatrix0to1[8][0] = -0.111017394186082;
   fWeightMatrix0to1[9][0] = -0.464235833787001;
   fWeightMatrix0to1[10][0] = -1.55754299662508;
   fWeightMatrix0to1[11][0] = 1.08177380556396;
   fWeightMatrix0to1[12][0] = 0.427820546217828;
   fWeightMatrix0to1[13][0] = -0.00096156729232802;
   fWeightMatrix0to1[14][0] = 0.818720190612208;
   fWeightMatrix0to1[15][0] = -0.0826220732011509;
   fWeightMatrix0to1[16][0] = -0.538249229820396;
   fWeightMatrix0to1[0][1] = 6.54424606525752;
   fWeightMatrix0to1[1][1] = 13.2127927397758;
   fWeightMatrix0to1[2][1] = 1.75086384569624;
   fWeightMatrix0to1[3][1] = 0.0587546192076868;
   fWeightMatrix0to1[4][1] = 3.63343101978135;
   fWeightMatrix0to1[5][1] = -2.59338174324647;
   fWeightMatrix0to1[6][1] = -8.84516696789669;
   fWeightMatrix0to1[7][1] = -0.260012781861232;
   fWeightMatrix0to1[8][1] = 0.595281885013279;
   fWeightMatrix0to1[9][1] = -0.813601742468197;
   fWeightMatrix0to1[10][1] = -0.81556371784502;
   fWeightMatrix0to1[11][1] = -6.80168037570905;
   fWeightMatrix0to1[12][1] = -0.108478844259068;
   fWeightMatrix0to1[13][1] = 1.1837098986228;
   fWeightMatrix0to1[14][1] = 2.53415995651723;
   fWeightMatrix0to1[15][1] = 0.981483125663883;
   fWeightMatrix0to1[16][1] = -2.07702121912231;
   fWeightMatrix0to1[0][2] = 0.710241446461506;
   fWeightMatrix0to1[1][2] = 0.0966813058324193;
   fWeightMatrix0to1[2][2] = -0.541649430660666;
   fWeightMatrix0to1[3][2] = 2.96488457018781;
   fWeightMatrix0to1[4][2] = -1.5647071100111;
   fWeightMatrix0to1[5][2] = 3.96991087055345;
   fWeightMatrix0to1[6][2] = 1.5503188353626;
   fWeightMatrix0to1[7][2] = 1.19670827311657;
   fWeightMatrix0to1[8][2] = 2.86797648084016;
   fWeightMatrix0to1[9][2] = 1.52719907037006;
   fWeightMatrix0to1[10][2] = 0.152412445279986;
   fWeightMatrix0to1[11][2] = -3.00178675493331;
   fWeightMatrix0to1[12][2] = 0.722178340802797;
   fWeightMatrix0to1[13][2] = 0.565254305667548;
   fWeightMatrix0to1[14][2] = 0.706086111793715;
   fWeightMatrix0to1[15][2] = 0.339960051070184;
   fWeightMatrix0to1[16][2] = 0.000233787670683847;
   fWeightMatrix0to1[0][3] = 0.139329699384086;
   fWeightMatrix0to1[1][3] = 3.11401002998885;
   fWeightMatrix0to1[2][3] = 2.24170552744276;
   fWeightMatrix0to1[3][3] = -1.27469838622321;
   fWeightMatrix0to1[4][3] = 0.542952409530418;
   fWeightMatrix0to1[5][3] = -3.42138852649843;
   fWeightMatrix0to1[6][3] = -5.94731000136625;
   fWeightMatrix0to1[7][3] = 1.92141710857326;
   fWeightMatrix0to1[8][3] = 0.762033923174538;
   fWeightMatrix0to1[9][3] = -3.0481507583288;
   fWeightMatrix0to1[10][3] = 0.69196261232123;
   fWeightMatrix0to1[11][3] = -2.10244298286555;
   fWeightMatrix0to1[12][3] = 0.463695078144157;
   fWeightMatrix0to1[13][3] = 2.23009577351435;
   fWeightMatrix0to1[14][3] = 0.107080550849294;
   fWeightMatrix0to1[15][3] = 0.897728537771469;
   fWeightMatrix0to1[16][3] = -1.47136238099611;
   fWeightMatrix0to1[0][4] = -0.299414896868411;
   fWeightMatrix0to1[1][4] = -0.178302265429241;
   fWeightMatrix0to1[2][4] = -0.706422472941404;
   fWeightMatrix0to1[3][4] = 0.690469321435864;
   fWeightMatrix0to1[4][4] = -0.415107814294569;
   fWeightMatrix0to1[5][4] = 0.196305682711577;
   fWeightMatrix0to1[6][4] = 0.249099363280206;
   fWeightMatrix0to1[7][4] = -0.740516559608369;
   fWeightMatrix0to1[8][4] = -2.8212157747968;
   fWeightMatrix0to1[9][4] = -0.0261692929672152;
   fWeightMatrix0to1[10][4] = 2.20167280001367;
   fWeightMatrix0to1[11][4] = -0.308962596086428;
   fWeightMatrix0to1[12][4] = 1.80382515349357;
   fWeightMatrix0to1[13][4] = 0.302032937420943;
   fWeightMatrix0to1[14][4] = -0.177502486159687;
   fWeightMatrix0to1[15][4] = -0.319459735666828;
   fWeightMatrix0to1[16][4] = -0.52992983479266;
   fWeightMatrix0to1[0][5] = 0.0209597763827185;
   fWeightMatrix0to1[1][5] = -0.0011843502531366;
   fWeightMatrix0to1[2][5] = -1.01068893885846;
   fWeightMatrix0to1[3][5] = 0.385171903022295;
   fWeightMatrix0to1[4][5] = -0.268399549054564;
   fWeightMatrix0to1[5][5] = 0.0671332545332679;
   fWeightMatrix0to1[6][5] = 0.104039264760952;
   fWeightMatrix0to1[7][5] = 2.23919676525042;
   fWeightMatrix0to1[8][5] = 0.17133570425233;
   fWeightMatrix0to1[9][5] = -0.27963221815517;
   fWeightMatrix0to1[10][5] = 0.384976222872958;
   fWeightMatrix0to1[11][5] = -0.0971814807024282;
   fWeightMatrix0to1[12][5] = -0.190723319083123;
   fWeightMatrix0to1[13][5] = -1.10825150496351;
   fWeightMatrix0to1[14][5] = -0.389910293156673;
   fWeightMatrix0to1[15][5] = -1.70545669603829;
   fWeightMatrix0to1[16][5] = -0.166677494845533;
   fWeightMatrix0to1[0][6] = 0.0463105195102044;
   fWeightMatrix0to1[1][6] = -0.326684598338379;
   fWeightMatrix0to1[2][6] = -0.774595468172434;
   fWeightMatrix0to1[3][6] = 0.888705014741333;
   fWeightMatrix0to1[4][6] = 1.02127128118598;
   fWeightMatrix0to1[5][6] = -0.932558090596249;
   fWeightMatrix0to1[6][6] = -1.29182350935084;
   fWeightMatrix0to1[7][6] = 0.158364602461881;
   fWeightMatrix0to1[8][6] = -2.19443362705697;
   fWeightMatrix0to1[9][6] = -0.217790564938043;
   fWeightMatrix0to1[10][6] = 0.740954765704774;
   fWeightMatrix0to1[11][6] = 1.7666377192004;
   fWeightMatrix0to1[12][6] = -1.69655232941194;
   fWeightMatrix0to1[13][6] = 2.37276523745117;
   fWeightMatrix0to1[14][6] = 0.508943411254993;
   fWeightMatrix0to1[15][6] = -0.773941120370965;
   fWeightMatrix0to1[16][6] = 1.17641221001677;
   fWeightMatrix0to1[0][7] = -0.541904346494964;
   fWeightMatrix0to1[1][7] = 0.151985902736133;
   fWeightMatrix0to1[2][7] = 0.373715502089738;
   fWeightMatrix0to1[3][7] = -1.11347230707014;
   fWeightMatrix0to1[4][7] = 2.36511748913919;
   fWeightMatrix0to1[5][7] = 2.4906282513436;
   fWeightMatrix0to1[6][7] = 0.260840431789981;
   fWeightMatrix0to1[7][7] = -0.215189411263963;
   fWeightMatrix0to1[8][7] = -0.0734390488370109;
   fWeightMatrix0to1[9][7] = -3.47808457108168;
   fWeightMatrix0to1[10][7] = 1.72305214617039;
   fWeightMatrix0to1[11][7] = 0.177909489687927;
   fWeightMatrix0to1[12][7] = -3.41887539532853;
   fWeightMatrix0to1[13][7] = 0.713492786239387;
   fWeightMatrix0to1[14][7] = 0.325403158297487;
   fWeightMatrix0to1[15][7] = -1.38995541884365;
   fWeightMatrix0to1[16][7] = 0.00285162760227615;
   fWeightMatrix0to1[0][8] = 3.67971831119039;
   fWeightMatrix0to1[1][8] = -0.216283625565823;
   fWeightMatrix0to1[2][8] = 0.81770031199268;
   fWeightMatrix0to1[3][8] = 1.72288414785784;
   fWeightMatrix0to1[4][8] = 1.7210763822435;
   fWeightMatrix0to1[5][8] = 0.230079071929961;
   fWeightMatrix0to1[6][8] = -0.00407049015801192;
   fWeightMatrix0to1[7][8] = -0.14578939415734;
   fWeightMatrix0to1[8][8] = -1.07451003175022;
   fWeightMatrix0to1[9][8] = 1.18410711906612;
   fWeightMatrix0to1[10][8] = 1.31176844961813;
   fWeightMatrix0to1[11][8] = -0.132703904179386;
   fWeightMatrix0to1[12][8] = -2.12418439687416;
   fWeightMatrix0to1[13][8] = 0.123497679780319;
   fWeightMatrix0to1[14][8] = -0.207496245106073;
   fWeightMatrix0to1[15][8] = 0.188352456367792;
   fWeightMatrix0to1[16][8] = -0.770963122438986;
   fWeightMatrix0to1[0][9] = -0.00679886772695323;
   fWeightMatrix0to1[1][9] = 0.0119254630586845;
   fWeightMatrix0to1[2][9] = 0.641066537059516;
   fWeightMatrix0to1[3][9] = 2.49454600147258;
   fWeightMatrix0to1[4][9] = -0.137985447997933;
   fWeightMatrix0to1[5][9] = 0.445270669766016;
   fWeightMatrix0to1[6][9] = 0.0550963505063713;
   fWeightMatrix0to1[7][9] = 1.10420312056378;
   fWeightMatrix0to1[8][9] = -0.898947488178856;
   fWeightMatrix0to1[9][9] = 0.0435564819966079;
   fWeightMatrix0to1[10][9] = -0.497030274763558;
   fWeightMatrix0to1[11][9] = 0.201477110212029;
   fWeightMatrix0to1[12][9] = -0.626926172757365;
   fWeightMatrix0to1[13][9] = -0.274140914996036;
   fWeightMatrix0to1[14][9] = 3.7419831959691;
   fWeightMatrix0to1[15][9] = -1.93644952803654;
   fWeightMatrix0to1[16][9] = -0.0759863389812192;
   fWeightMatrix0to1[0][10] = -0.328358354350331;
   fWeightMatrix0to1[1][10] = -1.20210345122414;
   fWeightMatrix0to1[2][10] = -0.0886244451095001;
   fWeightMatrix0to1[3][10] = 1.79903814293897;
   fWeightMatrix0to1[4][10] = 2.14218578316726;
   fWeightMatrix0to1[5][10] = -3.90767830160248;
   fWeightMatrix0to1[6][10] = -0.10923981291133;
   fWeightMatrix0to1[7][10] = 0.68788692136592;
   fWeightMatrix0to1[8][10] = 0.288890803183832;
   fWeightMatrix0to1[9][10] = -1.13994799222339;
   fWeightMatrix0to1[10][10] = -2.09951925502493;
   fWeightMatrix0to1[11][10] = 1.53654577511084;
   fWeightMatrix0to1[12][10] = 1.68606865532062;
   fWeightMatrix0to1[13][10] = -1.09268092517474;
   fWeightMatrix0to1[14][10] = -0.986441338845396;
   fWeightMatrix0to1[15][10] = -0.293545866536356;
   fWeightMatrix0to1[16][10] = -0.100650655848625;
   fWeightMatrix0to1[0][11] = 3.18279395868568;
   fWeightMatrix0to1[1][11] = 0.114331603820721;
   fWeightMatrix0to1[2][11] = -0.688005682260692;
   fWeightMatrix0to1[3][11] = -0.860882592693609;
   fWeightMatrix0to1[4][11] = -1.22927365962803;
   fWeightMatrix0to1[5][11] = 0.216508450249584;
   fWeightMatrix0to1[6][11] = 0.422388419134174;
   fWeightMatrix0to1[7][11] = 0.140106542968384;
   fWeightMatrix0to1[8][11] = -0.686994804452515;
   fWeightMatrix0to1[9][11] = -0.439035343559631;
   fWeightMatrix0to1[10][11] = 1.5737587459553;
   fWeightMatrix0to1[11][11] = -0.335428351122128;
   fWeightMatrix0to1[12][11] = -1.21588013456318;
   fWeightMatrix0to1[13][11] = -0.487957617887392;
   fWeightMatrix0to1[14][11] = -0.10567385441575;
   fWeightMatrix0to1[15][11] = -1.43486743149671;
   fWeightMatrix0to1[16][11] = -4.13593302593901;
   fWeightMatrix0to1[0][12] = 6.52803441067754;
   fWeightMatrix0to1[1][12] = 17.9755787377274;
   fWeightMatrix0to1[2][12] = -1.01861385307224;
   fWeightMatrix0to1[3][12] = -0.832353363122368;
   fWeightMatrix0to1[4][12] = 4.40546933011749;
   fWeightMatrix0to1[5][12] = -3.47299255887641;
   fWeightMatrix0to1[6][12] = -15.1452643585912;
   fWeightMatrix0to1[7][12] = -0.987775943675691;
   fWeightMatrix0to1[8][12] = -0.786172071038203;
   fWeightMatrix0to1[9][12] = -4.50723253694749;
   fWeightMatrix0to1[10][12] = 0.826467575321477;
   fWeightMatrix0to1[11][12] = -6.15202387125304;
   fWeightMatrix0to1[12][12] = 1.42387794523131;
   fWeightMatrix0to1[13][12] = -0.465442177262981;
   fWeightMatrix0to1[14][12] = 5.05704492076201;
   fWeightMatrix0to1[15][12] = -0.00296886204883379;
   fWeightMatrix0to1[16][12] = -6.20568982609174;
   // weight matrix from layer 1 to 2
   fWeightMatrix1to2[0][0] = -1.79487320102748;
   fWeightMatrix1to2[0][1] = -3.8528811150003;
   fWeightMatrix1to2[0][2] = 1.96697060737841;
   fWeightMatrix1to2[0][3] = 1.02728156597104;
   fWeightMatrix1to2[0][4] = -1.59174302292806;
   fWeightMatrix1to2[0][5] = 1.18262708820737;
   fWeightMatrix1to2[0][6] = 3.79616549052602;
   fWeightMatrix1to2[0][7] = -2.44387578557635;
   fWeightMatrix1to2[0][8] = -0.376367941695034;
   fWeightMatrix1to2[0][9] = 1.17067208547506;
   fWeightMatrix1to2[0][10] = 0.531368807412176;
   fWeightMatrix1to2[0][11] = 1.6325154956001;
   fWeightMatrix1to2[0][12] = -0.980716047560631;
   fWeightMatrix1to2[0][13] = 0.277859154656301;
   fWeightMatrix1to2[0][14] = -0.866818405766635;
   fWeightMatrix1to2[0][15] = -0.45898623294597;
   fWeightMatrix1to2[0][16] = 1.26355988055152;
   fWeightMatrix1to2[0][17] = -0.420924925700712;
}

inline double ReadMLP::GetMvaValue__( const std::vector<double>& inputValues ) const
{
   if (inputValues.size() != (unsigned int)fLayerSize[0]-1) {
      std::cout << "Input vector needs to be of size " << fLayerSize[0]-1 << std::endl;
      return 0;
   }

   for (int l=0; l<fLayers; l++)
      for (int i=0; i<fLayerSize[l]; i++) fWeights[l][i]=0;

   for (int l=0; l<fLayers-1; l++)
      fWeights[l][fLayerSize[l]-1]=1;

   for (int i=0; i<fLayerSize[0]-1; i++)
      fWeights[0][i]=inputValues[i];

   // layer 0 to 1
   for (int o=0; o<fLayerSize[1]-1; o++) {
      for (int i=0; i<fLayerSize[0]; i++) {
         double inputVal = fWeightMatrix0to1[o][i] * fWeights[0][i];
         fWeights[1][o] += inputVal;
      }
      fWeights[1][o] = ActivationFnc(fWeights[1][o]);
   }
   // layer 1 to 2
   for (int o=0; o<fLayerSize[2]; o++) {
      for (int i=0; i<fLayerSize[1]; i++) {
         double inputVal = fWeightMatrix1to2[o][i] * fWeights[1][i];
         fWeights[2][o] += inputVal;
      }
      fWeights[2][o] = OutputActivationFnc(fWeights[2][o]);
   }

   return fWeights[2][0];
}

double ReadMLP::ActivationFnc(double x) const {
   // hyperbolic tan
   return tanh(x);
}
double ReadMLP::OutputActivationFnc(double x) const {
   // sigmoid
   return 1.0/(1.0+exp(-x));
}
   
// Clean up
inline void ReadMLP::Clear() 
{
   // clean up the arrays
   for (int lIdx = 0; lIdx < 3; lIdx++) {
      delete[] fWeights[lIdx];
   }
}
   inline double ReadMLP::GetMvaValue( const std::vector<double>& inputValues ) const
   {
      // classifier response value
      double retval = 0;

      // classifier response, sanity check first
      if (!IsStatusClean()) {
         std::cout << "Problem in class \"" << fClassName << "\": cannot return classifier response"
                   << " because status is dirty" << std::endl;
         retval = 0;
      }
      else {
         if (IsNormalised()) {
            // normalise variables
            std::vector<double> iV;
            iV.reserve(inputValues.size());
            int ivar = 0;
            for (std::vector<double>::const_iterator varIt = inputValues.begin();
                 varIt != inputValues.end(); varIt++, ivar++) {
               iV.push_back(NormVariable( *varIt, fVmin[ivar], fVmax[ivar] ));
            }
            Transform( iV, -1 );
            retval = GetMvaValue__( iV );
         }
         else {
            std::vector<double> iV;
            int ivar = 0;
            for (std::vector<double>::const_iterator varIt = inputValues.begin();
                 varIt != inputValues.end(); varIt++, ivar++) {
               iV.push_back(*varIt);
            }
            Transform( iV, -1 );
            retval = GetMvaValue__( iV );
         }
      }

      return retval;
   }

//_______________________________________________________________________
inline void ReadMLP::InitTransform_1()
{
   // Normalization transformation, initialisation
   fMin_1[0][0] = 0.0049333521165;
   fMax_1[0][0] = 97.8059234619;
   fMin_1[1][0] = 0.0251949690282;
   fMax_1[1][0] = 99.0441207886;
   fMin_1[2][0] = 0.0049333521165;
   fMax_1[2][0] = 99.0441207886;
   fMin_1[0][1] = 0.000118228665087;
   fMax_1[0][1] = 0.384745001793;
   fMin_1[1][1] = 0.000358640885679;
   fMax_1[1][1] = 2.82946825027;
   fMin_1[2][1] = 0.000118228665087;
   fMax_1[2][1] = 2.82946825027;
   fMin_1[0][2] = 0.746487259865;
   fMax_1[0][2] = 16.0348377228;
   fMin_1[1][2] = 0.349412500858;
   fMax_1[1][2] = 21.5535106659;
   fMin_1[2][2] = 0.349412500858;
   fMax_1[2][2] = 21.5535106659;
   fMin_1[0][3] = 5.65686798096;
   fMax_1[0][3] = 1015.1350708;
   fMin_1[1][3] = 5.5704293251;
   fMax_1[1][3] = 3133.14746094;
   fMin_1[2][3] = 5.5704293251;
   fMax_1[2][3] = 3133.14746094;
   fMin_1[0][4] = 0;
   fMax_1[0][4] = 1;
   fMin_1[1][4] = 0;
   fMax_1[1][4] = 1;
   fMin_1[2][4] = 0;
   fMax_1[2][4] = 1;
   fMin_1[0][5] = 5.43361875316e-06;
   fMax_1[0][5] = 0.995340228081;
   fMin_1[1][5] = 1.1731815448e-06;
   fMax_1[1][5] = 0.994231462479;
   fMin_1[2][5] = 1.1731815448e-06;
   fMax_1[2][5] = 0.995340228081;
   fMin_1[0][6] = 0.000125885009766;
   fMax_1[0][6] = 5.06948852539;
   fMin_1[1][6] = 0.000748634338379;
   fMax_1[1][6] = 3.28829360008;
   fMin_1[2][6] = 0.000125885009766;
   fMax_1[2][6] = 5.06948852539;
   fMin_1[0][7] = 0.00179446791299;
   fMax_1[0][7] = 0.342327773571;
   fMin_1[1][7] = 0.00206668907776;
   fMax_1[1][7] = 0.370722234249;
   fMin_1[2][7] = 0.00179446791299;
   fMax_1[2][7] = 0.370722234249;
   fMin_1[0][8] = 0.0573893263936;
   fMax_1[0][8] = 0.805370807648;
   fMin_1[1][8] = 0.0696834176779;
   fMax_1[1][8] = 0.855991661549;
   fMin_1[2][8] = 0.0573893263936;
   fMax_1[2][8] = 0.855991661549;
   fMin_1[0][9] = -1.93126273155;
   fMax_1[0][9] = 1.94589090347;
   fMin_1[1][9] = -1.97041010857;
   fMax_1[1][9] = 1.96270954609;
   fMin_1[2][9] = -1.97041010857;
   fMax_1[2][9] = 1.96270954609;
   fMin_1[0][10] = 0;
   fMax_1[0][10] = 21.6617603302;
   fMin_1[1][10] = 0;
   fMax_1[1][10] = 26.7062072754;
   fMin_1[2][10] = 0;
   fMax_1[2][10] = 26.7062072754;
   fMin_1[0][11] = 0.00227414560504;
   fMax_1[0][11] = 0.321451663971;
   fMin_1[1][11] = 0.00222029048018;
   fMax_1[1][11] = 0.379105210304;
   fMin_1[2][11] = 0.00222029048018;
   fMax_1[2][11] = 0.379105210304;
}

//_______________________________________________________________________
inline void ReadMLP::Transform_1( std::vector<double>& iv, int cls) const
{
   // Normalization transformation
   if (cls < 0 || cls > 2) {
   if (2 > 1 ) cls = 2;
      else cls = 2;
   }
   const int nVar = 12;

   // get indices of used variables

   // define the indices of the variables which are transformed by this transformation
   static std::vector<int> indicesGet;
   static std::vector<int> indicesPut;

   if ( indicesGet.empty() ) { 
      indicesGet.reserve(fNvars);
      indicesGet.push_back( 0);
      indicesGet.push_back( 1);
      indicesGet.push_back( 2);
      indicesGet.push_back( 3);
      indicesGet.push_back( 4);
      indicesGet.push_back( 5);
      indicesGet.push_back( 6);
      indicesGet.push_back( 7);
      indicesGet.push_back( 8);
      indicesGet.push_back( 9);
      indicesGet.push_back( 10);
      indicesGet.push_back( 11);
   } 
   if ( indicesPut.empty() ) { 
      indicesPut.reserve(fNvars);
      indicesPut.push_back( 0);
      indicesPut.push_back( 1);
      indicesPut.push_back( 2);
      indicesPut.push_back( 3);
      indicesPut.push_back( 4);
      indicesPut.push_back( 5);
      indicesPut.push_back( 6);
      indicesPut.push_back( 7);
      indicesPut.push_back( 8);
      indicesPut.push_back( 9);
      indicesPut.push_back( 10);
      indicesPut.push_back( 11);
   } 

   static std::vector<double> dv;
   dv.resize(nVar);
   for (int ivar=0; ivar<nVar; ivar++) dv[ivar] = iv[indicesGet.at(ivar)];
   for (int ivar=0;ivar<12;ivar++) {
      double offset = fMin_1[cls][ivar];
      double scale  = 1.0/(fMax_1[cls][ivar]-fMin_1[cls][ivar]);
      iv[indicesPut.at(ivar)] = (dv[ivar]-offset)*scale * 2 - 1;
   }
}

//_______________________________________________________________________
inline void ReadMLP::InitTransform()
{
   InitTransform_1();
}

//_______________________________________________________________________
inline void ReadMLP::Transform( std::vector<double>& iv, int sigOrBgd ) const
{
   Transform_1( iv, sigOrBgd );
}