#TMVA Analysis

The events from background and the signal sample are categorized based on three muon mass resolution.  For training the BDT the background samples are further classified into two categories based on the sidebands region they are inF (right or left).  For training half the events are used and for testing the rest half is used. As an example, NewTMVAVars was given as a 
The following features were used to train the BDT:
var_vertexKFChi2 (chi sq of the fit of the secondary vertex)
var_svpvTauAngle (The angle between PV-SV vector and the tau vector)
var_flightLenSig (Flight length significance of the tau candidate)
var_sumMuTrkKinkChi2 (sum of chi sq of the kink of all three muons)
var_segCompMuMin (Minimum of the segment compatibility of the three muons)
var_MinMIPLikelihood (Minimum of the calorimeter compatibility of the three muons)

##How to run TMVA

Write a class with all the features and their corresponding histograms. Run the analysis from framework.
Go to macros and run add_files.py to combine all the TMVA trees from all the MC sets into one file and TMVA trees from all the data sets into another.
Example:
`python add_files.py —inputFile NewTMVAVarsInput —instanceName NewTMVAVars`
The two files will be stored as NewTMVAVars_MC_merged.root and NewTMVAVars_DATA_merged.root and the TMVA Tree with signal and background will be stored as NewTMVAVars_Input.root

Then go to CommonUtils/tmva and run makeTMVAClassifier.py
Example:
`Python makeTMVAClassifier.py  --inputFile NewTMVAVars_Input.root —outputFile NewTMVAVars_Output.root —category 1 —classifier NewTMVAClassification`
This will create a new classifier named NewTMVAClassification_1.cxx. The number of events from category 1 of signal and background will be split in half for training and training purpose. Use the bask script make classifiers for all the categories.

Then use runClassification.py to run the TMVA analysis by specifying the methods in the options.
Example:
`Python runClassification.py —classifier NewTMVAVarsClassification —methods BDT,DNN,Fisher`

The output, the xml file containing weights and the BDT node class corresponding to each method used in training, will be store in dataset/weights. 



