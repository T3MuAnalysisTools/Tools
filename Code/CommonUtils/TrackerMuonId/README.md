Run the MuonPionTree class in joshi to create an Ntuple containing variables related to there muons and those relevant to TrackerMuonId BDT on DsToTau(signal muons) and DsToPhiPi(background muons).
The class can be modified to include other quantities, but at the moment most of the necessary quantities are included. 
The phase-space of the background has to be reweighed to match the phase-space of the muons from the signal. To do so, follow the example below.

`cd plots\./plot_weights.py -i <input_ntuple> -o <path_to_weight_file>`


The weights are stored in a 2D histogram with pt range [0,15] (nbins = 30) and eta range [-2.5,2.5] (nbins = 25). The weights can then be added to the MuonPionTree using the following script.

`cd ../utils\root -q -b -l makeWeighedTree.C(<input_ntuple>, <path_to_weight_file>)`

For TMVA training run the following commands.

`cd ../TMVA\./makeClassifier.py -i <input_file> -o <output_file> -c <classifier_name> -v <variable_set>`
