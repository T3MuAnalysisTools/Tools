#python makeTMVAClassifier.py --inputFile ../../macros/TMVATrees_Run2017_MC.root --outputFile TMVATrees_Run2017_MC_Output_1.root --category 1 --classifier Run2017Classification_1
python makeTMVAClassifier.py --inputFile ../../macros/TMVATrees_Run2017_MC.root --outputFile TMVATrees_Run2017_MC_Output_2.root --category 2 --classifier Run2017Classification_2
python makeTMVAClassifier.py --inputFile ../../macros/TMVATrees_Run2017_MC.root --outputFile TMVATrees_Run2017_MC_Output_3.root --category 3 --classifier Run2017Classification_3
python makeTMVAClassifier.py --inputFile ../../macros/TMVATrees_Run2017_MC.root --outputFile TMVATrees_Run2017_MC_Output_4.root --category 4 --classifier Run2017Classification_4
python makeTMVAClassifier.py --inputFile ../../macros/TMVATrees_Run2017_MC.root --outputFile TMVATrees_Run2017_MC_Output_5.root --category 5 --classifier Run2017Classification_5
python makeTMVAClassifier.py --inputFile ../../macros/TMVATrees_Run2017_MC.root --outputFile TMVATrees_Run2017_MC_Output_6.root --category 6 --classifier Run2017Classification_6
python makeTMVAClassifier.py --inputFile ../../macros/TMVATrees_Run2017_MC.root --outputFile TMVATrees_Run2017_MC_Output_7.root --category 7 --classifier Run2017Classification_7
python makeTMVAClassifier.py --inputFile ../../macros/TMVATrees_Run2017_MC.root --outputFile TMVATrees_Run2017_MC_Output_8.root --category 8 --classifier Run2017Classification_8
python makeTMVAClassifier.py --inputFile ../../macros/TMVATrees_Run2017_MC.root --outputFile TMVATrees_Run2017_MC_Output_9.root --category 9 --classifier Run2017Classification_9
#root -q -b -l ./Run2017Classification_1.cxx\(\"BDT,MLP\"\)
root -q -b -l ./Run2017Classification_2.cxx\(\"BDT,MLP\"\)
root -q -b -l ./Run2017Classification_3.cxx\(\"BDT,MLP\"\)
root -q -b -l ./Run2017Classification_4.cxx\(\"BDT,MLP\"\)
root -q -b -l ./Run2017Classification_5.cxx\(\"BDT,MLP\"\)
root -q -b -l ./Run2017Classification_6.cxx\(\"BDT,MLP\"\)
root -q -b -l ./Run2017Classification_7.cxx\(\"BDT,MLP\"\)
root -q -b -l ./Run2017Classification_8.cxx\(\"BDT,MLP\"\)
root -q -b -l ./Run2017Classification_9.cxx\(\"BDT,MLP\"\)
