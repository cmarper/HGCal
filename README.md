# HGCal
Development of a L1 tau algorithm for CMS HGCal, C. Martin Perez (LLR), 2019

* *Step 1*: Postprocessing of raw ntuples

Instructions (based on Z->tautau samples):
1. Run tree_skimmer.C to remove gentaus out of the HGCal acceptance.
2. Run tree_matcher.C to match 3D-clusters to gentaus.
3. Run tree_sorter.C to sort 3D-clusters by pT and BDT score, as well as to get the input variables for the BDTs


* *Step 2*: Decay mode multiclassifier training/performance 

Instructions:
1. Select and dump input variables with tree_sorter.C.
2. Run training and performance evaluation with HGCal_BDT/multiclass_DM/training.py


* *Step 3*: Calibration regressor training/performance 

Instructions:
1. Select and dump input variables with tree_sorter.C.
2. Run training and performance evaluation with HGCal_BDT/regression_calib/training.py
