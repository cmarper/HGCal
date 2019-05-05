# HGCal
Development of a L1 tau algorithm for CMS HGCal
C. Martin Perez, 2019

Instructions (based on Z->tautau samples):
1. Run tree_skimmer.C to remove gentaus out of the HGCal acceptance.
2. Run tree_matcher.C to match 3D-clusters to gentaus.
3. Run tree_sorter.C to sort 3D-clusters by pT and BDT score.
