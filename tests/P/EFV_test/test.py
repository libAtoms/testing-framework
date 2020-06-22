from utilities import path_of_file, evaluate_file

filenamelist = ["P_test_set.xyz", "Hittorf_ActaB1969_isolated_layer.xyz", "Hittorf_Angew2020_isolated_layer.xyz", "phosphorene_ribbon_armchair.xyz", "phosphorene_ribbon_zigzag.xyz", "phosphorene_ribbon_single_layer.xyz"]

for filename in filenamelist: 
    evaluate_file(path_of_file(__file__)+"/"+filename)

properties = {"files": filenamelist}

 
