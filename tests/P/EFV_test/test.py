from utilities import path_of_file, evaluate_file

filenamelist = ["P_test_set.xyz"]

for filename in filenamelist: 
    evaluate_file(path_of_file(__file__)+"/"+filename)

properties = {"files": filenamelist}

 
