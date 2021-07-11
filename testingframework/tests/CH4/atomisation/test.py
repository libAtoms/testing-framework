from utilities import path_of_file, evaluate_file

filenamelist = ["atomisation_1.xyz", "atomisation_2.xyz", "atomisation_3.xyz", "atomisation_4.xyz"]

for filename in filenamelist: 
    evaluate_file(path_of_file(__file__)+"/"+filename)

properties = {"files": filenamelist}

 
