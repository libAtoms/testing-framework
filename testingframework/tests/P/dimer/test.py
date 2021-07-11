from utilities import path_of_file, evaluate_file

filenamelist = ["dimer_mbd_reference.xyz"]

for filename in filenamelist: 
    evaluate_file(path_of_file(__file__)+"/"+filename)

properties = {"files": filenamelist}
