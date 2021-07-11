from utilities import path_of_file, evaluate_file

filenamelist = ["Hittorf_scan_ActaB1969_mbd.xyz",   "Hittorf_scan_Angew2020_mbd.xyz", "Hittorf_ActaB1969.vasp_STEPS.xyz", "Hittorf_Angew2020.vasp_STEPS.xyz"]

for filename in filenamelist: 
    evaluate_file(path_of_file(__file__)+"/"+filename)

properties = {"files": filenamelist}

 
