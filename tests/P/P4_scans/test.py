from utilities import path_of_file, evaluate_file

filenamelist = ["P2_dimer_scan_C2v.xyz","P4_dimer_scan_C3v.xyz","P4_dimer_scan_D3h.xyz"]

for filename in filenamelist: 
    evaluate_file(path_of_file(__file__)+"/"+filename)

properties = {"files": filenamelist}
