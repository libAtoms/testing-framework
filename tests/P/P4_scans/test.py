from utilities import path_of_file, evaluate_file

filenamelist = ["P2_dimer_scan_C2v_mbd.xyz","P4_dimer_scan_C3v_mbd.xyz","P4_dimer_scan_D3h_mbd.xyz"]

for filename in filenamelist: 
    evaluate_file(path_of_file(__file__)+"/"+filename)

properties = {"files": filenamelist}
