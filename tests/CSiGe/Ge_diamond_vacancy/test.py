import os.path, vacancy

properties = vacancy.do_all_vacancies(os.path.abspath(os.path.dirname(__file__)), relax_radial=0.13, relax_symm_break=0.01, nn_cutoff=2.7)
