import os.path
import vacancy

properties = vacancy.do_all_vacancies(os.path.abspath(os.path.dirname(__file__)), nn_cutoff=2.7)
