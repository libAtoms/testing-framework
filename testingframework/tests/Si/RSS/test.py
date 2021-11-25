import os.path
from RSS import do_RSS

properties = do_RSS(os.path.join(os.path.abspath(os.path.dirname(__file__)),'random_structs.extxyz'),':400:2')
