import os.path, interstitial

properties = interstitial.do_interstitial(
    os.path.abspath(os.path.dirname(__file__)), nn_cutoff=2.7
)
