from phonons import *

# properties = do_phonons(['bulk_black', 'bulk_black_primitive', 'bulk_betaP4', 'bulk_Hittorf', 'bulk_fibrous'], n_supercell=[3,4,2,1,2], band_paths=['GX','GX','GX','GX','GX'])
properties = do_phonons(
    ["bulk_phosphorene"], n_supercell=[[6, 1, 6]], band_paths=["GX"]
)
