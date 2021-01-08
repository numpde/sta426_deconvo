# RA, 2020-12-22
# RA, 2021-01-08

"""
Load deconvolution results of the main R script.
"""

import twig
import datasource

try:
    fgcz_deco_by_bisque = datasource.read("code/**/20201226-AllDeco/output/fgcz_darm_bisque_data.csv")
    fgcz_deco_by_bisque.columns = fgcz_deco_by_bisque.columns.str.replace('.', '-')
except:
    twig.log.warning("Deconvolution results by Bisque not found.")
    fgcz_deco_by_bisque = None

try:
    fgcz_deco_by_music = datasource.read("code/**/20201226-AllDeco/output/fgcz_darm_music_data.csv")
    fgcz_deco_by_music.columns = fgcz_deco_by_music.columns.str.replace('.', '-')
except:
    twig.log.warning("Deconvolution results by MuSiC not found.")
    fgcz_deco_by_music = None

if __name__ == '__main__':
    print(fgcz_deco_by_bisque)
    print(fgcz_deco_by_music)
