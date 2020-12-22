# RA, 2020-12-22

import twig
import datasource

# Attempt to read the deconvolution results by Bisque
try:
    fgcz_deco_by_bisque = datasource.read("code/**/20201222-Bisque/b_*/fgcz_darm_data.csv")
    fgcz_deco_by_bisque.columns = fgcz_deco_by_bisque.columns.str.replace('.', '-')
except:
    twig.log.warning("Deconvolution results by Bisque not found.")
    fgcz_deco_by_bisque = None

if __name__ == '__main__':
    print(fgcz_deco_by_bisque)
