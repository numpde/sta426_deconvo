# RA, 2020-12-10

"""
Multiple deconvolutions in one triangle
on synthetic pseudobulk samples.
"""

from pathlib import Path

import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from progressbar import progressbar

from datasource import norm1
from tcga.utils import unlist1, mkdir

from a_deconvolution3 import scref
from nnls import deco3, corners

out_dir = mkdir(Path(__file__).with_suffix(""))

full = unlist1(scref(frac=1))


def qc(reco_prop):
    explained_fraction = sum(reco_prop)
    mode = max(reco_prop)
    return (explained_fraction > 0.9 > mode)


def experiment(i):
    rs = np.random.RandomState(i)

    igroup = (lambda s: s.groupby(s.index))
    to_rnd = (lambda x: rs.random())

    # Sample equally each cell type and within each cell type
    full_prop = pd.Series(index=full.columns, data=0).transform(to_rnd)
    full_prop = full_prop / igroup(full_prop).agg(np.sum)[full_prop.index]
    full_prop = full_prop * igroup(full_prop).agg(to_rnd)[full_prop.index]
    full_prop = norm1(full_prop)

    sample = full @ full_prop
    full_prop = igroup(full_prop).sum().sort_index()

    for frac in [0.3, 0.5, 0.8]:
        with deco3(bulk=sample, scref=scref(frac, rs=rs), qc=qc) as px:
            name = F"i={i:04}"
            px.a.plot(*(corners.T @ full_prop), 'kx', ms=18)
            px.f.savefig((mkdir(out_dir / F"frac={frac}") / name).with_suffix(".png"))


if __name__ == '__main__':
    Parallel(n_jobs=6)(delayed(experiment)(i) for i in progressbar(list(range(20))))
