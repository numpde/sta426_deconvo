# RA, 2020-12-07

"""
Multiple deconvolutions per bulk sample in one triangle.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from progressbar import progressbar

from tcga.utils import mkdir

from datasource import fgcz, darm, normalize, fgcz_meta, darm_celltypes as celltypes
from nnls import deco3

out_dir = mkdir(Path(__file__).with_suffix(""))

unwanted_celltypes = [celltypes.fetal_quiescent, celltypes.fetal_replicating, celltypes.hybrid]
darm = darm.drop(labels=unwanted_celltypes, axis=1)

darm = normalize(darm)
fgcz = normalize(fgcz)

# Reduce the number of reference cell types to /three/
darm.columns = [(i if i in ['astrocytes', 'neurons', 'others'] else 'others') for i in darm.columns]


def scref(frac=0.3, repeats=1000, rs=np.random.RandomState(43)):
    if (frac == 1):
        yield darm
    else:
        for __ in range(repeats):
            yield darm.sample(frac=frac, random_state=rs, axis=1)


def qc(reco_prop):
    explained_fraction = sum(reco_prop)
    mode = max(reco_prop)
    return (explained_fraction > 0.5) and (mode < 0.9)


def main():
    for (n, (kind, df)) in enumerate(fgcz.groupby(fgcz_meta.Condition, axis=1)):
        for (sample_id, sample) in progressbar(list(df.iteritems())):
            for frac in [0.3, 0.5, 0.8]:
                with deco3(bulk=sample, scref=scref(frac), qc=qc) as px:
                    name = "_".join(map(str, [
                        sample_id,
                        fgcz_meta.Condition[sample_id],
                        fgcz_meta.Source[sample_id],
                        fgcz_meta.Gender[sample_id],
                        fgcz_meta.Age[sample_id],
                    ]))

                    px.f.savefig((mkdir(out_dir / F"frac={frac}") / name).with_suffix(".png"))


if __name__ == '__main__':
    main()
