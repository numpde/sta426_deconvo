# RA, 2020-12-07

"""
Multiple deconvolutions per bulk sample:
bar chart with error bars.
"""

import numpy as np
import pandas as pd
from pathlib import Path

from plox import Plox

from tcga.utils import mkdir, First
from datasource import normalize, norm1
from datasource import darm, fgcz, fgcz_meta, darm_celltypes as celltypes
from nnls import nnls

rs = np.random.RandomState(9)
out_dir = mkdir(Path(__file__).with_suffix(""))

unwanted_celltypes = [celltypes.fetal_quiescent, celltypes.fetal_replicating, celltypes.hybrid]
darm = darm.drop(labels=unwanted_celltypes, axis=1)

darm = normalize(darm)
fgcz = normalize(fgcz)

f = First(norm1).then(lambda s: s.groupby(level=0).sum()).then(np.squeeze)

for (sample, fgcz) in fgcz.iteritems():
    deco = pd.DataFrame(
        data={
            r: f(nnls(bulk=pd.DataFrame(fgcz), scref=darm.sample(frac=0.5, random_state=rs, axis=1)))
            for r in range(20)
        }
    )

    m = deco.mean(axis=1)
    s = deco.std(axis=1)

    style = {
        'xtick.labelsize': "xx-small",
    }

    with Plox(style) as px:
        px.a.bar(m.index, height=m, yerr=s)
        px.a.set_ylim(0, 1)
        px.a.set_title(", ".join(fgcz_meta.loc[sample, ['Age', 'Condition', 'Source']].astype(str)))
        px.f.savefig(out_dir / F"{sample}.png")
