# RA, 2020-12-07

"""
Multiple deconvolutions per bulk sample in one triangle.
"""

import numpy as np
import pandas as pd
from pathlib import Path

from plox import Plox

from tcga.utils import mkdir, First
from datasource import fgcz, darm, normalize, norm1, fgcz_meta, darm_celltypes as celltypes
from nnls import nnls

from progressbar import progressbar

from scipy.stats import gaussian_kde
from scipy.spatial import ConvexHull

rs = np.random.RandomState(9)
out_dir = mkdir(Path(__file__).with_suffix(""))

unwanted_celltypes = [celltypes.fetal_quiescent, celltypes.fetal_replicating, celltypes.hybrid]
darm = darm.drop(labels=unwanted_celltypes, axis=1)

darm = normalize(darm)
fgcz = normalize(fgcz)

corners = np.array([[0, 0], [1, 0], [0.5, 1]])
cells = sorted(['astrocytes', 'neurons', 'others'])

# Reduce the number of reference cell types to /three/
darm.columns = [(i if i in cells else 'others') for i in darm.columns]

for frac in [0.3, 0.5, 0.8]:
    groupby = 'Condition'
    for (n, (kind, df)) in enumerate(fgcz.groupby(fgcz_meta[groupby], axis=1)):
        for (sample, df) in progressbar(list(df.iteritems())):
            with Plox() as px:
                f = First(norm1).then(lambda s: s.groupby(level=0).sum()).then(np.squeeze)

                mm = pd.DataFrame(
                    data={
                        r: f(nnls(bulk=pd.DataFrame(df), scref=darm.sample(frac=frac, random_state=rs, axis=1)))
                        for r in range(1000)
                    }
                )

                # Collapse to /three/ cell types
                mm = mm.groupby(mm.index).sum().sort_index()

                # pp = [xx, yy] - coordinates of the point cloud
                pp = (corners.T @ mm).values

                # Contours
                zz = gaussian_kde(pp).evaluate(pp)
                px.a.tricontourf(pp[0], pp[1], zz, levels=12, cmap="gray_r")
                # px.a.tricontour(pp[0], pp[1], zz, levels=5, colors='k', linewidths=0.4)

                # Convex hull
                for simplex in ConvexHull(pp.T).simplices:
                    px.a.plot(pp[0, simplex], pp[1, simplex], c='k', lw=0.3, alpha=0.1)

                # # Point cloud
                # px.a.plot(*pp, '.', c='k', alpha=0.01, mec='none', zorder=(-n * 1000))

                # Big triangle
                px.a.plot(list(corners.T[0]) * 2, list(corners.T[1]) * 2, c='k', lw=0.5)
                # Label the corners of the big triangle
                for (p, c) in zip(corners, cells):
                    px.a.text(*p, s=c, va=('top' if (p[1] == 0) else 'bottom'), ha='center')

                # px.a.legend(loc="upper right")
                px.a.axis('off')

                # Name for the figure file
                name = "_".join(map(str, [
                    sample,
                    fgcz_meta.Condition[sample],
                    fgcz_meta.Source[sample],
                    fgcz_meta.Gender[sample],
                    fgcz_meta.Age[sample],
                ]))

                px.f.savefig((mkdir(out_dir / F"frac={frac}") / name).with_suffix(".png"))
