# RA, 2020-12-07


import contextlib
import typing

import pandas as pd
import numpy as np

from scipy.stats import gaussian_kde
from scipy.spatial import ConvexHull
import scipy.optimize

from plox import Plox
from tcga.utils import First, first

from datasource import norm1


def nnls(*, bulk: pd.Series = None, scref: pd.DataFrame = None) -> pd.Series:
    """
    Perform deconvolution via NNLS.
    """

    if (bulk is None) and (scref is not None):
        return (lambda x: nnls(bulk=x, scref=scref))
    if (bulk is not None) and (scref is None):
        return (lambda x: nnls(bulk=bulk, scref=x))

    assert (bulk is not None)
    assert (scref is not None)

    x = pd.Series(data=first(scipy.optimize.nnls(scref.values, bulk.values)), index=scref.columns)
    x = norm1(x) * ((scref @ x).sum() / bulk.sum())

    return x


corners = np.array([[0, 0], [1, 0], [0.5, 1]])


@contextlib.contextmanager
def deco3(bulk: pd.Series, scref: typing.Iterable[pd.DataFrame], qc=None) -> Plox:
    with Plox() as px:
        f = First(nnls(bulk=bulk)).then(lambda s: s.groupby(s.index).sum().sort_index().squeeze())

        mm = pd.DataFrame(
            data={
                r: norm1(reco_prop)
                for (r, reco_prop) in enumerate(map(f, scref))
                if (qc is None) or qc(reco_prop)
            }
        )

        # Collapse to /three/ cell types
        mm = mm.groupby(mm.index).sum().sort_index()
        assert (3 == len(mm.index))

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
        for (p, c) in zip(corners, mm.index):
            px.a.text(*p, s=c, va=('top' if (p[1] == 0) else 'bottom'), ha='center')

        # px.a.legend(loc="upper right")
        px.a.axis('off')

        yield px
