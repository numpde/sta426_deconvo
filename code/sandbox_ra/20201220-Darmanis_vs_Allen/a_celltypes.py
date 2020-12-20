# RA, 2020-12-20

from pathlib import Path
from twig import log

import numpy as np
import pandas as pd

from progressbar import progressbar

from plox import Plox
from tcga.utils import at_most_n, mkdir, relpath

from datasource import abm1, darm, drop_zero_cols

out_dir = mkdir(Path(__file__).with_suffix(''))
out_csv = out_dir / "heatmap.csv"

if out_csv.is_file():
    log.info(F"Skipping making {relpath(out_csv)}.")
else:
    log.info(F"Making {relpath(out_csv)}.")

    def norm2_col(df: pd.DataFrame):
        return df * (df * df).sum(axis=0).transform(lambda x: (1 / np.sqrt(x or 1)))


    def dot(x: pd.Series, y: pd.Series):
        return np.dot(x, y)


    df = pd.DataFrame(data={
        b: {
            a: np.mean([
                dot(x, y)
                for (_, x) in at_most_n(X.iteritems(), n=-1)
                for (_, y) in at_most_n(Y.iteritems(), n=-1)
            ])
            for (a, X) in norm2_col(abm1).groupby(abm1.columns, axis=1)
        }
        for (b, Y) in progressbar(list(norm2_col(darm).groupby(darm.columns, axis=1)))
    })

    df.to_csv(out_csv, sep='\t')

if out_csv.is_file():
    log.info(F"Making figure/s.")

    df = pd.read_csv(out_csv, sep='\t', index_col=0)

    # https://matplotlib.org/tutorials/introductory/customizing.html
    style = {
        'legend.fontsize': "xx-small",
        'xtick.labelsize': "xx-small",
        'ytick.labelsize': "xx-small",

        'savefig.pad_inches': 0.01,
    }

    with Plox(style) as px:
        px.a.imshow(df.T)
        px.a.tick_params(length=0)
        px.a.set_xticks(range(0, len(df.index)))
        px.a.set_xticklabels(df.index)
        px.a.set_yticks(range(0, len(df.columns)))
        px.a.set_yticklabels(df.columns)
        px.f.savefig(out_csv.with_suffix(".png"))
