# RA, 2020-12-20

from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd

from progressbar import progressbar

from twig import log
from plox import Plox
from tcga.utils import at_most_n, mkdir, relpath

from datasource import abm1, abm1_meta, darm, darm_meta, norm2_col

# Don't need the sample IDs
darm.columns = darm_meta['cell type'][darm.columns]
abm1.columns = pd.Index(name="cell type", data=abm1_meta['celltype'][abm1.columns])

out_dir = mkdir(Path(__file__).with_suffix(''))


def compute_and_plot(agg):
    """
    `agg` is `np.mean` or `np.max` (or similar).

    Procedure:
     o) For each cell in Darmanis:
         - Compute the similarity with each cell in AllenBrain-M1,
         - Aggregate by cell type using `agg`.
     o) Average by cell type.
    """

    out_csv = out_dir / F"heatmap_agg={agg.__name__}.csv"

    if out_csv.is_file():
        log.info(F"Skipping {relpath(out_csv)}.")
    else:
        log.info(F"Making {relpath(out_csv)}.")

        df = pd.DataFrame(data={
            b: {
                a: pd.Series((Y.T @ X).apply(agg, axis=1)).mean()
                for (a, X) in norm2_col(abm1).groupby(abm1.columns, axis=1)
            }
            for (b, Y) in progressbar(list(norm2_col(darm).groupby(darm.columns, axis=1)))
        })

        df.to_csv(out_csv, sep='\t')

    if out_csv.is_file():
        log.info(F"Making figure/s.")

        df = pd.read_csv(out_csv, sep='\t', index_col=0)

        for i in [0, 1]:
            df = df.sort_index(axis=i, key=(lambda x: x.str.lower()))

        # https://matplotlib.org/tutorials/introductory/customizing.html
        style = {
            'legend.fontsize': "xx-small",
            'xtick.labelsize': "xx-small",
            'ytick.labelsize': "xx-small",

            'savefig.pad_inches': 0.01,
        }

        with Plox(style) as px:
            px.a.imshow(df.T, cmap='Greys')
            px.a.tick_params(length=0)
            px.a.set_xticks(range(0, len(df.index)))
            px.a.set_xticklabels([F"{c}\n{Counter(abm1.columns)[c]}" for c in df.index])
            px.a.set_xlabel(F"Allen Brain M1")
            px.a.set_yticks(range(0, len(df.columns)))
            px.a.set_yticklabels([F"{c}\n{Counter(darm.columns)[c]}" for c in df.columns])
            px.a.set_ylabel(F"Darmanis")
            px.f.savefig(out_csv.with_suffix(".png"))


def qt90(s: pd.Series):
    return s.quantile(0.9)


def qt95(s: pd.Series):
    return s.quantile(0.95)


def main():
    for agg in [np.mean, np.max, qt90, qt95]:
        compute_and_plot(agg)


if __name__ == '__main__':
    main()
