# RA, 2020-12-21

"""
Visualize each cell from Darmanis:
 - similarity to cell types from [AB]
 - on the t-SNE map of [AB]

[AB] Allen-Brain-M1
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats.kde import gaussian_kde
from matplotlib.colors import to_rgb, TABLEAU_COLORS as COLORS

from progressbar import progressbar as lazy_progressbar

from tcga.utils import unlist1, mkdir
from plox import Plox

from datasource import darm, darm_meta, abm1, abm1_meta, datapath, norm2_col

progressbar = (lambda X: lazy_progressbar(list(X)))
COLORS = list(COLORS.values())

tsne_map = pd.read_csv(unlist1(datapath.glob("*/*M1/c_*/tsne.csv.gz")), sep='\t', index_col=0)
assert tsne_map.index.equals(abm1.columns)

out_dir = mkdir(Path(__file__).with_suffix(''))

darm = norm2_col(darm)
abm1 = norm2_col(abm1)

# https://matplotlib.org/tutorials/introductory/customizing.html
style = {
    'legend.fontsize': "large",

    'axes.labelsize': "large",

    'xtick.labelsize': "small",

    'axes.spines.bottom': False,
    'axes.spines.left': False,
    'axes.spines.right': False,
    'axes.spines.top': False,
}

for (celltype, (m, expr)) in zip(progressbar(darm_meta['cell type'][darm.columns]), darm.iteritems()):
    name = F"{m}_{celltype}.png"

    similarities = expr @ abm1
    assert all((0 <= similarities) & (similarities <= 1))

    # Histograms
    with Plox(style) as px:
        for first_pass in [True, False]:
            type_to_sims = similarities.groupby(abm1_meta.celltype)

            SCORE_TYPES_BY = np.mean
            type_to_score = type_to_sims.aggregate(SCORE_TYPES_BY)
            type_to_score = type_to_score / (max(type_to_score) or 1) * len(type_to_score)

            # For each candidate cell
            # and each reference cell type
            # only keep similarity scores
            # below the quantile:
            SIMILARITY_CUTOFF_HIGH = 0.95

            for ((t, sim), c) in zip(type_to_sims, COLORS):
                sim = sim[sim <= sim.quantile(SIMILARITY_CUTOFF_HIGH)]

                xx = np.linspace(1e-2, 1, 201)
                if sim.std() and (len(sim) >= 2):
                    yy = gaussian_kde(sim).evaluate(xx)
                else:
                    yy = 0 * xx

                # Line width
                lw = 1 + type_to_score[t]

                if first_pass:
                    px.a.plot(xx[xx > 0.1], yy[xx > 0.1], '-', c=c, label=t, lw=lw)
                    ylim = px.a.get_ylim()
                else:
                    (s, a) = ('-', 0.9) if (max(yy[xx <= 0.1]) < max(ylim)) else ('--', 0.3)
                    px.a.plot(xx[xx <= 0.11], yy[xx <= 0.11], s, c=c, alpha=a, lw=lw)
                    px.a.set_ylim(*ylim)

        # px.a.set_xlabel("Similarity")
        px.a.legend(loc="upper right")
        px.a.axis('off')
        px.a.set_xticks([])
        px.a.set_yticks([])

        px.f.savefig(mkdir(out_dir / "hist") / name, dpi=60)

    # t-SNE plots
    with Plox(style) as px:
        for (c, (t, tsne)) in zip(COLORS, tsne_map.groupby(abm1_meta.celltype, axis=0)):
            px.a.scatter(tsne.x.mean(), tsne.y.mean(), marker='o', s=20, edgecolors=c, c='none', label=t)
            c = pd.DataFrame(data=[to_rgb(c)], index=tsne.index).assign(a=similarities)
            px.a.scatter(tsne.x, tsne.y, s=(4e3 / len(tsne)), c=c, edgecolors='none')
            px.a.axis('off')

        px.a.legend(loc="upper left")

        px.f.savefig(mkdir(out_dir / "tsne") / name, dpi=60)
