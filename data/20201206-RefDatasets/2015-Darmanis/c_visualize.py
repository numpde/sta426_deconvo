# RA, 2020-12-06

import numpy as np
import pandas as pd
from plox import Plox
from sklearn.manifold import TSNE

from pathlib import Path

out_dir = Path(__file__).parent / "c_visualized"

datafile = list(Path(__file__).parent.glob("b*/data*.csv.gz")).pop()
metafile = list(Path(__file__).parent.glob("b*/meta*.csv.gz")).pop()

df = pd.read_csv(datafile, sep='\t', index_col=0).astype(int)
df = df.drop(index=["ambiguous", "alignment_not_unique", "no_feature"])

df_meta = pd.read_csv(metafile, sep='\t', index_col=0)['cell type']

# Normalize
# df = df / df.sum(axis=0)
df = np.log(df + 1)

X = pd.DataFrame(
    index=df.columns,
    columns=['x', 'y'],
    data=TSNE(random_state=2).fit_transform(df.T),
)

# https://matplotlib.org/tutorials/introductory/customizing.html
style = {
    'legend.fontsize': "xx-small",
    'legend.framealpha': 0.5,
}

with Plox(style) as px:
    for (key, df) in X.groupby(df_meta):
        px.a.plot(df.x, df.y, '.', label=key)
    px.a.legend()
    px.a.axis('off')
    px.f.savefig(out_dir / "tsne.png")
