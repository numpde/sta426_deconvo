# RA, 2020-11-30

"""
Exploratory analysis of the FGCZ dataset.

RA, 2020-11-30
"""

from collections import Counter

import pandas as pd
import numpy as np
from pathlib import Path
from plox import Plox
from collections import OrderedDict
from tcga.utils import mkdir, from_iterable

rs = np.random.RandomState(seed=0)

# Marker genes from
# https://github.com/sta426hs2020/material/blob/8c57e3b/week13-07dec2020/workflow.Rmd#L152
markers = list(map(str.upper, from_iterable(dict(
    neuronal=["Snap25", "Stmn2", "Syn1", "Rbfox3", "Dlg4"],
    neuronal_excitatory=["Slc17a7", "Camk2a", "Grin2b"],
    neuronal_inhibitory=["Gad1", "Lhx6"],
    astrocytes=["Aqp4", "Gfap", "Fgfr3", "Dio2"],
    endothelial=["Cldn5", "Nostrin", "Flt1"],
    microglia=["C1qb", "Tyrobp", "P2ry12", "Csf1r", "Irf8"],
    oligodendrocyte=["Opalin", "Plp1", "Mag", "Mog"],
    OPC=["Pdgfra", "Sox6", "Bcan"],
).values())))

out_dir = mkdir(Path(__file__).with_suffix(''))

data_root = (Path(__file__).parent.parent.parent.parent / "data")
assert data_root.is_dir()

fgcz_data = list(data_root.glob("*FGCZ/*.zip")).pop()
fgcz_meta = list(data_root.glob("*FGCZ/*_infos.tsv")).pop()

fgcz_data = pd.read_csv(fgcz_data, sep='\t')
fgcz_meta = pd.read_csv(fgcz_meta, sep='\t')

print("Sample condition:", dict(Counter(fgcz_meta.Condition)))
print("Sample source:", dict(Counter(fgcz_meta.Source)))

assert fgcz_data["Feature ID"].is_unique

fgcz_data = fgcz_data.drop(columns="Feature ID").set_index("gene_name")

df = fgcz_data  #.transform(lambda x: np.log(1 + x))
df = df / df.sum(axis=0)

from sklearn.manifold import TSNE
X = pd.DataFrame(
    index=fgcz_meta.index,
    columns=['x', 'y'],
    data=TSNE(random_state=rs).fit_transform(df[df.index.isin(markers)].T),
)

# Round the age column
fgcz_meta.Age = fgcz_meta.Age.apply(lambda x: F"{np.round(x, -1)}s")
fgcz_meta = fgcz_meta.sort_values(by='Age')

style = {'savefig.pad_inches': 0.01}

for col in ['Condition', 'Gender', 'Age', 'Source']:
    with Plox(style) as px:
        markerstyles = ['^', 'v', 'p', 's', 'o']
        for ((key, df), m) in zip(X.groupby(fgcz_meta[col]), markerstyles):
            px.a.plot(df.x, df.y, ls='none', label=key, marker=m)
            px.a.legend()
        px.a.axis('off')
        px.f.savefig(out_dir / F"tsne_by_{col}.png")

for col in ['Gender', 'Age', 'pmDelay', 'LibConc_100_800bp', 'RIN', 'Source']:
    with Plox(style) as px:
        data = OrderedDict(sorted(fgcz_meta[col].groupby(fgcz_meta.Condition)))
        px.a.hist(data.values(), label=list(data.keys()))
        px.a.set_xlabel(col)
        px.a.legend()
        px.f.savefig(out_dir / F"hist_by_{col}.png")
