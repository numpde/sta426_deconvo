# RA, 2020-11-30

"""
Preliminary exploratory analysis
of the FGCZ dataset.

RA, 2020-11-30
"""

from collections import Counter

import pandas as pd
import numpy as np
from pathlib import Path
from plox import Plox
from collections import OrderedDict

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

df = fgcz_data.transform(lambda x: np.log(1 + x))
df = df / df.sum(axis=0)

# genes = ["C9orf72", "PRPH", "UNC13A", "SLC1A2", "SNCG", "NEFH", "SOD2", "APOE", "LOX", "RNASE2"]
# df = df[df.index.isin(genes)]

from sklearn.manifold import TSNE
X = pd.DataFrame(
    index=fgcz_meta.index,
    columns=['x', 'y'],
    data=TSNE().fit_transform(df.T),
)

with Plox() as px:
    for (key, df) in X.groupby(fgcz_meta.Condition):
        px.a.plot(df.x, df.y, '.', label=key)
        px.a.legend()
    px.show()

for col in ['Gender', 'Age', 'pmDelay', 'LibConc_100_800bp', 'RIN']:
    with Plox() as px:
        data = OrderedDict(sorted(fgcz_meta[col].groupby(fgcz_meta.Condition)))
        px.a.hist(data.values(), label=list(data.keys()))
        px.a.set_xlabel(col)
        px.a.legend()
        px.show()
