# RA, 2020-12-06

import pandas as pd
import numpy as np
from plox import Plox
from pathlib import Path
from tcga.utils import unlist1, assert_exists, first, mkdir, from_iterable, First
from scipy.optimize import nnls as NNLS

# Output folder
out_dir = mkdir(Path(__file__).with_suffix(""))

from datasource import fgcz, darm, subset_to_markers

if True:
    unwanted_celltypes = ['fetal_quiescent', 'fetal_replicating', 'hybrid']
    # unwanted_celltypes = ['hybrid']
    ref_df = darm.drop(labels=unwanted_celltypes, axis=1)

if True:
    ref_df = subset_to_markers(ref_df)
    fgcz_df = subset_to_markers(fgcz)

# # Venn diagram of genes in bulk and ref datasets
# a = set(df.index) - set(ref_df.index)
# b = set(ref_df.index) - set(df.index)
# c = set(ref_df.index) & set(df.index)  # Common genes
# from matplotlib_venn import venn2
# from plox import Plox
# with Plox() as px:
#     venn2(subsets=[len(a), len(b), len(c)], ax=px.a)
#     px.show()

#

# Align on the common genes
(ref_df, df) = ref_df.align(df, join='inner', axis=0)
assert (list(ref_df.index) == list(df.index))

# Drop the reference cells that have zero expression
ref_df = ref_df.loc[:, (ref_df.sum(axis=0) > 0)]

#

# # Collapse cell types # TODO
# ref_df = ref_df.groupby(level=0, axis=1).mean()

# Normalization of the reference dataset
ref_df = ref_df / ref_df.sum(axis=0)

# with Plox() as px:
#     px.a.imshow(ref_df.groupby(ref_df.columns, axis=1).mean())
#     px.show()
#     exit()

# # Normalize ...
# ref_df = ref_df / ref_df.sum(axis=0)
# ref_df = ref_df / (
#         ref_df.sum(axis=0).groupby(ref_df.columns).sum()[ref_df.columns]
#         *
#         len(set(ref_df.columns))
# )
# # ... such that:
# # - The total sums to one
# assert np.isclose(sum(ref_df.sum(axis=0)), 1)
# # - Each cell type,
# for c in set(ref_df.columns):
#     # has equal contribution (1 / #celltypes),
#     assert np.isclose(ref_df[c].sum().sum(), 1 / len(set(ref_df.columns)))
#     # and each cell contributes equally there
#     assert np.isclose(0, ref_df[c].sum(axis=0).std().max())

#

# Normalize the bulk samples
df = df / df.sum(axis=0)

# Subset of bulk samples (for debugging)
df = df.iloc[:, 0:1]

# # Pseudobulk
# print(ref_df)
# df = ref_df @ pd.DataFrame(data={'pseudobulk1': [1, 2, 3, 4, 5, 6]}, index=ref_df.columns)

# print(ref_df.sum(axis=0))
# print(df.sum(axis=0))

# Perform deconvolution via NNLS
dec_df: pd.DataFrame
dec_df = pd.DataFrame(
    data={
        col: first(NNLS(ref_df.values, df[col]))
        for col in df.columns
    },
    index=ref_df.columns,
)

# # Residual
# print((df - ref_df @ dec_df).sum(axis=0))

# Cell type proportions
dec_df = dec_df / dec_df.sum(axis=0)

# Save to disk
dec_df.to_csv(out_dir / "celltypes.csv", compression=None, sep='\t')

