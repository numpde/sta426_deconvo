# RA, 2020-12-07
# RA, 2020-12-20

from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from tcga.utils import unlist1, assert_exists, from_iterable, First

# Data folder
datapath = (Path(__file__).parent.parent.parent.parent / "data").resolve()
assert datapath.is_dir()

# Data reader
read_csv = (lambda f: pd.read_csv(assert_exists(f), sep='\t', index_col=0))

# Read all data
datasets = {
    '2020-FGCZ': {
        'data': read_csv(unlist1(datapath.glob("*FGCZ/*.zip"))),
        'meta': read_csv(unlist1(datapath.glob("*FGCZ/*.tsv"))),
    },
    '2015-Darmanis': {
        'data': read_csv(unlist1(datapath.glob("*/*Darmanis/b_*/data.csv.gz"))),
        'meta': read_csv(unlist1(datapath.glob("*/*Darmanis/b_*/meta.csv.gz"))),
    },
    '2019-AllenBrain-M1': {
        'data': read_csv(unlist1(datapath.glob("*/*M1/b_*/data.csv.gz"))),
        'meta': read_csv(unlist1(datapath.glob("*/*M1/b_*/meta.csv.gz"))),
    },
}

# # [genes] x [bulk samples]
# fgcz_meta: pd.DataFrame = datasets['2020-FGCZ']['meta']
# fgcz: pd.DataFrame = datasets['2020-FGCZ']['data']
# fgcz = fgcz.groupby('gene_name').sum()
# fgcz.index = fgcz.index.str.upper()

# [genes] x [single cells]
darm_meta: pd.DataFrame = datasets['2015-Darmanis']['meta']
darm: pd.DataFrame = datasets['2015-Darmanis']['data']
darm.index = darm.index.str.upper()
darm = darm.sort_index(axis=1)  # why?

# [genes] x [single cells],  preprocessed to a marker gene set
abm1_meta: pd.DataFrame = datasets['2019-AllenBrain-M1']['meta']
abm1: pd.DataFrame = datasets['2019-AllenBrain-M1']['data']
assert (str(abm1.index.name) in ["sample_name", "gene_name"])
abm1 = abm1.T if (str(abm1.index.name) == "sample_name") else abm1
abm1.index = pd.Series(name="gene_name", data=abm1.index.str.upper())
abm1 = abm1.sort_index()
# Note: don't sort the columns

assert (141 == len(set(abm1.index)))
assert (118 == len(set(abm1.index) & set(darm.index)))

# Align on the common genes
(darm, abm1) = darm.align(abm1, join='inner', axis=0)
assert (len(darm.index) == 118 == len(abm1.index))


# Cell types in the Darmanis dataset
class darm_celltypes:
    oligodendrocytes = 'oligodendrocytes'
    hybrid = 'hybrid'
    endothelial = 'endothelial'
    neurons = 'neurons'
    microglia = 'microglia'
    astrocytes = 'astrocytes'
    fetal_replicating = 'fetal_replicating'
    OPC = 'OPC'
    fetal_quiescent = 'fetal_quiescent'


# Cell types in the AllenBrain-M1 dataset
class abm1_celltypes:
    Exc = 'Exc'
    Inh = 'Inh'
    Oligo = 'Oligo'
    Astro = 'Astro'
    OPC = 'OPC'
    Micro = 'Micro'
    Endo = 'Endo'
    VLMC = 'VLMC'


def drop_zero_cols(df) -> pd.DataFrame:
    return df.loc[:, (df.sum(axis=0) > 0)]


def norm2_col(df: pd.DataFrame):
    return df * (df * df).sum(axis=0).transform(lambda x: (1 / np.sqrt(x or 1)))


def norm1_col(df) -> pd.DataFrame:
    return df / df.sum(axis=0)


normalize = First(drop_zero_cols).then(norm1_col)
