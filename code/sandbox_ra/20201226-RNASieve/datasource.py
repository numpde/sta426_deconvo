# RA, 2020-12-07

from twig import log

import pandas as pd

from pathlib import Path
from tcga.utils import unlist1, assert_exists, from_iterable, First

log.info(F"Datasource: {Path(__file__).parent.name}")

# Data folder
datapath = Path(__file__).parent.parent.parent.parent.resolve()
assert datapath.is_dir()

# CSV reader
read = First(datapath.glob).then(unlist1).then(lambda f: pd.read_csv(f, sep='\t', index_col=0))

# Read data
datasets = {
    '2020-FGCZ': {
        'data': read("data/**/20201128-FGCZ/*.zip"),
        'meta': read("data/**/20201128-FGCZ/*.tsv"),
    },
    '2015-Darmanis': {
        'data': read("data/**/2015-Darmanis/b_*/data.csv.gz"),
        'meta': read("data/**/2015-Darmanis/b_*/meta.csv.gz"),
    },
    '2019-AllenBrain-M1': {
        # 'data': read("*/*M1/b_*/data.csv.gz"),
        # 'meta': read("*/*M1/b_*/meta.csv.gz"),
        'mark': read("data/**/*M1/b_*/marker_genes.csv"),
    },
}

# [genes] x [bulk samples]
fgcz_meta: pd.DataFrame = datasets['2020-FGCZ']['meta']
fgcz: pd.DataFrame = datasets['2020-FGCZ']['data']
fgcz = fgcz.groupby('gene_name').sum()
fgcz.index = fgcz.index.str.upper()

# [genes] x [single cells]
darm_meta: pd.DataFrame = datasets['2015-Darmanis']['meta']
darm: pd.DataFrame = datasets['2015-Darmanis']['data']
darm.index = darm.index.str.upper()
darm.columns = datasets['2015-Darmanis']['meta']['cell type']
darm = darm.sort_index(axis=1, key=(lambda s: s.str.lower()))

# Align on the common genes
(darm, fgcz) = darm.align(fgcz, join='inner', axis=0)
assert (list(darm.index) == list(fgcz.index))


# Cell types in the reference dataset
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


markers = sorted(datasets['2019-AllenBrain-M1']['mark'].index.str.upper())


def subset_to_markers(df) -> pd.DataFrame:
    return df.reindex(df.index.intersection(markers))


def drop_zero_cols(df) -> pd.DataFrame:
    return df.loc[:, (df.sum(axis=0) > 0)]


def norm1(df) -> pd.DataFrame:
    return df / df.sum(axis=0)


normalize = First(drop_zero_cols).then(norm1).then(subset_to_markers).then(drop_zero_cols)

if __name__ == '__main__':
    print(fgcz)
    print(darm)
