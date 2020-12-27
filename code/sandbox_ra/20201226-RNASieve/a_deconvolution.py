# RA, 2020-12-27

"""
Note:
At the time of writing there is a bug in RNA-Sieve, see
    https://github.com/songlab-cal/rna-sieve/pull/2
To install with a patch, use
    pip3 install git+https://github.com/numpde/rna-sieve.git@92862b1
"""

from pathlib import Path

import pandas as pd
import numpy as np

from tcga.utils import mkdir
from datasource import fgcz, darm, subset_to_markers, darm_celltypes as celltypes

from rnasieve.preprocessing import model_from_raw_counts

out_dir = mkdir(Path(__file__).with_suffix(''))

"""
About `model_from_raw_counts` (from rnasieve/preprocessing.py):

    Given raw counts in the form of a dictionary {label: matrix} and 
    a matrix of bulks, produce a RNASieveModel and corresponding psis.

    Matrices should be sorted by genes over the same set of genes.
    
My mini example:
    ngenes = 10
    rng = np.random.default_rng()
    
    scref = pd.DataFrame([
        pd.Series(name=celltype, data=rng.integers(0, 5, size=ngenes))
        for celltype in ("AB" * 8)
    ]).T
    
    bulk = pd.DataFrame(scref.sum(axis=1)).values
    scref = dict((k, X.values) for (k, X) in scref.groupby(scref.columns, axis=1))
    
    (model, cleaned) = model_from_raw_counts(raw_counts=scref, bulks=bulk)
    x = model.predict(cleaned)
"""

assert darm.index.equals(fgcz.index)
assert not darm.isna().any().any()
assert not fgcz.isna().any().any()

unwanted_celltypes = [celltypes.fetal_quiescent, celltypes.fetal_replicating, celltypes.hybrid]
darm = darm.drop(labels=unwanted_celltypes, axis=1)

# The following leads to:
#   UserWarning: Solution may be inaccurate.
#   Try another solver, [...]
# darm = subset_to_markers(darm)
# fgcz = subset_to_markers(fgcz)

# This is why we use this adhockery:
ii = (darm / darm.sum(axis=0)).std(axis=1).nlargest(30).index
darm = darm.loc[ii]
fgcz = fgcz.loc[ii]

scref = dict((k, X.values) for (k, X) in darm.groupby(darm.columns, axis=1))

(model, cleaned) = model_from_raw_counts(raw_counts=scref, bulks=fgcz.values)

proportions = model.predict(cleaned).T
proportions.columns = fgcz.columns
proportions.index.name = "celltype"

# Sort columns
proportions = proportions.sort_index(axis=0, key=(lambda s: s.str.lower()))

proportions.to_csv(out_dir / F"celltypes.csv", index=True, sep='\t')
