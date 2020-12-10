# RA, 2020-12-06


import pandas as pd
import numpy as np

import statsmodels.api as sm
from plox import Plox
from pathlib import Path
from tcga.utils import unlist1, assert_exists, first, mkdir

# Output folder
out_dir = mkdir(Path(__file__).with_suffix(""))

# Data reader
read_csv = (lambda f: pd.read_csv(assert_exists(f), sep='\t', index_col=0))

# [bulk samples] x [cell types]
df_cell: pd.DataFrame = read_csv(unlist1(Path(__file__).parent.glob("a_deco*/*.csv*"))).T
df_cell = df_cell.groupby(df_cell.columns, axis=1).sum()

# Meta
df_meta: pd.DataFrame = read_csv(unlist1(Path(__file__).parent.glob("../../../*/*FGCZ/*infos.tsv")))
assert df_cell.index.equals(df_meta.index)

# Combine into one dataframe [Samples] x [Attributes]
df = df_cell.merge(df_meta, left_index=True, right_index=True)

for x in ['Age', 'Source', 'neurons']:
    for cell in df_cell:
        with Plox() as px:
            (x, y) = (x, cell)
            for (n, (k, grp)) in enumerate(df.groupby(df.Condition)):
                # xx = np.linspace(min(grp[x]), max(grp[x]), 101)
                # yy = sm.OLS(grp[y], sm.add_constant(grp[x])).fit().predict(sm.add_constant(xx))

                # # Alternative
                # from sklearn.linear_model import LinearRegression
                # linreg = LinearRegression().fit(grp[x].values.reshape([-1, 1]), grp[y])
                # yy = linreg.predict(xx.reshape([-1, 1]))

                c = F"C{n}"
                px.a.plot(grp[x], grp[y], '.', c=c, label=k)
                # px.a.plot(xx, yy, '--', c=c, lw=0.5)

            px.a.set_xlabel(x)
            px.a.set_ylabel(y)
            px.a.set_ylim(0, max(px.a.get_ylim()))
            px.a.legend()
            px.f.savefig((out_dir / F"{x}__{y}".lower()).with_suffix('.png'))
