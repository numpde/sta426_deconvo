# RA, 2020-12-07

import pandas


def nnls(*, bulk: pandas.DataFrame, scref: pandas.DataFrame) -> pandas.DataFrame:
    """
    Perform deconvolution via NNLS.
    """
    from scipy.optimize import nnls as NNLS
    return pandas.DataFrame(
        data={n: next(iter(NNLS(scref.values, b))) for (n, b) in bulk.iteritems()},
        index=scref.columns,
    )
