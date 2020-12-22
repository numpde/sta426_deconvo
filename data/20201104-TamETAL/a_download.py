# RA, 2020-11-04

"""
References:

[1]
Postmortem cortex samples identify molecular subtypes of ALS:
retrotransposon activation, oxidative stress, and activated glia.
Tam et al. incl. The NYGC ALS Consortium.
Cell Rep, 2019.
https://www.sciencedirect.com/science/article/pii/S221112471931263X

[GSE124439]
RNA-seq analysis of ALS patient samples and
individuals without neurological disorders
Contributor: NYGC ALS Consortium
https://www.nygenome.org/als-consortium/
"""

from tcga.utils import download
from pathlib import Path

download = download.to(abs_path=Path(__file__).with_suffix(''))

URLS = {
    # GSE124439_RAW.tar
    'GSE124439': "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE124439&format=file",
    # GSE124439_series_matrix.txt.gz
    'GSE124439_series': "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124439/matrix/GSE124439_series_matrix.txt.gz",
}

for url in URLS.values():
    data = download(url).now
    print(data.meta)

import tarfile
import pandas as pd

# How to load all samples
with download(URLS['GSE124439']).now.open(mode='rb') as fd:
    with tarfile.TarFile(fileobj=fd, mode='r') as tar:
        df = pd.concat(
            axis=1,
            ignore_index=False,
            verify_integrity=True,
            join="outer",
            objs=(
                pd.read_csv(tar.extractfile(f), compression="gzip", sep='\t', quotechar='"', index_col=0).astype(int)
                for f in tar
            )
        ).sort_index()

print(df)
