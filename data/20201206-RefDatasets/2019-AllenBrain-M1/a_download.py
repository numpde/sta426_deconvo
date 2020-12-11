# RA, 2020-12-10

"""
For references, see
http://portal.brain-map.org/atlases-and-data/rnaseq/protocols-human-cortex
[https://archive.ph/RQTwT]
"""

import pandas as pd
from pathlib import Path

from tcga.utils import download

URLS = {
    'expr': "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv",
    'meta': "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv",
}

download = download.to(abs_path=Path(__file__).with_suffix(''))

for (k, url) in URLS.items():
    print(download(url).now.meta)

with download(URLS['expr']).now.open() as fd:
    df = pd.read_csv(fd, sep=',', nrows=10, index_col=0).astype(int)
    assert (df.shape == (len(df), 50281))

with download(URLS['meta']).now.open() as fd:
    df = pd.read_csv(fd, sep=',', index_col=0)
    assert (df.shape == (len(df), 38))

    print(list(df.columns))
    print(df.cluster_label[df.class_label == 'Non-Neuronal'].value_counts())
    exit()
    ct = df.cell_type_designation_label
    print(dict(ct[~ct.str.startswith("Neuron")].value_counts()))

# TODO: save a reduced dataset?
