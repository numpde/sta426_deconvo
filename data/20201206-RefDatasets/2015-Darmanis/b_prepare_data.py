# RA, 2020-12-06


URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE67835&format=file"

import typing
import contextlib
import datetime
import tarfile
import pandas as pd
from pathlib import Path
from tcga.utils import download

path_original = Path(__file__).parent / "a_original"
path_prepared = Path(__file__).parent / "b_prepared"

download = download.to(abs_path=path_original)


def get_data() -> typing.Iterable[pd.Series]:
    with download(URL).now.open(mode='rb') as fd:
        with tarfile.open(fileobj=fd) as tf:
            for file in tf.getmembers():
                with tf.extractfile(file) as fd:
                    s = pd.read_csv(
                        fd,
                        compression="gzip", sep='\t',
                        names=["gene_name", file.name],
                        index_col=0,
                        squeeze=True,
                    )

                    s = s.astype(int)
                    s.index = s.index.str.strip()

                    assert isinstance(s, pd.Series)
                    assert s.index.is_unique
                    assert (s.dtype == int)

                    # Process a name like
                    # GSM1657871_1772078217.C03.csv.gz
                    s.name = str(s.name).split('_')[0]

                    yield s


if __name__ == '__main__':
    df = pd.DataFrame(data={s.name: s for s in get_data()})
    assert (466 == len(df.columns))

    assert df.index.is_unique
    assert df.columns.is_unique

    df.to_csv(path_prepared / "data.csv.gz", compression="gzip", sep='\t')

    ts = datetime.datetime.now(tz=datetime.timezone.utc).strftime("%Z-%Y%m%d-%H%M%S")
    with open(path_prepared / "data_readme.txt", mode='w') as fd:
        with contextlib.redirect_stdout(fd):
            print("Source:    ", URL)
            print("Local copy:", download(URL).now.local_file.name)
            print("Datetime:  ", ts)
