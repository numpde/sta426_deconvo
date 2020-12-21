# RA, 2020-12-06

URL = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE67nnn/GSE67835/miniml/GSE67835_family.xml.tgz"

import contextlib
import datetime
import typing
import tarfile
import pandas as pd

from pathlib import Path
from tcga.utils import download

path_original = Path(__file__).parent / "a_original"
path_prepared = Path(__file__).parent / "b_prepared"

download = download.to(abs_path=path_original)


def get_data() -> typing.Iterable[pd.Series]:
    with download(URL).time.open(mode='rb') as fd:
        with tarfile.open(fileobj=fd) as tf:
            for file in tf.getmembers():
                with tf.extractfile(file) as fd:
                    import xmltodict
                    d = xmltodict.parse(fd.read())
                    d = d['MINiML']
                    # print(d['Contributor'])
                    # print(d['Database'])
                    # print(d['Platform'])
                    for sample in d['Sample']:
                        assert sample['Channel-Count'] == '1'
                        yield pd.Series(
                            data={
                                **dict(pd.DataFrame(sample['Channel']['Characteristics']).values),
                                **dict(pd.DataFrame(sample['Relation']).values),
                            },
                            name=sample['@iid'],
                        )


if __name__ == '__main__':
    df = pd.DataFrame(get_data())
    assert (466 == len(df.index))

    assert df.index.is_unique
    assert df.columns.is_unique

    df.to_csv(path_prepared / "meta.csv.gz", compression="gzip", sep='\t')

    ts = datetime.datetime.now(tz=datetime.timezone.utc).strftime("%Z-%Y%m%d-%H%M%S")
    with open(path_prepared / "meta_readme.txt", mode='w') as fd:
        with contextlib.redirect_stdout(fd):
            print("Source:    ", URL)
            print("Local copy:", download(URL).time.local_file.name)
            print("Datetime:  ", ts)
