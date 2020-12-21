# RA, 2020-12-21

from twig import log

from datetime import datetime, timezone
from contextlib import redirect_stdout
from pathlib import Path

from tcga.utils import relpath, First, unlist1

from datasource import darm, darm_meta, datapath

time = datetime.now(tz=timezone.utc).strftime("%Z-%Y%m%d-%H%M%S")
this = Path(__file__)
find = First(datapath.glob).then(unlist1).then(relpath)


preambulations = F"""
We compare the scRNA expression from
Darmanis et al [[1]]({find("**/2015-Darmanis")})
with that from
Allen Brain M1 [[2]]({find("**/2019-AllenBrain-M1")}).
The datasets are reduced to the marker genes
that occur in [2], and then further to 
the common set of genes. 
The similarity between two cells is then measured
as the cosine similarity.

This file was generated by [{this.name}]({this.name}) ({time}).

### Correspondence by cell type

The first heatmap shows the average similarity
grouped by cell types.
It is computed for each cell type pair as follows.
First, for each cell (of a given cell type) in [1], 
the similarity to each cell (of a given cell type) in [2]
is computed. 
Then it is averaged over the cells in [2].
Then it is averaged over those in [1].

This has the issue that there could be many
cells in [2] that dilute the similarity.
To counteract this, in the second heatmap,
the 95% quantile is taken instead of the first average.

![](a_celltypes/heatmap_agg=mean.png)
![](a_celltypes/heatmap_agg=qt95.png)


### Correspondence by individual cell

For each cell from [1],
as quoted in the subtitle:

- The first figure shows the distribution
of the cosine similarity measure against cells from [2]. 
The similarities above 95% quantile have been removed.
The thickness of the line scales with 
the average similarity (before removal). 

- The second figure shows a t-SNE plot 
of the cells from [2] colored by their type
(the empty circles indicate the barycenters).
The size of the dots is such that the total
marker area is the same across cell types.
The intensity of the color is proportional
to the similarity.
"""
with redirect_stdout((this.parent / "readme.md").open('w')):
    print(preambulations)
    print("")

    for (t, ii) in sorted(darm_meta.index.groupby(darm_meta['cell type'].str.lower()).items()):
        for i in ii:
            f = First((this.parent / "b_cellwise").glob).then(unlist1).then(relpath).then(Path)

            (hist, tsne) = (f(F"hist/{i}*.png"), f(F"tsne/{i}*.png"))
            assert (hist.stem == tsne.stem)

            log.info(F"{i} ({t})")

            print(F"#### {i} ({t})")
            print("")
            print(F"![{hist.stem}]({hist})")
            print(F"![{tsne.stem}]({tsne})")
            print("")
