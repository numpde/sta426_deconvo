# RA, 2020-12-06

"""
Barebones baseline implementation
of deconvolution by NNLS.

Reference dataset: Darmanis et al (2015).
"""

from pathlib import Path
from tcga.utils import mkdir
from datasource import fgcz, darm, normalize, norm1, darm_celltypes as celltypes
from nnls import nnls

unwanted_celltypes = [celltypes.fetal_quiescent, celltypes.fetal_replicating, celltypes.hybrid]
darm = darm.drop(labels=unwanted_celltypes, axis=1)

darm = normalize(darm)
fgcz = normalize(fgcz)

# Deconvolution
deco = nnls(bulk=fgcz, scref=darm)

# Cell type proportions
deco = norm1(deco)

# Save to disk
out_dir = mkdir(Path(__file__).with_suffix(""))
deco.to_csv(out_dir / "celltypes.csv", sep='\t')
