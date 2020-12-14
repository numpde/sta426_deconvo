#!/usr/bin/env bash
# RA, 2020-12-14

out_dir=a_original
mkdir -p $out_dir

wget \
  https://github.com/ellispatrick/CortexCellDeconv/raw/master/CellTypeDeconvAnalysis/Data/GSE67835.RData \
  -o $out_dir/readme.txt \
  -O $out_dir/GSE67835.RData

