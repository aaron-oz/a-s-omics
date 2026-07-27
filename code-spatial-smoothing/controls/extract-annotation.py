#!/usr/bin/env python
# Control C: extract MOSTA E16.5_E1S3 per-bin annotations + spatial coords, clip
# to our FOV, write a compact CSV aligned to the paper's bin-50 grid.
# Reads backed (does not load the 5.3 GB expression matrix into memory).
import sys, numpy as np, pandas as pd, anndata as ad

path = "/var/home/aoz/data/mosta/E16.5_E1S3.MOSTA.h5ad"
out  = "/var/home/aoz/data/mosta/e16-fov-annotation.csv"
# our sub.bound FOV (same coordinate system as the GEM_bin1 raw data)
XLO, XHI, YLO, YHI = 15000, 25000, 11500, 18000

A = ad.read_h5ad(path, backed="r")
print("n_obs:", A.n_obs, "n_vars:", A.n_vars)
print("obs columns:", list(A.obs.columns))
print("obsm keys:", list(A.obsm.keys()))

# find annotation column
ann_col = next((c for c in ["annotation","Annotation","celltype","cell_type","region","tissue","domain"]
                if c in A.obs.columns), None)
print("annotation column ->", ann_col)
if ann_col:
    print("annotation categories:", list(pd.Index(A.obs[ann_col]).astype(str).unique())[:40])

# spatial coords: prefer obsm['spatial'], else obs x/y
if "spatial" in A.obsm:
    xy = np.asarray(A.obsm["spatial"])[:, :2]
    x, y = xy[:,0], xy[:,1]
else:
    xcol = next(c for c in ["x","X","array_col","imagecol"] if c in A.obs.columns)
    ycol = next(c for c in ["y","Y","array_row","imagerow"] if c in A.obs.columns)
    x, y = A.obs[xcol].to_numpy(), A.obs[ycol].to_numpy()
print("coord ranges: x [%.0f, %.0f]  y [%.0f, %.0f]" % (x.min(), x.max(), y.min(), y.max()))

df = pd.DataFrame({"x": np.asarray(x), "y": np.asarray(y),
                   "annotation": pd.Index(A.obs[ann_col]).astype(str) if ann_col else "NA"})
fov = df[(df.x>=XLO)&(df.x<=XHI)&(df.y>=YLO)&(df.y<=YHI)].copy()
print("bins in FOV window:", len(fov), "of", len(df))
print("FOV annotation counts:\n", fov.annotation.value_counts().head(25))
fov.to_csv(out, index=False)
print("wrote", out)
