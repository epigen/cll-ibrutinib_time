import magic
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

scdata = magic.mg.SCData.from_csv(os.path.join(os.environ['PROCESSED'], 'cll-time_course/results/single_cell_RNA/15_Magic/matrix.csv.gz'),data_type='sc-seq', normalize=False)
scdata.run_magic(n_pca_components=20, random_pca=True, t=6, k=30, ka=10, epsilon=1, rescale_percent=99)

scdata.magic.data.to_csv('~/projects_shared/cll-time_course/results/single_cell_RNA/15_Magic/magic.matrix.csv')
