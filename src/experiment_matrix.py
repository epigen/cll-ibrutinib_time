import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)

df = pd.read_csv("metadata/annotation.csv")

df = df[df['library'] == "ATAC-seq"]

# replace failed parameters with -1 and fill not-assessed with 0
df.loc[:, ['pass_qc', 'pass_qc_TK', 'pass_counts', "pass_corr"]] = df.loc[:, ['pass_qc', 'pass_qc_TK', 'pass_counts', "pass_corr"]].replace(0, -1).fillna(0)
# tag samples as failed in containing any failure: 1) pass; 0) fail
df['pass'] = (~(df.loc[:, ['pass_qc', 'pass_qc_TK', 'pass_counts', "pass_corr"]] == -1).any(axis=1)).astype(int)
# pivot and marked non-existing/non-processed samples with -1
df_pivot = pd.pivot_table(data=df, columns='patient_id', index=["timepoint", "cell_type"], values="pass", aggfunc=sum).fillna(-1)

ratio = df_pivot.shape[0] / float(df_pivot.shape[1])
fig, axis = plt.subplots(1, figsize=(4 * ratio, 4))
sns.heatmap(df_pivot.T, cmap="RdYlGn", ax=axis, square=False)
axis.set_xticklabels(axis.get_xticklabels(), rotation=90)
fig.savefig("experiment_matrix.svg", bbox_inches="tight")
