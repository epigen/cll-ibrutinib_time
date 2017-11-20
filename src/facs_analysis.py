import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


# Setup
df = pd.read_csv("metadata/facs.cell_type.quantification.csv").set_index(["patient_id", "timepoint"])
patients = df.index.levels[0]
markers = df.columns[~df.columns.str.islower()]
pallettes = sns.color_palette("colorblind") + sns.color_palette("Set1")[::-1] + sns.color_palette("colorblind") + sns.color_palette("Set1")[::-1]


df = df[df.index.get_level_values('timepoint') != 2]

for variable, units, ylog in [("absolute", "Cells / uL", True), ("percentage", "% live cells", False)]:
    df2 = df[df['value_type'] == variable]

    # For each marker, see the abundance of each cell type dependent on time
    df_melted = pd.melt(df2[markers].reset_index(), id_vars=["patient_id", "timepoint"], var_name="cell_type", value_name=variable).rename(columns={variable: units})

    g = sns.factorplot(
        data=df_melted, x="timepoint", y=units, col="cell_type",
        size=2, sharey=False, ci=75, col_wrap=3, legend_out=True, palette="viridis", # color=sns.color_palette("colorblind")[0],
        errwidth=2)
    g.savefig(os.path.join("results", "facs.cell_type.marker_per_time.{}.line.svg".format(variable)), bbox_inches="tight")

    g = sns.factorplot(
        data=df_melted, x="timepoint", y=units, col="cell_type",
        size=2, sharey=False, ci=75, col_wrap=3, legend_out=True, palette="viridis", # color=sns.color_palette("colorblind")[0],
        kind="bar", errwidth=2)
    g.savefig(os.path.join("results", "facs.cell_type.marker_per_time.{}.bar.svg".format(variable)), bbox_inches="tight")

    # For each marker, plot patient evolution along timepoints
    n_rows = n_cols = int(np.ceil(np.sqrt(len(markers))))
    fig, axis = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 2), sharex=True, sharey=False)
    axis = axis.flatten()
    colors = dict(zip(patients, pallettes))
    for i, marker in enumerate(markers):
        for j, patient in enumerate(patients):
            series = df2.loc[(patient,), marker]
            axis[i].scatter(
                series.index, series,
                s=25, label=patient, color=colors[patient])
            axis[i].plot(
                series.index,
                series.fillna(method='pad'),
                linestyle="--", label=patient, color=colors[patient], lw=1)
            axis[i].set_xscale('symlog', basex=2)
            if ylog:
                axis[i].set_yscale('symlog', basex=10)
            axis[i].set_title(marker)
            axis[i].set_xlim(0, series.index.max())
    for i in [0, 3, 6]:
        axis[i].set_ylabel(units)
    for ax in axis[-n_cols:]:
        ax.set_xlabel("Time (days)")
    axis[5].legend(bbox_to_anchor=(1.2, 0.5))
    sns.despine(fig)
    fig.savefig(os.path.join("results", "facs.cell_type.marker_per_patient.{}.svg".format(variable)), bbox_inches="tight")

    # For each patient, plot marker evolution along timepoints
    n_rows = n_cols = int(np.ceil(np.sqrt(len(patients))))
    fig, axis = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 2), sharex=True, sharey=False)
    axis = axis.flatten()
    colors = dict(zip(markers, pallettes))
    for i, patient in enumerate(patients):
        for j, marker in enumerate(markers[~markers.str.contains("CM|EM|TEMRA|PD1|Treg|naive")]):
            series = df2.loc[(patient,), marker]
            axis[i].scatter(
                series.index, series,
                s=25, label=marker, color=colors[marker])
            axis[i].plot(
                series.index,
                series.fillna(method='pad'),
                linestyle="--", label=marker, color=colors[marker], lw=1)
            axis[i].set_xscale('symlog', basex=2)
            axis[i].set_xlim(0, series.index.max())
            if ylog:
                axis[i].set_yscale('symlog', basex=10)
            axis[i].set_title(patient)
    for i in [0, 3, 6]:
        axis[i].set_ylabel(units)
    for ax in axis[-3:]:
        ax.set_xlabel("Time (days)")
    axis[5].legend(bbox_to_anchor=(1.2, 0.5))
    sns.despine(fig)
    fig.savefig(os.path.join("results", "facs.cell_type.patient_per_marker.cell_types.{}.svg".format(variable)), bbox_inches="tight")

    # For each patient, plot marker evolution along timepoints
    n_rows = n_cols = int(np.ceil(np.sqrt(len(patients))))
    fig, axis = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 2), sharex=True, sharey=False)
    axis = axis.flatten()
    colors = dict(zip(markers, pallettes))
    for i, patient in enumerate(patients):
        for j, marker in enumerate(markers[markers.str.contains("CM|EM|TEMRA|PD1|Treg|naive")]):
            series = df2.loc[(patient,), marker]
            axis[i].scatter(
                series.index, series,
                s=25, label=marker, color=colors[marker])
            axis[i].plot(
                series.index,
                series.fillna(method='pad'),
                linestyle="--", label=marker, color=colors[marker], lw=1)
            axis[i].set_xscale('symlog', basex=2)
            axis[i].set_xlim(0, series.index.max())
            if ylog:
                axis[i].set_yscale('symlog', basex=10)
            axis[i].set_title(patient)
    for i in [0, 3, 6]:
        axis[i].set_ylabel(units)
    for ax in axis[-3:]:
        ax.set_xlabel("Time (days)")
    axis[5].legend(bbox_to_anchor=(1.2, 0.5))
    sns.despine(fig)
    fig.savefig(os.path.join("results", "facs.cell_type.patient_per_marker.Tcell_subtypes.{}.svg".format(variable)), bbox_inches="tight")


# Expression
df = pd.read_csv(os.path.join("metadata", "facs.expression.quantification.csv"))
df = df.rename(columns={"Nkcells_CD5": "NKcells_CD5"})

df2 = pd.melt(df, id_vars=['timepoint', 'patient_id'])
df2 = df2.join(df2['variable'].str.split("_").apply(pd.Series))
df2 = df2.rename(columns={0: "cell_type", 1: "marker"})

df2.loc[df2['value'] < 0, 'value'] = 0

df3 = pd.pivot_table(df2, index=['cell_type', 'timepoint', 'patient_id'], columns='marker', values="value")


# Heatmaps
from scipy.stats import zscore
zscore = lambda x, axis=None: (x - x.mean(axis=axis)) / x.std(axis=axis, ddof=0)

# per patient

fig, axis = plt.subplots(1, 3, figsize=(3 * 6, 1 * 6))
sns.heatmap(np.log2(1 + df3.dropna(how="all")), ax=axis[0], cbar_kws={"label": "Expression (log2)"}, square=False)
sns.heatmap(zscore(np.log2(1 + df3.dropna(how="all")), axis=0), ax=axis[1], cbar_kws={"label": "Expression (log2)"}, square=False)
g = zscore(np.log2(1 + df3.dropna(how="all")), axis=0).groupby(level=["cell_type", "timepoint"])
sns.heatmap(g.mean(), ax=axis[2], mask=g.count() < 2, square=True)
for ax in axis:
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize="xx-small")
fig.savefig(os.path.join("results", "facs.expression.heatmap.all_cell_types.all_markers.svg"), dpi=300, bbox_inches="tight")

for cell_type in df3.index.get_level_values('cell_type').unique():
    fig, axis = plt.subplots(1, 3, figsize=(3 * 6, 1 * 6))
    sns.heatmap(np.log2(1 + df3.loc[df3.index.get_level_values('cell_type') == cell_type].dropna(how="all")), ax=axis[0], square=True)
    sns.heatmap(zscore(np.log2(1 + df3.loc[df3.index.get_level_values('cell_type') == cell_type].dropna(how="all")), axis=0), ax=axis[1], square=True)
    g = zscore(np.log2(1 + df3.loc[df3.index.get_level_values('cell_type') == cell_type].dropna(how="all")), axis=0).groupby(level="timepoint")
    sns.heatmap(g.mean(), ax=axis[2], mask=g.count() < 2, square=True)
    for ax in axis:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize="xx-small")
    fig.savefig(os.path.join("results", "facs.expression.heatmap.{}.all_markers.svg".format(cell_type)), dpi=300, bbox_inches="tight")




    # # 
    # q = np.log2(1 + df3.loc[df3.index.get_level_values('cell_type') == cell_type].dropna(how="all"))
    # q = q.loc[~q.isnull().all(1), ~q.isnull().all()]

    # g = sns.clustermap(q, mask=q.isnull())
    # g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    # g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    # g.savefig(os.path.join("results", "facs.expression.heatmap.{}.all_markers.clustermap.svg".format(cell_type)), dpi=300, bbox_inches="tight")


# patient mean
df4 = df3.groupby(level=["cell_type", "timepoint"]).mean()
mask = df3.groupby(level=["cell_type", "timepoint"]).count() < 2

fig, axis = plt.subplots(1, 2, figsize=(2 * 6, 1 * 6))
sns.heatmap(np.log2(1 + df4), ax=axis[0], mask=mask)
axis[0].set_xticklabels(axis[0].get_xticklabels(), rotation=90)
axis[0].set_yticklabels(axis[0].get_yticklabels(), rotation=0, fontsize="xx-small")
sns.heatmap(zscore(np.log2(1 + df4[~mask]), axis=0), ax=axis[1], mask=mask)
axis[1].set_xticklabels(axis[1].get_xticklabels(), rotation=90)
axis[1].set_yticklabels(axis[1].get_yticklabels(), rotation=0, fontsize="xx-small")
fig.savefig(os.path.join("results", "facs.expression.heatmap.all_cell_types.all_markers.patient_mean.svg"), dpi=300, bbox_inches="tight")


n_cell_types = len(df4.index.get_level_values('cell_type').unique())
fig, axis = plt.subplots(2, n_cell_types, figsize=(n_cell_types * 4, 2 * 4), sharey=True, sharex=False)

for i, cell_type in enumerate(df4.index.get_level_values('cell_type').unique()):
    # fig, axis = plt.subplots(1)
    # g = sns.heatmap(df4.loc[df4.index.get_level_values('cell_type') == cell_type], ax=axis)
    # axis.set_xticklabels(axis.get_xticklabels(), rotation=90)
    # axis.set_yticklabels(axis.get_yticklabels(), rotation=0)
    # fig.savefig(os.path.join("results", "facs.expression.heatmap.{}.all_markers.patient_mean.svg".format(cell_type)), dpi=300, bbox_inches="tight")

    # fig, axis = plt.subplots(1)
    # g = sns.heatmap(zscore(df4.loc[df4.index.get_level_values('cell_type') == cell_type], axis=0), ax=axis)
    # axis.set_xticklabels(axis.get_xticklabels(), rotation=90)
    # axis.set_yticklabels(axis.get_yticklabels(), rotation=0)
    # fig.savefig(os.path.join("results", "facs.expression.heatmap.{}.all_markers.patient_mean.zscore.svg".format(cell_type)), dpi=300, bbox_inches="tight")

    # barplot
    tmp_df = pd.melt(df3.loc[df3.index.get_level_values('cell_type') == cell_type].reset_index(), id_vars=["timepoint", "cell_type", "patient_id"])
    tmp_df['value_log'] = np.log2(1 + tmp_df['value'])
    tmp_df_zscore = pd.melt(zscore(df3.loc[df3.index.get_level_values('cell_type') == cell_type], axis=0).reset_index(), id_vars=["timepoint", "cell_type", "patient_id"])
    tmp_df_zscore_log = pd.melt(zscore(df3.loc[df3.index.get_level_values('cell_type') == cell_type].apply(np.log2), axis=0).reset_index(), id_vars=["timepoint", "cell_type", "patient_id"])

    g = sns.barplot(data=tmp_df, x="marker", y="value_log", hue="timepoint", ax=axis[0, i], palette="inferno")
    axis[0, i].set_ylabel("Expression (log2)")
    g = sns.barplot(data=tmp_df_zscore_log, x="marker", y="value", hue="timepoint", ax=axis[1, i], palette="inferno")
    axis[1, i].set_ylabel("Expression (log2 Z-score)")
    axis[0, i].set_title(cell_type)
    axis[1, i].set_title(cell_type)
sns.despine(fig)
fig.savefig(os.path.join("results", "facs.expression.heatmap.all_cell_types.all_markers.patient_mean.barplot.svg".format(cell_type)), dpi=300, bbox_inches="tight")


# Test differences
from scipy.stats import ks_2samp, ttest_ind, mannwhitneyu

results = pd.DataFrame()
for i, cell_type in enumerate(df3.index.get_level_values('cell_type').unique()):
    for timepoint in df3.index.levels[1]:
        for marker in df3.columns:
            a = np.log2(1 + df3.loc[(df3.index.get_level_values("cell_type") == cell_type) & (df3.index.get_level_values("timepoint") == timepoint), marker])
            b = np.log2(1 + df3.loc[(df3.index.get_level_values("cell_type") == cell_type) & (df3.index.get_level_values("timepoint") == 0), marker])

            results = results.append(pd.Series(
                [cell_type, timepoint, marker, np.log2(a.mean() / b.mean()), *ks_2samp(a, b), *ttest_ind(a, b), *mannwhitneyu(a, b)],
                index=["cell_type", "timepoint", "marker", "fold_change", "KS_stat", "KS_p_value", "T_stat", "T_p_value", "MW_stat", "MW_p_value"]), ignore_index=True)

results.to_csv(os.path.join("results", "facs.expression.diff_test.csv"), index=False)
