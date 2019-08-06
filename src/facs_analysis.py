import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import dabest


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
df = df[df.index.get_level_values('patient_id') != "CLL3"]

for variable, units, ylog in [("absolute", "Cells / uL", True), ("percentage", "% live cells", False)]:
    df2 = df[df['value_type'] == variable]

    # For each marker, see the abundance of each cell type dependent on time
    df_melted = pd.melt(df2[markers].reset_index(), id_vars=["patient_id", "timepoint"], var_name="cell_type", value_name=variable).rename(columns={variable: units})

    # Test against timepoint 0
    import scipy
    res = pd.DataFrame()
    for cell_type in df_melted['cell_type'].drop_duplicates():
        for timepoint in df_melted.loc[df_melted['cell_type'] == cell_type, "timepoint"].drop_duplicates():
            if (timepoint == 0) or pd.isnull(timepoint):
                continue
            t = df_melted.loc[(df_melted['cell_type'] == cell_type) & (df_melted['timepoint'] == timepoint), units]
            t.index = range(t.shape[0])
            t.name = timepoint
            t0 = df_melted.loc[(df_melted['cell_type'] == cell_type) & (df_melted['timepoint'] == 0), units]
            t0.index = range(t0.shape[0])
            t0.name = 0
            t2 = t0.to_frame().join(t).dropna()
            ts, tp = scipy.stats.ttest_ind(t.dropna(), t0.dropna())
            pts, ptp = scipy.stats.ttest_rel(t2.loc[:, timepoint], t2.loc[:, 0])
            ws, wp = scipy.stats.mannwhitneyu(t.dropna(), t0.dropna())
            res = res.append(pd.Series(
                [cell_type, timepoint, t.dropna().shape[0], t0.dropna().shape[0], t2.shape[0],
                 t.mean(), t0.mean(), ts, tp, pts, ptp, ws, wp]),
                ignore_index=True)
    res.columns = ['cell_type', 'timepoint', 'n_t', 'n_0', 'pairs', 't_mean', 't0_mean', 'T_stat',
                   'T_p_value', 'paired_T_stat', 'paired_T_p_value', 'MW_stat', 'MW_p_value']
    res.to_csv(os.path.join("results", "facs.cell_type.marker_per_time.{}.statistics.csv".format(variable)), index=False)

    df_melted = df_melted.dropna()
    df_melted['timepoint'] = df_melted['timepoint'].astype(int)

    g = sns.catplot(
        data=df_melted, x="timepoint", y=units, col="cell_type",
        height=2, sharey=False, ci=75, col_wrap=3, legend_out=True, palette="inferno",
        kind="bar", errwidth=2)
    g.savefig(os.path.join("results", "facs.cell_type.marker_per_time.{}.bar.svg".format(variable)), bbox_inches="tight")
    for ax in g.axes:
        ax.set_yscale("log", basey=10)
    g.savefig(os.path.join("results", "facs.cell_type.marker_per_time.{}.bar.log.svg".format(variable)), bbox_inches="tight")

    g = sns.catplot(
        data=df_melted, x="timepoint", y=units, col="cell_type",
        height=2, sharey=False, ci=75, col_wrap=3, legend_out=True, palette="inferno",
        kind="point", markers="^", errwidth=2)
    g.map(sns.swarmplot, "timepoint", units, palette="inferno")
    g.savefig(os.path.join("results", "facs.cell_type.marker_per_time.{}.line.svg".format(variable)), bbox_inches="tight")
    for ax in g.axes:
        ax.set_yscale("log", basey=10)
    g.savefig(os.path.join("results", "facs.cell_type.marker_per_time.{}.line.log.svg".format(variable)), bbox_inches="tight")

    g = sns.catplot(
        data=df_melted, x="timepoint", y=units, col="cell_type",
        height=2, sharey=False, col_wrap=3, legend_out=True, palette="inferno",
        kind="swarm")
    g.savefig(os.path.join("results", "facs.cell_type.marker_per_time.{}.swarm.svg".format(variable)), bbox_inches="tight")
    for ax in g.axes:
        ax.set_yscale("log", basey=10)
    g.savefig(os.path.join("results", "facs.cell_type.marker_per_time.{}.swarm.log.svg".format(variable)), bbox_inches="tight")

    for cell_type in df_melted['cell_type'].unique():
        df3 = df_melted[df_melted['cell_type'] == cell_type]
        x_melted = pd.pivot_table(df3, columns='timepoint', index='patient_id')
        x_melted.columns = x_melted.columns.droplevel(0).astype(int).astype(str)
        fig, stats = dabest.plot(
            x_melted.reset_index(),
            idx=tuple(x_melted.columns),
            color_col='patient_id',
            paired=True, show_pairs=True,
            # custom_palette=sns.color_palette("Set1", n_colors=8, desat=.5)
            )
        fig.axes[0].set_ylim((0, 100))
        fig.savefig(os.path.join("results", "facs.cell_type.marker_per_time.{}.{}.dabest.all_timepoints.svg"
                    .format(variable, cell_type)), bbox_inches="tight")
        stats.to_csv(os.path.join("results", "facs.cell_type.marker_per_time.{}.{}.dabest.all_timepoints.csv"
                     .format(variable, cell_type)), index=False)

    # For each marker, plot patient evolution along timepoints
    n_rows = n_cols = int(np.ceil(np.sqrt(len(markers))))
    fig, axis = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 2), sharex=True, sharey=False)
    axis = axis.flatten()
    colors = dict(zip(patients, pallettes))
    for i, marker in enumerate(markers):
        for j, patient in enumerate(patients):
            series = df2.loc[(patient,), marker].dropna()
            if series.empty:
                continue
            axis[i].scatter(
                series.index, series,
                s=25, label=patient, color=colors[patient])
            axis[i].plot(
                series.index,
                series.fillna(method='pad'),
                linestyle="--", label=patient, color=colors[patient], lw=1)
            axis[i].set_xscale('symlog', basex=2)
            # if ylog:
            #     axis[i].set_yscale('symlog', basex=10)
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

    # g = sns.catplot(
    #     data=df_melted, x="timepoint", y=units, col="patient_id", hue="cell_type",
    #     height=2, sharey=False, sharex=False, ci=75, col_wrap=3, legend_out=True, dodge=False,
    #     kind="point", errwidth=2)
    # for ax in g.axes:
    #     ax.set_yscale("log", basey=10)
    # g.savefig(os.path.join("results", "facs.cell_type.marker_per_time.{}.bar.svg".format(variable)), bbox_inches="tight")

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


# # CD4/CD8 ratios
fig, axis = plt.subplots(1, 2, figsize=(3 * 2, 3), tight_layout=True)
for i, (variable, units, log) in enumerate([
        ("percentage", "% live cells", False),
        ("absolute", "Cells / uL", True)]):

    df2 = df.loc[df['value_type'] == variable, :]
    cs = axis[i].scatter(
        df2['CD4+ T cells'], df2['CD8+ T cells'],
        c=np.log2(df2.index.get_level_values("timepoint")), cmap="inferno")
    vmin, vmax = (
        min(df2['CD4+ T cells'].min(), df2['CD4+ T cells'].min()),
        max(df2['CD4+ T cells'].max(), df2['CD4+ T cells'].max()))
    axis[i].plot((vmin, vmax), (vmin, vmax), linestyle="--", color="grey")
    axis[i].set_xlabel("CD4+ ({})".format(units))
    axis[i].set_ylabel("CD8+ ({})".format(units))
    if log:
        axis[i].set_xscale("log")
        axis[i].set_yscale("log")
fig.savefig(
    os.path.join(
        "results",
        "facs.CD4_CD8_ratio.all_patients_timepoints.scatter.svg"),
    bbox_inches="tight")


tps = df.reset_index()['timepoint'].dropna().unique()
fig, axis = plt.subplots(2, len(tps), figsize=(3 * len(tps), 3 * 2), tight_layout=True)
for i, (variable, units, log) in enumerate([
        ("percentage", "% live cells", False),
        ("absolute", "Cells / uL", True)]):

    for j, timepoint in enumerate(tps):
        df2 = df.loc[
            (df.index.get_level_values("timepoint") == timepoint) &
            (df['value_type'] == variable), :]
        cs = axis[i, j].scatter(
            df2['CD4+ T cells'], df2['CD8+ T cells'])
        vmin, vmax = (
            min(df2['CD4+ T cells'].min(), df2['CD4+ T cells'].min()),
            max(df2['CD4+ T cells'].max(), df2['CD4+ T cells'].max()))
        axis[i, j].plot((vmin, vmax), (vmin, vmax), linestyle="--", color="grey")
        axis[i, j].set_xlabel("CD4+ ({})".format(units))
        axis[i, j].set_ylabel("CD8+ ({})".format(units))
        if log:
            axis[i, j].set_xscale("log")
            axis[i, j].set_yscale("log")
        axis[i, j].set_title("Day {:.0f}".format(timepoint))
fig.savefig(
    os.path.join(
        "results",
        "facs.CD4_CD8_ratio.all_patients.per_timepoint.scatter.svg"),
    bbox_inches="tight")

ratio = (
    df.reset_index()
    .groupby(['value_type', 'timepoint'])
    .apply(lambda x: pd.Series(x['CD8+ T cells'] / x['CD4+ T cells'], name="ratio")))
ratio = np.log2(ratio)
fig, axis = plt.subplots(1, 1, figsize=(3, 3), tight_layout=True)
sns.swarmplot(data=ratio.reset_index(), x="timepoint", y='ratio', hue="value_type", ax=axis)
sns.pointplot(data=ratio.reset_index(), x="timepoint", y='ratio', hue="value_type", ax=axis)
axis.axhline(0, linestyle="--", color="grey")
axis.set_ylabel("log2(CD8+ T cells / CD4+ T cells)")
fig.savefig(
    os.path.join(
        "results",
        "facs.CD4_CD8_ratio.ratios.swarmplot.svg"),
    bbox_inches="tight")


# Expression
df = pd.read_csv(os.path.join("metadata", "facs.expression.quantification.csv"))
df = df.loc[df['patient_id'] != "KI", :]
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
sns.heatmap(np.log2(1 + df3.dropna(how="all")), ax=axis[0], cbar_kws={"label": "Expression (log2)"}, square=False, yticklabels=True)
sns.heatmap(zscore(np.log2(1 + df3.dropna(how="all")), axis=0), ax=axis[1], cbar_kws={"label": "Expression (log2)"}, square=False, cmap="coolwarm", yticklabels=True, center=0)
g = zscore(np.log2(1 + df3.dropna(how="all")), axis=0).groupby(level=["cell_type", "timepoint"])
sns.heatmap(g.mean(), ax=axis[2], mask=g.count() < 2, square=True, cmap="coolwarm", yticklabels=True, center=0)
for ax in axis:
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize="xx-small")
fig.savefig(os.path.join("results", "facs.expression.heatmap.all_cell_types.all_markers.svg"), dpi=300, bbox_inches="tight")

for cell_type in df3.index.get_level_values('cell_type').unique():
    fig, axis = plt.subplots(1, 4, figsize=(4 * 6, 1 * 6))
    q = np.log2(1 + df3.loc[df3.index.get_level_values('cell_type') == cell_type].dropna(how="all")).sort_index(level=['patient_id', 'timepoint']).drop(['PD1', 'PD1%'], axis=1)
    sns.heatmap(q, ax=axis[0], square=True, yticklabels=True)
    sns.heatmap(zscore(q, axis=0), ax=axis[1], square=True, cmap="coolwarm", yticklabels=True, center=0)
    g = q.groupby(level='timepoint')
    sns.heatmap(g.mean(), ax=axis[2], mask=g.count() < 2, square=True, yticklabels=True)
    g = zscore(q, axis=0).groupby(level="timepoint")
    sns.heatmap(g.mean(), ax=axis[3], mask=g.count() < 2, square=True, cmap="coolwarm", yticklabels=True, center=0)
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
from scipy.stats import ks_2samp, ttest_ind, mannwhitneyu, ttest_rel

results = pd.DataFrame()
for i, cell_type in enumerate(df3.index.get_level_values('cell_type').unique()):
    for timepoint in df3.index.levels[1]:
        if timepoint == 0:
            continue
        for marker in df3.columns:
            a = np.log2(1 + df3.loc[
                (df3.index.get_level_values("cell_type") == cell_type) &
                (df3.index.get_level_values("timepoint") == timepoint), marker]).dropna()
            b = np.log2(1 + df3.loc[
                (df3.index.get_level_values("cell_type") == cell_type) &
                (df3.index.get_level_values("timepoint") == 0), marker]).dropna()

            if a.empty and b.empty:
                continue
            a.name = timepoint
            b.name = 0
            j = a.reset_index(drop=True, level=[0, 1]).to_frame().join(b.reset_index(drop=True, level=[0, 1])).dropna()

            results = results.append(pd.Series(
                [cell_type, timepoint, marker, a.dropna().shape[0], b.dropna().shape[0], a.mean(), b.mean(), np.log2(a.mean() / b.mean()),
                 *ks_2samp(a, b), *ttest_rel(j[timepoint], j[0]), *ttest_ind(a, b), *mannwhitneyu(a, b)]), ignore_index=True)

results.columns = ["cell_type", "timepoint", "marker", "N_T", "N_0", "T_mean", "0_mean", "fold_change",
                   "KS_stat", "KS_p_value", "pairedT_stat", "pairedT_p_value", "T_stat", "T_p_value", "MW_stat", "MW_p_value"]
results.to_csv(os.path.join("results", "facs.expression.diff_test.csv"), index=False)


# Observe change in CLL counts from day 0 to day 30
# # Setup
df = pd.read_csv("metadata/facs.cell_type.quantification.csv").set_index(["patient_id", "timepoint"])
df = df[(df.index.get_level_values('timepoint') != 2) & (df.index.get_level_values('patient_id') != "CLL3")]
variable, units = "percentage", "% live cells"
df2 = df[df['value_type'] == variable]
df_melted = pd.melt(df2[markers].reset_index(), id_vars=["patient_id", "timepoint"], var_name="cell_type", value_name=variable).rename(columns={variable: units})
df2 = df_melted.pivot_table(index=["patient_id", "cell_type"], columns="timepoint", values=units)
