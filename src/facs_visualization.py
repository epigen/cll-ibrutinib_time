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
df = pd.read_csv("metadata/facs_quantification.csv").set_index(["patient_id", "time"])
patients = df.index.levels[0]
markers = df.columns[~df.columns.str.islower()]
pallettes = sns.color_palette("colorblind") + sns.color_palette("Set1")[::-1]


for variable, units, ylog in [("absolute", "Cells / uL", True), ("percentage", "% live cells", False)]:
    df2 = df[df['variable'] == variable]

    # For each marker, see the abundance of each cell type dependent on time
    df_melted = pd.melt(df2[markers].reset_index(), id_vars=["patient_id", "time"], var_name="cell_type", value_name=variable).rename(columns={variable: units})

    g = sns.factorplot(
        data=df_melted, x="time", y=units, col="cell_type",
        size=2, sharey=False, ci=75, col_wrap=3, legend_out=True, palette="viridis", # color=sns.color_palette("colorblind")[0],
        errwidth=2)
    g.savefig(os.path.join("results", "facs.marker_per_time.{}.line.svg".format(variable)), bbox_inches="tight")

    g = sns.factorplot(
        data=df_melted, x="time", y=units, col="cell_type",
        size=2, sharey=False, ci=75, col_wrap=3, legend_out=True, palette="viridis", # color=sns.color_palette("colorblind")[0],
        kind="bar", errwidth=2)
    g.savefig(os.path.join("results", "facs.marker_per_time.{}.bar.svg".format(variable)), bbox_inches="tight")

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
    for ax in axis[-3:]:
        ax.set_xlabel("time (days)")
    axis[5].legend(bbox_to_anchor=(1.2, 0.5))
    sns.despine(fig)
    fig.savefig(os.path.join("results", "facs.marker_per_patient.{}.svg".format(variable)), bbox_inches="tight")

    # For each patient, plot marker evolution along timepoints
    n_rows = n_cols = int(np.ceil(np.sqrt(len(patients))))
    fig, axis = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 2), sharex=True, sharey=False)
    axis = axis.flatten()
    colors = dict(zip(markers, pallettes))
    for i, patient in enumerate(patients):
        for j, marker in enumerate(markers):
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
        ax.set_xlabel("time (days)")
    axis[5].legend(bbox_to_anchor=(1.2, 0.5))
    sns.despine(fig)
    fig.savefig(os.path.join("results", "facs.patient_per_marker.{}.svg".format(variable)), bbox_inches="tight")


    # Observe evolution regarding timepoint 0
    n_rows = n_cols = int(np.ceil(np.sqrt(len(patients))))
    fig, axis = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 2), sharex=True, sharey=True)
    axis = axis.flatten()
    colors = dict(zip(patients, pallettes))

    for i, marker in enumerate(markers):
        dfm = df2.reset_index().set_index("time").groupby(['patient_id']).apply(lambda x: x.loc[:, marker])

        for j, patient in enumerate(patients):
            series = np.log2(dfm.loc[patient].T / dfm.loc[patient][0].T)

            axis[i].axhline(0, color="black", linestyle="--", alpha=0.6)
            axis[i].scatter(
                series.index, series,
                s=25, label=patient, color=colors[patient])
            axis[i].plot(
                series.index,
                series.fillna(method='pad'),
                linestyle="--", label=patient, color=colors[patient], lw=1)
            axis[i].set_xscale('symlog', basex=2)
            axis[i].set_xlim(0, series.index.max())
            axis[i].set_title(marker)
    for i in [0, 3, 6]:
        axis[i].set_ylabel("log2(fold-chcange over day 0)")
    for ax in axis[-3:]:
        ax.set_xlabel("time (days)")
    axis[5].legend(bbox_to_anchor=(1.2, 0.5))
    sns.despine(fig)
    fig.savefig(os.path.join("results", "facs.marker_per_patient.change.{}.svg".format(variable)), bbox_inches="tight")


    n_rows = n_cols = int(np.ceil(np.sqrt(len(patients))))
    fig, axis = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 2), sharex=True, sharey=True)
    axis = axis.flatten()
    colors = dict(zip(markers, pallettes))

    for j, marker in enumerate(markers):
        dfm = df2.reset_index().set_index("time").groupby(['patient_id']).apply(lambda x: x.loc[:, marker])

        for i, patient in enumerate(patients):
            series = np.log2(dfm.loc[patient].T / dfm.loc[patient][0].T)

            axis[i].axhline(0, color="black", linestyle="--", alpha=0.6)
            axis[i].scatter(
                series.index, series,
                s=25, label=marker, color=colors[marker])
            axis[i].plot(
                series.index,
                series.fillna(method='pad'),
                linestyle="--", label=marker, color=colors[marker], lw=1)
            axis[i].set_xscale('symlog', basex=2)
            axis[i].set_xlim(0, series.index.max())
            axis[i].set_title(patient)
    for i in [0, 3, 6]:
        axis[i].set_ylabel("log2(fold-chcange over day 0)")
    for ax in axis[-3:]:
        ax.set_xlabel("time (days)")
    axis[5].legend(bbox_to_anchor=(1.2, 0.5))
    sns.despine(fig)
    fig.savefig(os.path.join("results", "facs.patient_per_marker.change.{}.svg".format(variable)), bbox_inches="tight")

