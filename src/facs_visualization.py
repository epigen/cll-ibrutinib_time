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
df = pd.read_csv("metadata/facs_quantification.csv").set_index(["patient_id", "timepoint"])
patients = df.index.levels[0]
markers = df.columns
pallettes = sns.color_palette("colorblind") + sns.color_palette("Set1")[::-1]
fold_reduct_threshold = -0.5


# For each marker, plot patient evolution along timepoints
n_rows = n_cols = int(np.ceil(np.sqrt(len(markers))))
fig, axis = plt.subplots(n_rows, n_cols, figsize=(n_cols * 6, n_rows * 4), sharex=True, sharey=False)
axis = axis.flatten()
colors = dict(zip(patients, pallettes))
for i, marker in enumerate(markers):
    for j, patient in enumerate(patients):
        series = df.loc[patient, marker]
        axis[i].scatter(
            series.index, series,
            s=25, label=patient, color=colors[patient])
        axis[i].plot(
            series.index,
            series.fillna(method='pad'),
            linestyle="--", label=patient, color=colors[patient])
        axis[i].set_xscale('symlog', basex=2)
        axis[i].set_title(marker)
        # axis[i].set_xticklabels(l.set_text(str(l.get_position()[0] ** 2)) for l in axis[i].get_xticklabels())
        axis[i].set_xlim(0, series.index.max())
for i in [0, 3, 6]:
    axis[i].set_ylabel("% live cells")
for ax in axis[-3:]:
    ax.set_xlabel("time (days)")
axis[5].legend(bbox_to_anchor=(1.2, 0.5))
sns.despine(fig)
fig.savefig(os.path.join("results", "facs.marker_per_patient.svg"), bbox_inches="tight")


# For each patient, plot marker evolution along timepoints
n_rows = n_cols = int(np.ceil(np.sqrt(len(patients))))
fig, axis = plt.subplots(n_rows, n_cols, figsize=(n_cols * 6, n_rows * 4), sharex=True, sharey=False)
axis = axis.flatten()
colors = dict(zip(markers, pallettes))
for i, patient in enumerate(patients):
    for j, marker in enumerate(df.columns):
        series = df.loc[patient, marker]
        axis[i].scatter(
            series.index, series,
            s=25, label=marker, color=colors[marker])
        axis[i].plot(
            series.index,
            series.fillna(method='pad'),
            linestyle="--", label=marker, color=colors[marker])
        axis[i].set_xscale('symlog', basex=2)
        axis[i].set_xlim(0, series.index.max())
        axis[i].set_title(patient)
for i in [0, 3, 6]:
    axis[i].set_ylabel("% live cells")
for ax in axis[-3:]:
    ax.set_xlabel("time (days)")
axis[5].legend(bbox_to_anchor=(1.2, 0.5))
sns.despine(fig)
fig.savefig(os.path.join("results", "facs.patient_per_marker.svg"), bbox_inches="tight")


# Observe evolution regarding timepoint 0
n_rows = n_cols = int(np.ceil(np.sqrt(len(patients))))
fig, axis = plt.subplots(n_rows, n_cols, figsize=(n_cols * 6, n_rows * 4), sharex=True, sharey=True)
axis = axis.flatten()
colors = dict(zip(patients, pallettes))

for i, marker in enumerate(df.columns):
    dfm = df.reset_index().set_index("timepoint").groupby(['patient_id']).apply(lambda x: x.loc[:, marker])
    dfm = np.log2(dfm.T / dfm[0].T)

    for j, patient in enumerate(patients):
        series = dfm[patient]
        axis[i].axhline(0, color="black", linestyle="--", alpha=0.6)
        axis[i].scatter(
            series.index, series,
            s=25, label=patient, color=colors[patient])
        axis[i].plot(
            series.index,
            series.fillna(method='pad'),
            linestyle="--", label=patient, color=colors[patient])
        axis[i].set_xscale('symlog', basex=2)
        axis[i].set_xlim(0, series.index.max())
        axis[i].set_title(marker)
for i in [0, 3, 6]:
    axis[i].set_ylabel("log2(fold-chcange over day 0)")
for ax in axis[-3:]:
    ax.set_xlabel("time (days)")
axis[5].legend(bbox_to_anchor=(1.2, 0.5))
sns.despine(fig)
fig.savefig(os.path.join("results", "facs.marker_per_patient.change.svg"), bbox_inches="tight")


n_rows = n_cols = int(np.ceil(np.sqrt(len(patients))))
fig, axis = plt.subplots(n_rows, n_cols, figsize=(n_cols * 6, n_rows * 4), sharex=True, sharey=True)
axis = axis.flatten()
colors = dict(zip(markers, pallettes))

for j, marker in enumerate(df.columns):
    dfm = df.reset_index().set_index("timepoint").groupby(['patient_id']).apply(lambda x: x.loc[:, marker])
    dfm = np.log2(dfm.T / dfm[0].T)

    for i, patient in enumerate(patients):
        series = dfm[patient]
        axis[i].axhline(0, color="black", linestyle="--", alpha=0.6)
        axis[i].scatter(
            series.index, series,
            s=25, label=marker, color=colors[marker])
        axis[i].plot(
            series.index,
            series.fillna(method='pad'),
            linestyle="--", label=marker, color=colors[marker])
        axis[i].set_xscale('symlog', basex=2)
        axis[i].set_xlim(0, series.index.max())
        axis[i].set_title(patient)
for i in [0, 3, 6]:
    axis[i].set_ylabel("log2(fold-chcange over day 0)")
for ax in axis[-3:]:
    ax.set_xlabel("time (days)")
axis[5].legend(bbox_to_anchor=(1.2, 0.5))
sns.despine(fig)
fig.savefig(os.path.join("results", "facs.patient_per_marker.change.svg"), bbox_inches="tight")


# Stratify patients
# get patients which reach X reduction of CLL at the last timepoint (either 120 or 240 days)
dfm = df.reset_index().set_index("timepoint").groupby(['patient_id']).apply(lambda x: x.loc[:, "CLL"])
dfm = np.log2(dfm.T / dfm[0].T)

resp = dfm.columns[dfm.ix[["120d", "240d"]].min() < fold_reduct_threshold]

facs_classification = os.path.join("metadata", "facs.response_classification.txt")
with open(facs_classification, "w") as handle:
    handle.writelines("\n".join(resp))
