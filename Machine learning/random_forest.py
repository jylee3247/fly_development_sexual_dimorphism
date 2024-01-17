#!/usr/bin/env python
# coding: utf-8
##### Gut commensal microbe regulated by host Imd signaling  pathway accelerates sexual maturation in female fruit fly #####

# Random forest classifier using scikit-learn version 1.1.3 (https://scikit-learn.org/stable/index.html).
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import argparse
from sklearn.metrics import auc, RocCurveDisplay, ConfusionMatrixDisplay
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.ensemble import RandomForestClassifier

def load_data(X_path, y_path):
    """
    Load data table and metadata.
    """
    X_df = pd.read_csv(X_path, sep = "\t", index_col = 0)
    y_df = pd.read_csv(y_path, sep = "\t", index_col = 0)

    if not all(X_df.index == y_df.index):
        X_df = X_df.reindex(index = y_df.index)
        if not all(X_df.index == y_df.index):
            raise Exception('Failed to match index')

    return X_df, y_df

def run_randomforest(X_df, y_df, out_dir, cv_splits, cv_repeats, n_jobs, random_state):
    """
    Run Random Forest classifier with repated cross-validation.
    Save feature importance values.
    Draw ROC curve, confusion matrix.
    """
    
    os.makedirs(out_dir, exist_ok=True)
    X = X_df.to_numpy()
    y = y_df.to_numpy().ravel()

    clf = RandomForestClassifier(n_jobs = n_jobs,  random_state = random_state)
    cv = RepeatedStratifiedKFold(n_splits = cv_splits, n_repeats = cv_repeats)
    cv_num = cv_splits*cv_repeats

    print("Current settings:")
    print(f"\tRandomForest:\n\t\t- n_jobs:\t{n_jobs}")
    print(f"\tCV:\n\t\t- n_splits:\t{cv_splits}\n\t\t- n_repeates:\t{cv_repeats}")

    # ROC curve results
    aucs = []
    tprs = []
    mean_fpr = np.linspace(0, 1, 100)
    # Confusion matrix results
    confmat_res = []
    confmat_norm_res = []
    # Feature importance (Mean decrease in impurity[MDI])
    feature_importance = np.zeros([X.shape[1], cv_num])

    # Run CV trials
    for i, (train, test) in enumerate(cv.split(X, y)):
        # Training classifier
        clf.fit(X[train], y[train])
        
        # Get feature importance (Mean decrease in impurity, such as Gini)
        feature_importance[:,i] = clf.feature_importances_

        # Get ROC curve values
        viz = RocCurveDisplay.from_estimator(clf, X[test], y[test], name=f"ROC fold {i}")
        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)

        # Get Confusion matrix values
        confmat_viz = ConfusionMatrixDisplay.from_estimator(clf, X[test], y[test])
        cur_matrix = confmat_viz.confusion_matrix
        
        # Normalize confusion matrix
        cur_norm_matrix = np.zeros([2, 2])
        for i in range(cur_matrix.shape[0]):
            cur_norm_matrix[i] = cur_matrix[i] / cur_matrix.sum(axis = 1)[i]
        cur_res = cur_matrix.ravel().tolist()
        cur_norm_res = cur_norm_matrix.ravel().tolist()
        confmat_res.append(cur_res)
        confmat_norm_res.append(cur_norm_res)

    # Save importance
    importance = pd.DataFrame(feature_importance)
    importance.to_csv(os.path.join(out_dir, "RF_importance.tsv"), sep='\t')

    # Save ROC curve AUC values
    aucs_series = pd.Series(aucs, name = "AUC_values")
    aucs_series.to_csv(os.path.join(out_dir, "RF_aucs"), sep = "\t")

    # Save Confusion matrix values
    df_confmat = pd.DataFrame(confmat_res)
    df_confmat_norm = pd.DataFrame(confmat_norm_res)

    df_confmat.to_csv(os.path.join(out_dir, "Confusion_matrix.tsv"), sep = "\t")
    df_confmat_norm.to_csv(os.path.join(out_dir, "Confusion_matrix_normalized.tsv"), sep = "\t")

    # Draw ROC curves
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)

    ax.plot(mean_fpr, mean_tpr, color="b", label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc), lw=2, alpha=0.8)
    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color="grey", alpha=0.2, label=r"$\pm$ 1 std. dev.")
    ax.set(xlim=[0, 1], ylim=[0, 1], title=f"ROC curve", xlabel = "False positive rate", ylabel = "True positive rate")
    ax.legend(loc="lower right")

    # plt.show()
    plt.savefig(os.path.join(out_dir, "ROC_curve.svg"), dpi=600)

    # Draw confusion matrix
    first_label = confmat_viz.display_labels[0]
    second_label = confmat_viz.display_labels[1]
    mean_confmat = df_confmat_norm.mean(axis = 0)
    first_dim = mean_confmat[0:2]
    second_dim = mean_confmat[2:4]
    input_mean_confmat = np.array([first_dim, second_dim])
    conf_disp = ConfusionMatrixDisplay(input_mean_confmat, display_labels=(first_label, second_label)).plot(cmap = plt.cm.Blues)
    conf_disp.plot(cmap=plt.cm.Blues)
    plt.title(f"Confusion matrix")
    plt.savefig(os.path.join(out_dir, "Confusion_matrix.svg"), dpi = 600)

# Arguments
parser = argparse.ArgumentParser(description="Run Random Forest Classifier.")
parser.add_argument("-i", "--input_table", type=str, required=True, help="Path to the feature table.")
parser.add_argument("-m", "--metadata", type=str, required=True, help="Path to the metadata.")
parser.add_argument("-o", "--outdir", type=str, required=True, help="Output directory path.")
parser.add_argument("-cv", "--cv_splits", type=int, default=5, help="Number of splits for cross-validation.")
parser.add_argument("-r", "--cv_repeats", type=int, default=1, help="Number of repeats for cross-validation.")
parser.add_argument("-t", "--n_jobs", type=int, default=3, help="Number of jobs to run in parallel.")
parser.add_argument("-rs", "--random_state", type=int, default=111, help="Random state for reproducibility.")
args = parser.parse_args()

# Load data
X_df, y_df = load_data(args.input_table, args.metadata)

# Run Random Forest Classifier
run_randomforest(X_df, y_df, args.outdir, args.cv_splits, args.cv_repeats, args.n_jobs, args.random_state)