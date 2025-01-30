
import pandas as pd
import numpy as np
import os 
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

mpl.rcParams['font.size'] = 10
mpl.rcParams['text.color'] = 'k'
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['legend.handletextpad'] = '0.8'
mpl.rcParams['legend.labelspacing'] = '0.4'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['axes.labelsize'] = 10

def make_mut_sigs_plot(path_mut_sigs, ax):
    """
    Makes the classic bar charts showing mut sigs
    """
    sigs = pd.read_csv(path_mut_sigs)
    sigs.columns = ["Base_change", "value"]
    
    # COSMIC colors
    color_dict = {
        "C>A": "#01bcf1",
        "C>G": "#000000",
        "C>T": "#e52a25",
        "T>A": "#cac9c9",
        "T>C": "#a3ce62",
        "T>G": "#ecc6c5"
        }
    
    sigs["Base_change_extracted"] = sigs["Base_change"].str.extract(r'\[(.*?)\]')
    sigs["Color"] = sigs["Base_change_extracted"].map(color_dict)
    
    # Plotting
    ax.bar(sigs.index, sigs["value"], color=sigs["Color"], align='center')
    ax.set_xlim((-1, 96))
    ax.set_xticks(sigs.index)
    ax.set_xticklabels(sigs["Base_change"], rotation = 90, fontsize = 5)
    ax.tick_params(axis='x', bottom=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim((0, 0.15))
    ax.set_yticks([0, 0.05, 0.10, 0.15])
    ax.set_yticklabels(["0", "0.05", "0.10", "0.15"])
    ax.set_ylabel("Enrichment")
    
    return(ax)

def add_annotations(sigs, ax):
    """
    Mainly aes.
    """
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.tick_params(axis='x', bottom=False, labelbottom = False)
    for i, group in sigs.groupby("Base_change_extracted"):
        x_start = group.index.min()
        x_end = group.index.max()
        bar_size = x_end - x_start + 1
        bar_midpoint = (x_end - x_start)/2 + x_start
        bar_color = group["Color"].unique()
        
        ax.bar(bar_midpoint, 1, color = bar_color, width = bar_size)
        
        # For text annotations
        text = group["Base_change_extracted"].unique()[0]
        text_color = "black"
        if text == "C>G":
            text_color = "white"
        ax.text(bar_midpoint, 0.5, text, ha='center', va='center', color=text_color, fontsize = 9)
    return(ax)



fig, ax = plt.subplots(figsize = (8, 3))

gs = gridspec.GridSpec(2, 1, height_ratios=[0.07, 1], hspace = 0, wspace = 0)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1], sharex = ax1)


path_mut_sigs = "/groups/wyattgrp/users/amunzur/pipeline/results/mut_sigs/chip_sigs_product.csv"

ax1 = add_annotations(sigs, ax1)
ax2 = make_mut_sigs_plot(path_mut_sigs, ax2)
fig.tight_layout()
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_mut_sigs.png")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_mut_sigs.pdf")



