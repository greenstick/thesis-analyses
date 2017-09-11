#! /home/users/cordier/.linuxbrew/bin/python3

# Imports & Methods
from operator import itemgetter
from collections import OrderedDict
import operator, hashlib
import vcf as VCF
import argparse
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
import pandas as pd
import numpy as np
import re

def loadJSON (path):
    import json
    with open(path, "r") as data:
        loaded = json.load(data)
    return loaded

if __name__ == "__main__":

    metadata = {

        "datasets" : [
            {
                "name"  : "Set 1",
                "key"   : "set1"
            },
            {
                "name"  : "Set 2",
                "key"   : "set2"
            },
            {
                "name"  : "Set 3",
                "key"   : "set3"
            }
        ],

        "conditions" : [
            {
                "name"  : "BFC",
                "key"   : "bfc"
            },
            {
                "name"  : "Bloocoo",
                "key"   : "bloocoo"
            },
            {
                "name"  : "Lighter",
                "key"   : "lighter"
            },
            {
                "name"  : "No Model",
                "key"   : "nomodel"
            }
        ],
        "subconditions" : [
            {
                "name"  : "No BQSR",
                "key"   : "nobqsr"
            },
            {
                "name"  : "BQSR",
                "key"   : "bqsr"
            }
        ],
        "allranges"  : [("%.2f" % ((i / 2.0) / 10), "%.2f" % (((i + 1) / 2.0 ) / 10)) for i in range(20)],
        "lfvranges"  : [("%.2f" % (i / 100.0), "%.2f" % ((i + 1) / 100.0)) for i in range(20)]

    }

    #
    # Generate Heat Maps
    #

    # Statistics to Generate Heat Maps & Their Colors
    statistics = ["sensitivity", "specificity", "tpr", "tnr", "fpr", "fnr"]
    statcolors = ["#108070", "#108070", "#108070", "#108070", "#ED303C", "#ED303C"]
        

    # Structure Data into Dataframes
    for variantRange in ["lfvranges", "allranges"]:
        conditions = loadJSON("js/json/%s-heatmaps.json" % variantRange)
        matrices = {}
        for statistic, color in zip(statistics, statcolors):
            xLabels = []
            matrices[statistic] = {}
            for condition in conditions:
                xLabels.append(condition["label"])
                matrices[statistic][condition["metakey"]] = []
                column = matrices[statistic][condition["metakey"]]
                rangedAnalyses = sorted(condition["maxAnalyses"], key = lambda d: d['vafrange'][0]) 
                # if "set3-bloocoo-nobqsr" == condition["metakey"]:
                    # print(rangedAnalyses)
                    # exit()
                for rangedAnalysis in rangedAnalyses:
                    value = rangedAnalysis["analysis"][statistic]
                    column.append(value)

            matrix = matrices[statistic]

            #
            # Landscape
            #

            # Coerce to Data Frame
            df = pd.DataFrame(matrix)

            # Colors
            linear_cm = LinearSegmentedColormap.from_list("custom", ["#FFFFFF", color], N = 256, gamma = 1.0)
            colors = cm.get_cmap(linear_cm)

            # Plot it out
            fig, ax = plt.subplots()
            heatmap = ax.pcolor(df, cmap = colors, alpha = 0.8)

            # Format
            fig = plt.gcf()
            fig.set_size_inches(12, 8)

            # turn off the frame
            ax.set_frame_on(False)

            # put the major ticks at the middle of each cell
            ax.set_yticks(np.arange(df.shape[0]) + 0.5, minor = False)
            ax.set_xticks(np.arange(df.shape[1]) + 0.1, minor = False)

            # want a more natural, table-like display
            ax.invert_yaxis()
            ax.xaxis.tick_top()

            # Set Labels
            ax.set_xticklabels(xLabels, minor = False, fontsize = 5)
            ax.set_yticklabels(["%s - %s" % tup for tup in metadata[variantRange]], minor = False, fontsize = 6)

            # Text Styling
            plt.xticks(ha = "left")
            ax.grid(False)

            # Turn off all the ticks
            ax = plt.gca()

            for t in ax.xaxis.get_major_ticks():
                t.tick1On = False
                t.tick2On = False
            for t in ax.yaxis.get_major_ticks():
                t.tick1On = False
                t.tick2On = False
            fig.savefig("output/%s-%s-heatmap-landscape.pdf" % (statistic, variantRange[0:3]), bbox_inches = 'tight', transparent = True)
            # plt.show()
            plt.clf()

            #
            # Portrait
            #

            # Coerce to Data Frame
            df = pd.DataFrame(matrix).T

            # Colors
            linear_cm = LinearSegmentedColormap.from_list("custom", ["#FFFFFF", color], N = 256, gamma = 1.0)
            colors = cm.get_cmap(linear_cm)

            # Plot it out
            fig, ax = plt.subplots()
            heatmap = ax.pcolor(df, cmap = colors, alpha=0.8)

            # Format
            fig = plt.gcf()
            fig.set_size_inches(8, 11)

            # turn off the frame
            ax.set_frame_on(False)

            # put the major ticks at the middle of each cell
            ax.set_yticks(np.arange(df.shape[0]) + .5, minor = False)
            ax.set_xticks(np.arange(df.shape[1]) + .5, minor = False)

            # want a more natural, table-like display
            ax.invert_yaxis()
            ax.xaxis.tick_top()

            # Set Labels
            ax.set_yticklabels(xLabels, minor = False, fontsize = 6)
            ax.set_xticklabels(["%s - %s" % tup for tup in metadata[variantRange]], minor = False, fontsize = 6)

            # Text Styling
            plt.xticks(rotation = 45)
            ax.grid(False)

            # Turn off all the ticks
            ax = plt.gca()

            for t in ax.xaxis.get_major_ticks():
                t.tick1On = False
                t.tick2On = False
            for t in ax.yaxis.get_major_ticks():
                t.tick1On = False
                t.tick2On = False
            fig.savefig("output/%s-%s-heatmap-portrait.pdf" % (statistic, variantRange[0:3]), bbox_inches = 'tight', transparent = True)
            # plt.show()
            plt.clf()
else:
    pass