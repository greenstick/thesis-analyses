#! /home/users/cordier/.linuxbrew/bin/python3

# Imports & Methods
from collections import OrderedDict
from collections import defaultdict
from operator import itemgetter
from statistics import mean, median, mode, stdev, variance
from sklearn import metrics
import operator, argparse, matplotlib, hashlib
import vcf as VCF

# Allows Matplotlib to Render on Machine With No Display
matplotlib.use('Agg')
# Set matplotlib style
matplotlib.style.use('ggplot')

#
# Data Structuring
#

# Variant Hashing Function
def hashVariants (records, props, algorithm = "md5", collisionCheck = False):
    """
    
    Hash VCF Records and Check for Collisions.
    Accepts Reader object (records) and an algorithm string (algorithm).
    
    @Params:
        PARAMETER    USAGE      DEFAULT   DESCRIPTION
        records    - Required -         - Iterable records object
        props      - Required -         - Properties of record in records iterable
        algorithm  - Optional - "md5"   - String of hashing algorithm 
    @Returns
        hashMap    - OrderdDict
        collisions - Boolean
        
    """
    hashMap, count, collisions = OrderedDict(), 0, False
    for record in records:
        variantStr = "".join([str(getattr(record, prop)) for prop in props])
        variantHash = getattr(hashlib, algorithm)(variantStr.encode()).hexdigest()
        hashMap[variantHash] = record
        count += 1
    if collisionCheck:
        assert (len(hashMap.keys()) == count), print("Warning: Collision occurred (%s). Try using a difference hashing function.\nAvailable algorithms: %s\n" % (algorithm, ", ".join(hashlib.algorithms_available)))
    return hashMap


# Generate Analysis From Variant Hashmaps
def generateAnalysis (conditionHashMap, truthSet, metricKey = "TLOD", operation = "<", threshold = 100):
    
    # Map Operators
    operators = {
        '>': operator.gt,
        '<': operator.lt,
        '>=': operator.ge,
        '<=': operator.le,
        '=': operator.eq
    }

    # Sort Variants into TP/TN/FP/FN Sets
    tpSet, fpSet, tnSet, fnSet = set(), set(), set(), set()
    for hashKey, record in conditionHashMap.items():
        if operators[operation](float(record.INFO[metricKey]), threshold):
            if hashKey in truthSet:
                fnSet.add(hashKey)
            else:
                tnSet.add(hashKey)
        else:
            if hashKey in truthSet:
                tpSet.add(hashKey)
            else:
                fpSet.add(hashKey)

    # Analysis Object
    analysis = {}

    # True Positives
    tp = len(tpSet)
    analysis["tp"] = tp
    del tpSet
    # True Negatives
    tn = len(tnSet)
    analysis["tn"] = tn
    del tnSet
    # False Positives
    fp = len(fpSet)
    analysis["fp"] = fp
    del fpSet
    # False Negatives
    fn = len(fnSet)
    analysis["fn"] = fn
    del fnSet

    # Statistics (Some Computations Replicated for Naming)
    analysis["tp"] = tp
    analysis["tn"] = tn
    analysis["fp"] = fp
    analysis["fn"] = fn
    # Ternary Operations to Avoid Function Limits for Statistics (due to 0 values)
    analysis["sensitivity"]   = tp / (tp + fn) if tp > 0 else 0.0
    analysis["specificity"]   = tn / (tn + fp) if tn > 0 else 0.0
    analysis["precision"]     = tp / (tp + fp) if tp > 0 else 0.0
    analysis["recall"]        = tp / (tp + fn) if tp > 0 else 0.0
    analysis["tpr"]           = tp / (tp + fn) if tp > 0 else 0.0
    analysis["tnr"]           = tn / (tn + fp) if tn > 0 else 0.0
    analysis["fpr"]           = fp / (fp + tn) if fp > 0 else 0.0
    analysis["fnr"]           = fn / (fn + tp) if fn > 0 else 0.0
    analysis["accuracy"]      = (tp + tn) / (tp + tn + fp + fn) if (tp + tn) > 0 else 0.0
    analysis["f1score"]       = (2 * tp) / ((2 * tp) + (fp + fn)) if tp > 0 else 0.0
    analysis["fdr"]           = fp / (fp + tp) if fp > 0 else 0.0
    analysis["for"]           = fn / (fn + tn) if fn > 0 else 0.0
    analysis["ppv"]           = tp / (tp + fp) if tp > 0 else 0.0
    analysis["npv"]           = tn / (tn + fn) if tn > 0 else 0.0
    # Numerator / Denominator for Diagnostic Odds Ratio
    dorNum                    = tp / fp if fp > 0 else 1.0
    dorDen                    = tn / tn if tn > 0 else 1.0
    analysis["diagnosticOR"]  = dorNum / dorDen
    analysis["n"]             = tp + fp + tn + fn
    return analysis

# Compute Lots of Statistics n' Stuff
def generateAnalyses (conditionHashMap, truthHashMap, metricKey = "TLOD", thresholds = 100, skip = 1):
    # Setup Data Structures
    analyses   = OrderedDict()
    # Use Thresholds to Classify Positive or Negative 
    for threshold in range(0, thresholds, skip):
        analyses[threshold] = generateAnalysis(conditionHashMap, truthHashMap, metricKey = metricKey, threshold = threshold)
    return analyses

#
# Variant HashMap Subsetting Functions
#


def subsetHashMapByKeys (hashMap, keys, isfound = True):
    subset = {}
    for md5, record in hashMap.items():
        if isfound:
            if md5 in keys:
                subset[md5] = record
        else:
            if md5 not in keys:
                subset[md5] = record
    return subset

# Get Info Data From Records in Records HashMap
def subsetHashMapByInfoThreshold (hashMap, threshold, operation = ">", key = "TLOD"):
    subset, excluded = {}, {}
    operators = {
        '>': operator.gt,
        '<': operator.lt,
        '>=': operator.ge,
        '<=': operator.le,
        '=': operator.eq
    }
    for md5, record in hashMap.items():
        if key in record.INFO.keys() and operators[operation](record.INFO[key], threshold):
            subset[md5] = record
        else:
            excluded[md5] = record
    return subset, excluded

# Get Format Data From Records in Records HashMap
def subsetHashMapByFormatThreshold (hashMap, threshold, operation = ">", key = "AF", sampleIdx = 0):
    subset, excluded = {}, {}
    operators = {
        '>': operator.gt,
        '<': operator.lt,
        '>=': operator.ge,
        '<=': operator.le,
        '=': operator.eq
    }
    for md5, record in hashMap.items():
        sample = record.samples[sampleIdx]
        try:
            if type(sample[key]) is list and sum(operators[operation](sample[key], threshold)):
                subset[md5] = record
            elif operators[operation](sample[key], threshold):
                subset[md5] = record
            else: 
                excluded[md5] = record
        except:
            excluded[md5] = record
    return subset, excluded

# Get Info Data From Records in Records HashMap
def subsetHashMapByInfoRange (hashMap, rangeTup, key = "TLOD"):
    subset, excluded = {}, {}
    for md5, record in hashMap.items():
        if key in record.INFO.keys() and record.INFO[key] > rangeTup[0] and record.INFO[key] <= rangeTup[1]:
            subset[md5] = record
        else:
            excluded[md5] = record
    return subset, excluded

# Get Format Data From Records in Records HashMap
def subsetHashMapByFormatRange (hashMap, rangeTup, key = "AF", sampleIdx = 0):
    subset, excluded = {}, {}
    for md5, record in hashMap.items():
        sample = record.samples[sampleIdx]
        if sample[key]:
            if type(sample[key]) is list:
                add = sum(sample[key])
                if rangeTup[0] == 0:
                    subset[md5] = record
                elif add > rangeTup[0] and add <= rangeTup[1]:
                    subset[md5] = record
                else:
                    excluded[md5] = record
            else:
                if sample[key] > rangeTup[0] and sample[key] <= rangeTup[1]:
                    subset[md5] = record
                else:
                    excluded[md5] = record
        else:
            excluded[md5] = record
    return subset, excluded

#
# VCF Parsing
#

# Get Info Data From Records in Records HashMap
def getInfoData (hashMap, key = "TLOD"):
    data, excluded = [], {}
    for md5, record in hashMap.items():
        try:
            data.append(record.INFO[key])
        except:
            excluded[md5] = record
    return data, excluded

# Get Format Data From Records in Records HashMap
def getFormatData (hashMap, key = "AF", sampleIdx = 0):
    data, excluded = [], {}
    types = set()
    for md5, record in hashMap.items():
        try:
            sample = record.samples[sampleIdx]
            data.append(sample[key])
        except:
            excluded[md5] = record
    return data, excluded

# Retrieve Passing Variants (i.e. Those Which Have No Failing Filters)
def passingVariants (hashMap):
    passing = {}
    for md5, record in hashMap.items():
        if len(record.FILTER) == 0:
            passing[md5] = record
    return passing

#
# Statistics Computations & Reports
#

# Compute Quartiles of a Vector
def quartiles(vector):
    from statistics import median
    sort = sorted(vector)
    mid = len(sort) // 2
    if (len(sort) % 2 == 0):
        lowerQ = median(sort[:mid])
        upperQ = median(sort[mid:])
    else:
        lowerQ = median(sort[:mid])
        upperQ = median(sort[mid+1:])
    return lowerQ, upperQ    

# Generate Basic Statistics Summary
def stats(vector, title = "Title", output = ""):
    lowerQ, upperQ = quartiles(vector)
    lines = [
        "%s\n" % title,
        "Measures of Central Tendency",
        "\tMean          %s"  % mean(vector),
        "\tMedian        %s"  % median(vector),
        "\tMode          %s"  % mode(vector),
        "Measures of Spread",
        "\tstdev         %s"  % stdev(vector),
        "\tvariance      %s"  % variance(vector),
        "Distribution",
        "\tMin           %s"  % min(vector),
        "\tLower Quatile %s"  % lowerQ,
        "\tMedian        %s"  % median(vector),
        "\tUpper Quatile %s"  % upperQ,
        "\tMax           %s\n"% max(vector)
    ]
    for line in lines:
        print(line)
    if len(output) > 0:
        with open(output, "w") as file:
            for line in lines:
                file.write(line)
    
def printStatistics (analysis, output = ""):
    lines = []
    for threshold, records in analysis.items():
        report = [
            "--------------------------------------",
            "Tumor LOD Threshold: %d\n" % threshold,
            "\tTrue Positive:     %d" % records["tp"],
            "\tTrue Negative:     %d" % records["tn"],
            "\tFalse Positive:    %d" % records["fp"],
            "\tFalse Negative:    %d\n" % records["fn"],
            "\tTPR:               %0.2f%%" % (records["tpr"] * 100),
            "\tTNR:               %0.2f%%" % (records["tnr"] * 100),
            "\tFPR:               %0.2f%%" % (records["fpr"] * 100),
            "\tFNR:               %0.2f%%\n" % (records["fnr"] * 100),
            "\tSensitivity:       %0.2f%%" % (records["sensitivity"] * 100),
            "\tSpecificity:       %0.2f%%\n" % (records["specificity"] * 100),
            "\tPrecision:         %0.2f%%" % (records["precision"] * 100),
            "\tRecall:            %0.2f%%\n" % (records["recall"] * 100),
            "\tAccuracy:          %0.2f%%" % (records["accuracy"] * 100),
            "\tF1-Score:          %0.2f%%" % (records["f1score"] * 100),
            "\tFDR:               %0.2f%%\n" % (records["fdr"] * 100),
            "--------------------------------------"
        ]
        lines.extend(report)
        for line in report:
            print(line)
    # Write Report if Output Provided
    if len(output) > 0:
        with open(output, "w") as file:
            for line in lines:
                file.write(line)

# Print Best by Some Metric
def printBest (analysis, metric = "accuracy", output = ""):
    threshold, best = max(enumerate([records[metric] for threshold, records in analysis.items()]), key = itemgetter(1))
    records = analysis[threshold]
    lines = [
        "Best Threshold by '%s'\n" % metric,
        "Tumor LOD Threshold: %d\n" % threshold,
        "\tTrue Positive:     %d" % records["tp"],
        "\tTrue Negative:     %d" % records["tn"],
        "\tFalse Positive:    %d" % records["fp"],
        "\tFalse Negative:    %d\n" % records["fn"],
        "\tTPR:               %0.2f%%" % (records["tpr"] * 100),
        "\tTNR:               %0.2f%%" % (records["tnr"] * 100),
        "\tFPR:               %0.2f%%" % (records["fpr"] * 100),
        "\tFNR:               %0.2f%%\n" % (records["fnr"] * 100),
        "\tSensitivity:       %0.2f%%" % (records["sensitivity"] * 100),
        "\tSpecificity:       %0.2f%%\n" % (records["specificity"] * 100),
        "\tPrecision:         %0.2f%%" % (records["precision"] * 100),
        "\tRecall:            %0.2f%%\n" % (records["recall"] * 100),
        "\tAccuracy:          %0.2f%%" % (records["accuracy"] * 100),
        "\tF1-Score:          %0.2f%%\n" % (records["f1score"] * 100),
        "\tFDR:               %0.2f%%" % (records["fdr"] * 100),
        "\tFOR:               %0.2f%%\n" % (records["for"] * 100),
        "\tNPV:               %0.2f%%" % (records["npv"] * 100),
        "\tPPV:               %0.2f%%\n" % (records["ppv"] * 100),
        "\tDiagnostic OR:     %0.2f\n" % (records["diagnosticOR"] * 100)
    ]
    for line in lines:
        print(line)
    if len(output) > 0:
        with open(output, "w") as file:
            for line in lines:
                file.write(line)
    
#
# Repeated Computations
#

def computeQuals (hashMap):
    # Compute Quality Score Approximations
    QS = getFormatData(hashMap, "QSS")[0]
    AD = getFormatData(hashMap, "AD")[0]
    # Compute Mean Quality Scores for Reference & Alternate Alleles
    refQS = []
    altQS = []
    for qsPair, adPair in zip(QS, AD):
        if adPair[0] == 0 or adPair[1] == 0:
            pass
        else:
            refQS.append((qsPair[0] / adPair[0]))
            altQS.append((qsPair[1] / adPair[1]))
    return refQS, altQS

# Split Vector by Index of Value Ranges - Default by Value
def generateSplits (vector, n = 10, median = False):
    if median:
        # Split Vectors n Times by Their Median (Attempt to Evenly Distribute Vector by Rounding, Similar to a Split by Index)
        splits    = [vector[round(i * (len(vector) / n)): round((i + 1) * (len(vector) / n))] for i in range(n)]
        vecRanges = [(min(split), max(split)) for split in splits]
    else:
        # Sum Vector and Distribute Evenly by Value
        vecMin, vecMax = min(vector), max(vector)
        vecRanges = [(vecMin + (x * ((vecMax - vecMin) / n)), vecMin + ((x + 1) * ((vecMax - vecMin) / n))) for x in range(n)]
        splits    = [[x for x in vector if (x >= vecRange[0] and x < vecRange[1])] for vecRange in vecRanges]
    return splits, vecRanges

def getRanges (splits):
    ranges = []
    for split in splits:
        ranges.append([min(split), max(split)])
    return ranges    

#
# Interface Management & Output Formatting
#
    
# Print a Break
def brk(char = "-", n = 50):
    print("\n" + ("-" * n) + "\n")    

#
# Plotting
#

# Generate a Simple Inline Plot
def simpleInlinePlot (x, y, title = "", xlabel = "", ylabel = "", color = "green", plotType = "plot", size = (8, 8), dpi = 300, domain = [], range = [], opacity = 0.1, bins = 1000):
    # Setup
    fig = matplotlib.pylot.figure(1, size, dpi = dpi)
    domain = domain if len(domain) == 2 else [min(x), max(x)]
    range  = range if len(range) == 2 else [min(y), max(y)]
    # Plot
    matplotlib.pylot.subplot(1,1,1)
    if plotType == "hist":
        getattr(Plot, plotType)(y, bins = bins, color=color, lw=2, alpha = opacity)
    if plotType == "plot":
        matplotlib.pylot.xlim(domain)
        matplotlib.pylot.ylim(range)
        getattr(Plot, plotType)(x, y, color=color, lw=2, alpha = opacity)
        matplotlib.pylot.plot([0, 1], [0, 95000], color='navy', lw=2, linestyle='--')
    if plotType == "scatter":
        matplotlib.pylot.xlim(domain)
        matplotlib.pylot.ylim(range)
        getattr(Plot, plotType)(x, y, color = color, alpha = opacity)
    axis = fig.add_subplot(1,1,1)
    axis.set_title(title, fontsize = 12)
    axis.set_xlabel(xlabel, fontsize = 10)
    axis.set_ylabel(ylabel, fontsize = 10)
    if plotType == "bar":
        assert len(x) == len(y)
        xs = [n for n in range(len(x))]
        rects = axis.bar(xs, y, color = color, align = "center")
        i = 0
        for rect in rects:
            h = rect.get_height()
            label = x[i]
            axis.text(rect.get_x() + (rect.get_width() / 2), 1.05 * h, label, ha = 'center', va = 'bottom')
            i += 1
    
    return fig
    
# Generates a Plot Scatter
def scatterXY (analysis, x, y, color = "green", opacity = 0.5, n = 100, skip = 1, method = "scatter"):
    X = [analysis[threshold][x] for threshold in range(0, n, skip)]
    Y = [analysis[threshold][y] for threshold in range(0, n, skip)]
    return getattr(Plot, method)(X, Y, color = color, alpha = opacity), X, Y
    
def plotCurve (conditions, labels, xKey = "", yKey = "", xLabel = "", yLabel = "", title = "Title", colors = ["green", "orange", "red", "purple", "blue"], diag = "tl", method = "plot", legendLoc = "lower left"):
    # Setup
    fig         = matplotlib.pylot.figure(1, (8, 8), dpi = 300)
    axis        = fig.add_subplot(111)
    axis.set_xlabel(xLabel, fontsize = 10)
    axis.set_ylabel(yLabel, fontsize = 10)
    # Scatter
    legendKeys  = []
    i           = 0
    nLabels     = []
    assert len(conditions) == len(labels)
    # Plot Conditions
    for condition in conditions:
        colorIdx = i % len(colors)
        scattered, X, Y = scatterXY(condition, xKey, yKey, colors[colorIdx], method = method)
        auc = metrics.auc(X, Y)
        legendKeys.append(scattered)
        label = "%s (AUC = %.4f)" % (str(labels[i]), auc) 
        nLabels.append(label)
        i += 1
        
    # Plot Diagonal
    if diag is "tl":
        matplotlib.pylot.plot((1, 0), color='navy', lw=1, linestyle='--')
    if diag is "tr":
        matplotlib.pylot.plot((0, 1), color='navy', lw=1, linestyle='--')
    matplotlib.pylot.ylim([0, 1])
    matplotlib.pylot.xlim([0, 1])
    matplotlib.pylot.title(title, fontsize = 12)
    # Plot Legend
    if method is "scatter":
        matplotlib.pylot.legend(legendKeys, nLabels, loc = legendLoc, fontsize = 10)
    else:
        handles = [matplotlib.patches.Patch(color = key[0]._color, label = label) for key, label in zip(legendKeys, nLabels)]
        matplotlib.pylot.legend(handles=handles, loc = legendLoc, fontsize = 10)
    # Save & Render
    
    return fig
    
# Plot Maximum of Two Combined Values (e.g. Maximize TPR & TNR)
def plotMaximized (analysis, x, y, color = "red", n = 100, skip = 1):
    combinations = [analysis[threshold][x] + analysis[threshold][y] for threshold in range(0, n, skip)]
    index, combination = max(enumerate(combinations), key = itemgetter(1))
    return matplotlib.pylot.scatter(analysis[index][x], analysis[index][y], s=200, facecolors = 'none', edgecolors = color)

# Get Analysis of Threshold with Maximum Combined Metrics (e.g. sensitivity and specificity or tpr and fpr)
def getMaximizedAnalysis (analysis, x, y, skip = 1):
    combinations = [analysis[threshold][x] + analysis[threshold][y] for threshold in range(0, len(analysis), skip)]
    index, combination = max(enumerate(combinations), key = itemgetter(1))
    threshold = index + 1
    return threshold, analysis[threshold]

# Quick & Dirty Function to Plot Multiple Curves & Stop me Polluting the Global Scope. Tsk Tsk. 
def plotCurves (c1Tups, c2Tups, truthSet, c1Map, c2Map, subsetKey = "", xKey = "specificity", yKey = "sensitivity", diag = "tl", legendLoc = "lower left", c1Label = "", c2Label = "", titleStr = ""):
    i = 1
    plots = []
    for c1, c2 in zip(c1Tups, c2Tups):
        # Subset Hashmaps
        c1Hm, c1Exc = subsetHashMapByFormatRange(c1Map, c1, key = subsetKey)
        c2Hm, c2Exc = subsetHashMapByFormatRange(c2Map, c2, key = subsetKey)
        # Regenerate Analysis on Subset
        c1A = generateAnalyses(c1Hm, truthSet, "TLOD")
        c2A = generateAnalyses(c2Hm, truthSet, "TLOD")
        c1Title = c1Label + " Range (%.2f, %.2f)" % tuple(c1)
        c2Title = c2Label + " Range (%.2f, %.2f)" % tuple(c2)
        c1Legend = str(c1Label) + ", n = %d" % len(c1Hm.keys())
        c2Legend = str(c2Label) + ", n = %d" % len(c2Hm.keys())
        xLab = xKey[0].upper() + xKey[1:]
        yLab = yKey[0].upper() + yKey[1:]
        title = titleStr + " (Curve Pair " + str(i) + ")\n" +  c1Title + " / " + c2Title
        # Call CurvePlotting Function
        plot = plotCurve([c1A, c2A], (c1Legend, c2Legend), xKey = xKey, yKey = yKey, xLabel = xLab, yLabel = yLab, title = title, method = "plot", legendLoc = legendLoc,  diag = diag)
        plots.append(plot)
        i += 1
        print("Records Excluded From Analysis: %d (%s) – %d (%s)" % (len(c1Exc.keys()), c1Label, len(c2Exc.keys()), c2Label))
    return plots
        
# Plot Multiple Curves on One Canvas
def multiCurve (c1Tups, c2Tups, truthSet, c1Map, c2Map, subsetKey = "", c1Label = "", c2Label = "", titleStr = "", xKey = "sensitivity", yKey = "specificity", xLabel = "Sensitivity", yLabel = "Specificity", title = "", colors = [["#004400", "#005300", "#006200", "#007100", "#008000", "#008F00", "#009E00", "#00AD00", "#00BC00", "#00CC00"], ["#601D00", "#712300", "#822A00", "#933100", "#A43800", "#B53E00", "#C64500", "#D74C00", "#E85300", "#FA5A00"]], diag = "tl", method = "plot", legendLoc = "lower left"):
    # Compute Curves
    analyses = []
    for c1, c2 in zip(c1Tups, c2Tups):
        # Subset Hashmaps
        c1Hm, excluded = subsetHashMapByFormatRange(c1Map, c1, key = subsetKey)
        c2Hm, excluded = subsetHashMapByFormatRange(c2Map, c2, key = subsetKey)
        # Regenerate Analysis on Subset
        c1Analysis = generateAnalyses(c1Hm, truthSet, "TLOD")
        c2Analysis = generateAnalyses(c2Hm, truthSet, "TLOD")
        analyses.append([c1Analysis, c2Analysis])
    # Plot Setup
    fig         = matplotlib.pylot.figure(1, (8, 8), dpi = 300)
    axis        = fig.add_subplot(111)
    axis.set_xlabel(xLabel, fontsize = 10)
    axis.set_ylabel(yLabel, fontsize = 10)
    alphaSteps = len(analyses) + 1
    i = 0
    # For Each Analysis
    for analysis in analyses:
        legendKeys = []
        j = 0
        alpha = ((i + 1) / alphaSteps * 0.8)
        # For Each Condition
        for condition in analysis:
            colorIdx = i % len(colors[j])
            # Call Scatter Function, 
            legendKeys.append(scatterXY(condition, xKey, yKey, colors[j][colorIdx], opacity = alpha, method = method))
            j += 1
        i += 1
    if diag is "tl":
        matplotlib.pylot.plot((1, 0), color='navy', lw=1, linestyle='--')
    if diag is "tr":
        matplotlib.pylot.plot((0, 1), color='navy', lw=1, linestyle='--')
    matplotlib.pylot.ylim([0, 1])
    matplotlib.pylot.xlim([0, 1])
    matplotlib.pylot.title(title, fontsize = 12)
    if method is "scatter":
        matplotlib.pylot.legend(legendKeys, labels, loc = legendLoc, fontsize = 10)
    else:
        labels = [c1Label, c2Label]
        handles = [matplotlib.patches.Patch(color = key[0]._color, label = label) for key, label in zip(legendKeys, labels)]
        matplotlib.pylot.legend(handles=handles, loc = legendLoc, fontsize = 10)
    
    return fig

def plotQualityScores (ref, alt, title = "Title"):
    # Setup
    fig         = matplotlib.pylot.figure(1, (16, 12), dpi = 300)
    axis        = fig.add_subplot(111)
    axis.set_ylabel('Quality Score', fontsize = 10)
    axis.set_xlabel('Variant Allele', fontsize = 10)
    # Scatter
    r = matplotlib.pylot.scatter(range(len(ref)), ref, color = "orange", alpha = 0.2)
    a = matplotlib.pylot.scatter(range(len(alt)), alt, color = "blue", alpha = 0.2)
    # Plot
    matplotlib.pylot.ylim([0, 40])
    matplotlib.pylot.xlim([0, len(ref)])
    matplotlib.pylot.title(title, fontsize = 12)
    matplotlib.pylot.legend((r, a), ('Reference', 'Alternate'), loc = 'lower right', fontsize = 10)
    
    return fig
    
# Plot Mutect Filters by Frequency
def plotFilterFrequencies (hashmap, title = "", xlabel = "", ylabel = ""):
    # Which Filters are Catching Variants?
    counts = defaultdict(int)
    for md5hash, record in hashmap.items():
        if len(record.FILTER) > 0:
            for mutectFilter in record.FILTER:
                counts[mutectFilter] += 1
    # Get Frequencies
    frequencies = {key : (value / len(hashmap)) for key, value in counts.items()}
    # Plot Bar Plot
    return simpleInlinePlot(list(frequencies.keys()), list(frequencies.values()), title = title, xlabel = xlabel, ylabel = ylabel, plotType = "bar")
    
def plotMetrics(analysis, metrics = ["sensitivity", "specificity"], reduce = lambda x: len(x), title = "Title", colors = ["red", "green", "yellow", "orange"]):
    metricsDict = { metric : [] for metric in metrics }
    thresholds  = []
    for threshold, records in analysis.items():
        thresholds.append(threshold)
        for metric in metrics:
            if hasattr(records[metric], '__iter__'): 
                metricsDict[metric].append(len(records[metric]))
            else:
                metricsDict[metric].append(records[metric])
    index = 0
    figure = matplotlib.pylot.figure(1, (8, 8), dpi = 300) 
    axis = figure.add_subplot(1,1,1)
    axis.set_title(title, fontsize = 12)
    axis.set_xlabel("Threshold", fontsize = 10)
    axis.set_ylabel("Performance", fontsize = 10)
    box = axis.get_position()
    axis.set_position([box.x0, box.y0, box.width, box.height * 0.9])
    for metric, values in metricsDict.items():
        color = colors[index % len(colors)]
        label = metric
        x = thresholds
        y = values
        matplotlib.pylot.subplot(1,1,1)
        matplotlib.pylot.plot(x, y, label = label, color = color, alpha = 0.6, lw = 2)
        matplotlib.pylot.legend(fontsize = 10, loc = "upper center", ncol = 5, bbox_to_anchor=(0.5, -0.05))
        index += 1
    return fig

#
# Main Routine
#

if __name__ == "__main__":

    # Imports
    import argparse

    # Parse Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--set", type = str, help = "Data Set (e.g. 'set1')")
    parser.add_argument("-c", "--condition", type = str, help = "Error Model Condition")
    parser.add_argument("-o", "--output", type = str, help = "Output Prefix")
    argsDict = vars(parser.parse_args())
    dataset = argsDict["set"]
    condition = argsDict["condition"]
    outputdir = argsDict["output"]

    assert dataset is not None, "Set Argument Required (-s | --set)"
    assert condition is not None, "Condition Argument Required (-c | --condition)"

    if outputdir is None:
        outputdir = "output/%s/%s" % (dataset, condition)

    # Selecte Appropriate Truth Set
    truthsets = [
        "synthetic.challenge.set1.tumor.all.truth.vcf",
        "synthetic.challenge.set2.tumor.all.truth.vcf",
        "synthetic.challenge.set3.tumor.20pctmasked.truth.vcf"
    ]
    selected = [truthset for truthset in truthsets if dataset in truthset][0]

    # Print Configuration
    print("%s – %s – %s" % (dataset, condition, selected))

    # VCF Locations
    condition1VCF  = "vcf/synthetic.challenge.%s.%s.default.%s.raw.snps.indels.vcf" % (dataset, condition, subcondition)
    condition2VCF  = "vcf/synthetic.challenge.%s.%s.default.bqsr.raw.snps.indels.vcf" % (dataset, condition)
    truthVCF       = "vcf/truth-set/%s" % selected

    # Init VCF Readers
    c1Reader       = VCF.Reader(open(condition1VCF, 'r'))
    c2Reader       = VCF.Reader(open(condition2VCF, 'r'))
    trReader       = VCF.Reader(open(truthVCF, 'r'))

    # Variant Record Properties to Define Uniqueness
    uniqueProps    = ["CHROM", "POS", "REF", "ALT"]

    #
    # Store Metadata References
    #

    # INFO Terms
    c1Infos        = c1Reader.infos
    c2Infos        = c2Reader.infos
    trInfos        = trReader.infos

    # Metadata
    c1Metadata     = c1Reader.metadata
    c2Metadata     = c2Reader.metadata
    trMetadata     = trReader.metadata

    # Filters
    c1Filters      = c1Reader.filters
    c2Filters      = c2Reader.filters
    trFilters      = trReader.filters

    # Formats
    c1Formats      = c1Reader.formats
    c2Formats      = c2Reader.formats
    trFormats      = trReader.formats

#
# Generate Variant HashMaps
#

    # Hash VCF Records for Comparison - Use Values of Unique Props to Generate Hash
    print("Hashing Condition 1...")
    c1HashMap = hashVariants(c1Reader, uniqueProps, algorithm = "sha1")
    print("Hashing Condition 2...")
    c2HashMap = hashVariants(c2Reader, uniqueProps, algorithm = "sha1")
    print("Hashing Truth Set...")
    trHashMap = hashVariants(trReader, uniqueProps, algorithm = "sha1")
    print("Done")

    # Retrieve Passing Variants
    c1Passing = passingVariants(c1HashMap)
    c2Passing = passingVariants(c2HashMap)

    # Get Sets      
    trSet   = set(trHashMap.keys())
    c1Set   = set(c1HashMap.keys())
    c2Set   = set(c2HashMap.keys())

    # Lengths
    tr      = len(trSet)
    c1      = len(c1Set)
    c2      = len(c2Set)

    #
    # BQSR / No BQSR Interval Performance Comparison - All Variants
    #

    # Plot Across 5% Ranges
    for i in range(20):
        lower = (i / 2) / 10
        upper = ((i + 1) / 2) / 10
        lowerStr = str(lower)
        upperStr = str(upper)
        figure = plotCurves([(lower, upper)], [(lower, upper)], trSet, c1HashMap, c2HashMap, subsetKey = "AF", c1Label = "No BQSR", c2Label = "BQSR", titleStr = "Sensitivity-Specificity by All Variant Allele Frequency Interval")[0]
        figure.savefig(outputdir + "/senspec-" + lowerStr + "-" + upperStr + "-all-variants.pdf", bbox_inches='tight')
        figure.clf()

    # Plot Across 0.1% LFV Ranges
    for i in range(15):
        lower = i / 100
        upper = (i + 1) / 100
        lowerStr = str(lower)
        upperStr = str(upper)
        figure = plotCurves([(lower, upper)], [(lower, upper)], trSet, c1HashMap, c2HashMap, subsetKey = "AF", c1Label = "No BQSR", c2Label = "BQSR", titleStr = "Sensitivity-Specificity by All Variant Allele Frequency Interval")[0]
        figure.savefig(outputdir + "/senspec-" + lowerStr + "-" + upperStr + "-all-variants.pdf", bbox_inches='tight')
        figure.clf()
    
    #
    # BQSR / No BQSR Interval Performance Comparison - Passing Variants
    #

    # Plot Across 5% Ranges
    for i in range(20):
        lower = (i / 2) / 10
        upper = ((i + 1) / 2) / 10
        lowerStr = str(lower)
        upperStr = str(upper)
        figure = plotCurves([(lower, upper)], [(lower, upper)], trSet, c1Passing, c2Passing, subsetKey = "AF", c1Label = "No BQSR", c2Label = "BQSR", titleStr = "Sensitivity-Specificity by Passing Variant Allele Frequency Interval")[0]
        figure.savefig(outputdir + "/senspec-" + lowerStr + "-" + upperStr + "-passing-variants.pdf", bbox_inches='tight')
        figure.clf()
    
    # Plot Across 0.1% LFV Ranges
    for i in range(15):
        lower = i / 100
        upper = (i + 1) / 100
        lowerStr = str(lower)
        upperStr = str(upper)
        figure = plotCurves([(lower, upper)], [(lower, upper)], trSet, c1Passing, c2Passing, subsetKey = "AF", c1Label = "No BQSR", c2Label = "BQSR", titleStr = "Sensitivity-Specificity by Passing Variant Allele Frequency Interval")[0]
        figure.savefig(outputdir + "/senspec-" + lowerStr + "-" + upperStr + "-passing-variants.pdf", bbox_inches='tight')
        figure.clf()

#
# Which Mutect2 Filters Are Creating False Negatives?
#

    #
    # Truth Set
    #

    # Subset HashMaps by VAF Range
    trSub, excluded = subsetHashMapByInfoThreshold(trHashMap, 0.15, operation = "<=", key = "VAF")
    # Truth Set 0.10 - 0.15
    trSubGt10, trExcludedLt10 = subsetHashMapByInfoThreshold(trSub, 0.10, operation = "<=", key = "VAF")
    # Truth Set 0.00 - 0.05; 0.05 - 0.10
    trSubGt00, trSubGt05 = subsetHashMapByInfoThreshold(trExcludedLt10, 0.05, operation = "<=", key = "VAF")

    # Print / Write Analyses
    with open(outputdir + "/true-lfv-counts", "w") as output:
        output.write("n True Variants <= 0.05 VAF: %d\n" % len(trSubGt00))
        output.write("n True Variants <= 0.10 VAF: %d\n" % (len(trSubGt00) + len(trSubGt05)))
        output.write("n True Variants <= 0.15 VAF: %d\n" % (len(trSubGt00) + len(trSubGt05) + len(trSubGt10)))
        output.write("n True Variants in range (0.00 - <= 0.05: %d\n" % len(trSubGt00))
        output.write("n True Variants in range (0.05 - <= 0.10: %d\n" % len(trSubGt05))
        output.write("n True Variants in range (0.10 - <= 0.15: %d\n" % len(trSubGt10))

    #
    # No BQSR - All
    #

    c1Sub, excluded = subsetHashMapByFormatThreshold(c1HashMap, 0.15, operation = "<=", key = "AF")
    # All - No BQSR Set 0.10 - 0.15
    c1SubGt10, c1ExcludedLt10 = subsetHashMapByFormatThreshold(c1Sub, 0.10, operation = "<=", key = "AF")
    # All - No BQSR Set 0.00 - 0.05; 0.05 - 0.10
    c1SubGt00, c1SubGt05 = subsetHashMapByFormatThreshold(c1ExcludedLt10, 0.05, operation = "<=", key = "AF")

    #
    # No BQSR - Passing
    #

    c1PassSub, excluded = subsetHashMapByFormatThreshold(c1Passing, 0.15, operation = "<=", key = "AF")
    # Passing - No BQSR Set 0.10 - 0.15
    c1PassSubGt10, c1PassExcludedLt10 = subsetHashMapByFormatThreshold(c1PassSub, 0.10, operation = "<=", key = "AF")
    # Passing - No BQSR Set 0.00 - 0.05; 0.05 - 0.10
    c1PassSubGt00, c1PassSubGt05 = subsetHashMapByFormatThreshold(c1PassExcludedLt10, 0.05, operation = "<=", key = "AF")

    #
    # BQSR - All
    #

    c2Sub, excluded = subsetHashMapByFormatThreshold(c2HashMap, 0.15, operation = "<=", key = "AF")
    # All - BQSR Set 0.10 - 0.15
    c2SubGt10, c2ExcludedLt10 = subsetHashMapByFormatThreshold(c2Sub, 0.10, operation = "<=", key = "AF")
    # All - BQSR Set 0.00 - 0.05; 0.05 - 0.10
    c2SubGt00, c2SubGt05 = subsetHashMapByFormatThreshold(c2ExcludedLt10, 0.05, operation = "<=", key = "AF")

    #
    # BQSR - Passing
    #

    c2PassSub, excluded = subsetHashMapByFormatThreshold(c2Passing, 0.15, operation = "<=", key = "AF")
    # Passing - BQSR Set 0.10 - 0.15
    c2PassSubGt10, c2PassExcludedLt10 = subsetHashMapByFormatThreshold(c2PassSub, 0.10, operation = "<=", key = "AF")
    # Passing - BQSR Set 0.00 - 0.05; 0.05 - 0.10
    c2PassSubGt00, c2PassSubGt05 = subsetHashMapByFormatThreshold(c2PassExcludedLt10, 0.05, operation = "<=", key = "AF")

    # Plot Filter Frequencies – False Negatives - No BQSR Set 0.10 - 0.15

    # Get True Keys For False Negatives
    diff = set(c1SubGt10.keys()) - set(c1PassSubGt10.keys())
    c1SubGt10Fn = subsetHashMapByKeys(c1SubGt10, diff)
    if len(c1SubGt10Fn) > 0:
        # Generate Filter Frequencies Plot
        figure = plotFilterFrequencies(c1SubGt10Fn, title = "Filter Frequencies for False Negatives (No BQSR) - VAF [0.10 - 0.15]", xlabel = "Filters", ylabel = "Proportion")
        figure.savefig(outputdir + "/fn-filtfreq-0.10-0.15-nobqsr-variants.pdf", bbox_inches='tight')
        figure.clf()

    # Plot Filter Frequencies – False Negatives - BQSR Set 0.10 - 0.15

    # Get True Keys For False Negatives
    diff = set(c2SubGt10.keys()) - set(c2PassSubGt10.keys())
    c2SubGt10Fn = subsetHashMapByKeys(c2SubGt10, diff)
    if len(c2SubGt10Fn) > 0:
        # Generate Filter Frequencies Plot
        figure = plotFilterFrequencies(c2SubGt10Fn, title = "Filter Frequencies for False Negatives (BQSR) - VAF [0.10 - 0.15]", xlabel = "Filters", ylabel = "Proportion")
        figure.savefig(outputdir + "/fn-filtfreq-0.10-0.15-bqsr-variants.pdf", bbox_inches='tight')
        figure.clf()

    # Plot Filter Frequencies – False Negatives - No BQSR Set 0.05 - 0.10

    # Get True Keys For False Negatives
    diff = set(c1SubGt05.keys()) - set(c1PassSubGt05.keys())
    c1SubGt05Fn = subsetHashMapByKeys(c1SubGt05, diff)
    if len(c1SubGt05Fn) > 0:
        # Generate Filter Frequencies Plot
        figure = plotFilterFrequencies(c1SubGt05Fn, title = "Filter Frequencies for False Negatives (No BQSR) - VAF [0.05 - 0.10]", xlabel = "Filters", ylabel = "Proportion")
        figure.savefig(outputdir + "/fn-filtfreq-0.05-0.10-nobqsr-variants.pdf", bbox_inches='tight')
        figure.clf()

    # Plot Filter Frequencies – False Negatives - BQSR Set 0.05 - 0.10

    # Get True Keys For False Negatives
    diff = set(c2SubGt05.keys()) - set(c2PassSubGt05.keys())
    c2SubGt05Fn = subsetHashMapByKeys(c2SubGt05, diff)
    if len(c2SubGt05Fn) > 0:
        # Generate Filter Frequencies Plot
        figure = plotFilterFrequencies(c2SubGt05Fn, title = "Filter Frequencies for False Negatives (BQSR) - VAF [0.05 - 0.10]", xlabel = "Filters", ylabel = "Proportion")
        figure.savefig(outputdir + "/fn-filtfreq-0.05-0.10-bqsr-variants.pdf", bbox_inches='tight')
        figure.clf()

    # Plot Filter Frequencies – False Negatives - No BQSR Set < 0.01 - 0.05

    # Get True Keys For False Negatives
    diff = set(c1SubGt00.keys()) - set(c1PassSubGt00.keys())
    c1SubGt00Fn = subsetHashMapByKeys(c1SubGt00, diff)
    if len(c1SubGt00Fn) > 0:
        # Generate Filter Frequencies Plot
        figure = plotFilterFrequencies(c1SubGt00Fn, title = "Filter Frequencies for False Negatives (No BQSR) - VAF [0.00 - 0.05]", xlabel = "Filters", ylabel = "Proportion")
        figure.savefig(outputdir + "/fn-filtfreq-0.00-0.05-nobqsr-variants.pdf", bbox_inches='tight')
        figure.clf()

    # Plot Filter Frequencies – False Negatives - BQSR Set < 0.01 - 0.05

    # Get True Keys For False Negatives
    diff = set(c2SubGt00.keys()) - set(c2PassSubGt00.keys())
    c2SubGt00Fn = subsetHashMapByKeys(c2SubGt00, diff)
    if len(c2SubGt00Fn) > 0:
        # Generate Filter Frequencies Plot
        figure = plotFilterFrequencies(c2SubGt00Fn, title = "Filter Frequencies for False Negatives (BQSR) - VAF [0.00 - 0.05]", xlabel = "Filters", ylabel = "Proportion")
        figure.savefig(outputdir + "/fn-filtfreq-0.00-0.05-bqsr-variants.pdf", bbox_inches='tight')
        figure.clf()
    # Mutect2 Performance at Low VAF
    # Called / Not Called by MuTect2 by VAF Range [> 0.00 - 0.05]; [0.05 - 0.10]; [0.10 - 0.15]

    #
    # Map Analysis
    #

    trCombinations = [trSubGt10, trSubGt05, trSubGt00]
    combinations = [
        [
            {"label": "No BQSR - All - [0.10 - 0.15] VAF", "data": c1SubGt10}, 
            {"label": "No BQSR - All - [0.05 - 0.10] VAF", "data": c1SubGt05},
            {"label": "No BQSR - All - [> 0.00 - 0.05] VAF", "data": c1SubGt00}
        ],
        [
            {"label": "BQSR - All - [0.10 - 0.15] VAF", "data": c2SubGt10}, 
            {"label": "BQSR - All - [0.05 - 0.10] VAF", "data": c2SubGt05}, 
            {"label": "BQSR - All - [> 0.00 - 0.05] VAF", "data": c2SubGt00}
        ],
        [
            {"label": "No BQSR - Passing - [0.10 - 0.15] VAF", "data": c1PassSubGt10}, 
            {"label": "No BQSR - Passing - [0.05 - 0.10] VAF", "data": c1PassSubGt05}, 
            {"label": "No BQSR - Passing - [> 0.00 - 0.05] VAF", "data": c1PassSubGt00}
        ],
        [
            {"label": "BQSR - Passing - [0.10 - 0.15] VAF", "data": c2PassSubGt10}, 
            {"label": "BQSR - Passing - [0.05 - 0.10] VAF", "data": c2PassSubGt05}, 
            {"label": "BQSR - Passing - [> 0.00 - 0.05] VAF", "data": c2PassSubGt00}
        ]
    ]

#
# Curve Analysis
#

    # Individual Plots
    figure = plotCurves([(0.0, 1.0)], [(0.0, 1.0)], trSet, c1HashMap, c2HashMap, subsetKey = "AF", xKey = "fpr", yKey = "tpr", diag = "tr", legendLoc = "lower right", c1Label = "No BQSR", c2Label = "BQSR", titleStr = "ROC Curve by Tumor LOD Thresholds for All Variants (No BQSR v. BQSR)")[0]
    figure.savefig(outputdir + "/roc-all-variants.pdf", bbox_inches='tight')   
    figure.clf()

    # Individual Plots
    figure = plotCurves([(0.0, 1.0)], [(0.0, 1.0)], trSet, c1Passing, c2Passing, subsetKey = "AF", xKey = "fpr", yKey = "tpr", diag = "tr", legendLoc = "lower right", c1Label = "No BQSR", c2Label = "BQSR", titleStr = "ROC Curve by Tumor LOD Thresholds for Passing Variants (TLOD >= 10.0) (No BQSR v. BQSR)")[0]
    figure.savefig(outputdir + "/roc-passing-variants.pdf", bbox_inches='tight')
    figure.clf()

    # Sensitivity Specificity Curve (No BQSR Condition v. BQSR Condition)

    #
    # Generate Sensitivity-Specificity Curve
    #

    # Individual Plots
    figure = plotCurves([(0.0, 1.0)], [(0.0, 1.0)], trSet, c1HashMap, c2HashMap, subsetKey = "AF", c1Label = "No BQSR", c2Label = "BQSR", titleStr = "Sensitivity Specificity Curve – Tumor LOD Thresholds – All Variants (No BQSR v. BQSR)")[0]
    figure.savefig(outputdir + "/senspec-all-variants.pdf", bbox_inches='tight')
    figure.clf()

    # Individual Plots
    figure = plotCurves([(0.0, 1.0)], [(0.0, 1.0)], trSet, c1Passing, c2Passing, subsetKey = "AF", c1Label = "No BQSR", c2Label = "BQSR", titleStr = "Sensitivity Specificity Curve – Tumor LOD Thresholds – Passing Variants (No BQSR v. BQSR)")[0]
    figure.savefig(outputdir + "/senspec-passing-variants.pdf", bbox_inches='tight')
    figure.clf()

#
# Additional Passing Analysis
#

    # Second Look: What Threshold Results in The Max Combined Sensitivity-Specificity?

    # Generate Analyses - With Tumor LOD As Metric to Compute Statistics On
    condition1Analysis = generateAnalyses(c1Passing, trSet, "TLOD")
    del c1Passing
    condition2Analysis = generateAnalyses(c2Passing, trSet, "TLOD")
    del c2Passing

    #
    # Generate Sensitivity-Specificity Curve - With Maximum Combination
    #

    # Setup
    fig         = matplotlib.pylot.figure(1, (8, 8), dpi = 300)
    axis        = fig.add_subplot(111)
    axis.set_ylabel('Sensitivity', fontsize = 10)
    axis.set_xlabel('Specificity', fontsize = 10)
    # Scatter
    c1, c1X, c1Y = scatterXY(condition1Analysis, "specificity", "sensitivity", "green")
    c2, c2X, c2Y = scatterXY(condition2Analysis, "specificity", "sensitivity", "orange")
    c1Auc = metrics.auc(c1X, c1Y)
    c2Auc = metrics.auc(c2X, c2Y)
    c1Label = "No BQSR (AUC = %.4f)" % c1Auc
    c2Label = "BQSR (AUC = %.4f)" % c2Auc
    # Annotate Maximized Threshold
    m  = plotMaximized(condition1Analysis, "specificity", "sensitivity")
    # Plot
    matplotlib.pylot.plot((1, 0), color='navy', lw=1, linestyle='--')
    matplotlib.pylot.ylim([0, 1])
    matplotlib.pylot.xlim([0, 1])
    matplotlib.pylot.title("Sensitivity Specificity Curve – Tumor LOD Thresholds – Passing Variants (No BQSR v. BQSR)", fontsize = 12)
    matplotlib.pylot.legend((c1, c2), (c1Label, c2Label), loc = 'lower left', fontsize = 10)
    fig.savefig(outputdir + "/senspec-passing-optimal-thresh.pdf", bbox_inches='tight')
    fig.clf()    
        
    #        
    # Optimal Threshold Analysis
    #

    # No BQSR Condition Statistics
    printBest(condition1Analysis, output=outputdir + "/c1-passing-analysis-accuracy") # Accuracy
    printBest(condition1Analysis, "f1score", output=outputdir + "/c1-passing-analysis-f1score")
    printBest(condition1Analysis, "ppv", output=outputdir + "/c1-passing-analysis-ppv")
    c1MaxThresh, c1MaxAnalysis = getMaximizedAnalysis(condition1Analysis, "sensitivity", "specificity")
    # Print / Write Analyses
    with open(outputdir + "/max-sens-spec-passing-nobqsr", "w") as output:
        output.write("==================================================\n")
        output.write("Threshold: %d\n" % c1MaxThresh)
        for key, value in c1MaxAnalysis.items():
            line = "%s = %s\n" % (str(key), str(value))
            output.write(line)
        output.write("==================================================\n")

    # BQSR Condition Statistics
    printBest(condition2Analysis, output=outputdir + "/c2-passing-analysis-accuracy") # Accuracy
    printBest(condition2Analysis, "f1score", output=outputdir + "/c2-passing-analysis-f1score")
    printBest(condition2Analysis, "ppv", output=outputdir + "/c2-passing-analysis-ppv")
    c2MaxThresh, c2MaxAnalysis = getMaximizedAnalysis(condition2Analysis, "sensitivity", "specificity")
    # Print / Write Analyses
    with open(outputdir + "/max-sens-spec-passing-bqsr", "w") as output:
        output.write("==================================================\n")
        output.write("Threshold: %d\n" % c2MaxThresh)
        for key, value in c2MaxAnalysis.items():
            line = "%s = %s\n" % (str(key), str(value))
            output.write(line)
        output.write("==================================================\n")

    del condition1Analysis
    del condition2Analysis

#
# Additional All Analysis
#

    # Second Look: What Threshold Results in The Max Combined Sensitivity-Specificity?

    # Generate Analyses - With Tumor LOD As Metric to Compute Statistics On
    condition1Analysis = generateAnalyses(c1HashMap, trSet, "TLOD")
    del c1HashMap
    condition2Analysis = generateAnalyses(c2HashMap, trSet, "TLOD")
    del c2HashMap

    #
    # Generate Sensitivity-Specificity Curve - With Maximum Combination
    #

    # Setup
    fig         = matplotlib.pylot.figure(1, (8, 8), dpi = 300)
    axis        = fig.add_subplot(111)
    axis.set_ylabel('Sensitivity', fontsize = 10)
    axis.set_xlabel('Specificity', fontsize = 10)
    # Scatter
    c1, c1X, c1Y = scatterXY(condition1Analysis, "specificity", "sensitivity", "green")
    c2, c2X, c2Y = scatterXY(condition2Analysis, "specificity", "sensitivity", "orange")
    c1Auc = metrics.auc(c1X, c1Y)
    c2Auc = metrics.auc(c2X, c2Y)
    c1Label = "No BQSR (AUC = %.4f)" % c1Auc
    c2Label = "BQSR (AUC = %.4f)" % c2Auc
    # Annotate Maximized Threshold
    m  = plotMaximized(condition1Analysis, "specificity", "sensitivity")
    # Plot
    matplotlib.pylot.plot((1, 0), color='navy', lw=1, linestyle='--')
    matplotlib.pylot.ylim([0, 1])
    matplotlib.pylot.xlim([0, 1])
    matplotlib.pylot.title("Sensitivity Specificity Curve – Tumor LOD Thresholds – All Variants (No BQSR v. BQSR)", fontsize = 12)
    matplotlib.pylot.legend((c1, c2), (c1Label, c2Label), loc = 'lower left', fontsize = 10)
    fig.savefig(outputdir + "/senspec-all-optimal-thresh.pdf", bbox_inches='tight')
    fig.clf()
        
    #        
    # Optimal Threshold Analysis
    #

    # No BQSR Condition Statistics
    printBest(condition1Analysis, output=outputdir + "/c1-all-analysis-accuracy") # Accuracy
    printBest(condition1Analysis, "f1score", output=outputdir + "/c1-all-analysis-f1score")
    printBest(condition1Analysis, "ppv", output=outputdir + "/c1-all-analysis-ppv")
    c1MaxThresh, c1MaxAnalysis = getMaximizedAnalysis(condition1Analysis, "sensitivity", "specificity")
    # Print / Write Analyses
    with open(outputdir + "/max-sens-spec-all-nobqsr", "w") as output:
        output.write("==================================================\n")
        output.write("Threshold: %d\n" % c1MaxThresh)
        for key, value in c1MaxAnalysis.items():
            line = "%s = %s\n" % (str(key), str(value))
            output.write(line)
        output.write("==================================================\n")

    # BQSR Condition Statistics
    printBest(condition2Analysis, output=outputdir + "/c2-all-analysis-accuracy") # Accuracy
    printBest(condition2Analysis, "f1score", output=outputdir + "/c2-all-analysis-f1score")
    printBest(condition2Analysis, "ppv", output=outputdir + "/c2-all-analysis-ppv")
    c2MaxThresh, c2MaxAnalysis = getMaximizedAnalysis(condition2Analysis, "sensitivity", "specificity")
    # Print / Write Analyses
    with open(outputdir + "/max-sens-spec-all-bqsr", "w") as output:
        output.write("==================================================\n")
        output.write("Threshold: %d\n" % c2MaxThresh)
        for key, value in c2MaxAnalysis.items():
            line = "%s = %s\n" % (str(key), str(value))
            output.write(line)
        output.write("==================================================\n")

else:
    pass
