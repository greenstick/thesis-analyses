#! /home/users/cordier/.linuxbrew/bin/python3

# Imports & Methods
from collections import OrderedDict
import operator, hashlib
import vcf as VCF

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
# File Writing
#

# Write Python List of Dictionaries to Directory as JSON File
def writeJSON (data, path="output.json"):
    """
    Write JSON output to file.
    @params:
        data        - Required  : data to write (Dict)
        path        - Optional  : output path (Str)
    Returns: path (Str)
    """
    import json
    with open(path, "w") as file:
        file.write(json.dumps(data, indent = 4, separators=(',', ': '), sort_keys = False))
        file.close()
    return path

#
# UI
#

# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()

#
# Main Routine
#

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

    truthsets = [
        "synthetic.challenge.set1.tumor.all.truth.vcf",
        "synthetic.challenge.set2.tumor.all.truth.vcf",
        "synthetic.challenge.set3.tumor.20pctmasked.truth.vcf"
    ]
    matrices = []
    uniqueProps = ["CHROM", "POS", "REF", "ALT"] # Variant Record Properties to Define Uniqueness

    for dataset in metadata["datasets"]:
        # Set Active Truth & Open VCF
        truthActive = [truthset for truthset in truthsets if dataset["key"] in truthset][0]
        with open("vcf/truth-set/%s" % truthActive, 'r') as truthHandle:
            print("Status: Hashing Truth Set – %s" % truthActive)
            truthReader = VCF.Reader(truthHandle)
            truthSet = hashVariants(truthReader, uniqueProps, algorithm = "sha1").keys()
            i, total = 0, len(metadata["conditions"]) * len(metadata["subconditions"])
            for condition in metadata["conditions"]:
                for subcondition in metadata["subconditions"]:
                    # Progress Bar
                    printProgressBar(i, total, prefix = "Computing Matrix:", suffix = "%s - %s - %s" % (dataset["name"], condition["name"], subcondition["name"]))
                    # Set Active Condition & Open VCF
                    conditionActive = "vcf/synthetic.challenge.%s.%s.default.%s.raw.snps.indels.vcf" % (dataset["key"], condition["key"], subcondition["key"])
                    with open(conditionActive) as conditionHandle:
                        conditionReader = VCF.Reader(conditionHandle)
                        conditionHashMap = hashVariants(conditionReader, uniqueProps, algorithm = "sha1")
                        analyses = generateAnalyses(conditionHashMap, truthSet, metricKey = "TLOD", thresholds = 100, skip = 1)
                        label = dataset["name"] + " \n " + condition["name"] + " \n " + subcondition["name"]
                        metakey = dataset["key"] + condition["key"] + subcondition["key"]
                        matrices.append({
                            "label"     : label,
                            "key"       : metakey,
                            "analyses"  : analyses
                        })
                    i += 1
                    printProgressBar(i, total, prefix = "Computing Matrices:", suffix = "%s - %s - %s" % (dataset["name"], condition["name"], subcondition["name"]))
    writeJSON(matrices, "heatmaps.json")

else:
    pass

