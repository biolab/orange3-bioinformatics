""" Documentation script

Check all data files from GEO, find those which include at least N
samples in all sample subsets of at least one sample type. Useful
when, for instance, filtering out the data sets that could be used for
supervised machine learning.
"""

from orangecontrib.bioinformatics.geo.dataset import GDSInfo, GDS


def valid(info, n=40):
    """Return a set of subset types containing more than n samples in every subset"""
    invalid = set()
    subsets = set([sinfo["type"] for sinfo in info["subsets"]])
    for sampleinfo in info["subsets"]:
        if len(sampleinfo["sample_id"]) < n:
            invalid.add(sampleinfo["type"])
    return subsets.difference(invalid)


def report(stypes, info):
    """Pretty-print GDS and valid susbset types"""
    for id, sts in stypes:
        print(id)
        for st in sts:
            gds = info[id]
            print("  %s:" % st +
                  ", ".join(["%s/%d" % (sinfo["description"], len(sinfo["sample_id"]))
                             for sinfo in gds["subsets"] if sinfo["type"] == st]))



gdsinfo = GDSInfo()
valid_subset_types = [(id, valid(info)) for id, info in sorted(gdsinfo.items()) if valid(info)]
report(valid_subset_types, gdsinfo)

print('datasets = ' + str(len(valid_subset_types)))
print('type subsets = ' + str(sum(len(b) for _, b in valid_subset_types)))
