import csv, collections

def load_targetscan(path):
    d = collections.defaultdict(set)
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            d[row["miRNA"]].add(row["Gene"])
    return dict(d)

def load_mirdb(path, score_min=80):
    d = collections.defaultdict(set)
    with open(path) as fh:
        for line in fh:
            mir, gene, score = line.rstrip("\n").split("\t")[:3]
            if float(score) >= score_min:
                d[mir].add(gene)
    return dict(d)

def target_set(mirna, targetscan, mirdb, mode="union"):
    a = targetscan.get(mirna, set())
    b = mirdb.get(mirna, set())
    return (a | b) if mode == "union" else (a & b)

def harmonize_mirna_names(name):
    base = name.split(".")[0]
    return base.lower()
