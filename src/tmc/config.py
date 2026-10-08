"""Loaders for the frozen experiment map (config/accessions.tsv)."""
import csv

CONTROL = "dmso_control"
# The nine PRTR Class-I chemicals (benzene first: the only known carcinogen).
CHEMICALS = ["benzene", "dehp", "chloroform", "p_cresol", "p_dcb",
             "phenol", "pyrocatechol", "tnbp", "tce"]


def load_accessions(path):
    """Return the accession<->chemical<->layer map as a list of dict rows.

    Columns: layer, run_accession, chemical, replicate, fastq_ftp, fastq_md5.
    """
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))
