"""Audit punchline: how much of a predicted-target-based GO/pathway claim is
actually supported by observed repression in the matched chemical's mRNA layer."""


def observed_support_rate(claim_genes, repressed_genes):
    """Fraction of the claim's genes that are actually repressed (observed down)
    in the matched chemical. NaN when the claim set is empty."""
    claim = set(claim_genes)
    if not claim:
        return float("nan")
    return len(claim & set(repressed_genes)) / len(claim)
