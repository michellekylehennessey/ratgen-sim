from itertools import product
from collections import Counter
import random

# Genes modeled: Agouti (A dominant over a), Pink-eyed dilution (P dominant over p)
GENES = ["A", "P"]

def parse_parent(genostr: str):
    """
    Parse 'A/a; P/p' or 'AA; Pp' etc. Returns dict like {'A': ('A','a'), 'P': ('P','p')}
    """
    geno = {}
    parts = [p.strip() for p in genostr.replace("|",";").split(";") if p.strip()]
    if len(parts) != len(GENES):
        raise ValueError(f"Expected {len(GENES)} loci separated by ';' (e.g. 'A/a; P/p').")
    for gene, part in zip(GENES, parts):
        part = part.replace("/", "").strip()
        if len(part) == 2:
            alleles = (part[0], part[1])
        else:
            a, b = [x.strip() for x in part.split("/")[:2]]
            alleles = (a, b)
        geno[gene] = alleles
    return geno

def gametes(parent):
    """Return all gametes (one allele per locus) with equal probability."""
    combos = list(product(*[parent[g] for g in GENES]))
    probs = 1 / len(combos)
    return [{"A": a, "P": p, "prob": probs} for a, p in combos]

def combine(g1, g2):
    """Combine two gametes -> offspring genotype dict like {'A':('A','a'),'P':('P','p')}."""
    return {g: tuple(sorted([g1[g], g2[g]])) for g in GENES}

def phenotype_from_genotype(gt):
    A = gt["A"]
    P = gt["P"]
    A_dom = "A" in A
    P_dom = "P" in P
    if A_dom and P_dom:
        return "Agouti"
    if (not A_dom) and P_dom:
        return "Black"
    if A_dom and (not P_dom):
        return "Silver Fawn"
    return "Beige"  # aa pp

def cross(parent1: str, parent2: str, litter_size: int | None = None, seed: int | None = None):
    """
    Returns:
      - 'table': dict of genotype -> probability
      - 'phenotypes': dict of phenotype -> probability
      - 'litter': list of phenotypes if litter_size provided
    """
    p1 = parse_parent(parent1)
    p2 = parse_parent(parent2)
    g1 = gametes(p1)
    g2 = gametes(p2)

    # Punnett probabilities
    geno_counter = Counter()
    for x in g1:
        for y in g2:
            gt = combine(x, y)
            key = f"A{''.join(gt['A'])};P{''.join(gt['P'])}"
            geno_counter[key] += x["prob"] * y["prob"]

    # Phenotype probabilities
    pheno_counter = Counter()
    rep_map = {}
    for key, prob in geno_counter.items():
        gt = {"A": tuple(key.split(";")[0][1:]), "P": tuple(key.split(";")[1][1:])}
        ph = phenotype_from_genotype(gt)
        pheno_counter[ph] += prob
        rep_map[key] = ph

    result = {
        "table": dict(sorted(geno_counter.items(), key=lambda x: -x[1])),
        "phenotypes": dict(sorted(pheno_counter.items(), key=lambda x: -x[1])),
    }

    # Optional litter simulation
    if litter_size:
        if seed is not None:
            random.seed(seed)
        phenos, weights = zip(*result["phenotypes"].items())
        cum = []
        total = 0
        for w in weights:
            total += w
            cum.append(total)
        litter = []
        for _ in range(litter_size):
            r = random.random()
            for ph, c in zip(phenos, cum):
                if r <= c:
                    litter.append(ph)
                    break
        result["litter"] = litter
    return result
