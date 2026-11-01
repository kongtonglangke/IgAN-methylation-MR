"""Stage 2a: type-I error simulation for candidate pathway-level statistics.

Question
--------
The 27 metabolite-target genes carry ~1 700 usable CpGs, but neighbouring
cis-SNPs sit in tight LD, so the per-CpG Wald z-scores are correlated. Which
aggregate statistic keeps a nominal 5 % error rate when that correlation is
taken into account only by the textbook independence assumption, and which one
collapses?

Statistics
----------
  minp  : Sidak-corrected minimum p            p = 1-(1-p_min)^k
  fisher: Fisher combination                   -2*sum ln p ~ chi2(2k)
  chisq : k * mean(z^2)                        ~ chi2(k)
  sumz  : sum(z)/sqrt(k)                       ~ N(0,1)
  pct   : Pathway Concordance Test -- the pathway mean z is standardised by a
          matched *control* CpG set instead of by sqrt(k):
              z_pct = (mean_p - mean_c) / sqrt(sd_p^2/k + sd_c^2/k)
          so whatever LD is present enters both arms and cancels.

Correlation model
-----------------
Compound symmetry inside each gene block, independence between blocks:
    z_i = sqrt(rho)*g_block + sqrt(1-rho)*e_i
which is the exact Cholesky form for equicorrelated z, so rho can be pushed to
0.95 without matrix algebra.
"""
import math
import random

random.seed(20261030)

BLOCK = 60          # fallback CpGs per gene block

# observed number of usable (deduplicated) CpGs per gene in the real data:
# sizes range from 2 (LIF) to 320 (MADCAM1), i.e. strongly unbalanced clusters
GENE_SIZES = [2, 10, 14, 21, 29, 30, 31, 34, 35, 36, 36, 37, 38, 49, 54, 61,
              77, 87, 93, 101, 104, 194, 206, 320]


# ---------------------------------------------------------------- distributions
def norm_sf(x):
    return 0.5 * math.erfc(x / math.sqrt(2))


def _gammainc_upper_reg(a, x):
    """Regularised upper incomplete gamma Q(a, x) (Numerical Recipes gammq)."""
    if x < 0 or a <= 0:
        return float("nan")
    if x == 0:
        return 1.0
    gln = math.lgamma(a)
    if x < a + 1.0:
        # series for P(a,x)
        ap, s, de = a, 1.0 / a, 1.0 / a
        for _ in range(1000):
            ap += 1.0
            de *= x / ap
            s += de
            if abs(de) < abs(s) * 1e-14:
                break
        return max(0.0, 1.0 - s * math.exp(-x + a * math.log(x) - gln))
    # continued fraction for Q(a,x)
    b = x + 1.0 - a
    c = 1e300
    d = 1.0 / b
    h = d
    for i in range(1, 1000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < 1e-300:
            d = 1e-300
        c = b + an / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < 1e-14:
            break
    return min(1.0, max(0.0, math.exp(-x + a * math.log(x) - gln) * h))


def chi2_sf(x, df):
    if x <= 0:
        return 1.0
    return _gammainc_upper_reg(df / 2.0, x / 2.0)


def _betacf(a, b, x):
    """Continued fraction for the incomplete beta function (Lentz)."""
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < 1e-14:
            break
    return h


def _betainc_reg(a, b, x):
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    front = math.exp(lbeta + a * math.log(x) + b * math.log(1 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def _t_sf2(t, df):
    """Two-sided p for Student t."""
    p = _betainc_reg(df / 2.0, 0.5, df / (df + t * t))
    return min(1.0, max(0.0, p))


# ------------------------------------------------------------------- simulation
def block(n, rho, shift=0.0):
    g = random.gauss(0.0, 1.0) + shift
    a = math.sqrt(max(rho, 0.0))
    b = math.sqrt(max(1.0 - rho, 1e-12))
    return [a * g + b * random.gauss(0.0, 1.0) for _ in range(n)]


def draw(k, rho, shift=0.0, sizes=None):
    """Draw k correlated z-scores. With `sizes`, each list entry is one gene
    cluster of that size (unbalanced, as in the real data); otherwise fixed
    blocks of BLOCK."""
    if sizes is None:
        out = []
        while len(out) < k:
            n = min(BLOCK, k - len(out))
            out.extend(block(n, rho, shift))
        return out
    out = []
    for n in sizes:
        out.extend(block(n, rho, shift))
    return out


def two_sided_p(z):
    if abs(z) > 8.5:
        return 0.0
    return 2.0 * norm_sf(abs(z))


def stats(z, ctrl, sizes=None):
    k = len(z)
    p = [two_sided_p(x) for x in z]
    p = [max(x, 1e-300) for x in p]
    out = {}
    out["minp"] = 1.0 - (1.0 - min(p)) ** k
    out["fisher"] = chi2_sf(-2.0 * sum(math.log(x) for x in p), 2 * k)
    out["chisq"] = chi2_sf(sum(x * x for x in z), k)
    out["sumz"] = two_sided_p(sum(z) / math.sqrt(k))
    mp = sum(z) / k
    mc = sum(ctrl) / len(ctrl)
    vp = sum((x - mp) ** 2 for x in z) / max(k - 1, 1)
    vc = sum((x - mc) ** 2 for x in ctrl) / max(len(ctrl) - 1, 1)
    se = math.sqrt(vp / k + vc / len(ctrl))
    out["pct"] = two_sided_p((mp - mc) / se) if se > 0 else 1.0
    # cluster-robust variant: gene is the cluster, SE taken from the spread of
    # the gene means, t reference with G-1 df
    gm = []
    if sizes is not None:
        pos = 0
        for n in sizes:
            seg = z[pos:pos + n]
            pos += n
            gm.append(sum(seg) / len(seg))
    else:
        for i in range(0, k, BLOCK):
            seg = z[i:i + BLOCK]
            gm.append(sum(seg) / len(seg))
    G = len(gm)
    m = sum(gm) / G
    if G > 2:
        sd = math.sqrt(sum((x - m) ** 2 for x in gm) / (G - 1))
        t = m / (sd / math.sqrt(G)) if sd > 0 else 0.0
        out["pctcr"] = _t_sf2(abs(t), G - 1)
    else:
        out["pctcr"] = 1.0
    return out


KEYS = ("minp", "fisher", "chisq", "sumz", "pct", "pctcr")


def run(rho, sizes, reps, alpha=(0.05, 0.01)):
    k = sum(sizes)
    rej = {s: {a: 0 for a in alpha} for s in KEYS}
    for _ in range(reps):
        z = draw(k, rho, sizes=sizes)
        c = draw(k, rho, sizes=sizes)
        r = stats(z, c, sizes=sizes)
        for s in KEYS:
            for a in alpha:
                if r[s] < a:
                    rej[s][a] += 1
    return {s: {a: rej[s][a] / reps for a in alpha} for s in KEYS}


def power(rho, sizes, frac, delta, reps, alpha=0.05):
    k = sum(sizes)
    G = len(sizes)
    sig = set(random.sample(range(G), max(1, int(round(G * frac)))))
    hit = {s: 0 for s in KEYS}
    for _ in range(reps):
        z = []
        for g, n in enumerate(sizes):
            z.extend(block(n, rho, delta if g in sig else 0.0))
        c = draw(k, rho, sizes=sizes)
        r = stats(z, c, sizes=sizes)
        for s in KEYS:
            if r[s] < alpha:
                hit[s] += 1
    return {s: hit[s] / reps for s in KEYS}


if __name__ == "__main__":
    print("=== A. type-I error rate, nominal 0.05 / 0.01 (2000 replicates) ===")
    print("  gene sizes taken from the real data (24 genes, 2-320 CpGs each),")
    print("  rho = within-gene correlation between cis-SNP instruments")
    print("%-5s %-6s | " % ("rho", "k") + " | ".join("%-13s" % s for s in KEYS))
    for scale in (1, 3):
        sizes = GENE_SIZES * scale
        for rho in (0.0, 0.3, 0.6, 0.9):
            r = run(rho, sizes, 2000)
            print("%-5.2f %-6d | " % (rho, sum(sizes)) + " | ".join(
                "%.3f / %.3f" % (r[s][0.05], r[s][0.01]) for s in KEYS), flush=True)

    print("\n=== B. power at alpha=0.05, rho=0.6, 24 genes (1000 replicates) ===")
    print("  a fraction of genes carries a per-CpG z shift delta")
    print("%-6s %-7s %-7s | " % ("frac", "delta", "OR") +
          " | ".join("%-9s" % s for s in KEYS))
    for frac in (0.3, 1.0):
        for delta in (0.02, 0.05, 0.10):
            pw = power(0.6, GENE_SIZES, frac, delta, 1000)
            print("%-6.2f %-7.2f %-7.3f | " % (frac, delta, math.exp(delta)) +
                  " | ".join("%-9.3f" % pw[s] for s in KEYS), flush=True)
