#!/usr/bin/env python3
"""
ComputeRkBound.py
=================
Explicit lower bound for R_k(n), k odd square-free, (n,k)=1, n >= 8e9, via the
divisor-sum form of Proposition 4.1 (the explicit Estermann estimate, for any
omega(k)):

  R_k(n)/n  >  C_Artin * prod_{q|k} lambda_q  -  (I) - (II) - (III),

  lambda_q = q(q-2)/(q^2-q-1),   and with c>=1, Z>c integers, A in (0,1/2):

  (I)   = (1/log n) sum_{a<=c,(a,k)=1} mu^2(a) sum_{j|k} c_theta(j a^2)
          [the j=a=1 modulus-1 term uses Broadbent: 0.375/(log n)^3]
  (II)  = ( prod_{q|k} (q-2)/(q-1) + 2/(1-2A) )
          ( sum_{c<a<=Z,(a,k)=1} mu^2(a)/phi(a^2) + 4/Z )
  (III) = ( k n^{-2A} + sqrt(k) n^{-A} + n^{A-1}/sqrt(k) ) log n.
          [the three pieces are the Sigma_4 tail (d > n^A/sqrt(k)) and the
           bad-gcd range (d,n)>1 contributing n^{A-1}/sqrt(k); see Prop 4.1 proof]

c_theta(m) is the SHARP per-modulus error bound for theta(n;m,a): for each modulus
the smallest available of
  - the raw Bennett et al. table entry   (bennett_c_theta.tsv; validity x0 in bennett_x0.tsv),
  - the Ramare-Rumely Table 1 entries    (rr_theta_table1.tsv; valid for n >= the listed x0),
  - the Ramare-Rumely Table 2 sqrt-bound (rr_theta_table2.tsv; valid for n <= 1e10),
  - the Bennett envelope  1/840 (m<=1e4), 1/160 (m<=1e5)   [fallback of last resort].
Using the raw tables, not the envelope, is what yields the constants in the paper
(e.g. R_15 > 0.09527 n, not 0.088 n). The modulus-1 term uses Broadbent et al.,
|theta(n) - n| <= 0.375 n / (log n)^3.

DATA FILES (REQUIRED).  Four tables are read from this script's own directory:
    bennett_c_theta.tsv   bennett_x0.tsv   rr_theta_table1.tsv   rr_theta_table2.tsv
They MUST be present. If any is missing or empty the script aborts, rather than
silently falling back to the (weaker) envelope -- which does not reproduce the
paper and can turn a valid bound negative (e.g. R_51: 0.05556 with tables, 0.0523
without).

Companion to ComputeBqBound.py. Pure Python 3 standard library; no third-party deps.

Usage:
    python3 ComputeRkBound.py --k 15                  # R_15(8e9)/n   (Corollary 4.2)
    python3 ComputeRkBound.py --k 105 --n 2.55e11     # evaluate at a chosen n
    python3 ComputeRkBound.py --k 429 --threshold     # least n with R_k(n) > 0
"""
import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from math import gcd

# SHARPEN: when True, Sigma_2 is bounded by Proposition 4.1's sharpened form --
#   (W) the Brun-Titchmarsh loss factor log n / log(n/(k d^2)) is kept INSIDE the
#       tail sum instead of being replaced by its endpoint value 1/(1-2A); and
#   (R) the coprimality condition (n-p,k)=1 is retained rather than dropped, so
#       Brun-Titchmarsh is applied to the modulus k d^2 over the prod_{q|k}(q-2)
#       admissible classes, contributing the factor alpha = prod (q-2)/(q-1),
#       at the cost of one exceptional term of size at most log k per d.
# Set SHARPEN = False to recover the published (unsharpened) constants exactly.
SHARPEN = True

C_ARTIN = 0.3739558136192023
BROADBENT_THETA_CONST = 0.375
BROADBENT_LOWER_X0 = math.exp(20.0)
RR_TABLE2_MAX_X = 1.0e10
DEFAULT_N = 8e9
MIN_VALID_N = 8e9
Z = 100_000
MAX_TABLE_MOD = 100_000
TABLE_FILES = ("bennett_c_theta.tsv", "bennett_x0.tsv", "rr_theta_table1.tsv", "rr_theta_table2.tsv")


def bennett_c0(modulus):
    return 1.0 / 840.0 if modulus <= 10_000 else 1.0 / 160.0


def normalize_header(header):
    return "".join(ch.lower() for ch in header.strip() if ch.isalnum())


def read_tsv_rows(path):
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8", errors="ignore") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        if not sample.strip():
            return []
        dialect = csv.excel_tab
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        except csv.Error:
            pass
        reader = csv.DictReader(fh, dialect=dialect)
        rows = []
        for raw in reader:
            row = {normalize_header(k): (v or "").strip() for k, v in raw.items() if k is not None}
            if any(row.values()):
                rows.append(row)
        return rows


def pick_value(row, aliases, required=True):
    for alias in aliases:
        if alias in row and row[alias] != "":
            return row[alias]
    if required:
        raise KeyError(f"Missing required column. Tried: {', '.join(aliases)}")
    return None


def parse_int(text):
    return int(float(text.replace("_", "")))


def parse_float(text):
    return float(text.replace("_", ""))


@dataclass(frozen=True)
class BennettBound:
    c_theta: float
    x_theta: int


@dataclass(frozen=True)
class RRTable2Bound:
    theta_sqrt_const: float


@dataclass(frozen=True)
class CandidateBound:
    source: str
    modulus: int
    contribution: float
    lower_x0: float


def require_tables(here):
    """Abort loudly if any required table is missing or empty, instead of
    silently degrading to the envelope (which does not reproduce the paper)."""
    problems = []
    for name in TABLE_FILES:
        path = here / name
        if not path.exists():
            problems.append(f"  missing: {name}")
        elif not read_tsv_rows(path):
            problems.append(f"  empty:   {name}")
    if problems:
        raise SystemExit(
            "ComputeRkBound.py: required c_theta table(s) not found next to the script:\n"
            + "\n".join(problems)
            + f"\n\nExpected all of {', '.join(TABLE_FILES)} in:\n  {here}\n"
            "Without them the bound falls back to the Bennett envelope, which is\n"
            "weaker and does not reproduce the paper's constants. Aborting."
        )


def load_bennett_c_theta(path):
    out = {}
    for row in read_tsv_rows(path):
        modulus = parse_int(pick_value(row, ["m", "modulus", "k", "q"]))
        c_theta = parse_float(pick_value(row, ["ctheta", "c_theta", "rawctheta", "raw_c_theta", "thetabound", "thetaepsilon"]))
        out[modulus] = c_theta
    return out


def load_bennett_x0(path):
    out = {}
    for row in read_tsv_rows(path):
        modulus = parse_int(pick_value(row, ["m", "modulus", "k", "q"]))
        x_theta = parse_int(pick_value(row, ["xtheta", "x_theta", "x0theta", "x0_theta", "xthetamin", "xthetah"]))
        out[modulus] = x_theta
    return out


def load_rr_table1(path):
    threshold_aliases = {
        1.0e10: ["eps1e10", "rr1e10", "1e10", "1010", "10pow10", "10to10"],
        1.0e13: ["eps1e13", "rr1e13", "1e13", "1013", "10pow13", "10to13"],
        1.0e30: ["eps1e30", "rr1e30", "1e30", "1030", "10pow30", "10to30"],
        1.0e100: ["eps1e100", "rr1e100", "1e100", "10100", "10pow100", "10to100"],
    }
    out = {}
    for row in read_tsv_rows(path):
        modulus = parse_int(pick_value(row, ["k", "modulus", "m", "q"]))
        entries = sorted(
            [(x0, parse_float(v)) for x0, aliases in threshold_aliases.items()
             if (v := pick_value(row, aliases, required=False)) is not None],
            key=lambda t: t[0],
        )
        if entries:
            out[modulus] = entries
    return out


def load_rr_table2(path):
    out = {}
    for row in read_tsv_rows(path):
        modulus = parse_int(pick_value(row, ["k", "modulus", "m", "q"]))
        theta = parse_float(pick_value(row, ["theta", "rrtheta", "thetasqrt", "thetaconst"]))
        out[modulus] = RRTable2Bound(theta_sqrt_const=theta)
    return out


def mobius_sieve(n):
    mu = [1] * (n + 1)
    mu[0] = 0
    primes, is_comp = [], [False] * (n + 1)
    for i in range(2, n + 1):
        if not is_comp[i]:
            primes.append(i)
            mu[i] = -1
        for p in primes:
            ip = i * p
            if ip > n:
                break
            is_comp[ip] = True
            if i % p == 0:
                mu[ip] = 0
                break
            mu[ip] = -mu[i]
    return mu


def spf_sieve(n):
    spf = list(range(n + 1))
    for i in range(2, int(n**0.5) + 1):
        if spf[i] == i:
            for j in range(i*i, n+1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def phi_sieve(n):
    phi = list(range(n + 1))
    for p in range(2, n + 1):
        if phi[p] == p:
            for j in range(p, n + 1, p):
                phi[j] -= phi[j] // p
    return phi


def phi_square(a, spf):
    res, x, last = a * a, a, 0
    while x > 1:
        p = spf[x]
        if p != last:
            res -= res // p
            last = p
        while x % p == 0:
            x //= p
    return res


def factor_squarefree(k):
    if k % 2 == 0:
        raise ValueError(f"k={k} must be odd.")
    primes, x, p = [], k, 3
    while p * p <= x:
        if x % p == 0:
            primes.append(p)
            x //= p
            if x % p == 0:
                raise ValueError(f"k={k} is not squarefree.")
        p += 2
    if x > 1:
        primes.append(x)
    return sorted(primes)


def divisors_from_primes(primes):
    divs = [1]
    for p in primes:
        divs += [d * p for d in divs]
    return sorted(divs)


def alpha_coeff(primes):
    a = 1.0
    for p in primes:
        a *= (p - 2) / (p - 1)
    return a


def beta_coeff(primes):
    b = 1.0
    for p in primes:
        b *= (p * (p - 2)) / (p*p - p - 1)
    return b


def collect_required_terms(k, c, mu):
    divs = divisors_from_primes(factor_squarefree(k))
    return [
        (a, d, d * a * a)
        for a in range(1, c + 1) if gcd(a, k) == 1 and mu[a] != 0
        for d in divs
    ]


def rr_table1_candidate(n, modulus, phi_m, rr_table1):
    entries = rr_table1.get(modulus)
    if not entries:
        return None
    valid = [(x0, eps) for x0, eps in entries if n >= x0]
    if not valid:
        return None
    x0, eps = max(valid, key=lambda t: t[0])
    return CandidateBound(source=f"RR-Table1(x0={x0:.0e})", modulus=modulus, contribution=eps / phi_m, lower_x0=x0)


def rr_table2_candidate(n, modulus, rr_table2):
    entry = rr_table2.get(modulus)
    if entry is None or n > RR_TABLE2_MAX_X:
        return None
    return CandidateBound(source="RR-Table2", modulus=modulus, contribution=entry.theta_sqrt_const / math.sqrt(n), lower_x0=0.0)


def bennett_candidate(n, modulus, bennett):
    entry = bennett.get(modulus)
    if entry is None or n < entry.x_theta:
        return None
    return CandidateBound(source="Bennett-raw", modulus=modulus, contribution=entry.c_theta / math.log(n), lower_x0=float(entry.x_theta))


def bennett_c0_candidate(n, modulus):
    if not (3 <= modulus <= MAX_TABLE_MOD) or n < MIN_VALID_N:
        return None
    return CandidateBound(source="Bennett-c0", modulus=modulus, contribution=bennett_c0(modulus) / math.log(n), lower_x0=MIN_VALID_N)


def choose_best_bound(n, modulus, phi_m, bennett, rr_table1, rr_table2):
    candidates = [c for c in [
        bennett_candidate(n, modulus, bennett),
        rr_table2_candidate(n, modulus, rr_table2),
        rr_table1_candidate(n, modulus, phi_m, rr_table1),
        bennett_c0_candidate(n, modulus),
    ] if c is not None]
    if not candidates:
        raise ValueError(f"No valid bound for modulus m={modulus} at n={n:.0f}.")
    return min(candidates, key=lambda c: c.contribution)


def compute_first_sum(n, k, c, mu, phi, bennett, rr_table1, rr_table2):
    total = broadbent_total = non_m1_total = max_x0 = 0.0
    for a, d, modulus in collect_required_terms(k, c, mu):
        if modulus == 1:
            contrib = BROADBENT_THETA_CONST / math.log(n)**3
            total += contrib
            broadbent_total += contrib
            max_x0 = max(max_x0, BROADBENT_LOWER_X0)
        else:
            best = choose_best_bound(n, modulus, phi[modulus], bennett, rr_table1, rr_table2)
            total += best.contribution
            non_m1_total += best.contribution
            max_x0 = max(max_x0, best.lower_x0)
    return total, broadbent_total, non_m1_total, max_x0


def tail_constant(k, c, mu, spf):
    return sum(
        1.0 / phi_square(a, spf)
        for a in range(c + 1, Z + 1) if gcd(a, k) == 1 and mu[a] != 0
    ) + 4.0 / Z


def endpoint_term(n, k, c, mu):
    """The artificial p = n endpoint of theta(n; e d^2, n), in units of n.

    For k > 1 the endpoint contributes log n * sum_{e|k} mu(e) = 0 to the signed
    combination and cancels exactly.  For k = 1 there is nothing to cancel
    against, so one term of weight at most log n survives for each admissible
    d <= c and must be carried.  It is tiny (about 5e-7 at n = 8e9) but the
    comparison device of Section 5 subtracts a k = 1 bound, so it is stated
    rather than assumed away.
    """
    if k > 1:
        return 0.0
    cnt = sum(1 for d in range(1, c + 1) if mu[d])
    return cnt * math.log(n) / n


_SUFFIX = {}


def sieve_suffix(k, mu, spf):
    """S[y] = sum_{y < d <= Z, (d,k)=1} mu^2(d)/phi(d^2), cached per k.

    Used to bound the Brun-Titchmarsh range ABOVE the evaluation point sharply,
    instead of leaning on Ramare's 4/y across the whole tail: the exact sum is
    available for d <= Z and only the genuine remainder beyond Z needs 4/Z.
    """
    if k not in _SUFFIX:
        S = [0.0] * (Z + 2)
        run = 0.0
        for d in range(Z, 0, -1):
            S[d] = run
            if mu[d] and gcd(d, k) == 1:
                run += 1.0 / phi_square(d, spf)
        S[0] = run
        _SUFFIX[k] = S
    return _SUFFIX[k]


def bt_prefix(n, k, c, mu, spf):
    """Prefix sums of the d-dependent Brun-Titchmarsh weight over the Sigma_2 range.

    P[D] = sum_{c < d <= D, (d,k)=1, mu(d) != 0} w(d) / phi(d^2),
        w(d) = log n / log(n / (k d^2)),
    which is the exact loss factor from applying Brun-Titchmarsh to the modulus
    k d^2 at height n.  At the endpoint d = n^A / sqrt(k) one has k d^2 = n^{2A}
    and w(d) = 1/(1-2A), the uniform factor used in the unsharpened bound; for
    every smaller d the weight is strictly smaller.  Entries with k d^2 >= n are
    left at infinity: they are never read, since the Sigma_2 range stops at
    d = n^A / sqrt(k) and A < 1/2.

    Precomputing this once per (n, k) keeps the search over A an O(1) lookup.
    """
    logn = math.log(n)
    P = [0.0] * (Z + 2)
    run = 0.0
    for d in range(1, Z + 1):
        if d > c and mu[d] and gcd(d, k) == 1:
            denom = logn - math.log(k * d * d)
            if denom <= 0.0:
                run = math.inf
            else:
                run += (logn / denom) / phi_square(d, spf)
        P[d] = run
    P[Z + 1] = run
    return P


def err_bound(n, A, k, alpha, first_sum_total, tail, prefix=None, suffix=None):
    logn = math.log(n)
    if SHARPEN:
        # (II): the completed principal tail (coefficient alpha) plus the
        # Brun-Titchmarsh range, which after (R) also carries the factor alpha.
        # Sigma_2 runs over c < d <= D(n) = n^A/sqrt(k), which GROWS with n.
        # The bound must hold for every n in [n0, infinity), so the sum is
        # evaluated at the interval's left endpoint n0 (= n here) and the range
        # above it is absorbed by a Ramare tail:
        #     sum_{c<d<=D(n)} mu^2 w_n / phi   <=   sum_{c<d<=Y} mu^2 w_{n0}/phi
        #                                          + (1/(1-2A)) * 4/Y,
        # where Y = D(n0).  This is legitimate because w_n(d) is DECREASING in n
        # at fixed d, so w_{n0} dominates on d <= Y, while on Y < d <= D(n) the
        # weight never exceeds its endpoint value w(D) = 1/(1-2A).  Without this
        # remainder the certification would not propagate past the n at which
        # D(n) overtakes the tabulated range.
        D = int(math.floor(n ** A / math.sqrt(k)))
        if D < 1 or D > Z or not math.isfinite(prefix[D]):
            return math.inf                      # parameter choice inadmissible
        ramare_bt = (suffix[D] + 4.0 / Z) / (1.0 - 2.0 * A)
        term_II = alpha * (tail + 2.0 * (prefix[D] + ramare_bt))
        # (III): Sigma_4 as before; the bad-gcd range now also absorbs the
        # exceptional-prime term of (R), at most sum_{q|k} log q = log k per d.
        sigma4 = (k * n ** (-2.0 * A) + math.sqrt(k) * n ** (-A)) * logn
        E_gcd = (n ** (A - 1.0) / math.sqrt(k)) * (logn + math.log(k))
        return first_sum_total + term_II + sigma4 + E_gcd
    # (III): the two trivial ranges of the Moebius sieve.
    #   * Sigma_4 (d > n^A/sqrt(k)):    (k n^{-2A} + sqrt(k) n^{-A}) log n
    #   * Sigma_3 (bad gcd, (d,n)>1):   n^{A-1}/sqrt(k) log n  =: E_gcd
    # The Sigma_3 contribution is bounded in the proof by |Sigma_3| <= (n^A/sqrt(k)) log n
    # and must be carried explicitly; it is not absorbed into the Sigma_4 bound (the
    # -1/sqrt(n) and +sqrt(n) boundary pieces there cancel, leaving Sigma_4 alone).
    sigma4 = (k * n**(-2.0 * A) + math.sqrt(k) * n**(-A)) * logn
    E_gcd = (n**(A - 1.0) / math.sqrt(k)) * logn
    return (
        first_sum_total
        + (alpha + 2.0 / (1.0 - 2.0 * A)) * tail
        + sigma4
        + E_gcd
    )


def minimize_A(n, k, alpha, first_sum_total, tail, prefix=None, suffix=None):
    f = lambda A: err_bound(n, A, k, alpha, first_sum_total, tail, prefix, suffix)
    best_i = min(range(500, 4900), key=lambda i: f(i/10000))
    bestA = best_i / 10000.0
    lo, hi = max(1e-6, bestA - 0.002), min(0.499999, bestA + 0.002)
    bestA = min((lo + t*(hi-lo)/8000 for t in range(8001)), key=f)
    return bestA, f(bestA)


def evaluate(n, k, c, mu, phi, bennett, rr_table1, rr_table2, tail, alpha, beta, spf=None):
    total, broadbent_total, non_m1_total, x0_req = compute_first_sum(n, k, c, mu, phi, bennett, rr_table1, rr_table2)
    total += endpoint_term(n, k, c, mu)
    prefix = bt_prefix(n, k, c, mu, spf) if SHARPEN else None
    suffix = sieve_suffix(k, mu, spf) if SHARPEN else None
    A_star, err = minimize_A(n, k, alpha, total, tail, prefix, suffix)
    # Proposition 4.1 requires Z >= n^{A}/sqrt(k) (the tail truncation must reach the
    # trivial-range cutoff). Z is fixed at MAX_TABLE_MOD; assert the hypothesis holds
    # at the optimised A rather than return a bound whose derivation is unjustified.
    # Raising Z would only decrease (II), so a failure here means Z must be increased,
    # never that the bound is unsafe -- but we refuse to certify silently either way.
    z_required = n ** A_star / math.sqrt(k)
    if z_required > Z:
        raise AssertionError(
            f"Prop 4.1 hypothesis violated for k={k}, n={n:.3e}: need Z >= "
            f"n^A/sqrt(k) = {z_required:.3e} but Z = {Z}. Increase Z (MAX_TABLE_MOD).")
    return {
        "bound": beta * C_ARTIN - err,
        "A_star": A_star, "err": err, "first_sum": total,
        "broadbent": broadbent_total, "non_m1": non_m1_total, "x0_req": x0_req,
    }


def find_threshold(k, c, mu, phi, bennett, rr_table1, rr_table2, tail, alpha, beta, spf=None, nmax=1e13):
    f = lambda n: evaluate(n, k, c, mu, phi, bennett, rr_table1, rr_table2, tail, alpha, beta, spf)["bound"]
    lo, hi = MIN_VALID_N, nmax
    if f(lo) > 0:
        return lo
    if f(hi) <= 0:
        return None
    for _ in range(200):
        if hi / lo <= 1.0005:
            break
        mid = math.sqrt(lo * hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    return hi


def main():
    parser = argparse.ArgumentParser(description="Explicit lower bound for R_k(n) (Proposition 4.1)")
    parser.add_argument("--k", type=int, required=True, help="odd square-free modulus k")
    parser.add_argument("--n", type=float, default=DEFAULT_N, help="evaluation point (default 8e9)")
    parser.add_argument("--threshold", action="store_true", help="also report the least n with R_k(n) > 0")
    args = parser.parse_args()

    k, n = args.k, args.n
    if n < MIN_VALID_N:
        raise SystemExit(f"n = {n:.3e} is below the validity floor {MIN_VALID_N:.0e}; the imported bounds require n >= 8e9.")

    here = Path(__file__).resolve().parent
    require_tables(here)

    bennett_c_theta = load_bennett_c_theta(here / "bennett_c_theta.tsv")
    bennett_x0 = load_bennett_x0(here / "bennett_x0.tsv")
    bennett = {m: BennettBound(c_theta=bennett_c_theta[m], x_theta=bennett_x0[m])
               for m in bennett_c_theta if m in bennett_x0}
    rr_table1 = load_rr_table1(here / "rr_theta_table1.tsv")
    rr_table2 = load_rr_table2(here / "rr_theta_table2.tsv")

    mu = mobius_sieve(Z)
    spf = spf_sieve(Z)
    phi = phi_sieve(MAX_TABLE_MOD)

    primes = factor_squarefree(k)
    alpha = alpha_coeff(primes)
    beta = beta_coeff(primes)
    c = int(math.floor(math.sqrt(MAX_TABLE_MOD / k)))
    tail = tail_constant(k, c, mu, spf)

    r = evaluate(n, k, c, mu, phi, bennett, rr_table1, rr_table2, tail, alpha, beta, spf)

    print(f"k={k}, primes={primes}, c={c}, n={n:.3e}, x0_req={r['x0_req']:.0e}")
    print(f"alpha = {alpha:.12f},  beta = {beta:.12f}")
    print(f"first_sum = {r['first_sum']:.12f}  (broadbent = {r['broadbent']:.12f},  non_m1 = {r['non_m1']:.12f})")
    print(f"tail  = {tail:.12f}")
    print(f"A*    = {r['A_star']:.7f},  err = {r['err']:.12f}")
    print(f"R_k(n)/n >= {r['bound']:.12f}   {'>0' if r['bound'] > 0 else '<=0'}")

    if args.threshold:
        t = find_threshold(k, c, mu, phi, bennett, rr_table1, rr_table2, tail, alpha, beta, spf)
        if t is None:
            print("threshold: R_k(n) > 0 not established below 1e13")
        else:
            print(f"threshold: R_k(n) > 0 for n >= {t:.3e}  (log10 {math.log10(t):.2f};"
                  f" within 2.55e11: {t <= 2.55e11})")


if __name__ == "__main__":
    main()
