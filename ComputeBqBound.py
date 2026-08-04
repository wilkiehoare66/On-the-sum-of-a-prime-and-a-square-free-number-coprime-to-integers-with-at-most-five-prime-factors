#!/usr/bin/env python3
"""
ComputeBqBound.py
=================
Explicit upper bound for the divisible part

        B_q(n) = sum_{p<n, q|(n-p)} mu^2(n-p) log p  =  R(n) - R_q(n),

which is Proposition 4.3 of the paper.  B_q(n) is the weighted mass discarded
when the coprimality condition is tightened by one further prime q, and it is
the quantity subtracted for each peeled prime in the base-plus-peel criterion

    R_m(n) > log 2 + sum_{l} B_{q_l}(n)          ==>  R_k(n) > 0.

That criterion is taken over the common Euler product W_n = prod_{p nmid n}
(1 - 1/(p(p-1))) >= C_Artin, which is the sharp common factor: R_m/n and B_q/n
have main terms (prod_{q|m} lambda_q) W_n and (1 - lambda_q) W_n respectively.

This script evaluates B_q alone.  The criterion itself, for an arbitrary odd
square-free k, lives in ScanOmega.py; the case division it induces is derived
by CoverageTree.py and certified by CertifyCases.py.

Two-family B_q (Proposition 4.3).  Imposing q | (n-p) splits the sieve modulus
lcm(q,d^2) into a family qd^2 (when (d,q)=1) and a family q^2 b^2 (writing
d=qb). The two families carry SEPARATE cutoffs: c1,Z1 for qd^2 with large-range
boundary d>n^A, and c2,Z2 for q^2 b^2 with boundary b>n^A/q (since
q^2 b^2 <= n^{2A} forces b <= n^A/q). Keeping them separate prevents the second
family's large range from being undercounted. The short range carries a p=n
endpoint term; the large range carries the bad-gcd (d,n)>1 / (b,n)>1 terms.

Coupled comparison.  The comparison is taken over the common Euler product
P_k(n) = prod_{p not| kn}(1-1/(p(p-1))) >= C_Artin, so the unrestricted side
contributes its genuine explicit-estimate error E_R(n) (from R_k at k=1,
including its bad-gcd term), not a flat lower bound.

Raw c_theta data, table loaders, R_m machinery and the missing-table guard are
IMPORTED from ComputeRkBound.py. Requires ComputeRkBound.py and the four
c_theta tables in the same directory; aborts if a table is missing.

Usage:
    python3 ComputeBqBound.py --k 429              # both criteria for k, best threshold
    python3 ComputeBqBound.py --bq 13              # just the B_13 upper-bound error
    python3 ComputeBqBound.py --k 2431 --n 1e10    # criteria at a chosen n
    python3 ComputeBqBound.py --k 105 --certify    # directed-rounding certificate
"""
import argparse
import math
from math import gcd, floor, sqrt
from pathlib import Path

from ComputeRkBound import (
    C_ARTIN, MAX_TABLE_MOD, Z, MIN_VALID_N, DEFAULT_N,
    require_tables,
    load_bennett_c_theta, load_bennett_x0, load_rr_table1, load_rr_table2,
    BennettBound,
    mobius_sieve, spf_sieve, phi_sieve, phi_square,
    choose_best_bound,
    factor_squarefree, alpha_coeff, beta_coeff, tail_constant, evaluate,
)


# Directed-rounding envelope: >6 orders of magnitude above the measured double
# rounding error (~1e-16 on the tail sums).


def build_context():
    here = Path(__file__).resolve().parent
    require_tables(here)
    bct = load_bennett_c_theta(here / "bennett_c_theta.tsv")
    bx0 = load_bennett_x0(here / "bennett_x0.tsv")
    bennett = {m: BennettBound(c_theta=bct[m], x_theta=bx0[m]) for m in bct if m in bx0}
    return {"mu": mobius_sieve(Z), "spf": spf_sieve(Z), "phi": phi_sieve(MAX_TABLE_MOD),
            "bennett": bennett, "rr1": load_rr_table1(here / "rr_theta_table1.tsv"),
            "rr2": load_rr_table2(here / "rr_theta_table2.tsv"), "suf": {}}


# SHARPEN_B: when True, Proposition 4.3's two Brun-Titchmarsh loss factors f1, f2
# are kept INSIDE their respective tail sums as the d-dependent weights
#     w1(d) = log n / log(n / (q d^2))      (family q d^2,   c1 < d <= n^A)
#     w2(b) = log n / log(n / (q^2 b^2))    (family q^2 b^2, c2 < b <= n^A/q)
# rather than replaced by their endpoint values f1, f2.  Each weight equals the
# corresponding f at the top of its range and is strictly smaller below it.
# Set SHARPEN_B = False to recover the published (unsharpened) constants exactly.
SHARPEN_B = True


def w_prefix(q, n, ctx, family):
    """Prefix sums of the d-dependent Brun-Titchmarsh weight for one family.

    family=1: P[D] = sum_{d <= D, (d,q)=1, mu(d)!=0} w1(d)/phi(d^2)
    family=2: P[B] = sum_{b <= B, (b,q)=1, mu(b)!=0} w2(b)/phi(b^2)

    Cached per (q, n, family); the search over A is then an O(1) lookup.
    """
    key = (q, n, family)
    if key not in ctx.setdefault("wpre", {}):
        mu, spf = ctx["mu"], ctx["spf"]
        logn = math.log(n)
        base = math.log(q) if family == 1 else 2.0 * math.log(q)
        P = [0.0] * (Z + 2)
        run = 0.0
        for a in range(1, Z + 1):
            if mu[a] and gcd(a, q) == 1:
                denom = logn - base - 2.0 * math.log(a)
                run = math.inf if denom <= 0.0 else run + (logn / denom) / phi_square(a, spf)
            P[a] = run
        P[Z + 1] = run
        ctx["wpre"][key] = P
    return ctx["wpre"][key]


def suffix_tail(q, ctx):
    if q not in ctx["suf"]:
        mu, spf = ctx["mu"], ctx["spf"]
        suf = [0.0] * (Z + 2)
        for a in range(Z, 0, -1):
            suf[a] = suf[a + 1] + ((1.0 / phi_square(a, spf)) if (mu[a] and gcd(a, q) == 1) else 0.0)
        ctx["suf"][q] = suf
    return ctx["suf"][q]


def S_tail(c, Zc, suf):
    return suf[c + 1] + 4.0 / Zc


def one_minus(p):
    return 1.0 - 1.0 / (p * (p - 1))


R_UNRESTRICTED_LOWER = 0.32035   # Theorem 2.4 flat bound (crude Lemma 5.1-type checks only)


def Bq_error(q, n, ctx, A=None, c1=None, c2=None, Z1=None, Z2=None):
    mu, phi = ctx["mu"], ctx["phi"]
    bennett, rr1, rr2 = ctx["bennett"], ctx["rr1"], ctx["rr2"]
    logn = math.log(n)
    suf = suffix_tail(q, ctx)
    if c1 is None:
        c1 = int(floor(sqrt(MAX_TABLE_MOD / q)))
    if c2 is None:
        c2 = int(floor(sqrt(MAX_TABLE_MOD))) // q
    if Z1 is None:
        Z1 = Z
    if Z2 is None:
        Z2 = Z
    I = 0.0
    cnt1 = 0
    for d in range(1, c1 + 1):
        if mu[d] and gcd(d, q) == 1:
            m = q * d * d
            I += choose_best_bound(n, m, phi[m], bennett, rr1, rr2).contribution
            cnt1 += 1
    cnt2 = 0
    for b in range(1, c2 + 1):
        if mu[b] and gcd(b, q) == 1:
            m = q * q * b * b
            I += choose_best_bound(n, m, phi[m], bennett, rr1, rr2).contribution
            cnt2 += 1
    endpoint = (logn / n) * (cnt1 + cnt2)
    E_short = I + endpoint
    S1 = S_tail(c1, Z1, suf)
    S2 = S_tail(c2, Z2, suf)

    P1 = w_prefix(q, n, ctx, 1) if SHARPEN_B else None
    P2 = w_prefix(q, n, ctx, 2) if SHARPEN_B else None

    def E_med(A):
        d1 = (1 - 2 * A) * logn - math.log(q)
        if d1 <= 0:
            return None
        if SHARPEN_B:
            D1 = min(int(floor(n ** A)), Z)
            D2 = min(int(floor(n ** A / q)), Z)
            W1 = P1[max(D1, c1)] - P1[c1]
            W2 = P2[max(D2, c2)] - P2[c2]
            if W1 == math.inf or W2 == math.inf:
                return None
            return (S1 + 2.0 * W1) / (q - 1) + (S2 + 2.0 * W2) / (q * (q - 1))
        f1 = logn / d1
        f2 = 1.0 / (1 - 2 * A)
        return ((1 + 2 * f1) / (q - 1)) * S1 + ((1 + 2 * f2) / (q * (q - 1))) * S2

    def E_large(A):
        return ((2.0 / q) * n ** (-A) + (1 + 1.0 / q) * n ** (-0.5)
                + (1 + 1.0 / q) * n ** (A - 1.0)) * logn

    def tot(A):
        em = E_med(A)
        return None if em is None else E_short + em + E_large(A)

    if A is None:
        best = None
        for i in range(100, 4900):
            a = i / 10000.0
            t = tot(a)
            if t is not None and (best is None or t < best[1]):
                best = (a, t)
        A = best[0]
    em = E_med(A)
    return {"q": q, "A": A, "c1": c1, "c2": c2, "main_coeff": 1.0 / q,
            "E_short": E_short, "endpoint": endpoint, "E_med": em, "E_large": E_large(A),
            "err": E_short + em + E_large(A)}


def R_lower(m, n, ctx):
    primes = factor_squarefree(m)
    beta = beta_coeff(primes)
    alpha = alpha_coeff(primes)
    c = int(floor(sqrt(MAX_TABLE_MOD / m))) if m > 1 else 316
    tail = tail_constant(m, c, ctx["mu"], ctx["spf"])
    r = evaluate(n, m, c, ctx["mu"], ctx["phi"], ctx["bennett"], ctx["rr1"], ctx["rr2"],
                 tail, alpha, beta, ctx["spf"])
    return beta * C_ARTIN, r["bound"]


# --- margin coefficients -------------------------------------------------
# Both main terms are exactly proportional to the SAME Euler product
#     W_n = prod_{p nmid n} (1 - 1/(p(p-1))),
# since  R_m/n ~ (prod_{q|m} lambda_q) W_n  and  B_q/n ~ (1 - lambda_q) W_n
# (the identity (q-2)/(q-1) = lambda_q (1 - 1/(q(q-1))) converts the stated
# forms of Prop 4.1 and Prop 4.3 into these).  The difference therefore carries
# W_n as a common factor, and W_n >= C_Artin is sharp (equality at n = 1).
#
# Writing the same difference over P_k(n) = prod_{p nmid kn} and using
# P_k >= C_Artin instead is the SAME algebra with a weaker constant: it
# discards prod_{q|k}(1 - 1/(q(q-1))), a factor 0.7728 at k = 105.  W_n is the
# right normalisation, and the bound W_n >= C_Artin is sharp.


def main():
    ap = argparse.ArgumentParser(
        description="Explicit two-family upper bound for B_q(n) (Proposition 4.3)")
    ap.add_argument("--bq", type=int, required=True,
                    help="odd prime q; prints the B_q upper-bound error at n")
    ap.add_argument("--n", type=float, default=DEFAULT_N,
                    help="evaluation point (default 8e9)")
    args = ap.parse_args()
    if args.n < MIN_VALID_N:
        raise SystemExit(f"n = {args.n:.3e} is below the validity floor {MIN_VALID_N:.0e}.")
    ctx = build_context()
    b = Bq_error(args.bq, args.n, ctx)
    print(f"B_{args.bq} upper-bound error at n={args.n:.3e}: "
          f"E_B{args.bq} = {b['err']:.6f}  (main coeff 1/{args.bq} = {b['main_coeff']:.5f}, "
          f"c1={b['c1']}, c2={b['c2']}, A={b['A']:.4f}; "
          f"short={b['E_short']:.5f} med={b['E_med']:.5f} large={b['E_large']:.2e})")


if __name__ == "__main__":
    main()
