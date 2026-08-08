#!/usr/bin/env python3
"""
ScanOmega.py
============
Evaluates the base-plus-peel criterion for an arbitrary odd square-free modulus
k = q_1 ... q_r (q_1 < ... < q_r), at a given n, using the sharpened
Theorems 4.1 and 4.4 (ComputeRkBound.py / ComputeBqBound.py with
SHARPEN / SHARPEN_B on).

THE CRITERION.  Choose a base m = q_1...q_s (an initial segment) and peel the
remaining primes.  Both main terms are proportional to the same Euler product
W_n = prod_{p nmid n}(1 - 1/(p(p-1))):

    R_m(n)/n   >=  beta(m) W_n - E_R(m,n),        beta(m) = prod_{q|m} lambda_q,
    B_q(n)/n   <=  (1 - lambda_q) W_n + E_B(q,n).

Peeled primes may be costed in either of two currencies:

  (W-currency)  via Theorem 4.4, cost  (1 - lambda_q) C_Artin + E_B(q,n).
                Requires the bracket  beta(m) - sum_{q in W} (1 - lambda_q) >= 0,
                since only then does W_n >= C_Artin apply in the right direction.
  (crude)       via the elementary bound (10) / Lemma 5.4, cost
                    1/q + (q-1)/(840 log n)                for q <= 100,
                    1/(q-1) + 1/(160 log n)                for 100 < q <= 1e5,
                    0.00026                                for q > 1e5.
                Absolute, so it carries no bracket constraint.

Both costs are absolute (in units of n), so each peeled prime independently takes
the cheaper of the two, and

    margin = beta(m) C_Artin - E_R(m,n) - sum_q cost(q) - log2/n.

s = r is the direct bound (no peel); s = 0 is the comparison of Section 3, with
E_R(1,n) the genuine explicit error at k = 1 rather than a flat lower bound.

Usage:
    python3 ScanOmega.py --k 15015                 # one modulus
    python3 ScanOmega.py --omega 5 --qmax 101      # exhaustive scan
    python3 ScanOmega.py --omega 5 --worst         # worst-case search (see below)
"""
import argparse
import math
from math import gcd, floor, sqrt
from itertools import combinations

import ComputeBqBound as BQ
import ComputeRkBound as RK

C_ARTIN = RK.C_ARTIN
LOG2 = math.log(2.0)
MAXMOD = RK.MAX_TABLE_MOD          # 1e5: the Bennett modulus cap
BASE_CAP = MAXMOD // 4             # largest base admitting c >= 2

_ctx = None
_cacheB, _cacheR = {}, {}


def ctx():
    global _ctx
    if _ctx is None:
        _ctx = BQ.build_context()
    return _ctx


def lam(q):
    return q * (q - 2) / (q * q - q - 1.0)


def beta(primes):
    b = 1.0
    for q in primes:
        b *= lam(q)
    return b


def E_B(q, n):
    """Theorem 4.4 error for B_q, or None if the modulus cap forbids it."""
    if q > MAXMOD:
        return None
    key = (q, n)
    if key not in _cacheB:
        try:
            _cacheB[key] = BQ.Bq_error(q, n, ctx())["err"]
        except Exception:
            _cacheB[key] = None
    return _cacheB[key]


def E_R(m, n):
    """Theorem 4.1 error for R_m, or None if the modulus cap forbids it.

    Theorem 4.1 needs an admissible c_theta(e d^2) for every e | m and
    d <= c, so c = floor(sqrt(MAXMOD/m)); c >= 1 requires m <= MAXMOD, and the
    tail term decreases like 1/c, so in practice m <= MAXMOD/4 (c >= 2) is the
    working range.
    """
    if m > MAXMOD:
        return None
    key = (m, n)
    if key not in _cacheR:
        try:
            _, bound = BQ.R_lower(m, n, ctx())
            _cacheR[key] = beta(BQ.factor_squarefree(m) if m > 1 else []) * C_ARTIN - bound
        except Exception:
            _cacheR[key] = None
    return _cacheR[key]


def crude_cost(q, n):
    """The elementary upper bound for B_q(n)/n: equation (10) and Lemma 5.4."""
    logn = math.log(n)
    if q <= 100:
        return 1.0 / q + (q - 1.0) / (840.0 * logn)
    if q <= MAXMOD:
        return 1.0 / (q - 1.0) + 1.0 / (160.0 * logn)
    return 0.00026


Q_ENV = 1000          # cost(q) is computed exactly for every prime q <= Q_ENV
_ENV = {}


def cost_envelope(n):
    """A certified NON-INCREASING majorant for the peel cost.

    cost(q) itself is NOT monotone: at n = 8e9 it rises from 0.0095116 at q = 97
    to 0.0095721 at q = 101, because c_1 = floor(sqrt(M/q)) drops from 32 to 31
    there and the tabulated c_theta(qd^2) entering the short-range term of
    Theorem 4.4 are not monotone in the modulus.  Any argument that reduces
    infinitely many completions to one worst case therefore cannot appeal to
    monotonicity of cost.

    Instead we build the running supremum from the right,

        cost_bar(q) = sup { cost(q') : q' >= q prime },

    which is non-increasing by construction and needs no monotonicity
    assumption whatever.  It is computable because cost(q) is evaluated exactly
    for the finitely many primes q <= Q_ENV, while for q > Q_ENV we use
    cost(q) <= crude(q) together with crude(q) <= crude(Q_ENV^+) for every
    prime q > Q_ENV: the second case of (10) is decreasing in q, and the third
    is the constant 0.00026 < 1/1008 <= crude(1009).  So the supremum over that
    tail is bounded by crude at the first prime above Q_ENV.

    Note that crude is NOT globally non-increasing, and the argument does not
    claim it is: at the junction q = M the second case of (10) decays like
    1/(160 log n) while the third is constant, so the third is the larger once
    n > e^25 ~ 7.2e10.  Only the domination above Q_ENV is used, and that holds
    for every n >= 8e9 with a factor of nearly four to spare.  This is the tail
    step in the proof of Lemma 5.3.

    Returns (env, tail): env[q] = cost_bar(q) for primes q <= Q_ENV, and tail is
    the uniform bound for q > Q_ENV.  Where cost happens to be monotone the
    envelope coincides with it, so nothing is lost.
    """
    key = n
    if key in _ENV:
        return _ENV[key]
    ps = [p for p in range(3, Q_ENV + 1)
          if all(p % d for d in range(2, int(sqrt(p)) + 1))]
    t = Q_ENV + 1
    while not all(t % d for d in range(2, int(sqrt(t)) + 1)):
        t += 1
    tail = crude_cost(t, n)
    env, run = {}, tail
    for p in reversed(ps):
        eb = E_B(p, n)
        cw = (1.0 - lam(p)) * C_ARTIN + eb if eb is not None else float("inf")
        run = max(run, min(cw, crude_cost(p, n)))
        env[p] = run
    _ENV[key] = (env, tail)
    return env, tail


def cost_bar(q, n):
    env, tail = cost_envelope(n)
    return env[q] if q <= Q_ENV else min(tail, crude_cost(q, n))


def margin(primes, n, s):
    """Margin for base = primes[:s], peeling primes[s:].  None if inadmissible."""
    base, peeled = list(primes[:s]), list(primes[s:])
    m = 1
    for q in base:
        m *= q
    er = E_R(m, n)
    if er is None:
        return None
    # Each peeled prime is charged the ENVELOPE cost_bar, not cost itself, so
    # that the evaluation at the worst completion bounds every larger
    # completion without assuming cost is monotone.  The bracket is taken in
    # its strictest form, debiting (1 - lambda_q) for EVERY peeled prime: the
    # constraint is only needed for those actually bounded by Theorem 4.4,
    # and dropping a prime from that set only increases the bracket, so this is
    # valid whichever currency achieves the minimum at each prime.
    total, bracket = 0.0, beta(base)
    for q in peeled:
        total += cost_bar(q, n)
        bracket -= (1.0 - lam(q))
    if bracket < 0.0:
        return None
    return beta(base) * C_ARTIN - er - total - LOG2 / n


def best_device(primes, n):
    """Best (device-name, margin) over all initial-segment bases."""
    r = len(primes)
    out = []
    for s in range(r, -1, -1):
        v = margin(primes, n, s)
        if v is not None:
            out.append(("direct" if s == r else ("cmp" if s == 0 else f"peel{r - s}"), v))
    return max(out, key=lambda t: t[1]) if out else ("none", float("-inf"))


def main():
    ap = argparse.ArgumentParser(description="Base-plus-peel coverage scan")
    ap.add_argument("--k", type=int)
    ap.add_argument("--omega", type=int)
    ap.add_argument("--qmax", type=int, default=101)
    ap.add_argument("--n", type=float, default=8e9)
    ap.add_argument("--show", type=int, default=8, help="how many worst cases to print")
    args = ap.parse_args()

    if args.k:
        ps = BQ.factor_squarefree(args.k)
        dev, s = best_device(ps, args.n)
        print(f"k={args.k} {ps}  device={dev}  margin={s:+.6f}  {'CLEARS' if s > 0 else 'FAILS'}")
        for j in range(len(ps), -1, -1):
            v = margin(ps, args.n, j)
            tag = "direct" if j == len(ps) else ("cmp" if j == 0 else f"peel{len(ps)-j}")
            print(f"    {tag:8s} base={ps[:j]!s:24s} " +
                  ("inadmissible" if v is None else f"margin={v:+.6f}"))
        return

    if args.omega:
        pool = [p for p in range(3, args.qmax + 1) if all(p % d for d in range(2, int(sqrt(p)) + 1))]
        fails, tested, results = [], 0, []
        for combo in combinations(pool, args.omega):
            tested += 1
            dev, s = best_device(list(combo), args.n)
            results.append((s, combo, dev))
            if s <= 0:
                fails.append((s, combo, dev))
        results.sort()
        print(f"omega={args.omega}, primes<={args.qmax}: tested {tested}, "
              f"FAILURES {len(fails)}")
        for s, combo, dev in results[:args.show]:
            k = 1
            for q in combo:
                k *= q
            print(f"   {'FAIL' if s <= 0 else 'ok  '} k={k:<12d} {combo}  {dev:8s} {s:+.6f}")


if __name__ == "__main__":
    main()
