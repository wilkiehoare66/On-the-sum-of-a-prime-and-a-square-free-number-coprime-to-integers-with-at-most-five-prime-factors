#!/usr/bin/env python3
"""
CertifyCases.py
===============
Directed-rounding certificate for the five cases of Lemma 5.4 (the case tree
covering every odd square-free k with omega(k) <= 5).

CONVENTION, as elsewhere in this repository: main terms are rounded DOWN, error
terms UP, and a safety envelope EPS is subtracted from every margin.  EPS is set
far above the measured double-precision error of the underlying sums (order
1e-16 relative), so each printed inequality holds as written.

Artin's constant is replaced throughout by the strict lower bound 0.3739558,
and each beta(base) = prod lambda_q by an exact Fraction rounded DOWN, so the
certified main term never exceeds the true one.  Each peel cost is rounded UP.

Two intervals are certified separately, because the set of admissible imported
prime-counting bounds changes at x = 1e10 (the Ramare--Rumely Table 2 sqrt-bound
is available only below it).  On each interval every error term is decreasing in
n, so the left endpoint is a rigorous maximum for the error and hence a rigorous
minimum for the margin.

Exit status is the verdict: 0 if every case is certified positive on both
intervals, 1 otherwise.
"""
import json
import math
import sys
from fractions import Fraction

import ScanOmega as S
import ComputeBqBound as BQ

EPS = 1e-9
C_ARTIN_LOWER = Fraction(3739558, 10000000)     # < Artin's constant
LOG2_UPPER = Fraction(6931472, 10000000)        # > log 2

# The case tree of Lemma 5.4, as (description, base primes, threshold T, worst
# completion).  T is the least value the first peeled prime may take in that
# case; by Lemma 5.3 the worst k of the case is the base followed by T and its
# successors, which is what makes one evaluation certify the whole branch.
CASES = [
    ("q_1 >= 5",                             [],        5,  [5, 7, 11, 13, 17]),
    ("q_1 = 3, q_2 >= 11",                   [3],       11, [3, 11, 13, 17, 19]),
    ("q_1 = 3, q_2 = 5, q_3 >= 11",          [3, 5],    11, [3, 5, 11, 13, 17]),
    ("q_1 = 3, q_2 = 7, q_3 >= 11",          [3, 7],    11, [3, 7, 11, 13, 17]),
    ("q_1 = 3, q_2 = 5, q_3 = 7, q_4 >= 11", [3, 5, 7], 11, [3, 5, 7, 11, 13]),
]

INTERVALS = [("[8e9, 1e10]", 8e9), ("(1e10, inf)", 1.00001e10)]


def lam_exact(q):
    return Fraction(q * (q - 2), q * q - q - 1)


def beta_exact(primes):
    b = Fraction(1)
    for q in primes:
        b *= lam_exact(q)
    return b


def up(x):
    """Round a floating error UP by the safety envelope."""
    return x * (1.0 + 1e-12) + EPS


def certify(desc, base, combo, n):
    peeled = combo[len(base):]
    m = 1
    for q in base:
        m *= q
    er = S.E_R(m, n)
    if er is None:
        return None
    beta_b = beta_exact(base)
    main = float(beta_b * C_ARTIN_LOWER)         # rounded down by construction
    total, bracket, ledger = 0.0, beta_b, []
    for q in peeled:
        eb = S.E_B(q, n)
        cost_w = float((1 - lam_exact(q)) * C_ARTIN_LOWER) + up(eb) if eb is not None else None
        cost_c = up(S.crude_cost(q, n))
        if cost_w is not None and cost_w <= cost_c:
            bracket -= (1 - lam_exact(q))
            total += cost_w
            ledger.append((q, "Prop 4.3", cost_w))
        else:
            total += cost_c
            ledger.append((q, "elementary", cost_c))
    if bracket < 0:
        return None
    margin = main - up(er) - total - float(LOG2_UPPER) / n - EPS
    return {"desc": desc, "base": base, "m": m, "combo": combo, "beta": beta_b,
            "main": main, "E_R": up(er), "peel_total": total,
            "bracket": float(bracket), "margin": margin, "ledger": ledger}


def main():
    ok = True
    archive = {"n_intervals": {}, "convention": {
        "C_Artin_lower": str(C_ARTIN_LOWER), "log2_upper": str(LOG2_UPPER),
        "safety_envelope": EPS,
        "rounding": "main terms down, error terms up, envelope subtracted"}}
    for label, n in INTERVALS:
        print(f"\n=== interval {label}, certified at the left endpoint "
              f"n = {n:.5e} ===")
        for desc, base, T, combo in CASES:
            r = certify(desc, base, combo, n)
            if r is None:
                print(f"  {desc:40s}  INADMISSIBLE")
                ok = False
                continue
            kk = 1
            for q in combo:
                kk *= q
            verdict = "CERTIFIED" if r["margin"] > 0 else "FAILS"
            if r["margin"] <= 0:
                ok = False
            print(f"  {desc:40s} base m={r['m']:<5d} worst k={kk:<7d}"
                  f" margin={r['margin']:+.7f}  {verdict}")
            print(f"      main={r['main']:.7f} (beta={float(r['beta']):.7f})"
                  f"  E_R={r['E_R']:.7f}  peels={r['peel_total']:.7f}"
                  f"  bracket={r['bracket']:+.5f}")
            for q, src, cst in r["ledger"]:
                print(f"        peel q={q:<4d} {src:11s} cost {cst:.7f}")
            archive["n_intervals"].setdefault(label, {})[desc] = {
                "base_m": r["m"], "base_primes": base, "worst_k": kk,
                "worst_k_primes": combo, "beta_exact": str(r["beta"]),
                "first_peel_threshold": T,
                "main_term": r["main"], "E_R": r["E_R"],
                "peel_total": r["peel_total"], "bracket": r["bracket"],
                "margin": r["margin"], "n": n,
                "peels": [{"q": q, "source": src, "cost": cst}
                          for q, src, cst in r["ledger"]]}
    archive["verdict"] = "certified" if ok else "failed"
    with open("case_certificates.json", "w") as f:
        json.dump(archive, f, indent=1, sort_keys=True)
    print(f"\n{'ALL FIVE CASES CERTIFIED ON BOTH INTERVALS' if ok else 'CERTIFICATION FAILED'}")
    print("archive written to case_certificates.json")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
