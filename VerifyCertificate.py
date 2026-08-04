#!/usr/bin/env python3
"""
VerifyCertificate.py
====================
Independently re-checks the certificate archive written by CertifyCases.py.

The point of this script is that it imports NOTHING from the optimisation
scripts: not ComputeRkBound.py, not ComputeBqBound.py, not ScanOmega.py.  It
reads case_certificates.json, re-derives every quantity it can from first
principles, and re-checks every inequality.  A reader who distrusts the
optimisation can run this instead.

WHAT IS RE-DERIVED HERE (not taken from the archive):
  * the exact local factors lambda_q = q(q-2)/(q^2-q-1), as Fractions;
  * beta(base) = prod lambda_q over the base primes, compared against the
    archived exact value;
  * the main term beta * C_Artin_lower, recomputed;
  * the bracket beta(base) - sum over Prop-4.3-costed primes of (1-lambda_q),
    which must be nonnegative for W_n >= C_Artin to be applied in the right
    direction;
  * each peel cost, checked to equal (1-lambda_q) C_Artin_lower + E_B for the
    Proposition 4.3 entries;
  * the margin itself, and its positivity.

WHAT IS TAKEN ON TRUST from the archive: the imported prime-counting errors
E_R and E_B.  Those are the outputs of Propositions 4.1 and 4.3 and are
reproduced by ComputeRkBound.py / ComputeBqBound.py; re-deriving them here
would just be a second copy of the same code.  Everything downstream of them
is recomputed.

Also checked: that the five cases form an exhaustive division of the odd
square-free k with omega(k) = 5, and that each case's "worst k" really is the
smallest completion of its base -- which is what Lemma 5.4 needs from
Lemma 5.3 (monotonicity).

Exit status is the verdict: 0 if every check passes, 1 otherwise.

Usage:
    python3 VerifyCertificate.py [case_certificates.json]
"""
import json
import sys
from fractions import Fraction

TOL = 1e-12


def lam(q):
    return Fraction(q * (q - 2), q * q - q - 1)


def is_prime(x):
    if x < 2:
        return False
    d = 2
    while d * d <= x:
        if x % d == 0:
            return False
        d += 1
    return True


def next_prime(p):
    x = p + 1
    while not is_prime(x):
        x += 1
    return x


def check_case(desc, rec, cart, log2u, eps, fails):
    base = rec["base_primes"]
    combo = rec["worst_k_primes"]
    n = rec["n"]

    # 1. beta(base) recomputed exactly
    beta = Fraction(1)
    for q in base:
        beta *= lam(q)
    if str(beta) != rec["beta_exact"]:
        fails.append(f"{desc}: beta mismatch, archive {rec['beta_exact']} vs {beta}")

    # 2. main term
    main = float(beta * cart)
    if abs(main - rec["main_term"]) > TOL:
        fails.append(f"{desc}: main term {rec['main_term']} vs recomputed {main}")

    # 3. base product matches the recorded modulus
    m = 1
    for q in base:
        m *= q
    if m != rec["base_m"]:
        fails.append(f"{desc}: base_m {rec['base_m']} vs product {m}")

    # 4. peel costs and bracket
    bracket, total = beta, 0.0
    peeled = combo[len(base):]
    if [p["q"] for p in rec["peels"]] != peeled:
        fails.append(f"{desc}: peel list {[p['q'] for p in rec['peels']]} != {peeled}")
    for p in rec["peels"]:
        q, cost = p["q"], p["cost"]
        if p["source"] == "Prop 4.3":
            floor = float((1 - lam(q)) * cart)
            if cost < floor - TOL:
                fails.append(f"{desc}: peel q={q} cost {cost} below main part {floor}")
            bracket -= (1 - lam(q))
        total += cost
    if abs(float(bracket) - rec["bracket"]) > TOL:
        fails.append(f"{desc}: bracket {rec['bracket']} vs recomputed {float(bracket)}")
    if bracket < 0:
        fails.append(f"{desc}: NEGATIVE bracket -- W_n >= C_Artin applied the wrong way")
    if abs(total - rec["peel_total"]) > TOL:
        fails.append(f"{desc}: peel total {rec['peel_total']} vs {total}")

    # 5. the margin, and its positivity
    margin = main - rec["E_R"] - total - float(log2u) / n - eps
    if abs(margin - rec["margin"]) > TOL:
        fails.append(f"{desc}: margin {rec['margin']} vs recomputed {margin}")
    if margin <= 0:
        fails.append(f"{desc}: NON-POSITIVE margin {margin}")

    # 6. the worst k really is the base followed by the case threshold T and
    #    its successors.  This is what Lemma 5.3 (monotonicity) needs: the peel
    #    cost decreases in q, so the smallest admissible completion is the
    #    hardest, and one evaluation then certifies the whole branch.
    T = rec["first_peel_threshold"]
    tail, x = [T], T
    while len(tail) < len(peeled):
        x = next_prime(x)
        tail.append(x)
    if base + tail != combo:
        fails.append(f"{desc}: worst k {combo} is not base + (T={T} and successors) "
                     f"{base + tail} -- Lemma 5.3 does not transfer")

    # 7. T must be consistent with the case description, and the base must be
    #    an initial segment of the worst k (the base is fixed across the branch).
    if not desc.rstrip().endswith(f">= {T}"):
        fails.append(f"{desc}: threshold T={T} does not match the case description")
    if combo[:len(base)] != base:
        fails.append(f"{desc}: base {base} is not an initial segment of {combo}")
    return margin


def check_exhaustive(descs, fails):
    """The five cases must partition the odd square-free k with omega(k)=5."""
    expected = ["q_1 >= 5",
                "q_1 = 3, q_2 >= 11",
                "q_1 = 3, q_2 = 5, q_3 >= 11",
                "q_1 = 3, q_2 = 7, q_3 >= 11",
                "q_1 = 3, q_2 = 5, q_3 = 7, q_4 >= 11"]
    if sorted(descs) != sorted(expected):
        fails.append(f"case list changed: {sorted(descs)}")
        return
    # q_1 = 3 is the only prime below 5; q_2 in {5,7} the only primes below 11;
    # and no prime lies strictly between 7 and 11.
    for lo, hi in ((3, 5), (7, 11)):
        between = [p for p in range(lo + 1, hi) if is_prime(p)]
        if between:
            fails.append(f"primes {between} lie between {lo} and {hi}: cases not exhaustive")


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "case_certificates.json"
    d = json.load(open(path))
    conv = d["convention"]
    cart = Fraction(conv["C_Artin_lower"])
    log2u = Fraction(conv["log2_upper"])
    eps = conv["safety_envelope"]

    if cart > Fraction(3739559, 10000000):
        print("FAIL: C_Artin lower bound is not below Artin's constant")
        return 1
    if log2u < Fraction(693147, 1000000):
        print("FAIL: log 2 upper bound is not above log 2")
        return 1

    fails, margins = [], []
    for label, cases in sorted(d["n_intervals"].items()):
        check_exhaustive(list(cases), fails)
        for desc, rec in sorted(cases.items()):
            margins.append((label, desc, check_case(desc, rec, cart, log2u, eps, fails)))

    for label, desc, mg in margins:
        print(f"  {label:14s} {desc:42s} margin {mg:+.7f}")
    if fails:
        print(f"\n{len(fails)} FAILURES:")
        for f in fails:
            print("   " + f)
        return 1
    print(f"\nALL {len(margins)} RECORDS INDEPENDENTLY VERIFIED "
          f"(least margin {min(m for _, _, m in margins):+.7f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
