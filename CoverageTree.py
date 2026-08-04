#!/usr/bin/env python3
"""
CoverageTree.py
===============
Produces the finite case division that establishes R_k(n) > 0 for EVERY odd
square-free k with omega(k) = r and (n,k) = 1, at a given n -- not merely for k
built from small primes.

THE MONOTONICITY THAT MAKES THIS FINITE
---------------------------------------
Fix a base m = q_1...q_s (an initial segment of k's primes) and peel the rest.
Each peeled prime q costs

    min( (1 - lambda_q) C_Artin + E_B(q,n),   crude(q,n) ),

and both expressions are strictly decreasing in q.  Hence, at a fixed base, the
margin is smallest when the peeled primes are as small as they are allowed to
be, i.e. the primes immediately following q_s.  So a single evaluation at that
worst completion certifies EVERY k sharing the prefix q_1...q_s.

THE CASE TREE
-------------
Starting from the empty prefix (base = 1, the Theorem 2.4 comparison), at each
node with prefix P = (q_1,...,q_s):

  * evaluate the worst completion of P; if the margin is positive, the whole
    branch is settled by a base that is an initial segment of P, and the node
    is a LEAF;
  * otherwise find the least prime T > q_s such that the prefix P + (T) is
    settled for every continuation -- by the monotonicity above, so is
    P + (q) for every prime q >= T -- and recurse on the finitely many
    primes q_s < q < T.

The tree is finite exactly when every branch terminates, i.e. when no prefix
survives to length r with product exceeding the modulus cap of Theorem 2.2.
A branch that reaches length r with product > 1e5 has no admissible device and
is reported as UNCOVERED.

Usage:
    python3 CoverageTree.py --omega 5
    python3 CoverageTree.py --omega 6 --n 1e10
"""
import argparse
import math
from math import sqrt

import ScanOmega as S

MAXMOD = S.MAXMOD


def is_prime(x):
    if x < 2:
        return False
    if x % 2 == 0:
        return x == 2
    d = 3
    while d * d <= x:
        if x % d == 0:
            return False
        d += 2
    return True


def next_prime(p):
    x = p + 2 if p > 2 else 3
    while not is_prime(x):
        x += 2
    return x


def successors(p, count):
    out, x = [], p
    for _ in range(count):
        x = next_prime(x)
        out.append(x)
    return out


def worst_margin(prefix, r, n, first_peel=None):
    """Margin for a base that is an INITIAL SEGMENT OF `prefix` -- hence fixed
    throughout the branch -- against the worst admissible completion.

    `first_peel` is the smallest value the next prime factor is allowed to take;
    the peeled primes are then first_peel and its successors.  Because every
    peel cost is decreasing in the prime, the margin is increasing in
    first_peel, so a positive value here certifies the whole branch
    q_{s+1} >= first_peel.

    Returns (margin, s_used, worst_k).
    """
    last = prefix[-1] if prefix else 1
    if first_peel is None:
        first_peel = next_prime(last)
    tail = [first_peel] + successors(first_peel, r - len(prefix) - 1) if r > len(prefix) else []
    combo = list(prefix) + tail
    best = (float("-inf"), None)
    for s in range(len(prefix), -1, -1):
        v = S.margin(combo, n, s)
        if v is not None and v > best[0]:
            best = (v, s)
    return best[0], best[1], combo


def settle(prefix, r, n, out, qcap=2000):
    """Settle the branch whose first len(prefix) primes are exactly `prefix`.

    At a leaf (len(prefix) == r) the base is all of k and the device is the
    direct bound of Proposition 4.1.  Otherwise we find the least prime T such
    that base = prefix clears every k with q_{s+1} >= T, and recurse on the
    finitely many primes strictly between max(prefix) and T.
    """
    if len(prefix) == r:
        v, s, combo = worst_margin(prefix, r, n)
        dev = "none" if s is None else ("direct" if s == r else
              ("cmp" if s == 0 else f"peel{r - s}"))
        out.append((tuple(prefix), s, v, tuple(combo), dev))
        return v > 0
    last = prefix[-1] if prefix else 1
    T, q, probes = None, next_prime(last), []
    while q <= qcap:
        vv, ss, cc = worst_margin(prefix, r, n, first_peel=q)
        probes.append((q, vv, ss, cc))
        if vv > 0:
            T = q
            break
        q = next_prime(q)
    ok = True
    for q, vv, ss, cc in (probes[:-1] if T is not None else probes):
        if not settle(prefix + [q], r, n, out, qcap):
            ok = False
    if T is None:
        out.append((tuple(prefix) + ("(no T)",), None, float("-inf"), (), "none"))
        return False
    q, vv, ss, cc = probes[-1]
    out.append((tuple(prefix) + (f">={T}",), ss, vv, tuple(cc),
                "cmp" if ss == 0 else f"peel{r - ss}"))
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--omega", type=int, required=True)
    ap.add_argument("--n", type=float, default=8e9)
    args = ap.parse_args()
    r, n = args.omega, args.n

    out = []
    ok = settle([], r, n, out)
    bad = [x for x in out if x[2] <= 0]
    print(f"omega(k) = {r}, n = {n:.3e}:  "
          f"{'COVERED' if ok else 'NOT COVERED'}   "
          f"({len(out)} cases, {len(bad)} failing)")
    print(f"{'case (first prime factors of k)':34s} {'device':8s} {'margin':>11s}   worst k")
    for pre, s, v, combo, dev in sorted(out, key=lambda t: (len(t[0]), [str(z) for z in t[0]])):
        kk = 1
        for q in combo:
            kk *= q
        flag = "  <-- FAILS" if v <= 0 else ""
        print(f"  {str(pre):32s} {dev:8s} {v:+11.6f}   "
              f"{str(list(combo)):26s} {'k='+str(kk) if combo else ''}{flag}")


if __name__ == "__main__":
    main()
