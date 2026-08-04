#!/usr/bin/env python3
"""
Threshold.py
============
Computes the SHARP threshold n_0(r): the least N such that every n >= N (of the
appropriate parity) is a sum of a prime and a square-free number coprime to k,
for every k > 1 with omega(k) <= r.

WHY NOT THE GREEDY SEARCH.  The reference script of Lee and O'Clarey
(Lemma4_2.py) certifies n by exhibiting r+1 PAIRWISE COPRIME square-free
witnesses.  That is sufficient but not necessary: an n at which the greedy
search fails may still admit a representation for every admissible k, since a
representation need only avoid the primes of the single k in hand.  A threshold
printed in a theorem statement should be the exact one, so this script decides
the question directly.

THE EXACT CRITERION.  Put
    W(n) = { n - p : p < n prime, mu^2(n-p) = 1 },
and for w in W(n) let supp(w) be the set of ODD primes dividing w.  A modulus k
defeats n precisely when every w in W(n) shares a prime factor with k, i.e. when
the odd prime factors of k form a HITTING SET for the family { supp(w) }.  Hence

    n is good for every odd k with omega(k) <= r    <=>    tau(n) > r,

tau(n) being the least size of a hitting set for that family.  If some w has
supp(w) empty (that is, n-1 or n-2 is prime) then no k defeats n.

For EVEN k = 2m the square-free summand must also be odd, so the relevant family
is W_odd(n) = { w in W(n) : w odd } and, omega(k) <= r allowing omega(m) <= r-1,
the criterion is tau_odd(n) > r-1.  For odd n the even-k case is the genuine
parity obstruction and is excluded from the theorem, exactly as at omega(k)<=3.

EARLY EXIT.  Hitting-set size is monotone under adding constraints: if some
SUBFAMILY already needs more than r primes, so does the whole family.  Witnesses
are therefore generated in increasing order of w -- smallest first, hence fewest
prime factors, hence most constraining -- and the search stops as soon as the
partial family forces tau > r.  Only genuinely failing n are examined in full.

Usage:
    python3 Threshold.py --rmax 6 --nmax 20000
"""
import argparse


def sieve(nmax):
    """Smallest-prime-factor table and primality flags up to nmax."""
    spf = list(range(nmax + 1))
    i = 2
    while i * i <= nmax:
        if spf[i] == i:
            for j in range(i * i, nmax + 1, i):
                if spf[j] == j:
                    spf[j] = i
        i += 1
    isp = [False] * (nmax + 1)
    for m in range(2, nmax + 1):
        isp[m] = spf[m] == m
    return spf, isp


def odd_support(w, spf):
    """frozenset of odd primes dividing w, or None if w is not square-free."""
    s, x = [], w
    while x > 1:
        p = spf[x]
        x //= p
        if x % p == 0:
            return None
        if p > 2:
            s.append(p)
    return frozenset(s)


def min_hitting_set(fam, cap):
    """Least hitting-set size, or cap+1 as soon as it is known to exceed cap."""
    if not fam:
        return 0
    reduced = []
    for s in sorted(set(fam), key=len):
        if not any(t <= s for t in reduced):
            reduced.append(s)
    best = [cap + 1]

    def rec(rem, chosen):
        if not rem:
            best[0] = min(best[0], chosen)
            return
        if chosen + 1 > cap or chosen + 1 >= best[0]:
            return
        pivot = min(rem, key=len)
        for q in pivot:
            rec([s for s in rem if q not in s], chosen + 1)

    rec(reduced, 0)
    return best[0]


def family(n, spf, isp, odd_only):
    """Odd-prime supports of the square-free witnesses, smallest witness first.
    Returns None if some witness has empty support (n then defeats every k)."""
    fam = []
    for w in range(1, n - 1):
        if odd_only and w % 2 == 0:
            continue
        p = n - w
        if p < 2 or not isp[p]:
            continue
        s = odd_support(w, spf)
        if s is None:
            continue
        if not s:
            return None
        fam.append(s)
        if len(fam) >= 6 and len(fam) % 6 == 0 and min_hitting_set(fam, 8) > 8:
            break
    return fam


def tau(n, r, spf, isp, odd_only):
    fam = family(n, spf, isp, odd_only)
    return r + 1 if fam is None else min_hitting_set(fam, r)


def defeating_modulus(n, r, spf, isp, odd_only):
    fam = family(n, spf, isp, odd_only)
    if fam is None:
        return None
    reduced = []
    for s in sorted(set(fam), key=len):
        if not any(t <= s for t in reduced):
            reduced.append(s)
    best = [None]

    def rec(rem, chosen):
        if best[0] is not None and len(chosen) >= len(best[0]):
            return
        if not rem:
            best[0] = list(chosen)
            return
        if len(chosen) + 1 > r:
            return
        for q in min(rem, key=len):
            rec([s for s in rem if q not in s], chosen + [q])

    rec(reduced, [])
    return best[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rmin", type=int, default=2)
    ap.add_argument("--rmax", type=int, default=6)
    ap.add_argument("--nmax", type=int, default=20000)
    ap.add_argument("--certify", nargs=2, type=int, metavar=("LO", "HI"),
                    help="certify every n in [LO,HI] for omega(k) <= rmax")
    args = ap.parse_args()

    if args.certify:
        lo, hi, r = args.certify[0], args.certify[1], args.rmax
        spf, isp = sieve(hi + 2)
        bad = []
        for n in range(lo, hi + 1):
            if n % 2:
                ok = tau(n, r, spf, isp, False) > r
            else:
                ok = (tau(n, r, spf, isp, False) > r
                      and tau(n, r - 1, spf, isp, True) > r - 1)
            if not ok:
                bad.append(n)
        print(f"exact criterion, omega(k) <= {r}, every n in [{lo}, {hi}]:")
        if bad:
            print(f"  FAILURES ({len(bad)}): {bad}")
        else:
            print(f"  ALL {hi - lo + 1} VALUES CERTIFIED "
                  f"(odd n: every odd k; even n: also every even k = 2m)")
        raise SystemExit(1 if bad else 0)

    spf, isp = sieve(args.nmax + 2)
    print(f"exact sharp thresholds, every n <= {args.nmax} examined")
    print(f"{'r':>3s} {'largest failing n':>18s} {'n_0':>8s}   defeating modulus")
    for r in range(args.rmin, args.rmax + 1):
        worst, desc = None, ""
        for n in range(3, args.nmax + 1):
            if n % 2:
                hit = tau(n, r, spf, isp, False) <= r
                mode = (False, r)
            else:
                b1 = tau(n, r, spf, isp, False) <= r
                b2 = (r >= 1) and tau(n, r - 1, spf, isp, True) <= r - 1
                hit = b1 or b2
                mode = (False, r) if b1 else (True, r - 1)
            if hit:
                worst = n
                qs = defeating_modulus(n, mode[1], spf, isp, mode[0])
                if qs:
                    kk = 2 if mode[0] else 1
                    for q in qs:
                        kk *= q
                    desc = (f"n = {n}, k = {kk} = " + ("2*" if mode[0] else "")
                            + "*".join(str(q) for q in sorted(qs)))
        print(f"{r:>3d} {worst:>18d} {worst + 1:>8d}   {desc}")


if __name__ == "__main__":
    main()
