#!/usr/bin/env python3
"""
verify_R_classes.py
===================
Brute-force check of the residue-class structure underlying the (R) step of
Proposition 4.1 (the retention of the coprimality condition in Sigma_2).

CLAIM.  Let k be odd and square-free, (n,k)=1, and let d satisfy (d,kn)=1.
Consider the residues r mod k d^2 with

        r == n (mod d^2)   and   r != n (mod q) for every q | k,

i.e. exactly the classes in which a prime p counted by
Sigma_2's inner sum can lie.  Then

  (a) there are exactly phi(k) = prod_{q|k} (q-1) such classes;
  (b) exactly prod_{q|k} (q-2) of them are coprime to k d^2;
  (c) every one of the remaining phi(k) - prod (q-2) classes consists entirely
      of integers divisible by some q | k, hence contains at most one prime,
      and that prime is q itself.

(a) and (b) give the factor alpha = prod (q-2)/(q-1) after Brun-Titchmarsh on the
modulus k d^2; (c) bounds the discarded classes by sum_{q|k} log q = log k in
total, since each prime q | k lies in exactly one class mod k d^2.

Run: python3 verify_R_classes.py
Exits 0 on success, 1 on any mismatch.
"""
from math import gcd
from itertools import combinations

ODD_PRIMES = [3, 5, 7, 11, 13, 17, 19, 23]


def squarefree_odd(x):
    if x % 2 == 0:
        return False
    m, p = x, 3
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return False
        p += 2
    return True


def primes_of(k):
    out, m, p = [], k, 3
    while p * p <= m:
        if m % p == 0:
            out.append(p)
            m //= p
        p += 2
    if m > 1:
        out.append(m)
    return out


def check(k, n, d):
    qs = primes_of(k)
    M = k * d * d
    classes = [r for r in range(M)
               if r % (d * d) == n % (d * d) and all((r - n) % q != 0 for q in qs)]
    expect_total = 1
    for q in qs:
        expect_total *= (q - 1)
    coprime = [r for r in classes if gcd(r, M) == 1]
    expect_coprime = 1
    for q in qs:
        expect_coprime *= (q - 2)
    if len(classes) != expect_total or len(coprime) != expect_coprime:
        return False, f"count: {len(classes)}/{expect_total}, {len(coprime)}/{expect_coprime}"
    # (c): every non-coprime class is entirely divisible by some q | k
    for r in classes:
        if gcd(r, M) == 1:
            continue
        g = gcd(r, M)
        if not any(g % q == 0 for q in qs):
            return False, f"class {r} not divisible by any q|k"
        # every member of the class shares that divisor, so the only prime is q
        q = next(q for q in qs if r % q == 0)
        if any((r + t * M) % q != 0 for t in range(4)):
            return False, f"class {r} not uniformly divisible by {q}"
    return True, None


def main():
    cases = 0
    for rlen in (1, 2, 3):
        for qs in combinations(ODD_PRIMES[:6], rlen):
            k = 1
            for q in qs:
                k *= q
            if k > 400:
                continue
            for d in (1, 2, 4, 5, 8, 11):
                if gcd(d, k) != 1:
                    continue
                for n in range(1, 60):
                    if gcd(n, k * d) != 1:
                        continue
                    ok, msg = check(k, n, d)
                    cases += 1
                    if not ok:
                        print(f"MISMATCH k={k} n={n} d={d}: {msg}")
                        return 1
    print(f"verify_R_classes: {cases} cases, 0 mismatches -- (a),(b),(c) all hold.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
