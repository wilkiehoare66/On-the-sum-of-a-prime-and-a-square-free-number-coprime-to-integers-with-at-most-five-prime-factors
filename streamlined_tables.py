#!/usr/bin/env python3
"""Tables 1 and 2 of streamlined.tex from Proposition 2.5 (E_1) and
Proposition 2.7 (E_2) exactly as printed, at n0 = 2e12, from the Bennett et
al. c_theta table.

Table 2 is computed with Z = 1e5 for the base and Z_1 = Z_2 = floor(n0^{A_q}/q)
for each peeled prime q -- the largest value the hypothesis q Z_i <= n0^A of
Proposition 2.7 allows.  Each A (base) and A_q (peeled prime) is optimised on
its own, since every error term depends on its own parameter only.  Lower
bounds are rounded down and upper bounds up, at the last printed digit.

Usage: python3 -I streamlined_tables.py PATH/TO/bennett_c_theta.tsv
"""
import math
import sys
from math import log, sqrt, floor
from fractions import Fraction

N0 = 2e12
L = log(N0)
CARTIN = 0.3739558
MMAX = 10**5
ZBIG = 10**5

ctheta = {}
with open(sys.argv[1]) as fh:
    next(fh)
    for line in fh:
        a, b = line.split()
        ctheta[int(a)] = float(b)

spf = list(range(ZBIG + 1))
for i in range(2, int(ZBIG ** 0.5) + 1):
    if spf[i] == i:
        for j in range(i * i, ZBIG + 1, i):
            if spf[j] == j:
                spf[j] = i
SQF, PHI2, PR = [False] * (ZBIG + 1), [0] * (ZBIG + 1), [()] * (ZBIG + 1)
for d in range(1, ZBIG + 1):
    x, ps, sf = d, [], True
    while x > 1:
        p = spf[x]
        x //= p
        if x % p == 0:
            sf = False
        if not ps or ps[-1] != p:
            ps.append(p)
    ph = d
    for p in ps:
        ph = ph // p * (p - 1)
    SQF[d], PHI2[d], PR[d] = sf, d * ph, tuple(ps)


def coprime(d, kp):
    return all(p not in kp for p in PR[d])


def prefix(kp, kval):
    """Prefix sums over square-free d coprime to kp of 1/phi(d^2) and of the
    Brun--Titchmarsh-weighted 1/phi(d^2) * (1 - log(kval d^2)/log n0)^{-1}."""
    P, W = [0.0] * (ZBIG + 1), [0.0] * (ZBIG + 1)
    s = w = 0.0
    for d in range(1, ZBIG + 1):
        if SQF[d] and coprime(d, kp):
            t = 1.0 / PHI2[d]
            s += t
            if kval * d * d < N0:
                w += t / (1 - log(kval * d * d) / L)
            else:
                w = float("nan")   # weight undefined: such Z is inadmissible
        P[d], W[d] = s, w
    return P, W


def E(m):
    return 0.375 / L ** 2 if m == 1 else ctheta[m]


def sqfree_divisors(kp):
    out = [1]
    for p in kp:
        out += [e * p for e in out]
    return out


# ---------------- Proposition 2.5 ----------------
class Base:
    def __init__(self, kp):
        self.kp = kp
        self.k = math.prod(kp) if kp else 1
        self.c = int(floor(sqrt(MMAX / self.k)))
        self.alpha = math.prod((q - 2) / (q - 1) for q in kp)
        self.main = CARTIN * math.prod((q * q - 2 * q) / (q * q - q - 1) for q in kp)
        divs = sqfree_divisors(kp)
        self.s1 = math.fsum(E(e * d * d) for d in range(1, self.c + 1)
                            if SQF[d] and coprime(d, kp) for e in divs)
        P, W = prefix(kp, self.k)
        self.s2 = P[ZBIG] - P[self.c] + 4.0 / ZBIG
        self.s3 = W[ZBIG] - W[self.c]
        assert self.k * ZBIG ** 2 < N0          # hypothesis k Z^2 < n0

    def E1(self, A):
        k, w = self.k, len(self.kp)
        return (self.s1 / L + N0 ** (A - 1) * log(k) / sqrt(k)
                + (self.s2 + 2 * self.s3 + 8 / (ZBIG * (1 - 2 * A))) * self.alpha
                + (2 ** w * N0 ** (A - 1) / sqrt(k)
                   + math.prod((q + 1) / q for q in self.kp) / (N0 ** A / sqrt(k) - 1)
                   + (2 ** w + (1 if k == 1 else 0)) / sqrt(N0)) * L)


# ---------------- Proposition 2.7 ----------------
class Peel:
    def __init__(self, q):
        self.q = q
        self.c1 = int(floor(sqrt(MMAX / q)))
        self.c2 = int(floor(10 ** 2.5 / q))
        assert q * self.c1 ** 2 <= MMAX and q * q * self.c2 ** 2 <= MMAX
        self.sh1 = math.fsum(E(q * d * d) for d in range(1, self.c1 + 1) if SQF[d] and d % q)
        self.sh2 = math.fsum(E(q * q * d * d) for d in range(1, self.c2 + 1) if SQF[d] and d % q)
        self.P, self.W1 = prefix((q,), q)
        _, self.W2 = prefix((q,), q * q)

    def Z_of(self, A):
        return int(floor(N0 ** A / self.q))

    def E2(self, A, Z=None):
        q = self.q
        Z1 = Z2 = (self.Z_of(A) if Z is None else Z)
        if Z1 > ZBIG:
            return float("inf")       # outside the precomputed range; never optimal
        assert Z1 > self.c1 and Z2 > self.c2
        P, W1, W2 = self.P, self.W1, self.W2
        s2_1 = P[Z1] - P[self.c1] + 4.0 / Z1
        s2_2 = P[Z2] - P[self.c2] + 4.0 / Z2
        s3_1 = W1[Z1] - W1[self.c1]
        s3_2 = W2[Z2] - W2[self.c2]
        e1 = self.sh1 / L + (2 / (q - 1)) * (s2_1 / 2 + s3_1 + 4 / (Z1 * (1 - 2 * A)))
        e2 = self.sh2 / L + (2 / (q * (q - 1))) * (s2_2 / 2 + s3_2 + 4 / (Z2 * (1 - 2 * A)))
        return e1 + e2 + (q * N0 ** (-2 * A) + N0 ** (-A) + N0 ** (A - 1) / q) * L


def argmin(f, lo, hi, steps=4000):
    xs = [lo + (hi - lo) * i / steps for i in range(steps + 1)]
    best = min(xs, key=f)
    h = (hi - lo) / steps
    a, b = max(lo, best - 2 * h), min(hi, best + 2 * h)
    for _ in range(200):          # ternary refinement inside the bracket
        m1, m2 = a + (b - a) / 3, b - (b - a) / 3
        if f(m1) <= f(m2):
            b = m2
        else:
            a = m1
    x = round((a + b) / 2, 6)     # the printed value is the one used
    return x, f(x)


def down(x, dp):
    return math.floor(x * 10 ** dp) / 10 ** dp


def up(x, dp):
    return math.ceil(x * 10 ** dp) / 10 ** dp


def primes_of(n):
    out, x, p = [], n, 2
    while p * p <= x:
        if x % p == 0:
            out.append(p)
            x //= p
        else:
            p += 1
    if x > 1:
        out.append(x)
    return tuple(out)


if __name__ == "__main__":
    print("Table 1 (Theorem 2.3), n0 = 2e12, c = floor(10^{5/2}/sqrt k), Z = 1e5")
    t1 = {}
    for k in (1, 3, 5, 7, 15, 21, 35, 105):
        b = Base(primes_of(k) if k > 1 else ())
        A, e = argmin(b.E1, 0.30, 0.4999)
        t1[k] = (A, b, b.main - e)
        print(f"  k={k:4d}  A={A:.6f}  Sigma1={up(b.s1, 5):.5f}  Sigma2={up(b.s2, 5):.5f}  "
              f"Sigma3={up(b.s3, 5):.5f}  alpha_k={down(b.main - e, 6):.6f}")

    rows = [(3, 138567), (3, 289731), (5, 230945), (7, 323323), (15, 36465), (15, 62985),
            (21, 51051), (21, 88179), (35, 85085), (105, 15015), (105, 23205)]
    peels = {q: Peel(q) for q in (11, 13, 17, 19, 23)}
    opt = {}
    print("\nPeeled primes: optimal A_q with Z_1 = Z_2 = floor(n0^{A_q}/q)")
    for q, pq in peels.items():
        A, e = argmin(pq.E2, 0.30, 0.4999)
        opt[q] = (A, e)
        print(f"  q={q}: A_q={A:.6f}  Z_1=Z_2={pq.Z_of(A)}  E_2={e:.6f}"
              f"   [same A_q, Z=1e5 as before: {pq.E2(A, ZBIG):.6f}]")

    print("\nTable 2 (corrected)")
    print("  m & k & A & (A_q) & LHS & RHS")
    for m, k in rows:
        mp = primes_of(m)
        qs = [q for q in primes_of(k) if m % q]
        A, e1 = t1[m][0], t1[m][1].E1(t1[m][0])
        lhs = (Fraction(1)
               * math.prod(Fraction(q * q - 2 * q, q * q - q - 1) for q in mp)
               - sum(Fraction(q - 1, q * q - q - 1) for q in qs)) * Fraction("0.3739558")
        rhs = log(2) / N0 + e1 + sum(opt[q][1] for q in qs)
        print(f"  ${m}$ & ${k}$ & ${A:.6f}$ & {', '.join(f'{opt[q][0]:.6f}' for q in qs)} & "
              f"${down(float(lhs), 6):.6f}$ & ${up(rhs, 6):.6f}$ \\\\")
