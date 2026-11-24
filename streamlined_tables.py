#!/usr/bin/env python3
"""Tables 1 and 2 of the paper, and the domination check in the proof of
Corollary 3.6, at n0 = 2e12, from the Bennett--Martin--O'Bryant--Rechnitzer
c_theta table.

Method (as described in the remark after Table 2):

  * The parameters are chosen with scipy's differential_evolution, working in
    double precision, and then rounded to six decimal places:
      - Table 1: A for each k, minimising E_1 (Proposition 2.5).
      - Table 2: for each row, (A, A_q for q | k/m) jointly, minimising the
        right-hand side of (eqn:suff_cond_lr).
      - Domination check: A_q for each prime 11 <= q <= 37 on its own,
        minimising E_2 (Proposition 2.7).
  * Every bound is then recomputed at the rounded parameters with mpmath to
    at least 120 decimal places (mp.dps = 135; no quantity used exceeds
    n0 = 2e12 in size). Lower bounds are rounded down and upper bounds up, at
    the last printed digit.

The script then checks the values printed in the paper directly: Table 1 and
Table 2 at the parameters printed there, and the domination step with the
A_q printed in Table 2.

Parameters: n0 = 2e12;  c = floor(sqrt(1e5/m)), Z = 1e5 for the base m;
c_1 = floor(sqrt(1e5/q)), c_2 = floor(10^{5/2}/q),
Z_1 = Z_2 = min(floor(n0^{A_q}/q), 1e5) for each peeled prime q.

Usage:  python3 streamlined_tables.py bennett_c_theta.tsv [--x0 bennett_x0.tsv] [--seed 0]
Needs:  Python 3.8+, scipy, mpmath.
"""
import argparse
import math
import os
import sys
from math import isqrt

import mpmath
import scipy
from mpmath import mp, mpf
from scipy.optimize import differential_evolution

mp.dps = 135   # >= 120 decimal places for every quantity (all are below 1e13)

N0 = 2 * 10**12
MMAX = 10**5          # Bennett et al. table covers moduli 3 <= m <= 1e5
ZBIG = 10**5          # Z for the base, and cap on Z_1, Z_2
ABOUNDS = (0.30, 0.4999)
CARTIN = mpf("0.3739558")  # lower bound for Artin's constant 0.37395581...

LF = math.log(N0)
LM = mp.log(N0)

# ----------------------------------------------------------------------------
# Elementary arithmetic up to ZBIG
# ----------------------------------------------------------------------------
spf = list(range(ZBIG + 1))
for i in range(2, isqrt(ZBIG) + 1):
    if spf[i] == i:
        for j in range(i * i, ZBIG + 1, i):
            if spf[j] == j:
                spf[j] = i
SQF = [False] * (ZBIG + 1)
PHI2 = [0] * (ZBIG + 1)      # phi(d^2) = d * phi(d)
PR = [()] * (ZBIG + 1)       # distinct prime factors
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


def is_prime(n):
    if n <= ZBIG:
        return n >= 2 and spf[n] == n
    return all(n % p for p in range(2, isqrt(n) + 1))


def coprime(d, kp):
    return all(p not in kp for p in PR[d])


def primes_of(n):
    out, x, p = [], n, 2
    while p * p <= x:
        if x % p == 0:
            out.append(p)
            while x % p == 0:
                x //= p
        p += 1
    if x > 1:
        out.append(x)
    return tuple(out)


def sqfree_divisors(kp):
    out = [1]
    for p in kp:
        out += [e * p for e in out]
    return out


# ----------------------------------------------------------------------------
# Bennett et al. constants
# ----------------------------------------------------------------------------
CT_F, CT_M = {}, {}


def load_ctheta(path):
    with open(path) as fh:
        next(fh)
        for line in fh:
            a, b = line.split()
            CT_F[int(a)] = float(b)
            CT_M[int(a)] = mpf(b)


def check_x0(path):
    """Theorem 2.2 is applied for n >= n0; check x_theta(m) <= n0 for every m."""
    worst = 0
    with open(path) as fh:
        next(fh)
        for line in fh:
            a, b = line.split()
            worst = max(worst, float(b))
    assert worst <= N0, f"x_theta exceeds n0: {worst}"
    return worst


def E_f(m):
    assert m <= MMAX
    return 0.375 / LF**2 if m == 1 else CT_F[m]


def E_m(m):
    assert m <= MMAX
    return mpf("0.375") / LM**2 if m == 1 else CT_M[m]


# ----------------------------------------------------------------------------
# High-precision sums over square-free d in (lo, hi]
# ----------------------------------------------------------------------------
_INV, _LOGD, _SUMS = {}, {}, {}


def _inv(d):
    if d not in _INV:
        _INV[d] = mpf(1) / PHI2[d]
    return _INV[d]


def _logd(d):
    if d not in _LOGD:
        _LOGD[d] = mp.log(d)
    return _LOGD[d]


def mp_sums(kp, kval, lo, hi):
    """(sum 1/phi(d^2), sum 1/phi(d^2) * (1 - log(kval d^2)/log n0)^{-1}) over
    square-free d in (lo, hi] coprime to the primes kp."""
    key = (kp, kval, lo, hi)
    if key in _SUMS:
        return _SUMS[key]
    assert kval * hi * hi < N0, "Brun--Titchmarsh weight undefined"
    lk = mp.log(kval)
    s = w = mpf(0)
    for d in range(lo + 1, hi + 1):
        if SQF[d] and coprime(d, kp):
            t = _inv(d)
            s += t
            w += t / (1 - (lk + 2 * _logd(d)) / LM)
    _SUMS[key] = (s, w)
    return s, w


def float_prefix(kp, kval):
    """Double-precision prefix sums used only inside the optimiser."""
    P, W = [0.0] * (ZBIG + 1), [0.0] * (ZBIG + 1)
    s = w = 0.0
    lk = math.log(kval)
    for d in range(1, ZBIG + 1):
        if SQF[d] and coprime(d, kp):
            t = 1.0 / PHI2[d]
            s += t
            if kval * d * d < N0:
                w += t / (1 - (lk + 2 * math.log(d)) / LF)
            else:
                w = float("nan")
        P[d], W[d] = s, w
    return P, W


# ----------------------------------------------------------------------------
# Proposition 2.5: E_1(m, n0, c, A, Z)
# ----------------------------------------------------------------------------
class Base:
    def __init__(self, m):
        self.m = m
        self.kp = primes_of(m) if m > 1 else ()
        self.w = len(self.kp)
        self.c = isqrt(MMAX // m)                       # floor(sqrt(1e5/m))
        self.Z = ZBIG
        assert self.c >= 1 and self.Z > self.c
        assert m * self.Z**2 < N0                       # k Z^2 < n0
        divs = sqfree_divisors(self.kp)
        ds = [d for d in range(1, self.c + 1) if SQF[d] and coprime(d, self.kp)]
        assert all(e * d * d <= MMAX for d in ds for e in divs)
        self.s1_f = math.fsum(E_f(e * d * d) for d in ds for e in divs)
        self.s1_m = mp.fsum(E_m(e * d * d) for d in ds for e in divs)
        P, W = float_prefix(self.kp, m)
        self.s2_f = P[self.Z] - P[self.c] + 4.0 / self.Z
        self.s3_f = W[self.Z] - W[self.c]
        s, w = mp_sums(self.kp, m, self.c, self.Z)
        self.s2_m = s + mpf(4) / self.Z
        self.s3_m = w
        self.alpha_m = mp.fprod([mpf(q - 2) / (q - 1) for q in self.kp])
        self.alpha_f = float(self.alpha_m)
        self.main_m = CARTIN * mp.fprod([mpf(q * q - 2 * q) / (q * q - q - 1) for q in self.kp])
        self.up_m = mp.fprod([mpf(q + 1) / q for q in self.kp])
        self.up_f = float(self.up_m)

    def E1_f(self, A):
        k, w = self.m, self.w
        rk = math.sqrt(k)
        return (self.s1_f / LF + N0 ** (A - 1) * math.log(k) / rk
                + (self.s2_f + 2 * self.s3_f + 8 / (self.Z * (1 - 2 * A))) * self.alpha_f
                + (2**w * N0 ** (A - 1) / rk
                   + self.up_f / (N0**A / rk - 1)
                   + (2**w + (k == 1)) / math.sqrt(N0)) * LF)

    def E1_m(self, A):
        A = mpf(A)
        k, w = self.m, self.w
        rk, n = mp.sqrt(k), mpf(N0)
        assert self.c * rk <= n**A                      # c sqrt(k) <= n0^A
        return (self.s1_m / LM + n ** (A - 1) * mp.log(k) / rk
                + (self.s2_m + 2 * self.s3_m + 8 / (self.Z * (1 - 2 * A))) * self.alpha_m
                + (2**w * n ** (A - 1) / rk
                   + self.up_m / (n**A / rk - 1)
                   + (2**w + (1 if k == 1 else 0)) / mp.sqrt(n)) * LM)


# ----------------------------------------------------------------------------
# Proposition 2.7: E_2(q, n0, c_1, c_2, A, Z_1, Z_2)
# ----------------------------------------------------------------------------
class Peel:
    def __init__(self, q):
        assert q % 2 and is_prime(q)
        self.q = q
        self.c1 = isqrt(MMAX // q)                      # floor(sqrt(1e5/q))
        self.c2 = isqrt(MMAX // (q * q))                # floor(10^{5/2}/q)
        assert self.c1 >= 1 and self.c2 >= 1
        assert q * self.c1**2 <= MMAX and q * q * self.c2**2 <= MMAX
        d1 = [d for d in range(1, self.c1 + 1) if SQF[d] and d % q]
        d2 = [d for d in range(1, self.c2 + 1) if SQF[d] and d % q]
        self.sh1_f = math.fsum(E_f(q * d * d) for d in d1)
        self.sh2_f = math.fsum(E_f(q * q * d * d) for d in d2)
        self.sh1_m = mp.fsum(E_m(q * d * d) for d in d1)
        self.sh2_m = mp.fsum(E_m(q * q * d * d) for d in d2)
        self.P, self.W1 = float_prefix((q,), q)
        _, self.W2 = float_prefix((q,), q * q)

    def Z_f(self, A):
        return min(int(math.floor(N0**A / self.q)), ZBIG)

    def Z_m(self, A):
        return min(int(mp.floor(mpf(N0) ** mpf(A) / self.q)), ZBIG)

    def E2_f(self, A):
        q, Z = self.q, self.Z_f(A)
        if Z <= max(self.c1, self.c2):
            return float("inf")
        s2_1 = self.P[Z] - self.P[self.c1] + 4.0 / Z
        s2_2 = self.P[Z] - self.P[self.c2] + 4.0 / Z
        s3_1 = self.W1[Z] - self.W1[self.c1]
        s3_2 = self.W2[Z] - self.W2[self.c2]
        if math.isnan(s3_1) or math.isnan(s3_2):
            return float("inf")
        e1 = self.sh1_f / LF + (2 / (q - 1)) * (s2_1 / 2 + s3_1 + 4 / (Z * (1 - 2 * A)))
        e2 = self.sh2_f / LF + (2 / (q * (q - 1))) * (s2_2 / 2 + s3_2 + 4 / (Z * (1 - 2 * A)))
        return e1 + e2 + (q * N0 ** (-2 * A) + N0 ** (-A) + N0 ** (A - 1) / q) * LF

    def E2_m(self, A):
        q, Z = self.q, self.Z_m(A)
        A, n = mpf(A), mpf(N0)
        assert Z > self.c1 and Z > self.c2
        assert q * Z <= n**A                            # q Z_i <= n0^A
        s1, w1 = mp_sums((q,), q, self.c1, Z)
        s2, w2 = mp_sums((q,), q * q, self.c2, Z)
        s2_1, s2_2 = s1 + mpf(4) / Z, s2 + mpf(4) / Z
        e1 = self.sh1_m / LM + (mpf(2) / (q - 1)) * (s2_1 / 2 + w1 + 4 / (Z * (1 - 2 * A)))
        e2 = self.sh2_m / LM + (mpf(2) / (q * (q - 1))) * (s2_2 / 2 + w2 + 4 / (Z * (1 - 2 * A)))
        return e1 + e2 + (q * n ** (-2 * A) + n ** (-A) + n ** (A - 1) / q) * LM

    def h_m(self, A):
        """Contribution of q in the proof of Corollary 3.6."""
        q = self.q
        return CARTIN * mpf(q - 1) / (q * q - q - 1) + self.E2_m(A)


# ----------------------------------------------------------------------------
# Lemma 2.6 / Corollary 3.5: c_q(n0)
# ----------------------------------------------------------------------------
def c_q(q):
    if q <= 100:
        return mpf(1) / q + 1 / (420 * LM)
    if q <= 10**5:
        return mpf(1) / (q - 1) + 1 / (160 * LM)
    return mpf("0.000041")


def c_q_f(q):
    if q <= 100:
        return 1 / q + 1 / (420 * LF)
    if q <= 10**5:
        return 1 / (q - 1) + 1 / (160 * LF)
    return 0.000041


# ----------------------------------------------------------------------------
# Optimisation and rounding
# ----------------------------------------------------------------------------
def de_min(f, dim, seed):
    res = differential_evolution(lambda x: f(x), [ABOUNDS] * dim, seed=seed,
                                 tol=1e-14, atol=0, maxiter=5000, popsize=40,
                                 mutation=(0.5, 1.0), recombination=0.9,
                                 polish=False, init="sobol")
    return [round(float(v), 6) for v in res.x]


def down(x, dp):
    v = int(mp.floor(x * 10**dp))
    return f"{v // 10**dp}.{v % 10**dp:0{dp}d}"


def up(x, dp):
    v = int(mp.ceil(x * 10**dp))
    return f"{v // 10**dp}.{v % 10**dp:0{dp}d}"


GAPS = []


def grid_check(f, A_de, steps=20000):
    """Sanity check on the optimiser, for the one-dimensional searches only:
    how much a 20001-point grid improves on the rounded DE optimum. Optimality
    is not needed for any bound to be valid; this only confirms the search."""
    lo, hi = ABOUNDS
    best = min(f(lo + (hi - lo) * i / steps) for i in range(steps + 1))
    GAPS.append(max(0.0, f(A_de) - best))


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
TABLE1_K = (1, 3, 5, 7, 15, 21, 35, 105)

# Values as printed in the paper.
# Theorem 2.3 and Table 1 (tab:detailed_comps): k -> (A, Sigma_1, Sigma_2, Sigma_3, alpha_k)
PAPER_T1 = {
    1: ("0.417301", "0.59123", "0.00322", "0.00619", "0.336775"),
    3: ("0.443722", "0.40654", "0.00367", "0.0071", "0.200473"),
    5: ("0.439804", "0.32124", "0.00561", "0.01092", "0.262482"),
    7: ("0.440077", "0.27167", "0.00714", "0.01392", "0.279591"),
    15: ("0.461668", "0.20577", "0.00669", "0.013", "0.15678"),
    21: ("0.461901", "0.1867", "0.00836", "0.01632", "0.166938"),
    35: ("0.458742", "0.15419", "0.01368", "0.02661", "0.203566"),
    105: ("0.475597", "0.10578", "0.0167", "0.0323", "0.120652"),
}
# Table 2 (tab:remaining_cases): (m, k, A, (A_q for q | k/m, increasing q), LHS, RHS)
PAPER_T2 = [
    (3, 138567, "0.443717", ("0.467435", "0.468303", "0.468351", "0.468681"), "0.119296", "0.070475"),
    (3, 289731, "0.443601", ("0.468152", "0.468516", "0.468748", "0.468869"), "0.137312", "0.066038"),
    (5, 230945, "0.439830", ("0.467916", "0.468116", "0.468559", "0.468827"), "0.190150", "0.079321"),
    (7, 323323, "0.440037", ("0.467856", "0.468062", "0.468666", "0.468659"), "0.214153", "0.086214"),
    (15, 36465, "0.461550", ("0.467945", "0.468108", "0.468582"), "0.091799", "0.056884"),
    (15, 62985, "0.461715", ("0.468068", "0.468585", "0.468776"), "0.106367", "0.053387"),
    (21, 51051, "0.461936", ("0.467830", "0.468094", "0.468672"), "0.106200", "0.061128"),
    (21, 88179, "0.461906", ("0.468284", "0.468527", "0.468819"), "0.120768", "0.057630"),
    (35, 85085, "0.458677", ("0.468142", "0.468405", "0.468491"), "0.166686", "0.084985"),
    (105, 15015, "0.475572", ("0.468015", "0.468364"), "0.087955", "0.056593"),
    (105, 23205, "0.475571", ("0.468355", "0.468543"), "0.100184", "0.053545"),
]
TABLE2_ROWS = [(3, 138567), (3, 289731), (5, 230945), (7, 323323), (15, 36465), (15, 62985),
               (21, 51051), (21, 88179), (35, 85085), (105, 15015), (105, 23205)]
SMALL_Q = [q for q in range(11, 38) if is_prime(q)]     # q <= 37: through h(q)
Q_LARGE = 41                                            # q >= 41: through c_q(n0)


def domination(rows, peels, hopt, label):
    """Proof of Corollary 3.6: for each row and each prime q* of k/m other than
    the least, every prime q > q* contributes at most h(q*) (taken at the row's
    A_q): through h(q) at its own optimal A_q for q <= 37, and through
    c_q(n0) <= c_41(n0) for q >= 41."""
    ok, worst = True, None
    for m, k, qs, A, Aqs, lhs, rhs in rows:
        hrow = {q: peels[q].h_m(aq) for q, aq in zip(qs, Aqs)}
        for qstar in qs[1:]:
            cands = [(q, hopt[q]) for q in SMALL_Q if q > qstar] + [(Q_LARGE, c_q(Q_LARGE))]
            for q, val in cands:
                margin = hrow[qstar] - val
                if worst is None or margin < worst[0]:
                    worst = (margin, m, k, qstar, q)
                if margin < 0:
                    ok = False
                    print(f"  FAILS: row m={m}, k={k}: q={q} exceeds q*={qstar}")
    margin, m, k, qstar, q = worst
    print(f"  every prime q > q* contributes at most h(q*), {label}: "
          f"{'OK' if ok else 'FAILS'} (tightest: row m={m}, k={k}, "
          f"q* = {qstar} against q = {q}, margin {mp.nstr(margin, 6)})")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ctheta", help="bennett_c_theta.tsv")
    ap.add_argument("--x0", default=None, help="bennett_x0.tsv (default: beside ctheta)")
    ap.add_argument("--seed", type=int, default=0, help="seed for differential_evolution")
    args = ap.parse_args()

    load_ctheta(args.ctheta)
    x0path = args.x0 or os.path.join(os.path.dirname(os.path.abspath(args.ctheta)), "bennett_x0.tsv")
    print(f"python {sys.version.split()[0]}, scipy {scipy.__version__}, mpmath {mpmath.__version__}, "
          f"mp.dps = {mp.dps}, seed = {args.seed}, n0 = {N0:.0e}")
    if os.path.exists(x0path):
        print(f"max x_theta(m) over 3 <= m <= 1e5: {check_x0(x0path):.0f} <= n0  [OK]")
    ok = True

    # ------------------------------------------------------------- Table 1
    print("\nTable 1 (Theorem 2.3 and tab:detailed_comps): c = floor(sqrt(1e5/k)), Z = 1e5")
    print("  k & A & Sigma_1 & Sigma_2 & Sigma_3 & alpha_k")
    bases = {}
    for k in TABLE1_K:
        b = Base(k)
        (A,) = de_min(lambda x: b.E1_f(x[0]), 1, args.seed)
        grid_check(b.E1_f, A)
        alpha = b.main_m - b.E1_m(A)
        bases[k] = (b, A, alpha)
        print(f"  ${k}$ & ${A:.6f}$ & ${up(b.s1_m, 5)}$ & ${up(b.s2_m, 5)}$ & ${up(b.s3_m, 5)}$ "
              f"& ${down(alpha, 6)}$ \\\\")

    # ------------------------------------------------------------- Table 2
    print("\nTable 2: (A, A_q) jointly by differential evolution for each row, "
          "Z_1 = Z_2 = min(floor(n0^{A_q}/q), 1e5)")
    print("  m & k & A & (A_q) & LHS & RHS")
    peels = {}
    rowdata = []
    for idx, (m, k) in enumerate(TABLE2_ROWS):
        b = bases[m][0]
        qs = [q for q in primes_of(k) if m % q]
        for q in qs:
            peels.setdefault(q, Peel(q))
        f = lambda x: b.E1_f(x[0]) + sum(peels[q].E2_f(a) for q, a in zip(qs, x[1:]))
        xs = de_min(f, 1 + len(qs), args.seed + idx)
        A, Aqs = xs[0], xs[1:]
        a = mp.fprod([mpf(q * q - 2 * q) / (q * q - q - 1) for q in b.kp])
        s = mp.fsum([mpf(q - 1) / (q * q - q - 1) for q in qs])
        lhs = (a - s) * CARTIN
        rhs = mp.log(2) / N0 + b.E1_m(A) + mp.fsum([peels[q].E2_m(aq) for q, aq in zip(qs, Aqs)])
        good = (a > s) and (lhs > rhs)
        ok &= good
        rowdata.append((m, k, qs, A, Aqs, lhs, rhs))
        print(f"  ${m}$ & ${k}$ & ${A:.6f}$ & {', '.join(f'{x:.6f}' for x in Aqs)} & "
              f"${down(lhs, 6)}$ & ${up(rhs, 6)}$ \\\\" + ("" if good else "   [FAILS]"))

    # ------------------------------------------------- Domination check
    print("\nDomination check (proof of Corollary 3.6): h(q) = C(q-1)/(q^2-q-1) + E_2(q) "
          "for 11 <= q <= 37 at its own optimal A_q; c_q(n0) for q >= 41")
    hopt = {}
    for q in SMALL_Q:
        pq = peels.setdefault(q, Peel(q))
        (Aq,) = de_min(lambda x: pq.E2_f(x[0]), 1, args.seed)
        grid_check(pq.E2_f, Aq)
        hopt[q] = pq.h_m(Aq)
        print(f"  q = {q:2d}: A_q = {Aq:.6f}  Z = {pq.Z_m(Aq):6d}  h(q) = {up(hopt[q], 6)}")

    # c_q(n0) is non-increasing in q over all odd primes: checked up to 1e5 in
    # double precision (consecutive gaps are far above rounding error), and
    # c_q = 0.000041 beyond 1e5 is below every earlier value.
    plist = [p for p in range(3, ZBIG + 1) if is_prime(p)]
    mono = all(c_q_f(plist[i]) > c_q_f(plist[i + 1]) for i in range(len(plist) - 1))
    mono &= c_q(plist[-1]) > c_q(ZBIG + 3)
    ok &= mono
    print(f"  c_q(n0) non-increasing over all odd primes q: {'OK' if mono else 'FAILS'}; "
          f"c_41(n0) = {up(c_q(Q_LARGE), 6)}")

    ok &= domination(rowdata, peels, hopt, "with the A_q found above")

    # ------------------------------------------ values printed in the paper
    print("\nValues printed in the paper")
    print("  Table 1 and Theorem 2.3, at the printed A: printed Sigma_i >= computed, "
          "printed alpha_k <= computed")
    for k in TABLE1_K:
        b = bases[k][0]
        A, S1, S2, S3, al = PAPER_T1[k]
        comp = (b.s1_m, b.s2_m, b.s3_m)
        alpha = b.main_m - b.E1_m(A)
        good = all(mpf(x) >= c for x, c in zip((S1, S2, S3), comp)) and mpf(al) <= alpha
        ok &= good
        print(f"    k = {k:3d}: {'OK' if good else 'FAILS'}  (computed alpha_k = {down(alpha, 6)}, "
              f"printed {al})")
    print("  Table 2, at the printed A and A_q: printed LHS <= computed LHS, "
          "printed RHS >= computed RHS, printed LHS > printed RHS")
    prow = []
    for m, k, A, Aqs, L_p, R_p in PAPER_T2:
        b = bases[m][0]
        qs = [q for q in primes_of(k) if m % q]
        assert len(qs) == len(Aqs)
        a = mp.fprod([mpf(q * q - 2 * q) / (q * q - q - 1) for q in b.kp])
        s = mp.fsum([mpf(q - 1) / (q * q - q - 1) for q in qs])
        lhs = (a - s) * CARTIN
        rhs = mp.log(2) / N0 + b.E1_m(A) + mp.fsum([peels[q].E2_m(aq) for q, aq in zip(qs, Aqs)])
        good = (a > s) and mpf(L_p) <= lhs and mpf(R_p) >= rhs and mpf(L_p) > mpf(R_p)
        ok &= good
        same = down(lhs, 6) == L_p and up(rhs, 6) == R_p
        prow.append((m, k, qs, A, Aqs, lhs, rhs))
        print(f"    m = {m:3d}, k = {k:6d}: {'OK' if good else 'FAILS'}  "
              f"(LHS {down(lhs, 6)}, RHS {up(rhs, 6)}"
              f"{', identical to printed' if same else ''})")
    ok &= domination(prow, peels, hopt, "with the A_q printed in Table 2")

    print(f"\nOptimiser sanity check (1-D searches): a 20001-point grid improves on the rounded "
          f"DE optimum by at most {max(GAPS):.1e} (below half a unit in the sixth decimal: "
          f"{'yes' if max(GAPS) < 5e-7 else 'NO'})")
    print("\nALL CHECKS PASSED" if ok else "\nSOME CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
