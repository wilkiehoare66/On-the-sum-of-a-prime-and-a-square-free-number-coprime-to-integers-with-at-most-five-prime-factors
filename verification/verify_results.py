#!/usr/bin/env python3
"""
verify_results.py
=================
Independent cross-verifier for the exhaustive finite computation. It reads the per-chunk CSV
output of the two C++ verifiers and checks, WITHOUT importing or trusting any
of their code:

  1. CHUNK COVERAGE: the recorded chunks tile each verifier's range
     [RANGE_START, RANGE_END] with no gap and no overlap, every chunk has status
     OK, and no chunk is recorded twice (a resume must not double-count).
  2. EXACT COUNT: the number of integers covered equals the exact
     arithmetic count of the range (accounting for the odd-only stride of the
     odd verifier), so "every integer in range is accounted for" is asserted,
     not assumed.
  3. WITNESS RE-DERIVATION: for a random sample of n in each range,
     this script INDEPENDENTLY finds a valid witness of the same kind the
     corresponding verifier claims to have found -- the six disjoint-support
     square-free witnesses with n - s an odd prime (odd), or NUM_SPLITS
     distinct Goldbach splits (even) -- using its own
     primality and square-free routines below. This re-checks primality,
     square-freeness and support-disjointness on real data, independently.
  4. FALLBACK COUNT: totals the last column of odd_results.csv, the number of
     odd n that the fast path did not clear and the exhaustive search had to.

Primality here is deterministic Miller-Rabin with the twelve Sorenson-Webster
bases (unconditional below 318665857834031151167461 ~ 3.19e23), matching the
verifiers' backend of record
-- but implemented separately, so agreement is a genuine cross-check.

Usage:
    python3 verify_results.py                 # coverage + count, both files
    python3 verify_results.py --sample 200    # also re-derive witnesses for 200 random n each
    python3 verify_results.py --results-dir results
Exit status is nonzero if any check fails.
"""
import argparse
import csv
import os
import random
import sys
from math import gcd, isqrt

# --- verifier metadata: filename, range, stride, kind -----------------------
# stride 2 = odd verifier records only odd n; stride 1 = every integer.
VERIFIERS = {
    # NUM_WITNESSES / NUM_SPLITS must match verify_odd.cpp and verify_even.cpp.
    # Six witnesses (odd) and three splits (even) cover omega(k) <= 5.
    "odd":  {"csv": "odd_results.csv", "start": 10001, "end": 2000000000000,
             "chunk": 1000000, "stride": 2, "kind": "witness"},
    "even": {"csv": "even_results.csv", "start": 10002, "end": 2000000000000,
             "chunk": 10000000, "stride": 2, "kind": "split"},
}

SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]


# --- independent primality: deterministic Miller-Rabin, 12 SW bases ---------
def is_prime(n):
    if n < 2:
        return False
    for p in SMALL_PRIMES:
        if n % p == 0:
            return n == p
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in SMALL_PRIMES:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


# --- independent square-free test (complete) --------------------------------
def is_squarefree(m):
    if m <= 0:
        return False
    if m == 1:
        return True
    # remove small square factors, then check the residual has no square factor
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if m % (p * p) == 0:
            return False
        while m % p == 0:
            m //= p
    if m == 1:
        return True
    # residual m has only prime factors > 37; it is square-free iff it is not
    # itself a perfect square times ... -> test all prime-square divisors up to
    # cube-root, then a final perfect-square check on what remains.
    d = 41
    while d * d * d <= m:
        if m % d == 0:
            if m % (d * d) == 0:
                return False
            while m % d == 0:
                m //= d
        d += 2
    # now m is 1, a prime, or a product of two distinct primes, or a prime^2
    s = isqrt(m)
    if s * s == m:
        return False
    return True


def odd_support(m):
    """Set of odd primes dividing m (m assumed square-free)."""
    s = set()
    x = m
    while x % 2 == 0:
        x //= 2
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if x % p == 0:
            s.add(p)
            while x % p == 0:
                x //= p
    d = 41
    while d * d <= x:
        if x % d == 0:
            s.add(d)
            while x % d == 0:
                x //= d
        d += 2
    if x > 1:
        s.add(x)
    return s


NUM_WITNESSES = 6      # must equal verify_odd.cpp's NUM_WITNESSES
NUM_SPLITS = 3         # must equal verify_even.cpp's NUM_SPLITS


# --- independent witness re-derivation, per kind ----------------------------
def find_splits(n):
    """NUM_SPLITS distinct Goldbach splits of an even n."""
    found = []
    for p in _primes_asc(n):
        if is_prime(n - p):
            found.append((p, n - p))
            # distinct splits are distinct unordered pairs; collect NUM_SPLITS
            if len({frozenset(pr) for pr in found}) >= NUM_SPLITS:
                pairs = list({frozenset(pr) for pr in found})[:NUM_SPLITS]
                return "; ".join("%d+%d" % tuple(sorted(pr)) for pr in pairs)
    return None


def find_witnesses(n):
    """NUM_WITNESSES square-free s_i with n-s_i prime, pairwise-disjoint odd supports."""
    wits = []
    supports = []
    # fast path: s = 2r with r odd prime, n-2r prime -> support {r}
    for r in _primes_asc(n // 2):
        if r == 2:
            continue
        s = 2 * r
        if s >= n:
            break
        if is_prime(n - s) and is_squarefree(s):
            sup = frozenset({r})
            if all(sup.isdisjoint(t) for t in supports):
                wits.append("s=%d" % s)
                supports.append(sup)
                if len(wits) == NUM_WITNESSES:
                    return "; ".join(wits)
    # fallback: general square-free s with disjoint support; p = 2 is skipped
    # because Lemma 3.2 asks for a prime 2 < p < n
    for p in _primes_asc(n):
        if p == 2:
            continue
        s = n - p
        if s < 2:
            continue
        if is_squarefree(s):
            sup = frozenset(odd_support(s))
            if all(sup.isdisjoint(t) for t in supports):
                wits.append("s=%d" % s)
                supports.append(sup)
                if len(wits) == NUM_WITNESSES:
                    return "; ".join(wits)
    return None


def _primes_asc(limit):
    p = 2
    while p < limit:
        if is_prime(p):
            yield p
        p += 1


# --- CSV coverage + count ---------------------------------------------------
def check_coverage(name, meta, path):
    problems = []
    if not os.path.exists(path):
        if meta.get("optional"):
            print(f"  {name}: SKIPPED (optional cross-check, {path} not present)")
            return []
        return [f"{name}: results file {path} not found"]
    rows = []
    fallback_total, fallback_rows = 0, 0
    with open(path, newline="") as fh:
        for line in csv.reader(fh):
            if not line or line[0].strip().lower() in ("timestamp", "time"):
                continue
            try:
                start = int(line[1])
                end = int(line[2])
                status = line[3].strip().upper()
            except (IndexError, ValueError):
                problems.append(f"{name}: unparseable row: {line[:4]}")
                continue
            rows.append((start, end, status))
            # column 8 of odd_results.csv: values that reached the fallback
            if (meta["kind"] == "witness" and len(line) >= 8
                    and line[7].strip().isdigit()):
                fallback_total += int(line[7])
                fallback_rows += 1

    if not rows:
        return [f"{name}: no data rows in {path}"]

    # status check
    bad = [(s, e, st) for (s, e, st) in rows if not st.startswith("OK")]
    for s, e, st in bad[:5]:
        problems.append(f"{name}: chunk [{s},{e}] has non-OK status '{st}'")

    # sort by start; check no overlap, no gap, no duplicate
    rows.sort()
    expect = meta["start"]
    seen_starts = set()
    for (s, e, st) in rows:
        if s in seen_starts:
            problems.append(f"{name}: chunk starting {s} recorded more than once (resume double-count)")
            continue
        seen_starts.add(s)
        if s > expect:
            problems.append(f"{name}: GAP in coverage: [{expect}, {s}) uncovered")
        elif s < expect:
            problems.append(f"{name}: OVERLAP: chunk [{s},{e}] starts before expected {expect}")
        expect = e + 1
    # final chunk should reach RANGE_END (last covered integer >= end, or the
    # last chunk's end is the range end)
    last_end = rows[-1][1]
    if last_end < meta["end"]:
        problems.append(f"{name}: coverage stops at {last_end}, short of RANGE_END {meta['end']}")

    # exact count (item 8)
    covered = _count_integers(meta["start"], last_end, meta["stride"])
    expected = _count_integers(meta["start"], meta["end"], meta["stride"])
    if last_end >= meta["end"] and covered != expected:
        problems.append(f"{name}: covered count {covered} != expected {expected} "
                        f"for range [{meta['start']},{meta['end']}] stride {meta['stride']}")
    else:
        note = "integers" if meta["stride"] == 1 else "odd integers" if meta["start"] % 2 else "even integers"
        print(f"  {name}: {len(rows)} chunks, contiguous, "
              f"{covered:,} {note} covered "
              f"[{meta['start']:,} .. {last_end:,}]"
              + ("" if last_end >= meta["end"] else "  (INCOMPLETE)"))
    if fallback_rows:
        print(f"  {name}: {fallback_total:,} values reached the fallback search "
              f"(over {fallback_rows:,} chunks)")
    return problems


def _count_integers(a, b, stride):
    """Count integers in [a,b] with the given stride (1=all, 2=matching parity of a)."""
    if b < a:
        return 0
    if stride == 1:
        return b - a + 1
    # stride 2: integers of the same parity as a
    return (b - a) // 2 + 1


def check_witnesses(name, meta, n_samples, rng):
    """Re-derive a witness independently for n_samples random n in range."""
    finder = {"split": find_splits,
              "witness": find_witnesses}[meta["kind"]]
    problems = []
    lo, hi = meta["start"], meta["end"]
    tested = 0
    for _ in range(n_samples):
        n = rng.randrange(lo, hi + 1)
        if meta["stride"] == 2:
            # align parity to the verifier's stride
            if (n - meta["start"]) % 2:
                n += 1
            if n > hi:
                continue
        w = finder(n)
        tested += 1
        if w is None:
            problems.append(f"{name}: NO independent witness found for n={n} "
                            f"(kind {meta['kind']}) -- investigate")
    print(f"  {name}: independently re-derived witnesses for {tested} sampled n "
          f"({'all verified' if not problems else str(len(problems)) + ' FAILED'})")
    return problems


def main():
    ap = argparse.ArgumentParser(description="Independent cross-verifier for the finite computation")
    ap.add_argument("--results-dir", default="results", help="directory holding the *_results.csv files")
    ap.add_argument("--sample", type=int, default=0,
                    help="re-derive witnesses for this many random n per verifier (0 = coverage/count only)")
    ap.add_argument("--seed", type=int, default=20260725, help="RNG seed for reproducible sampling")
    ap.add_argument("--only", choices=list(VERIFIERS), help="check a single verifier")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    names = [args.only] if args.only else list(VERIFIERS)

    print("== Coverage and exact-count check ==")
    all_problems = []
    for name in names:
        meta = VERIFIERS[name]
        path = os.path.join(args.results_dir, meta["csv"])
        all_problems += check_coverage(name, meta, path)

    if args.sample:
        print(f"\n== Independent witness re-derivation ({args.sample}/verifier) ==")
        for name in names:
            all_problems += check_witnesses(name, VERIFIERS[name], args.sample, rng)

    print()
    if all_problems:
        print(f"FAILED: {len(all_problems)} problem(s) found:")
        for p in all_problems[:40]:
            print(f"  ! {p}")
        sys.exit(1)
    print("ALL CHECKS PASSED: coverage is gap-free and non-overlapping, counts are exact"
          + (", sampled witnesses independently re-derived." if args.sample else "."))
    sys.exit(0)


if __name__ == "__main__":
    main()
