# On the sum of a prime and a square-free number coprime to integers with at most five prime factors

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21834718.svg)](https://doi.org/10.5281/zenodo.21813035)

Code and data for the computational verification in the paper *"On the sum of
a prime and a square-free number coprime to integers with at most five prime
factors"* by W. Hoare, E. S. Lee and A. Pearce-Crump.

## What is verified

Lemma 3.2 of the paper: for every odd square-free `k` with `omega(k) <= 5`,
every integer `108 <= n <= 2e12`, of either parity, is a sum `p + s` of a prime
`2 < p < n` and a square-free number `s` coprime to `k`. The analytic part of
the paper covers `n > 2e12`. The prime is required to exceed 2 because Lemma
3.1 of the paper derives the even-`k` case from this one.

The range splits at `1e4`:

* `108 <= n <= 1e4`: `Threshold.py` decides the question exactly, by a
  hitting-set criterion (Section 6, *The sharp threshold*).
* `1e4 < n <= 2e12`: two C++ verifiers exhibit witnesses. `verify_odd` finds,
  for each odd `n`, six square-free `s` with `n - s` prime and pairwise coprime
  odd parts; `verify_even` finds, for each even `n`, three distinct Goldbach
  splits. Either way at most five of the six survive the prime factors of `k`.

`Threshold.py` also computes the sharp thresholds quoted in the remark after
Theorem 1.1: 36, 60, 84 and 108 for `omega(k) <= 2, 3, 4, 5`.

## Requirements

* Python 3.8+ (standard library only)
* A C++17 compiler with OpenMP (`g++ -fopenmp`)

## Layout

```
.
├── Threshold.py              exact criterion: 108 <= n <= 1e4, and the thresholds
├── README.md
└── verification/
    ├── verify_odd.cpp         odd n: six square-free witnesses, pairwise coprime odd parts
    ├── verify_even.cpp        even n: three distinct Goldbach splits
    ├── prime64.hpp            Baillie-PSW backend for the optional -DUSE_BPSW builds
    ├── verify_results.py      independent cross-check of the run's output
    ├── run_manifest.sh        build and run provenance
    └── results/               output of the run
        ├── odd_results.csv
        ├── even_results.csv
        ├── even_deep_checks.csv
        ├── odd_run.log.gz
        ├── even_run.log.gz
        └── run_manifest.txt
```

## What corresponds to what

| paper | produced by | command |
|---|---|---|
| Remark after Theorem 1.1 (thresholds) | `Threshold.py` | `--rmin 2 --rmax 5 --nmax 3000` |
| Lemma 3.2, `108 <= n <= 1e4` | `Threshold.py` | `--rmax 5 --certify 108 10000` |
| Lemma 3.2, odd `n` in `(1e4, 2e12]` | `verify_odd.cpp` | see *The C++ run* |
| Lemma 3.2, even `n` in `(1e4, 2e12]` | `verify_even.cpp` | see *The C++ run* |
| Section 6, *Reproducibility* | `verify_results.py` | `--sample 200` |

## The range 108 <= n <= 1e4

For each `n`, put `W(n) = { n - p : p < n prime, n - p square-free }` and let
`supp(w)` be the set of odd primes dividing `w`. An odd `k` admits no
representation of `n` exactly when its prime factors meet every `supp(w)`, so
`n` is good for every odd `k` with `omega(k) <= r` exactly when the least
hitting set of `{ supp(w) }` has more than `r` elements.

```
python3 Threshold.py --rmax 5 --certify 108 10000     # seconds
```

applies this to Lemma 3.2 itself: odd `k`, both parities of `n`, and only
primes `2 < p < n`. It should print

```
exact criterion, omega(k) <= 5, every n in [108, 10000]:
  ALL 9893 VALUES CERTIFIED (every odd k, both parities of n, primes 2 < p < n only)
```

The threshold mode admits `p = 2`, as Theorem 1.1 does, and also handles even
`k` (the summand must then be odd). Excluding `p = 2` leaves every threshold
unchanged.

```
python3 Threshold.py --rmin 2 --rmax 5 --nmax 3000
```

```
exact sharp thresholds, every n <= 3000 examined
  r  largest failing n      n_0   defeating modulus
  2                 35       36   n = 35, k = 33 = 3*11
  3                 59       60   n = 59, k = 759 = 3*11*23
  4                 83       84   n = 83, k = 3795 = 3*5*11*23
  5                107      108   n = 107, k = 275655 = 3*5*17*23*47
```

The `r = 2` row reproduces the threshold 36 of Lee and O'Clarey, with the same
integer 33 attaining it, which is the external check on the script.

The C++ device takes over above `1e4` because it is sufficient but not
necessary, and below `1e4` the gap shows: at `n = 167` no six pairwise coprime
witnesses exist, yet `n = 167` is good for every `omega(k) <= 5`.

## The C++ run

Build without `-DUSE_BPSW`: the results of record use the twelve-base
Miller-Rabin test of Sorenson and Webster, which is unconditional below
`3.19e23`. The Baillie-PSW variant is faster but exhaustive only conjecturally.

```
cd verification
g++ -O3 -fopenmp -std=c++17 -o verify_odd  verify_odd.cpp
g++ -O3 -fopenmp -std=c++17 -o verify_even verify_even.cpp

./verify_odd  --selftest && ./verify_odd  --test
./verify_even --selftest && ./verify_even --test

export BUILD_FLAGS="-O3 -fopenmp -std=c++17"
bash run_manifest.sh                                   # provenance, before the run
./verify_odd  2>&1 | tee results/odd_run.log  && gzip -f results/odd_run.log
./verify_even 2>&1 | tee results/even_run.log && gzip -f results/even_run.log
bash run_manifest.sh                                   # provenance, after the run
python3 verify_results.py --sample 200
```

* **Threads.** `verify_odd` uses OpenMP's thread count (`OMP_NUM_THREADS` if
  set, otherwise every logical CPU); `verify_even` runs on a fixed six threads
  (`NUM_THREADS`). The results do not depend on either.
* **Resuming.** Both are resumable with `--resume` and write their results once
  per batch of chunks, so a killed run loses at most the current batch. Without
  `--resume`, `verify_odd` starts its results file afresh, while `verify_even`
  appends to an existing one, so move the old even files aside before a fresh
  even run.
* **The fallback count.** Each row of `odd_results.csv` ends with
  `fallback_count`: how many odd `n` in that chunk the fast path did not clear,
  so the exhaustive search had to. `verify_results.py` reports the total.
* **Counts in one place.** The witness and split counts appear once each,
  `NUM_WITNESSES` in `verify_odd.cpp` and `NUM_SPLITS` in `verify_even.cpp`;
  `verify_results.py` carries copies that must be kept in step with them.
* **Provenance.** `run_manifest.sh` appends a snapshot to
  `results/run_manifest.txt`: commit, compiler and flags, and checksums of the
  sources, the executables and the results. The executable checksum depends on
  the exact compiler flags, including their order, so record them in
  `BUILD_FLAGS` as built.

## Citing this repository

Archived on Zenodo. The DOI below is the *concept* DOI: it always resolves to
the latest release, so it stays correct as the repository is updated.

> Hoare, W. *Reproducibility code and data for "On the sum of a prime and a
> square-free number coprime to integers with at most five prime factors"*.
> Zenodo. https://doi.org/10.5281/zenodo.21813035

```bibtex
@misc{hoare2026code,
  author       = {Hoare, Wilkie},
  title        = {Reproducibility code and data for {On the sum of a prime and
                  a square-free number coprime to integers with at most five
                  prime factors}},
  year         = {2026},
  howpublished = {Zenodo archived release},
  doi          = {10.5281/zenodo.21813035}
}
```

## Results of record

The default (Miller-Rabin) builds of `verify_odd` and `verify_even`, together
with `Threshold.py --certify` for `108 <= n <= 1e4`.
