# ETBF Longline Data Simulator — Documentation

Documents two files:

- `simdatfun.r` — the core simulation function `simdatfun()`.
- `ETBF_LL_simulator.R` — a driver script that sets parameters, calls the function, and writes output.

Author: Simon Hoyle, Jan 2016, building on Henning Winker's code (per header).
Purpose (inferred): generate synthetic catch-and-effort records for an
Eastern Tuna and Billfish Fishery (ETBF) longline fishery, with a known
"true" biomass trajectory, to test CPUE standardisation / spatial methods.

> Note on confidence. I could not run the code (no R available here). Claims
> about *what the code does* are from reading + R language semantics. Claims
> about *intent* are marked **(inferred)**.

---

## 1. `simdatfun.r` — `simdatfun()`

### 1.1 Conceptual model

Catch per set is built from the relationship (stated in the file's comments):

```
catch = effort × q × density
      = effort × q_tactic(sp) × q_vessel(vessel) × q_season(month) × B(sp,year) × relden(sp,location)
```

In the code `effort` is a constant `1` per set, so each set carries one unit
of effort and the per-set expected catch reduces to
`q × B × relative_density`. Observed catch is then drawn from a Tweedie or
zero-inflated negative binomial distribution.

### 1.2 Sampling hierarchy

```
vessel ─┬─ fixed base tactic (cycles 1:ntactics); lower-half tactics switch up at year floor(nyr/2)
        └─ vessel-year ─┬─ ntrip trips  (sample of tripmode ± 2)
                        └─ trip ─┬─ nsets sets   (sample 5–13)
                                 └─ set ─ assigned a day, month, lat, lon, tactic
```

Each set is the unit of observation in the returned data.

### 1.3 Parameters

| Arg | Meaning | Notes |
|---|---|---|
| `nclus` | (declared) number of clusters | **Unused** (kept for call compatibility) |
| `nspec` | number of species | General; `nspec==10` reproduces the original exactly |
| `spnames` | species short names | Length `nspec` |
| `nyr` | number of years | |
| `p` | Tweedie power parameter | passed to `rTweedie` |
| `disp` | dispersion | Tweedie `phi`, or NB `size` |
| `qq` | tactic × species catchability matrix | `ntactics × nspec` |
| `nvess` | number of vessels | |
| `nlat`, `nlon` | spatial grid dimensions | `nlon > 1` for `ew_trend`/`ew_rndom` |
| `ntactics` | number of fishing tactics | honoured (vessels cycle `1:ntactics`) |
| `tripmode` | central value for trips/vessel/year | actual = `tripmode ± 2` |
| `spdist` | spatial distribution model | `"rndom"`, `"ew_trend"`, `"ew_rndom"` |
| `rtype` | observation distribution | `"Tweedie"` or `"zinb"` |

### 1.4 Processing pipeline

1. **Biomass.** Per-species growth rate `r` drawn `U(0,0.1)` with random sign.
   Start biomass `B1` drawn lognormally (abundant vs rare split; `nspec==10`
   keeps the original 3/1/1/5 pattern). Trajectory
   `B[year, sp] = B1 × exp((year−1)·r)` — deterministic exponential.
2. **Seasonal catchability** `qmon`: a 12-month table, hardcoded for the first
   5 species; any further species are set to 1 (no seasonality).
3. **Vessels** `vess`: tactics assigned `rep(1:ntactics, length.out=nvess)`;
   the first `floor(ntactics/2)` tactics switch to a higher tactic at year
   `floor(nyr/2)`; vessel catchability `vq` lognormal, log-mean 0 (median 1).
4. **Vessel-years** `vessyr`: trips per year sampled; tactic copied, with the
   mid-series tactic change applied.
5. **Trips / sets** built per vessel-year by a vectorised expansion: `rep()`
   repeats each parent's fields and `sequence()` numbers the children within
   each parent (replacing the earlier per-row `data.frame()` + `rbind`).
6. **Calendar placement.** Sets per vessel-year are summed **by key and joined
   back** (order-safe); each year has `360 − nsets` "spare" days; trip start
   days are sampled and sets run consecutively, giving `day` and
   `mon = min(floor(day/30)+1, 12)` per set.
7. **Tactic per set.** Trip tactic is **joined** from the vessel-year on the
   key (`match`), then copied to each set.
8. **Spatial density** `relden`: relative density per species per cell, from
   `do_rndom`, `do_ewtrend` (W→E gradient along longitude, constant across
   latitude), or a logit blend (`ew_rndom`); then **normalised to mean 1 per
   species** across cells.
9. **Spatial effort.** First set of each trip placed at random; later sets
   follow a constrained random walk (reflecting at the grid edges). Computed on
   integer vectors with step increments drawn up front, then a sequential pass
   applies them (the per-step boundary clamp prevents a pure `cumsum`).
10. **Realised catchability** `qsp = q_tactic × q_vessel × q_month`.
11. **Expected catch** `catch_M = B × relden × qsp` (effort = 1).
12. **Observed catch** `catchobs` drawn via `rTweedie` or `rzinegbin`, rounded.

### 1.5 Output

A list with:

- `sets` — one row per set: identifiers (`vess`, `yr`, `trip`, `set`, keys),
  `day`, `mon`, `tactic`, `lat`, `lon`, `loc`, then expected catch columns
  (`<sp>_M`) and observed catch columns (`<sp>_R`).
- `B` — the `nyr × nspec` true biomass matrix.

### 1.6 Dependencies

| Function | Package | Status |
|---|---|---|
| `rTweedie` | `mgcv` | loaded by `simdatfun.R` |
| `rzinegbin` | `VGAM` | loaded by `simdatfun.R` (needed for `rtype="zinb"`) |
| `logit`, `inv.logit` | `boot` (or `gtools`) | loaded by `simdatfun.R` (needed for `spdist="ew_rndom"`) |
| `mat.or.vec` | base | n/a |

`simdatfun.R` now `requireNamespace`-checks and loads all three at source time.

### 1.7 Performance & reproducibility

Profiling (profvis/Rprof) of a representative run showed almost all time in the
per-set random walk (~92% of in-function samples, with data-frame row indexing
and ~12% GC), plus a smaller cost in the trip/set construction. Two changes
address this:

- **Trip/set expansion** built once with `rep()` + `sequence()` instead of a
  per-row `data.frame()` + `do.call(rbind, …)` (avoids repeated allocation and
  `rbind` name-matching).
- **Random walk** run on integer vectors with all step increments drawn up
  front, instead of per-row data-frame assignment and a `sample()` call per
  set.

Reproducibility caveat: both changes alter the order in which random numbers
are drawn, so a `set.seed()` run produces a **different** data set than the
pre-optimisation code. The statistical properties are unchanged; only exact
byte-for-byte continuity with previously generated data is lost.

---

## 2. `ETBF_LL_simulator.R` — driver (revised, portable)

### 2.1 What it does

1. **Configuration block.** `base_dir` (a single output root, default
   `file.path(getwd(),"sim_output")`), `fun_path` (location of the function
   file), and `set.seed(42)` for reproducibility. All output is written under
   `base_dir` via `file.path()` — no absolute or Windows-specific paths.
2. `source()`s `simdatfun.R`, which defines `simdatfun()` and loads its
   dependencies (`mgcv`, `VGAM`, `boot`).
3. Defines run constants (`nspec`, `nyr`, `ntactics`, `p`, `disp`, `spnames`)
   and the `qq` catchability matrix.
4. Runs scenarios:
   - `sets1`/`sets2`/`sets3` — single illustrative runs (varying grid size,
     dispersion, `spdist`, `rtype`); writes `smallsim.csv` / `smallsimB.csv`.
   - `run_replicates()` — a helper that writes 100 replicate datasets for a
     given `spdist`, called for `ew_trend` (`SDHsim3/`) and `ew_rndom`
     (`SDHsim4c/`), producing `simData<i>.csv` / `simB<i>.csv`.
   - A `memtests/` loop (`sim2000`…`sim16000`) of increasing size for
     memory/scaling tests.
5. **Diagnostic** `plot_cpue_vs_B()`: overlays normalised year-mean catch
   (expected `_M` or observed `_R`) on the normalised true biomass for each
   species, and writes a PNG. Cross-platform (`png()`/`dev.off()`), and the
   panel grid adapts to `nspec`. Called on `sets1` for both `_M` and `_R`.

### 2.2 The `qq` matrix

Rows = tactics (1 deep-LL BET, 2 shallow-LL SWO, 3 YFT+ALB, 4 other);
columns = species. Must have `ntactics` rows and `nspec` columns. The unused
alternative value block from the original has been dropped.

### 2.3 Notes on the driver

- Run constants are defined once and passed explicitly into each call. To run
  elsewhere, set `base_dir` and `fun_path`.
- `nclus` is still accepted by `simdatfun()` but unused; it is passed for
  backward compatibility.
- The earlier inert `nlat=5; nlon=6` assignments and the broken list-vs-data
  frame "First check" block have been removed/replaced by the diagnostic above.

---

## 3. Data model summary (quick reference)

| Object | Grain | Key fields |
|---|---|---|
| `B` | year × species | true biomass |
| `vess` | vessel | tactic, tactic-change, `vq` |
| `vessyr` | vessel × year | `ntrip`, `tac`, `nsets` |
| `trip` | trip | `nsets`, `tripst`, `tripnd`, `tactic` |
| `relden` | grid cell × species | relative density |
| `sets` | **set** | day, mon, lat, lon, tactic, `_M`, `_R` |
