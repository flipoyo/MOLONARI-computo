# Handoff: PyHeatMy in the MOLONARI-computo tree

*13/11/2026 · MOLONARI team, PyHeatMy work package (Yunus) · state described: the MOLONARI-computo tree 1.0.1 (pyheatmy 0.12.2, `bayesian_inversion` 0.4.3), branch `CGS-MOLO` of its five repositories*

This document is for the team that takes over PyHeatMy. It brings together, in a single text, the observation report (B1), the corrections report (B2) and the report on the split into two packages. It ends with where to start.

---

## 0. In brief

PyHeatMy estimates the properties of the river bed and the water flux between the river and the aquifer. It starts from temperatures measured at several depths and from the head difference. It has two roles:

- the **direct model**, which simulates the temperatures for given parameters;
- the **Bayesian inversion** (MCMC), which goes back from the measurements to the parameters.

This year, the work went through three stages.

1. **A verified base.** Before restructuring, the code was examined. The review produced **36 findings**, all handled with the same method (§1). Today:
   - 28 are fixed and verified (§2.2);
   - 2 are partly fixed (B-11, B-29);
   - 6 are documented with their analysis and await a decision (B-09, B-18, B-24, B-27, B-30, B-31; §2.3).
2. **Two independent packages.** The direct model (`pyheatmy`) and the inversion (`bayesian_inversion`) do not import each other. The only link is a function, `Column.simulate`, which the user passes to `run_mcmc` (§3).
3. **What you need to start.** There are three tutorials, tests in each package run at every push, and integration tests in the tree (§4).

The project manager then reorganised the tree (9–10/11): new repositories, one pixi environment per repository, and a fifth repository for processing raw data (§3.2). The code of the two packages did not change, apart from their version number and the removal of `lora_emitter.py`, so the measurements in this document still hold. The computo [README](README.md), for users, and the [developer guide](docs/README_DEV.md) describe the installation and the organisation of the tree; this document adds the scientific state of the code.

| Path in the tree | GitHub repository | Content | Version | Commit (13/11) |
|---|---|---|---|---|
| `.` | `flipoyo/MOLONARI-computo` | root: `integration_tests/`, `data/`, `research/`, `docs/`, pixi environment, tree description (`molo-computo.cgs`) | 1.0.1 | `211b2db` |
| `pyheatmy` | `flipoyo/MOLO-pyheatmy` | `pyheatmy`, direct model | 0.12.2 | `5d649bb` |
| `bayesian-inversion` | `flipoyo/MOLO-bayesian` | `bayesian_inversion`, MCMC (DREAM, random walk Metropolis) | 0.4.3 | `2ec31fe` |
| `data-pre-proc` | `flipoyo/MOLO-data-process` | `molonari-pipeline` (`raw_to_processed`), processing of raw measurements | 0.2.0 | `80fa923` |
| `TUTO` | `flipoyo/MOLO-tuto-cmpt` | tutorials 01 to 03 | 1.0.1 | `cf8d284` |

**Tags.** `MOLONARI-computo-1` (release of the tree on 30/10/2026) now exists only in MOLONARI-computo: the four other repositories were created later, without tags. The state described here is that of the commits in the table above.

**How to read this document:**

- to start quickly: §0, then §4;
- to pick up a finding: its row in §1.3 (what was seen), then §2.2 or §2.3 (what was done, or what remains to be decided);
- for the code: §3.

**Sources.** The detailed working reports are not in the repositories. This document keeps everything from them that is needed to continue. Commit identifiers refer to three histories:

- **before 7/10/2026, the origin of the findings:** the history of MOLONARI1D, folder `pyheatmy/pyheatmy/`;
- **fixes up to version 0.9.7:** the history of MOLONARI-computo, in the `pyheatmy/` folder, since removed but present at the tag `MOLONARI-computo-1`;
- **pyheatmy 0.10 to 0.12.1 and `bayesian_inversion` 0.1 to 0.4.2:** the repositories `YunusEmre4418/MOLONARI-pyheatmy`, `YunusEmre4418/MOLONARI-bayesian-inversion` and `YunusEmre4418/MOLONARI-tuto`, which are no longer updated. The `MOLO-*` repositories of the current tree start from a first commit, without this history.

---

## 1. Findings

### 1.1 One method for every finding

PyHeatMy has been developed since 2023, and each cohort builds on the work of the previous one. Behind every line, there is someone's choice. Each finding was handled in four steps:

1. **Compare expected and observed behaviour.** The expected behaviour comes from physics, mathematics or the documentation. The observed behaviour comes from running the code, whenever possible.
2. **Read the history of the code** to find out whether it was a deliberate choice. Sometimes a part that was originally right had changed behaviour during a merge. A deliberate choice is not "fixed": it is assessed and documented. Single precision (B-26) and constant viscosity (B-18) are examples.
3. **Fix** following the formulation in the literature (sources in Appendix A).
4. **Compare the outputs before and after**, to check that only what had to change has changed (§2.1).

**Levels of evidence** used in the table of §1.3:

| Level | Meaning |
|---|---|
| run | the code was run and the behaviour observed |
| test | a test written for the purpose fails on the old code |
| calculation | a value computed by hand differs from the code's |
| code | reading only, or search (`grep`) |

**Reading the history.** Searches are done in MOLONARI1D, on the main line:

- `git log ed4fee8 -S'<text>' -- '*<file>.py'` gives the commits that added or removed a text;
- `git blame ed4fee8 -- pyheatmy/pyheatmy/<file>.py` gives the origin of each line of the examined state.

Avoid `--all`, which also scans branches that were never merged. If `blame` points to `df5a87d` (October 2024), that is a general reformatting: the real origin is older and is found with `-S`.

### 1.2 The code at the start

The examined state is the release `MOLONARI1D-0.1.0`, that is pyheatmy 0.4.0 in MOLONARI-computo (`a464ba2`). On Point034, the direct model was reproducible to the bit. The MCMC was reproducible only if both random generators were seeded (B-11).

- **Inputs.** `Column` receives:
  - the head difference and the river temperature;
  - the temperatures of the shaft.

  The last sensor of the shaft is not an observation: it is the bottom boundary condition. Temperatures are in kelvins. The conversion of a voltage into head (sensor calibration) is done in the notebook, not in the library.
- **Direct model** (`compute_solve_transi`):
  1. a regular mesh of the column;
  2. the head H, solved first;
  3. the Darcy flux, derived from H;
  4. the temperature, with this flux.

  The solver is tridiagonal and compiled with `numba`. One branch handles a single layer, another several layers.
- **Parameters.** A `Layer` holds five parameters: IntrinK, n, λs, ρs·cs and q_s. Each can have a `Prior`, whose scale (log, symlog or linear) is chosen automatically from its range. The MCMC moves in this transformed space.
- **Inversion** (`Column.compute_mcmc`):
  - with a single chain, the algorithm is a random walk Metropolis (RWM);
  - with several chains, it is DREAM.

  The energy (minus the log-likelihood) is computed at the sensor cells. In `"no sigma"` mode, σ² = 1 is fixed. In `"estimesigma"` mode, σ² is estimated.
- **Randomness.** The code draws numbers with two generators: `random` (priors, RWM acceptance) and `np.random` (DREAM, synthetic noise).

The current state of the code is described in §3.3.

### 1.3 The 36 findings

The findings are grouped by area. The "Origin" column summarises step 2 of the method; "—" means it was not investigated. The "Status" column refers to §2.

#### A. Direct model: physics and numerics

| | Finding | Evidence | Origin | Status |
|---|---|---|---|---|
| B-01 | **Two conventions for the intrinsic permeability.** The solver uses IntrinK in m². The interface rows between layers and the initial head profile (multilayer) used it in −log10. With `IntrinK = 11`, half the cells oscillate at ±8.7·10⁷ °C, while the sensor cells stay plausible (RMSE 0.53 to 1.08 K): **neither the RMSE nor the energy sees it**. With a value in m², the interface cell is computed with k ≈ 1 m² | run | Switch from −log10 to m² (11/2025, `9287f72`): these two lines were not converted. The conversion of k into K was written in six places | fixed |
| B-12 | `get_dt` divides the total duration by the number of points instead of the number of steps: 899.3 s instead of 900 s on Point034. Only the time axis of the plots is affected | run | — | fixed |
| B-16 | With a variable time step, the interface rows always used the first step | test | Block moved into `correct_numerical_schema` (10/2024, `1468224`); the 2023 version was right | fixed |
| B-17 | At an interface between two layers, the second cell used the storage Ss of the first | test | Since it was first written (11/2023, `35756c0`) | fixed |
| B-18 | **Updating the viscosity μ(T) has no effect.** `compute_Mu` returns a constant, and K is not recomputed. Between 5 and 25 °C, μ varies by a factor 1.65, and so does K; on Point034 (13 to 18 °C), the difference is about 10 %. With constant μ, the MCMC absorbs this difference into k | code, calculation | The μ(T) law was commented out on 15/11/2024 (`fca431a`) after an *overflow* warning. That warning is reproduced with temperatures in °C, not in kelvins. The time coupling T → μ → K → H was never written | to decide |
| B-19 | **Two meshes.** `getListParameters` built its own mesh from the `zLow` of the last layer, not from the bottom of the column. As soon as the two differ (for instance `offset ≠ 0`), the interface shifts: 0.26 to 0.29 K at the sensors on Point034, and a factor of about 10 on the flux in the shifted cells | run | Two meshes since the first commit (06/2023) | fixed |
| B-26 | **Single-precision computation.** In pure conduction, the gap to the analytical solution reaches 0.17 K after 10 days and grows with time. It is proportional to the mean temperature (−0.170 K at 285 K, −0.015 K at 28.5 K): it is rounding, not physics | run | Design choice since the first version (11/2023, `0a585dc`), probably for memory | fixed |
| B-27 | **The lateral source term q_s does not contain the temperature of the incoming water.** The matrix adds a term proportional to T, in kelvins. A column at 285 K everywhere, without flow, warms by 0.13 K in 4 days with q_s = 10⁻⁸ s⁻¹, and by 1.3 K with 10⁻⁷. The diagnostic `source_heat_flux` computes the same term in °C | run | Reported by issue #2 of MOLONARI-computo (10/2026) | to decide |

#### B. Inversion (MCMC)

| | Finding | Evidence | Origin | Status |
|---|---|---|---|---|
| B-02 | **Probability of γ = 1 in DREAM.** In `np.random.choice([γ, 1], 1, [0.8, 0.2])`, the third argument is `replace`, not `p`: γ = 1 is drawn in 51 % of cases instead of 20 % | run | Since the first version (11/2023) | fixed |
| B-03 | **DREAM states shared the same memory** (a numpy view): all states of a chain held its last position. Energy and parameters of a state matched for only 42 states out of 296 | run | Manual merge of 11/2024 (`4bf008b`); the code was right on both sides of the merge | fixed |
| B-04 | **Initial DREAM energy** computed before the direct model, against an array of zeros (about 1.6·10⁸): the first proposals were accepted unconditionally | run | Same merge (`4bf008b`) | fixed |
| B-05 | **RWM states stored in physical space** and read back as if they were in MCMC space. The best state of the initial sweep was lost, and the parameter history gave k ≈ 1 m² | run | Introduction of the two spaces (11/2025): two lines not adapted | fixed |
| B-06 | **Acceptance rate divided twice.** RWM displayed 0.0033 for a real rate of 0.33, DREAM a rate 35 times too small | run | Merge of 11/2024 (`4bf008b`) | fixed |
| B-08 | **Quantiles included rejected proposals** (with three twin cases). The 95 % temperature quantile was shifted by 1.4 K, the flux median by 44 % | test | — | fixed |
| B-15 | The σ² proposal in DREAM always started from the first state of the chain | code, test | Merge of 11/2024 (`4bf008b`) | fixed |
| B-21 | In DREAM, the end-of-burn-in position was not in the history. A chain that rejected its first proposal repeated its state drawn **before** the burn-in | test | Merge of 11/2024 (`4bf008b`) | fixed |
| B-22 | **`estimesigma` mode.** σ² was drawn in log10 but used as a variance (negative σ², no proposal accepted), and the Jeffreys prior was counted twice | run | Automatic log scale of the priors (11/2025, `9287f72`): the σ² lines were not adapted | fixed |
| B-23 | **Likelihood.** The coefficient of log σ² was n instead of n/2. A hand-computed example (σ² = 2, 300 measurements, 1 K offset) expects 178.97; the code gave 282.94. The difference is invisible at σ² = 1 | calculation | Merge of two energy functions (12/2024, `3f81cb3`); the original function did divide by 2 | fixed |
| B-24 | **The scale chosen for a prior depends on the sign of its range.** `(1e-9, 1e-6)` is on a log scale, but `(-1e-6, -1e-9)` on a linear scale. No current prior is affected | run | — | to decide |
| B-29 | **DREAM burn-in.** The Gelman–Rubin criterion used all samples since the initial population. A fixed parameter received R = 2 and prevented convergence. The threshold depended on the unit, and no warning was given at the limit of 25 iterations. On Point034, the old criterion "converged" after 47 iterations while R̂ was still up to 1.41 | run | Criterion present since the first version of DREAM (10/2023); blocking by fixed parameters since those parameters were introduced (10/2025) | partly fixed |
| B-30 | **DREAM does not handle outlier chains.** One chain out of ten stayed at k ≈ 3·10⁻¹⁵ m² for 300 iterations and entered the posterior | code, observation | Step never written | to decide |
| B-32 | `compute_mcmc` stopped (`IndexError`) when the number of iterations was not a multiple of the quantile subsampling step | run | Storage condition changed without the array size (12/2024, `3f81cb3`) | fixed |
| B-34 | With fewer than `2·delta + 1` chains, DREAM stopped on an unclear numpy error | run | — | fixed |
| B-35 | RWM stopped as soon as the quantiles were subsampled in space or time | run | — | fixed |

#### C. Data and inputs

| | Finding | Evidence | Origin | Status |
|---|---|---|---|---|
| B-09 | `sigma_meas_T` and `sigma_meas_P` are accepted but never used. In `"no sigma"` mode, the uncertainty therefore does not depend on the sensor error | code | — | to decide |
| B-10 | **The periodic analytical solution (`Analy_Sol`) no longer worked.** The code read an attribute that no longer existed and used k instead of K. The velocity was divided by one cell and had the opposite sign. The signals were read in the old format, and the RMSE was computed incorrectly | run | Switch from −log10 K to −log10 k (10/2024, `b04d191`) | fixed |
| B-11 | **Reproducibility.** The demo only set `np.random.seed`, while the code uses two generators: two runs of the same notebook gave two different MCMCs | run | — | partly fixed |
| B-13 | **Input checks.** The code raised a `NameError` and checked neither dates nor depths. A single missing temperature (NaN) at a sensor made every energy NaN, and the RMSE divided by the total number of times | run | Checks since 06/2023; norm instead of `nansum` in the energy (12/2024, `3f81cb3`) | fixed |
| B-28 | **Noise of the synthetic data added in place**: the noise-free series were lost, and the bottom sensor received the noise twice (standard deviation √2·σ) | run | `+` became `+=` (10/2024, `cfdf9c4`) | fixed |
| B-31 | **`demoInversion.ipynb`.** Temperatures were converted to kelvins twice (561 to 572 K passed to the inversion). The P508 sensor calibration was applied to a head already in metres. The "true" values differed from those of the generator, with only 20 iterations. The same calibration error is in `integration_tests/conftest.py` (§2.3) | run | Conversions since 11/2024 (`2c7550c`); reference values changed in 10–11/2025 | to decide |

#### D. Interface and robustness of the code

| | Finding | Evidence | Origin | Status |
|---|---|---|---|---|
| B-07 | Three accessors (`sample_param`, `get_best_param`, `get_all_energy`) raised `AttributeError` | run | — | fixed |
| B-14 | Dead and duplicated code: unused `Gelman_Rubin.py`, dead attributes of `Prior`, ignored argument of `Layer.perturb`, commented-out blocks. The inventory has 20 items, each with its origin | code | — | fixed |
| B-20 | `Column(all_layers=…)` was ignored, although the API documentation recommended this usage | run | Argument added in 11/2024 (`5fd7bc6`), never wired | fixed |
| B-25 | `set_priors_from_dict` read the stored values back through the **new** priors: IntrinK = 10⁻¹² m² became 1.0 m² in the call order of the demo | run | Separation of the two spaces (11/2025, `9287f72`): the comment shows the intent, but the order of reading and assignment is reversed | fixed |
| B-33 | With DREAM, two accessors required a prior call to `compute_solve_transi` | run | DREAM worked on copies of the column | fixed |
| B-36 | **A solved column could no longer be serialised** (`pickle`), and its copy solved the original column. This prevented parallel execution | run | `checker` set on the instance a method bound to it | fixed |

**Observations that are not errors:**

- With four free parameters (k, n, λs, ρs·cs), DREAM accepts 0.8 to 3.2 % of the proposals, far from the theoretical optimum of about 23 % (Roberts, Gelman and Gilks, 1997). The temperatures constrain combinations of n, λs and ρs·cs rather than each of them.
- The demo prior of λs goes up to 10 W/m/K, while literature values do not exceed 4.5 W/m/K (Anderson, 2005).
- In the synthetic data, the generating column and the inversion column receive the same noisy river temperature. The inversion therefore never sees a measurement error at the boundary.
- **Specific storage.** Since the first commit (06/2023), the code uses `Ss = n / column height`, that is 0.25 m⁻¹ in the demo. The definition through compressibility, Ss = ρg(α + nβ), gives 10⁻⁶ to 10⁻³ m⁻¹ for sands and gravels. This order of magnitude was given from memory, after Freeze and Cherry (1979), and remains to be checked in the source. The response time of the head would then go from about 4000 s to less than 20 s (§2.3).

---

## 2. Fixes

### 2.1 Verifying a fix: two complementary checks

Every fix followed the same sequence:

1. analysis of the code and its history (§1.1);
2. a "before" test, which must fail on the old code;
3. the change;
4. a review of the diff;
5. two checks, described below.

One fix fits in one commit; a few fixes that could not be separated were grouped.

**(a) References ("golden"): a change detector.** A reference records what the code produces for a given input.

- **Starting point.** The outputs of `MOLONARI1D-0.1.0` on five scenarios, identical to the bit over two runs:

  | Scenario | Content |
  |---|---|
  | S1 | direct model, one layer, IntrinK = 10⁻¹² m² |
  | S2 | IntrinK variants and two-layer column |
  | S3 | RWM, one chain, 20 + 100 iterations |
  | S4 | DREAM, ten chains, 20 + 35 iterations |
  | S5 | one month of synthetic data, σ_T = 0.1 K |

  The demonstration notebook is added to these.
- **Chain.** After each fix, all scenarios are run again and compared **with the previous step**. Each array is compared to the bit: type, shape and values, with NaN equal to NaN. Each difference can thus be traced to a single fix.
- **Predicted signature.** Before each fix, we write which outputs must change and which must not. The check consists in finding exactly this signature. Any unexpected difference stops the work until it is explained.

  *Example:* for B-01, we expected the acceptance decisions of the MCMC to change. They stayed identical, because the effect in a single layer is too small for that. The prediction was corrected in writing, with this explanation.
- **Limit.** A reference records what the code produces, not what is right:
  - those of `MOLONARI1D-0.1.0` contained the wrong energies of B-04;
  - B-23 is invisible at σ² = 1, hence in every reference.

  Check (b) is what establishes that a fix is right.

**(b) Correctness tests: the expected value does not come from the code.** It comes from a hand calculation, an analytical solution, a known posterior distribution or a conservation law.

- Each test fails on the old code (marked `xfail(strict=True)` until the fix is made), then passes after it.
- When one error hid others, we used **mutations**: each part of the fix is undone, one at a time, in a fixed copy, and the tests must detect each undoing. Number of mutations detected:

  | Findings | Mutations detected |
  |---|---|
  | B-10 | 7 of 7 |
  | B-22 and B-23 | 4 of 4 |
  | B-29 | 8 of 8 |
  | B-13 | 13 of 13 |

The main correctness tests are listed below. They are now in the repositories.

| Expected value | File |
|---|---|
| Steady flux through two layers in series: ΔH / (L₁/K₁ + L₂/K₂), the same in every cell | `pyheatmy/tests/integration/test_steady_state.py` |
| Periodic analytical solution (Stallman): no flow, downward flow, upward flow, with a phase; tolerance 0.01 K | `pyheatmy/tests/integration/test_analytical_periodic.py` |
| Energy: half the sum of squares at σ² = 1; minus the Gaussian log-likelihood; NaN left out | `bayesian-inversion/tests/test_energy.py` |
| Samplers on a 2D Gaussian: mean and variance within the Monte Carlo error | `bayesian-inversion/tests/test_samplers.py` |
| Gelman–Rubin R̂ compared with hand-computed values | `bayesian-inversion/tests/test_gelman_rubin.py` |
| Posterior mean of σ² with fixed parameters: inverse gamma distribution (shape n/2, scale S/2), for DREAM and RWM | `integration_tests/test_sigma2_estimation.py` |
| MCMC with fixed parameters: quantiles equal to the direct model | `integration_tests/test_mcmc_subsampling.py` |
| Every recorded state: its energy is that of its parameters, replayed through the direct model | `integration_tests/test_dream_records.py`, `test_rwm_records.py` |

Today, `pixi run test-all` runs 68 tests for `pyheatmy`, 66 for `bayesian_inversion`, 120 for `molonari-pipeline` (one skipped) and 31 integration tests. The tests of each package also run on GitHub Actions, in pixi, with Python 3.10 and 3.13. The reference scripts and the references themselves (49 arrays, plus the notebook) stayed outside the repositories. For a new change, protocol (a) has to be set up again: record the outputs before, then compare to the bit after (§4.2).

### 2.2 The 28 fixed findings

Versions are those of `pyheatmy`, unless `bayesian_inversion` ("bi") is stated. The "Test" column gives the file that protects the fix against a regression today. Each file path is relative to `pyheatmy/tests`, `bayesian-inversion/tests` or `integration_tests/`.

#### A. Direct model

| | Fix | Version | Measured effect | Test |
|---|---|---|---|---|
| B-01 | A single function converts k into K (`calc_K(k, viscosity)`) and is used everywhere. `Layer` refuses an IntrinK outside 0 < k ≤ 10⁻⁴ m² | 0.8.5, 0.9.0 | **One layer:** ≤ 4·10⁻⁴ K. **Two layers** (10⁻¹² / 10⁻¹³ m²): up to 1.4 K; in steady state, the ratio of maximum to median flux goes from 5.57 to 1.015, so the flux becomes the same throughout the column, as mass conservation requires. **Old notation:** it raises an explicit `ValueError` | `test_steady_state.py`, `test_layer_permeability.py`, `test_intrink_prior.py` |
| B-12 | The mean step is the mean of the actual intervals | 0.5.0 | 899.31 s → 900 s; no other output changes | no dedicated test |
| B-16 | The current time step is passed as an argument to the interface rows | 0.8.4 | References unchanged (constant step in Point034) | `test_head_interface.py`: two steps at once = two one-step computations |
| B-17 | Each row uses the Ss of its own cell | 0.8.3 | References unchanged (same porosity in all layers) | `test_head_interface.py` |
| B-19 | The parameters are placed on the mesh of the column | 0.8.0 | References unchanged (`offset = 0` everywhere); inconsistent geometries now give the same result to the bit | dedicated reference scenarios; no test in the repositories |
| B-26 | Double-precision computation; only the samples stored for the quantiles stay in single precision (memory) | 0.9.5 | **Gap to the analytical solution:** 0.012–0.025 K → 0.0035–0.0047 K. **Computing time:** 7 to 8 % shorter. **References:** all change, by ≤ 0.0024 K for the direct model and ≤ 0.022 K for the MCMC quantiles | `test_analytical_periodic.py` (tolerance lowered from 0.05 to 0.01 K) |

#### B. Inversion

| | Fix | Version | Measured effect | Test |
|---|---|---|---|---|
| B-02 | Argument `p=`. The probability and the value of the jump (γ = 1) become constants of `config.py` | 0.8.2 | γ = 1 in 49.9 % → 19.8 % of the proposals (target 20 %) | `test_dream_proposal.py` |
| B-03 | Each state keeps its own copy of the parameters | 0.8.0 | Energy–parameter agreement: 42/296 → 296/296 states | `test_dream_records.py` |
| B-04 | Initial energy computed after the direct model | 0.8.2 | Agreement over all states: 296 → 360/360 | `test_dream_records.py` |
| B-05 | States stored in MCMC space | 0.8.1 | Contract respected: 0 → 101/101 states | `test_rwm_records.py` |
| B-06 | Acceptance rate computed once, in the branch of each algorithm | 0.8.0 | RWM: 0.0033 → 0.33; DREAM: ×35 | `test_dream_records.py`, `test_rwm_records.py` |
| B-08 | Outputs recorded only in the acceptance branch, at the four places concerned | 0.8.0 | Quantiles now come from accepted states only | `test_dream_records.py`, `test_rwm_records.py` |
| B-15 | Each chain keeps its current σ² | 0.8.2 | σ² contract: 126 violations → 0 | `test_dream_records.py` |
| B-21 | States are appended at the end of the burn-in | 0.8.2 | S4 quantiles identical to the replay: 0 → 6/6 | `test_dream_records.py` |
| B-22 | Chains move in log10 σ², and states keep the physical σ². A uniform distribution in log σ² is exactly the Jeffreys prior 1/σ²: the extra 1/σ² density is removed. RWM is treated the same way | 0.9.3 | References unchanged (all in `"no sigma"`). Posterior mean of σ²: gap to the inverse gamma distribution of 0.4 % (DREAM) and 1.1 % (RWM); before, σ² was negative (−1.29 and −0.47) | `test_sigma2_estimation.py` |
| B-23 | Coefficient n/2 | 0.9.3 | Without the `/ 2`, the estimate of σ² is exactly halved (0.00495 instead of 0.00993) | `test_energy.py`, `test_sigma2_estimation.py` |
| B-32 | Size of the quantile array rounded up | 0.9.4 | References unchanged | `test_mcmc_subsampling.py` |
| B-34 | `DREAM` gives the minimum number of chains as soon as it is built, before any model run | bi 0.3.2 | References unchanged | `test_samplers.py` |
| B-35 | Arrays allocated at the subsampled size | bi 0.3.3 | References unchanged | `test_mcmc_subsampling.py` |

#### C. Data and inputs

| | Fix | Version | Measured effect | Test |
|---|---|---|---|---|
| B-10 | Parameters read through the current API; K instead of k; velocity +α·ΔH/L; signals with their phase; RMSE | 0.9.1 | Gap between the direct model and the analytical solution: 0.012–0.025 K for an amplitude of 2 K, then 0.0035–0.0047 K after B-26. With the velocity of opposite sign, the gap would be 0.19 K | `test_analytical_periodic.py` |
| B-13 | NaN at the sensors are left out of the initial profile, the energy (n = number of measured values) and the RMSE, with a single warning. Structural errors raise an explicit `ValueError` before any computation: number of times, identical and strictly increasing dates, number of depths, NaN at the boundaries | 0.9.7 | References unchanged. Without NaN, energy and RMSE are identical to the bit; with NaN, DREAM and RWM accept proposals again | `test_column_inputs.py`, `test_missing_measurements.py`, `test_energy.py`, `test_missing_measurements_mcmc.py` |
| B-28 | The noise is added to a copy | 0.9.2 | Only S5 changes (bottom sensor: up to 0.37 K); the noisy series are identical to the bit | `test_synthetic_noise.py` |

#### D. Interface and robustness

| | Fix | Version | Measured effect | Test |
|---|---|---|---|---|
| B-07 | Accessors fixed, then replaced by the result object `MCMCResult` | 0.5.0, bi 0.4.0 | No computed output changes | `test_mcmc.py` |
| B-14 | **Inventory of 20 items**, each with its origin. All are removed, except two classes related to viscosity (`HTK_stratified`, `dK_list`), kept for the B-18 decision. `linear_system_ADIM`, a non-dimensional solver that is not wired in, is not dead code and stays | 0.6.0 | No output changes | references |
| B-20 | The argument goes through `set_layers` | 0.8.0 | No output changes | no dedicated test |
| B-25 | Values are read before the priors are replaced. Since 0.12.0, priors are no longer in `Layer` and the method no longer exists | 0.8.6 | References unchanged | not applicable |
| B-33 | The MCMC goes through `Column.simulate` and no longer depends on the state of the column | 0.12.0, bi 0.4.0 | References unchanged | `test_missing_measurements_mcmc.py`, without solving the column first |
| B-36 | `checker` only sets a flag | 0.11.1 | References unchanged | `test_checker.py` |

### 2.3 What remains to be decided

**Rule followed since 9/11.** The digital twin (Molonaviz 2, repository `flipoyo/MOLONARI-API`) has its own direct model and its own DREAM, written in JAX. It compares its direct model with `pyheatmy`, with a relative tolerance of 10⁻⁹. We changed nothing that would force it to modify its code: neither the physics of the direct model, nor the energy, nor the DREAM algorithm. Each decision below that touches one of these three points is therefore made **with the twin team** and requires new references.

#### Physics of the direct model

**B-27: the lateral source term q_s.**

- **What we know.** The term added to the diagonal is proportional to T, in kelvins. There is no inflow temperature. A column at uniform temperature drifts (+0.13 K in 4 days for q_s = 10⁻⁸ s⁻¹).
- **The project manager's assumption.** The water brought by q_s is at the temperature of the column. In that case, the heat advection term should vanish.
- **Options:**
  - (a) T_inflow = local T: the heat term disappears, and only the water balance remains;
  - (b) a term q_s·ρw·cw·(T_inflow − T), with a T_inflow given by the user, for instance that of the river or of the aquifer.
- **In both cases:** the diagnostic `source_heat_flux`, which is in °C, must be made consistent with the matrix.
- **Points of attention.**
  - The twin has no q_s.
  - The tutorials use `q_s = 0` by default and explain how to switch it on.
  - With q_s estimated, the 100-iteration burn-in does not converge (R̂ up to 9).

**B-18: the viscosity μ(T).**

- **What we know.** The result is right for a constant-viscosity model, but the code suggests otherwise.
- **Options:**
  - (a) remove the update loop and document constant μ;
  - (b) reactivate the μ(T) law, in kelvins, then recompute K and H at every update, with the coupling T → μ → K → H.
- **Preparation already done.** `calc_K` already takes the viscosity as an argument.
- **Cost of option (b).** All references change, and the twin must do the same.
- **What is at stake.** Without μ(T), the estimated k depends on the season of the measurement (≈ 10 % on Point034).

**Specific storage `Ss = n / height`** (observation of §1.3).

- **To decide:** keep this definition, with its justification, or switch to Ss = ρg(α + nβ).
- **Expected effect.** The response time of the head would go from about 4000 s to less than 20 s.
- **Before deciding:** check the order of magnitude in Freeze and Cherry (1979).

#### MCMC algorithm: towards a single DREAM

The following three points go together. Together, they would let the twin use `bayesian_inversion` with its JAX model, instead of maintaining its own DREAM. Each one changes the random draws or the results: it requires new references **and** a statistical test, for instance on the 2D Gaussian of `test_samplers.py`.

**B-30: outlier chains.**

- **What the source says** (Vrugt, 2016, p. 284 and 286):
  - the "fitness" of each chain is the mean log-density over the **second half** of its samples;
  - an outlier test compares these N values: IQR by default (Upton and Cook, 1996), or the Grubbs, Peirce or Chauvenet tests;
  - a chain declared an outlier takes the position of **another chain drawn at random**;
  - according to the author, this correction "violates detailed balance (irreversible transition)", but it is sometimes necessary to reach convergence. As long as one chain stays in a disjoint region, R̂ cannot fall below 1.2.
- **What the source does not say.** It does not state that the correction is limited to the burn-in: the details are left to Vrugt et al. (2009a), which is not among our sources. As the correction is irreversible, the cautious reading is to apply it only during the burn-in, like the adaptation of the crossover probabilities, which the same source restricts to the burn-in "to not destroy reversibility".
- **Today:** tutorial 02 sets aside after the run, with an explicit rule, the chains whose mean energy exceeds Q3 + 2·IQR. It does not move them.

**B-29, continued: maximum length of the burn-in.**

- **What is done.** The criterion and the warning are fixed.
- **What remains.** The default of `run_mcmc` is `nitmaxburning = 25`. Tutorials 02 and 03 pass 100, and their burn-in stops by itself at about 50 and 60 iterations.
- **Cost of a change.** Changing the default would change the results of every caller, the twin included.
- **Good to know.** When the burn-in converges, `nb_burn_in_iter` counts one iteration less. This behaviour is documented and not fixed, because the twin relies on it.

**B-11, continued: a single random generator.**

- **What is done.** `run_mcmc(..., seed=...)` seeds both generators (`random` and `np.random`), and two runs with the same seed are identical (`test_mcmc.py`, `test_samplers.py`).
- **What remains.** A single `np.random.Generator`, passed explicitly. This is the direction compatible with the explicit keys of JAX.
- **Cost.** All random sequences change.

**Parallel DREAM.** Today, DREAM updates its chains one after the other (documented in `MCMC_README.md`). The standard version makes proposals for all chains from the same population. With an entry point to evaluate the model in batches (a `map` function), it would allow parallel computation. This has been possible since B-36: a column can be serialised.

#### Interface

**B-09: `sigma_meas_T` and `sigma_meas_P`.**

- **Today.** The user passes the measurement variance directly to `run_mcmc` (`sigma2 = SIGMA_T**2` in tutorial 03). The arguments of `Column` still have no effect.
- **Options:**
  - (a) remove them from `Column`;
  - (b) expose them to build `sigma2`.
- **Point of attention:** the twin reads these arguments.

**B-24: scale of the priors depending on the sign.**

- **Options:**
  - (a) a mirrored log scale for ranges that are entirely negative;
  - (b) extend the symlog scale.
- **Today:** no effect, since no current prior is entirely negative.
- **Where to decide:** in `StateSpace` and `Prior` of `bayesian_inversion`.

#### Data

**B-31, continued: the P508 calibration in `integration_tests/conftest.py`.**

- **What we know.** In the processed data of Point034 (and of Point021), the head difference is **already in metres** (column "Pression différentielle (m)"). The integration tests still apply the P508 sensor calibration to it as if it were a voltage. They get +0.85 to +0.87 m instead of −0.031 to −0.006 m: 45 times too much, and with the opposite sign. The notebook where the error came from, `demoInversion.ipynb`, has been removed (tutorial 02 replaces it).
- **Why it is not fixed yet.** These tests use the demo column as it was, to stay comparable with the references of the refactoring.
- **What has to be done:** remove the calibration, then check each integration test again with the new inputs.
- **The twin is not affected:** it calibrates the raw voltages.

**Estimating σ² on Point034.**

- **Trial** with the settings of tutorial 03: the burn-in reaches its limit (100) without converging, with about 4 % acceptance and σ ≈ 0.035 K.
- **Choice of tutorial 03:** σ is fixed at 0.1 K.
- **To try:** a longer chain, or other settings.

#### Small known points

- **Public accessors of `synthetic_MOLONARI`.** They are missing: tutorial 01 reads private attributes (`_dates`, `_dH`, …).
- **NaN in the boundary conditions.** They raise an explicit error. Interpolating them in time remains to be done.
- **`Point021-1` data.** One sensor is entirely NaN: the data load with a warning.
- **Subsampling in space.** The rows of the quantiles correspond to one cell out of n. The sensor indices in the full grid therefore do not apply directly to these rows.
- **Fixes without a dedicated test.** B-12, B-19 and B-20 have no test in the repositories; they were verified through the references.

---

## 3. Split into two packages (refactoring)

### 3.1 Goal and principles

**Goal.** At the start, the inversion was a method of the model (`Column.compute_mcmc`). We wanted:

- the inversion to work with any model, including the JAX model of the twin;
- each part to evolve and be tested on its own.

The project defines three access levels: *public* (documentation and tutorials), *user* (direct model) and *developer* (state space and MCMC). The goal was set in two stages:

- on 30/10: the inversion may use the direct model, never the reverse;
- on 6/11: **full independence**, neither package knows the other. The project manager considered this independence valuable.

**Principles followed at every step:**

- **Move before changing.** A step that moves code does not change its behaviour. The references (§2.1) stay identical to the bit, demonstration notebook included, unless the purpose of the step is to change them.
- **Change the API only once.** Steps B.3 and B.4 were designed together so that users (tutorials, twin) have a single migration to make.
- **Check every copy.** A rewritten file is transferred by copying, then compared with its reference (`git diff`, `cmp`, `diff -rq`). This check caught, in B.2, a forgotten `description` line.
- **Check by running that the packages are independent.** A log of the imported modules shows that no test process of one package loads the other.
- **Start from new repositories.** The package repositories were created without history (decision of 21/10), with the single branch `CGS-MOLO`. The history stays in MOLONARI-computo and MOLONARI1D.

**Two routes to full independence**, compared on 30/10:

- **(A)** a third package holds the glue code and depends on the other two;
- **(B)** a contract: `bayesian_inversion` expects a "parameters → outputs" function, and `pyheatmy` provides one, `Column.simulate`, without importing the inversion.

Route (B) was chosen (2/11). It avoids a third repository, and the MCMC becomes generic (numpy only).

### 3.2 The steps

| Step (date) | What changes | Versions | Verification |
|---|---|---|---|
| **A** (29/10) | **Pure move** to `bayesian_inversion`, which does not import pyheatmy: `Prior`, Gelman–Rubin, the energy, the acceptance probability, the DREAM proposal (a method of `Column` turned into the function `perturbation_DREAM`, since it did not use the column), their constants and their unit tests. `compute_mcmc` stays in pyheatmy | pyheatmy 0.10.0, bi 0.1.0 | 140 tests, as before; references 49/49 and notebook 9/9 identical to the bit |
| **A+** (29/10) | `bayesian_inversion` becomes an optional dependency of pyheatmy (`pyheatmy[inversion]`) | pyheatmy 0.10.1 | Without bi: the direct model runs (60 tests passed, 21 skipped), direct-model references 31/31 |
| **B.1** (30/10) | **Pure move, through a subclass**: `compute_mcmc`, its accessors and its plots (851 lines) go into `InversionColumn(Column)`, in bi. The MCMC integration tests and the Point034 data follow | pyheatmy 0.11.0, bi 0.2.0 | 140 tests; references identical; without bi, no inversion module loaded |
| Release (30/10) | `MOLONARI-computo-1` on the four repositories of the time, with tutorials 01 and 02 | — | Tutorials: same results and figures as with pyheatmy 0.10.1 |
| **B.2** (6/11) | **Sampler classes** that do not know the model (`Sampler`, `DREAM`, `RandomWalk`). They see it only through `model(x)` and `energy(output, sigma2)`, and advance one iteration at a time (`step`). The order of the random draws is kept | bi 0.3.0, 0.3.1 | 147 tests; references identical; new tests on a 2D Gaussian, without pyheatmy |
| B-36 (6/11) | A column can be copied and serialised: a precondition of `simulate` | pyheatmy 0.11.1 | References identical |
| **B.3 + B.4** (6/11) | **The two packages no longer import each other:** <br>• `Layer` holds only the physics, and the priors move to `StateSpace`; <br>• `Column.simulate` is the direct model; <br>• `run_mcmc(model, observations, state_space, …, seed)` replaces `compute_mcmc` and returns an `MCMCResult`; <br>• `InversionColumn` disappears, with the dependencies of bi on pyheatmy, matplotlib and psutil; <br>• the shared tests move to `integration_tests/` in the tree | pyheatmy 0.12.0, bi 0.4.0 | MCMC and notebook references identical to the bit; tests 66 (bi) + 68 (pyheatmy) + 31 (tree); tutorial 02 identical to the byte |
| Last development day (9/11) | **Continuous integration** (GitHub Actions, Python 3.10 and 3.13, each package installed alone). **Clean-ups:** LF line endings; an unused output file and an unused constant removed; MCMC README up to date; RWM proposes σ² from the MCMC space. **Tree:** tutorial 03; a single deletion commit in computo; tree description; repositories moved to forks of the project manager | pyheatmy 0.12.1, bi 0.4.1, 0.4.2 | References identical to the bit; first continuous-integration runs passed |
| Reorganisation (9–10/11, project manager) | **Repositories:** new `MOLO-*` repositories, and `molonari-pipeline` (`raw_to_processed`) joins the tree; `TD4_analyse_donnees.ipynb` moves to `MOLO-data-process`. **Environment:** one `pixi.toml` per repository, a SemVer version number kept by `scripts/bump_version.py`, continuous integration in pixi. **Code:** `lora_emitter.py` removed from pyheatmy; otherwise, the code of the two packages and the tutorials are unchanged (file-by-file comparison, 13/11) | pyheatmy 0.12.2, bi 0.4.3, tree 1.0.1 | 68 + 66 + 120 + 31 tests passed (13/11) |

### 3.3 The architecture today

```
user ──► Column (pyheatmy)                              StateSpace (bayesian_inversion)
           │ data, layers, mesh                           │ one block per layer: {parameter: Prior}
           │                                              │
           └── column.simulate ─────────────────► run_mcmc(model, observations, state_space, seed=…)
               (parameters → outputs)                     │ DREAM (several chains) or RWM (one chain)
                                                          ▼
                                                     MCMCResult: states, energies, acceptance rates,
                                                                 quantiles, best state
```

The only link between the two packages is the function passed as the first argument. The user connects them in a few lines (excerpt from tutorial 02):

```python
from pyheatmy import Column, Layer
from bayesian_inversion import Prior, StateSpace, run_mcmc

column = Column(river_bed=..., depth_sensors=..., offset=..., dH_measures=..., T_measures=..., nb_cells=...)
column.set_layers(Layer.from_dict({"name": "river bed", "zLow": ..., "IntrinK": 1e-13, ...}))
space = StateSpace({"river bed": {name: Prior(*PRIORS[name]) for name in PARAMETERS}})

result = run_mcmc(column.simulate, column.get_measured_temperatures(), space,
                  typealgo="no sigma", sigma2=SIGMA_T**2, nb_iter=..., nb_chain=...,
                  nitmaxburning=..., seed=SEED)
result.chains()       # physical parameters, (iterations + 1, chains, layers, parameters)
result.energies()     # energy of each state
result.best()         # parameters of the lowest-energy state, a valid argument of column.simulate
```

**The contract.** What `run_mcmc` expects from a model:

- `model(params)` receives a list, with one dictionary `{name: physical value}` per block (a layer);
- it returns a dictionary of arrays:
  - `output["sensors"]` is compared with the observations (Gaussian errors of variance `sigma2`, NaN left out);
  - the outputs named in `stored` (by default `"temperatures"` and `"flows"`, of shape cells × times) are used for the quantiles.

Any model that follows this contract can be inverted, without pyheatmy.

**Main modules.**

- **`pyheatmy`** (direct model):
  - `core.py`: `Column`, `simulate`, `get_measured_temperatures`;
  - `layers.py`: `Layer`;
  - `params.py`: `Param`, `calc_K`;
  - `linear_system.py` and `solver.py`: solvers for H and T, with `numba`;
  - `synthetic_MOLONARI.py`: synthetic data;
  - `val_analy.py`: analytical solution;
  - `checker.py`: order of the computations;
  - `config.py`: constants.
- **Modules of `pyheatmy` outside the main path:**
  - `linear_system_ADIM.py`: non-dimensional solver, not wired in;
  - `frequency.py`, `time_series_multiperiodic.py`.
- **`bayesian_inversion`** (inversion):
  - `prior.py`: `Prior`, with the automatic choice of scale;
  - `state_space.py`: `StateSpace`;
  - `energy.py`: energy and acceptance;
  - `dream.py`: DREAM proposal;
  - `samplers.py`: `Sampler`, `DREAM`, `RandomWalk`;
  - `diagnostics.py`: R̂;
  - `mcmc.py`: `run_mcmc`, `MCMCResult`;
  - `MCMC_README.md`: the algorithm, step by step.

### 3.4 Design decisions and their reasons

| Decision | Reason |
|---|---|
| `simulate` receives the parameters **by name** and is **pure**: it solves a copy and does not modify the column | A pure function can be evaluated in parallel and has the same form as a JAX model. Before 6/11, a non-pure version modified the column |
| Column-specific plots leave the library; the tutorials draw their own figures | `bayesian_inversion` no longer depends on matplotlib; each user chooses their figures |
| Tests that use both packages live in the tree (`integration_tests/`), not in a package | Each package is tested alone, without the other |
| The seed is a parameter of `run_mcmc` and of the synthetic generation. At this stage, it seeds the two existing generators | Make runs reproducible without changing the bits; a single generator is the next step (§2.3) |
| B-33 is closed by the new interface, not by a separate patch | The MCMC no longer depends on the state of the column: the problem disappears, as a test shows |
| No automatic translation of the code to JAX (option discussed with the project manager) | The translation made by the twin is not mechanical (`lax.scan`, precomputed geometry), whereas our solver uses in-place assignments, Python loops and `numba`. Instead: <br>• the twin's comparison test, tied to our releases; <br>• every change of physics announced in the commit message; <br>• in the long term, a single source of the physics, independent of the array library (§4.3) |

### 3.5 Tree, repositories and publication

- **cgitsync tree.** The tree is described by `molo-computo.cgs`, at the root of MOLONARI-computo: five repositories, all on the branch `CGS-MOLO`. The file must be the only `.cgs` at the root, because the search for nested trees picks up every `*.cgs`.
  - **Rebuilding the tree:** `cgitsync bootstrap <path>/molo-computo.cgs` rebuilds the tree in an empty directory (checked on 13/11: `repos=5`, `errors=0`).
  - **Link with MOLONARI1D:** the trees of MOLONARI1D (`.config/*.cgs`) mount computo in their `computo/` folder.
- **Repositories.** They are the project manager's repositories (§0), with continuous integration enabled. The `YunusEmre4418/MOLONARI-*` repositories keep the history from 0.10 to 0.12.1 and are no longer updated. After a change of address, two operations are needed: `git remote set-url` in each repository, then `cgitsync pull <.cgs>`, which updates the recorded state of the tree.
- **Publication.** `cgitsync release freeze "<message>"` puts the same tag on every repository of the tree. Careful: the command runs `git add --all` in each repository, so all of them must be clean. Run it first with `--dry-run`.
- **Content of MOLONARI-computo.** A single deletion commit (`b2ce91f`, reversible) removed:
  - the former `pyheatmy/` package 0.9.7;
  - the demonstration notebooks, replaced by the tutorials;
  - `configuration/` and the LoRa files;
  - the former projects of `research/`.

  `research/` keeps six notebooks from 2023 to 2025, as a historical record that is not maintained. `TD4_analyse_donnees.ipynb` has since moved to `MOLO-data-process`.
- **The twin** (Molonaviz 2, `flipoyo/MOLONARI-API`, branch `CGS-MOLO`):
  - it uses pyheatmy as the reference of its tests and pinned it to version 0.9.7 (`7ebcc57`, the former computo folder);
  - with the current version, its tests fail only on the import of the energy;
  - if the energy is taken from `bayesian_inversion`, its 12 direct-model tests pass within 10⁻⁹ (checked on 6/11 with pyheatmy 0.11.0).

  For pinning, it can use the commits of §0.

---

## 4. Where to start

### 4.1 The first day

1. **Rebuild the tree.** You need git, [pixi](https://pixi.sh) and a clone of [ComplexGitSync](https://github.com/flipoyo/ComplexGitSync). From that clone:

   ```sh
   pixi run cgitsync bootstrap <path>/molo-computo.cgs
   ```

   The file `molo-computo.cgs` can be downloaded from the root of MOLONARI-computo, branch `CGS-MOLO`. You need an up-to-date ComplexGitSync (`git pull`, then `pixi install` in its clone): the file uses a recent syntax. The five repositories are cloned at their place, on the branch `CGS-MOLO`, and `bootstrap` prints the `export CGSHOME=...` line of the tree. Then, from inside the tree: `pixi run --manifest-path <ComplexGitSync>/pixi.toml cgitsync status`. A healthy state shows `repos=5`, `dirty=0`, `recorded_mismatch=0` and `errors=0`.
2. **Install and test.** At the root of the tree:

   ```sh
   pixi install
   pixi run test-all
   ```

   Expected result: 68 (pyheatmy) + 66 (`bayesian_inversion`) + 120 (`molonari-pipeline`, one skipped) + 31 (integration) tests passed.

   `pixi install` may rewrite `pixi.lock` at the root: this modified file must not end up in a commit by accident.
3. **Run the tutorials**, in order: `pixi run lab` opens JupyterLab in `TUTO/`, and `pixi run test-tuto` runs the three notebooks. The seeds are fixed: on the same machine, two runs give the same outputs.

   | Tutorial | Content | Checkpoints for the installation |
   |---|---|---|
   | `01_generate_data` | synthetic data and their true values (`truth.json`); pyheatmy only | writes the same files as those in the repository |
   | `02_inversion` | inversion of the synthetic data, compared with the true values; MCMC of about 4 min | the burn-in stops by itself at about 50 iterations; R̂ of the kept iterations 1.07 and 1.09 |
   | `03_real_data` | inversion of the real data of Point034; MCMC of about 3 to 4 min | see below |

   Checkpoints for tutorial 03:
   - burn-in of about 60 iterations; R̂ 1.04 and 1.12;
   - λs from 4.28 to 4.41 W/m/K;
   - k bounded only from above (95 %: 9.3·10⁻¹⁴ m²);
   - RMSE from 0.03 to 0.04 K.
4. **Read, in this order:**
   - the text of the tutorials;
   - `bayesian-inversion/bayesian_inversion/MCMC_README.md`, the algorithm step by step;
   - the [developer guide](docs/README_DEV.md);
   - §1.3 and §2.3 of this document.

**About the tutorials.** They are committed without their outputs. Tutorial 03 fixes σ at 0.1 K, because the estimation of σ² does not converge with its settings (§2.3). k is constrained there only from above, for a physical reason: the flux is upward and slow (at most about 3 mm/day), and the Péclet number stays below 0.02. Heat is therefore carried by conduction, and the temperatures determine the thermal diffusivity rather than k.

### 4.2 How we worked

These rules made it possible to fix 28 findings without a regression. We recommend keeping them.

- **References before and after every change:**
  1. before modifying the code, record the outputs of a fixed set of scenarios (direct model, RWM, DREAM, synthetic data, tutorials), with a fixed seed;
  2. write the predicted signature;
  3. modify;
  4. compare to the bit (§2.1).

  An unexpected difference stops the work until it is explained.
- **One correctness test per fix.** Its expected value comes from theory, not from the code, and it must fail before the fix.
- **No number hard-coded in the code** (rule of the project manager). A physical or tunable value becomes a named constant in `config.py`. A structural number disappears if the expression can be rewritten.
- **One version per change**, with semantic versioning, for each repository: `pixi run bump-version patch` (or `minor`, `major`), checked by `pixi run check-version` and by continuous integration. The commit message names the package and its version, and follows the rule of ComplexGitSync (`CommitMessagePolicy`).
- **A change of physics, energy or algorithm is flagged** in the commit message and coordinated with the twin team (§2.3).
- **For the tree:**
  - always the branch `CGS-MOLO`, never `main`;
  - rehearse every operation on the tree (`freeze`, change of address, deletion) on a copy before doing it for real.
- **Language.** The code, the docstrings, the READMEs and the tutorials are in English.
- **Before deleting code, read its history** (§1.1). Code that carries the trace of a physical choice is assessed before it is deleted (example: B-18).

### 4.3 Work items, by priority

The first four are ready to be picked up. The following ones need a decision, with the project manager and, for physics, with the twin team.

| # | Work item | Why now | Starting point |
|---|---|---|---|
| 1 | **Bring the reference protocol back into the tree**: scenario scripts and bit-level comparison, for instance in `integration_tests/`. Add the missing tests for B-12, B-19 and B-20 | Cheap. It protects all the following items | §2.1; scenarios S1–S5 |
| 2 | **Remove the P508 calibration** from `integration_tests/conftest.py`, then check each integration test again | The head difference of Point034 is already in metres | §2.3, B-31 |
| 3 | **Public accessors** for `synthetic_MOLONARI`, then tutorial 01 uses them | The tutorial reads private attributes | §2.3 |
| 4 | **Tutorial 04**: multilayer column, or frequency analysis | Multilayer is the real case, and B-01, B-17 and B-19 made it consistent | `research/multilayer/`; notebooks at the tag `MOLONARI-computo-1` |
| 5 | **A single DREAM**: outlier chains (B-30), burn-in length (B-29), single generator (B-11), parallel DREAM with a `map` function | The twin could use `bayesian_inversion` with its JAX model instead of maintaining its own DREAM | §2.3; Vrugt (2016); first, statistical tests shared by the two DREAMs |
| 6 | **Physics decisions**: B-27 (q_s), B-18 (μ(T)), Ss | These are the scientific choices that remain open. Each one changes all references and the twin's model | §2.3; `calc_K` already takes the viscosity |
| 7 | **Interface**: B-09 (`sigma_meas_T`), B-24 (scale of negative priors) | Low urgency | §2.3 |
| 8 | **Content of the tree**: <br>• add MOLONARI-data, the twin and the firmware to the `.cgs`; <br>• bring the firmware into computo (`raw_to_processed` is there since the reorganisation); <br>• licences to complete (`TUTO` and `data-pre-proc` have no `LICENSE` file; the others are MIT) | Organisation of the project, decision of the project manager | `docs/README_DEV.md` |
| 9 | **Performance**: compile the time loop with `numba`, **after measuring**. One call of the direct model takes about 60 ms, that is about 40 µs per time step (measured on 6/11) | Useful for longer MCMC runs | `linear_system.py` |
| 10 | **Long term**: a single source of the physics, independent of the array library (numpy or JAX) | Avoids maintaining two direct models | §3.4 |

**Molonaviz (the first one).** The decision of 19/10 is to keep it: the pyheatmy code it uses is not dead code. But it has not worked with pyheatmy since November 2024, and Molonaviz 2 takes its place. Whether to adapt it to the current interface is a question to raise again with the project manager.

### 4.4 Known pitfalls

- **Units:**
  - temperatures in kelvins in the library;
  - IntrinK in m² (otherwise `ValueError`);
  - head difference of the processed data already in metres, hence without calibration.
- **DREAM settings:**
  - at least `2·delta + 1` chains are needed, that is 7 with `delta = 3`;
  - the maximum burn-in length is 25 iterations by default, which is short: the tutorials pass 100;
  - a warning flags a burn-in stopped at its limit without converging.
- **Output folder.** `Column` and `synthetic_MOLONARI` create an output folder (`rac`). In scripts and tests, pass a temporary folder.
- **Freeze.** `cgitsync release freeze` runs `git add --all` in each repository: check beforehand with `git status --ignored`. The root of the tree must contain only one `.cgs`.
- **Searching the code.** Search tools such as `rg` or `ugrep` ignore the files listed in `.gitignore`. In the MOLONARI1D tree, the root `.gitignore` thus hides the nested repositories. Use `grep` without a filter, or switch the filter off.
- **macOS, pixi environment.** If the Python of `.pixi/envs/default` dies at once (exit code 137), even with `python -V`, look at `~/Library/Logs/DiagnosticReports/` before suspecting memory. The message `Code Signature Invalid` there points to a signature cached by the system. The remedy is to give the file a new inode, with the same content: `mv python3.13 python3.13.bak && cp -p python3.13.bak python3.13`.

---

## Appendix A: sources

- Anderson, M. P. (2005). Heat as a ground water tracer. *Ground Water*, 43(6), 951–968.
- Brooks, S. P., Gelman, A. (1998). General methods for monitoring convergence of iterative simulations. *Journal of Computational and Graphical Statistics*, 7(4), 434–455.
- Chib, S., Greenberg, E. (1995). Understanding the Metropolis–Hastings algorithm. *The American Statistician*, 49(4), 327–335.
- Freeze, R. A., Cherry, J. A. (1979). *Groundwater*. Prentice-Hall. Order of magnitude of Ss quoted from memory, to be checked (§1.3).
- Gelman, A., Rubin, D. B. (1992). Inference from iterative simulation using multiple sequences. *Statistical Science*, 7(4), 457–472.
- Roberts, G. O., Gelman, A., Gilks, W. R. (1997). Weak convergence and optimal scaling of random walk Metropolis algorithms. *Annals of Applied Probability*, 7(1), 110–120.
- Romary, T. *Bayesian inverse problems*. Lecture notes, Mines Paris (slide 39: likelihood in (σ²)^(−n/2); slide 63: monitoring of the energy).
- Stallman, R. W. (1965). Steady one-dimensional fluid flow in a semi-infinite porous medium with sinusoidal surface temperature. *Journal of Geophysical Research*, 70(12), 2821–2827. Reference for the analytical solution of `val_analy.py`, cited without being checked again in the article.
- Upton, G., Cook, I. (1996). *Understanding Statistics*. Oxford University Press. IQR test for outlier chains.
- Vrugt, J. A. (2016). Markov chain Monte Carlo simulation using the DREAM software package: theory, concepts, and MATLAB implementation. *Environmental Modelling & Software*, 75, 273–316. R̂: eq. 32–35; outlier chains: p. 284 and 286.
- Vrugt, J. A., ter Braak, C. J. F. (2011). DREAM(D): an adaptive MCMC simulation algorithm to solve discrete, noncontinuous, and combinatorial posterior parameter estimation problems. *Hydrology and Earth System Sciences*, 15, 3701–3713.
- Vrugt, J. A., ter Braak, C. J. F., Diks, C. G. H., Higdon, D., Robinson, B. A., Hyman, J. M. (2009a). Accelerating Markov chain Monte Carlo simulation by differential evolution with self-adaptive randomized subspace sampling. *International Journal of Nonlinear Sciences and Numerical Simulation*, 10(3), 273–290. Not consulted; it is the reference for the details of the outlier-chain correction (B-30).

## Appendix B: vocabulary

| Term | Meaning in this document |
|---|---|
| Reference (golden) | Recorded outputs of the code for a fixed input; detects a change, does not prove correctness |
| Predicted signature | List, written before a fix, of the outputs that must change and of those that must not |
| Correctness test | Test whose expected value comes from theory (hand calculation, analytical solution, known distribution) |
| Mutation | Deliberate undoing of part of a fix, to check that a test detects it |
| Physical space / MCMC space | Value of a parameter, and its transform (log, symlog or linear) in which the chains move |
| Block | Group of parameters of `StateSpace`; for a column, a layer |
| Energy | Minus the log posterior density (up to a constant) |
| Burn-in | Beginning of the chains, discarded; DREAM adapts its crossover probabilities there and stops when R̂ ≤ 1.2 for every free parameter |
| Twin | Molonaviz 2 (`flipoyo/MOLONARI-API`): direct model and DREAM in JAX, with an HTTP API |
| Tree, release | Set of repositories managed by cgitsync; common tag put on all of them by `release freeze` |

---

Good luck with what comes next!

Yunus Emre Aydogan, EV26
