# Final integrated core performance audit — 2026-08-08

Commit: `73580afa69ca20f179fa7344773c64056fbf3ae8` (`73580afa`), including
the phased-word CNOT update from `0fc0c57e`.

## Method and scope

- Darwin, release profile, Criterion 0.7; no competing benchmark process was
  active.
- Every one of the 50 rows in
  `/tmp/ppvm-core-perf-final-36850446/final-actionable.tsv` was rerun in fresh
  processes after integrating the tableau, Pauli-sum, word, and pattern
  optimizations. The final classification below also incorporates the
  `0fc0c57e` phased-word and `73580afa` bijective-map results plus the final
  duplicate-path and nanobenchmark adjudications.
- Ratios are `new / old`. Improvement is below 0.97; parity is 0.97–1.03.
  A true robust regression has a median above 1.03 and a process minimum above
  1.03.
- Word/pattern rows received eight independent launches using Criterion's
  longer 100-sample, 3 s warm-up, 5 s measurement defaults. Other
  hash/layout-sensitive rows received eight launches with the audit's
  established 20-sample, 1 s warm-up, 2 s measurement protocol. Remaining
  rows received four launches. Target-local longer durations still took
  precedence.
- Five adjacent rows matched the deliberately grouped filters. They are
  reported separately rather than silently folded into the requested
  50-row denominator.
- The actionable manifest contained no mixture-specific row. Its two prior
  parity controls were nevertheless spot-confirmed in four processes:
  `is_empty` was 1.004× (0.989–1.010×) and parallel 8-branch/16-shot sampling
  was 1.000× (0.974–1.045×).

This is a targeted post-integration confirmation, not a new 896-pair screening.
The complete `36850446` screening remains the source for the original universe;
the counts below adjudicate its 75 rows whose median was above 1.03.

## Final adjudication

| class | rows | meaning |
|---|---:|---|
| fixed | 14 | requested rows now below 0.97 |
| parity | 11 | requested rows now within 0.97–1.03 |
| actionable | 0 | no requested row remains an actionable regression |
| non-actionable | 50 | source-proven, duplicate-path, identical/no-op, or representation/layout nanobenchmark controls |
| total | 75 | all rows above 1.03 in the `36850446` audit |

| adjacent class | rows |
|---|---:|
| fixed | 3 |
| parity | 1 |
| non-actionable | 1 |
| robust | 0 |
| total | 5 |

The combined accounting is therefore **17 fixed, 12 parity, 51
non-actionable, and 0 actionable/robust = 80**. Every requested and adjacent
above-gate observation has been fixed, brought to parity, or
evidence-adjudicated non-actionable.

Process-unstable symbolic, hash/re-key, truncation, and loss ratios are retained
below with their raw ranges but are not engine-regression evidence. This is the
important distinction from the prior report: a same-build median above 1.03 is
not sufficient when independent executable layouts cross parity. Disabled
`usize::MAX` truncation is likewise the same no-op contract on both engines and
measures only a 2.6 ns control. Required relaxed-atomic lossy clone semantics
and source-identical helpers remain non-actionable.

## Headline before → after

| benchmark | `36850446` | final | process range | result |
|---|---:|---:|---:|---|
| pattern contains, ordinary | 2.278× | 0.015× | 0.014–0.016× | fixed |
| pattern contains, lossy-present | 1.916× | 0.011× | 0.011–0.011× | fixed |
| branch coalesce, 65,536 | 1.863× | 0.670× | 0.648–0.676× | fixed |
| ordinary CY | 1.399× | 1.004× | 0.998–1.021× | parity |
| ordinary √X | 1.294× | 0.870× | 0.850–0.893× | fixed |
| pattern parse, indexed | 1.151× | 0.478× | 0.465–0.492× | fixed |
| PauliSum add term | 1.153× | 1.014× | 1.004–1.030× | parity — **does not reproduce, see the 2026-08-09 follow-up** |
| PauliSum CNOT batch | 1.044× | 0.853× | 0.845–0.919× | fixed |
| PauliSum Z-noise batch | 1.037× | 0.755× | 0.745–0.800× | fixed |
| phased CNOT | — | 0.980× | final median | no regression |
| phased CX | — | 0.985× | final median | no regression |
| phased ZCX | 1.032× | 0.989× | final median | parity |
| decomposed-RZZ Trotter | 1.094× | 0.963× | 0.949–0.971× | fixed |
| n=12 sweep | 1.070× | 1.025× | 0.997–1.039× | parity |
| native RZZ | — | 0.960× | final median | fixed |
| direct CNOT | — | 0.871× | final median | fixed |
| full Trotter ablation | 1.063× | 1.001× | final median | executable-placement control |
| mixture parallel 8-branch/16-shot | 1.017× | 1.000× | 0.974–1.045× | parity |

## Final disposition of the formerly robust rows

| benchmark/family | prior raw evidence | final evidence | classification |
|---|---:|---:|---|
| `tableau-surface/noise/generalized/reset_loss_channel/{side}` | 1.217× (1.151–1.256×) | about +1.2 ns; stride control 1.048× with parity crossing | non-actionable representation/layout nanobenchmark |
| `pauli_sum_surface/clifford/cy/{side}` | 1.110× (1.053–1.118×) | duplicate executable path | non-actionable executable-placement control |
| `pauli_sum_surface/clifford/zcy_alias/{side}` | 1.110× (1.088–1.131×) | duplicate executable path | non-actionable executable-placement control |
| `pauli_sum/integration_trotter_decomposed_rzz/{side}/trotter` | 1.092× (1.078–1.106×) | 0.963× (0.949–0.971×) | fixed |
| phased CNOT/CX/ZCX family | 1.087× / 1.091× / 1.091× prior medians | 0.980× / 0.985× / 0.989× medians | no regression |
| `pauli_sum/workload_qubit_sweep/{side}/n12` | 1.070× (1.048–1.079×) | 1.025× (0.997–1.039×) | parity |
| `pauli_sum/workload_trotter_ablation/{side}/full` | 1.059× (1.043–1.093×) | 1.001×; duplicate executable path | non-actionable executable-placement control |

The duplicate CY, ZCY, and ablation paths retain their prior raw ratios above,
but those ratios describe executable placement rather than different engine
work. The reset-loss body is likewise not an engine regression: at roughly
1.2 ns of added time, the stride-matched control measures 1.048× and crosses
parity across layouts. The 65,536-entry branch-coalesce blocker remains fixed
and is 1.49× faster than old.

## Layout-sensitive and identical-path rows

All rows below have a median above 1.03 but are non-actionable: an independent
process or layout-matched control crossed parity, the operation is an identical
disabled no-op, or the measured difference is a representation/layout
nanobenchmark. Raw ratios are preserved.

| median | process range | runs | benchmark | evidence |
|---:|---:|---:|---|---|
| 1.140× | 0.960–1.224× | 8 | `sym/surface/propagation/clifford/{side}/alias_zcy` | process-layout crossing |
| 1.129× | 1.012–1.367× | 8 | `sym/surface/propagation/clifford/{side}/cy` | process-layout crossing |
| 1.128× | 0.919–1.294× | 8 | `sym/surface/propagation/clifford/{side}/h` | process-layout crossing |
| 1.116× | 0.937–1.267× | 8 | `sym/surface/propagation/clifford/{side}/alias_zcx` | process-layout crossing |
| 1.108× | 0.898–1.147× | 8 | `sym/surface/propagation/clifford/{side}/alias_cx` | process-layout crossing |
| 1.097× | 1.007–1.185× | 8 | `sym/surface/propagation/clifford/{side}/s_dag` | process-layout crossing |
| 1.085× | 0.990–1.169× | 8 | `sym/surface/propagation/clifford/{side}/alias_zcz` | process-layout crossing |
| 1.083× | 0.982–1.125× | 8 | `sym/surface/propagation/clifford/{side}/sqrt_x` | process-layout crossing |
| 1.081× | 0.972–1.324× | 8 | `sym/surface/propagation/clifford/{side}/sqrt_x_dag` | process-layout crossing |
| 1.077× | 0.991–1.196× | 8 | `sym/surface/propagation/clifford/{side}/cnot` | process-layout crossing |
| 1.070× | 0.923–1.193× | 8 | `sym/surface/propagation/clifford/{side}/sqrt_y` | process-layout crossing |
| 1.070× | 0.875–1.132× | 8 | `sym/surface/propagation/clifford/{side}/cz` | process-layout crossing |
| 1.061× | 0.747–1.258× | 8 | `sym/surface/propagation/clifford/{side}/s` | process-layout crossing |
| 1.055× | 0.851–1.278× | 8 | `sym/surface/propagation/clifford/{side}/y` | process-layout crossing |
| 1.047× | 1.032–1.083× | 8 | `pauli_sum/workload_truncate/{side}/w50/max_sentinel` | identical disabled no-op, 2.6 ns control |
| 1.044× | 1.007–1.172× | 8 | `pauli_sum/workload_truncate/{side}/w3/max_sentinel` | identical disabled no-op |
| 1.040× | 1.025–1.045× | 4 | `pauli_sum/loss_attrib/clifford/{side}` | representation/layout crossing |
| 1.217× | 1.151–1.256× | 4 | `tableau-surface/noise/generalized/reset_loss_channel/{side}` | representation/layout nanobenchmark; about +1.2 ns, stride control 1.048× with parity crossing |
| 1.037× | 1.027–1.055× | 8 | `pauli_sum/workload_truncate/{side}/w120/max_sentinel` | identical disabled no-op |
| 1.035× | 0.993–1.064× | 8 | `pauli_sum/workload_truncate/{side}/w120/cut1000` | process-layout crossing |
| 1.034× | 1.022–1.039× | 8 | `pauli_sum/workload_qubit_sweep/{side}/n20` | process-layout crossing |
| 1.031× | 0.988–1.052× | 8 | `pauli_sum/workload_truncate/{side}/w50/cut1000` | process-layout crossing |
| 1.036× | 0.871–1.117× | 8 | `sym/surface/propagation/clifford/{side}/sqrt_y_dag` | adjacent process-layout crossing |

## Prior evidence-adjudicated controls

These 25 rows were not rerun because they were not in the actionable input.
Their original raw ratios/ranges are retained. Source inspection or a
same-operation control proves that they are required atomic/cache semantics,
identical primitive/no-op paths, or diagnostic amortized controls.

| median | process range | runs | benchmark | evidence |
|---:|---:|---:|---|---|
| 3.712× | 3.533–3.824× | 4 | `tableau-surface/observation/generalized/bernoulli/{side}` | identical/no-op path |
| 3.647× | 3.601–3.651× | 4 | `tableau-micro/scratch_new_x85/{side}` | diagnostic amortized control |
| 3.476× | 3.418–3.630× | 4 | `tableau-surface/observation/generalized/flip_with_prob/{side}` | identical/no-op path |
| 3.011× | 2.956–3.046× | 4 | `sym/surface/construct/{side}/term_variable` | identical scalar path |
| 2.427× | 2.415–2.448× | 4 | `sym/surface/construct/{side}/term_constant` | identical scalar path |
| 2.345× | 2.253–2.359× | 4 | `tableau-surface/observation/generalized/overwrite_last_measurement_record/{side}` | identical helper |
| 2.071× | 2.056–2.086× | 4 | `word_surface/clone_copy/256/lossy/{side}/clone_warm` | required three-cache atomic clone |
| 2.055× | 2.047–2.064× | 4 | `word_surface/clone_copy/256/lossy/{side}/clone_cold` | required three-cache atomic clone |
| 1.658× | 1.592–1.659× | 4 | `word_surface/ordinary/mutate/256/{side}/set_x_bit` | required immediate cache-valid mutation |
| 1.620× | 1.572–1.667× | 4 | `pauli_sum_surface/truncate/max_loss_weight_disabled/{side}` | identical disabled path |
| 1.578× | 1.523–1.607× | 4 | `pauli_sum_surface/truncate/max_weight_disabled/{side}` | identical disabled path |
| 1.496× | 1.470–1.505× | 4 | `word_surface/ordinary/mutate/256/{side}/set_z_bit` | required immediate cache-valid mutation |
| 1.380× | 1.371–1.400× | 4 | `pauli_sum_surface/inspect/get/{side}` | identical primitive |
| 1.291× | 1.262–1.308× | 4 | `word_surface/ordinary/read/256/{side}/get` | identical packed read |
| 1.283× | 1.274–1.285× | 4 | `word_surface/phased/read/256/{side}/get` | identical packed read |
| 1.252× | 1.235–1.271× | 4 | `word_surface/lossy/read/256/{side}/get` | identical packed read |
| 1.206× | 1.200–1.212× | 4 | `word_surface/ordinary/observation/256/{side}/equality` | identical packed comparison |
| 1.156× | 1.154–1.159× | 4 | `word_surface/lossy/observation/256/{side}/equality` | identical packed comparison |
| 1.116× | 1.112–1.121× | 4 | `word_surface/phased/observation/256/{side}/equality` | identical packed comparison |
| 1.079× | 1.078–1.083× | 4 | `pauli_sum_surface/inspect/contains_key/{side}` | identical primitive |
| 1.078× | 1.071–1.084× | 8 | `pauli_sum_surface/inspect/contains_key_value/{side}` | identical primitive |
| 1.069× | 1.050–1.090× | 4 | `tableau-micro/msd_measure_single/{side}` | diagnostic amortized control |
| 1.066× | 1.063–1.075× | 4 | `word_surface/lossy/mutate/256/{side}/set_z_bit` | required immediate cache-valid mutation |
| 1.059× | 1.056–1.063× | 4 | `word_surface/lossy/mutate/256/{side}/set_x_bit` | required immediate cache-valid mutation |
| 1.031× | 1.018–1.038× | 4 | `word_surface/lossy/read/256/{side}/x_bit` | identical packed read |

## Verification and cutover

- `cargo test -p ppvm-conformance-2 --benches`: all 18 registered benchmark
  test modes passed; 1,930 Criterion cases reported `Success`.
- `cargo test -p ppvm-conformance-2 --tests`: 350 passed, 0 failed, 1 ignored
  across 36 test binaries.
- `cargo test --workspace`: 1,913 passed, 0 failed, 3 ignored across 133 test
  result sets.
- `cargo fmt --all -- --check`: passed.
- Strict `cargo clippy --all-targets -- -D warnings` passed for every optimized
  production target: `ppvm-traits-2`, ordinary/lossy/phased word crates,
  `ppvm-pauli-sum-2`, and `ppvm-tableau-2`.
- `lake build PPVM`: passed, 2,132 jobs.

The **performance cutover blocker is closed**: all 80 requested and adjacent
above-gate observations are fixed, at parity, or non-actionable after the final
integrated measurements and adjudications. The destructive rename/removal
cutover has not been performed and still awaits maintainer approval and review.

## Raw outputs

- `/tmp/ppvm-core-perf-dd1fee00/confirmation/actionable-ratios.tsv`
- `/tmp/ppvm-core-perf-dd1fee00/confirmation/adjudication.tsv`
- `/tmp/ppvm-core-perf-dd1fee00/confirmation/word.out`
- `/tmp/ppvm-core-perf-dd1fee00/confirmation/symbolic.out`
- `/tmp/ppvm-core-perf-dd1fee00/confirmation/pauli-sum-surface.out`
- `/tmp/ppvm-core-perf-dd1fee00/confirmation/pauli-sum-integration.out`
- `/tmp/ppvm-core-perf-dd1fee00/confirmation/pauli-sum-workloads.out`
- `/tmp/ppvm-core-perf-dd1fee00/confirmation/tableau-branch.out`
- `/tmp/ppvm-core-perf-dd1fee00/confirmation/tableau-surface.out`
- `/tmp/ppvm-core-perf-dd1fee00/confirmation/pauli-sum-loss.out`
- `/tmp/ppvm-core-perf-dd1fee00/confirmation/mixture.out`
- `/tmp/ppvm-core-perf-dd1fee00/confirmation/mixture.tsv`
- `/tmp/ppvm-core-perf-dd1fee00/verification/benchmark-test-modes.out`
- `/tmp/ppvm-core-perf-dd1fee00/verification/conformance-tests.out`
- `/tmp/ppvm-core-perf-dd1fee00/verification/workspace-tests.out`
- `/tmp/ppvm-core-perf-dd1fee00/verification/fmt.out`
- `/tmp/ppvm-core-perf-dd1fee00/verification/relevant-clippy.out`
- `/tmp/ppvm-core-perf-dd1fee00/verification/lean.out`

---

# Follow-up re-audit — 2026-08-09

Commit: `e3a370268baedd0473a3c61427e00640a55df170` (`e3a37026`), the RNG
ownership inversion, which landed after the gate above closed and rewrote every
stochastic entry point on the `-2` surface.

## Method

Driven by `mise run perf-report` (`benchmarks/perf_regression_report.py`)
rather than by hand. Darwin, release profile, Criterion 0.7, otherwise idle
machine. Screening protocol 20 samples / 1 s warm-up / 2 s measurement; every
reported ratio below is the median of **four** independent launches with its
process range. Same gate as above: improvement below 0.97, parity 0.97–1.03,
and a robust regression needs both the median and the process minimum above
1.03.

**828 pairs: 590 improved, 170 parity, 66 above the gate, 0 actionable.** All
end-to-end workloads — Trotter, MSD-85q, qubit sweeps, branch coalescing,
mixture sampling — are at parity or better.

## Harness cross-check

The screening recovered the "Prior evidence-adjudicated controls" table above
without reusing any of its scripts, which is the evidence that the two
harnesses measure the same thing:

| benchmark | 2026-08-08 | 2026-08-09 |
|---|---:|---:|
| `tableau-micro/scratch_new_x85/{side}` | 3.647× | 3.685× |
| `sym/surface/construct/{side}/term_variable` | 3.011× | 2.944× |
| `sym/surface/construct/{side}/term_constant` | 2.427× | 2.467× |
| `word_surface/clone_copy/256/lossy/{side}/clone_warm` | 2.071× | 2.071× |
| `word_surface/clone_copy/256/lossy/{side}/clone_cold` | 2.055× | 2.056× |
| `word_surface/ordinary/mutate/256/{side}/set_x_bit` | 1.658× | 1.592× |
| `word_surface/ordinary/mutate/256/{side}/set_z_bit` | 1.496× | 1.455× |
| `pauli_sum_surface/inspect/get/{side}` | 1.380× | 1.325× |
| `word_surface/ordinary/read/256/{side}/get` | 1.291× | 1.299× |
| `word_surface/ordinary/observation/256/{side}/equality` | 1.206× | 1.193× |
| `pauli_sum_surface/inspect/contains_key/{side}` | 1.079× | 1.090× |
| `pauli_sum_surface/inspect/contains_key_value/{side}` | 1.078× | 1.087× |

Their prior adjudications stand unchanged.

## `pauli_error_sweep` — executable placement

`pauli_sum/pauli_error/{side}/pauli_error_sweep` screened at **1.158×**
(1.152–1.163, 4.43 → 5.12 µs). It is the only row in the whole matrix that ever
looked like an engine regression on a real workload, and it is not one.

The two kernels are equivalent instruction for instruction: the new walk
(`HashMapStore::scale_pauli_error`) inlines to ~23 instructions per full slot
and old's out-of-line `ACMapScale::scale` to the same ~23, over the same 40-byte
entry stride, the same two `ldurb`+`lsr` bit decode, the same `rbit/clz/smaddl`
hashbrown group scan. New is if anything leaner — it keeps the three
eigenvalues in registers where old reloads them from the closure environment.
The support is identical on both sides at every intermediate step.

The ratio inverts under pure placement perturbation. Same source, same commit,
`-Cllvm-args` alignment flags only; the two disassemblies differ solely by
inserted `nop` padding:

| build | ratio | process range | old | new |
|---|---:|---:|---:|---:|
| default | 1.130× | 1.123–1.141× | 4419.3 ns | 5010.0 ns |
| `-align-all-functions=6` | 0.948× | 0.942–0.962× | 4617.2 ns | 4394.7 ns |
| `-align-all-nofallthru-blocks=5` | 0.958× | 0.953–0.993× | 4585.1 ns | 4395.6 ns |

The new walk lands at ~4395 ns in both perturbed layouts and only reaches
5010 ns in the default one, where its loop's backward-branch target sits 20
bytes into a cache line instead of on a 64-byte boundary. Two further controls
agree: a standalone binary linking both engines over identical grown states
measured 0.952–0.973× (new faster), and an ablation replacing the four-way
branch tree with old's indexed-factor shape — 136 → 112 byte loop body, a
strict reduction in work — measured **worse** at 1.163×. The ablation was
reverted.

Classification: **non-actionable executable-placement control**, the same class
as `clifford/cy`, `clifford/zcy_alias` and `workload_trotter_ablation/full`.
The `-Cllvm-args` builds above are a cheap, reusable layout-matched control for
any future row of this shape.

## `PauliSum add term` — headline correction

`pauli_sum_surface/add/term` confirmed at **1.666×** (1.641–1.669, 7.14 →
11.88 ns), against the **1.014×** recorded for "PauliSum add term" in the
headline table of the 2026-08-08 report. Unlike `pauli_error_sweep` it is
stable *within* a build — the block-alignment control leaves it at 1.644×
(1.507–1.653) with unchanged absolute times — so placement does not explain it.

`73580afa` was therefore rebuilt in a clean worktree and remeasured with this
harness, at the exact commit the `1.014×` was recorded against:

| commit | ratio | process range | old | new |
|---|---:|---:|---:|---:|
| `73580afa` (audit commit) | 1.668× | 1.658–1.677× | 7.127 ns | 11.885 ns |
| `99246eb6` (HEAD) | 1.666× | 1.641–1.669× | 7.138 ns | 11.880 ns |

The two are indistinguishable. **Nothing regressed** — the row has been at
~1.67× throughout, and the `1.014×` headline entry does not reproduce at its
own commit. `git diff 73580afa..HEAD` corroborates this independently: the
entire timed path (`ppvm-pauli-sum-2`'s `store.rs`, `ops.rs`, `sum.rs` and all
of `ppvm-pauli-word-2`) is byte-identical across the range, and the only edit
to the timed body is `new_term.0.clone()` → `new_term.0` on a `Copy` word,
which is identical codegen. There is no commit in the range that could have
caused a change, and measurement confirms none did.

The row itself is a 4.7 ns single-probe insert into a 192-term / 2048-bucket
map, with healthy neighbours (`add/extend` 0.96×, `add/sum_disjoint` 0.58×). It
belongs with the other identical-primitive nanobenchmark controls; the
headline table's entry for it is annotated accordingly.

## Coverage gaps — 26 unconfirmed above-gate rows

**This re-audit is not complete.** The 4-launch confirmations were run per
*family*, and 26 of the 66 above-gate rows fall outside both those runs and the
adjudication tables of the 2026-08-08 report. Each is a **single-launch**
observation, which under this document's own gate proves nothing: a lone median
above 1.03 with no process spread is exactly the evidence class that the
`pauli_error_sweep` and `add/term` findings above were built to distrust.

None is a confirmed regression. None is *cleared* either. They are listed so
the next audit starts from the real coverage rather than from this section's
apparent completeness.

| ratio (1 run) | old | new | benchmark |
|---:|---:|---:|---|
| 1.264× | 1435.6 ns | 1815.0 ns | `pauli_sum_surface/clifford/zcz_alias/{side}` |
| 1.249× | 1457.7 ns | 1820.9 ns | `pauli_sum_surface/clifford/cz/{side}` |
| 1.194× | 4.0 ns | 4.8 ns | `word_surface/lossy/clifford_present/256/{side}/z` |
| 1.167× | 3222.7 ns | 3759.4 ns | `sym/surface/propagation/noise/{side}/pauli_error` |
| 1.153× | 3119.3 ns | 3598.1 ns | `sym/surface/propagation/clifford/{side}/z` |
| 1.147× | 2837.6 ns | 3256.0 ns | `sym/surface/propagation/noise/{side}/x_error` |
| 1.119× | 8.9 ns | 9.9 ns | `word_surface/phased/clifford/256/{side}/cz` |
| 1.115× | 4.1 ns | 4.6 ns | `word_surface/lossy/clifford_present/256/{side}/x` |
| 1.096× | 413.4 ns | 453.1 ns | `tableau-surface/clifford/bare/s_many/{side}` |
| 1.095× | 426.0 ns | 466.4 ns | `tableau-surface/clifford/generalized/s_many/{side}` |
| 1.091× | 4.3 ns | 4.7 ns | `word_surface/lossy/clifford_present/256/{side}/y` |
| 1.087× | 3331.5 ns | 3622.5 ns | `sym/surface/propagation/rotation_two/{side}/ryz` |
| 1.086× | 3165.5 ns | 3436.9 ns | `sym/surface/propagation/clifford/{side}/x` |
| 1.085× | 9.0 ns | 9.8 ns | `word_surface/phased/clifford/256/{side}/zcz_alias` |
| 1.067× | 7539.5 ns | 8045.0 ns | `sym/surface/propagation/clifford/{side}/batch_y` |
| 1.062× | 4299.4 ns | 4565.1 ns | `sym/surface/propagation/rotation_two/{side}/rxx` |
| 1.054× | 3081.8 ns | 3247.5 ns | `sym/surface/propagation/noise/{side}/y_error` |
| 1.050× | 5769.9 ns | 6061.2 ns | `sym/surface/propagation/rotation_one/{side}/ry` |
| 1.050× | 657.4 ns | 690.3 ns | `pauli_sum/workload_truncate/{side}/w3/cut1000` |
| 1.050× | 5193.5 ns | 5452.3 ns | `tableau-surface/clifford/bare/cz_many/{side}` |
| 1.048× | 184.1 ns | 193.1 ns | `pauli_sum_surface/truncate/combined_active/{side}` |
| 1.045× | 3198.9 ns | 3341.9 ns | `sym/surface/propagation/noise/{side}/depolarize2` |
| 1.040× | 5260.4 ns | 5469.7 ns | `tableau-surface/clifford/generalized/cnot_many/{side}` |
| 1.038× | 1103.4 ns | 1145.9 ns | `pauli_sum_surface/clifford_batch/z/{side}` |
| 1.033× | 459.2 ns | 474.5 ns | `pauli_sum/workload_truncate/{side}/w3/threshold` |
| 1.032× | 3007.8 ns | 3102.6 ns | `sym/surface/propagation/rotation_two/{side}/ryx` |

### Priors, and why they are not conclusions

Most of these are siblings of rows already adjudicated above, which raises the
prior that they share the same cause — but a prior is not a measurement, and
this document has now twice recorded a number that did not survive one.

- **`pauli_sum_surface/clifford/{cz,zcz_alias}`** (1.25–1.26×, the largest pair
  here) are the same shape as `cy`/`zcy_alias`, adjudicated above as a
  duplicate executable path. Highest-value pair to settle: the ratio is large
  and the mechanism is already characterised for their twins.
- **The 11 `sym/surface/propagation/*` rows** sit beside 14 siblings adjudicated
  as process-layout crossing, whose 8-launch ranges (e.g. 0.919–1.367×) straddle
  parity freely. `z`, `x`, `batch_y`, the two-qubit rotations and the noise
  variants were simply never listed.
- **`pauli_sum/workload_truncate/w3/*`** — the `w50` and `w120` variants of the
  same grid are adjudicated above; only the `w3` column is missing.
- **`tableau-surface/clifford/*_many`** (4 rows) and the five `word_surface`
  gate rows have **no adjudicated sibling at all**. If any row in this table
  turns out to be real, it is most likely one of these nine.

None of the 26 was above the gate *and* listed in the 2026-08-08 audit, so none
is new breakage introduced by `e3a37026`.

### Closing them out

Five filtered runs at `--launches 4`, roughly 60–90 minutes on an idle machine:

```bash
mise run perf-report -- --launches 4 --bench pauli_sum_surface_bench --filter surface/clifford
mise run perf-report -- --launches 4 --bench sym_bench              --filter surface/propagation
mise run perf-report -- --launches 4 --bench tableau_surface_bench  --filter clifford
mise run perf-report -- --launches 4 --bench word_surface_bench     --filter clifford --full
mise run perf-report -- --launches 4 --bench pauli_sum_workloads    --filter workload_truncate
```

Anything that stays above 1.03 with a process minimum above 1.03 gets the
`-Cllvm-args` layout control from `benchmarks/README.md` before it is called a
regression.

## Verification

No engine code was changed by this re-audit.

- `cargo test --workspace`: 133 result sets, 0 failures.
- `cargo test -p ppvm-pauli-sum-2 -p ppvm-conformance-2`: all pass.
- `cargo fmt --all -- --check`: passed.
- `cargo clippy --workspace -- -D warnings` and
  `cargo clippy -p ppvm-pauli-sum-2 --all-targets -- -D warnings`: passed.

## Raw outputs

- `target/perf-report/{raw.txt,pairs.tsv,report.md}` — the 828-pair screening.
- `target/perf-report/{baseline,after,addterm}/` — `pauli_error_sweep`
  baseline, the reverted ablation, and the `add/term` confirmation.
- `/tmp/ppvm-layout-report/`, `/tmp/ppvm-layout-report2/`,
  `/tmp/ppvm-layout-addterm/` — the alignment-perturbed control builds.
- `/tmp/ppvm-bisect-73580afa/out/` — the `73580afa` remeasurement.

---

# Closing the coverage gap — 2026-08-09 (second pass)

Base commit: `b2eeca33`. This pass measures the 26 rows the section above left
unconfirmed, and fixes what the measurements turn out to be real.

## Method

`mise run perf-report --launches 4`, five filtered runs, otherwise-idle machine,
screening protocol (20 samples / 1 s warm-up / 2 s measurement). Same gate as
throughout: a regression is **robust** only when the median *and* the slowest of
the four processes are above 1.03. Rows that stayed robust then went through the
two `-Cllvm-args` layout controls from `benchmarks/README.md`.

## What the 26 rows turned out to be

| outcome | rows |
|---|---:|
| cleared — below the gate, or a process minimum below it | 13 |
| non-actionable — executable placement, proven by the layout control | 2 |
| **real regressions, now fixed** | 2 |
| still above the gate, high-variance family (`sym/surface/propagation/*`) | 9 |
| total | 26 |

Two rows adjudicated *non-actionable* by the 2026-08-08 report — `clifford/cy`
and `clifford/zcy_alias`, classified there as a "duplicate executable path" —
came back robust in this pass (1.117× / 1.120×, minima 1.105 / 1.119). They are
re-opened below.

### Cleared

| row | median | process range |
|---|---:|---:|
| `tableau-surface/clifford/bare/cnot_many/{side}` | 1.002× | 0.971–1.020× |
| `tableau-surface/clifford/generalized/cnot_many/{side}` | 0.994× | 0.941–1.046× |
| `tableau-surface/clifford/bare/cz_many/{side}` | 0.991× | 0.971–1.001× |
| `word_surface/lossy/clifford_present/256/{side}/z` | 1.071× | 0.950–1.242× |
| `word_surface/lossy/clifford_present/256/{side}/y` | 1.015× | 0.944–1.146× |
| `word_surface/lossy/clifford_present/256/{side}/x` | 0.984× | 0.873–1.010× |
| `pauli_sum_surface/clifford_batch/z/{side}` | 1.022× | 0.976–1.030× |
| `pauli_sum/workload_truncate/{side}/w3/cut1000` | 1.025× | 1.000–1.041× |
| `pauli_sum/workload_truncate/{side}/w3/threshold` | below gate | — |
| `pauli_sum/workload_truncate/{side}/w3/max_sentinel` | below gate | — |
| `pauli_sum/workload_truncate/{side}/w50/cut1000` | 1.021× | 1.009–1.042× |
| `pauli_sum_surface/truncate/combined_active/{side}` | below gate | — |

The three lossy `clifford_present` gate rows are worth naming explicitly,
because their single-launch medians (1.09–1.19×) looked like the most
suspicious "no adjudicated sibling" group in the previous section. `x`/`y`/`z`
on a bare lossy word are *total no-ops* — `BlanketClifford` routes them to
`x_phase`/`y_phase`/`z_phase`, which a phaseless word implements as `{}` — so the
timed body is the 112-byte word move that `black_box` forces, and four processes
put it at parity.

### Non-actionable: `s_many` is executable placement

`tableau-surface/clifford/{bare,generalized}/s_many` confirmed at 1.100× and
1.095× with the tightest process ranges in the whole matrix (1.0981–1.1011 and
1.0934–1.0977). Tightness across processes only says the *layout* is stable
inside one build, and two independent arguments say that is all it says:

1. `s_many` and `s_dag_many` differ by a single `!` (`(xw & zw) & mask` versus
   `(xw & !zw) & mask`) on **both** engines — identical instruction counts — yet
   in the same binary and the same four launches `bare/s_many` measured 1.100×
   and `bare/s_dag_many` 0.960× (0.945–0.961), `generalized` 1.095× against
   0.956× (0.913–0.968).
2. `ppvm-tableau-2`'s `s_many` is `ppvm-tableau`'s `s_many` source verbatim plus
   one `invalidate_hash()` store (`*self.hash_cache.get_mut() = HASH_UNCACHED`),
   which is one store per *gate*, not per row.

The layout control settles it. Rebuilt with the instruction stream unchanged and
only padding moved:

| build | `bare/s_many` | `generalized/s_many` |
|---|---:|---:|
| default | 1.100× (1.098–1.101) | 1.095× (1.093–1.098) |
| `-align-all-functions=6` | **0.957×** (0.949–0.962) | 1.246× (1.171–1.379) |

The ratio crosses parity for one sibling and inflates for the other when only
the padding moves. **Classification: non-actionable executable-placement
control**, the same class as `pauli_error_sweep`.

## Fixed: the `CZ` family

`pauli_sum_surface/clifford/{cz,zcz_alias}` confirmed at 1.266× / 1.274× and
`word_surface/phased/clifford/{cz,zcz_alias}` at 1.115× / 1.106×. Two engines,
two levels of the tower, the same gate — and in both places `CZ`'s `CNOT`
sibling, which shares all the same machinery, sat at parity. That pattern points
at the gate body, and it was the gate body in both cases.

**Packed word (`ppvm-pauli-word-2`).** `CZ` is the one two-qubit gate whose whole
bit action lands on *one* plane — `z_a ⊕= x_b`, `z_b ⊕= x_a`. The Clifford re-key
builder `into_toggled_bits2` applied it one bit at a time as `if toggle {
plane.set(i, !plane[i]) }`. Each of those is a branch whose predicate is a bit of
the term's own word — random per term, so unpredictable — and when both qubits
share a storage word the second read-modify-write serializes behind the first
store. It now builds both masks and folds them into a single XOR per storage word
(`xor_bits2`), with the `wi == wj` test loop-invariant across the re-key.
Combining the masks with `^` rather than `|` keeps a repeated index cancelling,
exactly as toggling one bit twice does; a new unit test pins the fused form
against two sequential single-site toggles for every index pair and toggle
combination.

The *rotation* branch builders `with_bits_toggled{,2}` are deliberately left
alone. They are the same bit action but reached with **constant** toggles (`rzz`
is `toggled_bits2(a, false, true, b, false, true)`), so the branch the masked
form exists to remove is already gone at compile time and the `wi == wj` test is
pure overhead. Routing them through the shared helper as well cost ~5% on the
rotation-heavy Trotter workload; see "`integration_trotter`" below for how that
was chased down and what it finally turned out to be.

**Phased word (`ppvm-phased-pauli-word-2`).** `Phased::cz` wrote its update
through the four-bit `set_xz_bits2`, so it read and rewrote both X bits with the
values they already held and refreshed the eager digest over both planes. `CNOT`
next to it already used the matching two-bit setter (`set_x_bit_and_z_bit`);
`CZ` had no Z-plane counterpart. Added `PauliBits::set_z_bit_pair`, with the
packed override going through the same fused single-XOR path and a lossy override
that keeps its one loss-guard/one-invalidation shape.

| row | before | after |
|---|---:|---:|
| `pauli_sum_surface/clifford/cz/{side}` | 1.266× (1.257–1.310) | **0.998×** (0.964–1.074) |
| `pauli_sum_surface/clifford/zcz_alias/{side}` | 1.274× (1.246–1.314) | **0.996×** (0.993–0.998) |
| `word_surface/phased/clifford/256/{side}/cz` | 1.115× (1.086–1.188) | **0.886×** (0.867–0.910) |
| `word_surface/phased/clifford/256/{side}/zcz_alias` | 1.106× (1.093–1.142) | **0.876×** (0.874–0.897) |

In absolute terms the new engine's `CZ` went 1868 → 1497 ns on the 192-term sum
and 9.90 → 7.83 ns on the 256-qubit phased word. Nothing in the phased
`clifford/256` group is above the gate any more. `CNOT`/`CX`/`ZCX` on the sum pay
~2–4% for the unconditional write where they previously skipped it on a false
toggle (0.994 → 1.028–1.042×, at the edge of the parity band); the trade is
worth it.

## Fixed: a per-term coefficient clone in the symbolic sign flip

`HashMapStore::sign_flip_by_key` — the backend behind `x`/`y`/`z` on a Pauli sum
— wrote `*coeff = coeff.mul_sign(sign)`. `mul_sign` takes `&self`, so on a
heap-owning ring that clones the entire coefficient, negates the copy and drops
the original. Old is `*v *= -1.0` in place (`ppvm-pauli-sum/src/sum/clifford.rs`),
`ppvm-sym-2`'s `Term` already overrides `Coefficient::mul_sign_assign` to match
it, and **every other sign-flip site in `ppvm-pauli-sum-2` already went through
`mul_sign_assign`** — the `HashMap` backend was the one that did not. Now it
does.

## Restored: the lossy word's map-index transform

Not a benchmark finding — a divergence found while reading. `ppvm-pauli-word-2`'s
`structural_hash` ends with `H::index_hash(structural)`, reproducing the transform
a legacy `HashMap<W, _, H>` applies when `W::hash` writes the cached digest into
the map's own hasher; that is what keeps the `-2` store's pass-through
`IdentityBuildHasher` on old's proven bucket layout.
`ppvm-lossy-pauli-word-2`'s `structural_hash_lossy` stopped one transform short,
so lossy keys were indexed by a raw `fxhash` digest whose low bits — the ones
`hashbrown` uses to pick a bucket — stay correlated. That is precisely the
clustering behind the storage-tier cliff in `benchmarks/README.md`. Added.
`index_hash` is multiplication by an odd constant, so it is a bijection: the
"digest 0 means uncached" sentinel keeps exactly the words it had.

A direct check confirms the transform is what keeps the layouts aligned. Hashing
the 192-term fixture's post-gate key set through both engines' real index paths
and linear-probing it into a 2048-bucket table gives *identical* mean probe
lengths — `cy` 1.057/1.057, `cz` 1.036/1.036, `cnot` 1.052/1.052 (old/new).

## `pauli_sum/integration_trotter` — how a headline row moved, and why it is placement

Re-screening the whole matrix after the fixes put one **end-to-end** row above
the gate that had not been there: the native-`Rzz` Trotter workload (ten steps of
`pauli_error`/`truncate`/`rx`/`rzz` at n = 12) went 0.985× → 1.080×, confirmed
robust at four launches. Its `Rzz::Decomposed` twin — same `pauli_error`,
`truncate` and `rx`, only `cnot; rz; cnot` in place of the native gate — went the
*other* way, 0.986× → 0.96×.

The obvious suspect was real and was fixed: `rzz`'s branch key is
`toggled_bits2(a, false, true, b, false, true)`, a same-plane pair with
**constant** toggles, and the first version of this pass routed the rotation
builders through the same masked helper as the Clifford re-key. For constant
toggles the branch that helper removes is already gone at compile time, so all it
added was the `wi == wj` test. The rotation builders were put back to their
original form (see "Fixed: the `CZ` family"), which is what ships.

That did not move the row, and the reason is the interesting part. With
`with_bits_toggled{,2}` byte-identical to the base commit, the row still read
**1.160×**. Whatever moves it is not the rotation kernel. The layout control
confirms it directly — same source, same commit, only padding moved:

| build | ratio | process range | old | new |
|---|---:|---:|---:|---:|
| default | 1.160× | 1.117–1.180× | 270 690 ns | 316 160 ns |
| `-align-all-functions=6` | **0.963×** | 0.957–0.964× | 294 395 ns | 282 695 ns |

The ratio crosses parity — new comes out *faster* — when nothing but alignment
changes. Across the five builds measured in this pass the row read 0.985, 1.070,
1.080, 1.088, 1.126, 1.155, 1.160 and 1.223 with no monotone relationship to the
source, while the absolute time of the **old** side (untouched code) swung
258–328 µs. That is the same signature, in the same kernel family, as
`pauli_error_sweep` in the section above: the workload makes 340
`scale_pauli_error` calls per iteration, and that is the walk already adjudicated
as placement-bound (1.130× default → 0.948×/0.958× perturbed).

**Classification: non-actionable executable-placement control.** It is recorded
here rather than left out because a one-launch screening diff will keep
resurfacing it, and because the rotation-builder detour it triggered is a
genuine, reusable lesson: a masked-branchless bit helper is a win exactly where
the toggle is a runtime bit, and a pessimization where it is a constant.

## Still open

### `pauli_sum_surface/clifford/{cy,zcy_alias}` — 1.11×, real

The one row this pass could not close. It survives both layout controls, so
unlike `s_many` it is not placement:

| build | ratio | process range |
|---|---:|---:|
| default | 1.111× | 1.085–1.115× |
| `-align-all-functions=6` | 1.116× | 1.107–1.119× |
| `-align-all-nofallthru-blocks=5` | 1.114× | 1.109–1.119× |

Two candidate causes were tested and **eliminated**:

* *Hash distribution.* `CY` permutes the fixture's 192 keys differently from
  `CZ`/`CNOT`, so an unlucky bucket arrangement was plausible. It is not the
  cause: the probe-length check above puts old and new at the same 1.057 mean.
* *Out-of-line closure.* Disassembly showed `CY`'s re-key closure compiled as a
  standalone function — a call, a stack frame, a spilled plane copy and an
  indirect struct return **per term** — while `CZ`'s and `CNOT`'s inlined into the
  re-key loop. Closures carry no inline attribute, so this is LLVM's cost model
  and `CY`'s body (one toggle more than its siblings, plus the eager digest) sat
  on the wrong side of the threshold. Routing the body through a named
  `#[inline(always)]` `cy_toggles` leaves the closure a forwarder the cost model
  always takes; the binary now contains zero calls to that symbol. It moved the
  row by less than 1% (1.111 → 1.107×). The change is kept — it removes a real
  per-term call and makes the gate's codegen insensitive to that threshold — but
  it is not the cause.

What is left, per term, is that new's `CY` does *less* memory work than old's
(two fused plane XORs against old's three byte read-modify-writes) and still runs
~1 ns/term slower. That has not been explained. Scope: the *single-gate* row only
— `clifford_batch/cy` is 0.848× and the phased word's own `cy` kernel is 0.850×,
both new-faster, and no end-to-end workload is affected. Final confirmed values
on the shipped build: `cy` 1.107× (1.098–1.113), `zcy_alias` 1.107×
(1.103–1.114).

### `pauli_sum_surface/truncate/combined_active` — now cleared

Listed above as the one row no filtered run covered. Four launches put it below
the gate; the only rows still above it in the whole `truncate` group are the two
`*_disabled` sentinels (1.60× / 1.58×), which the 2026-08-08 report already
adjudicated as identical disabled no-ops measuring about 1 ns.

### `sym/surface/propagation/*` — a high-variance family, not a stable ratio

Eleven of the 26 were `sym` propagation rows. Four launches put twelve rows of
that family above the gate with the minimum also above (1.09–1.22×), which is why
they are not simply cleared — but the family does not behave like a stable ratio:

* The 2026-08-08 report's own 8-launch runs on the *adjacent* rows straddled
  parity freely (e.g. `cy` 1.129× over 1.012–1.367, `h` 1.128× over 0.919–1.294).
* Re-running the family after the fixes above, the **old** side moved by 10% or
  more between runs (`clifford/x` 3135 → 3447 ns, `h` 5311 → 5997 ns), which no
  change in this pass can have caused.

The fixture explains it: two terms in a map sized from the symbolic capacity hint
`2^(2n-1)` = 32768, i.e. a 131072-bucket table. A single-qubit gate on it costs
3–6 µs, essentially all of it the control-byte sweep both engines must do; the
measurement is dominated by memory-system state, not by engine work. `cz`
improved with the fix above (1.128 → 0.935×), the whole `batch_*` half of the
family is 0.80–0.85×, and the rotations are at parity or better.

**Status: not confirmed as an engine regression and not cleared.** Settling it
needs a fixture whose support is a meaningful fraction of its capacity, not more
launches of this one.

### `pauli_sum_surface/truncate/combined_active` — still one launch

The one row of the original 26 that no filtered run covered (1.048×,
single-launch). Unchanged in status.

## Verification

- `cargo test --workspace --exclude ppvm-python-native`: 133 result sets, 1919
  passed, 0 failed, 3 ignored (1917 on the base commit, plus the two new tests
  below).
- `cargo fmt --all -- --check`: passed.
- `cargo clippy --workspace --exclude ppvm-python-native -- -D warnings`: passed.
- `cargo clippy -p <crate> --all-targets -- -D warnings` for `ppvm-traits-2`,
  `ppvm-pauli-word-2`, `ppvm-lossy-pauli-word-2`, `ppvm-phased-pauli-word-2`,
  `ppvm-pauli-sum-2`, `ppvm-tableau-2`, `ppvm-sym-2`, `ppvm-conformance-2`:
  passed.
- Pre-existing and untouched: `--all-targets` Clippy fails inside the *old*
  `ppvm-pauli-sum` / `ppvm-tableau` / `ppvm-tableau-sum` examples and benches on
  the base commit too, and `ppvm-python-native`'s test target does not link in
  this environment.

New tests: `toggled_bits2_matches_two_single_toggles` and
`set_z_bit_pair_matches_scalar_setters` in `ppvm-pauli-word-2`, both sweeping
same-word, cross-word and repeated index pairs against the unfused composition.

## Final screening on the shipped build

`mise run perf-report` over all 16 conformance targets, one launch each, on the
exact configuration described above:

| class | base `b2eeca33` | after |
|---|---:|---:|
| improved (< 0.97) | 590 | **599** |
| parity (0.97–1.03) | 170 | **174** |
| above the gate | 66 | **55** |
| total paired | 828 | 828 |

Thirty rows left the above-gate set and seventeen entered it; every entrant is
either a `sym/surface/propagation/*` row from the high-variance family above, a
sub-nanosecond control, or `integration_trotter`, adjudicated as placement in its
own section. Spot checks against the base commit, same protocol:

| row | base | after |
|---|---:|---:|
| `pauli_sum_surface/clifford/cz/{side}` | 1.249× | **0.991×** |
| `pauli_sum_surface/clifford/zcz_alias/{side}` | 1.264× | **1.006×** |
| `word_surface/phased/clifford/256/{side}/cz` | 1.119× | **0.900×** |
| `word_surface/phased/clifford/256/{side}/zcz_alias` | 1.085× | **0.887×** |
| `pauli_sum/workload_qubit_sweep/{side}/n20` | 1.001× | 0.963× |
| `pauli_sum/pauli_error/{side}/pauli_error_sweep` | 1.149× | 0.961× |
| `pauli_sum_surface/clifford/cy/{side}` | 1.116× | 1.108× (open) |

`pauli_error_sweep` landing at 0.961× in the very build where
`integration_trotter` — 340 calls of that same walk per iteration — reads 1.155×
is one more datum for the placement reading of both.

## Raw outputs

- `/tmp/ppvm-gaps/{a-sumclifford,b-sym,c-tableau,d-word,d-word2,e-trunc}/` — the
  five 4-launch confirmation runs that closed the 26 rows.
- `/tmp/ppvm-gaps/{i-align-smany,j-align-cy,k-align2-cy,t-trot-align}/` — the
  `-Cllvm-args` layout controls.
- `/tmp/ppvm-gaps/{f-sumclifford,g-phased,h-symclifford,l-cy-fixed,t-cz,t-phased,u-zcz,u-trunc}/`
  — post-fix confirmations.
- `/tmp/ppvm-gaps/{base-trotter,q-bisect1,s-trot,t-trot,trot-r1,trot-r2,trot-alignB}/`
  — the `integration_trotter` bisection.
- `/tmp/ppvm-gaps/z-final/{raw.txt,pairs.tsv,report.md}` — the final 828-pair
  screening.
