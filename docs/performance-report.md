# Core old/new benchmark report — 2026-08-07

This report compares the latest `-2` core crates with the old reference crates.
Every comparable public operation is represented; operations without a valid old
semantic twin are listed separately rather than assigned an artificial ratio.

## Method

- Platform: Darwin, release profile, Criterion 0.7.
- Full screening: every comparative target, 10 samples, 0.2 s warm-up, 0.5 s measurement.
- Confirmation: targeted screening regressions and complete surface suites were
  rerun in fresh Cargo-launched processes, 20 samples, 1 s warm-up,
  2 s measurement. The third complete-surface pass was stopped during the
  pathological mixture branch-scaling case; no fourth full pass began.
  Rows with fewer than three measurements remain provisional.
- Ratio is `new / old`; below 0.97 is improvement, 0.97–1.03 is parity.
  Above 1.03 is confirmed with at least three processes and provisional otherwise.
- Setup is excluded unless construction or cloning is itself the target.

## Summary

- Comparable benchmark pairs: **901**
- Improvements: **587**
- Parity: **162**
- Confirmed regressions: **102**
- Provisional regressions: **50**

## Regressions

| status | ratio | process range | runs | old | new | benchmark |
|---|---:|---:|---:|---:|---:|---|
| confirmed regression | **3.948×** | 3.429–4.073× | 7 | 2.012 ns | 7.952 ns | `tableau-surface/observation/generalized/flip_with_prob/{side}` |
| confirmed regression | **3.881×** | 3.481–4.011× | 7 | 2.016 ns | 7.987 ns | `tableau-surface/observation/generalized/bernoulli/{side}` |
| provisional regression | **3.506×** | 3.506–3.506× | 1 | 74.001 ns | 259.460 ns | `tableau-micro/scratch_new_x85/{side}` |
| confirmed regression | **2.560×** | 2.263–2.669× | 7 | 1.733 ns | 4.522 ns | `tableau-surface/observation/generalized/overwrite_last_measurement_record/{side}` |
| provisional regression | **2.213×** | 2.213–2.213× | 1 | 3.783 ns | 8.371 ns | `word_surface/clone_copy/256/lossy/{side}/clone_warm` |
| provisional regression | **2.202×** | 2.202–2.202× | 1 | 3.763 ns | 8.284 ns | `word_surface/clone_copy/256/lossy/{side}/clone_cold` |
| provisional regression | **1.999×** | 1.975–2.023× | 2 | 0.837 ns | 1.674 ns | `sym/surface/sum/{side}/add_const` |
| provisional regression | **1.962×** | 1.962–1.962× | 1 | 2.835 ns | 5.561 ns | `word_surface/clone_copy/256/ordinary/{side}/clone_cold` |
| provisional regression | **1.945×** | 1.924–1.966× | 2 | 0.798 ns | 1.551 ns | `sym/surface/operator_add/{side}/sum_add_coefficient` |
| provisional regression | **1.924×** | 1.924–1.924× | 1 | 2.869 ns | 5.520 ns | `word_surface/clone_copy/256/ordinary/{side}/clone_warm` |
| provisional regression | **1.902×** | 1.894–1.911× | 2 | 256.305 ns | 487.575 ns | `word_surface/pattern/match_contains/256/{side}/ordinary` |
| provisional regression | **1.801×** | 1.801–1.801× | 1 | 3.227 ns | 5.813 ns | `word_surface/clone_copy/256/phased/{side}/clone_cold` |
| confirmed regression | **1.799×** | 1.793–1.809× | 7 | 6.818 ns | 12.298 ns | `word_surface/lossy/branch_key/256/{side}/one_site_clone_then_bits` |
| provisional regression | **1.787×** | 1.787–1.787× | 1 | 3.258 ns | 5.821 ns | `word_surface/clone_copy/256/phased/{side}/clone_warm` |
| provisional regression | **1.774×** | 1.709–1.840× | 2 | 309.025 ns | 548.515 ns | `word_surface/pattern/match_contains/256/{side}/lossy_present` |
| provisional regression | **1.713×** | 1.644–1.781× | 2 | 2.611 µs | 4.475 µs | `sym/surface/propagation/{side}/h` |
| confirmed regression | **1.668×** | 1.597–1.691× | 7 | 1.574 ns | 2.609 ns | `pauli_sum_surface/truncate/max_weight_disabled/{side}` |
| confirmed regression | **1.668×** | 1.624–1.690× | 7 | 1.574 ns | 2.632 ns | `pauli_sum_surface/truncate/max_loss_weight_disabled/{side}` |
| provisional regression | **1.630×** | 1.516–1.745× | 2 | 2.661 µs | 4.353 µs | `sym/surface/propagation/{side}/s` |
| provisional regression | **1.604×** | 1.521–1.687× | 2 | 2.598 µs | 4.151 µs | `sym/surface/propagation/{side}/cz` |
| confirmed regression | **1.598×** | 1.577–1.620× | 7 | 7.166 ns | 11.556 ns | `word_surface/lossy/branch_key/256/{side}/two_site_clone_then_bits` |
| provisional regression | **1.558×** | 1.558–1.558× | 1 | 6.600 ns | 10.283 ns | `word_surface/lossy/mutate/256/{side}/set_x_bit` |
| provisional regression | **1.541×** | 1.541–1.541× | 1 | 6.709 ns | 10.337 ns | `word_surface/lossy/mutate/256/{side}/set_z_bit` |
| confirmed regression | **1.533×** | 1.455–1.538× | 4 | 1.765 ns | 2.659 ns | `pauli_sum/workload_truncate/{side}/w120/max_sentinel` |
| confirmed regression | **1.533×** | 1.470–1.635× | 4 | 1.740 ns | 2.675 ns | `pauli_sum/workload_truncate/{side}/w3/max_sentinel` |
| confirmed regression | **1.519×** | 1.492–1.550× | 4 | 1.758 ns | 2.686 ns | `pauli_sum/workload_truncate/{side}/w50/max_sentinel` |
| confirmed regression | **1.477×** | 1.433–1.560× | 6 | 5.961 ns | 8.598 ns | `sym/surface/construct/{side}/fold_cos_constant` |
| provisional regression | **1.448×** | 1.334–1.561× | 2 | 2.952 µs | 4.248 µs | `sym/surface/propagation/{side}/rx` |
| confirmed regression | **1.414×** | 1.407–1.443× | 3 | 300.820 µs | 425.420 µs | `pauli_sum/loss_attrib/clifford/{side}` |
| provisional regression | **1.410×** | 1.352–1.469× | 2 | 2.803 µs | 3.947 µs | `sym/surface/propagation/{side}/cnot` |
| provisional regression | **1.392×** | 1.392–1.392× | 1 | 0.510 ns | 0.710 ns | `word_surface/lossy/hash_protocol/256/{side}/warm` |
| confirmed regression | **1.388×** | 1.360–1.428× | 6 | 5.867 ns | 8.059 ns | `sym/surface/construct/{side}/fold_sin_constant` |
| confirmed regression | **1.374×** | 1.315–1.411× | 7 | 1.134 ns | 1.579 ns | `pauli_sum_surface/inspect/get/{side}` |
| confirmed regression | **1.339×** | 1.325–1.352× | 7 | 211.620 ns | 283.030 ns | `pauli_sum_surface/construct/empty/{side}` |
| provisional regression | **1.334×** | 1.223–1.444× | 2 | 3.120 µs | 4.135 µs | `sym/surface/propagation/{side}/ry` |
| confirmed regression | **1.329×** | 1.274–1.402× | 7 | 211.310 ns | 285.450 ns | `mixture/measure/lost/{side}` |
| provisional regression | **1.305×** | 1.305–1.305× | 1 | 374.110 ns | 488.040 ns | `pauli_sum_surface/compare/abs_diff_eq_near/{side}` |
| confirmed regression | **1.304×** | 1.258–1.372× | 6 | 81.417 ns | 108.250 ns | `sym/surface/eval/{side}/sum` |
| provisional regression | **1.294×** | 1.287–1.300× | 2 | 763.795 µs | 988.305 µs | `sym/trace_readout_k3/{side}/trace` |
| confirmed regression | **1.292×** | 1.235–1.314× | 7 | 16.145 ms | 20.861 ms | `pauli_sum/loss_attrib/loss/{side}` |
| confirmed regression | **1.291×** | 1.287–1.293× | 3 | 4.063 µs | 5.253 µs | `pauli_sum_surface/clifford_batch/s_dag/{side}` |
| provisional regression | **1.291×** | 1.147–1.436× | 2 | 2.002 µs | 2.552 µs | `sym/surface/propagation/{side}/rz` |
| confirmed regression | **1.290×** | 1.286–1.294× | 3 | 4.031 µs | 5.205 µs | `pauli_sum_surface/clifford_batch/s/{side}` |
| confirmed regression | **1.289×** | 1.236–1.300× | 6 | 1.502 ns | 1.936 ns | `sym/surface/construct/{side}/sum_new` |
| confirmed regression | **1.288×** | 1.283–1.289× | 3 | 0.941 ns | 1.212 ns | `word_surface/ordinary/read/256/{side}/get` |
| provisional regression | **1.284×** | 1.284–1.284× | 1 | 376.950 ns | 484.110 ns | `pauli_sum_surface/compare/abs_diff_eq_equal/{side}` |
| confirmed regression | **1.284×** | 1.281–1.294× | 3 | 0.944 ns | 1.211 ns | `word_surface/phased/read/256/{side}/get` |
| provisional regression | **1.280×** | 1.272–1.288× | 2 | 767.460 µs | 982.520 µs | `sym/trace_readout_k2/{side}/trace` |
| confirmed regression | **1.276×** | 1.268–1.295× | 3 | 5.216 µs | 6.708 µs | `pauli_sum_surface/clifford_batch/cy/{side}` |
| confirmed regression | **1.272×** | 1.175–1.288× | 7 | 417.320 ns | 532.110 ns | `pauli_sum_surface/construct/clone/{side}` |
| confirmed regression | **1.271×** | 1.243–1.319× | 3 | 5.176 µs | 6.581 µs | `pauli_sum/clifford_h/{side}/h` |
| provisional regression | **1.271×** | 1.261–1.282× | 2 | 783.155 µs | 995.660 µs | `sym/trace_readout_k4/{side}/trace` |
| confirmed regression | **1.265×** | 1.256–1.301× | 3 | 0.973 ns | 1.248 ns | `pauli_sum_surface/inspect/contains_key/{side}` |
| confirmed regression | **1.264×** | 1.260–1.286× | 3 | 4.167 µs | 5.275 µs | `pauli_sum_surface/clifford_batch/sqrt_x/{side}` |
| confirmed regression | **1.260×** | 1.257–1.270× | 3 | 5.206 µs | 6.589 µs | `pauli_sum/clifford_cnot/{side}/cnot` |
| confirmed regression | **1.256×** | 1.255–1.270× | 3 | 4.154 µs | 5.234 µs | `pauli_sum_surface/clifford_batch/sqrt_x_dag/{side}` |
| provisional regression | **1.248×** | 1.205–1.292× | 2 | 130.135 µs | 162.545 µs | `sym/surface/readout/{side}/trace` |
| confirmed regression | **1.237×** | 1.231–1.243× | 3 | 4.314 µs | 5.343 µs | `pauli_sum_surface/clifford_batch/h/{side}` |
| confirmed regression | **1.235×** | 1.228–1.237× | 3 | 4.327 µs | 5.346 µs | `pauli_sum_surface/clifford_batch/sqrt_y/{side}` |
| confirmed regression | **1.235×** | 1.232–1.238× | 3 | 1.230 µs | 1.522 µs | `pauli_sum_surface/clifford/s/{side}` |
| confirmed regression | **1.235×** | 1.230–1.235× | 3 | 1.222 µs | 1.510 µs | `pauli_sum_surface/clifford/s_dag/{side}` |
| confirmed regression | **1.233×** | 1.228–1.240× | 3 | 4.329 µs | 5.367 µs | `pauli_sum_surface/clifford_batch/sqrt_y_dag/{side}` |
| confirmed regression | **1.230×** | 1.227–1.233× | 7 | 1.027 ns | 1.261 ns | `word_surface/lossy/read/256/{side}/get` |
| confirmed regression | **1.225×** | 1.213–1.249× | 7 | 16.945 ms | 20.751 ms | `pauli_sum/loss_attrib/rotation/{side}` |
| confirmed regression | **1.221×** | 1.216–1.221× | 3 | 1.519 µs | 1.854 µs | `pauli_sum_surface/clifford/zcy_alias/{side}` |
| confirmed regression | **1.209×** | 1.202–1.216× | 3 | 1.528 µs | 1.854 µs | `pauli_sum_surface/clifford/cy/{side}` |
| provisional regression | **1.208×** | 1.108–1.309× | 2 | 4.449 µs | 5.359 µs | `sym/surface/propagation/clifford/{side}/s_dag` |
| confirmed regression | **1.207×** | 1.202–1.214× | 3 | 35.660 µs | 43.029 µs | `pauli_sum/loss_attrib/reset/{side}` |
| confirmed regression | **1.205×** | 1.193–1.242× | 4 | 487.960 ns | 589.720 ns | `pauli_sum/workload_truncate/{side}/w120/threshold` |
| confirmed regression | **1.203×** | 1.197–1.205× | 3 | 1.246 µs | 1.500 µs | `pauli_sum_surface/clifford/sqrt_x_dag/{side}` |
| confirmed regression | **1.198×** | 1.196–1.206× | 3 | 4.775 µs | 5.755 µs | `pauli_sum_surface/clifford_batch/cnot/{side}` |
| confirmed regression | **1.195×** | 1.178–1.198× | 3 | 1.274 µs | 1.501 µs | `pauli_sum_surface/clifford/sqrt_x/{side}` |
| confirmed regression | **1.194×** | 1.180–1.212× | 6 | 962.690 µs | 1.155 ms | `sym/trace_readout_k5/{side}/trace` |
| confirmed regression | **1.193×** | 1.179–1.193× | 3 | 1.307 µs | 1.559 µs | `pauli_sum_surface/clifford/sqrt_y_dag/{side}` |
| confirmed regression | **1.191×** | 1.189–1.212× | 3 | 1.307 µs | 1.567 µs | `pauli_sum_surface/clifford/h/{side}` |
| confirmed regression | **1.188×** | 1.155–1.208× | 3 | 1.049 µs | 1.217 µs | `pauli_sum_surface/noise_batch/depolarize2/{side}` |
| confirmed regression | **1.187×** | 1.157–1.210× | 7 | 4.479 ns | 5.341 ns | `tableau-surface/noise/generalized/reset_loss_channel/{side}` |
| confirmed regression | **1.183×** | 1.175–1.185× | 3 | 1.319 µs | 1.553 µs | `pauli_sum_surface/clifford/sqrt_y/{side}` |
| confirmed regression | **1.174×** | 1.128–1.216× | 7 | 7.124 ns | 8.123 ns | `pauli_sum_surface/add/term/{side}` |
| provisional regression | **1.173×** | 1.171–1.176× | 2 | 1.471 ns | 1.726 ns | `pauli_sum_surface/inspect/contains_key_value/{side}` |
| confirmed regression | **1.172×** | 1.158–1.174× | 3 | 1.410 µs | 1.655 µs | `pauli_sum_surface/clifford/cnot/{side}` |
| confirmed regression | **1.169×** | 1.150–1.181× | 3 | 1.401 µs | 1.638 µs | `pauli_sum_surface/clifford/zcx_alias/{side}` |
| confirmed regression | **1.164×** | 1.153–1.166× | 7 | 3.885 ns | 4.481 ns | `word_surface/construct/lossy/{side}/new_identity/8` |
| confirmed regression | **1.163×** | 1.146–1.170× | 7 | 3.904 ns | 4.472 ns | `word_surface/construct/lossy/{side}/new_identity/64` |
| confirmed regression | **1.163×** | 1.163–1.171× | 3 | 1.397 µs | 1.631 µs | `pauli_sum_surface/clifford/cx_alias/{side}` |
| provisional regression | **1.163×** | 1.145–1.180× | 2 | 63.692 ns | 74.065 ns | `word_surface/pattern/parse/{side}/indexed` |
| confirmed regression | **1.160×** | 1.151–1.167× | 7 | 3.900 ns | 4.487 ns | `word_surface/construct/lossy/{side}/new_identity/256` |
| confirmed regression | **1.156×** | 1.148–1.160× | 3 | 4.977 µs | 5.753 µs | `pauli_sum_surface/clifford_batch/cz/{side}` |
| provisional regression | **1.155×** | 1.155–1.155× | 1 | 1.174 ns | 1.356 ns | `word_surface/ordinary/observation/256/{side}/equality` |
| provisional regression | **1.151×** | 1.028–1.274× | 2 | 4.164 µs | 4.809 µs | `sym/surface/propagation/clifford/{side}/sqrt_x` |
| confirmed regression | **1.150×** | 1.111–1.155× | 4 | 478.530 ns | 546.435 ns | `pauli_sum/workload_truncate/{side}/w50/threshold` |
| confirmed regression | **1.149×** | 1.133–1.150× | 3 | 1.429 µs | 1.642 µs | `pauli_sum_surface/clifford/zcz_alias/{side}` |
| confirmed regression | **1.148×** | 1.146–1.169× | 3 | 3.999 µs | 4.642 µs | `tableau-surface/noise/generalized/loss_channel/{side}` |
| confirmed regression | **1.148×** | 1.140–1.150× | 3 | 1.446 µs | 1.663 µs | `pauli_sum_surface/clifford/cz/{side}` |
| provisional regression | **1.145×** | 1.145–1.145× | 1 | 1.525 ns | 1.746 ns | `word_surface/lossy/observation/256/{side}/equality` |
| confirmed regression | **1.143×** | 1.121–1.177× | 3 | 990.980 ns | 1.142 µs | `pauli_sum_surface/noise_batch/depolarize1/{side}` |
| confirmed regression | **1.143×** | 1.141–1.145× | 3 | 58.694 µs | 67.178 µs | `mixture/noise/depolarize2_many/{side}` |
| confirmed regression | **1.140×** | 1.131–1.160× | 7 | 3.939 ns | 4.499 ns | `word_surface/lossy/clifford_present/256/{side}/z` |
| confirmed regression | **1.135×** | 0.968–1.267× | 7 | 3.996 ns | 4.505 ns | `word_surface/lossy/clifford_present/256/{side}/x` |
| provisional regression | **1.134×** | 1.134–1.134× | 1 | 704.040 µs | 798.280 µs | `pauli_sum/workload_trotter_ablation/{side}/full` |
| confirmed regression | **1.134×** | 1.130–1.155× | 3 | 58.661 µs | 66.502 µs | `mixture/noise/two_qubit_pauli_error_many/{side}` |
| confirmed regression | **1.126×** | 1.109–1.157× | 7 | 698.320 ns | 795.270 ns | `mixture/measure/case_b/{side}` |
| confirmed regression | **1.122×** | 0.817–1.132× | 7 | 4.022 ns | 4.487 ns | `word_surface/lossy/clifford_present/256/{side}/y` |
| confirmed regression | **1.120×** | 1.080–1.152× | 4 | 681.085 ns | 763.820 ns | `tableau-micro/msd_measure_single/{side}` |
| confirmed regression | **1.119×** | 1.109–1.127× | 7 | 391.960 ns | 439.280 ns | `tableau-surface/clifford/bare/y_many/{side}` |
| provisional regression | **1.116×** | 1.093–1.139× | 2 | 4.307 µs | 4.801 µs | `sym/surface/propagation/clifford/{side}/cy` |
| confirmed regression | **1.113×** | 1.083–1.146× | 7 | 173.030 ns | 195.110 ns | `pauli_sum_surface/add/extend/{side}` |
| confirmed regression | **1.112×** | 1.107–1.113× | 4 | 3.499 ms | 3.887 ms | `pauli_sum/workload_qubit_sweep/{side}/n12` |
| provisional regression | **1.110×** | 1.025–1.194× | 2 | 3.782 µs | 4.202 µs | `sym/surface/propagation/rotation_two/{side}/rxx` |
| confirmed regression | **1.108×** | 1.102–1.131× | 7 | 316.790 ns | 349.030 ns | `mixture/lifecycle/clone/{side}` |
| confirmed regression | **1.108×** | 1.094–1.117× | 7 | 407.280 ns | 451.060 ns | `tableau-surface/clifford/generalized/y_many/{side}` |
| confirmed regression | **1.105×** | 1.102–1.110× | 7 | 429.710 ns | 475.110 ns | `tableau-surface/clifford/bare/s_dag_many/{side}` |
| confirmed regression | **1.100×** | 1.071–1.105× | 7 | 444.360 ns | 486.070 ns | `tableau-surface/clifford/generalized/s_dag_many/{side}` |
| confirmed regression | **1.098×** | 1.089–1.114× | 7 | 423.040 ns | 463.920 ns | `tableau-surface/clifford/bare/s_many/{side}` |
| provisional regression | **1.098×** | 1.098–1.098× | 1 | 1.361 ns | 1.493 ns | `word_surface/phased/observation/256/{side}/equality` |
| provisional regression | **1.096×** | 0.973–1.219× | 2 | 4.880 µs | 5.319 µs | `sym/surface/propagation/rotation_two/{side}/ryy` |
| confirmed regression | **1.094×** | 1.065–1.129× | 4 | 488.060 ns | 533.920 ns | `pauli_sum/workload_truncate/{side}/w3/threshold` |
| confirmed regression | **1.092×** | 1.090–1.101× | 7 | 435.670 ns | 475.670 ns | `tableau-surface/clifford/generalized/s_many/{side}` |
| confirmed regression | **1.092×** | 1.013–1.178× | 7 | 273.640 ns | 295.870 ns | `pauli_sum_surface/noise/depolarize2/{side}` |
| confirmed regression | **1.089×** | 1.064–1.096× | 3 | 5.541 µs | 6.023 µs | `tableau-surface/noise/generalized/asymmetric_loss_channel/{side}` |
| confirmed regression | **1.084×** | 1.055–1.089× | 7 | 371.750 ns | 404.750 ns | `tableau-surface/clifford/generalized/x_many/{side}` |
| provisional regression | **1.082×** | 1.082–1.082× | 1 | 5.073 ns | 5.487 ns | `word_surface/ordinary/mutate/256/{side}/set_z_bit` |
| provisional regression | **1.081×** | 1.081–1.081× | 1 | 691.310 µs | 747.260 µs | `pauli_sum/integration_trotter_decomposed_rzz/{side}/trotter` |
| confirmed regression | **1.076×** | 1.058–1.082× | 7 | 363.660 ns | 393.650 ns | `tableau-surface/clifford/bare/z_many/{side}` |
| confirmed regression | **1.075×** | 1.010–1.113× | 3 | 232.340 ns | 258.530 ns | `pauli_sum_surface/clifford/z/{side}` |
| provisional regression | **1.075×** | 1.067–1.083× | 2 | 356.335 ns | 383.070 ns | `pauli_sum_surface/inspect/equality_equal_support/{side}` |
| confirmed regression | **1.072×** | 1.034–1.124× | 7 | 378.120 ns | 405.590 ns | `tableau-surface/clifford/generalized/z_many/{side}` |
| confirmed regression | **1.067×** | 1.063–1.071× | 4 | 10.162 ms | 10.845 ms | `pauli_sum/workload_qubit_sweep/{side}/n20` |
| provisional regression | **1.067×** | 1.067–1.067× | 1 | 4.974 ns | 5.305 ns | `word_surface/ordinary/mutate/256/{side}/set_x_bit` |
| provisional regression | **1.065×** | 0.995–1.135× | 2 | 35.542 ns | 37.908 ns | `tableau-surface/projection/project_case_b/{side}` |
| provisional regression | **1.065×** | 1.027–1.103× | 2 | 4.183 µs | 4.461 µs | `sym/surface/propagation/clifford/{side}/alias_zcy` |
| confirmed regression | **1.062×** | 1.043–1.064× | 3 | 3.839 µs | 4.067 µs | `mixture/noise/depolarize2/{side}` |
| confirmed regression | **1.061×** | 1.045–1.074× | 6 | 171.675 ns | 181.385 ns | `sym/surface/eval/{side}/term` |
| confirmed regression | **1.061×** | 0.996–1.063× | 3 | 5.267 µs | 5.477 µs | `tableau-surface/clifford/bare/cz_many/{side}` |
| provisional regression | **1.060×** | 1.057–1.064× | 2 | 137.810 ns | 146.115 ns | `word_surface/ordinary/read/256/{side}/iter_traverse` |
| confirmed regression | **1.059×** | 1.023–1.095× | 7 | 363.820 ns | 389.460 ns | `tableau-surface/clifford/bare/x_many/{side}` |
| confirmed regression | **1.059×** | 1.053–1.068× | 3 | 3.796 µs | 4.055 µs | `mixture/noise/two_qubit_pauli_error/{side}` |
| confirmed regression | **1.055×** | 1.045–1.067× | 3 | 193.490 ns | 206.480 ns | `pauli_sum_surface/algebra/scale/{side}` |
| provisional regression | **1.051×** | 1.044–1.058× | 2 | 139.725 ns | 146.830 ns | `word_surface/phased/read/256/{side}/iter_traverse` |
| confirmed regression | **1.050×** | 1.047–1.053× | 4 | 19.241 ms | 20.203 ms | `pauli_sum/workload_qubit_sweep/{side}/n28` |
| provisional regression | **1.048×** | 1.044–1.052× | 2 | 6.339 ns | 6.643 ns | `sym/surface/operator_add/{side}/term_add_coefficient` |
| confirmed regression | **1.046×** | 1.037–1.056× | 3 | 5.727 µs | 6.022 µs | `tableau-surface/noise/generalized/correlated_loss_channel/{side}` |
| provisional regression | **1.045×** | 1.004–1.085× | 2 | 3.660 µs | 3.831 µs | `sym/surface/propagation/rotation_one/{side}/rz` |
| confirmed regression | **1.038×** | 0.968–1.140× | 7 | 257.580 ns | 260.690 ns | `pauli_sum_surface/truncate/max_loss_weight_active/{side}` |
| confirmed regression | **1.037×** | 1.030–1.044× | 3 | 954.150 ns | 982.640 ns | `pauli_sum/clifford_x/{side}/x` |
| confirmed regression | **1.036×** | 1.033–1.037× | 4 | 65.793 µs | 68.066 µs | `pauli_sum/workload_qubit_sweep/{side}/n4` |
| confirmed regression | **1.035×** | 1.018–1.077× | 7 | 8.068 ns | 8.363 ns | `word_surface/ordinary/clifford/256/{side}/cy` |
| provisional regression | **1.035×** | 1.033–1.037× | 2 | 4.077 µs | 4.221 µs | `sym/surface/propagation/clifford/{side}/s` |
| confirmed regression | **1.033×** | 1.027–1.048× | 3 | 0.796 ns | 0.822 ns | `word_surface/ordinary/read/256/{side}/weight` |
| provisional regression | **1.032×** | 0.994–1.071× | 2 | 3.674 µs | 3.798 µs | `sym/surface/propagation/rotation_two/{side}/rxz` |
| confirmed regression | **1.031×** | 0.998–1.083× | 7 | 8.102 ns | 8.369 ns | `word_surface/ordinary/clifford/256/{side}/zcy_alias` |
| confirmed regression | **1.031×** | 1.029–1.032× | 4 | 45.846 ms | 47.232 ms | `pauli_sum/workload_qubit_sweep/{side}/n44` |

## Attribution summary

- **Lossy branch keys (1.6–1.8×):** new clone-then-toggle copies three
  atomic caches and performs guarded invalidation; old copies one plain
  digest and uses unchecked setters. The benchmark excludes the later hash.
- **Disabled truncation (1.5–1.7×):** both sides return immediately;
  the roughly 1 ns absolute delta is confirmed but remains unattributed.
- **Pauli-sum construction/clone:** new preallocates persistent scratch and
  clones atomic key caches; this explains much of the 1.27–1.34× delta.
- **Pauli-sum Clifford families and qubit scaling:** drift is localized to
  repeated bijective re-keying. Existing controls exclude reserve and show
  only a small drain/clone contribution; the remaining mechanism is unknown.
- **Lossy-sum integration:** loss, reset, Clifford and rotation stages remain
  slower; the component-cache representation is a source difference, but no
  controlled ablation proves it is the cause.
- **Symbolic propagation/trace:** matched fixtures confirm the regression,
  but the remaining engine/trace mechanism is unattributed. Symbolic eval
  is slower on one-use variables because the new angle cache initializes 32
  entries and computes both sine and cosine on each miss.
- **Tiny tableau observation helpers (2.3–4.0×):** implementations are
  effectively identical and absolute deltas are 3–6 ns; this is attributed
  to inlining/code placement or struct-layout effects, not algorithmic work.
- **Tableau many-gate batches (about 1.06–1.12×):** new adds lazy-hash
  invalidation and uses a different inner loop shape; full scaling workloads
  are otherwise parity or faster.
- **Mixture clone/measurement:** new persistently clones fingerprint buckets
  and measurement eagerly dirties/rebuilds fingerprints. Two-qubit mixture
  noise repeatedly scans rows for each of 15 Pauli branches.

## Complete comparison table

| status | ratio | runs | old | new | benchmark |
|---|---:|---:|---:|---:|---|
| parity | 1.008× | 1 | 1.716 ns | 1.730 ns | `lossy_pauli_word/cnot/{side}/cnot` |
| parity | 1.011× | 1 | 0.540 ns | 0.545 ns | `lossy_pauli_word/key_hash/{side}/warm` |
| parity | 1.007× | 1 | 0.647 ns | 0.651 ns | `lossy_pauli_word/loss_weight/{side}/loss_weight` |
| parity | 1.008× | 1 | 0.746 ns | 0.752 ns | `lossy_pauli_word/weight/{side}/weight` |
| improvement | 0.817× | 3 | 314.480 ns | 256.930 ns | `mixture/gate/cnot_many/{side}` |
| improvement | 0.860× | 3 | 356.290 ns | 306.350 ns | `mixture/gate/cy_many/{side}` |
| improvement | 0.831× | 3 | 314.560 ns | 262.000 ns | `mixture/gate/cz_many/{side}` |
| improvement | 0.771× | 3 | 228.620 ns | 176.340 ns | `mixture/gate/s_dag_many/{side}` |
| improvement | 0.549× | 3 | 116.810 ns | 64.059 ns | `mixture/gate/sqrt_x_dag/{side}` |
| improvement | 0.749× | 3 | 221.300 ns | 164.470 ns | `mixture/gate/sqrt_x_dag_many/{side}` |
| improvement | 0.756× | 3 | 234.250 ns | 177.020 ns | `mixture/gate/sqrt_x_many/{side}` |
| improvement | 0.579× | 3 | 132.270 ns | 76.305 ns | `mixture/gate/sqrt_y_dag/{side}` |
| improvement | 0.815× | 3 | 289.680 ns | 237.240 ns | `mixture/gate/sqrt_y_dag_many/{side}` |
| improvement | 0.810× | 3 | 291.530 ns | 236.050 ns | `mixture/gate/sqrt_y_many/{side}` |
| improvement | 0.884× | 2 | 414.635 µs | 366.595 µs | `mixture/integration/noisy_build/{side}` |
| confirmed regression | 1.108× | 7 | 316.790 ns | 349.030 ns | `mixture/lifecycle/clone/{side}` |
| parity | 1.003× | 3 | 0.415 ns | 0.416 ns | `mixture/lifecycle/is_empty/{side}` |
| parity | 1.001× | 3 | 0.419 ns | 0.418 ns | `mixture/lifecycle/len/{side}` |
| improvement | 0.738× | 3 | 327.080 ns | 239.600 ns | `mixture/lifecycle/new/{side}` |
| improvement | 0.712× | 3 | 315.690 ns | 224.880 ns | `mixture/lifecycle/new_with_seed/{side}` |
| improvement | 0.679× | 3 | 24.733 ns | 17.477 ns | `mixture/lifecycle/normalize_probabilities/{side}` |
| improvement | 0.342× | 3 | 97.868 ns | 33.497 ns | `mixture/lifecycle/truncate_cutoff/{side}` |
| improvement | 0.715× | 3 | 2.035 µs | 1.444 µs | `mixture/measure/case_a/{side}` |
| confirmed regression | 1.126× | 7 | 698.320 ns | 795.270 ns | `mixture/measure/case_b/{side}` |
| confirmed regression | 1.329× | 7 | 211.310 ns | 285.450 ns | `mixture/measure/lost/{side}` |
| improvement | 0.656× | 3 | 923.720 ns | 605.710 ns | `mixture/noise/correlated_loss_channel/{side}` |
| improvement | 0.694× | 3 | 1.039 µs | 723.950 ns | `mixture/noise/depolarize1/{side}` |
| improvement | 0.799× | 3 | 15.366 µs | 12.376 µs | `mixture/noise/depolarize1_many/{side}` |
| confirmed regression | 1.062× | 3 | 3.839 µs | 4.067 µs | `mixture/noise/depolarize2/{side}` |
| confirmed regression | 1.143× | 3 | 58.694 µs | 67.178 µs | `mixture/noise/depolarize2_many/{side}` |
| improvement | 0.622× | 3 | 372.980 ns | 231.070 ns | `mixture/noise/loss_channel/{side}` |
| improvement | 0.701× | 3 | 1.023 µs | 720.980 ns | `mixture/noise/pauli_error/{side}` |
| improvement | 0.805× | 3 | 15.395 µs | 12.398 µs | `mixture/noise/pauli_error_many/{side}` |
| improvement | 0.917× | 3 | 314.180 ns | 287.710 ns | `mixture/noise/reset_loss_channel/{side}` |
| confirmed regression | 1.059× | 3 | 3.796 µs | 4.055 µs | `mixture/noise/two_qubit_pauli_error/{side}` |
| confirmed regression | 1.134× | 3 | 58.661 µs | 66.502 µs | `mixture/noise/two_qubit_pauli_error_many/{side}` |
| improvement | 0.625× | 3 | 511.230 ns | 323.000 ns | `mixture/noise/x_error/{side}` |
| improvement | 0.720× | 3 | 2.661 µs | 1.916 µs | `mixture/noise/x_error_many/{side}` |
| improvement | 0.630× | 3 | 518.080 ns | 325.440 ns | `mixture/noise/y_error/{side}` |
| improvement | 0.732× | 3 | 2.666 µs | 1.962 µs | `mixture/noise/y_error_many/{side}` |
| improvement | 0.633× | 3 | 513.240 ns | 325.130 ns | `mixture/noise/z_error/{side}` |
| improvement | 0.742× | 3 | 2.688 µs | 1.998 µs | `mixture/noise/z_error_many/{side}` |
| improvement | 0.821× | 3 | 3.633 µs | 2.984 µs | `mixture/reset/reset_many/{side}` |
| improvement | 0.839× | 3 | 765.790 ns | 635.770 ns | `mixture/reset/reset_x/{side}` |
| improvement | 0.740× | 3 | 4.282 µs | 3.148 µs | `mixture/reset/reset_x_many/{side}` |
| improvement | 0.816× | 3 | 814.590 ns | 670.130 ns | `mixture/reset/reset_y/{side}` |
| improvement | 0.758× | 3 | 4.318 µs | 3.257 µs | `mixture/reset/reset_y_many/{side}` |
| improvement | 0.830× | 3 | 703.530 ns | 582.930 ns | `mixture/reset/reset_z/{side}` |
| improvement | 0.817× | 3 | 3.691 µs | 2.971 µs | `mixture/reset/reset_z_many/{side}` |
| improvement | 0.881× | 3 | 348.220 ns | 309.920 ns | `mixture/rotation/rotate_1/{side}` |
| improvement | 0.769× | 3 | 704.580 ns | 536.200 ns | `mixture/rotation/rotate_2/{side}` |
| improvement | 0.946× | 3 | 1.037 µs | 983.470 ns | `mixture/rotation/rx_many/{side}` |
| improvement | 0.784× | 3 | 691.930 ns | 539.980 ns | `mixture/rotation/rxx/{side}` |
| improvement | 0.515× | 3 | 3.412 µs | 1.757 µs | `mixture/rotation/rxx_many/{side}` |
| improvement | 0.778× | 3 | 705.870 ns | 549.800 ns | `mixture/rotation/rxy/{side}` |
| improvement | 0.520× | 3 | 3.451 µs | 1.793 µs | `mixture/rotation/rxy_many/{side}` |
| improvement | 0.764× | 3 | 706.110 ns | 539.590 ns | `mixture/rotation/rxz/{side}` |
| improvement | 0.515× | 3 | 3.411 µs | 1.754 µs | `mixture/rotation/rxz_many/{side}` |
| improvement | 0.962× | 3 | 1.067 µs | 1.016 µs | `mixture/rotation/ry_many/{side}` |
| improvement | 0.772× | 3 | 705.740 ns | 548.060 ns | `mixture/rotation/ryx/{side}` |
| improvement | 0.517× | 3 | 3.458 µs | 1.783 µs | `mixture/rotation/ryx_many/{side}` |
| improvement | 0.777× | 3 | 717.230 ns | 557.170 ns | `mixture/rotation/ryy/{side}` |
| improvement | 0.519× | 3 | 3.481 µs | 1.800 µs | `mixture/rotation/ryy_many/{side}` |
| improvement | 0.775× | 3 | 707.590 ns | 547.910 ns | `mixture/rotation/ryz/{side}` |
| improvement | 0.514× | 3 | 3.449 µs | 1.773 µs | `mixture/rotation/ryz_many/{side}` |
| improvement | 0.941× | 3 | 928.680 ns | 873.930 ns | `mixture/rotation/rz_many/{side}` |
| improvement | 0.776× | 3 | 697.940 ns | 541.260 ns | `mixture/rotation/rzx/{side}` |
| improvement | 0.516× | 3 | 3.422 µs | 1.760 µs | `mixture/rotation/rzx_many/{side}` |
| improvement | 0.772× | 3 | 708.770 ns | 550.730 ns | `mixture/rotation/rzy/{side}` |
| improvement | 0.514× | 3 | 3.434 µs | 1.766 µs | `mixture/rotation/rzy_many/{side}` |
| improvement | 0.742× | 3 | 724.970 ns | 533.620 ns | `mixture/rotation/rzz/{side}` |
| improvement | 0.737× | 3 | 2.169 µs | 1.606 µs | `mixture/rotation/rzz_many/{side}` |
| improvement | 0.877× | 3 | 344.570 ns | 300.160 ns | `mixture/rotation/t_dag/{side}` |
| improvement | 0.943× | 3 | 917.940 ns | 866.070 ns | `mixture/rotation/t_dag_many/{side}` |
| improvement | 0.937× | 3 | 929.970 ns | 872.620 ns | `mixture/rotation/t_many/{side}` |
| improvement | 0.883× | 3 | 39.527 µs | 34.993 µs | `mixture/sampler/adaptive_sample_shots_128/{side}` |
| improvement | 0.745× | 3 | 803.190 ns | 598.180 ns | `mixture/sampler/construction/{side}` |
| improvement | 0.970× | 1 | 122.690 µs | 118.970 µs | `mixture/sampler/parallel_branch_shot_scaling/{side}/1_branches/1024_shots` |
| improvement | 0.967× | 1 | 61.020 µs | 59.036 µs | `mixture/sampler/parallel_branch_shot_scaling/{side}/1_branches/128_shots` |
| parity | 1.000× | 1 | 35.631 µs | 35.635 µs | `mixture/sampler/parallel_branch_shot_scaling/{side}/1_branches/16_shots` |
| improvement | 0.933× | 1 | 393.460 ns | 367.240 ns | `mixture/sampler/parallel_branch_shot_scaling/{side}/1_branches/1_shots` |
| parity | 0.978× | 1 | 130.730 µs | 127.900 µs | `mixture/sampler/parallel_branch_shot_scaling/{side}/64_branches/1024_shots` |
| improvement | 0.958× | 1 | 64.683 µs | 61.949 µs | `mixture/sampler/parallel_branch_shot_scaling/{side}/64_branches/128_shots` |
| improvement | 0.965× | 1 | 40.207 µs | 38.783 µs | `mixture/sampler/parallel_branch_shot_scaling/{side}/64_branches/16_shots` |
| improvement | 0.612× | 1 | 3.460 µs | 2.116 µs | `mixture/sampler/parallel_branch_shot_scaling/{side}/64_branches/1_shots` |
| improvement | 0.959× | 1 | 123.280 µs | 118.280 µs | `mixture/sampler/parallel_branch_shot_scaling/{side}/8_branches/1024_shots` |
| parity | 1.003× | 1 | 59.123 µs | 59.289 µs | `mixture/sampler/parallel_branch_shot_scaling/{side}/8_branches/128_shots` |
| parity | 0.983× | 1 | 36.284 µs | 35.668 µs | `mixture/sampler/parallel_branch_shot_scaling/{side}/8_branches/16_shots` |
| improvement | 0.861× | 1 | 641.980 ns | 553.000 ns | `mixture/sampler/parallel_branch_shot_scaling/{side}/8_branches/1_shots` |
| improvement | 0.846× | 1 | 660.910 ns | 559.270 ns | `mixture/sampler/serial_parallel_crossover/{side}_parallel/1` |
| parity | 0.988× | 1 | 121.580 µs | 120.090 µs | `mixture/sampler/serial_parallel_crossover/{side}_parallel/1024` |
| improvement | 0.885× | 1 | 66.834 µs | 59.135 µs | `mixture/sampler/serial_parallel_crossover/{side}_parallel/128` |
| parity | 0.976× | 1 | 37.202 µs | 36.303 µs | `mixture/sampler/serial_parallel_crossover/{side}_parallel/16` |
| improvement | 0.957× | 1 | 23.904 µs | 22.875 µs | `mixture/sampler/serial_parallel_crossover/{side}_parallel/2` |
| parity | 0.972× | 1 | 74.490 µs | 72.401 µs | `mixture/sampler/serial_parallel_crossover/{side}_parallel/256` |
| improvement | 0.919× | 1 | 44.825 µs | 41.184 µs | `mixture/sampler/serial_parallel_crossover/{side}_parallel/32` |
| parity | 0.982× | 1 | 32.418 µs | 31.823 µs | `mixture/sampler/serial_parallel_crossover/{side}_parallel/4` |
| parity | 0.992× | 1 | 91.375 µs | 90.628 µs | `mixture/sampler/serial_parallel_crossover/{side}_parallel/512` |
| parity | 0.977× | 1 | 50.294 µs | 49.131 µs | `mixture/sampler/serial_parallel_crossover/{side}_parallel/64` |
| parity | 0.991× | 1 | 33.970 µs | 33.666 µs | `mixture/sampler/serial_parallel_crossover/{side}_parallel/8` |
| improvement | 0.827× | 1 | 648.520 ns | 536.580 ns | `mixture/sampler/serial_parallel_crossover/{side}_serial/1` |
| improvement | 0.916× | 1 | 314.470 µs | 288.080 µs | `mixture/sampler/serial_parallel_crossover/{side}_serial/1024` |
| improvement | 0.899× | 1 | 40.103 µs | 36.049 µs | `mixture/sampler/serial_parallel_crossover/{side}_serial/128` |
| improvement | 0.905× | 1 | 5.290 µs | 4.786 µs | `mixture/sampler/serial_parallel_crossover/{side}_serial/16` |
| improvement | 0.863× | 1 | 966.080 ns | 833.330 ns | `mixture/sampler/serial_parallel_crossover/{side}_serial/2` |
| improvement | 0.914× | 1 | 78.845 µs | 72.089 µs | `mixture/sampler/serial_parallel_crossover/{side}_serial/256` |
| improvement | 0.895× | 1 | 10.392 µs | 9.296 µs | `mixture/sampler/serial_parallel_crossover/{side}_serial/32` |
| improvement | 0.881× | 1 | 1.595 µs | 1.405 µs | `mixture/sampler/serial_parallel_crossover/{side}_serial/4` |
| improvement | 0.906× | 1 | 158.850 µs | 143.940 µs | `mixture/sampler/serial_parallel_crossover/{side}_serial/512` |
| improvement | 0.898× | 1 | 20.309 µs | 18.228 µs | `mixture/sampler/serial_parallel_crossover/{side}_serial/64` |
| improvement | 0.891× | 1 | 2.850 µs | 2.539 µs | `mixture/sampler/serial_parallel_crossover/{side}_serial/8` |
| improvement | 0.824× | 3 | 625.540 ns | 515.590 ns | `mixture/sampler/serial_shot_scaling/{side}/1` |
| improvement | 0.880× | 3 | 313.230 µs | 274.210 µs | `mixture/sampler/serial_shot_scaling/{side}/1024` |
| improvement | 0.873× | 3 | 39.506 µs | 34.809 µs | `mixture/sampler/serial_shot_scaling/{side}/128` |
| improvement | 0.875× | 3 | 5.262 µs | 4.611 µs | `mixture/sampler/serial_shot_scaling/{side}/16` |
| improvement | 0.818× | 3 | 625.950 ns | 507.930 ns | `mixture/sampler/single_sample/{side}` |
| improvement | 0.546× | 2 | 42.124 ns | 23.002 ns | `mixture/scaling/branches_h/{side}/1` |
| improvement | 0.495× | 2 | 531.520 ns | 263.005 ns | `mixture/scaling/branches_h/{side}/16` |
| improvement | 0.532× | 2 | 75.228 ns | 39.984 ns | `mixture/scaling/branches_h/{side}/2` |
| improvement | 0.529× | 2 | 141.610 ns | 74.846 ns | `mixture/scaling/branches_h/{side}/4` |
| improvement | 0.511× | 2 | 273.255 ns | 139.615 ns | `mixture/scaling/branches_h/{side}/8` |
| improvement | 0.546× | 2 | 240.455 ns | 131.295 ns | `mixture/scaling/branches_loss/{side}/1` |
| improvement | 0.590× | 2 | 2.055 µs | 1.212 µs | `mixture/scaling/branches_loss/{side}/16` |
| improvement | 0.597× | 2 | 389.450 ns | 232.685 ns | `mixture/scaling/branches_loss/{side}/2` |
| improvement | 0.599× | 2 | 658.475 ns | 394.420 ns | `mixture/scaling/branches_loss/{side}/4` |
| improvement | 0.574× | 2 | 1.135 µs | 651.430 ns | `mixture/scaling/branches_loss/{side}/8` |
| improvement | 0.858× | 2 | 130.490 ns | 111.970 ns | `mixture/scaling/branches_sampler_construction/{side}/1` |
| improvement | 0.718× | 2 | 1.625 µs | 1.166 µs | `mixture/scaling/branches_sampler_construction/{side}/16` |
| improvement | 0.846× | 2 | 227.205 ns | 192.120 ns | `mixture/scaling/branches_sampler_construction/{side}/2` |
| improvement | 0.810× | 2 | 423.925 ns | 343.340 ns | `mixture/scaling/branches_sampler_construction/{side}/4` |
| improvement | 0.729× | 2 | 828.700 ns | 603.615 ns | `mixture/scaling/branches_sampler_construction/{side}/8` |
| improvement | 0.623× | 3 | 5.839 µs | 3.629 µs | `pauli_sum/build_batch/{side}/build_add_assign` |
| confirmed regression | 1.260× | 3 | 5.206 µs | 6.589 µs | `pauli_sum/clifford_cnot/{side}/cnot` |
| confirmed regression | 1.271× | 3 | 5.176 µs | 6.581 µs | `pauli_sum/clifford_h/{side}/h` |
| confirmed regression | 1.037× | 3 | 954.150 ns | 982.640 ns | `pauli_sum/clifford_x/{side}/x` |
| improvement | 0.779× | 1 | 327.340 µs | 255.140 µs | `pauli_sum/integration_trotter/{side}/trotter` |
| provisional regression | 1.081× | 1 | 691.310 µs | 747.260 µs | `pauli_sum/integration_trotter_decomposed_rzz/{side}/trotter` |
| confirmed regression | 1.414× | 3 | 300.820 µs | 425.420 µs | `pauli_sum/loss_attrib/clifford/{side}` |
| improvement | 0.694× | 3 | 50.610 ms | 35.602 ms | `pauli_sum/loss_attrib/correlated/{side}` |
| confirmed regression | 1.292× | 7 | 16.145 ms | 20.861 ms | `pauli_sum/loss_attrib/loss/{side}` |
| confirmed regression | 1.207× | 3 | 35.660 µs | 43.029 µs | `pauli_sum/loss_attrib/reset/{side}` |
| confirmed regression | 1.225× | 7 | 16.945 ms | 20.751 ms | `pauli_sum/loss_attrib/rotation/{side}` |
| improvement | 0.936× | 3 | 236.320 ms | 221.100 ms | `pauli_sum/loss_interleaved_n12/{side}` |
| improvement | 0.800× | 1 | 4.508 µs | 3.606 µs | `pauli_sum/multiply_word/{side}/mul_word` |
| improvement | 0.001× | 3 | 1.101 ms | 1.520 µs | `pauli_sum/overlap/{side}/overlap` |
| improvement | 0.920× | 3 | 1.437 µs | 1.326 µs | `pauli_sum/pauli_error/{side}/pauli_error` |
| improvement | 0.840× | 1 | 5.299 µs | 4.449 µs | `pauli_sum/pauli_error/{side}/pauli_error_sweep` |
| parity | 1.009× | 1 | 52.965 µs | 53.436 µs | `pauli_sum/rekey_cnot/{side}/cnot_sweep` |
| improvement | 0.892× | 3 | 5.862 µs | 5.259 µs | `pauli_sum/rotation_rx/{side}/rx` |
| improvement | 0.941× | 1 | 6.021 µs | 5.664 µs | `pauli_sum/rx/{side}/rx_sweep` |
| improvement | 0.957× | 1 | 11.504 µs | 11.005 µs | `pauli_sum/rzz/{side}/rzz_sweep` |
| parity | 1.008× | 3 | 500.810 ns | 508.400 ns | `pauli_sum/scale/{side}/scale` |
| parity | 0.989× | 1 | 239.620 ns | 237.030 ns | `pauli_sum/truncate/{side}/truncate` |
| parity | 1.012× | 1 | 286.130 ns | 289.620 ns | `pauli_sum/truncate_active/{side}/truncate` |
| confirmed regression | 1.112× | 4 | 3.499 ms | 3.887 ms | `pauli_sum/workload_qubit_sweep/{side}/n12` |
| confirmed regression | 1.067× | 4 | 10.162 ms | 10.845 ms | `pauli_sum/workload_qubit_sweep/{side}/n20` |
| confirmed regression | 1.050× | 4 | 19.241 ms | 20.203 ms | `pauli_sum/workload_qubit_sweep/{side}/n28` |
| parity | 1.028× | 4 | 30.957 ms | 31.841 ms | `pauli_sum/workload_qubit_sweep/{side}/n36` |
| confirmed regression | 1.036× | 4 | 65.793 µs | 68.066 µs | `pauli_sum/workload_qubit_sweep/{side}/n4` |
| confirmed regression | 1.031× | 4 | 45.846 ms | 47.232 ms | `pauli_sum/workload_qubit_sweep/{side}/n44` |
| improvement | 0.900× | 1 | 11.493 ms | 10.340 ms | `pauli_sum/workload_random_circuit/{side}/circuit` |
| provisional regression | 1.134× | 1 | 704.040 µs | 798.280 µs | `pauli_sum/workload_trotter_ablation/{side}/full` |
| parity | 1.000× | 1 | 39.561 µs | 39.564 µs | `pauli_sum/workload_trotter_ablation/{side}/no_rekey` |
| improvement | 0.910× | 4 | 1.503 µs | 1.368 µs | `pauli_sum/workload_truncate/{side}/w120/combined` |
| improvement | 0.899× | 4 | 1.103 µs | 991.140 ns | `pauli_sum/workload_truncate/{side}/w120/cut10` |
| parity | 0.995× | 4 | 684.640 ns | 682.810 ns | `pauli_sum/workload_truncate/{side}/w120/cut1000` |
| confirmed regression | 1.533× | 4 | 1.765 ns | 2.659 ns | `pauli_sum/workload_truncate/{side}/w120/max_sentinel` |
| confirmed regression | 1.205× | 4 | 487.960 ns | 589.720 ns | `pauli_sum/workload_truncate/{side}/w120/threshold` |
| improvement | 0.968× | 4 | 1.067 µs | 1.033 µs | `pauli_sum/workload_truncate/{side}/w3/combined` |
| parity | 1.012× | 4 | 676.320 ns | 681.305 ns | `pauli_sum/workload_truncate/{side}/w3/cut10` |
| parity | 1.001× | 4 | 678.725 ns | 679.135 ns | `pauli_sum/workload_truncate/{side}/w3/cut1000` |
| confirmed regression | 1.533× | 4 | 1.740 ns | 2.675 ns | `pauli_sum/workload_truncate/{side}/w3/max_sentinel` |
| confirmed regression | 1.094× | 4 | 488.060 ns | 533.920 ns | `pauli_sum/workload_truncate/{side}/w3/threshold` |
| improvement | 0.903× | 4 | 1.477 µs | 1.334 µs | `pauli_sum/workload_truncate/{side}/w50/combined` |
| improvement | 0.903× | 4 | 1.112 µs | 1.006 µs | `pauli_sum/workload_truncate/{side}/w50/cut10` |
| parity | 1.030× | 4 | 681.365 ns | 698.915 ns | `pauli_sum/workload_truncate/{side}/w50/cut1000` |
| confirmed regression | 1.519× | 4 | 1.758 ns | 2.686 ns | `pauli_sum/workload_truncate/{side}/w50/max_sentinel` |
| confirmed regression | 1.150× | 4 | 478.530 ns | 546.435 ns | `pauli_sum/workload_truncate/{side}/w50/threshold` |
| improvement | 0.558× | 1 | 56.789 µs | 31.714 µs | `pauli_sum_indexmap/build/{side}` |
| improvement | 0.880× | 1 | 8.256 µs | 7.268 µs | `pauli_sum_indexmap/gates/{side}/cnot` |
| improvement | 0.924× | 1 | 6.997 µs | 6.468 µs | `pauli_sum_indexmap/gates/{side}/rx` |
| improvement | 0.689× | 1 | 100.460 µs | 69.231 µs | `pauli_sum_indexmap/ordered_terms/{side}` |
| confirmed regression | 1.113× | 7 | 173.030 ns | 195.110 ns | `pauli_sum_surface/add/extend/{side}` |
| improvement | 0.611× | 3 | 608.020 ns | 375.140 ns | `pauli_sum_surface/add/sum_disjoint/{side}` |
| confirmed regression | 1.174× | 7 | 7.124 ns | 8.123 ns | `pauli_sum_surface/add/term/{side}` |
| improvement | 0.527× | 3 | 3.464 µs | 1.826 µs | `pauli_sum_surface/algebra/mul_word/{side}` |
| improvement | 0.010× | 3 | 41.569 µs | 408.460 ns | `pauli_sum_surface/algebra/overlap/{side}` |
| confirmed regression | 1.055× | 3 | 193.490 ns | 206.480 ns | `pauli_sum_surface/algebra/scale/{side}` |
| confirmed regression | 1.172× | 3 | 1.410 µs | 1.655 µs | `pauli_sum_surface/clifford/cnot/{side}` |
| confirmed regression | 1.163× | 3 | 1.397 µs | 1.631 µs | `pauli_sum_surface/clifford/cx_alias/{side}` |
| confirmed regression | 1.209× | 3 | 1.528 µs | 1.854 µs | `pauli_sum_surface/clifford/cy/{side}` |
| confirmed regression | 1.148× | 3 | 1.446 µs | 1.663 µs | `pauli_sum_surface/clifford/cz/{side}` |
| confirmed regression | 1.191× | 3 | 1.307 µs | 1.567 µs | `pauli_sum_surface/clifford/h/{side}` |
| confirmed regression | 1.235× | 3 | 1.230 µs | 1.522 µs | `pauli_sum_surface/clifford/s/{side}` |
| confirmed regression | 1.235× | 3 | 1.222 µs | 1.510 µs | `pauli_sum_surface/clifford/s_dag/{side}` |
| confirmed regression | 1.195× | 3 | 1.274 µs | 1.501 µs | `pauli_sum_surface/clifford/sqrt_x/{side}` |
| confirmed regression | 1.203× | 3 | 1.246 µs | 1.500 µs | `pauli_sum_surface/clifford/sqrt_x_dag/{side}` |
| confirmed regression | 1.183× | 3 | 1.319 µs | 1.553 µs | `pauli_sum_surface/clifford/sqrt_y/{side}` |
| confirmed regression | 1.193× | 3 | 1.307 µs | 1.559 µs | `pauli_sum_surface/clifford/sqrt_y_dag/{side}` |
| parity | 1.022× | 3 | 253.180 ns | 257.560 ns | `pauli_sum_surface/clifford/x/{side}` |
| parity | 1.010× | 3 | 257.620 ns | 260.080 ns | `pauli_sum_surface/clifford/y/{side}` |
| confirmed regression | 1.075× | 3 | 232.340 ns | 258.530 ns | `pauli_sum_surface/clifford/z/{side}` |
| confirmed regression | 1.169× | 3 | 1.401 µs | 1.638 µs | `pauli_sum_surface/clifford/zcx_alias/{side}` |
| confirmed regression | 1.221× | 3 | 1.519 µs | 1.854 µs | `pauli_sum_surface/clifford/zcy_alias/{side}` |
| confirmed regression | 1.149× | 3 | 1.429 µs | 1.642 µs | `pauli_sum_surface/clifford/zcz_alias/{side}` |
| confirmed regression | 1.198× | 3 | 4.775 µs | 5.755 µs | `pauli_sum_surface/clifford_batch/cnot/{side}` |
| confirmed regression | 1.276× | 3 | 5.216 µs | 6.708 µs | `pauli_sum_surface/clifford_batch/cy/{side}` |
| confirmed regression | 1.156× | 3 | 4.977 µs | 5.753 µs | `pauli_sum_surface/clifford_batch/cz/{side}` |
| confirmed regression | 1.237× | 3 | 4.314 µs | 5.343 µs | `pauli_sum_surface/clifford_batch/h/{side}` |
| confirmed regression | 1.290× | 3 | 4.031 µs | 5.205 µs | `pauli_sum_surface/clifford_batch/s/{side}` |
| confirmed regression | 1.291× | 3 | 4.063 µs | 5.253 µs | `pauli_sum_surface/clifford_batch/s_dag/{side}` |
| confirmed regression | 1.264× | 3 | 4.167 µs | 5.275 µs | `pauli_sum_surface/clifford_batch/sqrt_x/{side}` |
| confirmed regression | 1.256× | 3 | 4.154 µs | 5.234 µs | `pauli_sum_surface/clifford_batch/sqrt_x_dag/{side}` |
| confirmed regression | 1.235× | 3 | 4.327 µs | 5.346 µs | `pauli_sum_surface/clifford_batch/sqrt_y/{side}` |
| confirmed regression | 1.233× | 3 | 4.329 µs | 5.367 µs | `pauli_sum_surface/clifford_batch/sqrt_y_dag/{side}` |
| parity | 0.993× | 3 | 1.091 µs | 1.079 µs | `pauli_sum_surface/clifford_batch/x/{side}` |
| parity | 0.989× | 3 | 1.113 µs | 1.088 µs | `pauli_sum_surface/clifford_batch/y/{side}` |
| parity | 0.995× | 3 | 1.101 µs | 1.101 µs | `pauli_sum_surface/clifford_batch/z/{side}` |
| provisional regression | 1.284× | 1 | 376.950 ns | 484.110 ns | `pauli_sum_surface/compare/abs_diff_eq_equal/{side}` |
| provisional regression | 1.305× | 1 | 374.110 ns | 488.040 ns | `pauli_sum_surface/compare/abs_diff_eq_near/{side}` |
| improvement | 0.935× | 1 | 443.670 ns | 414.900 ns | `pauli_sum_surface/compare/relative_eq_equal/{side}` |
| improvement | 0.964× | 1 | 441.920 ns | 426.010 ns | `pauli_sum_surface/compare/relative_eq_near/{side}` |
| improvement | 0.619× | 3 | 1.115 µs | 687.960 ns | `pauli_sum_surface/construct/build_support/{side}` |
| confirmed regression | 1.272× | 7 | 417.320 ns | 532.110 ns | `pauli_sum_surface/construct/clone/{side}` |
| confirmed regression | 1.339× | 7 | 211.620 ns | 283.030 ns | `pauli_sum_surface/construct/empty/{side}` |
| improvement | 0.320× | 3 | 21.009 ns | 6.726 ns | `pauli_sum_surface/construct/parse_word/{side}` |
| improvement | 0.900× | 1 | 1.413 µs | 1.272 µs | `pauli_sum_surface/format/debug/{side}` |
| improvement | 0.872× | 1 | 1.393 µs | 1.215 µs | `pauli_sum_surface/format/display/{side}` |
| parity | 0.985× | 3 | 147.970 ns | 145.750 ns | `pauli_sum_surface/inspect/borrowed_traversal/{side}` |
| confirmed regression | 1.265× | 3 | 0.973 ns | 1.248 ns | `pauli_sum_surface/inspect/contains_key/{side}` |
| provisional regression | 1.173× | 2 | 1.471 ns | 1.726 ns | `pauli_sum_surface/inspect/contains_key_value/{side}` |
| provisional regression | 1.075× | 2 | 356.335 ns | 383.070 ns | `pauli_sum_surface/inspect/equality_equal_support/{side}` |
| confirmed regression | 1.374× | 7 | 1.134 ns | 1.579 ns | `pauli_sum_surface/inspect/get/{side}` |
| parity | 1.001× | 3 | 1.176 ns | 1.177 ns | `pauli_sum_surface/inspect/metadata/{side}` |
| improvement | 0.731× | 3 | 627.650 ns | 458.820 ns | `pauli_sum_surface/loss/channel/{side}` |
| improvement | 0.416× | 3 | 2.205 µs | 917.060 ns | `pauli_sum_surface/loss/correlated/{side}` |
| improvement | 0.634× | 3 | 1.077 µs | 684.150 ns | `pauli_sum_surface/loss/reset/{side}` |
| improvement | 0.673× | 3 | 836.570 ns | 568.410 ns | `pauli_sum_surface/noise/amplitude_damping/{side}` |
| parity | 1.023× | 3 | 238.780 ns | 242.880 ns | `pauli_sum_surface/noise/depolarize1/{side}` |
| confirmed regression | 1.092× | 7 | 273.640 ns | 295.870 ns | `pauli_sum_surface/noise/depolarize2/{side}` |
| improvement | 0.939× | 3 | 327.330 ns | 307.320 ns | `pauli_sum_surface/noise/pauli_error/{side}` |
| improvement | 0.666× | 3 | 685.200 ns | 460.270 ns | `pauli_sum_surface/noise/two_qubit_pauli_error/{side}` |
| improvement | 0.962× | 3 | 328.520 ns | 316.310 ns | `pauli_sum_surface/noise/x_error/{side}` |
| improvement | 0.965× | 3 | 336.960 ns | 323.560 ns | `pauli_sum_surface/noise/y_error/{side}` |
| improvement | 0.932× | 3 | 332.510 ns | 309.900 ns | `pauli_sum_surface/noise/z_error/{side}` |
| confirmed regression | 1.143× | 3 | 990.980 ns | 1.142 µs | `pauli_sum_surface/noise_batch/depolarize1/{side}` |
| confirmed regression | 1.188× | 3 | 1.049 µs | 1.217 µs | `pauli_sum_surface/noise_batch/depolarize2/{side}` |
| parity | 0.999× | 3 | 1.474 µs | 1.509 µs | `pauli_sum_surface/noise_batch/pauli_error/{side}` |
| improvement | 0.850× | 3 | 2.994 µs | 2.565 µs | `pauli_sum_surface/noise_batch/two_qubit_pauli_error/{side}` |
| parity | 1.025× | 3 | 1.472 µs | 1.489 µs | `pauli_sum_surface/noise_batch/x_error/{side}` |
| parity | 1.008× | 3 | 1.497 µs | 1.514 µs | `pauli_sum_surface/noise_batch/y_error/{side}` |
| parity | 1.026× | 3 | 1.441 µs | 1.478 µs | `pauli_sum_surface/noise_batch/z_error/{side}` |
| improvement | 0.787× | 3 | 44.321 ns | 34.760 ns | `pauli_sum_surface/projection/p0_iz_unit/{side}` |
| improvement | 0.784× | 3 | 44.035 ns | 34.780 ns | `pauli_sum_surface/projection/p1_iz_unit/{side}` |
| improvement | 0.901× | 3 | 3.709 µs | 3.338 µs | `pauli_sum_surface/rotation_one/rot_xy_r/{side}` |
| improvement | 0.005× | 3 | 467.070 ns | 2.258 ns | `pauli_sum_surface/rotation_one/rotate_1_i/{side}` |
| improvement | 0.630× | 3 | 1.232 µs | 796.030 ns | `pauli_sum_surface/rotation_one/rotate_1_x/{side}` |
| improvement | 0.744× | 3 | 1.233 µs | 917.070 ns | `pauli_sum_surface/rotation_one/rotate_1_y/{side}` |
| improvement | 0.707× | 3 | 1.233 µs | 871.950 ns | `pauli_sum_surface/rotation_one/rotate_1_z/{side}` |
| improvement | 0.808× | 3 | 1.048 µs | 873.120 ns | `pauli_sum_surface/rotation_one/rx/{side}` |
| improvement | 0.842× | 3 | 1.104 µs | 936.310 ns | `pauli_sum_surface/rotation_one/ry/{side}` |
| improvement | 0.837× | 3 | 1.054 µs | 883.020 ns | `pauli_sum_surface/rotation_one/rz/{side}` |
| improvement | 0.690× | 3 | 6.298 µs | 4.351 µs | `pauli_sum_surface/rotation_one_batch/rx/{side}` |
| improvement | 0.701× | 3 | 6.699 µs | 4.696 µs | `pauli_sum_surface/rotation_one_batch/ry/{side}` |
| improvement | 0.712× | 3 | 6.771 µs | 4.805 µs | `pauli_sum_surface/rotation_one_batch/rz/{side}` |
| improvement | 0.908× | 3 | 1.611 µs | 1.448 µs | `pauli_sum_surface/rotation_two/rotate_2_generic_xz/{side}` |
| improvement | 0.798× | 3 | 1.179 µs | 956.330 ns | `pauli_sum_surface/rotation_two/rxx/{side}` |
| improvement | 0.939× | 3 | 1.569 µs | 1.474 µs | `pauli_sum_surface/rotation_two/rxy/{side}` |
| improvement | 0.905× | 3 | 1.616 µs | 1.463 µs | `pauli_sum_surface/rotation_two/rxz/{side}` |
| improvement | 0.926× | 3 | 1.600 µs | 1.481 µs | `pauli_sum_surface/rotation_two/ryx/{side}` |
| improvement | 0.876× | 3 | 1.300 µs | 1.166 µs | `pauli_sum_surface/rotation_two/ryy/{side}` |
| improvement | 0.930× | 3 | 1.567 µs | 1.457 µs | `pauli_sum_surface/rotation_two/ryz/{side}` |
| improvement | 0.923× | 3 | 1.582 µs | 1.461 µs | `pauli_sum_surface/rotation_two/rzx/{side}` |
| improvement | 0.948× | 3 | 1.651 µs | 1.566 µs | `pauli_sum_surface/rotation_two/rzy/{side}` |
| improvement | 0.782× | 3 | 1.157 µs | 917.700 ns | `pauli_sum_surface/rotation_two/rzz/{side}` |
| improvement | 0.650× | 3 | 7.371 µs | 4.814 µs | `pauli_sum_surface/rotation_two_batch/rxx/{side}` |
| parity | 1.016× | 3 | 9.454 µs | 9.650 µs | `pauli_sum_surface/rotation_two_batch/rxy/{side}` |
| parity | 0.990× | 3 | 9.578 µs | 9.459 µs | `pauli_sum_surface/rotation_two_batch/rxz/{side}` |
| parity | 1.019× | 3 | 9.548 µs | 9.669 µs | `pauli_sum_surface/rotation_two_batch/ryx/{side}` |
| improvement | 0.686× | 3 | 7.414 µs | 5.112 µs | `pauli_sum_surface/rotation_two_batch/ryy/{side}` |
| parity | 1.015× | 3 | 9.586 µs | 9.708 µs | `pauli_sum_surface/rotation_two_batch/ryz/{side}` |
| parity | 0.986× | 3 | 9.599 µs | 9.462 µs | `pauli_sum_surface/rotation_two_batch/rzx/{side}` |
| parity | 1.021× | 3 | 9.766 µs | 9.943 µs | `pauli_sum_surface/rotation_two_batch/rzy/{side}` |
| improvement | 0.671× | 3 | 7.135 µs | 4.812 µs | `pauli_sum_surface/rotation_two_batch/rzz/{side}` |
| parity | 0.988× | 7 | 109.430 ns | 106.550 ns | `pauli_sum_surface/truncate/coefficient_active/{side}` |
| improvement | 0.926× | 7 | 113.190 ns | 104.620 ns | `pauli_sum_surface/truncate/coefficient_disabled/{side}` |
| improvement | 0.892× | 7 | 176.310 ns | 158.120 ns | `pauli_sum_surface/truncate/combined_active/{side}` |
| improvement | 0.943× | 7 | 106.840 ns | 100.400 ns | `pauli_sum_surface/truncate/combined_disabled/{side}` |
| confirmed regression | 1.038× | 7 | 257.580 ns | 260.690 ns | `pauli_sum_surface/truncate/max_loss_weight_active/{side}` |
| confirmed regression | 1.668× | 7 | 1.574 ns | 2.632 ns | `pauli_sum_surface/truncate/max_loss_weight_disabled/{side}` |
| improvement | 0.901× | 7 | 123.140 ns | 111.690 ns | `pauli_sum_surface/truncate/max_weight_active/{side}` |
| confirmed regression | 1.668× | 7 | 1.574 ns | 2.609 ns | `pauli_sum_surface/truncate/max_weight_disabled/{side}` |
| parity | 0.982× | 7 | 110.370 ns | 108.040 ns | `pauli_sum_surface/truncate/preserve_empty/{side}` |
| improvement | 0.744× | 7 | 278.080 ns | 208.830 ns | `pauli_sum_surface/truncate/preserve_nonempty/{side}` |
| improvement | 0.735× | 1 | 1.607 ns | 1.182 ns | `pauli_word/cnot/{side}/cnot` |
| parity | 1.009× | 1 | 0.701 ns | 0.707 ns | `pauli_word/weight/{side}/weight` |
| improvement | 0.935× | 1 | 1.631 ns | 1.524 ns | `phased_pauli_word/cnot/{side}/cnot` |
| parity | 0.998× | 1 | 2.589 ns | 2.585 ns | `phased_pauli_word/product/{side}/phased_mul_assign` |
| improvement | 0.404× | 2 | 999.240 µs | 403.640 µs | `sym/expectation_eval/{side}/eval_grid_1000` |
| improvement | 0.327× | 2 | 312.425 µs | 102.064 µs | `sym/expectation_propagate/{side}/propagate_trace` |
| improvement | 0.700× | 2 | 1.910 µs | 1.336 µs | `sym/micro_mul_term/{side}/mul_term` |
| improvement | 0.661× | 2 | 74.941 ns | 49.524 ns | `sym/micro_prod_mul/{side}/prod_mul` |
| improvement | 0.649× | 2 | 978.385 ns | 634.930 ns | `sym/micro_term_add/{side}/sum_plus_sum` |
| improvement | 0.758× | 2 | 27.349 µs | 20.746 µs | `sym/micro_term_mul/{side}/sum_x_sum` |
| improvement | 0.731× | 2 | 324.030 ms | 236.905 ms | `sym/random_circuit/{side}/full_replay` |
| improvement | 0.703× | 2 | 293.920 µs | 206.650 µs | `sym/random_circuit_clifford/{side}/clifford_prefix` |
| confirmed regression | 1.477× | 6 | 5.961 ns | 8.598 ns | `sym/surface/construct/{side}/fold_cos_constant` |
| confirmed regression | 1.388× | 6 | 5.867 ns | 8.059 ns | `sym/surface/construct/{side}/fold_sin_constant` |
| improvement | 0.679× | 6 | 16.201 ns | 10.965 ns | `sym/surface/construct/{side}/prod_cos` |
| improvement | 0.127× | 6 | 3.974 ns | 0.505 ns | `sym/surface/construct/{side}/prod_new` |
| improvement | 0.670× | 6 | 16.387 ns | 10.918 ns | `sym/surface/construct/{side}/prod_sin` |
| improvement | 0.613× | 6 | 22.377 ns | 14.123 ns | `sym/surface/construct/{side}/promote_cos` |
| improvement | 0.604× | 6 | 22.133 ns | 13.358 ns | `sym/surface/construct/{side}/promote_sin` |
| confirmed regression | 1.289× | 6 | 1.502 ns | 1.936 ns | `sym/surface/construct/{side}/sum_new` |
| parity | 1.004× | 6 | 0.425 ns | 0.428 ns | `sym/surface/construct/{side}/term_constant` |
| parity | 0.995× | 6 | 0.334 ns | 0.333 ns | `sym/surface/construct/{side}/term_variable` |
| improvement | 0.837× | 6 | 95.439 ns | 78.454 ns | `sym/surface/eval/{side}/prod` |
| confirmed regression | 1.304× | 6 | 81.417 ns | 108.250 ns | `sym/surface/eval/{side}/sum` |
| confirmed regression | 1.061× | 6 | 171.675 ns | 181.385 ns | `sym/surface/eval/{side}/term` |
| improvement | 0.444× | 1 | 34.714 ns | 15.427 ns | `sym/surface/observable/product/{side}/clone` |
| parity | 1.005× | 1 | 0.247 ns | 0.248 ns | `sym/surface/observable/product/{side}/cos_pow` |
| parity | 0.972× | 1 | 218.520 ns | 212.450 ns | `sym/surface/observable/product/{side}/display` |
| improvement | 0.137× | 1 | 12.648 ns | 1.734 ns | `sym/surface/observable/product/{side}/equality` |
| improvement | 0.374× | 1 | 7.340 ns | 2.743 ns | `sym/surface/observable/product/{side}/hash` |
| parity | 1.004× | 1 | 0.248 ns | 0.249 ns | `sym/surface/observable/product/{side}/pow` |
| parity | 1.001× | 1 | 0.249 ns | 0.249 ns | `sym/surface/observable/product/{side}/sin_pow` |
| improvement | 0.608× | 1 | 98.175 ns | 59.699 ns | `sym/surface/observable/sum/{side}/clone` |
| parity | 0.992× | 1 | 584.330 ns | 579.730 ns | `sym/surface/observable/sum/{side}/display` |
| improvement | 0.401× | 1 | 37.804 ns | 15.170 ns | `sym/surface/observable/sum/{side}/equality` |
| improvement | 0.607× | 1 | 98.162 ns | 59.600 ns | `sym/surface/observable/term/{side}/clone` |
| parity | 0.984× | 1 | 587.000 ns | 577.330 ns | `sym/surface/observable/term/{side}/display` |
| improvement | 0.406× | 1 | 37.908 ns | 15.385 ns | `sym/surface/observable/term/{side}/equality` |
| provisional regression | 1.945× | 2 | 0.798 ns | 1.551 ns | `sym/surface/operator_add/{side}/sum_add_coefficient` |
| improvement | 0.853× | 2 | 30.362 ns | 25.891 ns | `sym/surface/operator_add/{side}/sum_add_term` |
| provisional regression | 1.048× | 2 | 6.339 ns | 6.643 ns | `sym/surface/operator_add/{side}/term_add_coefficient` |
| improvement | 0.464× | 2 | 31.441 ns | 14.517 ns | `sym/surface/operator_add/{side}/term_negate` |
| parity | 1.013× | 2 | 7.190 ns | 7.286 ns | `sym/surface/operator_add/{side}/term_subtract_coefficient` |
| improvement | 0.534× | 2 | 48.389 ns | 25.775 ns | `sym/surface/operator_add/{side}/term_subtract_term` |
| improvement | 0.672× | 2 | 182.015 ns | 122.175 ns | `sym/surface/pauli_sum/{side}/mul_coefficient` |
| improvement | 0.857× | 2 | 1.825 ns | 1.564 ns | `sym/surface/product/{side}/add_phase` |
| improvement | 0.756× | 2 | 8.720 ns | 6.594 ns | `sym/surface/product/{side}/mul_cos` |
| improvement | 0.566× | 2 | 11.736 ns | 6.639 ns | `sym/surface/product/{side}/mul_sin` |
| improvement | 0.436× | 2 | 31.893 ns | 13.889 ns | `sym/surface/product/{side}/multiply` |
| parity | 0.972× | 2 | 4.373 µs | 4.250 µs | `sym/surface/propagation/clifford/{side}/alias_cx` |
| parity | 0.980× | 2 | 4.303 µs | 4.217 µs | `sym/surface/propagation/clifford/{side}/alias_zcx` |
| provisional regression | 1.065× | 2 | 4.183 µs | 4.461 µs | `sym/surface/propagation/clifford/{side}/alias_zcy` |
| parity | 1.001× | 2 | 4.305 µs | 4.312 µs | `sym/surface/propagation/clifford/{side}/alias_zcz` |
| parity | 0.972× | 2 | 9.178 µs | 8.895 µs | `sym/surface/propagation/clifford/{side}/batch_cnot` |
| improvement | 0.777× | 2 | 11.568 µs | 8.978 µs | `sym/surface/propagation/clifford/{side}/batch_cy` |
| parity | 1.016× | 2 | 9.221 µs | 9.371 µs | `sym/surface/propagation/clifford/{side}/batch_cz` |
| improvement | 0.639× | 2 | 17.009 µs | 10.878 µs | `sym/surface/propagation/clifford/{side}/batch_h` |
| improvement | 0.738× | 2 | 16.718 µs | 12.352 µs | `sym/surface/propagation/clifford/{side}/batch_s` |
| improvement | 0.792× | 2 | 16.771 µs | 13.317 µs | `sym/surface/propagation/clifford/{side}/batch_s_dag` |
| improvement | 0.853× | 2 | 13.335 µs | 11.381 µs | `sym/surface/propagation/clifford/{side}/batch_sqrt_x` |
| improvement | 0.875× | 2 | 13.271 µs | 11.630 µs | `sym/surface/propagation/clifford/{side}/batch_sqrt_x_dag` |
| improvement | 0.670× | 2 | 16.265 µs | 10.908 µs | `sym/surface/propagation/clifford/{side}/batch_sqrt_y` |
| improvement | 0.710× | 2 | 16.616 µs | 11.826 µs | `sym/surface/propagation/clifford/{side}/batch_sqrt_y_dag` |
| improvement | 0.796× | 2 | 7.033 µs | 5.596 µs | `sym/surface/propagation/clifford/{side}/batch_x` |
| improvement | 0.803× | 2 | 7.016 µs | 5.633 µs | `sym/surface/propagation/clifford/{side}/batch_y` |
| improvement | 0.817× | 2 | 7.017 µs | 5.734 µs | `sym/surface/propagation/clifford/{side}/batch_z` |
| improvement | 0.969× | 2 | 4.351 µs | 4.217 µs | `sym/surface/propagation/clifford/{side}/cnot` |
| provisional regression | 1.116× | 2 | 4.307 µs | 4.801 µs | `sym/surface/propagation/clifford/{side}/cy` |
| parity | 0.989× | 2 | 4.412 µs | 4.367 µs | `sym/surface/propagation/clifford/{side}/cz` |
| improvement | 0.944× | 2 | 4.611 µs | 4.355 µs | `sym/surface/propagation/clifford/{side}/h` |
| provisional regression | 1.035× | 2 | 4.077 µs | 4.221 µs | `sym/surface/propagation/clifford/{side}/s` |
| provisional regression | 1.208× | 2 | 4.449 µs | 5.359 µs | `sym/surface/propagation/clifford/{side}/s_dag` |
| provisional regression | 1.151× | 2 | 4.164 µs | 4.809 µs | `sym/surface/propagation/clifford/{side}/sqrt_x` |
| parity | 1.023× | 2 | 4.177 µs | 4.269 µs | `sym/surface/propagation/clifford/{side}/sqrt_x_dag` |
| parity | 0.984× | 2 | 4.394 µs | 4.327 µs | `sym/surface/propagation/clifford/{side}/sqrt_y` |
| improvement | 0.896× | 2 | 4.757 µs | 4.215 µs | `sym/surface/propagation/clifford/{side}/sqrt_y_dag` |
| improvement | 0.792× | 2 | 2.678 µs | 2.119 µs | `sym/surface/propagation/clifford/{side}/x` |
| improvement | 0.826× | 2 | 2.570 µs | 2.122 µs | `sym/surface/propagation/clifford/{side}/y` |
| parity | 1.029× | 2 | 2.595 µs | 2.676 µs | `sym/surface/propagation/clifford/{side}/z` |
| improvement | 0.745× | 2 | 7.522 µs | 5.568 µs | `sym/surface/propagation/noise/{side}/batch_depolarize1` |
| improvement | 0.752× | 2 | 5.244 µs | 3.929 µs | `sym/surface/propagation/noise/{side}/batch_depolarize2` |
| improvement | 0.779× | 2 | 7.497 µs | 5.840 µs | `sym/surface/propagation/noise/{side}/batch_pauli_error` |
| improvement | 0.779× | 2 | 6.130 µs | 4.774 µs | `sym/surface/propagation/noise/{side}/batch_two_qubit_pauli_error` |
| improvement | 0.810× | 2 | 7.490 µs | 6.063 µs | `sym/surface/propagation/noise/{side}/batch_x_error` |
| improvement | 0.761× | 2 | 8.020 µs | 6.109 µs | `sym/surface/propagation/noise/{side}/batch_y_error` |
| improvement | 0.824× | 2 | 7.510 µs | 6.204 µs | `sym/surface/propagation/noise/{side}/batch_z_error` |
| improvement | 0.926× | 2 | 2.724 µs | 2.534 µs | `sym/surface/propagation/noise/{side}/depolarize1` |
| improvement | 0.714× | 2 | 3.104 µs | 2.215 µs | `sym/surface/propagation/noise/{side}/depolarize2` |
| improvement | 0.821× | 2 | 2.854 µs | 2.343 µs | `sym/surface/propagation/noise/{side}/pauli_error` |
| improvement | 0.831× | 2 | 3.632 µs | 2.960 µs | `sym/surface/propagation/noise/{side}/two_qubit_pauli_error` |
| improvement | 0.761× | 2 | 2.920 µs | 2.222 µs | `sym/surface/propagation/noise/{side}/x_error` |
| improvement | 0.746× | 2 | 2.965 µs | 2.204 µs | `sym/surface/propagation/noise/{side}/y_error` |
| improvement | 0.927× | 2 | 2.829 µs | 2.624 µs | `sym/surface/propagation/noise/{side}/z_error` |
| improvement | 0.839× | 2 | 11.983 µs | 10.054 µs | `sym/surface/propagation/rotation_one/{side}/batch_rx` |
| improvement | 0.855× | 2 | 21.368 µs | 18.289 µs | `sym/surface/propagation/rotation_one/{side}/batch_ry` |
| improvement | 0.952× | 2 | 11.753 µs | 11.223 µs | `sym/surface/propagation/rotation_one/{side}/batch_rz` |
| improvement | 0.924× | 2 | 4.476 µs | 4.089 µs | `sym/surface/propagation/rotation_one/{side}/generic_rotate_1` |
| improvement | 0.924× | 2 | 12.028 µs | 11.120 µs | `sym/surface/propagation/rotation_one/{side}/rot_xy_r` |
| parity | 0.998× | 2 | 4.398 µs | 4.390 µs | `sym/surface/propagation/rotation_one/{side}/rx` |
| parity | 1.027× | 2 | 4.993 µs | 5.133 µs | `sym/surface/propagation/rotation_one/{side}/ry` |
| provisional regression | 1.045× | 2 | 3.660 µs | 3.831 µs | `sym/surface/propagation/rotation_one/{side}/rz` |
| improvement | 0.888× | 2 | 8.266 µs | 7.288 µs | `sym/surface/propagation/rotation_two/{side}/batch_rxx` |
| improvement | 0.910× | 2 | 7.476 µs | 6.793 µs | `sym/surface/propagation/rotation_two/{side}/batch_rxy` |
| improvement | 0.941× | 2 | 7.607 µs | 7.170 µs | `sym/surface/propagation/rotation_two/{side}/batch_rxz` |
| improvement | 0.900× | 2 | 7.206 µs | 6.489 µs | `sym/surface/propagation/rotation_two/{side}/batch_ryx` |
| improvement | 0.892× | 2 | 11.735 µs | 10.462 µs | `sym/surface/propagation/rotation_two/{side}/batch_ryy` |
| improvement | 0.937× | 2 | 7.140 µs | 6.683 µs | `sym/surface/propagation/rotation_two/{side}/batch_ryz` |
| improvement | 0.850× | 2 | 8.961 µs | 7.622 µs | `sym/surface/propagation/rotation_two/{side}/batch_rzx` |
| improvement | 0.872× | 2 | 7.529 µs | 6.566 µs | `sym/surface/propagation/rotation_two/{side}/batch_rzy` |
| improvement | 0.821× | 2 | 9.176 µs | 7.530 µs | `sym/surface/propagation/rotation_two/{side}/batch_rzz` |
| improvement | 0.915× | 2 | 4.525 µs | 4.169 µs | `sym/surface/propagation/rotation_two/{side}/generic_rotate_2` |
| provisional regression | 1.110× | 2 | 3.782 µs | 4.202 µs | `sym/surface/propagation/rotation_two/{side}/rxx` |
| improvement | 0.877× | 2 | 4.578 µs | 4.003 µs | `sym/surface/propagation/rotation_two/{side}/rxy` |
| provisional regression | 1.032× | 2 | 3.674 µs | 3.798 µs | `sym/surface/propagation/rotation_two/{side}/rxz` |
| improvement | 0.816× | 2 | 2.808 µs | 2.294 µs | `sym/surface/propagation/rotation_two/{side}/ryx` |
| provisional regression | 1.096× | 2 | 4.880 µs | 5.319 µs | `sym/surface/propagation/rotation_two/{side}/ryy` |
| improvement | 0.779× | 2 | 2.707 µs | 2.108 µs | `sym/surface/propagation/rotation_two/{side}/ryz` |
| parity | 0.982× | 2 | 4.273 µs | 4.190 µs | `sym/surface/propagation/rotation_two/{side}/rzx` |
| improvement | 0.914× | 2 | 4.049 µs | 3.692 µs | `sym/surface/propagation/rotation_two/{side}/rzy` |
| improvement | 0.938× | 2 | 4.421 µs | 4.125 µs | `sym/surface/propagation/rotation_two/{side}/rzz` |
| provisional regression | 1.410× | 2 | 2.803 µs | 3.947 µs | `sym/surface/propagation/{side}/cnot` |
| provisional regression | 1.604× | 2 | 2.598 µs | 4.151 µs | `sym/surface/propagation/{side}/cz` |
| provisional regression | 1.713× | 2 | 2.611 µs | 4.475 µs | `sym/surface/propagation/{side}/h` |
| provisional regression | 1.448× | 2 | 2.952 µs | 4.248 µs | `sym/surface/propagation/{side}/rx` |
| provisional regression | 1.334× | 2 | 3.120 µs | 4.135 µs | `sym/surface/propagation/{side}/ry` |
| provisional regression | 1.291× | 2 | 2.002 µs | 2.552 µs | `sym/surface/propagation/{side}/rz` |
| provisional regression | 1.630× | 2 | 2.661 µs | 4.353 µs | `sym/surface/propagation/{side}/s` |
| provisional regression | 1.248× | 2 | 130.135 µs | 162.545 µs | `sym/surface/readout/{side}/trace` |
| provisional regression | 1.999× | 2 | 0.837 ns | 1.674 ns | `sym/surface/sum/{side}/add_const` |
| improvement | 0.784× | 2 | 62.398 ns | 48.946 ns | `sym/surface/sum/{side}/add_term` |
| improvement | 0.757× | 2 | 86.392 ns | 65.438 ns | `sym/surface/sum/{side}/mul_scalar` |
| improvement | 0.658× | 2 | 1.390 µs | 909.845 ns | `sym/surface/sum/{side}/mul_term` |
| improvement | 0.667× | 2 | 304.455 ns | 203.045 ns | `sym/surface/term/{side}/add` |
| improvement | 0.674× | 2 | 94.545 ns | 63.745 ns | `sym/surface/term/{side}/mul_scalar` |
| improvement | 0.781× | 2 | 2.366 µs | 1.848 µs | `sym/surface/term/{side}/multiply` |
| improvement | 0.859× | 2 | 2.349 ns | 2.019 ns | `sym/surface/term_setters/{side}/set_max_sin` |
| improvement | 0.956× | 2 | 2.347 ns | 2.246 ns | `sym/surface/term_setters/{side}/set_min_eps` |
| improvement | 0.306× | 2 | 4.076 ms | 1.247 ms | `sym/tfim_trotter_k3/{side}/trotter` |
| improvement | 0.211× | 2 | 12.324 ms | 2.599 ms | `sym/tfim_trotter_k4/{side}/trotter` |
| improvement | 0.501× | 2 | 7.256 µs | 3.633 µs | `sym/trace_parametric/{side}/build_propagate_trace_eval` |
| provisional regression | 1.280× | 2 | 767.460 µs | 982.520 µs | `sym/trace_readout_k2/{side}/trace` |
| provisional regression | 1.294× | 2 | 763.795 µs | 988.305 µs | `sym/trace_readout_k3/{side}/trace` |
| provisional regression | 1.271× | 2 | 783.155 µs | 995.660 µs | `sym/trace_readout_k4/{side}/trace` |
| confirmed regression | 1.194× | 6 | 962.690 µs | 1.155 ms | `sym/trace_readout_k5/{side}/trace` |
| improvement | 0.532× | 2 | 85.517 µs | 45.501 µs | `sym/truncation_sweep/{side}/k1` |
| improvement | 0.406× | 2 | 151.390 µs | 61.469 µs | `sym/truncation_sweep/{side}/k2` |
| improvement | 0.324× | 2 | 312.225 µs | 101.272 µs | `sym/truncation_sweep/{side}/k3` |
| improvement | 0.268× | 2 | 632.145 µs | 169.275 µs | `sym/truncation_sweep/{side}/k4` |
| improvement | 0.253× | 2 | 1.036 ms | 261.820 µs | `sym/truncation_sweep/{side}/k5` |
| improvement | 0.326× | 2 | 311.295 µs | 101.357 µs | `sym/truncation_sweep_eps/{side}/1e-12` |
| improvement | 0.327× | 2 | 311.620 µs | 101.757 µs | `sym/truncation_sweep_eps/{side}/1e-6` |
| improvement | 0.326× | 2 | 313.630 µs | 102.224 µs | `sym/truncation_sweep_eps/{side}/eps` |
| improvement | 0.864× | 1 | 36.562 µs | 31.591 µs | `tableau-attrib/measure-sweep/decomp_only/{side}/128` |
| improvement | 0.825× | 1 | 2.568 µs | 2.117 µs | `tableau-attrib/measure-sweep/decomp_only/{side}/32` |
| improvement | 0.922× | 1 | 48.025 µs | 44.297 µs | `tableau-attrib/measure-sweep/decomp_sweep/{side}/128` |
| improvement | 0.926× | 1 | 4.070 µs | 3.771 µs | `tableau-attrib/measure-sweep/decomp_sweep/{side}/32` |
| parity | 0.975× | 1 | 25.015 µs | 24.402 µs | `tableau-attrib/measure-sweep/frame_sweep/{side}/128` |
| parity | 0.972× | 1 | 1.707 µs | 1.659 µs | `tableau-attrib/measure-sweep/frame_sweep/{side}/32` |
| improvement | 0.922× | 1 | 36.767 µs | 33.896 µs | `tableau-attrib/measure-sweep/many_sweep/{side}/128` |
| improvement | 0.913× | 1 | 2.680 µs | 2.447 µs | `tableau-attrib/measure-sweep/many_sweep/{side}/32` |
| parity | 0.995× | 1 | 115.360 µs | 114.740 µs | `tableau-integration/branch-coalesce/doubling/{side}/16384` |
| parity | 1.003× | 1 | 14.298 µs | 14.348 µs | `tableau-integration/branch-coalesce/doubling/{side}/2048` |
| parity | 0.990× | 1 | 1.905 µs | 1.886 µs | `tableau-integration/branch-coalesce/doubling/{side}/256` |
| improvement | 0.926× | 1 | 398.420 ns | 368.790 ns | `tableau-integration/branch-coalesce/doubling/{side}/32` |
| improvement | 0.856× | 1 | 217.610 ns | 186.340 ns | `tableau-integration/branch-coalesce/doubling/{side}/4` |
| parity | 0.998× | 1 | 511.970 µs | 511.170 µs | `tableau-integration/branch-coalesce/doubling/{side}/65536` |
| parity | 1.000× | 1 | 176.940 µs | 176.880 µs | `tableau-integration/branch-coalesce/merge/{side}/16384` |
| parity | 1.013× | 1 | 19.543 µs | 19.799 µs | `tableau-integration/branch-coalesce/merge/{side}/2048` |
| parity | 0.990× | 1 | 2.273 µs | 2.252 µs | `tableau-integration/branch-coalesce/merge/{side}/256` |
| improvement | 0.950× | 1 | 412.110 ns | 391.600 ns | `tableau-integration/branch-coalesce/merge/{side}/32` |
| improvement | 0.846× | 1 | 215.760 ns | 182.460 ns | `tableau-integration/branch-coalesce/merge/{side}/4` |
| parity | 1.007× | 1 | 1.296 ms | 1.305 ms | `tableau-integration/branch-coalesce/merge/{side}/65536` |
| improvement | 0.969× | 1 | 1.292 ms | 1.252 ms | `tableau-integration/fused-tgate/{side}/12t-85q` |
| parity | 0.992× | 1 | 27.326 ms | 27.112 ms | `tableau-integration/fused-tgate/{side}/16t-85q` |
| improvement | 0.853× | 1 | 66.249 µs | 56.479 µs | `tableau-integration/fused-tgate/{side}/8t-85q` |
| improvement | 0.797× | 1 | 25.884 µs | 20.641 µs | `tableau-integration/measure-all-msd/measure_all/{side}` |
| improvement | 0.792× | 1 | 26.280 µs | 20.821 µs | `tableau-integration/measure-all-msd/measure_loop/{side}` |
| improvement | 0.796× | 1 | 26.065 µs | 20.758 µs | `tableau-integration/measure-all-msd/measure_many/{side}` |
| improvement | 0.875× | 1 | 64.335 µs | 56.277 µs | `tableau-integration/msd-85q/fused/{side}` |
| improvement | 0.938× | 1 | 101.310 µs | 95.054 µs | `tableau-integration/msd-85q/naive/{side}` |
| improvement | 0.913× | 1 | 509.480 µs | 465.010 µs | `tableau-integration/noisy-shots/{side}` |
| improvement | 0.288× | 1 | 3.031 ms | 871.460 µs | `tableau-integration/rot2-brickwork/{side}/n10_l4_m1024` |
| improvement | 0.361× | 1 | 11.216 ms | 4.050 ms | `tableau-integration/rot2-brickwork/{side}/n12_l3_m4096` |
| improvement | 0.284× | 1 | 565.000 µs | 160.550 µs | `tableau-integration/rot2-brickwork/{side}/n8_l4_m256` |
| parity | 0.991× | 1 | 30.157 µs | 29.882 µs | `tableau-integration/scaling/gates/{side}/128` |
| parity | 0.994× | 1 | 2.456 µs | 2.440 µs | `tableau-integration/scaling/gates/{side}/32` |
| parity | 0.992× | 1 | 8.269 µs | 8.200 µs | `tableau-integration/scaling/gates/{side}/64` |
| parity | 0.994× | 1 | 17.434 µs | 17.321 µs | `tableau-integration/scaling/gates/{side}/96` |
| improvement | 0.895× | 1 | 36.817 µs | 32.948 µs | `tableau-integration/scaling/measure-sweep/{side}/128` |
| improvement | 0.914× | 1 | 2.678 µs | 2.447 µs | `tableau-integration/scaling/measure-sweep/{side}/32` |
| improvement | 0.872× | 1 | 10.044 µs | 8.757 µs | `tableau-integration/scaling/measure-sweep/{side}/64` |
| improvement | 0.888× | 1 | 21.310 µs | 18.927 µs | `tableau-integration/scaling/measure-sweep/{side}/96` |
| improvement | 0.916× | 1 | 322.090 ns | 294.930 ns | `tableau-micro/case_a_m1/{side}` |
| parity | 1.006× | 1 | 723.340 ns | 727.570 ns | `tableau-micro/case_a_m32/{side}` |
| improvement | 0.847× | 1 | 555.420 ns | 470.280 ns | `tableau-micro/construct/{side}` |
| improvement | 0.863× | 1 | 244.770 ns | 211.330 ns | `tableau-micro/cz_block17/{side}` |
| parity | 0.987× | 1 | 2.632 µs | 2.598 µs | `tableau-micro/cz_loop17/{side}` |
| improvement | 0.660× | 1 | 198.590 ns | 131.080 ns | `tableau-micro/frame_project/{side}` |
| improvement | 0.735× | 1 | 189.660 ns | 139.320 ns | `tableau-micro/measure/{side}` |
| confirmed regression | 1.120× | 4 | 681.085 ns | 763.820 ns | `tableau-micro/msd_measure_single/{side}` |
| improvement | 0.785× | 1 | 26.495 µs | 20.811 µs | `tableau-micro/msd_sweep_all/{side}` |
| improvement | 0.792× | 1 | 26.916 µs | 21.313 µs | `tableau-micro/msd_sweep_loop/{side}` |
| improvement | 0.920× | 1 | 370.990 ns | 341.440 ns | `tableau-micro/msd_z_expectation/{side}` |
| provisional regression | 3.506× | 1 | 74.001 ns | 259.460 ns | `tableau-micro/scratch_new_x85/{side}` |
| improvement | 0.711× | 1 | 230.870 ns | 164.080 ns | `tableau-micro/sqrt_y/{side}` |
| parity | 0.982× | 1 | 1.750 µs | 1.718 µs | `tableau-micro/sqrt_y_loop17/{side}` |
| improvement | 0.792× | 1 | 239.230 ns | 189.360 ns | `tableau-micro/sqrt_y_many17/{side}` |
| improvement | 0.779× | 1 | 243.850 ns | 190.070 ns | `tableau-micro/t-gate/{side}` |
| improvement | 0.720× | 3 | 246.880 ns | 177.640 ns | `tableau-surface/clifford/bare/cnot/{side}` |
| improvement | 0.958× | 3 | 5.338 µs | 5.112 µs | `tableau-surface/clifford/bare/cnot_many/{side}` |
| improvement | 0.723× | 3 | 246.890 ns | 179.130 ns | `tableau-surface/clifford/bare/cx/{side}` |
| improvement | 0.755× | 3 | 254.010 ns | 191.470 ns | `tableau-surface/clifford/bare/cy/{side}` |
| parity | 1.014× | 3 | 6.913 µs | 6.956 µs | `tableau-surface/clifford/bare/cy_many/{side}` |
| improvement | 0.730× | 3 | 243.630 ns | 177.900 ns | `tableau-surface/clifford/bare/cz/{side}` |
| improvement | 0.740× | 3 | 248.610 ns | 182.120 ns | `tableau-surface/clifford/bare/cz_block_pairs/{side}` |
| improvement | 0.725× | 3 | 249.690 ns | 181.110 ns | `tableau-surface/clifford/bare/cz_block_pairs_cross_word/{side}` |
| confirmed regression | 1.061× | 3 | 5.267 µs | 5.477 µs | `tableau-surface/clifford/bare/cz_many/{side}` |
| improvement | 0.693× | 3 | 233.510 ns | 163.090 ns | `tableau-surface/clifford/bare/h/{side}` |
| improvement | 0.884× | 3 | 366.140 ns | 323.520 ns | `tableau-surface/clifford/bare/h_many/{side}` |
| improvement | 0.688× | 3 | 221.410 ns | 152.100 ns | `tableau-surface/clifford/bare/s/{side}` |
| improvement | 0.688× | 3 | 220.860 ns | 151.900 ns | `tableau-surface/clifford/bare/s_dag/{side}` |
| confirmed regression | 1.105× | 7 | 429.710 ns | 475.110 ns | `tableau-surface/clifford/bare/s_dag_many/{side}` |
| confirmed regression | 1.098× | 7 | 423.040 ns | 463.920 ns | `tableau-surface/clifford/bare/s_many/{side}` |
| improvement | 0.698× | 3 | 223.690 ns | 156.860 ns | `tableau-surface/clifford/bare/sqrt_x/{side}` |
| improvement | 0.712× | 3 | 221.010 ns | 157.460 ns | `tableau-surface/clifford/bare/sqrt_x_dag/{side}` |
| improvement | 0.831× | 3 | 279.520 ns | 232.230 ns | `tableau-surface/clifford/bare/sqrt_x_dag_many/{side}` |
| improvement | 0.850× | 3 | 344.390 ns | 291.970 ns | `tableau-surface/clifford/bare/sqrt_x_many/{side}` |
| improvement | 0.697× | 3 | 236.240 ns | 164.230 ns | `tableau-surface/clifford/bare/sqrt_y/{side}` |
| improvement | 0.699× | 3 | 236.490 ns | 165.020 ns | `tableau-surface/clifford/bare/sqrt_y_dag/{side}` |
| improvement | 0.875× | 3 | 370.790 ns | 324.490 ns | `tableau-surface/clifford/bare/sqrt_y_dag_many/{side}` |
| improvement | 0.878× | 3 | 371.330 ns | 325.830 ns | `tableau-surface/clifford/bare/sqrt_y_many/{side}` |
| improvement | 0.690× | 3 | 196.590 ns | 135.650 ns | `tableau-surface/clifford/bare/x/{side}` |
| confirmed regression | 1.059× | 7 | 363.820 ns | 389.460 ns | `tableau-surface/clifford/bare/x_many/{side}` |
| improvement | 0.688× | 3 | 207.360 ns | 142.490 ns | `tableau-surface/clifford/bare/y/{side}` |
| confirmed regression | 1.119× | 7 | 391.960 ns | 439.280 ns | `tableau-surface/clifford/bare/y_many/{side}` |
| improvement | 0.683× | 3 | 199.330 ns | 136.890 ns | `tableau-surface/clifford/bare/z/{side}` |
| confirmed regression | 1.076× | 7 | 363.660 ns | 393.650 ns | `tableau-surface/clifford/bare/z_many/{side}` |
| improvement | 0.721× | 3 | 246.600 ns | 178.250 ns | `tableau-surface/clifford/bare/zcx/{side}` |
| improvement | 0.756× | 3 | 253.770 ns | 192.150 ns | `tableau-surface/clifford/bare/zcy/{side}` |
| improvement | 0.731× | 3 | 245.530 ns | 179.040 ns | `tableau-surface/clifford/bare/zcz/{side}` |
| improvement | 0.747× | 3 | 246.720 ns | 183.260 ns | `tableau-surface/clifford/generalized/cnot/{side}` |
| parity | 0.982× | 3 | 5.439 µs | 5.275 µs | `tableau-surface/clifford/generalized/cnot_many/{side}` |
| improvement | 0.748× | 3 | 246.600 ns | 184.380 ns | `tableau-surface/clifford/generalized/cx/{side}` |
| improvement | 0.758× | 3 | 254.790 ns | 193.870 ns | `tableau-surface/clifford/generalized/cy/{side}` |
| parity | 1.024× | 3 | 6.779 µs | 6.945 µs | `tableau-surface/clifford/generalized/cy_many/{side}` |
| improvement | 0.747× | 3 | 244.360 ns | 182.960 ns | `tableau-surface/clifford/generalized/cz/{side}` |
| improvement | 0.836× | 3 | 295.450 ns | 247.160 ns | `tableau-surface/clifford/generalized/cz_block/{side}` |
| improvement | 0.823× | 3 | 274.400 ns | 225.850 ns | `tableau-surface/clifford/generalized/cz_block_pairs/{side}` |
| improvement | 0.837× | 3 | 294.890 ns | 247.180 ns | `tableau-surface/clifford/generalized/cz_block_pairs_cross_word/{side}` |
| improvement | 0.955× | 3 | 5.531 µs | 5.288 µs | `tableau-surface/clifford/generalized/cz_many/{side}` |
| improvement | 0.701× | 3 | 235.070 ns | 164.310 ns | `tableau-surface/clifford/generalized/h/{side}` |
| improvement | 0.890× | 3 | 375.940 ns | 333.800 ns | `tableau-surface/clifford/generalized/h_many/{side}` |
| improvement | 0.693× | 3 | 221.740 ns | 153.840 ns | `tableau-surface/clifford/generalized/s/{side}` |
| improvement | 0.670× | 3 | 228.170 ns | 152.900 ns | `tableau-surface/clifford/generalized/s_dag/{side}` |
| confirmed regression | 1.100× | 7 | 444.360 ns | 486.070 ns | `tableau-surface/clifford/generalized/s_dag_many/{side}` |
| confirmed regression | 1.092× | 7 | 435.670 ns | 475.670 ns | `tableau-surface/clifford/generalized/s_many/{side}` |
| improvement | 0.682× | 3 | 231.060 ns | 158.030 ns | `tableau-surface/clifford/generalized/sqrt_x/{side}` |
| improvement | 0.716× | 3 | 220.800 ns | 158.210 ns | `tableau-surface/clifford/generalized/sqrt_x_dag/{side}` |
| improvement | 0.835× | 3 | 286.040 ns | 239.200 ns | `tableau-surface/clifford/generalized/sqrt_x_dag_many/{side}` |
| improvement | 0.867× | 3 | 351.090 ns | 304.320 ns | `tableau-surface/clifford/generalized/sqrt_x_many/{side}` |
| improvement | 0.678× | 3 | 244.580 ns | 166.380 ns | `tableau-surface/clifford/generalized/sqrt_y/{side}` |
| improvement | 0.682× | 3 | 244.060 ns | 166.180 ns | `tableau-surface/clifford/generalized/sqrt_y_dag/{side}` |
| improvement | 0.888× | 3 | 378.090 ns | 335.710 ns | `tableau-surface/clifford/generalized/sqrt_y_dag_many/{side}` |
| improvement | 0.888× | 3 | 377.440 ns | 335.220 ns | `tableau-surface/clifford/generalized/sqrt_y_many/{side}` |
| improvement | 0.699× | 3 | 198.460 ns | 138.540 ns | `tableau-surface/clifford/generalized/x/{side}` |
| confirmed regression | 1.084× | 7 | 371.750 ns | 404.750 ns | `tableau-surface/clifford/generalized/x_many/{side}` |
| improvement | 0.688× | 3 | 207.650 ns | 143.680 ns | `tableau-surface/clifford/generalized/y/{side}` |
| confirmed regression | 1.108× | 7 | 407.280 ns | 451.060 ns | `tableau-surface/clifford/generalized/y_many/{side}` |
| improvement | 0.691× | 3 | 198.190 ns | 136.970 ns | `tableau-surface/clifford/generalized/z/{side}` |
| confirmed regression | 1.072× | 7 | 378.120 ns | 405.590 ns | `tableau-surface/clifford/generalized/z_many/{side}` |
| improvement | 0.742× | 3 | 248.470 ns | 184.250 ns | `tableau-surface/clifford/generalized/zcx/{side}` |
| improvement | 0.760× | 3 | 255.130 ns | 194.470 ns | `tableau-surface/clifford/generalized/zcy/{side}` |
| improvement | 0.748× | 3 | 245.320 ns | 182.980 ns | `tableau-surface/clifford/generalized/zcz/{side}` |
| improvement | 0.803× | 3 | 92.832 ns | 74.358 ns | `tableau-surface/clifford/width-32/bare/cnot-edge/{side}` |
| improvement | 0.741× | 3 | 82.648 ns | 61.087 ns | `tableau-surface/clifford/width-32/bare/h/{side}` |
| improvement | 0.820× | 3 | 115.280 ns | 94.596 ns | `tableau-surface/clifford/width-32/generalized/cnot-edge/{side}` |
| improvement | 0.802× | 3 | 102.960 ns | 82.312 ns | `tableau-surface/clifford/width-32/generalized/h/{side}` |
| improvement | 0.851× | 3 | 146.420 ns | 124.910 ns | `tableau-surface/clifford/width-63/bare/cnot-edge/{side}` |
| improvement | 0.784× | 3 | 129.940 ns | 101.290 ns | `tableau-surface/clifford/width-63/bare/h/{side}` |
| improvement | 0.859× | 3 | 147.020 ns | 126.490 ns | `tableau-surface/clifford/width-63/generalized/cnot-edge/{side}` |
| improvement | 0.785× | 3 | 132.470 ns | 103.990 ns | `tableau-surface/clifford/width-63/generalized/h/{side}` |
| improvement | 0.847× | 3 | 150.180 ns | 127.210 ns | `tableau-surface/clifford/width-64/bare/cnot-edge/{side}` |
| improvement | 0.785× | 3 | 133.300 ns | 104.770 ns | `tableau-surface/clifford/width-64/bare/h/{side}` |
| improvement | 0.855× | 3 | 151.340 ns | 129.690 ns | `tableau-surface/clifford/width-64/generalized/cnot-edge/{side}` |
| improvement | 0.794× | 3 | 135.790 ns | 107.750 ns | `tableau-surface/clifford/width-64/generalized/h/{side}` |
| improvement | 0.796× | 3 | 175.350 ns | 139.300 ns | `tableau-surface/clifford/width-65/bare/cnot-edge/{side}` |
| improvement | 0.735× | 3 | 158.840 ns | 116.680 ns | `tableau-surface/clifford/width-65/bare/h/{side}` |
| improvement | 0.421× | 3 | 345.900 ns | 145.260 ns | `tableau-surface/clifford/width-65/generalized/cnot-edge/{side}` |
| improvement | 0.405× | 3 | 295.470 ns | 119.290 ns | `tableau-surface/clifford/width-65/generalized/h/{side}` |
| improvement | 0.772× | 3 | 26.688 ns | 23.610 ns | `tableau-surface/clifford/width-8/bare/cnot-edge/{side}` |
| improvement | 0.703× | 3 | 28.496 ns | 20.040 ns | `tableau-surface/clifford/width-8/bare/h/{side}` |
| improvement | 0.758× | 3 | 31.216 ns | 24.330 ns | `tableau-surface/clifford/width-8/generalized/cnot-edge/{side}` |
| improvement | 0.735× | 3 | 28.262 ns | 21.245 ns | `tableau-surface/clifford/width-8/generalized/h/{side}` |
| improvement | 0.756× | 3 | 260.620 ns | 197.740 ns | `tableau-surface/clifford/width-96/bare/cnot-edge/{side}` |
| improvement | 0.703× | 3 | 241.250 ns | 169.370 ns | `tableau-surface/clifford/width-96/bare/h/{side}` |
| improvement | 0.770× | 3 | 260.810 ns | 199.800 ns | `tableau-surface/clifford/width-96/generalized/cnot-edge/{side}` |
| improvement | 0.697× | 3 | 243.150 ns | 169.190 ns | `tableau-surface/clifford/width-96/generalized/h/{side}` |
| improvement | 0.699× | 3 | 109.630 ns | 76.978 ns | `tableau-surface/construction/bare/clone/32/{side}` |
| improvement | 0.667× | 3 | 204.650 ns | 136.500 ns | `tableau-surface/construction/bare/clone/63/{side}` |
| improvement | 0.666× | 3 | 207.870 ns | 138.630 ns | `tableau-surface/construction/bare/clone/64/{side}` |
| improvement | 0.678× | 3 | 211.090 ns | 143.350 ns | `tableau-surface/construction/bare/clone/65/{side}` |
| improvement | 0.818× | 3 | 42.409 ns | 34.690 ns | `tableau-surface/construction/bare/clone/8/{side}` |
| improvement | 0.657× | 3 | 301.480 ns | 198.300 ns | `tableau-surface/construction/bare/clone/96/{side}` |
| improvement | 0.857× | 1 | 210.310 ns | 180.260 ns | `tableau-surface/construction/bare/new-entropy/{side}/32` |
| improvement | 0.850× | 1 | 394.410 ns | 335.410 ns | `tableau-surface/construction/bare/new-entropy/{side}/63` |
| improvement | 0.843× | 1 | 401.280 ns | 338.350 ns | `tableau-surface/construction/bare/new-entropy/{side}/64` |
| improvement | 0.847× | 1 | 407.500 ns | 345.280 ns | `tableau-surface/construction/bare/new-entropy/{side}/65` |
| improvement | 0.846× | 1 | 84.210 ns | 71.278 ns | `tableau-surface/construction/bare/new-entropy/{side}/8` |
| improvement | 0.842× | 1 | 586.560 ns | 494.100 ns | `tableau-surface/construction/bare/new-entropy/{side}/96` |
| improvement | 0.844× | 3 | 199.190 ns | 168.560 ns | `tableau-surface/construction/bare/new/{side}/32` |
| improvement | 0.836× | 3 | 374.060 ns | 314.450 ns | `tableau-surface/construction/bare/new/{side}/63` |
| improvement | 0.843× | 3 | 379.470 ns | 319.250 ns | `tableau-surface/construction/bare/new/{side}/64` |
| improvement | 0.848× | 3 | 380.330 ns | 322.970 ns | `tableau-surface/construction/bare/new/{side}/65` |
| improvement | 0.857× | 3 | 78.523 ns | 67.742 ns | `tableau-surface/construction/bare/new/{side}/8` |
| improvement | 0.847× | 3 | 547.740 ns | 463.730 ns | `tableau-surface/construction/bare/new/{side}/96` |
| improvement | 0.868× | 3 | 291.550 ns | 248.600 ns | `tableau-surface/construction/bare/reset_all/32/{side}` |
| improvement | 0.878× | 3 | 426.760 ns | 374.600 ns | `tableau-surface/construction/bare/reset_all/63/{side}` |
| improvement | 0.870× | 3 | 439.480 ns | 378.390 ns | `tableau-surface/construction/bare/reset_all/64/{side}` |
| improvement | 0.864× | 3 | 459.970 ns | 396.760 ns | `tableau-surface/construction/bare/reset_all/65/{side}` |
| improvement | 0.750× | 3 | 78.525 ns | 58.885 ns | `tableau-surface/construction/bare/reset_all/8/{side}` |
| improvement | 0.859× | 3 | 627.990 ns | 544.760 ns | `tableau-surface/construction/bare/reset_all/96/{side}` |
| improvement | 0.858× | 3 | 231.550 ns | 200.110 ns | `tableau-surface/construction/generalized/clone/32/{side}` |
| improvement | 0.793× | 3 | 337.510 ns | 267.860 ns | `tableau-surface/construction/generalized/clone/63/{side}` |
| improvement | 0.800× | 3 | 334.620 ns | 268.170 ns | `tableau-surface/construction/generalized/clone/64/{side}` |
| improvement | 0.814× | 3 | 336.540 ns | 272.730 ns | `tableau-surface/construction/generalized/clone/65/{side}` |
| improvement | 0.966× | 3 | 156.620 ns | 151.270 ns | `tableau-surface/construction/generalized/clone/8/{side}` |
| improvement | 0.745× | 3 | 441.810 ns | 331.100 ns | `tableau-surface/construction/generalized/clone/96/{side}` |
| improvement | 0.854× | 3 | 233.650 ns | 199.360 ns | `tableau-surface/construction/generalized/fork/32/{side}` |
| improvement | 0.793× | 3 | 339.000 ns | 269.840 ns | `tableau-surface/construction/generalized/fork/63/{side}` |
| improvement | 0.786× | 3 | 340.780 ns | 269.790 ns | `tableau-surface/construction/generalized/fork/64/{side}` |
| improvement | 0.832× | 3 | 326.220 ns | 271.100 ns | `tableau-surface/construction/generalized/fork/65/{side}` |
| improvement | 0.957× | 3 | 160.160 ns | 152.950 ns | `tableau-surface/construction/generalized/fork/8/{side}` |
| improvement | 0.753× | 3 | 444.090 ns | 331.730 ns | `tableau-surface/construction/generalized/fork/96/{side}` |
| improvement | 0.870× | 1 | 234.730 ns | 204.140 ns | `tableau-surface/construction/generalized/new-entropy/{side}/32` |
| improvement | 0.852× | 1 | 421.780 ns | 359.240 ns | `tableau-surface/construction/generalized/new-entropy/{side}/63` |
| improvement | 0.849× | 1 | 426.060 ns | 361.870 ns | `tableau-surface/construction/generalized/new-entropy/{side}/64` |
| improvement | 0.847× | 1 | 433.680 ns | 367.160 ns | `tableau-surface/construction/generalized/new-entropy/{side}/65` |
| improvement | 0.883× | 1 | 106.130 ns | 93.704 ns | `tableau-surface/construction/generalized/new-entropy/{side}/8` |
| improvement | 0.851× | 1 | 612.920 ns | 521.860 ns | `tableau-surface/construction/generalized/new-entropy/{side}/96` |
| improvement | 0.873× | 3 | 226.430 ns | 198.820 ns | `tableau-surface/construction/generalized/new/{side}/32` |
| improvement | 0.853× | 3 | 397.110 ns | 338.670 ns | `tableau-surface/construction/generalized/new/{side}/63` |
| improvement | 0.855× | 3 | 403.660 ns | 345.120 ns | `tableau-surface/construction/generalized/new/{side}/64` |
| improvement | 0.857× | 3 | 410.450 ns | 351.520 ns | `tableau-surface/construction/generalized/new/{side}/65` |
| improvement | 0.886× | 3 | 104.570 ns | 91.850 ns | `tableau-surface/construction/generalized/new/{side}/8` |
| improvement | 0.849× | 3 | 576.670 ns | 489.140 ns | `tableau-surface/construction/generalized/new/{side}/96` |
| improvement | 0.869× | 3 | 332.410 ns | 289.690 ns | `tableau-surface/construction/generalized/reset_all/32/{side}` |
| improvement | 0.869× | 3 | 484.000 ns | 421.660 ns | `tableau-surface/construction/generalized/reset_all/63/{side}` |
| improvement | 0.873× | 3 | 485.340 ns | 422.120 ns | `tableau-surface/construction/generalized/reset_all/64/{side}` |
| improvement | 0.867× | 3 | 499.600 ns | 438.400 ns | `tableau-surface/construction/generalized/reset_all/65/{side}` |
| improvement | 0.887× | 3 | 170.070 ns | 150.180 ns | `tableau-surface/construction/generalized/reset_all/8/{side}` |
| improvement | 0.868× | 3 | 681.750 ns | 607.230 ns | `tableau-surface/construction/generalized/reset_all/96/{side}` |
| improvement | 0.616× | 1 | 2.996 µs | 1.845 µs | `tableau-surface/display/bare/{side}` |
| improvement | 0.960× | 1 | 33.367 µs | 32.018 µs | `tableau-surface/display/generalized/{side}` |
| improvement | 0.652× | 3 | 226.870 ns | 147.970 ns | `tableau-surface/measurement/bare/measure-deterministic/{side}` |
| improvement | 0.632× | 3 | 210.990 ns | 133.630 ns | `tableau-surface/measurement/bare/measure-random/{side}` |
| improvement | 0.629× | 3 | 211.530 ns | 133.030 ns | `tableau-surface/measurement/bare/reset/{side}` |
| improvement | 0.866× | 3 | 1.044 µs | 904.120 ns | `tableau-surface/measurement/bare/reset_many/{side}` |
| improvement | 0.721× | 3 | 295.240 ns | 212.210 ns | `tableau-surface/measurement/bare/reset_x/{side}` |
| improvement | 0.948× | 3 | 2.741 µs | 2.600 µs | `tableau-surface/measurement/bare/reset_x_many/{side}` |
| improvement | 0.789× | 3 | 362.290 ns | 285.790 ns | `tableau-surface/measurement/bare/reset_y/{side}` |
| parity | 0.970× | 3 | 3.873 µs | 3.741 µs | `tableau-surface/measurement/bare/reset_y_many/{side}` |
| improvement | 0.633× | 3 | 211.450 ns | 133.790 ns | `tableau-surface/measurement/bare/reset_z/{side}` |
| improvement | 0.859× | 3 | 1.047 µs | 901.170 ns | `tableau-surface/measurement/bare/reset_z_many/{side}` |
| improvement | 0.710× | 3 | 209.150 ns | 148.410 ns | `tableau-surface/measurement/generalized/measure-deterministic/{side}` |
| improvement | 0.910× | 3 | 328.220 ns | 300.850 ns | `tableau-surface/measurement/generalized/measure-random/{side}` |
| improvement | 0.705× | 3 | 13.292 µs | 9.366 µs | `tableau-surface/measurement/generalized/measure_all/{side}` |
| improvement | 0.712× | 3 | 13.261 µs | 9.572 µs | `tableau-surface/measurement/generalized/measure_all_with_scratch/{side}` |
| improvement | 0.760× | 3 | 4.312 µs | 3.275 µs | `tableau-surface/measurement/generalized/measure_many/{side}` |
| improvement | 0.759× | 3 | 4.317 µs | 3.306 µs | `tableau-surface/measurement/generalized/measure_many_with_scratch/{side}` |
| improvement | 0.908× | 3 | 330.190 ns | 299.690 ns | `tableau-surface/measurement/generalized/measure_noisy/{side}` |
| improvement | 0.926× | 3 | 403.790 ns | 374.700 ns | `tableau-surface/measurement/generalized/reset/{side}` |
| improvement | 0.783× | 3 | 4.798 µs | 3.765 µs | `tableau-surface/measurement/generalized/reset_many/{side}` |
| improvement | 0.933× | 3 | 504.200 ns | 468.700 ns | `tableau-surface/measurement/generalized/reset_x/{side}` |
| improvement | 0.841× | 3 | 6.565 µs | 5.519 µs | `tableau-surface/measurement/generalized/reset_x_many/{side}` |
| improvement | 0.952× | 3 | 583.180 ns | 558.440 ns | `tableau-surface/measurement/generalized/reset_y/{side}` |
| improvement | 0.861× | 3 | 7.645 µs | 6.596 µs | `tableau-surface/measurement/generalized/reset_y_many/{side}` |
| improvement | 0.927× | 3 | 407.420 ns | 377.270 ns | `tableau-surface/measurement/generalized/reset_z/{side}` |
| improvement | 0.793× | 3 | 4.791 µs | 3.813 µs | `tableau-surface/measurement/generalized/reset_z_many/{side}` |
| improvement | 0.700× | 3 | 208.070 ns | 145.720 ns | `tableau-surface/noise/bare/depolarize1/{side}` |
| improvement | 0.790× | 3 | 247.480 ns | 196.110 ns | `tableau-surface/noise/bare/depolarize1_many/{side}` |
| improvement | 0.699× | 3 | 209.050 ns | 146.770 ns | `tableau-surface/noise/bare/depolarize2/{side}` |
| improvement | 0.862× | 3 | 294.790 ns | 253.540 ns | `tableau-surface/noise/bare/depolarize2_many/{side}` |
| improvement | 0.695× | 3 | 208.280 ns | 144.540 ns | `tableau-surface/noise/bare/pauli_error/{side}` |
| improvement | 0.879× | 3 | 348.650 ns | 306.400 ns | `tableau-surface/noise/bare/pauli_error_many/{side}` |
| improvement | 0.699× | 3 | 208.500 ns | 145.330 ns | `tableau-surface/noise/bare/two_qubit_pauli_error/{side}` |
| improvement | 0.858× | 3 | 295.060 ns | 252.850 ns | `tableau-surface/noise/bare/two_qubit_pauli_error_many/{side}` |
| improvement | 0.687× | 3 | 210.390 ns | 144.450 ns | `tableau-surface/noise/bare/x_error/{side}` |
| improvement | 0.781× | 3 | 242.110 ns | 188.970 ns | `tableau-surface/noise/bare/x_error_many/{side}` |
| improvement | 0.704× | 3 | 225.310 ns | 158.340 ns | `tableau-surface/noise/bare/y_error/{side}` |
| improvement | 0.814× | 3 | 267.020 ns | 216.390 ns | `tableau-surface/noise/bare/y_error_many/{side}` |
| improvement | 0.684× | 3 | 211.050 ns | 144.420 ns | `tableau-surface/noise/bare/z_error/{side}` |
| improvement | 0.780× | 3 | 242.970 ns | 189.520 ns | `tableau-surface/noise/bare/z_error_many/{side}` |
| confirmed regression | 1.089× | 3 | 5.541 µs | 6.023 µs | `tableau-surface/noise/generalized/asymmetric_loss_channel/{side}` |
| confirmed regression | 1.046× | 3 | 5.727 µs | 6.022 µs | `tableau-surface/noise/generalized/correlated_loss_channel/{side}` |
| improvement | 0.704× | 3 | 207.250 ns | 145.940 ns | `tableau-surface/noise/generalized/depolarize1/{side}` |
| improvement | 0.790× | 3 | 247.400 ns | 196.350 ns | `tableau-surface/noise/generalized/depolarize1_many/{side}` |
| improvement | 0.705× | 3 | 208.560 ns | 147.160 ns | `tableau-surface/noise/generalized/depolarize2/{side}` |
| improvement | 0.880× | 3 | 296.610 ns | 261.380 ns | `tableau-surface/noise/generalized/depolarize2_many/{side}` |
| confirmed regression | 1.148× | 3 | 3.999 µs | 4.642 µs | `tableau-surface/noise/generalized/loss_channel/{side}` |
| improvement | 0.706× | 3 | 207.130 ns | 146.190 ns | `tableau-surface/noise/generalized/pauli_error/{side}` |
| improvement | 0.884× | 3 | 342.560 ns | 302.820 ns | `tableau-surface/noise/generalized/pauli_error_many/{side}` |
| confirmed regression | 1.187× | 7 | 4.479 ns | 5.341 ns | `tableau-surface/noise/generalized/reset_loss_channel/{side}` |
| improvement | 0.703× | 3 | 207.900 ns | 146.090 ns | `tableau-surface/noise/generalized/two_qubit_pauli_error/{side}` |
| improvement | 0.877× | 3 | 293.550 ns | 257.310 ns | `tableau-surface/noise/generalized/two_qubit_pauli_error_many/{side}` |
| improvement | 0.689× | 3 | 209.760 ns | 144.980 ns | `tableau-surface/noise/generalized/x_error/{side}` |
| improvement | 0.771× | 3 | 240.220 ns | 185.150 ns | `tableau-surface/noise/generalized/x_error_many/{side}` |
| improvement | 0.698× | 3 | 226.170 ns | 157.970 ns | `tableau-surface/noise/generalized/y_error/{side}` |
| improvement | 0.808× | 3 | 268.290 ns | 217.220 ns | `tableau-surface/noise/generalized/y_error_many/{side}` |
| improvement | 0.689× | 3 | 210.220 ns | 145.020 ns | `tableau-surface/noise/generalized/z_error/{side}` |
| improvement | 0.779× | 3 | 243.940 ns | 189.020 ns | `tableau-surface/noise/generalized/z_error_many/{side}` |
| parity | 0.971× | 3 | 7.128 ns | 6.971 ns | `tableau-surface/observation/generalized/append_measurement_record/{side}` |
| confirmed regression | 3.881× | 7 | 2.016 ns | 7.987 ns | `tableau-surface/observation/generalized/bernoulli/{side}` |
| improvement | 0.958× | 3 | 235.580 ns | 227.710 ns | `tableau-surface/observation/generalized/compute_decomposition/{side}` |
| improvement | 0.663× | 3 | 0.422 ns | 0.280 ns | `tableau-surface/observation/generalized/current_measurement_record/{side}` |
| parity | 0.972× | 3 | 1.593 µs | 1.550 µs | `tableau-surface/observation/generalized/expectation/{side}` |
| confirmed regression | 3.948× | 7 | 2.012 ns | 7.952 ns | `tableau-surface/observation/generalized/flip_with_prob/{side}` |
| improvement | 0.663× | 3 | 0.423 ns | 0.281 ns | `tableau-surface/observation/generalized/n_qubits/{side}` |
| parity | 1.028× | 7 | 40.632 ns | 40.974 ns | `tableau-surface/observation/generalized/odd_phase_destabilizer_mask/{side}` |
| confirmed regression | 2.560× | 7 | 1.733 ns | 4.522 ns | `tableau-surface/observation/generalized/overwrite_last_measurement_record/{side}` |
| improvement | 0.846× | 3 | 32.363 µs | 27.556 µs | `tableau-surface/observation/generalized/trace-pattern/{side}` |
| parity | 0.991× | 3 | 1.561 µs | 1.547 µs | `tableau-surface/observation/generalized/z_expectation/{side}` |
| parity | 1.001× | 2 | 22.743 ns | 22.762 ns | `tableau-surface/projection/compute_overlap_case_a/{side}` |
| parity | 0.998× | 2 | 9.060 ns | 9.045 ns | `tableau-surface/projection/compute_overlap_case_b/{side}` |
| improvement | 0.889× | 2 | 141.895 ns | 126.225 ns | `tableau-surface/projection/project_case_a/{side}` |
| provisional regression | 1.065× | 2 | 35.542 ns | 37.908 ns | `tableau-surface/projection/project_case_b/{side}` |
| parity | 1.003× | 3 | 5.484 µs | 5.450 µs | `tableau-surface/rotation/one/r_xy/{side}` |
| parity | 0.997× | 3 | 1.564 µs | 1.559 µs | `tableau-surface/rotation/one/rotate_1_x/{side}` |
| parity | 0.998× | 3 | 1.562 µs | 1.560 µs | `tableau-surface/rotation/one/rx/{side}` |
| parity | 1.002× | 3 | 5.313 µs | 5.323 µs | `tableau-surface/rotation/one/rx_many/{side}` |
| parity | 1.001× | 3 | 2.220 µs | 2.226 µs | `tableau-surface/rotation/one/ry/{side}` |
| parity | 0.997× | 3 | 8.378 µs | 8.362 µs | `tableau-surface/rotation/one/ry_many/{side}` |
| parity | 1.001× | 3 | 2.215 µs | 2.218 µs | `tableau-surface/rotation/one/rz/{side}` |
| parity | 0.996× | 3 | 8.453 µs | 8.312 µs | `tableau-surface/rotation/one/rz_many/{side}` |
| parity | 1.004× | 3 | 6.102 µs | 6.125 µs | `tableau-surface/rotation/one/u3/{side}` |
| parity | 1.006× | 3 | 2.227 µs | 2.239 µs | `tableau-surface/rotation/t/{side}` |
| parity | 1.005× | 3 | 2.043 µs | 2.058 µs | `tableau-surface/rotation/t_dag/{side}` |
| parity | 0.993× | 3 | 7.812 µs | 7.713 µs | `tableau-surface/rotation/t_dag_many/{side}` |
| parity | 0.995× | 3 | 7.937 µs | 7.901 µs | `tableau-surface/rotation/t_many/{side}` |
| improvement | 0.191× | 3 | 53.898 µs | 10.309 µs | `tableau-surface/rotation/two-batch/rxx_many/{side}` |
| improvement | 0.289× | 3 | 54.491 µs | 15.757 µs | `tableau-surface/rotation/two-batch/rxy_many/{side}` |
| improvement | 0.319× | 3 | 47.493 µs | 15.170 µs | `tableau-surface/rotation/two-batch/rxz_many/{side}` |
| improvement | 0.320× | 3 | 53.402 µs | 17.098 µs | `tableau-surface/rotation/two-batch/ryx_many/{side}` |
| improvement | 0.297× | 3 | 53.932 µs | 15.981 µs | `tableau-surface/rotation/two-batch/ryy_many/{side}` |
| improvement | 0.329× | 3 | 47.204 µs | 15.685 µs | `tableau-surface/rotation/two-batch/ryz_many/{side}` |
| improvement | 0.316× | 3 | 53.456 µs | 16.936 µs | `tableau-surface/rotation/two-batch/rzx_many/{side}` |
| improvement | 0.297× | 3 | 53.398 µs | 15.900 µs | `tableau-surface/rotation/two-batch/rzy_many/{side}` |
| improvement | 0.324× | 3 | 47.348 µs | 15.355 µs | `tableau-surface/rotation/two-batch/rzz_many/{side}` |
| improvement | 0.387× | 3 | 11.960 µs | 4.636 µs | `tableau-surface/rotation/two/rotate_2_xz/{side}` |
| improvement | 0.221× | 3 | 11.967 µs | 2.640 µs | `tableau-surface/rotation/two/rxx/{side}` |
| improvement | 0.382× | 3 | 11.959 µs | 4.566 µs | `tableau-surface/rotation/two/rxy/{side}` |
| improvement | 0.385× | 3 | 12.069 µs | 4.641 µs | `tableau-surface/rotation/two/rxz/{side}` |
| improvement | 0.345× | 3 | 12.294 µs | 4.245 µs | `tableau-surface/rotation/two/ryx/{side}` |
| improvement | 0.373× | 3 | 12.098 µs | 4.537 µs | `tableau-surface/rotation/two/ryy/{side}` |
| improvement | 0.372× | 3 | 11.847 µs | 4.420 µs | `tableau-surface/rotation/two/ryz/{side}` |
| improvement | 0.349× | 3 | 11.819 µs | 4.148 µs | `tableau-surface/rotation/two/rzx/{side}` |
| improvement | 0.374× | 3 | 11.933 µs | 4.433 µs | `tableau-surface/rotation/two/rzy/{side}` |
| improvement | 0.372× | 3 | 12.057 µs | 4.487 µs | `tableau-surface/rotation/two/rzz/{side}` |
| parity | 1.000× | 3 | 106.540 ns | 105.900 ns | `tableau-surface/sparse-amplitudes/add_or_insert/hit/{side}` |
| parity | 1.009× | 3 | 272.530 ns | 273.680 ns | `tableau-surface/sparse-amplitudes/add_or_insert/miss/{side}` |
| parity | 0.984× | 1 | 60.509 ns | 59.538 ns | `tableau-surface/sparse-amplitudes/clone/{side}` |
| parity | 1.000× | 1 | 1.053 ns | 1.053 ns | `tableau-surface/sparse-amplitudes/default/{side}` |
| parity | 0.998× | 1 | 63.136 ns | 63.033 ns | `tableau-surface/sparse-amplitudes/entries-traversal/{side}` |
| parity | 0.998× | 1 | 69.539 ns | 69.423 ns | `tableau-surface/sparse-amplitudes/equality/{side}` |
| parity | 1.000× | 3 | 16.964 ns | 16.915 ns | `tableau-surface/sparse-amplitudes/get/hit/{side}` |
| parity | 1.000× | 3 | 36.895 ns | 37.155 ns | `tableau-surface/sparse-amplitudes/get/miss/{side}` |
| parity | 0.995× | 1 | 131.420 ns | 130.740 ns | `tableau-surface/sparse-amplitudes/into_iter/{side}` |
| parity | 0.997× | 1 | 0.456 ns | 0.454 ns | `tableau-surface/sparse-amplitudes/is_empty/{side}` |
| parity | 1.004× | 1 | 63.141 ns | 63.410 ns | `tableau-surface/sparse-amplitudes/iter-traversal/{side}` |
| parity | 1.000× | 1 | 0.455 ns | 0.455 ns | `tableau-surface/sparse-amplitudes/len/{side}` |
| parity | 1.002× | 3 | 82.728 ns | 83.466 ns | `tableau-surface/sparse-amplitudes/mul_by/{side}` |
| parity | 1.000× | 3 | 110.040 ns | 109.550 ns | `tableau-surface/sparse-amplitudes/mul_element_by/hit/{side}` |
| parity | 1.000× | 3 | 82.228 ns | 80.978 ns | `tableau-surface/sparse-amplitudes/mul_element_by/miss/{side}` |
| parity | 0.995× | 1 | 1.055 ns | 1.050 ns | `tableau-surface/sparse-amplitudes/new/{side}` |
| parity | 0.993× | 3 | 228.100 ns | 226.600 ns | `tableau-surface/sparse-amplitudes/normalize/{side}` |
| parity | 1.003× | 3 | 366.160 ns | 364.250 ns | `tableau-surface/sparse-amplitudes/reserve/{side}` |
| parity | 1.014× | 3 | 134.450 ns | 136.210 ns | `tableau-surface/sparse-amplitudes/retain/{side}` |
| parity | 1.000× | 3 | 97.823 ns | 98.994 ns | `tableau-surface/sparse-amplitudes/trim/{side}` |
| parity | 1.025× | 3 | 188.940 ns | 188.970 ns | `tableau-surface/sparse-amplitudes/unsafe_insert/{side}` |
| provisional regression | 2.202× | 1 | 3.763 ns | 8.284 ns | `word_surface/clone_copy/256/lossy/{side}/clone_cold` |
| provisional regression | 2.213× | 1 | 3.783 ns | 8.371 ns | `word_surface/clone_copy/256/lossy/{side}/clone_warm` |
| provisional regression | 1.962× | 1 | 2.835 ns | 5.561 ns | `word_surface/clone_copy/256/ordinary/{side}/clone_cold` |
| provisional regression | 1.924× | 1 | 2.869 ns | 5.520 ns | `word_surface/clone_copy/256/ordinary/{side}/clone_warm` |
| provisional regression | 1.801× | 1 | 3.227 ns | 5.813 ns | `word_surface/clone_copy/256/phased/{side}/clone_cold` |
| provisional regression | 1.787× | 1 | 3.258 ns | 5.821 ns | `word_surface/clone_copy/256/phased/{side}/clone_warm` |
| confirmed regression | 1.160× | 7 | 3.900 ns | 4.487 ns | `word_surface/construct/lossy/{side}/new_identity/256` |
| confirmed regression | 1.163× | 7 | 3.904 ns | 4.472 ns | `word_surface/construct/lossy/{side}/new_identity/64` |
| confirmed regression | 1.164× | 7 | 3.885 ns | 4.481 ns | `word_surface/construct/lossy/{side}/new_identity/8` |
| improvement | 0.404× | 7 | 488.070 ns | 192.860 ns | `word_surface/construct/lossy/{side}/parse/256` |
| improvement | 0.398× | 7 | 132.930 ns | 52.408 ns | `word_surface/construct/lossy/{side}/parse/64` |
| improvement | 0.345× | 7 | 27.869 ns | 9.621 ns | `word_surface/construct/lossy/{side}/parse/8` |
| parity | 1.002× | 3 | 2.736 ns | 2.740 ns | `word_surface/construct/ordinary/{side}/new_identity/256` |
| parity | 0.995× | 3 | 2.740 ns | 2.731 ns | `word_surface/construct/ordinary/{side}/new_identity/64` |
| parity | 0.999× | 3 | 2.735 ns | 2.731 ns | `word_surface/construct/ordinary/{side}/new_identity/8` |
| improvement | 0.405× | 3 | 401.760 ns | 164.610 ns | `word_surface/construct/ordinary/{side}/parse/256` |
| improvement | 0.421× | 3 | 107.080 ns | 45.240 ns | `word_surface/construct/ordinary/{side}/parse/64` |
| improvement | 0.366× | 3 | 21.709 ns | 7.996 ns | `word_surface/construct/ordinary/{side}/parse/8` |
| parity | 1.000× | 3 | 3.055 ns | 3.059 ns | `word_surface/construct/phased/{side}/new_identity/256` |
| parity | 1.002× | 3 | 3.054 ns | 3.057 ns | `word_surface/construct/phased/{side}/new_identity/64` |
| parity | 1.000× | 3 | 3.047 ns | 3.048 ns | `word_surface/construct/phased/{side}/new_identity/8` |
| improvement | 0.609× | 3 | 679.430 ns | 413.480 ns | `word_surface/construct/phased/{side}/parse/256` |
| improvement | 0.725× | 3 | 242.030 ns | 177.210 ns | `word_surface/construct/phased/{side}/parse/64` |
| improvement | 0.617× | 3 | 41.686 ns | 25.708 ns | `word_surface/construct/phased/{side}/parse/8` |
| parity | 1.008× | 1 | 4.943 ns | 4.984 ns | `word_surface/construct/phased_explicit_phase/256/{side}/from_word_and_phase` |
| parity | 1.003× | 2 | 4.588 ns | 4.604 ns | `word_surface/construct/phased_from_word/256/{side}/from_existing_word` |
| confirmed regression | 1.799× | 7 | 6.818 ns | 12.298 ns | `word_surface/lossy/branch_key/256/{side}/one_site_clone_then_bits` |
| confirmed regression | 1.598× | 7 | 7.166 ns | 11.556 ns | `word_surface/lossy/branch_key/256/{side}/two_site_clone_then_bits` |
| improvement | 0.850× | 7 | 8.763 ns | 7.414 ns | `word_surface/lossy/clifford_lost_guard/256/{side}/h_lost_noop` |
| improvement | 0.964× | 7 | 12.105 ns | 11.661 ns | `word_surface/lossy/clifford_present/256/{side}/cnot` |
| improvement | 0.959× | 7 | 12.127 ns | 11.671 ns | `word_surface/lossy/clifford_present/256/{side}/cx_alias` |
| improvement | 0.842× | 7 | 12.240 ns | 10.287 ns | `word_surface/lossy/clifford_present/256/{side}/cy` |
| improvement | 0.821× | 7 | 10.848 ns | 8.898 ns | `word_surface/lossy/clifford_present/256/{side}/cz` |
| improvement | 0.591× | 7 | 13.948 ns | 8.259 ns | `word_surface/lossy/clifford_present/256/{side}/h` |
| improvement | 0.669× | 7 | 12.012 ns | 8.035 ns | `word_surface/lossy/clifford_present/256/{side}/s` |
| improvement | 0.669× | 7 | 12.006 ns | 8.040 ns | `word_surface/lossy/clifford_present/256/{side}/s_dag` |
| improvement | 0.682× | 7 | 14.006 ns | 9.559 ns | `word_surface/lossy/clifford_present/256/{side}/sqrt_x` |
| improvement | 0.683× | 7 | 14.033 ns | 9.578 ns | `word_surface/lossy/clifford_present/256/{side}/sqrt_x_dag` |
| improvement | 0.590× | 7 | 13.996 ns | 8.266 ns | `word_surface/lossy/clifford_present/256/{side}/sqrt_y` |
| improvement | 0.591× | 7 | 13.985 ns | 8.273 ns | `word_surface/lossy/clifford_present/256/{side}/sqrt_y_dag` |
| confirmed regression | 1.135× | 7 | 3.996 ns | 4.505 ns | `word_surface/lossy/clifford_present/256/{side}/x` |
| confirmed regression | 1.122× | 7 | 4.022 ns | 4.487 ns | `word_surface/lossy/clifford_present/256/{side}/y` |
| confirmed regression | 1.140× | 7 | 3.939 ns | 4.499 ns | `word_surface/lossy/clifford_present/256/{side}/z` |
| improvement | 0.965× | 7 | 12.114 ns | 11.689 ns | `word_surface/lossy/clifford_present/256/{side}/zcx_alias` |
| improvement | 0.844× | 7 | 12.243 ns | 10.308 ns | `word_surface/lossy/clifford_present/256/{side}/zcy_alias` |
| improvement | 0.819× | 7 | 10.831 ns | 8.926 ns | `word_surface/lossy/clifford_present/256/{side}/zcz_alias` |
| provisional regression | 1.392× | 1 | 0.510 ns | 0.710 ns | `word_surface/lossy/hash_protocol/256/{side}/warm` |
| improvement | 0.552× | 1 | 13.856 ns | 7.651 ns | `word_surface/lossy/mutate/256/{side}/clear_loss` |
| improvement | 0.612× | 1 | 13.847 ns | 8.469 ns | `word_surface/lossy/mutate/256/{side}/set_lost` |
| improvement | 0.599× | 1 | 13.672 ns | 8.186 ns | `word_surface/lossy/mutate/256/{side}/set_present` |
| provisional regression | 1.558× | 1 | 6.600 ns | 10.283 ns | `word_surface/lossy/mutate/256/{side}/set_x_bit` |
| provisional regression | 1.541× | 1 | 6.709 ns | 10.337 ns | `word_surface/lossy/mutate/256/{side}/set_z_bit` |
| improvement | 0.591× | 1 | 1.335 µs | 788.950 ns | `word_surface/lossy/observation/256/{side}/display` |
| provisional regression | 1.145× | 1 | 1.525 ns | 1.746 ns | `word_surface/lossy/observation/256/{side}/equality` |
| confirmed regression | 1.230× | 7 | 1.027 ns | 1.261 ns | `word_surface/lossy/read/256/{side}/get` |
| parity | 1.012× | 1 | 0.917 ns | 0.928 ns | `word_surface/lossy/read/256/{side}/is_lost` |
| parity | 1.008× | 2 | 155.915 ns | 157.185 ns | `word_surface/lossy/read/256/{side}/iter_traverse` |
| parity | 0.999× | 7 | 0.708 ns | 0.707 ns | `word_surface/lossy/read/256/{side}/loss_weight` |
| parity | 1.003× | 7 | 0.874 ns | 0.879 ns | `word_surface/lossy/read/256/{side}/weight` |
| parity | 0.994× | 2 | 0.473 ns | 0.470 ns | `word_surface/lossy/read/256/{side}/width` |
| parity | 0.997× | 1 | 0.943 ns | 0.940 ns | `word_surface/lossy/read/256/{side}/x_bit` |
| parity | 1.002× | 1 | 0.919 ns | 0.921 ns | `word_surface/lossy/read/256/{side}/z_bit` |
| parity | 0.988× | 1 | 4.252 ns | 4.199 ns | `word_surface/ordinary/branch_key/256/{side}/one_site` |
| improvement | 0.963× | 1 | 4.918 ns | 4.735 ns | `word_surface/ordinary/branch_key/256/{side}/two_site` |
| improvement | 0.931× | 3 | 7.593 ns | 7.070 ns | `word_surface/ordinary/clifford/256/{side}/cnot` |
| improvement | 0.930× | 3 | 7.671 ns | 7.101 ns | `word_surface/ordinary/clifford/256/{side}/cx_alias` |
| confirmed regression | 1.035× | 7 | 8.068 ns | 8.363 ns | `word_surface/ordinary/clifford/256/{side}/cy` |
| improvement | 0.817× | 3 | 7.358 ns | 6.014 ns | `word_surface/ordinary/clifford/256/{side}/cz` |
| improvement | 0.621× | 3 | 8.267 ns | 5.086 ns | `word_surface/ordinary/clifford/256/{side}/h` |
| improvement | 0.687× | 3 | 7.614 ns | 5.188 ns | `word_surface/ordinary/clifford/256/{side}/s` |
| improvement | 0.682× | 3 | 7.687 ns | 5.151 ns | `word_surface/ordinary/clifford/256/{side}/s_dag` |
| improvement | 0.916× | 3 | 7.977 ns | 7.329 ns | `word_surface/ordinary/clifford/256/{side}/sqrt_x` |
| improvement | 0.916× | 3 | 8.001 ns | 7.386 ns | `word_surface/ordinary/clifford/256/{side}/sqrt_x_dag` |
| improvement | 0.619× | 3 | 8.143 ns | 5.044 ns | `word_surface/ordinary/clifford/256/{side}/sqrt_y` |
| improvement | 0.612× | 3 | 8.235 ns | 5.012 ns | `word_surface/ordinary/clifford/256/{side}/sqrt_y_dag` |
| parity | 1.005× | 3 | 2.599 ns | 2.664 ns | `word_surface/ordinary/clifford/256/{side}/x` |
| parity | 1.016× | 3 | 2.587 ns | 2.695 ns | `word_surface/ordinary/clifford/256/{side}/y` |
| improvement | 0.969× | 3 | 2.670 ns | 2.600 ns | `word_surface/ordinary/clifford/256/{side}/z` |
| improvement | 0.926× | 3 | 7.660 ns | 7.065 ns | `word_surface/ordinary/clifford/256/{side}/zcx_alias` |
| confirmed regression | 1.031× | 7 | 8.102 ns | 8.369 ns | `word_surface/ordinary/clifford/256/{side}/zcy_alias` |
| improvement | 0.829× | 3 | 7.386 ns | 6.077 ns | `word_surface/ordinary/clifford/256/{side}/zcz_alias` |
| parity | 1.015× | 1 | 0.508 ns | 0.515 ns | `word_surface/ordinary/hash_protocol/256/{side}/warm` |
| provisional regression | 1.067× | 1 | 4.974 ns | 5.305 ns | `word_surface/ordinary/mutate/256/{side}/set_x_bit` |
| provisional regression | 1.082× | 1 | 5.073 ns | 5.487 ns | `word_surface/ordinary/mutate/256/{side}/set_z_bit` |
| improvement | 0.616× | 1 | 1.218 µs | 749.440 ns | `word_surface/ordinary/observation/256/{side}/display` |
| provisional regression | 1.155× | 1 | 1.174 ns | 1.356 ns | `word_surface/ordinary/observation/256/{side}/equality` |
| improvement | 0.326× | 1 | 11.536 ns | 3.766 ns | `word_surface/ordinary/product/256/{side}/product` |
| confirmed regression | 1.288× | 3 | 0.941 ns | 1.212 ns | `word_surface/ordinary/read/256/{side}/get` |
| provisional regression | 1.060× | 2 | 137.810 ns | 146.115 ns | `word_surface/ordinary/read/256/{side}/iter_traverse` |
| parity | 0.997× | 1 | 0.533 ns | 0.531 ns | `word_surface/ordinary/read/256/{side}/loss_weight` |
| confirmed regression | 1.033× | 3 | 0.796 ns | 0.822 ns | `word_surface/ordinary/read/256/{side}/weight` |
| parity | 0.996× | 2 | 0.472 ns | 0.470 ns | `word_surface/ordinary/read/256/{side}/width` |
| parity | 1.025× | 1 | 0.916 ns | 0.938 ns | `word_surface/ordinary/read/256/{side}/x_bit` |
| parity | 0.986× | 1 | 0.933 ns | 0.920 ns | `word_surface/ordinary/read/256/{side}/z_bit` |
| improvement | 0.371× | 2 | 3.607 µs | 1.337 µs | `word_surface/pattern/bounded_enumeration/8/{side}/enumerate_all` |
| provisional regression | 1.774× | 2 | 309.025 ns | 548.515 ns | `word_surface/pattern/match_contains/256/{side}/lossy_present` |
| provisional regression | 1.902× | 2 | 256.305 ns | 487.575 ns | `word_surface/pattern/match_contains/256/{side}/ordinary` |
| parity | 0.992× | 1 | 84.062 ns | 83.399 ns | `word_surface/pattern/observation/{side}/display` |
| parity | 0.995× | 1 | 3.556 ns | 3.537 ns | `word_surface/pattern/observation/{side}/equality` |
| provisional regression | 1.163× | 2 | 63.692 ns | 74.065 ns | `word_surface/pattern/parse/{side}/indexed` |
| improvement | 0.617× | 2 | 61.591 ns | 38.007 ns | `word_surface/pattern/parse/{side}/optional_repeat` |
| improvement | 0.701× | 2 | 56.320 ns | 39.497 ns | `word_surface/pattern/parse/{side}/star` |
| parity | 0.992× | 1 | 3.669 ns | 3.638 ns | `word_surface/phased/add_phase/256/{side}/add_phase` |
| improvement | 0.808× | 3 | 8.540 ns | 6.900 ns | `word_surface/phased/clifford/256/{side}/cnot` |
| improvement | 0.800× | 3 | 8.488 ns | 6.794 ns | `word_surface/phased/clifford/256/{side}/cx_alias` |
| improvement | 0.846× | 3 | 8.823 ns | 7.599 ns | `word_surface/phased/clifford/256/{side}/cy` |
| improvement | 0.819× | 3 | 8.182 ns | 6.693 ns | `word_surface/phased/clifford/256/{side}/cz` |
| improvement | 0.736× | 3 | 8.299 ns | 6.003 ns | `word_surface/phased/clifford/256/{side}/h` |
| improvement | 0.674× | 3 | 8.740 ns | 5.890 ns | `word_surface/phased/clifford/256/{side}/s` |
| improvement | 0.883× | 3 | 7.286 ns | 6.393 ns | `word_surface/phased/clifford/256/{side}/s_dag` |
| improvement | 0.797× | 3 | 8.004 ns | 6.313 ns | `word_surface/phased/clifford/256/{side}/sqrt_x` |
| improvement | 0.744× | 3 | 7.909 ns | 5.818 ns | `word_surface/phased/clifford/256/{side}/sqrt_x_dag` |
| improvement | 0.684× | 3 | 8.302 ns | 5.678 ns | `word_surface/phased/clifford/256/{side}/sqrt_y` |
| improvement | 0.670× | 3 | 8.107 ns | 5.435 ns | `word_surface/phased/clifford/256/{side}/sqrt_y_dag` |
| parity | 0.988× | 3 | 5.553 ns | 5.486 ns | `word_surface/phased/clifford/256/{side}/x` |
| improvement | 0.895× | 3 | 5.565 ns | 5.065 ns | `word_surface/phased/clifford/256/{side}/y` |
| parity | 0.988× | 3 | 5.820 ns | 5.441 ns | `word_surface/phased/clifford/256/{side}/z` |
| improvement | 0.794× | 3 | 8.826 ns | 6.864 ns | `word_surface/phased/clifford/256/{side}/zcx_alias` |
| improvement | 0.850× | 3 | 8.791 ns | 7.704 ns | `word_surface/phased/clifford/256/{side}/zcy_alias` |
| improvement | 0.824× | 3 | 8.208 ns | 6.743 ns | `word_surface/phased/clifford/256/{side}/zcz_alias` |
| improvement | 0.635× | 1 | 1.274 µs | 809.670 ns | `word_surface/phased/observation/256/{side}/display` |
| provisional regression | 1.098× | 1 | 1.361 ns | 1.493 ns | `word_surface/phased/observation/256/{side}/equality` |
| improvement | 0.332× | 1 | 11.591 ns | 3.846 ns | `word_surface/phased/product/256/{side}/product` |
| confirmed regression | 1.284× | 3 | 0.944 ns | 1.211 ns | `word_surface/phased/read/256/{side}/get` |
| improvement | 0.925× | 3 | 0.514 ns | 0.477 ns | `word_surface/phased/read/256/{side}/is_positive` |
| provisional regression | 1.051× | 2 | 139.725 ns | 146.830 ns | `word_surface/phased/read/256/{side}/iter_traverse` |
| parity | 1.001× | 1 | 0.484 ns | 0.485 ns | `word_surface/phased/read/256/{side}/phase` |
| improvement | 0.969× | 3 | 0.827 ns | 0.802 ns | `word_surface/phased/read/256/{side}/weight_delegate` |
| parity | 0.999× | 2 | 0.472 ns | 0.472 ns | `word_surface/phased/read/256/{side}/width` |

## No-old-twin measurements and exclusions

- New-only Pauli-sum `reduce`, Hermitian overlap, and sum×sum multiplication:
  old sum-RHS multiplication is uninstantiable even for a singleton RHS;
  `multiply_into` has no old callable semantic twin.
- Mixed/non-unit projection is new-only because Lean adjudicated the old behavior
  as incorrect; only the common I/Z unit-coefficient subset is compared.
- Lossy Pauli words intentionally have no native product.
- Complex symbolic evaluation and exact Gaussian-ring operations are new-only.
- Symbolic projection is blocked because the new symbolic term deliberately
  lacks `Halvable`; symbolic amplitude damping lacks `Float`; neither engine
  implements `PauliErrorAll` for symbolic sums.
- Direct tableau sampling has no common old/new API; mixture sampling is compared.

Unpaired measured benchmark IDs in the screening output: **126**.
