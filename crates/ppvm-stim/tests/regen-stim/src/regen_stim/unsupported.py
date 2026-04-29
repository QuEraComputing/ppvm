"""unsupported/: one fixture per phase-1-unsupported Stim instruction.

Each fixture: prep → apply unsupported gate → measure. Stim reference is
pre-recorded so phase-2 lifting only needs ppvm-side fields added (via
`regen-stim refresh`).
"""

from __future__ import annotations

from . import core


# (instruction_name, source_template, fixture_name)
UNSUPPORTED_FIXTURES: list[tuple[str, str, str]] = [
    ("SWAP", "X 0\nSWAP 0 1\nM 0 1\n", "swap_unsupported"),
    ("ISWAP", "X 0\nISWAP 0 1\nM 0 1\n", "iswap_unsupported"),
    ("ISWAP_DAG", "X 0\nISWAP_DAG 0 1\nM 0 1\n", "iswap_dag_unsupported"),
    ("SQRT_XX", "H 0 1\nSQRT_XX 0 1\nM 0 1\n", "sqrt_xx_unsupported"),
    ("SQRT_YY", "H 0 1\nSQRT_YY 0 1\nM 0 1\n", "sqrt_yy_unsupported"),
    ("SQRT_ZZ", "H 0 1\nSQRT_ZZ 0 1\nM 0 1\n", "sqrt_zz_unsupported"),
    ("CXSWAP", "X 0\nCXSWAP 0 1\nM 0 1\n", "cxswap_unsupported"),
    ("SWAPCX", "X 0\nSWAPCX 0 1\nM 0 1\n", "swapcx_unsupported"),
    ("XCX", "X 0\nXCX 0 1\nM 0 1\n", "xcx_unsupported"),
    ("XCY", "X 0\nXCY 0 1\nM 0 1\n", "xcy_unsupported"),
    ("XCZ", "X 0\nXCZ 0 1\nM 0 1\n", "xcz_unsupported"),
    ("YCX", "X 0\nYCX 0 1\nM 0 1\n", "ycx_unsupported"),
    ("YCY", "X 0\nYCY 0 1\nM 0 1\n", "ycy_unsupported"),
    ("YCZ", "X 0\nYCZ 0 1\nM 0 1\n", "ycz_unsupported"),
    ("C_XYZ", "H 0\nC_XYZ 0\nM 0\n", "c_xyz_unsupported"),
    ("C_ZYX", "H 0\nC_ZYX 0\nM 0\n", "c_zyx_unsupported"),
    ("H_XY", "H 0\nH_XY 0\nM 0\n", "h_xy_unsupported"),
    ("H_YZ", "H 0\nH_YZ 0\nM 0\n", "h_yz_unsupported"),
    ("MX", "H 0\nMX 0\n", "mx_unsupported"),
    ("MY", "H 0\nMY 0\n", "my_unsupported"),
    ("MRX", "H 0\nMRX 0\n", "mrx_unsupported"),
    ("MRY", "H 0\nMRY 0\n", "mry_unsupported"),
    ("MXX", "H 0 1\nMXX 0 1\n", "mxx_unsupported"),
    ("MYY", "H 0 1\nMYY 0 1\n", "myy_unsupported"),
    ("MZZ", "H 0 1\nMZZ 0 1\n", "mzz_unsupported"),
    ("MPP", "H 0 1\nMPP X0*X1\n", "mpp_unsupported"),
    ("HERALDED_ERASE", "X 0\nHERALDED_ERASE(0.1) 0\nM 0\n", "heralded_erase_unsupported"),
    (
        "HERALDED_PAULI_CHANNEL_1",
        "X 0\nHERALDED_PAULI_CHANNEL_1(0.05, 0.05, 0.05, 0.05) 0\nM 0\n",
        "heralded_pauli_channel_1_unsupported",
    ),
    (
        "CORRELATED_ERROR",
        "CORRELATED_ERROR(0.1) X0 X1\nM 0 1\n",
        "correlated_error_unsupported",
    ),
    (
        "ELSE_CORRELATED_ERROR",
        "CORRELATED_ERROR(0.1) X0\nELSE_CORRELATED_ERROR(0.1) X1\nM 0 1\n",
        "else_correlated_error_unsupported",
    ),
]


def run() -> int:
    paths = core.CorpusPaths.default()
    failures: list[str] = []
    written = 0
    for instruction, source, name in UNSUPPORTED_FIXTURES:
        meta = core.FixtureMeta(
            name=name,
            category="unsupported",
            source=source,
            test_num_shots=0,
        )
        try:
            core.write_unsupported_fixture(
                meta, paths, awaiting_phase2_instruction=instruction
            )
            written += 1
        except Exception as e:
            failures.append(f"{name}: {e}")
    print(f"regen-stim unsupported: wrote {written} fixtures")
    if failures:
        print("regen-stim unsupported: failures:")
        for f in failures:
            print(f"  {f}")
        return 1
    return 0
