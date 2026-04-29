# unsupported/

One fixture per phase-1-unsupported Stim instruction, plus auto-routed
fixtures from `regen-stim codes` whose generated circuit contains a
phase-1 gap (e.g. `RX`, `MX`, `C_XYZ`, `MR` in some sub-tasks).

Each fixture is recorded with `mode: "unsupported"` and Stim's reference
means pre-recorded. When phase-2 lifts a given instruction, the fixture
flips to `mode: "distribution"` (or `deterministic`).

## Phase-2 flip workflow

```bash
cd crates/ppvm-stim/tests/regen-stim
uv run regen-stim refresh ../data/unsupported/<name>.stim
```

`refresh` reads the existing JSON, sees that phase-2 now supports the
instruction (because `ppvm_stim::normalize::to_tableau` succeeds on the
source), runs the seed-search loop against the pre-recorded
`stim_bit_means`, and writes a new JSON with `mode: "distribution"`. The
`.stim` source itself never changes. After flipping, edit
`regen-stim/src/regen_stim/codes.py`'s `PHASE1_SUPPORTED` set so the next
`codes` regen no longer auto-routes circuits using the now-supported
instruction here.

## Regenerating

```bash
cd crates/ppvm-stim/tests/regen-stim
uv run regen-stim unsupported
```

Idempotent: overwriting any existing fixture with the same content is
safe. The corpus harness accepts rejection from any `parse(...)` failure
or `NormalizeError::Unsupported`; the latter checks that the rejected
instruction name matches `awaiting_phase2_instruction` exactly.
