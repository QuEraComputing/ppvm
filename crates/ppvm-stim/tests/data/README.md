# Stim corpus

Curated `.stim` fixtures used by `tests/stim_corpus.rs`.

The harness asserts:

1. **Parse**: every file under this directory parses (or matches its
   declared parse failure, including `ExtendedParseError::InvalidTag` for
   invalid extended-dialect tags).
2. **Prepare/execute**: every parsable file either executes successfully
   **or** fails with the specific `ExecError::Unsupported` variant declared
   by the harness table.

Phase-2 work converts each "expected Unsupported" entry into "expected to
execute" as features land — free regression coverage.

## Provenance

- `ghz.stim`, `x_only.stim`, `bell_pair.stim`, `repeat_block.stim`,
  `depolarize_smoke.stim`, `swap_unsupported.stim`, `mx_unsupported.stim`,
  `repetition_code_d3_r3.stim`, `feedback_cx_unsupported.stim`:
  hand-written by the ppvm team.
- Future fixtures pulled from `quantumlib/Stim` should record the upstream
  commit SHA in this file when added.
