# edge_cases/

Hand-written corpus of corner cases. No regen tooling required.

Fixture inventory:

- Empty / whitespace-only / comment-only programs.
- Every `pi_expr` shape: bare `pi`, `<coeff>*pi`, plain f64.
- Every tag shape: bare ident, single named, multiple named, multiple positional, mixed.
- REPEAT shapes: single-line, multi-line, three-deep nested, single-instruction body.
- Annotation `rec[-k]` targets (parser tolerates and discards for annotations).
- Dense measurement (`M 0 1 ... 31`).
- Sparse measurement (single-target `M`s interleaved with gates).
- One fixture exercising every Phase-1 gate family.

Distribution-mode fixtures use `num_shots=256`. Bootstrapped at fixture-creation time
by running ppvm and recording per-bit means.

Migrating from the previous flat layout, fixtures `bell_pair.stim`, `ghz.stim`,
`x_only.stim`, `repeat_block.stim`, `repetition_code_d3_r3.stim`, and
`depolarize_smoke.stim` were relocated here.
