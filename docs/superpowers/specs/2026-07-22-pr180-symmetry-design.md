# Translation-Symmetry Module Design

## Status

This design applies to the translation-symmetry functionality introduced by
PR 180. Implementation should occur on a branch based on
`split/1-translation-symmetry`; follow-up PRs remain responsible for their own
Lindbladian, Python-binding, and benchmarking changes.

## Goals

- Split `symmetry.rs` along semantic boundaries without changing the public
  `ppvm_pauli_sum::symmetry::*` import surface.
- Make `TranslationGroup` reject malformed generator descriptions eagerly.
- Give real-space merging and momentum projection explicit, distinct
  semantics.
- Correct momentum projection for words with nontrivial stabilizers.
- Make momentum-sector validation detect missing orbit members and
  stabilizer-incompatible sectors.
- Centralize group traversal so canonicalization, shift recovery, and orbit
  iteration cannot diverge.

## Non-goals

- Moving functionality from PRs 181–183 into PR 180.
- Adding Python bindings or changing Python APIs; those belong to PR 182.
- Splitting individual lattice constructors into separate files.
- Generalizing beyond finite abelian permutation actions.
- Introducing a new fallible public constructor API in this change.

## Module structure

Replace the single source file with the following module:

```text
crates/ppvm-pauli-sum/src/symmetry/
├── mod.rs       # module documentation and public re-exports
├── group.rs     # TranslationGroup, construction, traversal, canonicalization
├── merge.rs     # unnormalized k=0 Vec and PauliSum merging
├── momentum.rs  # characters, normalized projection, sector errors/checking
└── tests.rs     # unit and regression tests for the complete module
```

`lib.rs` continues to declare `pub mod symmetry`. `mod.rs` re-exports the
existing public names, so consumers do not need to change imports. Internal
traversal helpers are `pub(super)` at most.

`canonicalize_with_shift` remains in `group.rs`: although momentum code is its
current consumer, it is fundamentally an orbit operation. `character` belongs
in `momentum.rs` and may be implemented as an additional inherent
`TranslationGroup` implementation using the public generator accessors.

## Group model and validation

`TranslationGroup` represents an abstract product
`G = C_orders[0] × ... × C_orders[r-1]` together with a commuting permutation
action on qubits. The action may have a kernel: different elements of `G` may
produce the same qubit permutation or fix a particular Pauli word. Consequently,
`order()` returns the order of the abstract parameter group, while a word's
distinct orbit can be smaller.

This definition is preferred over requiring a faithful action. Proving
faithfulness by enumerating every composed permutation costs
`O(|G| × n_qubits)` during construction and would reject useful descriptions
that the stabilizer-aware projection can handle correctly. Two alternatives
were considered and rejected:

1. Require all generated permutations to be unique. This makes `order()` equal
   the image-group order but introduces potentially prohibitive construction
   work.
2. Support only the built-in lattice constructors. This avoids arbitrary-group
   validation but removes an intentional public extension point.

`from_generators` retains its current infallible signature and validates:

- equal numbers of permutations and orders;
- nonzero orders;
- permutation length, range, and uniqueness;
- `perm^order == identity` for every generator;
- pairwise commutativity of the generator permutations; and
- checked multiplication of generator orders into a cached `usize` group
  order.

The built-in constructors additionally reject zero dimensions, use checked
dimension products, and reject values that cannot be represented by their
`u32` permutation/order storage. Assertions must have precise messages. A
future fallible constructor can be added separately if callers need to recover
from invalid input.

All public operations perform release-mode shape validation. In particular,
word width must match `n_qubits`, and `character` requires one momentum mode
and one counter per generator.

## Group traversal

`group.rs` owns one internal mixed-radix traversal primitive. It maintains the
current counter and Pauli word as an odometer. Incrementing a digit applies its
generator once; rollover applies the same generator once more, returning that
digit to the identity before carrying. Order-one generators are skipped.

The traversal yields `(word, counter)` for every element of the abstract group,
including duplicate words caused by stabilizers. It uses `usize` for radix
arithmetic and casts only a checked per-digit remainder to `u32`.

`canonicalize`, `canonicalize_with_shift`, and `orbit` all consume this
primitive. `canonicalize_with_shift` returns the inverse counter taking the
chosen representative back to the input word. When several counters describe
the same mapping, it returns the first in traversal order; momentum code must
not assume that this counter is unique.

With order-one generators removed from carry processing, traversal performs an
amortized constant number of generator applications per group element, giving
the documented `O(|G| × n_qubits)` time and `O(r)` counter state.

## Merge semantics

The APIs intentionally expose two different coefficient conventions:

- `canonicalize_pauli_sum` and `symmetry_merge_pauli_sum` perform an
  **unnormalized merge**. Coefficients of input entries with the same orbit
  representative are summed.
- `canonicalize_pauli_sum_complex` performs a **normalized momentum
  projection**. Its output coefficient is the projected state's coefficient
  on the chosen representative.

The complex projection first coalesces duplicate input Pauli entries. For each
represented orbit it enumerates the abstract group action and determines:

1. the distinct orbit members and their character phases; and
2. whether the requested character is trivial on the representative's
   stabilizer.

If the character is nontrivial on the stabilizer, the group projector vanishes
on that orbit and the orbit is omitted from the result. Otherwise the projected
representative coefficient is the average of the phase-adjusted coefficients
over the **distinct orbit**, with absent input members contributing zero and
normalization by the orbit size rather than by `|G|`.
This is equivalent to the full `1/|G|` group projector because every distinct
orbit member occurs once per stabilizer element.

For `k = 0`, the complex routine is an orbit average, not the unnormalized real
merge. Documentation and test names must state that distinction explicitly.

## Momentum-sector validation

`check_momentum_sector` validates the represented operator, not merely pairs of
entries that happen to be present. It therefore:

1. coalesces duplicate Pauli entries;
2. chooses each unvisited orbit representative and infers its coefficient from
   the first present nonzero member;
3. enumerates all group elements acting on that representative;
4. verifies that repeated occurrences of the same word have compatible
   character phases, thereby checking stabilizer compatibility; and
5. compares every distinct orbit member's expected coefficient with its actual
   coefficient, treating a missing entry as zero.

Relative comparison continues to use `tol * max(|reference|, 1)`. Exact-zero
orbits produced by coalescing are ignored.

`SectorCheckError` gains a public `kind: SectorCheckErrorKind` field, where the
kind is either `CoefficientMismatch` or `IncompatibleStabilizer`. The error
retains the representative, offending or missing word, expected coefficient,
actual coefficient, and shift. Its `Debug` and `Display` implementations print
the relevant Pauli words; the necessary hasher bounds are placed on those
formatting implementations.

## Tests

The split must preserve all existing tests and add focused regressions for:

- zero generator orders and zero lattice dimensions;
- incorrect generator periods and noncommuting generators;
- checked group-order and lattice-dimension multiplication;
- public word-width, momentum-length, and counter-length checks;
- mixed-radix decoding without pre-modulo `u32` truncation;
- agreement of the shared traversal with brute-force composition on small 1D,
  2D, 3D, and ladder groups;
- a period-two word on a four-site chain round-tripping in `k=0`;
- the same orbit round-tripping in a compatible nonzero momentum sector;
- an incompatible momentum character vanishing on projection;
- rejection of a basis with missing orbit members;
- rejection of a nonzero coefficient with an incompatible stabilizer;
- acceptance of existing complete momentum eigenstates; and
- the documented normalization difference between real merge and complex
  projection.

Run `cargo fmt --check`, `cargo test -p ppvm-pauli-sum`, and
`cargo test --workspace` before updating PR 180.

## Stacked-PR integration

PR 180 owns the Rust core module split and all correctness changes above. Once
it is updated, rebase descendants in this order:

```text
split/1-translation-symmetry
└── split/2-ctpp-core                         # PR 181
    └── split/3-symmetric-evolution           # PR 182
        ├── split/4-autotune-ledgers           # PR 183
        └── kossakowski-dissipator             # PR 179
```

PR 181 uses only the stable real-space API and should need no source changes.
PR 182 is the first momentum consumer and must update its Rust/Python comments,
wrappers, and tests that currently describe a fixed `1/|G|` prefactor. Add a
Python period-two regression there after rebasing. PRs 183 and 179 should then
be rebased without importing their features into PR 180.

Review-thread replies should consolidate duplicates: generator validation,
`k=0` documentation, mixed-radix decoding, word-width checks, and zero-sized
constructors each have multiple Copilot threads. Resolve them only after the
corresponding implementation and regression tests are present.

## Acceptance criteria

- Existing public Rust import paths compile unchanged.
- Each new source file has one clear responsibility as described above.
- All malformed group definitions covered by the constructor contract fail
  immediately with clear messages.
- Complex projection round-trips valid eigenstates with both full and reduced
  orbits.
- Sector validation rejects missing members and incompatible stabilizers.
- Real merge remains unnormalized and behavior-compatible.
- The full Rust workspace passes after the split.
- Descendant PRs can be rebased in order without pulling their functionality
  into PR 180.
