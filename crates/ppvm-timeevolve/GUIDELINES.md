# ppvm-timeevolve Development Guidelines

## 1. Code Simplicity and Reuse

Less is more. Prefer reusing existing types and operations from `ppvm-runtime` over writing
new infrastructure. Avoid premature abstraction — if something is only used once, don't
generalize it. A few clear lines beat a clever helper.

## 2. Efficiency

Performance matters, but not at the cost of unreasonable complexity. Prefer approaches that
are fast in the common case and readable. Avoid micro-optimizations that obscure intent.

## 3. No Changes to Other Crates

Do not modify `ppvm-runtime` or any other existing crate. If a required capability is
genuinely absent and cannot be worked around without unreasonable code, raise the issue
explicitly with the maintainer before proceeding. Do not make the decision unilaterally.

## 4. Use Commits as Checkpoints

Commit after each self-contained, working piece of implementation. Each commit should leave
the crate in a buildable state. This makes it easy to review progress and roll back if
needed.

## 5. Changes Must Be Verifiable

Every non-trivial piece of logic must have unit tests before it is considered done. Tests
should cover the happy path and at least one edge case. Do not move on to the next task
until tests pass and the change has been reviewed.
