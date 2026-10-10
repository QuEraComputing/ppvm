# Symbolic propagation benchmark plan

User-approved scope: compare PPVM trigonometric-polynomial coefficients with
PauliPropagation.jl surrogates on parametrized TFIM circuits, measure construction
and repeated evaluation separately, and compare against repeated numeric
propagation. Publish reproducible results in a benchmark PR, then add blue paper
text and a matching figure.

1. Add Rust and Julia runners using identical forward circuits, reversed
   Heisenberg propagation, a Z observable, and independent or shared angles.
   Use exact propagation (no structural or magnitude truncation beyond floating
   point roundoff). Warm up before timing; report medians of repeated trials.
2. Validate coefficient vectors across engines and against numeric propagation at
   deterministic parameter assignments. Include analytic one-qubit checks and
   fail the driver on discrepancies.
3. Record separate construction (including observable readout preparation),
   single-evaluation, and numeric-propagation times. Serialize raw samples,
   versions, revisions, and machine metadata. Sweep depth conservatively to
   avoid uncontrolled symbolic expression growth.
4. Plot construction and evaluation scaling using the paper's QuEra colors.
   Report crossover counts only when symbolic evaluation beats numeric propagation.
5. Run relevant Rust tests and benchmark verification, review the diff, push the
   benchmark branch, create and attach a PR. Add a short blue algorithm paragraph
   and appendix benchmark discussion to the paper and rebuild with mise.
