/// Trait for computing trace(self * RHS)
/// if type implements `Trace` an implementation of `TraceBy` is also provided.
pub trait Trace<'a, RHS: 'a, T> {
    fn trace(&'a self, value: &'a RHS) -> T;
}
