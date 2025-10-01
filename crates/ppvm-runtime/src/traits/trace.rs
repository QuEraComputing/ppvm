/// Trait for computing trace(self * RHS)
/// if type implements `Trace` an implementation of `TraceBy` is also provided.
pub trait Trace<'a, RHS: 'a> {
    type Output;
    fn trace(&'a self, value: &'a RHS) -> Self::Output;
}
