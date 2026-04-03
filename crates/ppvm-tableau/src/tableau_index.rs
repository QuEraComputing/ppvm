use num::PrimInt;
use std::hash::Hash;
use std::ops::{BitAnd, BitOrAssign, BitXor, Shl};

pub trait TableauIndex:
    PartialEq
    + Eq
    + Hash
    + Copy
    + From<u8>
    + Shl<usize, Output = Self>
    + BitOrAssign<Self>
    + BitAnd<Self, Output = Self>
    + BitXor<Output = Self>
    + PrimInt
{
}

impl<I> TableauIndex for I where
    I: PartialEq
        + Eq
        + Hash
        + Copy
        + From<u8>
        + Shl<usize, Output = Self>
        + BitOrAssign<Self>
        + BitAnd<I, Output = I>
        + BitXor<Output = I>
        + PrimInt
{
}
