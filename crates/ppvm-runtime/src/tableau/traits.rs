use crate::config::Config;
use crate::traits::Clifford;
use std::hash::Hash;
use std::ops::{BitAnd, BitOrAssign, BitXor, Shl};

pub trait TGate<T: Config> {
    fn t(&mut self, addr0: usize);
    fn t_adj(&mut self, addr0: usize);
}

pub trait Measure {
    fn measure(&mut self, addr0: usize) -> bool;
}

pub trait CliffordExtensions: Clifford {
    fn s_adj(&mut self, addr0: usize);
    fn sqrt_x(&mut self, addr0: usize) {
        self.h(addr0);
        self.s(addr0);
        self.h(addr0);
    }

    fn sqrt_x_adj(&mut self, addr0: usize) {
        self.h(addr0);
        self.s_adj(addr0);
        self.h(addr0);
    }

    fn sqrt_y(&mut self, addr0: usize) {
        // NOTE: in matmul (RHS applied first)
        // SqrtY == S * SqrtX * S'
        self.s_adj(addr0);
        self.sqrt_x(addr0);
        self.s(addr0);
    }

    fn sqrt_y_adj(&mut self, addr0: usize) {
        self.s_adj(addr0);
        self.sqrt_x_adj(addr0);
        self.s(addr0);
    }
}

pub trait TableauIndex:
    PartialEq
    + Eq
    + Hash
    + Copy
    + From<u8>
    + Shl<usize>
    + BitOrAssign<<Self as Shl<usize>>::Output>
    + BitAnd<<Self as Shl<usize>>::Output, Output = Self>
    + BitXor<Output = Self>
{
}

impl<I> TableauIndex for I where
    I: PartialEq
        + Eq
        + Hash
        + Copy
        + From<u8>
        + Shl<usize>
        + BitOrAssign<<I as Shl<usize>>::Output>
        + BitAnd<<I as Shl<usize>>::Output, Output = I>
        + BitXor<Output = I>
{
}
