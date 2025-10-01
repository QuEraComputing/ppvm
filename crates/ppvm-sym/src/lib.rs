mod add;
mod display;
mod mul;
mod term;

pub use term::{Item, Prod, Sum};

pub fn sin(u: u32) -> Item {
    Item::Sin(u)
}

pub fn cos(u: u32) -> Item {
    Item::Cos(u)
}
