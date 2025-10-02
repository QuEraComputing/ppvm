use proc_macro::TokenStream;
use quote::quote;
use syn::{DeriveInput, parse_macro_input};

#[proc_macro_derive(ACMap)]
pub fn derive_acmap(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let name = input.ident;

    match input.data {
        syn::Data::Struct(_) => {}
        _ => {
            return syn::Error::new_spanned(name, "ACMap can only be derived for structs")
                .to_compile_error()
                .into();
        }
    }

    let expanded = quote! {
        impl<'a, S, V, State> ppvm_runtime::traits::ACMap<S, V> for #name<ppvm_runtime::word::PauliWord<S>, V, State>
        where
            S: ppvm_runtime::traits::PauliStorage,
            V: ppvm_runtime::traits::Coefficient,
            State: Clone + std::hash::BuildHasher + Default,
        {
            fn with_capacity(capacity: usize) -> Self {
                Self::with_capacity_and_hasher(capacity, State::default())
            }

            fn len(&self) -> usize {
                self.len()
            }
        }
    };

    TokenStream::from(expanded)
}
