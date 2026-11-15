//! Internal record codec derivation for the versioned body handover format.

use proc_macro::TokenStream;
use quote::{format_ident, quote};
use syn::{Data, DeriveInput, Fields, parse_macro_input};

/// Derives the internal codec. Skipped fields are transient and reset to default.
#[proc_macro_derive(Transfer, attributes(transfer))]
pub fn derive(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let name = &input.ident;
    let mut types = Vec::new();
    let mut fields = |fields: &Fields| {
        let mut bindings = Vec::new();
        let mut writes = Vec::new();
        let mut reads = Vec::new();
        for (i, field) in fields.iter().enumerate() {
            let binding = format_ident!("field_{i}");
            let skipped = field.attrs.iter().any(|a| a.path().is_ident("transfer"));
            let ty = &field.ty;
            if skipped {
                bindings.push(quote!(_));
                reads.push(quote!(Default::default()));
            } else {
                bindings.push(quote!(#binding));
                writes.push(quote!(Codec::write(#binding, writer)?;));
                reads.push(quote!(<#ty as Codec>::read(reader)?));
                types.push(ty.clone());
            }
        }
        let (pattern, value) = match fields {
            Fields::Named(named) => {
                let names: Vec<_> = named.named.iter().map(|f| &f.ident).collect();
                (
                    quote!({ #(#names: #bindings),* }),
                    quote!({ #(#names: #reads),* }),
                )
            }
            Fields::Unnamed(_) => (quote!((#(#bindings),*)), quote!((#(#reads),*))),
            Fields::Unit => (quote!(), quote!()),
        };
        (pattern, value, writes)
    };
    let (write, read) = match &input.data {
        Data::Struct(data) => {
            let (pattern, value, writes) = fields(&data.fields);
            (
                quote!(let Self #pattern = self; #(#writes)*),
                quote!(Self #value),
            )
        }
        Data::Enum(data) => {
            let mut writes = Vec::new();
            let mut reads = Vec::new();
            for (i, variant) in data.variants.iter().enumerate() {
                let tag = i as u32;
                let name = &variant.ident;
                let (pattern, value, body) = fields(&variant.fields);
                writes.push(
                    quote!(Self::#name #pattern => { Codec::write(&#tag, writer)?; #(#body)* }),
                );
                reads.push(quote!(#tag => Self::#name #value));
            }
            (
                quote!(match self { #(#writes),* }),
                quote!(match <u32 as Codec>::read(reader)? {
                    #(#reads,)* _ => return Err(invalid("unknown record variant")),
                }),
            )
        }
        Data::Union(_) => {
            return syn::Error::new_spanned(name, "unions cannot be transferred")
                .to_compile_error()
                .into();
        }
    };
    let mut generics = input.generics.clone();
    for param in generics.type_params_mut() {
        param
            .bounds
            .push(syn::parse_quote!(catallax_numbers::transfer::Codec));
    }
    let (impl_generics, ty_generics, where_clause) = generics.split_for_impl();
    quote! {
        impl #impl_generics catallax_numbers::transfer::Codec for #name #ty_generics #where_clause {
            fn write(&self, writer: &mut catallax_numbers::transfer::Writer) -> Result<(), catallax_numbers::Refusal> {
                use catallax_numbers::transfer::Codec;
                #write
                Ok(())
            }
            fn read(reader: &mut catallax_numbers::transfer::Reader<'_>) -> Result<Self, catallax_numbers::Refusal> {
                use catallax_numbers::transfer::{Codec, invalid};
                let _ = invalid;
                Ok(#read)
            }
            fn register(registry: &mut catallax_numbers::transfer::Registry) {
                #(registry.add::<#types>();)*
            }
        }
    }.into()
}
