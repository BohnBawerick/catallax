//! Standalone wasm timing entry points, compiled with rustc as documented in docs.

use catallax_app::handover::{read_body, write_body};
use catallax_features::Built;
use catallax_numbers::{Decider, with_operation_budget};
use std::cell::RefCell;

#[derive(Default)]
struct State {
    source: Option<Built>,
    bytes: Vec<u8>,
    adopted: Option<Built>,
}

thread_local! {
    static STATE: RefCell<State> = RefCell::new(State::default());
}

#[unsafe(no_mangle)]
pub extern "C" fn input_buffer(length: usize) -> *mut u8 {
    STATE.with_borrow_mut(|state| {
        state.bytes.resize(length, 0);
        state.bytes.as_mut_ptr()
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn source_from_input() {
    STATE.with_borrow_mut(|state| {
        state.source = Some(read_body(&state.bytes, &mut Decider::default()).unwrap());
    });
}

#[unsafe(no_mangle)]
pub extern "C" fn clear() {
    STATE.with_borrow_mut(|state| state.adopted = None);
}

#[unsafe(no_mangle)]
pub extern "C" fn probe_read() -> f64 {
    STATE.with_borrow_mut(|state| {
        let (body, work) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
            read_body(&state.bytes, &mut Decider::default())
        });
        state.adopted = Some(body.unwrap());
        work.used as f64
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn probe_write() -> f64 {
    STATE.with_borrow_mut(|state| {
        let (bytes, work) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
            write_body(state.source.as_ref().unwrap(), &mut Decider::default())
        });
        state.bytes = bytes.unwrap();
        work.used as f64
    })
}
