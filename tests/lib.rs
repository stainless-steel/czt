extern crate assert;
extern crate complex;
extern crate czt;

use complex::c64;

mod fixtures;

#[test]
fn forward() {
    let data = czt::forward(&fixtures::TIME_DATA);
    assert::close(as_f64(&data), &fixtures::FREQUENCY_DATA[..], 1e-14);
}

fn as_f64<'l>(slice: &'l [c64]) -> &'l [f64] {
    unsafe {
        std::slice::from_raw_parts(slice.as_ptr() as *const _, 2 * slice.len())
    }
}
