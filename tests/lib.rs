extern crate assert;
extern crate czt;

use czt::{Transform, c64};

macro_rules! chirp(
    ($m:expr) => ({
        use std::f64::consts::PI;
        c64::from_polar(&1.0, &(-2.0 * PI / $m as f64))
    });
);

#[test]
fn forward_real() {
    let data = [
        1.576130816775483e-01,
        9.705927817606157e-01,
        9.571669482429456e-01,
        4.853756487228412e-01,
        8.002804688888001e-01,
        1.418863386272153e-01,
        4.217612826262750e-01,
        9.157355251890671e-01,
        7.922073295595544e-01,
        9.594924263929030e-01,
    ];
    let expected_result = [
         6.602111831687766e+00, -3.330669073875470e-16,
         6.961651498427756e-01,  2.339802401557220e-02,
        -1.275236051698197e+00, -4.839463190859927e-04,
        -4.848214873681843e-01, -5.265485614227395e-01,
        -1.277071313383780e+00,  7.821458553540981e-01,
        -3.440536096975195e-01, -3.521657193471262e-16,
        -1.277071313383781e+00, -7.821458553540985e-01,
        -4.848214873681850e-01,  5.265485614227390e-01,
        -1.275236051698199e+00,  4.839463190835502e-04,
         6.961651498427770e-01, -2.339802401557728e-02,
    ];

    let m = data.len();
    let result = data.transform(m, chirp!(m), c64::new(1.0, 0.0));
    assert::close(as_f64(&result), &expected_result[..], 1e-14);
}

#[test]
fn forward_complex() {
    let data = [
        8.147236863931789e-01, 1.576130816775483e-01,
        9.057919370756192e-01, 9.705927817606157e-01,
        1.269868162935061e-01, 9.571669482429456e-01,
        9.133758561390194e-01, 4.853756487228412e-01,
        6.323592462254095e-01, 8.002804688888001e-01,
        9.754040499940952e-02, 1.418863386272153e-01,
        2.784982188670484e-01, 4.217612826262750e-01,
        5.468815192049838e-01, 9.157355251890671e-01,
        9.575068354342976e-01, 7.922073295595544e-01,
        9.648885351992765e-01, 9.594924263929030e-01,
    ];
    let expected_result = [
         6.238553055831749e+00,  6.602111831687766e+00,
         1.354181005313841e+00,  9.642216142292481e-01,
        -2.864562964973690e-01, -1.789026257465937e-01,
         1.251129460280184e+00, -1.037906571249956e+00,
        -1.333821911972047e+00, -2.172769214747543e+00,
        -6.184034494048681e-01, -3.440536096975196e-01,
         2.304697987361479e-01, -3.813734120200171e-01,
         1.980323374347046e-01,  6.826359651358684e-02,
        -2.874241891355410e-01, -2.371569477649806e+00,
         1.400977053344992e+00,  4.281086854562995e-01,
    ];

    let m = data.len() / 2;
    let result = as_c64(&data).transform(m, chirp!(m), c64::new(1.0, 0.0));
    assert::close(as_f64(&result), &expected_result[..], 1e-14);
}

#[test]
fn forward_complex_small_m() {
    let data = [
        8.147236863931789e-01, 1.576130816775483e-01,
        9.057919370756192e-01, 9.705927817606157e-01,
        1.269868162935061e-01, 9.571669482429456e-01,
        9.133758561390194e-01, 4.853756487228412e-01,
        6.323592462254095e-01, 8.002804688888001e-01,
        9.754040499940952e-02, 1.418863386272153e-01,
        2.784982188670484e-01, 4.217612826262750e-01,
        5.468815192049838e-01, 9.157355251890671e-01,
        9.575068354342976e-01, 7.922073295595544e-01,
        9.648885351992765e-01, 9.594924263929030e-01,
    ];
    let expected_result = [
         6.238553055831749e+00,  6.602111831687765e+00,
        -2.864562964973671e-01, -1.789026257465932e-01,
        -1.333821911972045e+00, -2.172769214747541e+00,
         2.304697987361482e-01, -3.813734120200151e-01,
        -2.874241891355401e-01, -2.371569477649796e+00,
    ];

    let m = data.len() / 2 - 5;
    let result = as_c64(&data).transform(m, chirp!(m), c64::new(1.0, 0.0));
    assert::close(as_f64(&result), &expected_result[..], 1e-14);
}

#[test]
fn forward_complex_large_m() {
    let data = [
        8.147236863931789e-01, 1.576130816775483e-01,
        9.057919370756192e-01, 9.705927817606157e-01,
        1.269868162935061e-01, 9.571669482429456e-01,
        9.133758561390194e-01, 4.853756487228412e-01,
        6.323592462254095e-01, 8.002804688888001e-01,
        9.754040499940952e-02, 1.418863386272153e-01,
        2.784982188670484e-01, 4.217612826262750e-01,
        5.468815192049838e-01, 9.157355251890671e-01,
        9.575068354342976e-01, 7.922073295595544e-01,
        9.648885351992765e-01, 9.594924263929030e-01,
    ];
    let expected_result = [
         6.238553055831749e+00,  6.602111831687766e+00,
         1.613215736819029e+00, -2.663924556100892e+00,
         3.890265936019357e+00, -4.924428190758741e-01,
        -2.864562964973641e-01, -1.789026257465957e-01,
         2.385011013625547e+00, -1.100610328758628e+00,
         2.026744627883097e+00, -1.046712023541919e+00,
        -1.333821911972041e+00, -2.172769214747538e+00,
         1.673305247584364e-01, -1.296818180928050e-01,
        -3.530623771613099e-01,  9.353211460426021e-01,
         2.304697987361495e-01, -3.813734120200123e-01,
         6.491612060807203e-01,  5.173275101128529e-01,
         4.304078885610212e-01,  1.011461952140921e+00,
        -2.874241891355410e-01, -2.371569477649782e+00,
        -3.705651933107217e-01,  3.387876411534603e+00,
        -2.778974524340422e+00,  4.480836493785334e-01,
    ];

    let m = data.len() / 2 + 5;
    let result = as_c64(&data).transform(m, chirp!(m), c64::new(1.0, 0.0));
    assert::close(as_f64(&result), &expected_result[..], 1e-13);
}

#[test]
fn forward_complex_different_w() {
    let data = [
        8.147236863931789e-01, 1.576130816775483e-01,
        9.057919370756192e-01, 9.705927817606157e-01,
        1.269868162935061e-01, 9.571669482429456e-01,
        9.133758561390194e-01, 4.853756487228412e-01,
        6.323592462254095e-01, 8.002804688888001e-01,
        9.754040499940952e-02, 1.418863386272153e-01,
        2.784982188670484e-01, 4.217612826262750e-01,
        5.468815192049838e-01, 9.157355251890671e-01,
        9.575068354342976e-01, 7.922073295595544e-01,
        9.648885351992765e-01, 9.594924263929030e-01,
    ];
    let expected_result = [
         6.238553055831749e+00,  6.602111831687766e+00,
         5.704267697623152e-01,  1.008425689197299e-01,
         1.066071341216861e+00, -2.969241251163725e+00,
         3.793137728761690e+00, -3.281105764362009e+00,
         6.616178717582136e-01,  1.063766981694492e+00,
        -1.507929713128197e+00, -4.342486874795381e-01,
         2.056896326647222e+00,  1.188905648029616e+00,
         6.799674664045305e-01, -2.144181560077285e+00,
         8.633970898736677e-02, -3.928896271523249e-01,
         2.144248604875231e+00, -2.670035172794396e+00,
    ];

    let m = data.len() / 2;
    let result = as_c64(&data).transform(m, c64::from_polar(&1.0, &-42.0), c64::new(1.0, 0.0));
    assert::close(as_f64(&result), &expected_result[..], 1e-13);
}

#[test]
fn forward_complex_different_a() {
    let data = [
        8.147236863931789e-01, 1.576130816775483e-01,
        9.057919370756192e-01, 9.705927817606157e-01,
        1.269868162935061e-01, 9.571669482429456e-01,
        9.133758561390194e-01, 4.853756487228412e-01,
        6.323592462254095e-01, 8.002804688888001e-01,
        9.754040499940952e-02, 1.418863386272153e-01,
        2.784982188670484e-01, 4.217612826262750e-01,
        5.468815192049838e-01, 9.157355251890671e-01,
        9.575068354342976e-01, 7.922073295595544e-01,
        9.648885351992765e-01, 9.594924263929030e-01,
    ];
    let expected_result = [
         8.322877644021808e+00,  2.336374861036936e+00,
         3.159332838220827e+00,  8.743786054620886e-01,
         8.214456639585297e-01,  9.585207566877277e-01,
         1.557156328730404e+00, -7.222742107919995e-01,
        -1.608963320209646e+00, -7.310985452013501e-01,
        -7.785003810981718e-01,  3.885863506629339e-01,
        -2.509921566194973e-01,  1.536868581117541e-02,
         9.883612335520264e-03,  5.685381440874357e-01,
        -2.415596643839245e+00, -1.461132642318109e+00,
        -6.694067215687522e-01, -6.511311886613711e-01,
    ];

    let m = data.len() / 2;
    let result = as_c64(&data).transform(m, chirp!(m), c64::from_polar(&1.0, &-69.0));
    assert::close(as_f64(&result), &expected_result[..], 1e-13);
}

fn as_f64<'l>(slice: &'l [c64]) -> &'l [f64] {
    unsafe {
        std::slice::from_raw_parts(slice.as_ptr() as *const _, 2 * slice.len())
    }
}

fn as_c64<'l>(slice: &'l [f64]) -> &'l [c64] {
    unsafe {
        std::slice::from_raw_parts(slice.as_ptr() as *const _, slice.len() / 2)
    }
}
