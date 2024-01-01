const EPSILON: f64 = 0.00001;

pub fn eq(a: f64, b: f64) -> bool {
    (a - b).abs() < EPSILON
}
