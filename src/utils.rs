use std::sync::atomic::{AtomicUsize, Ordering};

/// Determines how close values must be such that they are approximately equal.
const EPSILON: f64 = 0.00001;

/// Implements approximate equality for f64.
pub fn eq(a: f64, b: f64) -> bool {
    (a - b).abs() < EPSILON
}

/// Returns a unique id.
pub fn get_id() -> usize {
    static COUNTER: AtomicUsize = AtomicUsize::new(1);
    COUNTER.fetch_add(1, Ordering::Relaxed)
}
