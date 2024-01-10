use std::{f64::consts::PI, path::Path};

use raytracer::{Canvas, Color, Matrix, Tuple};

// Chapter 4: Putting It Together
fn main() {
    let width = 300;
    let height = 300;
    // This panics if the canvas is too large for the point in the middle to be
    // represented by i32. If this is a problem, decrease canvas size or change to
    // i64.
    let (middle_x, middle_y): (i32, i32) = (
        (width / 2).try_into().unwrap(),
        (height / 2).try_into().unwrap(),
    );
    let radian_step = 2.0 * PI / 12.0;

    let mut canvas = Canvas::new(width, height);
    for i in 1..=12 {
        let position = Matrix::rotation_z(radian_step * (i as f64)).scale(120.0, 120.0, 0.0)
            * Tuple::point(1.0, 0.0, 0.0);

        // This panics if the canvas is too small for the point to be drawn. If this is a problem,
        // increase the canvas size, or move the point closer to the origin.
        let canvas_x = (middle_x + (position.x.round() as i32)).try_into().unwrap();
        let canvas_y = (middle_y - (position.y.round() as i32)).try_into().unwrap();
        canvas[[canvas_x, canvas_y]] = Color::new(1.0, 0.0, 0.0);
    }

    canvas
        .save(Path::new("./images/clock.ppm"))
        .expect("Unable to save canvas.")
}
