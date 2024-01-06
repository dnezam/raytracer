use std::{f64::consts::PI, path::Path};

use raytracer::{Canvas, Color, Matrix, Tuple};

fn main() {
    let width = 400;
    let height = 400;
    let (middle_x, middle_y) = (width / 2, height / 2);
    let radian_step = 2.0 * PI / 12.0;

    let mut canvas = Canvas::new(width, height);
    for i in 1..=12 {
        let position = Matrix::rotation_z(radian_step * (i as f64)).scale(150.0, 150.0, 0.0)
            * Tuple::point(1.0, 0.0, 0.0);
        canvas[[
            ((middle_x as i32) + (position.x.round() as i32)) as usize,
            ((middle_y as i32) - (position.y.round() as i32)) as usize,
        ]] = Color::new(1.0, 0.0, 0.0);
    }

    canvas
        .save(Path::new("./clock.ppm"))
        .expect("Unable to save canvas.")
}
