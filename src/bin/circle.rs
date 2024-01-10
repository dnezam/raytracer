use std::path::Path;

use raytracer::{Canvas, Color, Matrix, Ray, Sphere, Tuple};

// Chapter 5: Putting It Together
fn main() {
    // Set up canvas
    let canvas_size = 100;
    let mut canvas = Canvas::new(canvas_size, canvas_size);

    // Set up scene
    let mut sphere = Sphere::new();
    let sphere_size = (canvas_size as f64) * 0.4;
    sphere
        .set_transform(Matrix::<4>::scaling(sphere_size, sphere_size, sphere_size))
        .unwrap();
    let camera_origin = Tuple::point(0.0, 0.0, -4.0 * sphere_size);

    for (canvas_x, canvas_y, color) in canvas.iter_mut_with_pos() {
        // Canvas coordinates -> Scene coordinates
        let scene_x = (canvas_x as f64) - (canvas_size as f64) / 2.0;
        let scene_y = -(canvas_y as f64) + (canvas_size as f64) / 2.0;
        let direction = (Tuple::point(scene_x, scene_y, 0.0) - camera_origin)
            .normalize()
            .unwrap();
        let ray = Ray::new(camera_origin, direction).unwrap();

        if sphere.intersect(ray).hit().is_some() {
            *color = Color::new(1.0, 0.0, 0.0);
        }
    }

    canvas
        .save(Path::new("./circle.ppm"))
        .expect("Unable to save canvas.")
}
