use std::path::Path;

use raytracer::{Canvas, Color, Matrix, Ray, Sphere, Tuple};

// Chapter 5: Putting It Together (following hints)
fn main() {
    // The key difference to circle.rs is that we will consider a wall that is not
    // set at z = 0.
    let canvas_pixels = 150;
    let mut canvas = Canvas::new(canvas_pixels, canvas_pixels);
    let red = Color::new(1.0, 0.0, 0.0);
    let mut sphere = Sphere::new();
    let transform = Matrix::<4>::IDENTITY.scale(1.0, 0.5, 1.0);
    sphere.set_transform(transform).unwrap();

    let wall_z = 10.0;
    let wall_size = 7.0;
    // Camera is at the looking at the center of the sphere, so half the wall
    // will be on its left/right. Additionally, the camera will be around the origin,
    // so half corresponds to the minimum and maximum of the x and y coordinates of the
    // wall.
    let half = wall_size / 2.0;
    // Moving one pixel on the canvas corresponds to moving pixel_size on the wall.
    let pixel_size = wall_size / canvas_pixels as f64;

    let camera_origin = Tuple::point(0.0, 0.0, -5.0);

    for (canvas_x, canvas_y, color) in canvas.iter_mut_with_pos() {
        // Top left corner of the wall is at (-half, half).
        // When canvas_x increases, we are moving to the right. Hence, we increment world_x.
        let world_x = -half + pixel_size * (canvas_x as f64);
        // When canvas_y increases, we are moving down. Hence, we decrement world_y.
        let world_y = half - pixel_size * (canvas_y as f64);

        // Wall coordinates corresponding to current canvas coordinates.
        let position = Tuple::point(world_x, world_y, wall_z);

        let direction = (position - camera_origin).normalize().unwrap();
        let camera_ray = Ray::new(camera_origin, direction).unwrap();

        if sphere.intersect(camera_ray).hit().is_some() {
            *color = red;
        }
    }

    canvas
        .save(Path::new("./images/circle_wall_shrink_y.ppm"))
        .unwrap();
}
