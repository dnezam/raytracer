use raytracer::{Matrix, Tuple};

// Chapter 3: Putting It Together
fn main() {
    println!("Inverted Identity");
    println!("{:?}", Matrix::<4>::IDENTITY.inverse());
    println!();

    println!("Multiplying a Matrix by its inverse");
    let matrix1 = Matrix::<4>::new([
        [3.0, -9.0, 7.0, 3.0],
        [3.0, -8.0, 2.0, -9.0],
        [-4.0, 4.0, 4.0, 1.0],
        [-6.0, 5.0, -1.0, 1.0],
    ]);
    println!("{:?}", matrix1 * matrix1.inverse().unwrap());
    println!();

    println!("Transpose of inverse");
    println!("{:?}", matrix1.inverse().unwrap().transpose());
    println!();
    println!("Inverse of Transpose");
    println!("{:?}", matrix1.transpose().inverse());
    println!();

    println!("Multiplying a tuple with almost the identity");
    let tuple = Tuple::new(1.0, 2.0, 3.0, 4.0);
    let mut almost_identity = Matrix::<4>::IDENTITY;
    almost_identity[[0, 1]] = 2.0;
    println!("{:?}", almost_identity * tuple);
}
