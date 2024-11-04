// This code is used to implement interpolation functionalities for the
// camera path.

pub mod interpolation {

    use std::path::PathBuf;

    use nalgebra as na;
    use crate::csv::csv::load_path;

    // --- Interpolation between two vectors ---

    // Define a custom trait for interpolation
    pub trait InterpolatableVector: Sized {
        fn interpolate(frac: f64, v1: Self, v2: Self) -> Self;
    }

    // Implement the trait for na::Vector3<f64>
    impl InterpolatableVector for na::Vector3<f64> {
        fn interpolate(frac: f64, v1: Self, v2: Self) -> Self {
            v1 + frac * (v2 - v1)
        }
    }

    // Implement the trait for na::Vector4<f64>
    impl InterpolatableVector for na::Vector4<f64> {
        fn interpolate(frac: f64, v1: Self, v2: Self) -> Self {
            v1 + frac * (v2 - v1)
        }
    }

    // Generic function to interpolate Vector3 or Vector4 objects
    pub fn interpolate_vector<T: InterpolatableVector>(frac: f64, v1: T, v2: T) -> T {
        // Check that frac is in [0,1]
        if !(frac >= 0.0 && frac <= 1.0) {
            panic!("frac must be between 0 and 1")
        }
    
        // Interpolate
        T::interpolate(frac, v1, v2)
    }

    // --- Interpolator struct ---

    pub struct Interpolator {
        pub positions: Vec<na::Vector4<f64>>,
        pub forward_vectors: Vec<na::Vector3<f64>>,
        pub up_vectors: Vec<na::Vector3<f64>>,
    }

    impl Interpolator {

        pub fn min_time(&self) -> f64 {
            return self.positions[0][0]; // first index, first coordinate (time)
        }

        pub fn max_time(&self) -> f64 {
            return self.positions[self.positions.len() - 1][0]; // last index, first coordinate (time)
        }

        /// Given the time variable t, it returns the indexes of t1 and t2
        /// between which t is located and the fractional part of t.
        fn time_indexes_and_frac_from_time(&self, t: f64) -> (usize, usize, f64) {
            // if t is not in the range of the camera path, panics
            if t < self.min_time() {
                panic!("Interpolation time cannot be smaller than first time in positions[0].")
            }

            if t > self.max_time() {
                panic!("Interpolation time cannot be greater than last time in positions[0].")
            }

            // Find the two times between which t is located
            let mut t1 = self.min_time();
            let mut t2 = self.max_time();
            let mut i = 0;

            while t > self.positions[i][0] {
                t1 = self.positions[i][0];
                t2 = self.positions[i + 1][0];
                i += 1;
            }

            // Compute the fractional part of t
            let frac = (t - t1) / (t2 - t1);

            let index1 = i;
            let index2 = i + 1;

            (index1, index2, frac)
        }
        
        pub fn camera_position(&self, t: f64) -> na::Vector4<f64> {
            let (index1, index2, frac) = self.time_indexes_and_frac_from_time(t);
            interpolate_vector(frac, self.positions[index1], self.positions[index2])
        }

        pub fn camera_up(&self, t: f64) -> na::Vector3<f64> {
            let (index1, index2, frac) = self.time_indexes_and_frac_from_time(t);
            interpolate_vector(frac, self.up_vectors[index1], self.up_vectors[index2])
        }

        pub fn camera_forward(&self, t: f64) -> na::Vector3<f64> {
            let (index1, index2, frac) = self.time_indexes_and_frac_from_time(t);
            interpolate_vector(frac, self.forward_vectors[index1], self.forward_vectors[index2])
        }

        pub fn from_file(path_to_csv_file: &PathBuf) -> Self {
            let (positions, forward_vectors, up_vectors) = load_path(path_to_csv_file);
            Interpolator { positions, forward_vectors, up_vectors }
        }
    }

}