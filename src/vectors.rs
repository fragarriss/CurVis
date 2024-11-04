/// This module contains the definition of the RelativisticVector and RelativisticObject structs,
/// representing respectively a 4-vector and a 4-vector pair (position, momentum) in the context of
/// General Relativity.

pub mod vectors {

    use core::fmt;
    use std;
    use nalgebra as na;
    
    /// Represents the covariance of a RelativisticVector.
    #[derive(Debug, Copy, Clone)]
    pub enum Covariance {
        Covariant,
        Contravariant,
    }

    impl fmt::Display for Covariance {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            match self {
                Covariance::Covariant => write!(f, "Covariant"),
                Covariance::Contravariant => write!(f, "Contravariant"),
            }
        }
    }

    impl PartialEq for Covariance {
        fn eq(&self, other: &Self) -> bool {
            match (self, other) {
                (Covariance::Covariant, Covariance::Covariant) => true,
                (Covariance::Contravariant, Covariance::Contravariant) => true,
                _ => false,
            }
        }
    }

    /// Represents a 4-dimensional vector in space-time coordinates, with its covariance.
    /// First coordinate is time, the other three are space.
    #[derive(Debug, Copy, Clone)]
    pub struct RelativisticVector {
        pub vector: na::Vector4<f64>,
        pub covariance: Covariance,
    }

    impl fmt::Display for RelativisticVector {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "{} ({}, {}, {}, {})", self.covariance, self.vector[0], self.vector[1], self.vector[2], self.vector[3])
        }
    }
    
    impl RelativisticVector {

        pub fn new(vector: na::Vector4<f64>, covariance: Covariance) -> RelativisticVector {
            RelativisticVector { vector, covariance }
        }

        /// Returns the i-th component of the vector.
        pub fn v(&self, i:usize) -> f64 {
            self.vector[i]
        }
    }

    impl std::ops::Add<f64> for RelativisticVector {
        type Output = RelativisticVector;

        fn add(self, other: f64) -> RelativisticVector {
            let new_vec = na::Vector4::new(self.vector[0] + other, self.vector[1] + other, self.vector[2] + other, self.vector[3] + other);
            RelativisticVector { vector: new_vec, covariance: self.covariance }
        }
    }
    
    impl std::ops::Sub<f64> for RelativisticVector {
        type Output = RelativisticVector;

        fn sub(self, other: f64) -> RelativisticVector {
            let new_vec = na::Vector4::new(self.vector[0] - other, self.vector[1] - other, self.vector[2] - other, self.vector[3] - other);
            RelativisticVector { vector: new_vec, covariance: self.covariance }
        }
    }

    impl std::ops::Mul<f64> for RelativisticVector {
        type Output = RelativisticVector;

        fn mul(self, other: f64) -> RelativisticVector {
            let new_vec = na::Vector4::new(self.vector[0] * other, self.vector[1] * other, self.vector[2] * other, self.vector[3] * other);
            RelativisticVector { vector: new_vec, covariance: self.covariance }
        }
    }

    impl std::ops::Div<f64> for RelativisticVector {
        type Output = RelativisticVector;

        fn div(self, other: f64) -> RelativisticVector {

            if other == 0.0 {
                panic!("Division by zero");
            }

            let new_vec = na::Vector4::new(self.vector[0] / other, self.vector[1] / other, self.vector[2] / other, self.vector[3] / other);
            RelativisticVector { vector: new_vec, covariance: self.covariance }
        }
    }

    // Implementing standard operations for RelativisticVector and RelativisticVector

    impl std::ops::Add<RelativisticVector> for RelativisticVector {
        type Output = RelativisticVector;

        fn add(self, other: RelativisticVector) -> RelativisticVector {
            if self.covariance != other.covariance {
                panic!("Cannot add vectors with different covariance");
            }
            let new_vec = self.vector + other.vector;
            RelativisticVector { vector: new_vec, covariance: self.covariance }
        }
    }

    impl std::ops::Sub<RelativisticVector> for RelativisticVector {
        type Output = RelativisticVector;

        fn sub(self, other: RelativisticVector) -> RelativisticVector {
            if self.covariance != other.covariance {
                panic!("Cannot subtract vectors with different covariance");
            }
            let new_vec = self.vector - other.vector;
            RelativisticVector { vector: new_vec, covariance: self.covariance }
        }
    }

    
    /// A pair of 4-dimensional vectors in space-time coordinates, representing the
    /// position and momentum of a relativistic object.
    /// For both position and moentum, the first coordinate is time, the other three are space.
    #[derive(Debug, Copy, Clone)]
    pub struct RelativisticObject
    {
        pub position: RelativisticVector,
        pub momentum: RelativisticVector,
    }

    impl RelativisticObject {

        pub fn new(position: RelativisticVector, momentum: RelativisticVector) -> RelativisticObject {
            RelativisticObject { position, momentum }
        }

        // Convenience functions to access the components of 
        // the position and momentum vectors

        /// Returns the i-th component of the position vector.
        pub fn x(&self, i:usize) -> f64 {
            self.position.v(i)
        }

        /// Returns the i-th component of the momentum vector.
        pub fn p(&self, i:usize) -> f64 {
            self.momentum.v(i)
        }

        // Convenience functions to access the covariance of
        // the position and momentum vectors

        /// Returns the covariance of the position vector.
        pub fn covariance_x(&self) -> &Covariance {
            &self.position.covariance
        }

        /// Returns the covariance of the momentum vector.
        pub fn covariance_p(&self) -> &Covariance {
            &self.momentum.covariance
        }

    }
}


// Tests

#[cfg(test)]
mod tests_vectors {

    use super::vectors::*;
    use nalgebra as na;

    #[test]
    fn test_relativistic_vector_addition() {
        let v1 = RelativisticVector::new(na::Vector4::new(1.0, 2.0, 3.0, 4.0), Covariance::Covariant);
        let v2 = RelativisticVector::new(na::Vector4::new(5.0, 6.0, 7.0, 8.0), Covariance::Covariant);
        let v3 = v1 + v2;
        assert_eq!(v3.vector, na::Vector4::new(6.0, 8.0, 10.0, 12.0));
    }

    #[test]
    fn test_relativistic_vector_subtraction() {
        let v1 = RelativisticVector::new(na::Vector4::new(1.0, 2.0, 3.0, 4.0), Covariance::Covariant);
        let v2 = RelativisticVector::new(na::Vector4::new(5.0, 6.0, 7.0, 8.0), Covariance::Covariant);
        let v3 = v1 - v2;
        assert_eq!(v3.vector, na::Vector4::new(-4.0, -4.0, -4.0, -4.0));
    }

    #[test]
    fn test_relativistic_vector_scalar_addition() {
        let v1 = RelativisticVector::new(na::Vector4::new(1.0, 2.0, 3.0, 4.0), Covariance::Covariant);
        let v2 = v1 + 5.0;
        assert_eq!(v2.vector, na::Vector4::new(6.0, 7.0, 8.0, 9.0));
    }

    #[test]
    fn test_relativistic_vector_scalar_subtraction() {
        let v1 = RelativisticVector::new(na::Vector4::new(1.0, 2.0, 3.0, 4.0), Covariance::Covariant);
        let v2 = v1 - 5.0;

        assert_eq!(v2.vector, na::Vector4::new(-4.0, -3.0, -2.0, -1.0));
    }

    #[test]
    #[should_panic]
    fn test_relativistic_vector_scalar_division_by_zero() {
        let v1 = RelativisticVector::new(na::Vector4::new(1.0, 2.0, 3.0, 4.0), Covariance::Covariant);
        let _v2 = v1 / 0.0;
    }

    #[test]
    #[should_panic]
    fn test_relativistic_vector_addition_different_covariance() {
        let v1 = RelativisticVector::new(na::Vector4::new(1.0, 2.0, 3.0, 4.0), Covariance::Covariant);
        let v2 = RelativisticVector::new(na::Vector4::new(5.0, 6.0, 7.0, 8.0), Covariance::Contravariant);
        let _v3 = v1 + v2;
    }

    #[test]
    #[should_panic]
    fn test_relativistic_vector_subtraction_different_covariance() {
        let v1 = RelativisticVector::new(na::Vector4::new(1.0, 2.0, 3.0, 4.0), Covariance::Covariant);
        let v2 = RelativisticVector::new(na::Vector4::new(5.0, 6.0, 7.0, 8.0), Covariance::Contravariant);
        let _v3 = v1 - v2;
    }

}