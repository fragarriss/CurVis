/// This module implement metrics of wormhole space-times.
pub mod metrics {

    use nalgebra as na;
    use crate::vectors::vectors::{RelativisticVector, Covariance, RelativisticObject};
    use std::f64::consts::PI;

    /// Panics if the passed argument is not Covariance::Contravariant
    fn check_contravariance(vec: &RelativisticVector) {
        if vec.covariance != Covariance::Contravariant {
            panic!("The vector must be contravariant.");
        }
    }

    /* not used */
    // /// Panics if the passed argument is not Covariance::Covariant
    // fn check_covariance(vec: &RelativisticVector) {
    //     if vec.covariance != Covariance::Covariant {
    //         panic!("The vector must be covariant.");
    //     }
    // }

    /// Trait that defines the methods that a metric must implement to be used in the project.
    /// In particular, the crate currently support spherically symmetric space-times with
    /// metrics that are diagonal (i.e. the metric tensor has only non-zero diagonal components).
    /// 
    /// The metric is defined in terms of the time coordinate *t*, the radial coordinate *l* and
    /// the two polar coordinates *theta* and *phi*.
    /// 
    /// To implement this trait, a struct must define the methods `r(l)`, `r_squared(l)` and `r_derivative(l)`.
    /// See the general documentation of this crate for more information on the meaning `r(l)` for
    /// wormhole metrics.
    /// 
    /// The methods `object_position_diff_contr()` and `object_momentum_diff_cov()` provided by this
    /// trait are those that implement the differential equations for the evolution of a relativistic
    /// object in the metric, as described by O. James et al. in their paper (https://doi.org/10.1119/1.4916949).
    /// 
    /// The `update_relativistic_object()` method uses the previous two methods to update the position
    /// and momentum of a relativistic object over a small proper time intervall.
    pub trait DiagonalSphericalMetric {

        fn r(&self, l:f64) -> f64;
        fn r_squared(&self, l:f64) -> f64;
        fn r_derivative(&self, l:f64) -> f64;

        // Covariant metric components

        /// Returns the covariant (0,0) (time) component of the metric tensor at the given position.
        fn g00(&self, position_contr: &RelativisticVector) -> f64 { 
            check_contravariance(position_contr);
            -1.0 
        }

        /// Returns the covariant (1,1) (radial) component of the metric tensor at the given position.
        fn g11(&self, position_contr: &RelativisticVector) -> f64 {
            check_contravariance(position_contr);
            1.0
        }

        /// Returns the covariant (2,2) (theta) component of the metric tensor at the given position.
        fn g22(&self, position_contr: &RelativisticVector) -> f64 {
            check_contravariance(position_contr);
            self.r_squared(position_contr.v(1)) } // Optimized self.r(l).powi(2)

        /// Returns the covariant (3,3) (phi) component of the metric tensor at the given position.
        fn g33(&self, position_contr: &RelativisticVector) -> f64 {
            check_contravariance(position_contr);
            self.r_squared(position_contr.v(1)) * position_contr.v(2).sin().powi(2) } // Optimized self.r(l).powi(2) * theta.powi(2).sin()

        // Returns the (i,i) component of the covariant metric tensor at the given position.
        fn gii(&self, i:usize, position:&RelativisticVector) -> f64 {
            match i {
                0 => self.g00(position),
                1 => self.g11(position),
                2 => self.g22(position),
                3 => self.g33(position),
                _ => panic!("Invalid index for the covariant metric components."),
            }
        }

        // Contravariant metric components

        /// Returns the contravariant (0,0) (time) component of the metric tensor at the given position.
        fn g00_contr(&self, position_contr:&RelativisticVector) -> f64 { self.g00(position_contr).powi(-1) }

        /// Returns the contravariant (1,1) (radial) component of the metric tensor at the given position.
        fn g11_contr(&self, position_contr:&RelativisticVector) -> f64 { self.g11(position_contr).powi(-1) }

        /// Returns the contravariant (2,2) (theta) component of the metric tensor at the given position.
        fn g22_contr(&self, position_contr:&RelativisticVector) -> f64 { self.g22(position_contr).powi(-1) }

        /// Returns the contravariant (3,3) (phi) component of the metric tensor at the given position.
        fn g33_contr(&self, position_contr:&RelativisticVector) -> f64 { self.g33(position_contr).powi(-1) }

        // Returns the (i,i) component of the contravariant metric tensor at the given position.
        fn gii_contr(&self, i:usize, position_contr:&RelativisticVector) -> f64 {
            match i {
                0 => self.g00_contr(position_contr),
                1 => self.g11_contr(position_contr),
                2 => self.g22_contr(position_contr),
                3 => self.g33_contr(position_contr),
                _ => panic!("Invalid index for the contravariant metric components."),
            }
        }

        // Frame field components

        /// Returns the (0,0) (time) component of the frame field at the given position.
        fn frame_field_00(&self, position_contr: &RelativisticVector) -> f64 { 
            check_contravariance(position_contr);
            1.0 
        }

        /// Returns the (1,1) (radial) component of the frame field at the given position.
        fn frame_field_11(&self, position_contr: &RelativisticVector) -> f64 { 
            check_contravariance(position_contr);
            1.0 
        }

        /// Returns the (2,2) (theta) component of the frame field at the given position.
        fn frame_field_22(&self, position_contr: &RelativisticVector) -> f64 { 
            check_contravariance(position_contr);
            self.r(position_contr.v(1))
        }

        /// Returns the (3,3) (phi) component of the frame field at the given position.
        fn frame_field_33(&self, position_contr: &RelativisticVector) -> f64 { 
            check_contravariance(position_contr);
            self.r(position_contr.v(1)) * position_contr.v(2).sin()
        }

        // Inverse frame field components

        /// Returns the (0,0) (time) component of the inverse frame field at the given position.
        fn inverse_frame_field_00(&self, position_contr: &RelativisticVector) -> f64 { 
            check_contravariance(position_contr);
            1.0 
        }

        /// Returns the (1,1) (radial) component of the inverse frame field at the given position.
        fn inverse_frame_field_11(&self, position_contr: &RelativisticVector) -> f64 { 
            check_contravariance(position_contr);
            1.0 
        }

        /// Returns the (2,2) (theta) component of the inverse frame field at the given position.
        fn inverse_frame_field_22(&self, position_contr: &RelativisticVector) -> f64 { 
            check_contravariance(position_contr);
            1.0 / self.r(position_contr.v(1))
        }

        /// Returns the (3,3) (phi) component of the inverse frame field at the given position.
        fn inverse_frame_field_33(&self, position_contr: &RelativisticVector) -> f64 { 
            check_contravariance(position_contr);
            1.0 / (self.r(position_contr.v(1)) * position_contr.v(2).sin())
        }

        /// Returns the covariant `RelativisticVector` to its contravariant counterpart assuming 
        /// it is calculated at the given position.
        fn to_covariant(&self, position_contr: &RelativisticVector, vector_contr: &RelativisticVector) -> RelativisticVector {
            check_contravariance(&vector_contr);
            RelativisticVector::new(
                na::Vector4::new(
                    vector_contr.v(0)*self.g00(position_contr),
                    vector_contr.v(1)*self.g11(position_contr),
                    vector_contr.v(2)*self.g22(position_contr),
                    vector_contr.v(3)*self.g33(position_contr),
                ),
                Covariance::Covariant,
            )
        }

        /// Mutates the covariant `RelativisticVector` to its contravariant counterpart assuming it is calculated
        /// at the given position.
        fn to_covariant_mut(&self, position_contr: &RelativisticVector, vector_contr: &mut RelativisticVector) {
            if vector_contr.covariance == Covariance::Covariant {
                panic!("The vector is already covariant.");
            }
            vector_contr.vector = na::Vector4::new(
                vector_contr.v(0)*self.g00(position_contr),
                vector_contr.v(1)*self.g11(position_contr),
                vector_contr.v(2)*self.g22(position_contr),
                vector_contr.v(3)*self.g33(position_contr),
            );
            vector_contr.covariance = Covariance::Covariant;
        }

        /// Returns the contravariant `RelativisticVector` to its covariant counterpart assuming it is calculated
        /// at the given position.
        fn to_contravariant(&self, position_contr: &RelativisticVector, vector_cov: &RelativisticVector) -> RelativisticVector {
            if vector_cov.covariance == Covariance::Contravariant {
                panic!("The vector is already contravariant.");
            }
            RelativisticVector::new(
                na::Vector4::new(
                    vector_cov.v(0)*self.g00_contr(position_contr),
                    vector_cov.v(1)*self.g11_contr(position_contr),
                    vector_cov.v(2)*self.g22_contr(position_contr),
                    vector_cov.v(3)*self.g33_contr(position_contr),
                ),
                Covariance::Contravariant,
            )
        }

        /// Mutates the contravariant `RelativisticVector` to its covariant counterpart assuming it is calculated
        /// at the given position.
        fn to_contravariant_mut(&self, position_contr: &RelativisticVector, vector_cov: &mut RelativisticVector) {
            if vector_cov.covariance == Covariance::Contravariant {
                panic!("The vector is already contravariant.");
            }
            vector_cov.vector = na::Vector4::new(
                vector_cov.v(0)*self.g00_contr(position_contr),
                vector_cov.v(1)*self.g11_contr(position_contr),
                vector_cov.v(2)*self.g22_contr(position_contr),
                vector_cov.v(3)*self.g33_contr(position_contr),
            );
            vector_cov.covariance = Covariance::Contravariant;
        }

        /// Returns the differential of the object's position with contravariant components.
        /// This is given by multiplying the object's covariant momentum by the covariant
        /// metric components.
        fn object_position_diff_contr(&self, object: &RelativisticObject) -> RelativisticVector {    
            if object.position.covariance != Covariance::Contravariant {
                panic!("The position vector must be contravariant.");
            }

            let momentum_cov = if object.momentum.covariance == Covariance::Contravariant {
                self.to_covariant(&object.position, &object.momentum)
            }
            else {
                object.momentum
            };

            RelativisticVector::new(
                na::Vector4::new(
                    momentum_cov.v(0) * self.g00_contr(&object.position),
                    momentum_cov.v(1) * self.g11_contr(&object.position),
                    momentum_cov.v(2) * self.g22_contr(&object.position),
                    momentum_cov.v(3) * self.g33_contr(&object.position),
                ),
                Covariance::Contravariant,
            )
        }

        /// Returns the differential of the object's momentum with covariant components
        fn object_momentum_diff_cov(&self, object: &RelativisticObject) -> RelativisticVector {

            // Normalizing object such that momentum is covariant
            let momentum_cov = if object.momentum.covariance == Covariance::Contravariant {
                self.to_covariant(&object.position, &object.momentum)
            }
            else {
                object.momentum
            };

            let b_squared = momentum_cov.v(2).powi(2) + momentum_cov.v(3).powi(2)/(object.x(2).sin().powi(2));

            let dp = na::Vector4::new(
                0.0,
                b_squared * self.r_derivative(object.x(1)) / (self.r(object.x(1)).powi(3)),
                momentum_cov.v(3).powi(2) * (object.x(2).cos() / (self.r_squared(object.x(1)) * object.x(2).sin().powi(3))),
                0.0,
            );

            RelativisticVector::new(
                dp,
                Covariance::Covariant,
            )
        }
   
        /// This function updates the position and momentum of a relativistic object that
        /// is propagating over the metric for a small proper time differential `delta`.
        /// 
        /// Uses internally the `object_position_diff_contr()` and `object_momentum_diff_cov()`
        /// methods.
        /// 
        /// The `delta` factor can be negative, which corresponds to evolving the object
        /// back in time
        /// 
        /// [To be checked and tested if positive delta actually correspond to propagating forward
        /// in time.]
        fn update_relativistic_object(&self, object: &mut RelativisticObject, delta: f64) {

            // Normalizing object such that momentum is covariant
            if object.momentum.covariance == Covariance::Contravariant {
                self.to_covariant_mut(&object.position, &mut object.momentum);
            }

            let x_diff_contr: RelativisticVector = self.object_position_diff_contr(object);
        
            // Calculating momentum differential
            let p_diff_cov: RelativisticVector = self.object_momentum_diff_cov(object);

            object.position = object.position + x_diff_contr * delta;
            object.momentum = object.momentum + p_diff_cov * delta;
        }

        /// Creates a new `RelativisticObject` from the given position and direction such
        /// that the object is a photon (light-like object).
        fn new_photon (&self, position: &RelativisticVector, direction:na::Vector3<f64>) -> RelativisticObject {
            
            /* Code removed because I was wrongly convinced that a light-like object has
            vanishing momentum norm. It's not true: it's *world-line* has length zero, not the norm. */

            // let space_norm = norm(
            //     &RelativisticVector::new(
            //         na::Vector4::new(0.0, direction[0], direction[1], direction[2]),
            //         Covariance::Contravariant,
            //     ),
            //     position,
            //     metric,
            // );

            // let p_but_p0 = RelativisticVector::new(
            //     na::Vector4::new(0.0, direction[0]/space_norm, direction[1]/space_norm, direction[2]/space_norm),
            //     Covariance::Contravariant,
            // );

            let direction = direction.normalize();

            RelativisticObject::new(
                RelativisticVector::new(position.vector, Covariance::Contravariant),
                RelativisticVector::new(
                    na::Vector4::new(
                        1.0,
                        direction[0],
                        direction[1] * self.r(position.v(1)),
                        direction[2] * self.r(position.v(1)) * position.v(2).sin(),
                    ),
                    Covariance::Covariant
                )
            )
        }

        /// Given a `RelativisticVector` at position `position_contr`, returns the direction
        /// of the vector in the tangent space at that position (the basis of the tangent
        /// space is given by the frame field).
        fn relativistic_vector_to_direction(&self, vector: &RelativisticVector, position_contr: &RelativisticVector) -> na::Vector3<f64> {
            let vector_contr: &RelativisticVector = if vector.covariance == Covariance::Covariant {
                &self.to_contravariant(position_contr, vector)
            } else { vector };

            na::Vector3::new(
                vector_contr.v(1) * self.frame_field_11(position_contr),
                vector_contr.v(2) * self.frame_field_22(position_contr),
                vector_contr.v(3) * self.frame_field_22(position_contr),
            )
        }
        
    }

    /// Computes the dot product between two vector `v1` and `v2` defined at the given `position` on the metric.
    /// It internally takes into account the covariance of `v1` and `v2`.
    pub fn dot_product(v1: &RelativisticVector, v2: &RelativisticVector, position: &RelativisticVector, metric: &impl DiagonalSphericalMetric) -> f64 {
        
        // Should expand this to include the case of covariant vectors and contravariant metrics
        let v1 = if v1.covariance != Covariance::Contravariant {
            &metric.to_contravariant(position, v1)
        } else { v1 };

        let v2 = if v2.covariance != Covariance::Contravariant {
            &metric.to_contravariant(position, v2)
        } else { v2 };
        
        let mut result = 0.0;
        for i in 0..4 {
            result += v1.vector[i] * v2.vector[i] * metric.gii(i, position);
        }
        result
    }

    /// Computes the norm of vector `v` at the given `position` on the metric.
    /// It internally takes into account the covariance of `v`.
    pub fn norm(v: &RelativisticVector, position: &RelativisticVector, metric: &impl DiagonalSphericalMetric) -> f64 {
        dot_product(v, v, position, metric).sqrt()
    }

    /// Computes the square of the norm of vector `v` at the given `position` on the metric.
    /// It internally takes into account the covariance of `v`
    pub fn squared_norm(v: &RelativisticVector, position: &RelativisticVector, metric: &impl DiagonalSphericalMetric) -> f64 {
        dot_product(v, v, position, metric)
    }

    /// Computes the angle between two vectors `v1` and `v2` at the given `position` on the metric.
    /// It internally takes into account the covariance of `v1` and `v2`.
    pub fn angle(v1: &RelativisticVector, v2: &RelativisticVector, position: &RelativisticVector, metric: &impl DiagonalSphericalMetric) -> f64 {
        dot_product(v1, v2, position, metric) / (norm(v1, position, metric) * norm(v2, position, metric))
    }

    
    // Here below we define actual metrics that implement the traits above
    
    /// Ellis wormholes are one of the first examples of traversable wormholes that were
    /// studied in the literature. See https://en.wikipedia.org/wiki/Ellis_wormhole for
    /// more information.
    /// 
    /// Ellis wormholes are uniquely chracterized by the radius parameter *rho*.
    pub struct EllisMetric {
        rho: f64,
    }

    impl EllisMetric {

        /// Creates a new `EllisMetric` with the given radius parameter.
        /// `rho`: throat radius parameter
        pub fn new(rho: f64) -> EllisMetric {

            if rho <= 0.0 {
                panic!("The rho parameter for Ellis Metrics must be positive.");
            }

            EllisMetric { rho }
        }
    }

    impl DiagonalSphericalMetric for EllisMetric {
        fn r(&self, l:f64) -> f64 { (self.rho.powi(2) + l.powi(2)).sqrt() }
        fn r_squared(&self, l:f64) -> f64 { return self.rho.powi(2) + l.powi(2) }
        fn r_derivative(&self, l:f64) -> f64 { l/self.r(l) }
    }

    /// In the movie Interstellar, a wormhole is shown. The `InterstellarMetric` is used it
    /// to implement its metric in the crate.
    /// 
    /// Such a wormhole, as described in O. James et al. (https://doi.org/10.1119/1.4916949),
    /// is characterized by three parameters:
    /// - the mass parameter `m`
    /// - the throat length parameter `a`
    /// - the throat radius parameter `rho`
    pub struct InterstellarMetric {
        m: f64,
        a: f64,
        rho: f64,
    }
    
    impl InterstellarMetric {

        /// Creates a new Interstellar Metric with the given parameters.
        /// `m`: mass parameter
        /// `a`: throat length parameter
        /// `rho`: throat radius parameter
        pub fn new(
            m: f64, a: f64, rho: f64) -> InterstellarMetric {

            if m <= 0.0 {
                panic!("The mass parameter for Interstellar Metrics must be positive.");
            }

            if a <= 0.0 {
                panic!("The angular momentum parameter for Interstellar Metrics must be positive.");
            }

            if rho <= 0.0 {
                panic!("The rho parameter for Interstellar Metrics must be positive.");
            }

            InterstellarMetric { m, a, rho }
        }

        fn scaled_distance(&self, l:f64) -> f64 { 2.0*(l.abs() - self.a)/(PI*self.m)}

    }

    impl DiagonalSphericalMetric for InterstellarMetric {
        
        fn r(&self, l:f64) -> f64 {
            if l.abs() > self.a {
                let x = self.scaled_distance(l);
                self.rho + self.m*(x*x.atan() - (1.0 + x*x).ln()/2.0)
            } else {
                self.rho
            }
        }
        
        fn r_squared(&self, l:f64) -> f64 { self.r(l).powi(2) }
        
        fn r_derivative(&self, l:f64) -> f64 { 
            if l.abs() > self.a {
                let x = self.scaled_distance(l);
                (2.0/PI)*(l.signum())*(x.atan())
            } else {
                0.0
            }
        }

    }

    
    /// Standard 3d flat Euclidean metric, in polar coordinates
    /// [Not tested]
    pub struct FlatSphericalMetric {}

    impl FlatSphericalMetric {

        pub fn new() -> FlatSphericalMetric {
            FlatSphericalMetric {}
        }
    }

    impl DiagonalSphericalMetric for FlatSphericalMetric {
        fn r(&self, l:f64) -> f64 { l }
        fn r_squared(&self, l:f64) -> f64 { l.powi(2) }
        fn r_derivative(&self, _l:f64) -> f64 { 1.0 }
    }



    #[cfg(test)]
    mod tests_metrics {

        use super::*;
        use approx::assert_relative_eq;

        #[test]
        fn test_photon_normalization_and_direction_in_plane() {

            let metric = EllisMetric { rho: 1.0 };

            let position = RelativisticVector::new(
                na::Vector4::new(0.0, 5.0, PI/2.0, 0.0),
                Covariance::Contravariant,
            );

            let angle = PI/4.0;

            let x = angle.cos();
            let y = 0.0;
            let z = angle.sin();

            let direction = na::Vector3::new(x, y, z);

            let photon = EllisMetric::new_photon(&metric, &position, direction);
            let norm = squared_norm(&photon.momentum, &position, &metric);

            let dir_from_photon = metric.relativistic_vector_to_direction(&photon.momentum, &position);

            assert_relative_eq!(norm, 0.0);
            assert_relative_eq!(direction, dir_from_photon);

        }

        #[test]
        fn test_photon_propagation_in_plane() {

            let metric = EllisMetric { rho: 1.0 };

            let position = RelativisticVector::new(
                na::Vector4::new(0.0, 5.0, PI/2.0, 0.0),
                Covariance::Contravariant,
            );

            let angle = PI/4.0;

            let x = angle.cos();
            let y = 0.0;
            let z = angle.sin();

            let direction = na::Vector3::new(x, y, z);

            let mut photon = EllisMetric::new_photon(&metric, &position, direction);

            for _ in 0..100 {
                metric.update_relativistic_object(&mut photon, 0.01);
            }

            let norm = squared_norm(&photon.momentum, &position, &metric);
            assert_relative_eq!(norm, 0.0);

        }


    }

    pub enum Metric {
        Ellis(EllisMetric),
        Interstellar(InterstellarMetric),
    }

}