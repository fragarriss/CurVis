pub mod algebra {

    use nalgebra as na;
    use std::f64::consts::PI;

    /// This struct represents an orientation in 3D space, defined by a pair of
    /// forward and up vectors.
    pub struct Orientation {
        forward: na::Vector3<f64>,
        up: na::Vector3<f64>,
        rotation_matrix: na::Rotation3<f64>,
        inverse_rotation_matrix: na::Rotation3<f64>,
    }

    impl Orientation {
        pub fn new(forward: na::Vector3<f64>, up: na::Vector3<f64>) -> Orientation {

            // Checking if parallel
            if forward.cross(&up).norm() == 0.0 {
                panic!("Forward and up vectors must not be parallel");
            }

            let rotation = rotation_matrix_from_forward_up_pairs(
                na::Vector3::<f64>::x(),
                na::Vector3::<f64>::z(),
                forward,
                up);
            let inverse_rotation = rotation.inverse();

            let orthogonal_up = rotation*na::Vector3::<f64>::z();

            Orientation {
                forward: forward,
                up: orthogonal_up,
                rotation_matrix: rotation,
                inverse_rotation_matrix: inverse_rotation,
            }
        }

        pub fn forward(&self) -> na::Vector3<f64> {
            self.forward
        }

        pub fn up(&self) -> na::Vector3<f64> {
            self.up
        }

        pub fn rotation_matrix(&self) -> na::Rotation3<f64> {
            self.rotation_matrix
        }

        pub fn inverse_rotation_matrix(&self) -> na::Rotation3<f64> {
            self.inverse_rotation_matrix
        }

    }


    /// Evaluates the rotation matrix that rotates a pair of forward-up vectors
    /// to a new pair of such vectors.
    /// 
    /// Vectors need not be normalized, and up vectors must not be collinear
    /// with their respective forward vectors.
    pub fn rotation_matrix_from_forward_up_pairs(
            forward_previous: na::Vector3<f64>,
            up_previous: na::Vector3<f64>,
            forward_next: na::Vector3<f64>,
            up_next: na::Vector3<f64>,
        ) -> na::Rotation3<f64> {
        
            let rotation1 = na::Rotation3::face_towards(&forward_previous, &up_previous);
            let rotation2 = na::Rotation3::face_towards(&forward_next, &up_next);
            rotation2*rotation1.inverse()
        }

    /// Evaluates the rotation amtrix that rotates the x axis to the direction
    /// defined by the angles theta and phi.
    /// In particular, the rotation is equivalent to first rotating by an angle
    /// (PI/2 - theta) around the negative y axis, and then to rotating by
    /// and angle phi around the z axis.

    pub fn rotation_matrix_from_theta_phi(theta: f64, phi: f64) -> na::Rotation3<f64> {

        let (theta, phi) = normalize_theta_phi(theta, phi);
        
        let rotation1 = na::Rotation3::new(na::Vector3::new(0.0, theta - PI/2.0, 0.0) );
        let rotation2 = na::Rotation3::new(na::Vector3::new(0.0, 0.0, phi) );

        rotation2*rotation1
    }

    pub fn rotation_from_two_vectors(v1: &na::Vector3<f64>, v2: &na::Vector3<f64>) -> na::Rotation3<f64> {
        
        // Checking if parallel
        if v1.cross(&v2).norm() == 0.0 {
            panic!("v1 and v2 must not be parallel");
        }

        // Rotation that rotates v1 to v2
        na::Rotation3::rotation_between(v1, v2).unwrap()
    }

    /// Ensures theta is always between 0 and PI, and phi is always between 0 and 2*PI
    /// When theta is negative, it is made positive (abs) and phi is increased by PI;
    /// then phi is taken modulo 2*PI.
    pub fn normalize_theta_phi(theta: f64, phi: f64) -> (f64, f64) {

        let (theta, phi) = if theta < 0.0 {
            (theta.abs(), phi + PI)
        } else {
            (theta, phi)
        };

        let phi = phi.rem_euclid(2.0*PI);
        (theta, phi)
    }

    pub fn vector3_from_theta_phi(theta: f64, phi: f64) -> na::Vector3<f64> {

        let (theta, phi) = normalize_theta_phi(theta, phi);

        let x = theta.sin()*phi.cos();
        let y = theta.sin()*phi.sin();
        let z = theta.cos();
        na::Vector3::new(x, y, z)
    }

    pub fn theta_phi_from_vector3(v: na::Vector3<f64>) -> (f64, f64) {
        let r = v.norm();
        let theta = (v.z/r).acos();
        let phi = v.y.atan2(v.x);
        
        normalize_theta_phi(theta, phi)
    }

    #[cfg(test)]
    mod tests_algebra {

        use super::*;
        use nalgebra as na;
        use std::f64::consts::PI;
        use approx::assert_relative_eq;

        #[test]
        fn test_orientation_constructor() {
            let forward = na::Vector3::new(1.0, 0.0, 0.0);
            let up = na::Vector3::new(0.0, 0.0, 1.0);
            let orientation = Orientation::new(forward, up);

            assert_eq!(orientation.forward(), forward);
            assert_eq!(orientation.up(), up);
        }

        #[test]
        fn test_orientation_constructor_non_orthogonal() {
            let forward_1 = na::Vector3::new(1.0, 0.0, 0.0);
            let up_1 = na::Vector3::new(1.0, 0.0, 1.0);
            let orientation_1 = Orientation::new(forward_1, up_1);

            let forward_2 = na::Vector3::new(1.0, 1.0, 0.0);
            let up_2 = na::Vector3::new(-1.0, -1.0, 1.0);
            let orientation_2 = Orientation::new(forward_2, up_2);

            let forward_3 = na::Vector3::new(1.0, 0.0, 1.0);
            let up_3 = na::Vector3::new(1.0, 1.0, 1.0);
            let orientation_3 = Orientation::new(forward_3, up_3);

            assert_eq!(orientation_1.forward(), forward_1);
            assert_eq!(orientation_1.up(), na::Vector3::new(0.0, 0.0, 1.0));

            assert_eq!(orientation_2.forward(), forward_2);
            assert_eq!(orientation_2.up(), na::Vector3::new(0.0, 0.0, 1.0));

            assert_eq!(orientation_3.forward(), forward_3);
            assert_eq!(orientation_3.up(), na::Vector3::new(0.0, 1.0, 0.0));
        }

        #[test]
        #[should_panic]
        fn test_orientation_constructor_parallel() {
            let forward = na::Vector3::new(1.0, 0.0, 0.0);
            let up = na::Vector3::new(1.0, 0.0, 0.0);
            let _orientation = Orientation::new(forward, up);

            println!("up: {:?}", _orientation.up());
            println!("forward {:?}", _orientation.forward());
        }

        #[test]
        #[should_panic]
        fn test_orientation_constructor_antiparallel() {
            let forward = na::Vector3::new(1.0, 0.0, 0.0);
            let up = na::Vector3::new(-1.0, 0.0, 0.0);
            let _orientation = Orientation::new(forward, up);

            println!("up: {:?}", _orientation.up());
            println!("forward {:?}", _orientation.forward());
        }

        #[test]
        fn test_rotation_matrix_xz() {
            let forward = na::Vector3::new(1.0, 0.0, 0.0);
            let up = na::Vector3::new(0.0, 0.0, 1.0);
            let orientation = Orientation::new(forward, up);

            let rotation_matrix = orientation.rotation_matrix();

            assert_eq!(rotation_matrix, na::Rotation3::identity());
        }

        #[test]
        fn test_rotation_matrix_identity() {

            use rand;
            use rand::Rng;

            let forward_x = rand::thread_rng().gen_range(-1.0..1.0);
            let forward_y = rand::thread_rng().gen_range(-1.0..1.0);
            let forward_z = rand::thread_rng().gen_range(-1.0..1.0);

            let up_x = rand::thread_rng().gen_range(-1.0..1.0);
            let up_y = rand::thread_rng().gen_range(-1.0..1.0);
            let up_z = rand::thread_rng().gen_range(-1.0..1.0);

            let forward = na::Vector3::new(forward_x, forward_y, forward_z);
            let up = na::Vector3::new(up_x, up_y, up_z);
            let orientation = Orientation::new(forward, up);

            let rotation_matrix = orientation.rotation_matrix();
            let inverse_rotation_matrix = orientation.inverse_rotation_matrix();

            let identity = rotation_matrix*inverse_rotation_matrix;

            assert_relative_eq!(identity, na::Rotation3::identity());
        }

        #[test]
        fn test_rotation_matrix_from_theta_phi() {

            use rand;
            use rand::Rng;

            for _ in 0..1000 {
                let theta = rand::thread_rng().gen_range(0.0 .. PI);
                let phi = rand::thread_rng().gen_range(0.0 .. 2.0*PI);

                let rotation = rotation_matrix_from_theta_phi(theta, phi);

                let target = na::Vector3::new(theta.sin()*phi.cos(), theta.sin()*phi.sin(), theta.cos());

                let rotated = rotation*na::Vector3::x();

                let e = 2e12*f64::EPSILON;
                assert_relative_eq!(rotated, target, epsilon = e);
            }

        }

        #[test]
        fn test_vector3_from_theta_phi() {

            let inv_sqrt_2: f64 = 1.0/(2.0_f64).sqrt();

            assert_relative_eq!(vector3_from_theta_phi(0.0, 0.0), na::Vector3::new(0.0, 0.0, 1.0));
            assert_relative_eq!(vector3_from_theta_phi(PI/2.0, 0.0), na::Vector3::new(1.0, 0.0, 0.0));
            assert_relative_eq!(vector3_from_theta_phi(PI, 0.0), na::Vector3::new(0.0, 0.0, -1.0));

            assert_relative_eq!(vector3_from_theta_phi(PI/2.0, PI/4.0), na::Vector3::new(inv_sqrt_2, inv_sqrt_2, 0.0));
            assert_relative_eq!(vector3_from_theta_phi(-PI/2.0, PI/4.0), na::Vector3::new(-inv_sqrt_2, -inv_sqrt_2, 0.0));
            assert_relative_eq!(vector3_from_theta_phi(PI/2.0, -PI/4.0), na::Vector3::new(inv_sqrt_2, -inv_sqrt_2, 0.0));
            assert_relative_eq!(vector3_from_theta_phi(-PI/2.0, -PI/4.0), na::Vector3::new(-inv_sqrt_2, inv_sqrt_2, 0.0));

            assert_relative_eq!(vector3_from_theta_phi(PI/2.0, PI/2.0), na::Vector3::new(0.0, 1.0, 0.0));
            assert_relative_eq!(vector3_from_theta_phi(-PI/2.0, PI/2.0), na::Vector3::new(0.0, -1.0, 0.0));

            assert_relative_eq!(vector3_from_theta_phi(PI/2.0, 3.0*PI/4.0), na::Vector3::new(-inv_sqrt_2, inv_sqrt_2, 0.0));
            assert_relative_eq!(vector3_from_theta_phi(PI/2.0, PI), na::Vector3::new(-1.0, 0.0, 0.0));
            assert_relative_eq!(vector3_from_theta_phi(PI/2.0, 5.0*PI/4.0), na::Vector3::new(-inv_sqrt_2, -inv_sqrt_2, 0.0));
            assert_relative_eq!(vector3_from_theta_phi(PI/2.0, 3.0*PI/2.0), na::Vector3::new(0.0, -1.0, 0.0));
            assert_relative_eq!(vector3_from_theta_phi(PI/2.0, 7.0*PI/4.0), na::Vector3::new(inv_sqrt_2, -inv_sqrt_2, 0.0));

        }

        #[test]
        fn test_theta_phi_from_vector3() {

            use rand;
            use rand::Rng;

            for _ in 0..1000 {
                let theta = rand::thread_rng().gen_range(0.0 .. PI);
                let phi = rand::thread_rng().gen_range(0.0 .. 2.0*PI);

                let r = rand::thread_rng().gen_range(0.1 .. 5.0);
                let x = r*theta.sin()*phi.cos();
                let y = r*theta.sin()*phi.sin();
                let z = r*theta.cos();

                let v = na::Vector3::new(x, y, z);

                let (theta_, phi_) = theta_phi_from_vector3(v);

                let e = 2e12*f64::EPSILON;
                assert_relative_eq!(theta, theta_, epsilon = e);
                assert_relative_eq!(phi, phi_, epsilon = e);
            }


        }


    }

}
