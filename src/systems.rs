/// Within this module we find the implementation of the RelativisticSystem struct, which is the main
/// struct that performs the raytracing on the metric to render the images of the wormhole.
/// 
/// The `RelativisticSystem` struct essentially constructs a scene to be rendered by putting together
/// a `Camera` (the viewpoint), two 'SphericalImage's (the backgrounds for the two side of the wormhole)
/// and a `DiagonalSphericalMetric` (the metric that defines the wormhole).
/// Then, it uses the methods of the camera to construct light rays that are propagated over the metric
/// until they reach one of the two backgrounds.
/// 
/// In this module we also find the implementation of the `PhotonEscape` and `EscapeAngle` enums, which
/// are used to represent the state of a photon and its eventual escape from the spacetime.

pub mod systems {
    use core::f64;
    use std::f64::consts::PI;
    // use std::io::Write;

    use image::GenericImage;
    use nalgebra as na;
    
    use crate::images::spherical_images::SphericalImage;
    use crate::cameras::cameras::Camera;
    use crate::metrics::metrics::DiagonalSphericalMetric;
    use crate::vectors::vectors::{RelativisticObject, RelativisticVector, Covariance};
    use crate::algebra::algebra::{vector3_from_theta_phi, rotation_from_two_vectors};
    use crate::sampling::sampling::doubly_sample_function;


    /// Represents the state of a photon as it propagates in the metric.
    /// `PositiveEscape` and `NegativeEscape` represent the photon after it has
    /// escaped the metric, while `NotEscaped` represents the photon if it has not.
    /// 
    /// The `PhotonEscape` and `NegativeEscape` variants contain the photon itself
    /// because its state is needed to calculate the direction in which it is propagating
    /// in background space after the escape.
    /// 
    /// In the future the `NotEscaped` variant may contain the photon itself as well, but
    /// for now it is not needed.
    #[derive(Clone)]
    enum PhotonEscape {
        PositiveSpace(RelativisticObject),
        NegativeSpace(RelativisticObject),
        NotEscaped,
    }

    /// Given a photon that has propagated on the equatorial plane of a spherically symmetric
    /// space-time that is asymptotically flat, the escape angle is defined as the angle between the photon's
    /// direction and the x-axis in world space coordinates.
    /// 
    /// This enum keeps track of the escape angle of the photon, and whether it has escaped
    /// in the positive or negative spaces at the sides of the wormhole (or if it has not escaped,
    /// in which case the angle is not tracked).
    #[derive(Clone)]
    pub enum EscapeAngle {
        PositiveSpace(f64),
        NegativeSpace(f64),
        NotEscaped,
    }

    /// `RelativisticSystem` is a struct that contains all the information needed to render a scene
    /// in a relativistic spacetime.
    /// 
    /// # Fields
    /// - `metric`: the metric that defines the spacetime. Must implement the `DiagonalSphericalMetric` trait.
    /// - `background_positive`: the image that is used as background for the positive-l space.
    /// - `background_negative`: the image that is used as background for the negative-l space.
    /// - `camera`: the camera that represents the viewpoint from which the scene is rendered.
    pub struct RelativisticSystem<M: DiagonalSphericalMetric> {
        pub metric: M,
        pub background_positive: SphericalImage,
        pub background_negative: SphericalImage,
        pub camera: Camera,
    }
    
    /// Given a photon, it evaluates it future path on the metric and returns it as a vector
    /// of `RelativisticObject`s.
    pub fn compute_photon_trajectory (
        photon: &mut RelativisticObject,
        metric: &impl DiagonalSphericalMetric,
        iterations: u32,
        delta: f64,
    ) -> Vec<RelativisticObject> {
        let mut trajectory = Vec::new();

        for _ in 0..iterations {
            trajectory.push(photon.clone());
            metric.update_relativistic_object(photon, delta);
        }

        return trajectory;

    }

    /// Given a photon, this function propagates it in the metric until it escapes the spacetime
    /// region defined by the maximum radius.
    /// 
    /// More formally, the photon is considered escaped if the absolute value of its radial
    /// coordinate (l) exceeds the maximum radius passed as argument.
    /// 
    /// The function returns a `PhotonEscape` enum that contains the photon if it has escaped
    /// and in which space it has escaped (positive or negative).
    /// 
    /// If the photon has not escaped after the maximum number of iterations, the function
    /// returns `PhotonEscape::NotEscaped`.
    /// 
    /// # Arguments
    /// - `metric`: the metric in which the photon is propagating.
    /// - `photon`: the photon to propagate.
    /// - `delta`: the step size for the numerical integration (proper time).
    /// - `max_iterations`: the maximum number of iterations to perform.
    /// - `max_radius`: the radius beyond which a photon escapes from the metric.
    /// 
    /// # Panics
    /// - If the photon starting position is already boyond the maximum radius.
    fn escape_photon(
        metric: &impl DiagonalSphericalMetric,
        photon: &mut RelativisticObject,
        delta: f64,
        max_iterations: u32,
        max_radius: f64) -> PhotonEscape {

        if photon.x(1).abs() > max_radius {
            panic!("Photon already beyond the maximum radius. Cannot evaluate escape.")
        }

        for _ in 0..max_iterations {
            metric.update_relativistic_object(photon, delta);

            if photon.x(1) > max_radius {
                return PhotonEscape::PositiveSpace(photon.clone());
            }
            else if photon.x(1) < -max_radius {
                return PhotonEscape::NegativeSpace(photon.clone());
            }
        }

        return PhotonEscape::NotEscaped;

    }


    /// Given a photon that has escaped, this function returns a 3D vector representing
    /// the direction the photon is propagating in world space coordinates.
    fn escaped_photon_to_world_direction(metric: &impl DiagonalSphericalMetric, escape: &PhotonEscape) -> na::Vector3<f64> {
        
        /* The photon direction can be calcualted from the photon's momentum in the euclidean
        tangent space of the metric at the photon's position.
        
        The tangent euclidean space will be in general rotated with respect to the world
        euclidean space (the space in which the background images are defined, which is
        the euclidean appoximation of the curved metric at large radii).
        
        The world direction we are interested in can be calculated by rotating the photon's
        direction in the tangent space by the angle between the tangent space and the world
        space.
        */
        
        // WARNING: What if the photon is in the negative space? (l < 0) It may be that
        // the tangent space direction should be flipped in that case.

        let photon = match escape {
            PhotonEscape::PositiveSpace(photon) => {
                photon
            }
            PhotonEscape::NegativeSpace(photon) => {
                photon
            }
            PhotonEscape::NotEscaped => { panic!("Photon not escaped. Cannot evaluate world direction.") }
        };

        // Direction of the photon on its tangent space
        let photon_tangent_space_direction: na::Vector3<f64> = 
        DiagonalSphericalMetric::relativistic_vector_to_direction(metric, &photon.momentum, &photon.position);

        let photon_theta = photon.x(2);
        let photon_phi = photon.x(3);

        let world_position = vector3_from_theta_phi(photon_theta, photon_phi);

        let tangent_space_rotation = rotation_from_two_vectors(
            &na::Vector3::x(),
            &world_position
        );

        return tangent_space_rotation * photon_tangent_space_direction;

    }

    /// This function propagates a photon starting at position `(0, l, PI/2, 0)` on the metric,
    /// with direction `alpha` in the equatorial plane (on the tangent space of the metric)
    /// until it escapes; it then returns the angle at which the photon is escaped in the world space
    /// as an `EscapeAngle` enum (see the definition of `EscapeAngle` for more details).
    /// 
    /// The photon is considered escaped when the absolute value of its radial coordinate (`l`) exceeds
    /// the `max_radius` argument. The sign of `l` determines in which space the photon has escaped.
    /// 
    /// Arguments:
    /// - `l`: initial radial coordinate of the photon.
    /// - `alpha`: initial angle between the photon direction and the radial direction.
    /// - `delta`: step size for the numerical integration (proper time).
    /// - `max_iterations`: maximum number of iterations to evaluate the escape.
    /// - `max_radius`: the radius beyond which a photon escapes from the metric.
    pub fn compute_escape_angle(
        metric: &impl DiagonalSphericalMetric,
        l:f64,
        alpha:f64,
        delta:f64,
        max_iterations:u32,
        max_radius:f64) -> EscapeAngle {
        
        // Here we are going to propagate the photon in the equatorial plane of the metric.
        // In such plane, the tangential euclidean space is defined by the frame field of
        // the metric. In Curvis, we choose the three basis vectors x, y, and z to be aligned
        // respectively with the direction of increasing coordinates l, theta and phi.

        // Hence, the direction of an object moving in the equatorial plane will have a 
        // vanishing y component (y is the direction of increasing theta).

        // Here we define the initial direction of the photon in the equatorial plane.
        // (positive x -> increasing l; positive z -> increasing phi).
        let direction_world = na::Vector3::new(alpha.cos(), 0.0, alpha.sin());

        // Here we define the initial position of the photon in the equatorial plane (theta = PI/2.0).
        let position = RelativisticVector::new(
            na::Vector4::new(0.0, l, PI/2.0, 0.0),
            Covariance::Contravariant
        );

        // Here we instantiate the photon
        let mut photon = metric.new_photon(&position, direction_world);

        // And here we propagate until it is escaped (or until the maximum number of iterations is reached).
        let escape = escape_photon(metric, &mut photon, delta, max_iterations, max_radius);

        // Once escaped, we calculate the angle at which the photon is escaped, defined as
        // the angle of the photon's momentum in the tangential euclidean space.

        match escape {
            
            PhotonEscape::NotEscaped => {
                return EscapeAngle::NotEscaped;
            }
            _ => {}
        }

        let mut world_direction: na::Vector3<f64> = escaped_photon_to_world_direction(metric, &escape);
        world_direction.normalize_mut();

        let vx= world_direction.dot(&na::Vector3::x());
        let vy= world_direction.dot(&na::Vector3::y());

        let angle = if vy >= 0.0 { vx.acos() } else { 2.0*PI - vx.acos() };
        
        match escape {

            PhotonEscape::PositiveSpace(_) => { return EscapeAngle::PositiveSpace(angle) }
            PhotonEscape::NegativeSpace(_) => { return EscapeAngle::NegativeSpace(angle) }
            _ => panic!("Escape state should be PositiveSpace or NegativeSpace")
        }
        
    }

    /// Works like `compute_escape_angle()`, but for a range of initial angles.
    /// See `compute_escape_angle()` for more details about the arguments.
    pub fn compute_escape_angles_range(
        metric: &impl DiagonalSphericalMetric,
        l:f64,
        alphas: &Vec<f64>,
        delta:f64,
        max_iterations:u32,
        max_radius:f64) -> Vec<EscapeAngle> {

        let mut escapes: Vec<EscapeAngle> = Vec::new();

        for alpha in alphas {
            escapes.push(compute_escape_angle(metric, l, *alpha, delta, max_iterations, max_radius));
        }

        return escapes;

    }

    impl <M: DiagonalSphericalMetric> RelativisticSystem<M> 
    {

        pub fn new(metric: M, background_positive: SphericalImage, background_negative: SphericalImage, camera: Camera) -> RelativisticSystem<M> {
            RelativisticSystem { metric, background_positive, background_negative, camera }
        }

        pub fn metric(&self) -> &M {
            &self.metric
        }

        pub fn background_positive(&self) -> &SphericalImage {
            &self.background_positive
        }

        pub fn background_negative(&self) -> &SphericalImage {
            &self.background_negative
        }

        pub fn camera(&self) -> &Camera {
            &self.camera
        }

        /// The main method that renders the scene.
        pub fn render_image(
            &self,
            max_iterations: u32,
            max_radius: f64,
            delta: f64,
        ) -> image::DynamicImage {

            let mut img = image::DynamicImage::new_rgb8(self.camera.resolution_width(), self.camera.resolution_height());

            for i in 0..self.camera.resolution_width() {
                
                println!("Rendering column {} of {}", i, self.camera.resolution_width());
                
                for j in 0..self.camera.resolution_height() {
                    let mut photon = self.camera_pixels_x_y_to_photon(i, j);
                    let escape = escape_photon(&self.metric, &mut photon, delta, max_iterations, max_radius);
                    let pixel = self.photon_escape_to_pixel(escape);
                    img.put_pixel(i, j, pixel)
                }
            };

            return img;

        }

        /// This method renders the scene in a more efficient way compared to `render_image()`.
        pub fn render_image_efficient(
            &self,
            max_iterations_propagation: u32,
            max_radius: f64,
            delta: f64,
            alpha_nums: u32,
            max_iterations_sampling: u32,
            sampling_convergence_threshold_1: f64,
            sampling_convergence_threshold_2: f64,
            // data_output_path: Option<&str>,
        ) -> image::DynamicImage {

            /* This function expoints the fact that in spherically symmetric space-times, the propagation on the
            equatorial plane is equivalent to any other plane.

            This is an advantage for two reasons:
            - The polar coordinates chosen for the metric behave badly at the poles
            - Once we know the path of a photon that starts at a given angle alpha from the center of the metric, we can
            re-use this information for any other photon that starts at the same angle from the center.

            The first point is important to avoid image artifacts, while the second massively reduces the computational
            cost of the rendering. On a laptop render_image_efficient can reduce the rendering times of an HD frame
            from minutes to seconds.

            Furthermore, this method uses the sampling algorithm implemented in the sampling module to further reduce
            the computation time. In fact, this allows to evaluate more trajectories where they end up in regions
            of the metric where the escape angle varies more rapidly (e.g., when a photon passes very close to the
            wormhole-s edge).

            The function works as follows:

            Step 1. For each pixel in the camera image, we calculate the corresponding world vector direction.

            Step 2. For each of these directions, we calculate the initial alpha angles.
            This is done by rotating each vector on the equatorial plane and taking the angle between the rotated vector
            and the axis joining the camera to the center of the metric.

            Step 3 (sampling). We choose another range of alphas between 0 and PI and we use the sampling algorithm to evaluate
            the escape angle and escape space for each of them. This steps returns three vectors of f64: the alphas,
            the escape angles and the escape spaces. The length of these vectors can be significantly higher than the
            initial number of alphas.

            Step 4. We interpolate the vectors calculated at step 3 and apply them to the alphas calculated at step 2.

            At this point, we know the escape angle that each pixel in the image would have if it propagated on the
            equatorial plane of the metric.

            Step 5. We apply the appropriate rotations to determine the actual point in the background images where the
            photons end up.
            
            (I didn't bother splitting this function is sub-functions as so far it is the only place where this sequence
            of steps is applied. In the future I might split it if needed.) */

            // Initial empty image
            let mut img = image::DynamicImage::new_rgb8(self.camera.resolution_width(), self.camera.resolution_height());

            // Step 1. For each pixel in the camera image, we calculate the corresponding world vector
            // direction and store it in a linear vector of na::Vector3<f64> with dimensions (width * height).
            
            // Calculating the direction of the camera position on the background space
            let camera_theta = self.camera.position().v(2);
            let camera_phi = self.camera.position().v(3);

            // This joins the center of the metric to the initial positions of the photons
            let camera_position_on_bg_space = vector3_from_theta_phi(camera_theta, camera_phi);

            let mut start_dirs_on_tangent_space: Vec<na::Vector3<f64>> = Vec::new();
            let mut start_dirs_on_bg_space: Vec<na::Vector3<f64>> = Vec::new();
            let mut rotation_axes_on_bg_space: Vec<na::Vector3<f64>> = Vec::new();

            // Step 2

            for i in 0..self.camera.resolution_width() {
                for j in 0..self.camera.resolution_height() {
                    
                    // outward_vector_on_world_space_from_x_y() returns the direction of the pixel with coordinates
                    // specified on the camera space (euclidean tangential space to the metric on the camera position)
                    let outward_vector_on_tangent_space: na::Vector3<f64> = self.camera.outward_vector_on_world_space_from_x_y(i, j);
                    let outward_vector_on_bg_space: na::Vector3<f64> = rotation_from_two_vectors(&na::Vector3::<f64>::x(), &camera_position_on_bg_space)*outward_vector_on_tangent_space;
                    
                    
                    let rotation_axis = 
                        camera_position_on_bg_space
                        .cross(&outward_vector_on_bg_space) as na::Vector3<f64>;

                    start_dirs_on_tangent_space.push(outward_vector_on_tangent_space);
                    start_dirs_on_bg_space.push(outward_vector_on_bg_space);
                    rotation_axes_on_bg_space.push(rotation_axis);
                    
                }
            }
            

            // Then for each of them we calculate the angle between the vector and the local x axis
            let mut img_alphas: Vec<f64> = Vec::new();

            for direction in start_dirs_on_tangent_space {
                let x = na::Vector3::x();
                let angle = direction.dot(&x).acos();
                img_alphas.push(angle);
            }

            // Step 3. Sampling of alpha angles.

            let min_alpha = -0.1*PI;
            let max_alpha = 1.1*PI;

            /* Before choosing the -0.1*PI/1,1*PI range, I used to calculate the min and max alpha
            angles from the image alphas and then use them as the range for the sampling.
            See the commented lines below.

            However, for some reason I did not understnd,, the max alpha calculated this way
            was very unstable under small camera movements.
            The previous apprach would be preferable, because it would allow a better sampling
            performance with narrow field of view images. */

            // I calculate the min and max angles

            // let min_alpha: f64 = img_alphas.iter().fold(PI, |a, &b| a.min(b));
            // let max_alpha: f64 = img_alphas.iter().fold(0.0, |a, &b| a.max(b));

            // println!("min alpha: {}", &min_alpha);
            // println!("max alpha: {}", &max_alpha);

            // I consider here compute_escape_angle as a function that needs to be sampled.
            // `doubly_sample_function()` implements an algorithm that samples the function more
            // densely where it varies more rapidly. See the sampling module for more info
            // on the algorithm.
            
            let (
                alphas_sampling,
                escape_angles_sampling,
                escape_spaces_sampling) = doubly_sample_function(
                min_alpha,
                max_alpha,
                alpha_nums as usize,
                max_iterations_sampling as usize,
                sampling_convergence_threshold_1,
                sampling_convergence_threshold_2,

                |alpha| {
                    match compute_escape_angle(
                        &self.metric, 
                        self.camera.position().v(1),
                        alpha,
                        delta,
                        max_iterations_propagation,
                        max_radius) {
                        EscapeAngle::PositiveSpace(e) => (e, 1.0),
                        EscapeAngle::NegativeSpace(e) => (e, -1.0),
                        EscapeAngle::NotEscaped => (f64::NAN, f64::NAN)
                    }
                }
            );

            // Step 4. Interpolation of the escape angles and spaces to the image alphas

            // I map the image alphas to the image excape angles
            let img_escape_angles = interp::interp_slice(&alphas_sampling, &escape_angles_sampling, &img_alphas);
            // I map the image alphas to the image escape spaces
            let img_escape_spaces = interp::interp_slice(&alphas_sampling, &escape_spaces_sampling, &img_alphas);


            // Step 5. I evaluate the final directions of the photons after the escape

            let mut final_dirs_on_bg_space: Vec<na::Vector3<f64>> = Vec::new();

            for (img_escape_angle, rotation_axis)
                in img_escape_angles
                    .iter()
                    .zip(rotation_axes_on_bg_space.iter()) {
                let rotation = na::Rotation3::from_axis_angle(&na::Unit::new_normalize(*rotation_axis), *img_escape_angle);
                final_dirs_on_bg_space.push(rotation * camera_position_on_bg_space);
            }

            // I map the final directions to the background images to assemble the final image

            for (index, (final_direction, escape_space))
                in final_dirs_on_bg_space.iter().zip(img_escape_spaces.iter()).enumerate() {

                let i = index as u32 / self.camera.resolution_height();
                let j = index as u32 % self.camera.resolution_height();

                let pixel = match escape_space {
                    1.0 => self.background_positive.get_pixel_from_vector3(final_direction),
                    -1.0 => self.background_negative.get_pixel_from_vector3(final_direction),
                    _ => image::Rgba([0, 0, 0, 255]), // Black
                };
                // let pixel = self.background_positive.get_pixel_from_vector3(final_direction);
                img.put_pixel(i, j, pixel);
            }

            return img;

        }

        /// Given a pixel in the camera image, this function returns the corresponding
        /// photon on the metric.
        fn camera_pixels_x_y_to_photon(&self, pixel_x: u32, pixel_y: u32) -> RelativisticObject {
            let direction = self.camera.outward_vector_on_world_space_from_x_y(pixel_x, pixel_y);
            return M::new_photon(&self.metric, self.camera.position(), direction)
        }
        
        /// Given an escaped photon, this function returns the appropriate pixel from
        /// the background images defined for the system.
        /// 
        /// A not escaped photon is returned as a black pixel.
        fn photon_escape_to_pixel(&self, photon_escape: PhotonEscape) -> image::Rgba<u8> {

            match photon_escape {

                PhotonEscape::PositiveSpace(photon) => {
                    return self.background_positive.get_pixel_from_vector3(
                        &self.metric.relativistic_vector_to_direction(&photon.momentum, &photon.position)
                    );
                }

                PhotonEscape::NegativeSpace(photon) => {
                    return self.background_negative.get_pixel_from_vector3(
                        &self.metric.relativistic_vector_to_direction(&photon.momentum, &photon.position)
                    );
                }

                PhotonEscape::NotEscaped => {
                    return image::Rgba([0, 0, 0, 255]); // black
                }
            }

        }

    }

}