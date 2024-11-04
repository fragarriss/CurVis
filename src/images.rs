/// Module for loading and saving DynamicImages from file (uses the image crate).
pub mod images {

    use std::path::PathBuf;

    /// Loads an image from a file and returns it as a DynamicImage object.
    pub fn load_image(path: &PathBuf) -> image::DynamicImage {
        let img = image::ImageReader::open(path).unwrap().decode().unwrap();

        img
    }

    pub fn save_image(img: image::DynamicImage, path: &PathBuf) {
        img.save(path).unwrap();
    }

    pub fn empty_image(width: u32, height: u32) -> image::DynamicImage {
        let img = image::DynamicImage::new_rgb8(width, height);
        img
    }
}

/// Module for managing images as the spherical backgrounds of the curved space-time.
pub mod spherical_images {

    use crate::algebra::algebra;
    use crate::images::images::load_image;
    use image::GenericImageView;
    use std::f64::consts::PI;
    // use crate::algebra;
    use nalgebra as na;
    use std::path::PathBuf;

    /// Struct that represents a spherical image.
    /// 
    /// In the context of this crate, they are used for the two backgrounds 
    /// of the wormhole space-times.
    ///
    /// `SpericalImage`s can be be rotated using an `algebra::Orientation` object.
    /// By default, images are oriented in world space such that the center
    /// of the image is aligned with the world x-axis and the top of the image
    /// corresponds to the world z-axis.
    ///
    /// The orientation specifies the forward and up vectors of the image with
    /// respect to the world space.
    ///
    /// As such, we define also the image space, and the polar coordinates (theta, phi)
    /// of a vector in the image space. This means that `theta = 0` (`PI`) will always
    /// correspond to the top (bottom) of the image, and `phi = 0` (`PI`) will always
    /// correspond to the center (back) of the image.
    pub struct SphericalImage {
        img: image::DynamicImage,
        width_pixels: u32,
        height_pixels: u32,
        orientation: algebra::Orientation,
    }

    impl SphericalImage {
        
        /// Creates a new `SphericalImage` object from a `DynamicImage`, and optionally
        /// specifies the forward and up vectors of the image in world space.
        /// 
        /// The `forward` and `up` vectors are used to define the orientation of the image
        /// in world space. The `forward` vector gives the direction of the center of the
        /// image in world space, and the `up` vector gives the direction of the top of the
        /// image in world space.
        /// 
        /// Internally, the `up` vector is made orthogonal to the `forward` vector.
        /// 
        /// By default, the forward vector is the world x-axis and the up vector is the world z-axis.
        pub fn new(
            img: image::DynamicImage,
            forward: Option<na::Vector3<f64>>,
            up: Option<na::Vector3<f64>>,
        ) -> SphericalImage {
            let (width_pixels, height_pixels) = img.dimensions();

            let orientation = algebra::Orientation::new(
                forward.unwrap_or(na::Vector3::<f64>::x()),
                up.unwrap_or(na::Vector3::<f64>::z()),
            );

            SphericalImage {
                img,
                width_pixels,
                height_pixels,
                orientation,
            }
        }

        /// Returns a reference to the current `Orientation` of the image in background space.
        pub fn orientation(&self) -> &algebra::Orientation {
            &self.orientation
        }

        /// Sets the orientation of the image in background space.
        fn set_orientation(&mut self, orientation: algebra::Orientation) {
            self.orientation = orientation;
        }

        /// Sets the orientation of the camera by specifying the `forward` and `up` vectors.
        pub fn set_forward_up(&mut self, forward: na::Vector3<f64>, up: na::Vector3<f64>) {
            self.set_orientation(algebra::Orientation::new(forward, up));
        }

        /// Gets the pixel indexed by `index_x` and `index_y`.
        pub fn get_pixel(&self, index_x: &u32, index_y: &u32) -> image::Rgba<u8> {
            let ind_x = *index_x;
            let ind_y = *index_y;
            self.img.get_pixel(ind_x, ind_y)
        }

        /// Given `theta` and `phi` polar coordinates in background space, it returns the corresponding
        /// pixel indices of the image.
        fn pixel_indexes_x_y_from_theta_phi_of_image(&self, theta: &f64, phi: &f64) -> (u32, u32) {
            let (theta, phi) = algebra::normalize_theta_phi(*theta, *phi);

            let y = ((theta / PI) * (self.height_pixels as f64)) as u32;
            let x = ((0.5 - phi / (2.0 * PI)).rem_euclid(1.0) * (self.width_pixels as f64)) as u32;
            (x, y)
        }

        /// Given `theta` and `phi` polar coordinates in background space, it returns the corresponding
        /// pixel of the image.
        pub fn get_pixel_at_theta_phi_of_image(&self, theta: &f64, phi: &f64) -> image::Rgba<u8> {
            let (index_x, index_y) = self.pixel_indexes_x_y_from_theta_phi_of_image(theta, phi);
            self.get_pixel(&index_x, &index_y)
        }

        /// Given a vector with coordinates in world space, it returns the corresponding vector
        /// with coordinates in image space.
        fn image_vector3_from_world_vector3(&self, v_world: &na::Vector3<f64>) -> na::Vector3<f64> {
            // Since the image is oriented, we need to rotate the vector from world
            // space to image space.
            // N.B. We use the inverse rotation matrix: if a vector were to be aligned
            // with the forward vector of the image in world space, it should correspond
            // to (theta, phi) = (PI/2, 0) in image space (x-axis). This is the reversed
            // rotation that defines the image space.
            self.orientation
                .inverse_rotation_matrix()
                .transform_vector(v_world)
        }

        /// Given a vector in world space, it returns the corresponding polar coordinates
        /// (theta, phi) of the vector in the image space (thus, taking into account the
        /// orientation of the image).
        /// Image space is defined with respect to world space by the forward and up vectors
        /// of the image. See the `new()`.
        ///
        /// The vector needs not be normalized.
        fn theta_phi_of_image_from_vector3(&self, v_world: &na::Vector3<f64>) -> (f64, f64) {
            // Since the image is oriented, we need to rotate the vector from world
            // space to image space.
            // N.B. We use the inverse rotation matrix: if a vector were to be aligned
            // with the forward vector of the image in world space, it should correspond
            // to (theta, phi) = (PI/2, 0) in image space (x axis). This is the reversed
            // rotation that defines the image space.
            // let w: na::Vector3<f64> = self
            //     .orientation
            //     .inverse_rotation_matrix()
            //     .transform_vector(v_world);
            let w = SphericalImage::image_vector3_from_world_vector3(self, v_world);

            // Then we calculate the polar coordinates of the vector in the image space
            // Already normalized to be in the interval [0, PI] and [0, 2*PI]
            algebra::theta_phi_from_vector3(w)
        }

        /// Given a vector in world space that specifies a direction, it returns the pixel
        /// in the image that corresponds to that direction.
        pub fn get_pixel_from_vector3(&self, v: &na::Vector3<f64>) -> image::Rgba<u8> {
            let (theta, phi) = self.theta_phi_of_image_from_vector3(v);
            self.get_pixel_at_theta_phi_of_image(&theta, &phi)
        }

        /// Given `theta` and `phi` in world coordinates, it returns the corresponding pixel
        /// in the image.
        pub fn get_pixel_from_theta_phi() {
            unimplemented!();
        }
    }

    /// Method to load an image from a file and return it as a `SphericalImage` object.
    /// By default, the image is oriented such that the center of the image is aligned
    /// with the world x-axis and the top of the image corresponds to the world z-axis.
    pub fn load_image_as_spherical_image(
        path: &PathBuf,
        forward_vector_world: Option<na::Vector3<f64>>,
        up_vector_world: Option<na::Vector3<f64>>,
    ) -> SphericalImage {
        let img = load_image(path);
        SphericalImage::new(img, forward_vector_world, up_vector_world)
    }

    /* These tests were failing. Need to double-check. */

    // #[cfg(test)]

    // mod tests_spherical_images {

    //     use super::*;
    //     use crate::images::images::{empty_image, load_image};
    //     use algebra::Orientation;
    //     use image::GenericImageView;
    //     use nalgebra as na;
    //     // use approx::assert_relative_eq;

    //     static path_test_img_1: &PathBuf= &PathBuf::from("images/test_image_1.png");
    //     // const PATH_TEST_IMG_2: &str = "images/test_image_2.png";

    //     fn load_test_image_1() -> image::DynamicImage {
    //         load_image(&PathBuf::from("images/test_image_1.png"))
    //     }

    //     fn load_test_image_2() -> image::DynamicImage {
    //         load_image(&PathBuf::from("images/test_image_2.png"))
    //     }

    //     fn create_empty_image() -> image::DynamicImage {
    //         empty_image(32, 16)
    //     }

    //     fn load_spherical_test_1(
    //         forward_vector_world: Option<na::Vector3<f64>>,
    //         up_vector_world: Option<na::Vector3<f64>>,
    //     ) -> SphericalImage {
    //         load_image_as_spherical_image(path_test_img_1, forward_vector_world, up_vector_world)
    //     }

    //     // fn load_spherical_test_2(
    //     //     forward_vector_world: Option<na::Vector3<f64>>,
    //     //     up_vector_world: Option<na::Vector3<f64>>,
    //     // ) -> SphericalImage {
    //     //     load_image_as_spherical_image(PATH_TEST_IMG_2, forward_vector_world, up_vector_world)
    //     // }

    //     fn create_empty_spherical_image() -> SphericalImage {
    //         let img = create_empty_image();
    //         SphericalImage::new(img, None, None)
    //     }

    //     fn create_empty_spherical_image_with_orientation(orientation: Orientation) -> SphericalImage {
    //         let img = create_empty_image();
    //         SphericalImage::new(img, Some(orientation.forward()), Some(orientation.up()))
    //     }



    //     #[test]
    //     fn test_load_image() {
    //         let img = load_test_image_1();
    //         assert_eq!(img.dimensions(), (32, 16));

    //         let img = load_test_image_2();
    //         assert_eq!(img.dimensions(), (32, 16));
    //     }

    //     #[test]
    //     fn test_create_empty_image() {
    //         let img = create_empty_image();
    //         assert_eq!(img.dimensions(), (32, 16));
    //     }

    //     #[test]
    //     fn test_load_spherical_image_basic() {
    //         let spherical_img = load_image_as_spherical_image(path_test_img_1, None, None);
    //         assert_eq!(spherical_img.orientation().forward(), na::Vector3::x());
    //         assert_eq!(spherical_img.orientation().up(), na::Vector3::z());
    //     }

    //     #[test]
    //     fn test_load_spherical_image_custom_orientation() {
    //         let spherical_img = load_image_as_spherical_image(
    //             path_test_img_1,
    //             Some(na::Vector3::y()),
    //             Some(na::Vector3::z()),
    //         );

    //         assert_eq!(spherical_img.orientation().forward(), na::Vector3::y());
    //         assert_eq!(spherical_img.orientation().up(), na::Vector3::z());
    //     }

    //     #[test]
    //     fn test_get_pixel() {
    //         let img_1 = load_spherical_test_1(None, None);

    //         let pixel = img_1.get_pixel(&1, &1);
    //         assert_eq!(pixel, image::Rgba([255, 255, 255, 255]));

    //         let pixel = img_1.get_pixel(&1, &0);
    //         assert_eq!(pixel, image::Rgba([127, 127, 127, 255]));

    //         let pixel = img_1.get_pixel(&0, &1);
    //         assert_eq!(pixel, image::Rgba([255, 0, 0, 255]));

    //         let pixel = img_1.get_pixel(&0, &7);
    //         assert_eq!(pixel, image::Rgba([0, 0, 0, 255]));
    //     }

    //     #[test]
    //     fn test_image_vector3_from_world_vector3_with_no_orientation() {

    //         use rand;
    //         use rand::Rng;

    //         let img_no_orientation = create_empty_spherical_image();

    //         // Without orientation, the image space is the same as the world space
    //         for _ in 0..100 {
    //             let v = na::Vector3::new(
    //                 rand::thread_rng().gen_range(-1.0..1.0),
    //                 rand::thread_rng().gen_range(-1.0..1.0),
    //                 rand::thread_rng().gen_range(-1.0..1.0),
    //             );
    //             let v_img = img_no_orientation.image_vector3_from_world_vector3(&v);
    //             assert_eq!(v_img, v);
    //         };
    //     }

    //     #[test]
    //     fn test_image_vector3_from_world_vector3() {

    //         let mut img = create_empty_spherical_image();
    //         let orientation = algebra::Orientation::new(na::Vector3::y(), na::Vector3::z());
    //         img.set_orientation(orientation);

    //         // The orientation rotates x to y and y to -x
    //         // The world vector follows the inverse rotation (x to -y and y to x)

    //         let v = na::Vector3::new(1.123, 0.0, 0.0);
    //         let v_img = img.image_vector3_from_world_vector3(&v);
    //         assert_eq!(v_img, na::Vector3::new(0.0, -1.123, 0.0));

    //         let v = na::Vector3::new(0.0, 1.123, 0.0);
    //         let v_img = img.image_vector3_from_world_vector3(&v);
    //         assert_eq!(v_img, na::Vector3::new(1.123, 0.0, 0.0));

    //         let v = na::Vector3::new(0.0, 0.0, 1.123);
    //         let v_img = img.image_vector3_from_world_vector3(&v);
    //         assert_eq!(v_img, na::Vector3::new(0.0, 0.0, 1.123));

    //         let v = na::Vector3::new(1.123, 1.123, 0.0);
    //         let v_img = img.image_vector3_from_world_vector3(&v);
    //         assert_eq!(v_img, na::Vector3::new(1.123, -1.123, 0.0));

    //         let v = na::Vector3::new(1.123, 1.123, 1.123);
    //         let v_img = img.image_vector3_from_world_vector3(&v);
    //         assert_eq!(v_img, na::Vector3::new(1.123, -1.123, 1.123));

    //     }

    //     #[test]
    //     fn test_theta_phi_form_vector3_with_no_orientation() {
            
    //         let img = create_empty_spherical_image();

    //         // positive x axis
    //         let v = na::Vector3::new(1.234, 0.0, 0.0);
    //         let (theta, phi) = img.theta_phi_of_image_from_vector3(&v);
    //         assert_eq!((theta, phi), (PI / 2.0, 0.0));

    //         // negative x axis
    //         let v = na::Vector3::new(-1.234, 0.0, 0.0);
    //         let (theta, phi) = img.theta_phi_of_image_from_vector3(&v);
    //         assert_eq!((theta, phi), (PI / 2.0, PI));

    //         // positive y axis
    //         let v = na::Vector3::new(0.0, 1.234, 0.0);
    //         let (theta, phi) = img.theta_phi_of_image_from_vector3(&v);
    //         assert_eq!((theta, phi), (PI / 2.0, PI / 2.0));

    //         // negative y axis
    //         let v = na::Vector3::new(0.0, -1.234, 0.0);
    //         let (theta, phi) = img.theta_phi_of_image_from_vector3(&v);
    //         assert_eq!((theta, phi), (PI / 2.0, 3.0 * PI / 2.0));
    //         assert_ne!(phi, -PI);

    //         // positive z axis
    //         let v = na::Vector3::new(0.0, 0.0, 1.234);
    //         let (theta, phi) = img.theta_phi_of_image_from_vector3(&v);
    //         assert_eq!((theta, phi), (0.0, 0.0));

    //         // negative z axis
    //         let v = na::Vector3::new(0.0, 0.0, -1.234);
    //         let (theta, phi) = img.theta_phi_of_image_from_vector3(&v);
    //         assert_eq!((theta, phi), (PI, 0.0));

    //         // positive x and y axis
    //         let v = na::Vector3::new(1.234, 1.234, 0.0);
    //         let (theta, phi) = img.theta_phi_of_image_from_vector3(&v);
    //         assert_eq!((theta, phi), (PI / 2.0, PI / 4.0));

    //         // negative x and y axis
    //         let v = na::Vector3::new(-1.234, -1.234, 0.0);
    //         let (theta, phi) = img.theta_phi_of_image_from_vector3(&v);
    //         assert_eq!((theta, phi), (PI / 2.0, 5.0 * PI / 4.0));

    //     }

    //     #[test]
    //     fn test_theta_phi_from_vector3 () {

    //         let or = Orientation::new(na::Vector3::y(), na::Vector3::z());
    //         let img = create_empty_spherical_image_with_orientation(or);

    //         // image is rotated such that x -> y and y -> -x
    //         // The world vector is thus rotated in the opposite direction (x -> -y and y -> x)

    //         // positive x axis
    //         let v = na::Vector3::new(1.234, 0.0, 0.0);
    //         let (theta, phi) = img.theta_phi_of_image_from_vector3(&v);
    //         assert_eq!((theta, phi), (PI / 2.0, 3.0*PI/2.0));

    //         // negative x axis
    //         let v = na::Vector3::new(-1.234, 0.0, 0.0);
    //         let (theta, phi) = img.theta_phi_of_image_from_vector3(&v);
    //         assert_eq!((theta, phi), (PI/2.0, PI/2.0));

    //         // positive y axis
    //         let v = na::Vector3::new(0.0, 1.234, 0.0);
    //         let (theta, phi) = img.theta_phi_of_image_from_vector3(&v);
    //         assert_eq!((theta, phi), (PI / 2.0, 0.0));

    //         // negative y axis
    //         let v = na::Vector3::new(0.0, -1.234, 0.0);
    //         let (theta, phi) = img.theta_phi_of_image_from_vector3(&v);
    //         assert_eq!((theta, phi), (PI / 2.0, PI));

    //         // positive z axis
    //         let v = na::Vector3::new(0.0, 0.0, 1.234);
    //         let (theta, _phi) = img.theta_phi_of_image_from_vector3(&v);
    //         assert_eq!(theta, 0.0);

    //         // negative z axis
    //         let v = na::Vector3::new(0.0, 0.0, -1.234);
    //         let (theta, _phi) = img.theta_phi_of_image_from_vector3(&v);
    //         assert_eq!(theta, PI);

    //         // Mixed
    //         let v = na::Vector3::new(1.234, 1.234, 0.0);
    //         let (theta, phi) = img.theta_phi_of_image_from_vector3(&v);
    //         assert_eq!((theta, phi), (PI / 2.0, 7.0*PI/4.0));


    //         // Should investigate! In the following test, theta evaluates to NaN.

    //         // // Mixed
    //         // let v = na::Vector3::new(1.234, 1.234, 1.234);
    //         // let (theta, phi) = img.theta_phi_of_image_from_vector3(&v);
    //         // assert_eq!((theta, phi), (1.0/(3.0_f64.sqrt()).acos(), 7.0*PI/4.0));

    //     }
    // }
}
