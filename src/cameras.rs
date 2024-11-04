/// This module contains the definition of the `Camera` struct, which represents a viewpoint in
/// curved spacetime. `Camera`s also contain information to render images, such as the resolution
/// and focal length.
pub mod cameras {

    use nalgebra as na;
    use crate::algebra::algebra;
    use crate::vectors::vectors::{Covariance, RelativisticVector};

    /// This struct represents a camera in 3D space, defining the viewpoint
    /// and how an image should be rendered.
    /// 
    /// The camera position and orientation in space is defined with three
    /// vectors:
    /// - `position`: the position of the camera in the world space
    /// - `forward`: the direction the camera is pointing at
    /// - `up`: the direction that is considered "up" for the camera
    /// 
    /// The `position` vector is a contravariant `RelativisticVector`
    /// that lives on the metric, and thus has the same coordinates of the
    /// metric (for Ellis and Interstellar metrics, here we use time + 
    /// the three polar coordinates l, theta and phi).
    /// 
    /// The `forward` and `up` vectors are 3D vectors that live in the
    /// local tangent space of the camera. They are used to define the
    /// orientation of the camera and have no time coordinate.
    ///  
    /// See the `new` method for more information on the parameters that
    /// can be specified.
    pub struct Camera {
        position: RelativisticVector,
        // forward: na::Vector3<f64>,
        // up: na::Vector3<f64>,
        orientation: algebra::Orientation,
        focal_length: f64,
        sensor_width: f64,
        sensor_height: f64,
        // sensor_diagonal: f64,
        // aspect_ratio: f64,
        resolution_width: u32,
        resolution_height: u32,
        camera_to_world_rotation_matrix: na::Rotation3<f64>,
    }

    impl Camera {

        /// Creates a new camera.
        /// 
        /// # Arguments
        ///  - `position`: the position of the camera in the world space
        ///     (with the same coordinates of the metric).
        ///  - `forward`: the direction the camera is pointing at (with
        ///     coordinates on tangent space). The vector is internally
        ///     normalized.
        ///  - `up`: the direction that is considered "up" for the camera
        ///     (with coordinates on tangent space). The vector is internally
        ///     normalized and made orthogonal to the forward vector.
        ///  - `focal_length`: controls the zoom of the camera.
        ///  - `sensor_diagonal`: the diagonal of the sensor in the camera.
        ///  - `resolution_width`: the width in pixels of the images produced
        ///    using this camera.
        ///  - `resolution_height`: the height in pixels of the images produced
        ///    using this camera.
        /// 
        /// The `position` vector is a contravariant `RelativisticVector`
        /// that lives on the metric, and thus has the same coordinates of the
        /// metric (for Ellis and Interstellar metrics, here we use time + 
        /// the three polar coordinates l, theta and phi).
        /// 
        /// The `forward` and `up` vectors are 3D vectors that live in the
        /// local tangent space of the camera. They are used to define the
        /// orientation of the camera and have no time coordinate.
        /// 
        /// # Panics
        /// - If the `position` vector is not contravariant.
        /// - If `focal_length` is less than or equal to 0.
        /// - If `sensor_diagonal` is less than or equal to 0.
        /// - If `resolution_width` or `resolution_height` are 0.
        pub fn new(
            position: RelativisticVector,
            forward_world: na::Vector3<f64>,
            up_world: na::Vector3<f64>,
            focal_length: f64,
            sensor_diagonal: f64,
            resolution_width: u32,
            resolution_height: u32,
        ) -> Camera {

            if position.covariance != Covariance::Contravariant {
                panic!("The camera position vector must be contravariant");
            }
            
            // [Checks on arguments]
            if focal_length <= 0.0 {
                panic!("focal_length must be greater than 0");
            }
            if sensor_diagonal <= 0.0 {
                panic!("sensor_diagonal must be greater than 0");
            }
            if resolution_width == 0 || resolution_height == 0 {
                panic!("resolution_width and resolution_height must be greater than 0");
            }
            
            let orientation = algebra::Orientation::new(forward_world, up_world);
            let camera_to_world_rotation_matrix = orientation.rotation_matrix();

            let aspect_ratio = resolution_width as f64 / resolution_height as f64;
            let aspect_ratio_squared = aspect_ratio.powi(2);
            let sensor_height = (sensor_diagonal.powi(2) / (aspect_ratio_squared + 1.0)).sqrt();
            let sensor_width = aspect_ratio * sensor_height;

            Camera {
                position,
                orientation,
                focal_length,
                sensor_width,
                sensor_height,
                resolution_width,
                resolution_height,
                camera_to_world_rotation_matrix,
            }
        }

        /// Getter method for the position vector
        pub fn position(&self) -> &RelativisticVector {
            &self.position
        }

        /// Getter method for the orientation
        pub fn orientation(&self) -> &algebra::Orientation {
            &self.orientation
        }

        /// Getter method for the focal length
        pub fn update_position(&mut self, new_position: RelativisticVector) {
            if new_position.covariance != Covariance::Contravariant {
                panic!("The camera position vector must be contravariant");
            }
            self.position = new_position;
        }

        /// Getter method for the focal length
        pub fn update_orientation(&mut self, forward_world: na::Vector3<f64>, up_world: na::Vector3<f64>) {
            self.orientation = algebra::Orientation::new(forward_world, up_world);
            self.camera_to_world_rotation_matrix = self.orientation.rotation_matrix();
        }

        /// Returns the vector pointing out of the camera in the direction of the image,
        /// identified by the pixel indexes x and y. The vector is in the camera space.
        pub fn outward_vector_on_camera_space(&self, pixel_x: u32, pixel_y: u32) -> na::Vector3<f64> {
            
            let res_x = self.resolution_width as f64;
            let res_y = self.resolution_height as f64;
            
            let h = 0.5 - (pixel_y as f64 / res_y); // height
            let w = (pixel_x as f64 / res_x) - 0.5; // width
            let x = &self.focal_length * 1.0;
            let y = -&self.sensor_width * w;
            let z = &self.sensor_height * h;
            
            let vec = na::Vector3::new(x, y, z);

            vec.normalize()
        }

        /// [This doc is to be checked]
        /// Returns the vector pointing out of the camera in the direction of the image,
        /// identified by the pixel indexes x and y. The vector is in the world space, not on camera space.
        pub fn outward_vector_on_world_space_from_x_y(&self, pixel_x: u32, pixel_y: u32) -> na::Vector3<f64> {
            let vec = self.outward_vector_on_camera_space(pixel_x, pixel_y);
            self.camera_to_world_rotation_matrix.transform_vector(&vec)
        }

        /// Given a vector with coordinates in camera space, returns the corresponding
        /// vector in world space.
        pub fn outward_vector_on_world_space_from_camera_vector(&self, vec: &na::Vector3<f64>) -> na::Vector3<f64> {
            self.camera_to_world_rotation_matrix.transform_vector(vec)
        }

        /// Getter method for resolution_width
        pub fn resolution_width(&self) -> u32 {
            self.resolution_width
        }
        
        /// Getter method for resolution_height
        pub fn resolution_height(&self) -> u32 {
            self.resolution_height
        }

    }

}