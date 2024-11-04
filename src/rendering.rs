/// This module containes useful structs that can be used to quickly setup a `RelativisticSystem``
/// to render images and videos.
/// 
/// See the systems.rs module for more informations on how `RenderingSystem`s work.

pub mod rendering {

    use nalgebra as na;
    use serde::{Serialize, Deserialize};
    use std::path::PathBuf;
    use crate::{
        cameras::cameras::Camera, images::spherical_images::load_image_as_spherical_image, interpolation::interpolation::Interpolator, metrics::metrics::DiagonalSphericalMetric, systems::systems::RelativisticSystem, vectors::vectors::Covariance
    };

    /// Use this struct to quickly setup an object that can render images.
    pub struct ImageRenderingSystem<M: DiagonalSphericalMetric> {
        relativistic_system: RelativisticSystem<M>,
        image_rendering_settings: ImageRenderingSettings
        // pass
    }

    impl<M: DiagonalSphericalMetric> ImageRenderingSystem<M> {

        /// Returns an instance of `ImageRenderingSystem` from a chosen metric and corresponding
        /// settings.
        /// 
        /// Internally, it creates a `RelativisticSystem` that holds a `Camera` and two `SphericalImage`s
        /// that will be used to render the warped space-time.
        /// 
        /// # Arguments
        /// - `metric`: the chosen metric.
        /// - `image_rendering_settings`: the settings for the image rendering.
        pub fn new(metric: M, image_rendering_settings: ImageRenderingSettings) -> Self {
            
            // Loading background images
            let image_1 = load_image_as_spherical_image(
                &image_rendering_settings.path_to_background_image_1, None, None);
            let image_2 = load_image_as_spherical_image(
                &image_rendering_settings.path_to_background_image_2, None, None);

            // Creating camera
            let camera = Camera::new(
                crate::vectors::vectors::RelativisticVector {
                    vector: image_rendering_settings.camera_position,
                    covariance: Covariance::Contravariant
                },
                image_rendering_settings.camera_forward,
                image_rendering_settings.camera_up,
                image_rendering_settings.camera_focal_length,
                image_rendering_settings.camera_diagonal,
                image_rendering_settings.resolution_x,
                image_rendering_settings.resolution_y
            );

            // Creating relativistic system

            let rel_system = RelativisticSystem::new(
                metric,
                image_1,
                image_2,
                camera
            );

            ImageRenderingSystem {
                relativistic_system: rel_system,
                image_rendering_settings
            }

        }

        /// Generates an image given settings passed to the `ImageRenderingSystem` and saves
        /// it to file. See `ImageRenderingSystem` for more information on the settings that
        /// can be specified.
        /// 
        /// Internally, it uses a `RelativisticSystem` to render the image.
        /// 
        /// # Errors
        /// Issues with the rendering process are generally returned as a `String` error describing
        /// the problem.
        /// 
        /// # Panics
        /// Some functions called by this method may panic.
        /// 
        /// See the documentation of inner functions for more details about errors and panics.
        pub fn render(&self) -> Result<(), String> {

            // Creating folder for output image

            let output_folder = &self.image_rendering_settings.path_to_output_folder;
            if !output_folder.exists() {
                match std::fs::create_dir(output_folder) {
                    Ok(()) => {},
                    Err(err) => return Err(format!("Could not create video output folder {output_folder:?} due to error: {err}"))
                }
            }

            let image = self.relativistic_system.render_image_efficient(
                self.image_rendering_settings.max_iterations_propagation,
                self.image_rendering_settings.escape_radius,
                self.image_rendering_settings.ray_integration_step,
                self.image_rendering_settings.alphas_num,
                self.image_rendering_settings.max_iterations_sampling,
                self.image_rendering_settings.sampling_convergence_threshold_1,
                self.image_rendering_settings.sampling_convergence_threshold_2,
                // &self.image_rendering_settings.path_to_output_image
            );

            let path_of_image = output_folder.join(&self.image_rendering_settings.output_image_name).with_extension("png");

            match image.save(&path_of_image) {
                Ok(()) => {}
                Err(err) => {return Err(format!("Could not save image frame {path_of_image:?} due to error: {err}"))}
            }

            Ok(())
            
        }
    }

    /// Contains all the informations required to render an image of a warped space-time.
    /// 
    /// # Fields to specify image characteristics
    /// - `resolution_x`: the horizontal resolution of the image.
    /// - `resolution_y`: the vertical resolution of the image.
    /// - `camera_diagonal`: the diagonal of the camera sensor.
    /// - `camera_focal_length`: the focal length of the camera.
    /// - `camera_position`: the position of the camera with metric coordinates (4D vector: t, l, theta, phi).
    /// - `camera_forward`: the (x, y, z) vector specifiying where the camera is pointed at, with tangent space coordinates.
    ///     If you were the camera, this is where your eyes are looking at.
    /// - `camera_up`: the (x, y, z) vector specifiying the direction of the camera's up, with tangent space coordinates.
    ///     If you were the camera, this is the direction of the top of your head.
    /// 
    /// # Fields to specify input and output image files
    /// - `path_to_background_image_1`: the path to the image that will be used as the background of the warped space-time, on the positive-l side.
    /// - `path_to_background_image_2`: the path to the image that will be used as the background of the warped space-time, on the negative-l side.
    /// - `path_to_output_folder`: the path to the folder where the output image will be saved (must exist).
    /// - `output_image_name`: the name of the output image (wihtout extension).
    /// 
    /// # Fields to control the rendering process
    /// **REFERENCE TO UPDATE** See also the documentation on the rendering and sampling algorithms.
    /// - `escape_radius`: the (absolute) l coordinate past which the ray is considered to have escaped the metric.
    /// - `max_iterations_propagation`: the maximum number of iterations to propagate light rays.
    /// - `ray_integration_step`: the step size (proper time) to integrate light rays.
    /// - `alphas_num`: the initial number of alphas for the sampling algorithm.
    /// - `max_iterations_sampling`: the maximum number of iterations for the sampling algorithm.
    /// - `sampling_convergence_threshold_1`: the threshold for the angular convergence criterion of the sampling algorithm.
    /// - `sampling_convergence_threshold_2`: the threshold for the space-sign convergence criterion of the sampling algorithm.
    pub struct ImageRenderingSettings {
        pub resolution_x: u32,
        pub resolution_y: u32,
        pub camera_diagonal: f64,
        pub camera_focal_length: f64,
        pub camera_position: na::Vector4<f64>,
        pub camera_forward: na::Vector3<f64>,
        pub camera_up: na::Vector3<f64>,
        pub path_to_background_image_1: PathBuf,
        pub path_to_background_image_2: PathBuf,
        pub path_to_output_folder: PathBuf,
        pub output_image_name: String,
        pub escape_radius: f64,
        pub max_iterations_propagation: u32,
        pub ray_integration_step: f64,
        pub alphas_num: u32,
        pub max_iterations_sampling: u32,
        pub sampling_convergence_threshold_1: f64,
        pub sampling_convergence_threshold_2: f64,
    }

    /// Use this struct to quickly setup an object that can render videos.
    pub struct VideoRenderingSystem<M: DiagonalSphericalMetric> {
        relativistic_system: RelativisticSystem<M>,
        interpolator: Interpolator,
        video_rendering_settings: VideoRenderingSettings
    }

    impl <M: DiagonalSphericalMetric> VideoRenderingSystem<M> {

        /// Returns an instance of `VideoRenderingSystem` from a chosen metric and corresponding
        /// settings.
        /// 
        /// Internally, it creates a `RelativisticSystem` that holds a `Camera` and two `SphericalImage`s
        /// that will be used to render the warped space-time. It also creates an `Interpolator` object
        /// that will be used to interpolate the camera path over time.
        /// 
        /// # Arguments
        /// - `metric`: the chosen metric.
        /// - `video_rendering_settings`: the settings for the video rendering.
        pub fn new(metric: M, video_rendering_settings: VideoRenderingSettings) -> Self {
            let image_1 = load_image_as_spherical_image(
                &video_rendering_settings.filepath_to_background_image_1, None, None);
            let image_2 = load_image_as_spherical_image(
                &video_rendering_settings.filepath_to_background_image_2, None, None);
            
            let interpolator = Interpolator::from_file( &video_rendering_settings.filepath_to_camera_path);

            let camera = Camera::new(
                crate::vectors::vectors::RelativisticVector {
                    vector: interpolator.camera_position(interpolator.min_time()),
                    covariance: Covariance::Contravariant
                },
                interpolator.camera_forward(interpolator.min_time()),
                interpolator.camera_up(interpolator.min_time()),
                video_rendering_settings.camera_focal_length,
                video_rendering_settings.camera_diagonal,
                video_rendering_settings.resolution_x,
                video_rendering_settings.resolution_y
            );

            let rel_system = RelativisticSystem::new(
                metric,
                image_1,
                image_2,
                camera
            );

            VideoRenderingSystem {
                relativistic_system: rel_system,
                interpolator,
                video_rendering_settings
            }
        }

        /// Returns the time coordinates of all the frames that are to be rendered.
        fn times_of_frames(&self) -> Vec<f64> {
            let min_time = self.interpolator.min_time();
            let max_time = self.interpolator.max_time();
            let delta_time = 1.0/self.video_rendering_settings.frame_rate;

            let mut times = Vec::new();
            let mut t = min_time;

            while t < max_time {
                times.push(t);
                t += delta_time;
            }

            times
        }

        /// Given a time t, updates the camera position and orientation to the values
        /// interpolated from the camera path.
        fn update_camera(&mut self, t: f64) {
            
            self.relativistic_system.camera.update_position(crate::vectors::vectors::RelativisticVector {
                vector: self.interpolator.camera_position(t),
                covariance: Covariance::Contravariant
            });
            
            self.relativistic_system.camera.update_orientation(
                self.interpolator.camera_forward(t),
                self.interpolator.camera_up(t)
            );
        }

        /// Generates a video given settings passed to the `VideoRenderingSystem`. Currently,
        /// it saves all the frames to a temporary folder and does not actually create a single
        /// video file.
        pub fn render(&mut self) -> Result<(), String> {

            let times = self.times_of_frames();

            // Creating tmp folder in the same folder of path_to_output_video
        

            // Creating folder for output video

            let output_folder = &self.video_rendering_settings.filepath_to_output_folder;
            if !output_folder.exists() {
                match std::fs::create_dir(output_folder) {
                    Ok(()) => {},
                    Err(err) => return Err(format!("Could not create video output folder {output_folder:?} due to error: {err}"))
                }
            }

            // Creating tmp folder
            let tmp_folder = output_folder.join("tmp");
            if tmp_folder.exists() {
                match std::fs::remove_dir_all(&tmp_folder) {
                    Ok(()) => {},
                    Err(err) => return Err(format!("Could not remove pre-existing tmp folder {tmp_folder:?} due to error: {err}"))
                }
            }

            match std::fs::create_dir(&tmp_folder) {
                Ok(()) => {},
                Err(err) => return Err(format!("Could not create tmp output folder {tmp_folder:?} due to error: {err}"))
            }

            println!("Rendering {} frames...", times.len());

            for (index, time) in times.iter().enumerate() {

                // Printing progress
                println!("Rendering frame {}/{}...", index+1, times.len());

                let path_to_output_video_frame = tmp_folder.join(format!("frame_{}.png", index));

                self.update_camera(*time);
                let image_frame = self.relativistic_system.render_image_efficient(
                    self.video_rendering_settings.max_iterations_propagation,
                    self.video_rendering_settings.escape_radius,
                    self.video_rendering_settings.ray_integration_step,
                    self.video_rendering_settings.alphas_num,
                    self.video_rendering_settings.max_iterations_sampling,
                    self.video_rendering_settings.sampling_convergence_threshold_1,
                    self.video_rendering_settings.sampling_convergence_threshold_1,
                );
                

                // Saving image to file
                match image_frame.save(&path_to_output_video_frame) {
                    Ok(()) => {}
                    Err(err) => {return Err(format!("Could not save image frame {path_to_output_video_frame:?} due to error: {err}"))}
                }

            }

            // Creating video from frames
            // [ToDo]
            
            // Removing tmp folder and its content

            // std::fs::remove_dir_all(&tmp_folder).expect("Could not remove tmp folder");

            Ok(())

        }
    }

    /// Contains all the informations required to render a video of a warped space-time.
    /// 
    /// # Fields to specify video characteristics
    /// - `frame_rate`: the frame rate of the video.
    /// - `resolution_x`: the horizontal resolution of the video frames.
    /// - `resolution_y`: the vertical resolution of the video frames.
    /// - `camera_diagonal`: the diagonal of the camera sensor.
    /// - `camera_focal_length`: the focal length of the camera.
    /// 
    /// # Fields to specify input and output video files
    /// - `path_to_camera_path`: must be a csv file that contains the camera path over time. See the `Interpolator` struct for more informations.
    /// - `path_to_background_image_1`: the path to the image that will be used as the background of the warped space-time, on the positive-l side.
    /// - `path_to_background_image_2`: the path to the image that will be used as the background of the warped space-time, on the negative-l side.
    /// - `path_to_output_folder`: the path to the folder where the output video will be saved (must exist **CHECK**).
    /// - `output_video_name`: the name of the output video (wihtout extension).
    /// 
    /// # Fields to control the rendering process
    /// **REFERENCE TO UPDATE** See also the documentation on the rendering and sampling algorithms.
    /// - `escape_radius`: the (absolute) l coordinate past which the ray is considered to have escaped the metric.
    /// - `max_iterations_propagation`: the maximum number of iterations to propagate light rays.
    /// - `ray_integration_step`: the step size (proper time) to integrate light rays.
    /// - `alphas_num`: the initial number of alphas for the sampling algorithm.
    /// - `max_iterations_sampling`: the maximum number of iterations for the sampling algorithm.
    /// - `sampling_convergence_threshold_1`: the threshold for the angular convergence criterion of the sampling algorithm.
    /// - `sampling_convergence_threshold_2`: the threshold for the space-sign convergence criterion of the sampling algorithm.
    #[derive(Serialize, Deserialize)]
    pub struct VideoRenderingSettings {
        pub frame_rate: f64,
        pub resolution_x: u32,
        pub resolution_y: u32,
        pub camera_diagonal: f64,
        pub camera_focal_length: f64,
        pub filepath_to_camera_path: PathBuf,
        pub filepath_to_background_image_1: PathBuf,
        pub filepath_to_background_image_2: PathBuf,
        pub filepath_to_output_folder: PathBuf,
        pub output_video_name: String,
        pub escape_radius: f64,
        pub max_iterations_propagation: u32,
        pub ray_integration_step: f64,
        pub alphas_num: u32,
        pub max_iterations_sampling: u32,
        pub sampling_convergence_threshold_1: f64,
        pub sampling_convergence_threshold_2: f64,
    }
    
}