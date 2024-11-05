// Here we define the various structs containing the settings to be used
// for the rendering

pub mod settings {

    use crate::filepaths::filepaths::{self, check_extension, ensure_is_toml_and_exists, resolve_path};
    use duplicate::duplicate_item;
    use serde::Deserialize;
    use toml;
    use std::fs;
    use std::default::Default;
    use std::path::PathBuf;
    
    pub trait Validate {
        fn validate(&self) -> Result<(), String>;
    }

    pub trait Normalize {
        fn normalize(&mut self) {}
    }

    #[derive(Deserialize)]
    pub struct VideoSettings {
        pub video_name: String,
        pub frame_rate: f64,
        pub filepath_to_camera_path: PathBuf,
    }

    impl Normalize for VideoSettings {
        fn normalize(&mut self) {
            self.filepath_to_camera_path = resolve_path(&self.filepath_to_camera_path);
        }
    }

    impl Validate for VideoSettings {
        
        fn validate(&self) -> Result<(), String> {
            
            if self.video_name.is_empty() {
                return Err(format!("Video name cannot be an empty string."));
            }

            // Check that the camera path exists and is a .csv file

            if !check_extension(&self.filepath_to_camera_path, 
            "csv") {
                return Err(format!("The camera path {:?} is not a csv file.", self.filepath_to_camera_path));
            }

            if !std::path::Path::new(&self.filepath_to_camera_path).exists() {
                return Err(format!("The camera path {:?} does not exist.", self.filepath_to_camera_path));
            }

            Ok(())
        }
    }

    #[derive(Deserialize)]
    pub struct ImageSettings {
        pub image_name: String,
        pub t: f64,
        pub l: f64,
        pub theta: f64,
        pub phi: f64,
        pub forward_x: f64,
        pub forward_y: f64,
        pub forward_z: f64,
        pub up_x: f64,
        pub up_y: f64,
        pub up_z: f64,
    }

    impl Normalize for ImageSettings {}
    impl Validate for ImageSettings {
        fn validate(&self) -> Result<(), String> { 
            
            if self.image_name.is_empty() {
                return Err(format!("Image name cannot be an empty string."));
            }
            Ok(()) }
    }

    #[derive(Deserialize)]
    pub struct CameraSettings {
        pub resolution_x: u32,
        pub resolution_y: u32,
        pub diagonal: f64,
        pub focal_length: f64,
    }

    impl Normalize for CameraSettings {}
    impl Validate for CameraSettings {
        fn validate(&self) -> Result<(), String> { 

            // Check that all the values are positive

            if self.resolution_x <= 0 {
                return Err(String::from("The resolution in the x direction must be larger than zero."));
            }

            if self.resolution_y <= 0 {
                return Err(String::from("The resolution in the y direction must be larger than zero."));
            }

            if self.diagonal <= 0.0 {
                return Err(String::from("The diagonal of the camera must be larger than zero."));
            }

            if self.focal_length <= 0.0 {
                return Err(String::from("The focal length of the camera must be larger than zero."));
            }

            Ok(())

        }
    }

    #[derive(Deserialize)]
    pub struct SimulationSettings {
        pub escape_radius: f64,
        pub ray_integration_max_itarations: u32,
        pub ray_integration_step: f64,
        pub sampling_initial_nums: u32,
        pub sampling_max_iterations: u32,
        pub sampling_convergence_threshold_1: f64,
        pub sampling_convergence_threshold_2: f64,
    }

    impl Normalize for SimulationSettings {}
    impl Validate for SimulationSettings {
        fn validate(&self) -> Result<(), String> { 

            // Check that all the values are positive (apart from sampling_initial_nums, that must be larger than 1)

            if self.escape_radius <= 0.0 {
                return Err(String::from("The escape radius must be larger than zero."));
            }

            if self.ray_integration_max_itarations <= 0 {
                return Err(String::from("The maximum number of iterations for the ray integration must be larger than zero."));
            }

            if self.ray_integration_step <= 0.0 {
                return Err(String::from("The step for the ray integration must be larger than zero."));
            }

            if self.sampling_initial_nums <= 1 {
                return Err(String::from("The initial number of samples must be larger than two."));
            }

            if self.sampling_max_iterations <= 0 {
                return Err(String::from("The maximum number of iterations for the sampling must be larger than zero."));
            }

            if self.sampling_convergence_threshold_1 <= 0.0 {
                return Err(String::from("The first convergence threshold for the sampling must be larger than zero."));
            }

            if self.sampling_convergence_threshold_2 <= 0.0 {
                return Err(String::from("The second convergence threshold for the sampling must be larger than zero."));
            }

            Ok(())

        }
    }

    #[derive(Deserialize)]
    pub struct EllisMetricSettings {
        pub rho: f64
    }

    impl Normalize for EllisMetricSettings {}
    impl Validate for EllisMetricSettings {
        fn validate(&self) -> Result<(), String> { 

            // Check that rho is positive

            if self.rho <= 0.0 {
                return Err(String::from("The density parameter rho must be larger than zero."));
            }

            Ok(())

        }
    }

    #[derive(Deserialize)]
    pub struct InterstellarMetricSettings {
        pub m: f64,
        pub a: f64,
        pub rho: f64,
    }

    impl Normalize for InterstellarMetricSettings {}
    impl Validate for InterstellarMetricSettings {
        
        fn validate(&self) -> Result<(), String> { 

            // Check that all the values are positive

            if self.m <= 0.0 {
                return Err(String::from("The mass parameter m must be larger than zero."));
            }

            if self.a <= 0.0 {
                return Err(String::from("The spin parameter a must be larger than zero."));
            }

            if self.rho <= 0.0 {
                return Err(String::from("The density parameter rho must be larger than zero."));
            }

            Ok(())

        }
    }
    
    
    pub trait FromToml {
        fn from_toml_file(toml_file_path: &PathBuf) -> Result<Self, String> where Self:Sized;
    }

    #[duplicate_item(name;
        [VideoSettings];
        [ImageSettings];
        [CameraSettings];
        [SimulationSettings];
        [EllisMetricSettings];
        [InterstellarMetricSettings];
    )]
    impl FromToml for name {
        
        fn from_toml_file(toml_file_path: &PathBuf) -> Result<name, String> {
            
            // Checking if the file exists and if it is a .toml file

            ensure_is_toml_and_exists(toml_file_path)?;

            let string = match fs::read_to_string(toml_file_path) {
                Ok(string) => { string }
                Err(error) => { return Err(format!("Could not read file {toml_file_path:?} due to error: {error:?}")) }
            };

            let read = toml::from_str(&string);

            match read
            {
                Result::Ok(settings) => { Ok(settings) }
                Result::Err(error) => { Err(String::from(error.message())) }     
            }
        }
    }

    // Implementation of Default trait for all settings
    // Settings are loaded from default files in the curvis/settings/defaults folder

    // Default
    impl Default for VideoSettings {
        fn default() -> Self {
            VideoSettings::from_toml_file( &filepaths::get_video_settings_default_path()).unwrap()
        }
    }

    impl Default for ImageSettings {
        fn default() -> Self {
            ImageSettings::from_toml_file( &filepaths::get_image_settings_default_path()).unwrap()
        }
    }

    impl Default for EllisMetricSettings {
        fn default() -> Self {
            EllisMetricSettings::from_toml_file( &filepaths::get_ellis_metric_settings_default_path()).unwrap()
        }
    }

    impl Default for InterstellarMetricSettings {
        fn default() -> Self {
            InterstellarMetricSettings::from_toml_file( &filepaths::get_interstellar_metric_settings_default_path()).unwrap()
        }
    }

    impl Default for CameraSettings {
        fn default() -> Self {
            CameraSettings::from_toml_file( &filepaths::get_camera_settings_default_path()).unwrap()
        }
    }

    impl Default for SimulationSettings {
        fn default() -> Self {
            SimulationSettings::from_toml_file( &filepaths::get_simulation_settings_default_path()).unwrap()
        }
    }
}