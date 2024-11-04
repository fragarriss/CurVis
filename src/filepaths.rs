// Here we store file paths to various files used by curvis.exe

pub mod filepaths {

    use std::env;
    use std::path::PathBuf;

    // These paths are specified with respect to the package folder (the parent of src)
    const VIDEO_SETTINGS_DEFAULT: &str = "settings/defaults/video_settings.toml";
    const IMAGE_SETTINGS_DEFAULT: &str = "settings/defaults/image_settings.toml";
    const CAMERA_SETTINGS_DEFAULT: &str = "settings/defaults/camera_settings.toml";
    const ELLIS_METRIC_SETTINGS_DEFAULT: &str = "settings/defaults/ellis_metric_settings.toml";
    const INTERSTELLAR_METRIC_SETTINGS_DEFAULT: &str = "settings/defaults/interstellar_metric_settings.toml";
    const SIMULATION_SETTINGS_DEFAULT: &str = "settings/defaults/simulation_settings.toml";

    // The functions here are intentionally panicking when the paths are not
    // found, as the program cannot run if these are not setup correctly.

    /// Returns the full path to the curvis.exe file
    /// Should be either curvis/target/release/curvis.exe or curvis/target/debug/curvis.exe
    /// N.B. For security reason, paths obtained from this function, should NOT be used
    /// to call other executables. Here we use it to:
    ///  - Load default settings files
    fn get_exe_path() -> PathBuf {
        env::current_exe().unwrap()
    }

    /// Returns the full path to the curvis package folder
    pub fn get_package_path() -> PathBuf {
        
        let exe_path = get_exe_path();
        let pck_path = PathBuf::from(exe_path
            .parent().unwrap()
            .parent().unwrap()
            .parent().unwrap());

        pck_path
    }

    /// Joins the curvis package folder path with `file_path` if it is relative.
    /// If `file_path` is absolute, it is returned as is.
    pub fn resolve_path(file_path: &PathBuf) -> PathBuf {
        
        // Replaces with file_path if file_path is absolute
        return get_package_path().join(file_path)
        
    }

    /// Returns true if `file_path` has an extension that matches `extension`, false otherwise.
    /// `extension` is specified without leading point.
    pub fn check_extension(file_path: &PathBuf, extension: &str) -> bool {
        match file_path.extension() {
            None => { return false }
            Some(ext) => {
                if ext != extension {
                    return false ;
                };
            }
        };
        return true
    }

    /// Returns the full path to a setting file. It is used with the relative
    /// paths specified at the top of the `filepath` module.
    fn get_settings_file_path(relative_to_package: &PathBuf) -> Result<PathBuf, String> {

        let full_path = resolve_path(relative_to_package);

        match ensure_is_toml_and_exists(&full_path) {
            Ok(_) => {return Ok(full_path)}
            Err(string) => {return Err(string)}
        }
    }

    /// Checks whether the file_path corresponds to an existing toml file.
    pub fn ensure_is_toml_and_exists(file_path: &PathBuf) -> Result<&PathBuf, String> {

        if !file_path.exists() {
            return Result::Err(format!("The path {:?} does not exist", file_path));
        };

        if !file_path.is_file() {
            return Result::Err(format!("The path {:?} does not correspond to a file", file_path));
        };

        if !check_extension(file_path, "toml") {
            return Result::Err(format!("The path {:?} is not a toml file", file_path));
        }

        Ok(file_path)

    }

    /// Return the full path to the default video settings file
    pub fn get_video_settings_default_path() -> PathBuf {
        return get_settings_file_path(&PathBuf::from(VIDEO_SETTINGS_DEFAULT)).unwrap()
    }

    /// Return the full path to the default image settings file
    pub fn get_image_settings_default_path() -> PathBuf {
        return get_settings_file_path(&PathBuf::from(IMAGE_SETTINGS_DEFAULT)).unwrap()
    }

    /// Return the full path to the default camera settings file
    pub fn get_camera_settings_default_path() -> PathBuf {
        return get_settings_file_path(&PathBuf::from(CAMERA_SETTINGS_DEFAULT)).unwrap()
    }

    /// Return the full path to the default Ellis metric settings file
    pub fn get_ellis_metric_settings_default_path() -> PathBuf {
        return get_settings_file_path(&PathBuf::from(ELLIS_METRIC_SETTINGS_DEFAULT)).unwrap()
    }

    /// Return the full path to the default Interstellar metric setting file
    pub fn get_interstellar_metric_settings_default_path() -> PathBuf {
        return get_settings_file_path(&PathBuf::from(INTERSTELLAR_METRIC_SETTINGS_DEFAULT)).unwrap()
    }

    /// Return the full path to the default simulation settings file
    pub fn get_simulation_settings_default_path() -> PathBuf {
        return get_settings_file_path(&PathBuf::from(SIMULATION_SETTINGS_DEFAULT)).unwrap()
    }

}