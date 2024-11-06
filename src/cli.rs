/// Command line interface module

/// Curvis works with either of three command patterns:
/// 1) curvis image [arguments]
/// 2) curvis video [arguments]
/// 3) curvis custom
/// 
/// Arguments specific to image:
/// --image-settings
/// 
/// Arguments specific to video:
/// --video-settings
/// 
/// Common arguments
/// --background-image-1
/// --background-image-2
/// --metric-settings
/// --camera-settings
/// --simulation-settings
/// 
/// "background-image-1" and "background-image-2" are required arguments, and point to the images
/// to be used as background for the rendering.
/// 
/// The other arguments are optional, and point to toml files which contain the various settings
/// for the rendering. If they are not specified, default values are used.

pub mod cli {
    use clap::{command, Arg, Command, ArgMatches, value_parser};
    use either::Either;
    use crate::settings::settings::{
        CameraSettings, EllisMetricSettings, FromToml, ImageSettings, InterstellarMetricSettings, SimulationSettings, VideoSettings};
    use std::process::exit;
    use std::path::PathBuf;

    fn get_cli() -> Command {
        
        // Defining arguments

        let background1 = Arg::new("background_image_1")
            .value_parser(value_parser!(PathBuf))
            .required(true)
            .value_name("IMAGE FILE 1");

        let background2 = Arg::new("background_image_2")
            .value_parser(value_parser!(PathBuf))
            .required(true)
            .value_name("IMAGE FILE 2");

        let output_folder = Arg::new("output_folder")
            .value_parser(value_parser!(PathBuf))
            .required(false)
            .value_name("OUTPUT FOLDER")
            .default_value(None);

        let video = Arg::new("video_settings")
            .long("video-settings")
            .short('v')
            .required(false)
            // .value_parser(value_parser!(PathBuf))
            .value_name("TOML FILE"); // We require a file path

        let image = Arg::new("image_settings")
            .long("image-settings")
            .short('i')
            .required(false)
            // .value_parser(value_parser!(PathBuf))
            .value_name("TOML FILE"); // We require a file path

        let metric = Arg::new("metric_settings")
            .long("metric-settings")
            .short('m')
            .required(false)
            // .value_parser(value_parser!(PathBuf))
            .value_name("TOML FILE"); // We require a file path

        let camera = Arg::new("camera_settings")
            .long("camera-settings")
            .short('c')
            .required(false)
            // .value_parser(value_parser!(PathBuf))
            .value_name("TOML FILE"); // We require a file path

        let simulation = Arg::new("simulation_settings")
            .long("simulation-settings")
            .short('s')
            .required(false)
            // .value_parser(value_parser!(PathBuf))
            .value_name("TOML FILE"); // We require a file path

        // Defining command with its subcommands

        let command = command!("curvis")
            .subcommand(
                Command::new("image")
                    .about("renders a single image frame")
                    .arg(background1.clone())
                    .arg(background2.clone())
                    .arg(output_folder.clone())
                    .arg(image)
                    .arg(metric.clone())
                    .arg(camera.clone())
                    .arg(simulation.clone())
                )
            .subcommand(
                    Command::new("video")
                    .about("renders a video")
                    .arg(background1.clone())
                    .arg(background2.clone())
                    .arg(output_folder.clone())
                    .arg(video)
                    .arg(metric.clone())
                    .arg(camera.clone())
                    .arg(simulation.clone())
            )
            .subcommand(
                Command::new("custom")
                .about("runs the custom script")
            );

        command

    }

    fn get_global_matches(command: &mut Command) -> ArgMatches {
        command.get_matches_mut()
    }

    pub fn get_subcommand() -> String {
        let global_matches = get_global_matches(&mut get_cli());
        let sc = global_matches.subcommand();
        match sc {
            Some((sc, _)) => { return String::from(sc); }
            None => {
                eprintln!("Subcommand not found");
                exit(1)
            }
        }
    }

    // Extracting settings from the matches

    fn get_command_specific_matches(global_matches: &ArgMatches) -> (String, &ArgMatches) {
        match global_matches.subcommand() {
            Some(("image", specific_matches)) => { return (String::from("image"), specific_matches) } // image was used
            Some(("video", specific_matches)) => { return (String::from("video"), specific_matches) } // video was used
            Some(("custom", specific_matches)) => { return (String::from("custom"), specific_matches) } // custom was used
            _ => panic!("Unrecognized subcommand.")
        }
    }

    fn check_filebuf_exists(file_path: PathBuf) -> Result<PathBuf, String> {
        if !(file_path.exists()) {
            Err(format!("File {file_path:?} not found."))
        }
        else { Ok(file_path) }
    }
    
    fn check_filebuf_ref_exists(file_path: &PathBuf) -> Result<&PathBuf, String> {
        if !(file_path.exists()) {
            Err(format!("File {file_path:?} not found."))
        }
        else { Ok(file_path) }
    }

    fn check_filebuf_ref_is_folder(file_path: &PathBuf) -> Result<&PathBuf, String> {
        if !(file_path.is_dir()) {
            Err(format!("{file_path:?} is not a folder."))
        }
        else { Ok(file_path) }
    }

    fn filepath_to_bg_img_1_from_matches(arg_matches: &ArgMatches) -> Result<PathBuf, String> {
        return match arg_matches.get_one::<PathBuf>("background_image_1") {
            None => { Err(String::from("Missing file path to first background image")) }
            Some(file_path) => { 
                let file_path = PathBuf::from(file_path);
                check_filebuf_exists(file_path)
            }
        }
    }
    
    fn filepath_to_bg_img_2_from_matches(arg_matches: &ArgMatches) -> Result<PathBuf, String> {
        return match arg_matches.get_one::<PathBuf>("background_image_2") {
            None => { Err(String::from("Missing file path to first background image")) }
            Some(file_path) => {
                let file_path = PathBuf::from(file_path);
                check_filebuf_exists(file_path)
            }
        }
    }

    fn output_folder_path_from_matches(arg_matches: &ArgMatches) -> Result<PathBuf, String> {
        return match arg_matches.get_one::<String>("output_folder") {
            None => { 
                // Using current working directory
                let curr_dir = std::env::current_dir();
                match curr_dir {
                    Ok(curr_dir) => { Ok(curr_dir) }
                    Err(err) => { Err(format!("Could not get current working directory: {err}.")) }
                }
             }
            Some(file_path) => {
                let file_path = &PathBuf::from(file_path);
                check_filebuf_ref_exists(file_path)?;
                check_filebuf_ref_is_folder(file_path)?;
                Ok(PathBuf::from(file_path))
            }
        }
    }
    
    fn video_settings_from_matches(arg_matches: &ArgMatches) -> Result<VideoSettings, String> {
        return match arg_matches.get_one::<String>("video_settings") {
            None => { Ok(VideoSettings::default()) }
            Some(file_to_settings) => {
                VideoSettings::from_toml_file(
                    check_filebuf_ref_exists(&PathBuf::from(file_to_settings))?
                )
            }
        }
    }

    fn image_settings_from_matches(arg_matches: &ArgMatches) -> Result<ImageSettings, String> {
        return match arg_matches.get_one::<String>("image_settings") {
            None => { Ok(ImageSettings::default()) }
            Some(file_to_settings) => {
                ImageSettings::from_toml_file(
                    check_filebuf_ref_exists(&PathBuf::from(file_to_settings))?
                )
            }
        }
    }
    
    fn metric_settings_from_matches (arg_matches: &ArgMatches)
        -> Result<Either<EllisMetricSettings, InterstellarMetricSettings>, String> {
        
        return match arg_matches.get_one::<String>("metric_settings") {
            
            None => { Ok(Either::Left(EllisMetricSettings::default())) }

            Some(file_to_settings) => {

                // First I try reading the file as an ellis metric settings
                // If it fails, I try again with the Interstellar metric

                let file_path = PathBuf::from(file_to_settings);
                check_filebuf_ref_exists(&file_path)?;

                if let Ok(settings) = InterstellarMetricSettings::from_toml_file(&file_path) {
                    Ok(Either::Right(settings))
                }
                else if let Ok(settings) = EllisMetricSettings::from_toml_file(&file_path) {
                    Ok(Either::Left(settings))
                }
                else {
                    Err(String::from("Could not read the metric configuration file."))
                }

            }
        }

    }

    fn camera_settings_from_matches (arg_matches: &ArgMatches) -> Result<CameraSettings, String> {
        return match arg_matches.get_one::<String>("camera_settings") {
            None => { Ok(CameraSettings::default()) }
            Some(file_to_settings) => {
                CameraSettings::from_toml_file(
                check_filebuf_ref_exists(&PathBuf::from(file_to_settings))?)
            }
        }
    }

    fn simulation_settings_from_matches (arg_matches: &ArgMatches) -> Result<SimulationSettings, String> {
        return match arg_matches.get_one::<String>("simulation_settings") {
            None => { Ok(SimulationSettings::default()) }
            Some(file_to_settings) => {
                SimulationSettings::from_toml_file(
                check_filebuf_ref_exists(&PathBuf::from(file_to_settings))?)
            }
        }
    }

    /// Parses the command line arguments and returns the paths to the background images and the settings
    /// to be used for video rendering.
    pub fn parse_video_arguments() -> (
        PathBuf,
        PathBuf,
        PathBuf,
        VideoSettings,
        Either<EllisMetricSettings, InterstellarMetricSettings>,
        CameraSettings,
        SimulationSettings,
    ) {
        let mut command = get_cli();
        let global_matches = get_global_matches(&mut command);
        let (_, matches) = get_command_specific_matches(&global_matches);

        let bg_img_1 = filepath_to_bg_img_1_from_matches(&matches)
            .unwrap_or_else(|string| {
                eprintln!("Error with background image 1: {}", string);
                exit(1);
            });

        let bg_img_2 = filepath_to_bg_img_2_from_matches(&matches)
            .unwrap_or_else(|string| {
                eprintln!("Error with background image 2: {}", string);
                exit(1);
            });

        let output_folder = output_folder_path_from_matches(&matches)
            .unwrap_or_else(|string| {
                eprintln!("Error with output folder: {}", string);
                exit(1);
            });

        let video_settings = video_settings_from_matches(&matches)
            .unwrap_or_else(|string| {
                eprintln!("Error with video settings: {}", string);
                exit(1);
            });
        let metric_settings = metric_settings_from_matches(&matches)
            .unwrap_or_else(|string| {
                eprintln!("Error with metric settings: {}", string);
                exit(1);
            });
        let camera_settings = camera_settings_from_matches(&matches)
            .unwrap_or_else(|string| {
                eprintln!("Error with camera settings: {}", string);
                exit(1);
            });
        let simulation_settings = simulation_settings_from_matches(&matches)
            .unwrap_or_else(|string| {
                eprintln!("Error with simulation settings: {}", string);
                exit(1);
            });

        return (bg_img_1, bg_img_2, output_folder, video_settings, metric_settings, camera_settings, simulation_settings);
    }

    /// Parses the command line arguments and returns the paths to the background images and the settings
    /// to be used for image rendering.
    pub fn parse_image_arguments() -> (
        PathBuf,
        PathBuf,
        PathBuf,
        ImageSettings,
        Either<EllisMetricSettings, InterstellarMetricSettings>,
        CameraSettings,
        SimulationSettings,
    ){
        let mut command = get_cli();
        let global_matches = get_global_matches(&mut command);
        let (_, matches) = get_command_specific_matches(&global_matches);

        let bg_img_1 = filepath_to_bg_img_1_from_matches(&matches)
            .unwrap_or_else(|string| {
                eprintln!("Error with background image 1: {}", string);
                exit(1);
            });

        let bg_img_2 = filepath_to_bg_img_2_from_matches(&matches)  
            .unwrap_or_else(|string| {
                eprintln!("Error with background image 2: {}", string);
                exit(1);
            });

        let output_folder = output_folder_path_from_matches(&matches)
            .unwrap_or_else(|string| {
                eprintln!("Error with output folder: {}", string);
                exit(1);
            });

        let image_settings = image_settings_from_matches(&matches)
            .unwrap_or_else(|string| {
                eprintln!("Error with video or image settings: {}", string);
                exit(1);
            });
        let metric_settings = metric_settings_from_matches(&matches)
            .unwrap_or_else(|string| {
                eprintln!("Error with metric settings: {}", string);
                exit(1);
            });
        let camera_settings = camera_settings_from_matches(&matches)
            .unwrap_or_else(|string| {
                eprintln!("Error with camera settings: {}", string);
                exit(1);
            });
        let simulation_settings = simulation_settings_from_matches(&matches)
            .unwrap_or_else(|string| {
                eprintln!("Error with simulation settings: {}", string);
                exit(1);
            });

        return (bg_img_1, bg_img_2, output_folder, image_settings, metric_settings, camera_settings, simulation_settings);
    }

    pub fn parse_custom_arguments() {}

}