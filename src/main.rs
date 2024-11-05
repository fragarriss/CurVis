use std::env;
use std::process::exit;
use std::path::PathBuf;

use either::Either;
use curvis::cli::cli::{get_subcommand, parse_image_arguments, parse_video_arguments};
use curvis::custom::custom::custom_main;
use curvis::metrics::metrics::{EllisMetric, InterstellarMetric, Metric};
use curvis::settings::settings::{CameraSettings, EllisMetricSettings, ImageSettings, InterstellarMetricSettings, Normalize, SimulationSettings, Validate, VideoSettings};
use curvis::rendering::rendering::{ImageRenderingSystem, VideoRenderingSystem, ImageRenderingSettings, VideoRenderingSettings};
use nalgebra as na;

/// Returns a single VideoRenderingSettings instance from the passed arguments.
fn setup_video_rendering_settings(
    filepath_bg_1: PathBuf,
    filepath_bg_2: PathBuf,
    output_folder: PathBuf,
    video_settings: &mut VideoSettings,
    camera_settings: &mut CameraSettings,
    simulation_settings: &mut SimulationSettings,
) -> Result<VideoRenderingSettings, String> {

    video_settings.normalize();
    video_settings.validate()?;

    camera_settings.normalize();
    camera_settings.validate()?;

    simulation_settings.normalize();
    simulation_settings.validate()?;

    Ok(VideoRenderingSettings {
        frame_rate: video_settings.frame_rate,
        resolution_x: camera_settings.resolution_x,
        resolution_y: camera_settings.resolution_y,
        camera_diagonal: camera_settings.diagonal,
        camera_focal_length: camera_settings.focal_length,
        filepath_to_camera_path: video_settings.filepath_to_camera_path.to_path_buf(),
        filepath_to_background_image_1: filepath_bg_1,
        filepath_to_background_image_2: filepath_bg_2,
        filepath_to_output_folder: output_folder,
        output_video_name: video_settings.video_name.clone(),
        escape_radius: simulation_settings.escape_radius,
        max_iterations_propagation: simulation_settings.ray_integration_max_itarations,
        ray_integration_step: simulation_settings.ray_integration_step,
        alphas_num: simulation_settings.sampling_initial_nums,
        max_iterations_sampling: simulation_settings.sampling_initial_nums,
        sampling_convergence_threshold_1: simulation_settings.sampling_convergence_threshold_1,
        sampling_convergence_threshold_2: simulation_settings.sampling_convergence_threshold_2,
    })
}

/// Returns a single ImageRenderingSettings instance from the passed arguments
fn setup_image_rendering_system(
    filepath_bg_1: PathBuf,
    filepath_bg_2: PathBuf,
    output_folder: PathBuf,
    image_settings: &mut ImageSettings,
    camera_settings: &mut CameraSettings,
    simulation_settings: &mut SimulationSettings,
) -> Result<ImageRenderingSettings, String> {

    image_settings.normalize();
    image_settings.validate()?;

    camera_settings.normalize();
    camera_settings.validate()?;

    simulation_settings.normalize();
    simulation_settings.validate()?;

    let camera_pos = na::Vector4::new (
        image_settings.t,
        image_settings.l,
        image_settings.theta,
        image_settings.phi,
    );

    let camera_for = na::Vector3::new (
        image_settings.forward_x,
        image_settings.forward_y,
        image_settings.forward_z,
    );

    let camera_up = na::Vector3::new (
        image_settings.up_x,
        image_settings.up_y,
        image_settings.up_z,
    );

    Ok(ImageRenderingSettings {
        resolution_x: camera_settings.resolution_x,
        resolution_y: camera_settings.resolution_y,
        camera_diagonal: camera_settings.diagonal,
        camera_focal_length: camera_settings.focal_length,
        camera_position: camera_pos,
        camera_forward: camera_for,
        camera_up: camera_up,
        path_to_background_image_1: filepath_bg_1,
        path_to_background_image_2: filepath_bg_2,
        path_to_output_folder: output_folder,
        output_image_name: image_settings.image_name.clone(),
        escape_radius: simulation_settings.escape_radius,
        max_iterations_propagation: simulation_settings.ray_integration_max_itarations,
        ray_integration_step: simulation_settings.ray_integration_step,
        alphas_num: simulation_settings.sampling_initial_nums,
        max_iterations_sampling: simulation_settings.sampling_initial_nums,
        sampling_convergence_threshold_1: simulation_settings.sampling_convergence_threshold_1,
        sampling_convergence_threshold_2: simulation_settings.sampling_convergence_threshold_2,
    })
}


fn instantiate_metric(metric_settings: Either<EllisMetricSettings, InterstellarMetricSettings>) -> Metric {

    match metric_settings {
        Either::Left(metric_settings) => {
            let metric = EllisMetric::new(metric_settings.rho);
            
            return Metric::Ellis(metric);
        }
        Either::Right(metric_settings) => {
            let metric = InterstellarMetric::new(
                metric_settings.m,
                metric_settings.a,
                metric_settings.rho,
            );

            return Metric::Interstellar(metric);
        }
    }
}

/// Main script executed when calling `curvis video`
fn video_main() -> Result<(), String> {

    let (
        filepath_bg_1,
        filepath_bg_2,
        output_folder,
        mut video_settings,
        metric_settings,
        mut camera_settings,
        mut simulation_settings
    ) = parse_video_arguments();
    
    let metric = instantiate_metric(metric_settings);
    let video_rendering_settings = setup_video_rendering_settings(
        filepath_bg_1,
        filepath_bg_2,
        output_folder,
        &mut video_settings,
        &mut camera_settings,
        &mut simulation_settings)?;
    
    let outcome = match metric {
        Metric::Ellis(metric) => {
            let mut system = VideoRenderingSystem::new(metric, video_rendering_settings);
            system.render()
        }
        Metric::Interstellar(metric) => {
            let mut system = VideoRenderingSystem::new(metric, video_rendering_settings);
            system.render()
        }
    };

    return outcome;
}

/// Main script executed when calling `curvis image`
fn image_main() -> Result<(), String> {
    
    let (
        filepath_bg_1,
        filepath_bg_2,
        output_folder,
        mut image_settings,
         metric_settings,
         mut camera_settings,
         mut simulation_settings
    ) = parse_image_arguments();

    let metric = instantiate_metric(metric_settings);
    let image_rendering_settings = setup_image_rendering_system(
        filepath_bg_1,
        filepath_bg_2,
        output_folder,
        &mut image_settings,
        &mut camera_settings,
        &mut simulation_settings,
    )?;

    let outcome = match metric {
        Metric::Ellis(metric) => {
            let system = ImageRenderingSystem::new(metric, image_rendering_settings);
            system.render()
        }
        Metric::Interstellar(metric) => {
            let system = ImageRenderingSystem::new(metric, image_rendering_settings);
            system.render()
        }
    };

    return outcome;
}


fn main() {

    // env::set_var("RUST_BACKTRACE", "1");
    env::set_var("RUST_BACKTRACE", "full");

    // Dispatching based on the subcommand

    let subcommand = get_subcommand();

    if subcommand == "video" {
        println!("Video rendering");
        video_main().unwrap_or_else(|err| {eprintln!("Error in rendering video: {err}"); exit(1)});
    }
    else if subcommand == "image" {
        println!("Image rendering");
        image_main().unwrap_or_else(|err| {eprintln!("Error in rendering image: {err}"); exit(1)});
    }
    else if subcommand == "custom" {
        println!("Custom script");
        custom_main().unwrap_or_else(|err| {eprintln!("Error in curstom script: {err}"); exit(1)});
    }
    else {
        eprintln!("Unrecognized subcommand {subcommand:?}");
        exit(1);
    }

    exit(0);
}
