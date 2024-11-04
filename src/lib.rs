//! Curvis is a crate that let you use the equations of Einstein's General Relativity to render
//! image and videos of the curved space-time of wormholes.
//! 
//! It takes direct inspiration from this paper by O. James and collaborators (https://doi.org/10.1119/1.4916949It),
//! where you can find a detailed explanation of the equations that have been used within this project.
//! 
//! The project also was inspired by the great videos of Scott Manley, a Youtuber who already
//! explains very well what is going on with this code. I highly recommend a watch of those videos!


pub mod images;
// mod spherical_images;
pub mod algebra;
pub mod cameras;
pub mod metrics;
pub mod vectors;
pub mod systems;
pub mod csv;
pub mod interpolation;
pub mod rendering;
pub mod sampling;
pub mod settings;
pub mod cli;
pub mod filepaths;

pub mod custom;

// Re-exports

pub use crate::rendering::rendering::{ImageRenderingSystem, ImageRenderingSettings, VideoRenderingSystem, VideoRenderingSettings};
pub use crate::metrics::metrics::{DiagonalSphericalMetric, EllisMetric, InterstellarMetric};
pub use crate::settings::settings::{CameraSettings, VideoSettings, ImageSettings, InterstellarMetricSettings, EllisMetricSettings, SimulationSettings};
pub use crate::vectors::vectors::{RelativisticObject, RelativisticVector, Covariance};
pub use crate::images::spherical_images::{SphericalImage, load_image_as_spherical_image};
pub use crate::cameras::cameras::Camera;
pub use crate::algebra::algebra::Orientation;
pub use crate::systems::systems::{compute_photon_trajectory, compute_escape_angle};