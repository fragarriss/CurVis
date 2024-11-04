// This module contains code for reading a camera path from a csv file.

// Reads csv file and returns a list of camera positions.

// The csv format is as follows:
// First line is the header (ignored in the present implementation); then the
// following lines are comma-separated float values in the following order:
// x0,x1,x2,x3, fx,fy,fz, ux, uy, uz
// where (x0,x1,x2,x3) are the relativistic coordinates of the camera on the
// curved space, (fx,fy,fz) and (ux,uy,uz) are the forward and up vectors that
// define the camera orientation in the tangent space in (x0,x1,x2,x3).

pub mod csv {

    use nalgebra as na;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::PathBuf;

    /// Loads a csv file and returns the three lists corresponding to the camera
    /// relativistic positions, forward vectors and up vectors.
    /// At this stage, information is returned as a tuple of three lists of na:Vector3<f64>
    /// and na:Vector4<f64> objects.
    pub fn load_path(path_to_csv_file: &PathBuf) -> (Vec<na::Vector4<f64>>, Vec<na::Vector3<f64>>, Vec<na::Vector3<f64>>) {
        let mut positions: Vec<na::Vector4<f64>> = Vec::new();
        let mut forward_vectors: Vec<na::Vector3<f64>> = Vec::new();
        let mut up_vectors: Vec<na::Vector3<f64>> = Vec::new();

        let file = File::open(path_to_csv_file).expect("Could not open file");
        let reader = BufReader::new(file);

        for (index, line) in reader.lines().enumerate() {

            // Skipping the header
            if index == 0 {
                continue;
            }

            let line = line.expect("Could not read line");
            let mut values = line.split(",").map(|x| x.parse::<f64>().expect("Could not parse float"));

            let x0 = values.next().expect("Could not read x0");
            let x1 = values.next().expect("Could not read x1");
            let x2 = values.next().expect("Could not read x2");
            let x3 = values.next().expect("Could not read x3");

            let fx = values.next().expect("Could not read fx");
            let fy = values.next().expect("Could not read fy");
            let fz = values.next().expect("Could not read fz");

            let ux = values.next().expect("Could not read ux");
            let uy = values.next().expect("Could not read uy");
            let uz = values.next().expect("Could not read uz");

            positions.push(na::Vector4::new(x0, x1, x2, x3));
            forward_vectors.push(na::Vector3::new(fx, fy, fz));
            up_vectors.push(na::Vector3::new(ux, uy, uz));
        }

        (positions, forward_vectors, up_vectors)

    }

}