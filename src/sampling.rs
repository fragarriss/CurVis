/// In this module we implement a sampling algorithm that is used to evaluate
/// an expensive function more densely where it varies more rapidly.
/// In the context of this crate, the expensive function is the propagation function
/// that computes the trajectory of a photon and returns the escape angle and
/// the space where the photon has escaped (f64 -> (f64, f64)).
/// See `systems::render_image_efficient` for more details on how this algorithm
/// is used within the rendering process.
pub mod sampling {

    /// This structure actually contains the initial angle alpha, the eventual escape angle
    /// and a third float that is equal to either 1.0 or -1.0, depending on the
    /// space where the photon has escaped.
    #[derive(Debug, Copy, Clone)]
    struct BiPoint {
        a: f64, // stands for alpha
        e: f64, // stands for escape angle
        s: f64, // stands for sign or space
    }

    /// Returns the list with NAN and INFINITE values removed
    fn clean_bipoints(bipoints: Vec<BiPoint>) -> Vec<BiPoint> {

        let mut cleaned = Vec::new();

        for bipoint in bipoints {
            if bipoint.a.is_finite() && bipoint.e.is_finite() && bipoint.s.is_finite() {
                cleaned.push(bipoint);
            };
        }

        cleaned
    }

    /// Samples the expensive function such that more points are returned
    /// where the function varies more rapidly. See the source code for
    /// the implementation details.
    /// 
    /// # Arguments
    /// - `a_min`: the minimum value of the initial alpha angle range.
    /// - `a_max`: the maximum value of the initial alpha angle range.
    /// - `initial_points_number`: the number of initial range of alpha angles to evaluate.
    /// - `max_iterations`: the maximum number of iterations to perform.
    /// - `area_threshold_1`: the first area threshold parameter (for escape angles).
    /// - `area_threshold_2`: the second area threshold parameter (for escape spaces).
    /// - `expensive_function`: the expensive function to evaluate.
    pub fn doubly_sample_function<F>(
        a_min: f64,
        a_max: f64,
        initial_points_number: usize,
        max_iterations: usize,
        area_threshold_1: f64,
        area_threshold_2: f64,
        expensive_function: F,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>)
    where 
        F: Fn(f64) -> (f64, f64)
    {
        // At the beginning, we evaluate a range of uniformly spaced xs
        let uniform_xs = compute_uniform_range(
            a_min,
            a_max, 
            initial_points_number);
        
        // Then we compute the correponding bi-points
        let mut bipoints = uniform_xs.iter().map(|x| {
            let (y1, y2) = expensive_function(*x);
            BiPoint {
                a: *x,
                e: y1,
                s: y2,
            }
        }).collect::<Vec<BiPoint>>();

        bipoints = clean_bipoints(bipoints);


        // Now we are ready to iterate.
        // Each iteration computes the area of the triangle formed by three consecutive points
        // adding new points between if the area is larger than the threshold.

        // N.B. The current algorithm is not optimal and can be improved: in fact, here
        // we keep evaluating points that result in areas smaller than the threshold even if it was determined
        // to be the case in previous iterations.
        // Still, the area evaluation is way cheaper to compute than the expensive function, so
        // this is not a major issue.

        let mut amounts: Vec<usize> = Vec::new();
        amounts.push(bipoints.len());
        
        let mut iteration = 0;
        while iteration < max_iterations {

            let previous_amount = bipoints.len();
            
            bipoints = evaluate_denser_bipoints(bipoints, area_threshold_1, area_threshold_2, &expensive_function);
            amounts.push(bipoints.len());

            if bipoints.len() < previous_amount {
                // println!("Warning: sampling algorithm is failing. Iterations are reducing the number of points instead of increasing them (iteration {iteration}).");
                // let amounts_to_print: Vec<&usize> = amounts.iter().rev().take(10).rev().collect();
                // println!("Previous amounts: {amounts_to_print:?}");
                break;
            }

            // If no new points were added, we have reached the desired density
            if bipoints.len() == previous_amount {
                break;
            }

            iteration += 1;
        }

        // We issue a warning if the maximum number of iterations was reached
        if iteration == max_iterations {
            println!("Warning: maximum number of iterations ({max_iterations}) reached in sampling algorithm.");
        }

        // Finally, we unzip the bipoints to obtain the alphas, escape angles and signs
        let alphas = bipoints.iter().map(|bipoint| bipoint.a).collect::<Vec<f64>>();
        let escapes = bipoints.iter().map(|bipoint| bipoint.e).collect::<Vec<f64>>();
        let signs = bipoints.iter().map(|bipoint| bipoint.s).collect::<Vec<f64>>();

        (alphas, escapes, signs)
    }


    /// Returns a uniform range of alphas between `min_alpha` and `max_alpha`
    /// with `initial_points_number` points.
    fn compute_uniform_range(
        x_min: f64,
        x_max: f64,
        initial_points_number: usize,
    ) -> Vec<f64> {
        let mut uniform_xs = vec![];
        let step = (x_max - x_min) / ((initial_points_number-1) as f64);
        for i in 0..initial_points_number {
            uniform_xs.push(x_min + i as f64 * step);
        }
        uniform_xs
    }

    /// Adds points to the list of points passed if the area of the triangle formed
    /// by three consecutive points is larger than the area_threshold.
    fn evaluate_denser_bipoints<F>(
        bipoints:Vec<BiPoint>,
        threshold_parameter_1: f64,
        threshold_parameter_2: f64,
        expensive_function: &F) -> Vec<BiPoint>
    where
        F: Fn(f64) -> (f64, f64)
    {

        let mut new_bipoints: Vec<BiPoint> = vec![];
        let bipoints = clean_bipoints(bipoints);

        if bipoints.len() < 3 {
            panic!("bipoints list has length < 3. Cannot proceed to evaluate denser bipoints.")
        }

        let mut i: usize = 0;
        while i < (bipoints.len()-2) {

            let bipoint1 = &bipoints[i];
            let bipoint2 = &bipoints[(i+1) % bipoints.len()];
            let bipoint3 = &bipoints[(i+2) % bipoints.len()];

            let (score1, score2) = evaluate_convergence_scores(bipoint1, bipoint2, bipoint3);

            if !(score1 > threshold_parameter_1 || score2 > threshold_parameter_2) {
                // Within the bounds
                new_bipoints.push(BiPoint{a:bipoint1.a, e:bipoint1.e, s:bipoint1.s});
                i += 1;
                continue;
            }
            else {
                // We add new points between bipoint1 and bipoint2 and between bipoint2 and bipoint3
                let new_a_1 = (bipoint1.a + bipoint2.a) / 2.0;
                let new_a_2 = (bipoint2.a + bipoint3.a) / 2.0;

                let (new_e_1, new_s_1) = expensive_function(new_a_1);
                let (new_e_2, new_s_2) = expensive_function(new_a_2);

                new_bipoints.push(BiPoint{a:bipoint1.a, e:bipoint1.e, s:bipoint1.s});
                new_bipoints.push(BiPoint{a:new_a_1, e:new_e_1, s:new_s_1});
                new_bipoints.push(BiPoint{a:bipoint2.a, e:bipoint2.e, s:bipoint2.s});
                new_bipoints.push(BiPoint{a:new_a_2, e:new_e_2, s:new_s_2});
                
                // Skip the next iteration
                i += 2;
            }
        }

        clean_bipoints(new_bipoints)

    }

    /// Evaluates the convergence scores (area parameters) for three consecutive bipoints.
    fn evaluate_convergence_scores(bipoint1:&BiPoint, bipoint2:&BiPoint, bipoint3:&BiPoint) -> (f64, f64) {
        
        /* This is actually a placeholder implementation that could be improved.
    
        The points we are considering here have coordinate that live on the 0-2pi interval,
        thus we should take extra care in evaluating the area of the triangle formed by points
        living in such space.

        Another issue with this algorithm is that points are never added between
        the last and second-last point in the list.

        In terms of the final rendered image, if space-time varies rapidly
        close to the image (for instance because the wormhole edge is very close to the
        image border), the escape angle will vary rapidly in a point where alpha angles
        are the largest. In that situation, this alorithm will fail, introducing artifacts.

        The area of the triangle is evaluated using the Shoelace formula.
        The 1/2 factor is not included because we are not interested in the area
        per se, but how it compares to the threshold. */
    
        let area1 = (
            (
                bipoint1.a * bipoint2.e +
                bipoint2.a * bipoint3.e +
                bipoint3.a * bipoint1.e 
            ) -
            (
                bipoint1.e * bipoint2.a +
                bipoint2.e * bipoint3.a +
                bipoint3.e * bipoint1.a
            )
        ).abs();

        let area2 = (
            (
                bipoint1.a * bipoint2.s +
                bipoint2.a * bipoint3.s +
                bipoint3.a * bipoint1.s 
            ) -
            (
                bipoint1.s * bipoint2.a +
                bipoint2.s * bipoint3.a +
                bipoint3.s * bipoint1.a
            )
        ).abs();

        (area1, area2)
    }


}

