use crate::mask::Mask;
use crate::vector::Vector;
use log::debug;
use rand::Rng;
use rayon::prelude::*;
use std::sync::{Arc, RwLock};
use std::thread;
use strum::IntoEnumIterator;
use strum_macros::EnumIter;

#[derive(Debug, Clone)]
pub struct PatientBox {
    pub x_size: i64,
    pub y_size: i64,
    pub z_size: i64,
}

impl PatientBox {
    pub const fn grid_size(&self) -> i64 {
        self.x_size * self.y_size * self.z_size
    }
}

#[derive(Debug, Clone)]
pub struct TissueBox {
    pub x: i64,
    pub y: i64,
    pub z: i64,
    pub x_width: i64,
    pub y_width: i64,
    pub z_width: i64,
    pub tissue_type: Option<TissueType>,
}

#[derive(Debug, Clone)]
pub enum TissueType {
    Tumour,
    SerialOrgan,
    ParallelOrgan,
}

#[derive(Debug, EnumIter)]
pub enum PatientBoxSide {
    LeftFace,
    RightFace,
    FrontFace,
    BackFace,
    BottomFace,
    TopFace,
}

pub fn compute_beam_entry(face: &PatientBoxSide, patient_box: &PatientBox) -> Vector {
    let mut rng = rand::rng();
    let entry_point = match face {
        PatientBoxSide::LeftFace => Vector::new(
            0.0,
            rng.random_range(0.0..(patient_box.y_size as f32)),
            rng.random_range(0.0..(patient_box.z_size as f32)),
        ),
        PatientBoxSide::RightFace => Vector::new(
            patient_box.x_size as f32,
            rng.random_range(0.0..(patient_box.y_size as f32)),
            rng.random_range(0.0..(patient_box.z_size as f32)),
        ),
        PatientBoxSide::FrontFace => Vector::new(
            rng.random_range(0.0..(patient_box.x_size as f32)),
            0.0,
            rng.random_range(0.0..(patient_box.z_size as f32)),
        ),
        PatientBoxSide::BackFace => Vector::new(
            rng.random_range(0.0..(patient_box.x_size as f32)),
            patient_box.y_size as f32,
            rng.random_range(0.0..(patient_box.z_size as f32)),
        ),
        PatientBoxSide::BottomFace => Vector::new(
            rng.random_range(0.0..(patient_box.x_size as f32)),
            rng.random_range(0.0..(patient_box.y_size as f32)),
            0.0,
        ),
        PatientBoxSide::TopFace => Vector::new(
            rng.random_range(0.0..(patient_box.x_size as f32)),
            rng.random_range(0.0..(patient_box.y_size as f32)),
            patient_box.z_size as f32,
        ),
    };
    entry_point
}

pub fn generate_beam_entries(patient_box: &PatientBox) -> Vec<Vector> {
    let mut beams: Vec<Vector> = vec![];
    for face in PatientBoxSide::iter() {
        beams.push(compute_beam_entry(&face, patient_box));
    }
    beams
}

pub struct ComputeDoseParams<const N: usize> {
    pub patient_box: PatientBox,
    pub beams: Vec<Vector>,
    pub tumour: TissueBox,
    pub dose_matrix: Arc<RwLock<Box<[f32; N]>>>,
}

pub struct ComputeDoseParamsIter<const N: usize> {
    pub patient_box: PatientBox,
    pub beams: Vec<Vector>,
    pub tumour: TissueBox,
    pub dose_matrix: Box<[f32; N]>,
}

const BEAM_RADIUS: f32 = 1.5;
const E_DEPOSITED: f32 = 0.50;
const MU: f32 = 0.03;

pub fn compute_dose<const N: usize>(params: &mut ComputeDoseParams<{ N }>) {
    let beams_vec = params.beams.clone();
    for beam_entry in beams_vec {
        let tumour_vector = beam_entry.beam_direction(&params.tumour.clone());
        let self_dot = tumour_vector.dot(&tumour_vector);
        let mut handle_vec = Vec::new();
        let local_ymax = params.patient_box.y_size;
        let local_zmax = params.patient_box.z_size;
        let local_xmax = params.patient_box.x_size;
        for x in 0..local_xmax {
            let dose_matrix_clone = Arc::clone(&params.dose_matrix);
            handle_vec.push(thread::spawn(move || {
                let mut local_dose: Vec<f32> = vec![];
                for y in 0..local_ymax {
                    for z in 0..local_zmax {
                        let mut dose_val: f32 = 0.0;
                        let test_idx = to_index(
                            x as usize,
                            y as usize,
                            z as usize,
                            local_xmax as usize,
                            local_ymax as usize,
                        );
                        let coords = to_coords(test_idx as i64, local_xmax, local_ymax);
                        let mut vector =
                            Vector::new(coords.0 as f32, coords.1 as f32, coords.2 as f32);
                        vector.calculate_offset(&beam_entry);
                        let dist = vector.dist_to_beam();
                        let dot_prod = vector.dot(&tumour_vector);
                        let projection_point = vector.mult_vec(dot_prod / self_dot);
                        let project_dist = vector.dist_to_vector(&projection_point);
                        if project_dist <= BEAM_RADIUS {
                            dose_val = E_DEPOSITED * (dist * -MU).exp();
                        }
                        local_dose.push(dose_val);
                    }
                }
                let mut w = dose_matrix_clone.write().unwrap();
                let mut counter = 0usize;
                for y in 0..local_ymax {
                    for z in 0..local_zmax {
                        let test_idx = to_index(
                            x as usize,
                            y as usize,
                            z as usize,
                            local_xmax as usize,
                            local_ymax as usize,
                        );
                        w[test_idx] += local_dose[counter];
                        counter += 1;
                    }
                }
            }));
        }

        for handle in handle_vec {
            handle.join().unwrap();
        }
    }
}

fn to_coords(index: i64, x_max: i64, y_max: i64) -> (i64, i64, i64) {
    let z = index / (x_max * y_max);
    let y = (index % (x_max * y_max)) / x_max;
    let x = index % x_max;
    (x, y, z)
}

fn to_index(x: usize, y: usize, z: usize, x_max: usize, y_max: usize) -> usize {
    x + y * x_max + z * x_max * y_max
}

pub fn compute_dose_iter<const N: usize>(params: &mut ComputeDoseParamsIter<{ N }>) {
    let beams_vec = params.beams.clone();
    for beam_entry in beams_vec {
        let tumour_vector = beam_entry.beam_direction(&params.tumour.clone());
        let self_dot = tumour_vector.dot(&tumour_vector);
        let local_ymax = params.patient_box.y_size;
        let local_xmax = params.patient_box.x_size;
        params
            .dose_matrix
            .par_iter_mut()
            .enumerate()
            .for_each(|(idx, val)| {
                let mapped_idx = to_coords(idx.try_into().unwrap(), local_xmax, local_ymax);
                let mut vector = Vector::new(
                    mapped_idx.0 as f32,
                    mapped_idx.1 as f32,
                    mapped_idx.2 as f32,
                );
                vector.calculate_offset(&beam_entry);
                let dist = vector.dist_to_beam();
                let dot_prod = vector.dot(&tumour_vector);
                let projection_point = vector.mult_vec(dot_prod / self_dot);
                let project_dist = vector.dist_to_vector(&projection_point);
                if project_dist <= BEAM_RADIUS {
                    *val += E_DEPOSITED * (dist * -MU).exp();
                }
            });
    }
}

const WEIGHT_TUMOUR: f32 = 1.0;
const WEIGHT_SERIAL: f32 = 1.0;
const WEIGHT_PARALLEL: f32 = 1.0;
const WEIGHT_HEALTHY: f32 = 0.75;

const D_PERSCRIBED: f32 = 40.0;
const D_THRESHOLD_S: f32 = 0.0;
const D_THRESHOLD_P: f32 = 0.0;
const D_THRESHOLD_H: f32 = 0.375;

pub fn compute_cost<const N: usize>(
    dose_params: &mut ComputeDoseParams<{ N }>,
    masks: &Vec<Mask>,
) -> f32 {
    let mut tumour_cost: f32 = 0.0;
    let mut serial_oar_cost: f32 = 0.0;
    let mut parallel_oar_cost: f32 = 0.0;
    let mut parallel_oar_intersections: i64 = 0;
    let mut healthy_tissue_cost: f32 = 0.0;

    let dose_maxtrix_read = dose_params.dose_matrix.read().unwrap();
    for x in 0..dose_params.patient_box.x_size {
        for y in 0..dose_params.patient_box.y_size {
            for z in 0..dose_params.patient_box.z_size {
                let index = to_index(
                    x as usize,
                    y as usize,
                    z as usize,
                    dose_params.patient_box.x_size as usize,
                    dose_params.patient_box.y_size as usize,
                );
                let mut mask_hit: bool = false;
                for mask in masks {
                    // Fix issue with multiple organs and the calculations
                    // being for per organ
                    if mask.bound_check(x, y, z) {
                        match mask.t_type {
                            TissueType::Tumour => {
                                //println!("Hit Tumour");
                                mask_hit = true;
                                tumour_cost += dose_maxtrix_read[index];
                            }
                            TissueType::SerialOrgan => {
                                mask_hit = true;
                                serial_oar_cost +=
                                    (dose_maxtrix_read[index] - D_THRESHOLD_S).max(0.0);
                            }
                            TissueType::ParallelOrgan => {
                                mask_hit = true;
                                parallel_oar_cost += dose_maxtrix_read[index];
                                parallel_oar_intersections += 1;
                            }
                        }
                    }
                }

                if !mask_hit {
                    healthy_tissue_cost += (dose_maxtrix_read[index] - D_THRESHOLD_H).max(0.0);
                }
            }
        }
    }

    if tumour_cost > 0.0 {
        tumour_cost = (tumour_cost - D_PERSCRIBED).abs();
    } else {
        tumour_cost = 1e6;
    }
    serial_oar_cost = (serial_oar_cost - D_THRESHOLD_S).max(0.0);

    let mean_dose: f32 = parallel_oar_cost / parallel_oar_intersections as f32;
    if mean_dose > D_THRESHOLD_P {
        parallel_oar_cost = mean_dose - D_THRESHOLD_P;
    } else {
        parallel_oar_cost = 0.0;
    }

    debug!("Tumour Cost: {}", tumour_cost);
    debug!("Serial Cost: {}", serial_oar_cost);
    debug!("Parallel Cost: {}", parallel_oar_cost);
    debug!("Healthy Tissue Cost: {}", healthy_tissue_cost);

    let total_cost: f32 = WEIGHT_TUMOUR * tumour_cost
        + WEIGHT_SERIAL * serial_oar_cost
        + WEIGHT_PARALLEL * parallel_oar_cost
        + WEIGHT_HEALTHY * healthy_tissue_cost;

    debug!("Total Cost: {}", total_cost);
    total_cost
}

pub fn compute_cost_iter<const N: usize>(
    dose_params: &mut ComputeDoseParamsIter<{ N }>,
    masks: &Vec<Mask>,
) -> f32 {
    let mut tumour_cost: f32 = 0.0;
    let mut serial_oar_cost: f32 = 0.0;
    let mut parallel_oar_cost: f32 = 0.0;
    let mut parallel_oar_intersections: i64 = 0;
    let mut healthy_tissue_cost: f32 = 0.0;

    for x in 0..dose_params.patient_box.x_size {
        for y in 0..dose_params.patient_box.y_size {
            for z in 0..dose_params.patient_box.z_size {
                let index: usize = to_index(
                    x as usize,
                    y as usize,
                    z as usize,
                    dose_params.patient_box.x_size as usize,
                    dose_params.patient_box.y_size as usize,
                );
                let mut mask_hit: bool = false;
                for mask in masks {
                    // Fix issue with multiple organs and the calculations
                    // being for per organ
                    if mask.bound_check(x, y, z) {
                        match mask.t_type {
                            TissueType::Tumour => {
                                //println!("Hit Tumour");
                                mask_hit = true;
                                tumour_cost += dose_params.dose_matrix[index];
                            }
                            TissueType::SerialOrgan => {
                                mask_hit = true;
                                serial_oar_cost +=
                                    (dose_params.dose_matrix[index] - D_THRESHOLD_S).max(0.0);
                            }
                            TissueType::ParallelOrgan => {
                                mask_hit = true;
                                parallel_oar_cost += dose_params.dose_matrix[index];
                                parallel_oar_intersections += 1;
                            }
                        }
                    }
                }

                if !mask_hit {
                    healthy_tissue_cost +=
                        (dose_params.dose_matrix[index] - D_THRESHOLD_H).max(0.0);
                }
            }
        }
    }

    if tumour_cost > 0.0 {
        tumour_cost = (tumour_cost - D_PERSCRIBED).abs();
    } else {
        tumour_cost = 1e6;
    }
    serial_oar_cost = (serial_oar_cost - D_THRESHOLD_S).max(0.0);

    let mean_dose: f32 = parallel_oar_cost / parallel_oar_intersections as f32;
    if mean_dose > D_THRESHOLD_P {
        parallel_oar_cost = mean_dose - D_THRESHOLD_P;
    } else {
        parallel_oar_cost = 0.0;
    }

    debug!("Tumour Cost: {}", tumour_cost);
    debug!("Serial Cost: {}", serial_oar_cost);
    debug!("Parallel Cost: {}", parallel_oar_cost);
    debug!("Healthy Tissue Cost: {}", healthy_tissue_cost);

    let total_cost: f32 = WEIGHT_TUMOUR * tumour_cost
        + WEIGHT_SERIAL * serial_oar_cost
        + WEIGHT_PARALLEL * parallel_oar_cost
        + WEIGHT_HEALTHY * healthy_tissue_cost;

    debug!("Total Cost: {}", total_cost);
    total_cost
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Instant;

    #[test]
    fn test_cost_function() {
        const PATIENT: PatientBox = PatientBox {
            x_size: 100,
            y_size: 50,
            z_size: 30,
        };

        let tumour = TissueBox {
            x: 40,
            y: 35,
            z: 12,
            x_width: 5,
            y_width: 5,
            z_width: 5,
            tissue_type: Some(TissueType::Tumour),
        };

        let serial_organ = TissueBox {
            x: 42,
            y: 38,
            z: 18,
            x_width: 3,
            y_width: 3,
            z_width: 3,
            tissue_type: Some(TissueType::SerialOrgan),
        };

        let parallel_organ = TissueBox {
            x: 50,
            y: 20,
            z: 12,
            x_width: 30,
            y_width: 15,
            z_width: 10,
            tissue_type: Some(TissueType::ParallelOrgan),
        };

        let mut mask_holder: Vec<Mask> = vec![];
        mask_holder.push(Mask::from_tissue_box(&tumour, &PATIENT));
        mask_holder.push(Mask::from_tissue_box(&serial_organ, &PATIENT));
        mask_holder.push(Mask::from_tissue_box(&parallel_organ, &PATIENT));
        let mut beams: Vec<Vector> = vec![];
        for face in PatientBoxSide::iter() {
            let entry_point = match face {
                PatientBoxSide::LeftFace => Vector::new(0.0, 25.0, 15.0),
                PatientBoxSide::RightFace => Vector::new(PATIENT.x_size as f32, 20.0, 12.3),
                PatientBoxSide::FrontFace => Vector::new(64.6, 0.0, 1.2),
                PatientBoxSide::BackFace => Vector::new(24.3, PATIENT.y_size as f32, 22.0),
                PatientBoxSide::BottomFace => Vector::new(98.9, 34.2, 0.0),
                PatientBoxSide::TopFace => Vector::new(46.0, 44.2, PATIENT.z_size as f32),
            };
            beams.push(entry_point);
        }

        const N_SIZE: usize = PATIENT.grid_size() as usize;
        println!(
            "Rough Memory Size of Dose Matrix: {} MB",
            (N_SIZE * std::mem::size_of::<[f32; 1]>()) / 1024 / 1024
        );

        let now = Instant::now();
        let mut dose_params: ComputeDoseParams<{ N_SIZE }> = ComputeDoseParams {
            patient_box: PATIENT.clone(),
            beams: beams.clone(),
            tumour: tumour.clone(),
            dose_matrix: Arc::new(RwLock::new(vec![0f32; N_SIZE].try_into().unwrap())),
        };
        compute_dose(&mut dose_params);
        let fitness = compute_cost(&mut dose_params, &mask_holder);
        println!("Fitness Preset: {}", fitness);
        assert_eq!(fitness, 118.73086);
        println!(
            "Time Taken Compute cost and total: {} Miliseconds",
            now.elapsed().as_millis()
        );

        let now2 = Instant::now();
        let mut dose_params_iter: ComputeDoseParamsIter<{ N_SIZE }> = ComputeDoseParamsIter {
            patient_box: PATIENT.clone(),
            beams: beams.clone(),
            tumour: tumour.clone(),
            dose_matrix: vec![0f32; N_SIZE].try_into().unwrap(),
        };

        compute_dose_iter(&mut dose_params_iter);
        let fitness = compute_cost_iter(&mut dose_params_iter, &mask_holder);
        println!("Fitness Preset: {}", fitness);
        assert_eq!(fitness, 118.73086);
        println!(
            "Time Taken Compute cost and total for iter version: {} Miliseconds",
            now2.elapsed().as_millis()
        );
    }

    #[test]
    fn test_coordinate_conversion() {
        let x_max = 5i64;
        let y_max = 6i64;
        let _z_max = 5i64;

        let index1 = to_index(0, 0, 0, x_max as usize, y_max as usize);
        let index2 = to_index(0, 1, 0, x_max as usize, y_max as usize);
        let index3 = to_index(3, 4, 2, x_max as usize, y_max as usize);

        assert_eq!((0, 0, 0), to_coords(index1 as i64, x_max, y_max));
        assert_eq!((0, 1, 0), to_coords(index2 as i64, x_max, y_max));
        assert_eq!((3, 4, 2), to_coords(index3 as i64, x_max, y_max));
    }
}
