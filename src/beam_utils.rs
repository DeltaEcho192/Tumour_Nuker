use crate::mask::{Mask,  MaskHolder};
use crate::vector::Vector;
use rand::Rng;
use strum::IntoEnumIterator;
use strum_macros::EnumIter;

#[derive(Debug)]
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

#[derive(Debug)]
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
    pub dose_matrix: Box<[f32; N]>,
}

const BEAM_RADIUS: f32 = 1.5;
const E_DEPOSITED: f32 = 0.50;
const MU: f32 = 0.03;

pub fn compute_dose<const N: usize>(params: &mut ComputeDoseParams<{ N }>) {
    for beam_entry in &params.beams {
        let tumour_vector = beam_entry.beam_direction(&params.tumour);
        let self_dot = tumour_vector.dot(&tumour_vector);
        for x in 0..params.patient_box.x_size {
            for y in 0..params.patient_box.y_size {
                for z in 0..params.patient_box.z_size {
                    let mut dose_val: f32 = 0.0;
                    let mut vector = Vector::new(x as f32, y as f32, z as f32);
                    vector.calculate_offset(&beam_entry);
                    let dist = vector.dist_to_beam();
                    let dot_prod = vector.dot(&tumour_vector);
                    let projection_point = vector.mult_vec(dot_prod / self_dot);
                    let project_dist = vector.dist_to_vector(&projection_point);
                    if project_dist <= BEAM_RADIUS {
                        dose_val = E_DEPOSITED * (dist * -MU).exp();
                    }
                    let index: usize = (x + params.patient_box.y_size
                        * (y + params.patient_box.z_size * z))
                        as usize;
                    params.dose_matrix[index] += dose_val;
                }
            }
        }
    }
}

const WEIGHT_TUMOUR: f32 = 1.0;
const WEIGHT_SERIAL: f32 = 1.0;
const WEIGHT_PARALLEL: f32 = 1.0;
const WEIGHT_HEALTHY: f32 = 1.0;

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

    for x in 0..dose_params.patient_box.x_size {
        for y in 0..dose_params.patient_box.y_size {
            for z in 0..dose_params.patient_box.z_size {
                let index: usize = (x + dose_params.patient_box.y_size
                    * (y + dose_params.patient_box.z_size * z))
                    as usize;
                let mut mask_hit: bool = false;
                for mask in masks {
                    // Fix issue with multiple organs and the calculations
                    // being for per organ
                    if mask.bound_check(x, y, z) {
                        match mask.t_type {
                            TissueType::Tumour => {
                                mask_hit = true;
                                tumour_cost += dose_params.dose_matrix[index];
                            }
                            TissueType::SerialOrgan => {
                                mask_hit = true;
                                serial_oar_cost += (dose_params.dose_matrix[index] - D_THRESHOLD_S).max(0.0);
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
                    healthy_tissue_cost += (dose_params.dose_matrix[index] - D_THRESHOLD_H).max(0.0);
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

    println!("Tumour Cost: {}", tumour_cost);
    println!("Serial Cost: {}", serial_oar_cost);
    println!("Parallel Cost: {}", parallel_oar_cost);
    println!("Healthy Tissue Cost: {}", healthy_tissue_cost);


    let total_cost: f32 = WEIGHT_TUMOUR * tumour_cost
        + WEIGHT_SERIAL * serial_oar_cost
        + WEIGHT_PARALLEL * parallel_oar_cost
        + WEIGHT_HEALTHY * healthy_tissue_cost;

    println!("Total Cost: {}", total_cost);
    total_cost
}
