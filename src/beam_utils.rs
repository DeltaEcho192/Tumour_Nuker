use crate::vector::Vector;
use rand::Rng;

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

#[derive(Debug)]
pub enum TissueType {
    Tumour,
    SerialOrgan,
    ParallelOrgan,
    HealthyTissue,
}

#[derive(Debug)]
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


pub struct ComputeDoseParams {
    pub patient_box: PatientBox,
    pub beams: Vec<Vector>,
}

pub fn compute_dose(params: &ComputeDoseParams) {

    for beam in params.beams {
        let self_dot = beam.dot(&beam);

        for x in 0..params.patient_box.x_size {
            for y in 0..params.patient_box.y_size {
                for z in 0..params.patient_box.z_size {
                    let mut dose_val: f32 = 0.0;
                    let mut vector = Vector::new(x as f32, y as f32, z as f32);
                    vector.calculate_offset(&beam);
                    let dist = vector.dist_to_beam();
                    let dot_prod = vector.dot(&beam);
                    let projection_point = vector.mult_vec(dot_prod / self_dot);
                    let project_dist = vector.dist_to_vector(&projection_point);
                    if project_dist <= 0.75 {
                        dose_val = project_dist;
                    }
                    //let index: usize = (x + grid_y * (y + grid_z * z)) as usize;
                    dose.push(dose_val);
                }
            }
        }
    }
}
