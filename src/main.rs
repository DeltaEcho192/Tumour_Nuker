use tumour_nuker::beam_utils::{
    PatientBox, TissueBox, TissueType,
};
use tumour_nuker::ga::ga;
use tumour_nuker::mask::Mask;
use std::time::Instant;

fn main() {
    println!("Running Tumour Nuker Optimizer");
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

    const N_SIZE: usize = PATIENT.grid_size() as usize;
    println!(
        "Rough Memory Size of Dose Matrix: {} MB",
        (N_SIZE * std::mem::size_of::<[f32; 1]>()) / 1024 / 1024
    );

    let now = Instant::now();
    ga::<{ N_SIZE }>(20, 10, PATIENT, tumour, mask_holder, 5);
    println!(
        "Time Taken Compute cost and total: {} Miliseconds",
        now.elapsed().as_millis()
    );
}
