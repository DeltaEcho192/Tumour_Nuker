use Tumour_Nuker::beam_utils::{compute_dose, compute_cost, generate_beam_entries, ComputeDoseParams, PatientBox, TissueBox, TissueType};
use Tumour_Nuker::mask::Mask;
use std::time::Instant;
use Tumour_Nuker::vector::Vector;
use std::thread;

fn loop_test() {
    println!("Run Loop test");
    //    const grid_x: i32 = 200;
    //    const grid_y: i32 = 100;
    //    const grid_z: i32 = 80;
    //    const total_size: usize = (grid_x * grid_y * grid_z) as usize;
    //
    //    let mut inbeam_count: i64 = 0;
    //    let mut operationcounter: i64 = 0;
    //    let mut handle_vec = Vec::new();
    //    let beam_amt: i32 = 12;
    //    for i in 1..beam_amt {
    //        handle_vec.push(thread::spawn(|| {
    //            let mut dose: Vec<f32> = Vec::with_capacity(total_size);
    //            let beam_vec: Vector = Vector::new(0.8, 1.43, 6.3);
    //
    //            let beam_ent_x: i32 = 2;
    //            let beam_ent_y: i32 = 2;
    //            let beam_ent_z: i32 = 0;
    //
    //            let self_dot = beam_vec.dot(&beam_vec);
    //
    //            let mut inbeam_count: i64 = 0;
    //            let mut operationcounter: i64 = 0;
    //            for x in 0..grid_x {
    //                for y in 0..grid_y {
    //                    for z in 0..grid_z {
    //                        let mut dose_val: f32 = 0.0;
    //                        let mut vector = Vector::new(x as f32, y as f32, z as f32);
    //                        vector.calculate_offset(beam_ent_x, beam_ent_y, beam_ent_z);
    //                        let dist = vector.dist_to_beam();
    //                        let dot_prod = vector.dot(&beam_vec);
    //                        let projection_point = vector.mult_vec(dot_prod / self_dot);
    //                        let project_dist = vector.dist_to_vector(&projection_point);
    //                        if project_dist <= 0.75 {
    //                            inbeam_count += 1;
    //                            dose_val = project_dist;
    //                        }
    //                        //let index: usize = (x + grid_y * (y + grid_z * z)) as usize;
    //                        dose.push(dose_val);
    //                        operationcounter += 1;
    //                    }
    //                }
    //            }
    //        }));
    //    }
    //
    //    for handle in handle_vec {
    //        handle.join().unwrap();
    //    }
}

fn main() {
    println!("Running Tumour Nuker Optimizer");
    const PATIENT: PatientBox = PatientBox {
        x_size: 200,
        y_size: 100,
        z_size: 40,
    };

    let tumour = TissueBox {
        x: 35,
        y: 40,
        z: 12,
        x_width: 5,
        y_width: 5,
        z_width: 5,
        tissue_type: Some(TissueType::Tumour),
    };

    let serial_organ = TissueBox {
        x: 38,
        y: 42,
        z: 18,
        x_width: 3,
        y_width: 3,
        z_width: 3,
        tissue_type: Some(TissueType::SerialOrgan),
    };

    let parallel_organ = TissueBox {
        x: 20,
        y: 50,
        z: 12,
        x_width: 15,
        y_width: 30,
        z_width: 10,
        tissue_type: Some(TissueType::ParallelOrgan),
    };

    let mut mask_holder: Vec<Mask> = vec![]; 
    mask_holder.push(Mask::from_tissue_box(&tumour, &PATIENT));
    mask_holder.push(Mask::from_tissue_box(&serial_organ, &PATIENT));
    mask_holder.push(Mask::from_tissue_box(&parallel_organ, &PATIENT));

    let beams_i = generate_beam_entries(&PATIENT);

    const N_SIZE: usize = PATIENT.grid_size() as usize;
    let mut dose_params: ComputeDoseParams<{N_SIZE}> = ComputeDoseParams {
        patient_box: PATIENT,
        beams: beams_i,
        tumour: tumour,
        dose_matrix: vec![0f32; N_SIZE].try_into().unwrap(),
    };

    println!("Rough Memory Size of Dose Matrix: {} MB", (N_SIZE * std::mem::size_of::<[f32; 1]>())/1024/1024);

    let now = Instant::now();
    compute_dose(&mut dose_params);
    println!("Time Taken Compute Dose {} Miliseconds", now.elapsed().as_millis());
    println!();
    compute_cost(&mut dose_params, &mask_holder);
    println!("Time Taken Compute cost and total: {} Miliseconds", now.elapsed().as_millis());
}

